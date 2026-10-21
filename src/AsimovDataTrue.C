/*
Author: Yael Deniz
Phd Candidate, University of Idaho

Description:
Class for a simulation of the interacting rates and interacting events.
It implements a custom 3D model of the Earth  (definitions in Earth3D.h)
3D model includes LLPVs.

Possible bin variables are: E,logE ; \theta,\cos(\theta); \phi 

Fucntions: 
    GetIntRates computes the neutrino rates Rint -> Integrate to get Nint

IMPORTANT:

All angular inputs are in degrees (main files)

All angular variables are in radias (Here)


Convention:

    ibins - Total number of Energy [GeV] bins
    jbins - Total number of zenith("theta") [rad] bins | cos(theta)
    kbins - Total number of azimuth [rad] bins 



Notes:
3D code has a "Break memory violation" if theta is 90 \pm 0.1 -> Fixed with a patch
that changes to oscprob prem for that range in zenith. A potential issue appears
if we consider the relative position of the LLVPs-detector 

Neutrino flux is interpolated in the space log10(dPhi/dEdcth) and (log10E, cos(\theta))
    ** To obtain dPhi/dEdcth, do pow(log10(dPhi/dEdcth), 10) to the interpolation results
    **IF using \theta as angular variable, rember to adjust the result of the bilinear interpolation into dPhi/dEdth =  sth * dPhi/dEdcth

GetTrueEvents3Dthlog10E
    Energy bins: Even space bins in log10(E)
    Zenith bins: Even space bins in  \theta


*/
#include "GetTrueEvents.h"

//C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//Random number generator
#include <random>


//Cern Root
#include <math.h>
#include "TF1.h"
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include <Math/Interpolator.h>
#include "TGraph.h"
//#include "TGraphEventsD.h"

//OSCPROB
#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif

// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"
#include "Earth3DModel.h"
//Layer Updates
#include "LayerUpdate.h"


using namespace std;


// Some Constants i need
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds
# define mN 1.67492749804E-27 // Approximate mass of nucleon

void AsimovSimulation::ChangeLayerProperties(int n,
                                             double rhopct,
                                             double zoapct)
{
    // ChangeLayerProperties
    // ---------------------
    // Registers or updates a modification for a given PREM layer
    // in the AsimovSimulation internal list.
    //
    // Parameters:
    //   n       : 1-based index of the layer (row) in the PREM CSV file
    //   rhopct  : relative change in density (%) for this layer
    //   zoapct  : relative change in Z/A (%) for this layer
    //
    // The updates are stored in LayerUpdates_intsim as LayerUpdate structs.
    // They will later be forwarded to Earth3DModel inside GetIntEvents().

    // If this layer already has a modification, update it instead of duplicating
    for (auto& upd : LayerUpdates_intsim) {
        if (upd.tableNumber == n) {              // n is 1-based
            upd.densityContrast  = rhopct;
            upd.chemicalContrast = zoapct;
            return;
        }
    }

    // Otherwise add a new modification
    LayerUpdate NewUpdate;
    NewUpdate.tableNumber      = n;      // keep 1-based
    NewUpdate.densityContrast  = rhopct;
    NewUpdate.chemicalContrast = zoapct;

    LayerUpdates_intsim.push_back(NewUpdate);
}


void AsimovSimulation::AddUpdatesList(const std::vector<LayerUpdate>& updates)
{
    // AddUpdatesList
    // --------------
    // Registers a batch of PREM layer property updates in the AsimovSimulation.
    //
    // Parameters
    // ----------
    //   updates : vector of LayerUpdate objects. Each entry corresponds to a
    //             PREM layer (1-based row index in the CSV) and contains:
    //               - tableNumber      : 1-based index of the PREM layer
    //               - densityContrast  : relative change in density (%)
    //               - chemicalContrast : relative change in Z/A (%)
    //
    // For each element, this method forwards the values to ChangeLayerProperties(),
    // which stores or updates the corresponding entry in LayerUpdates_intsim.
    // These updates will later be passed to Earth3DModel in GetIntEvents().

    if (updates.empty()) {
        Error("AddUpdatesList",
              "No layer updates were provided. At least one LayerUpdate "
              "is required to modify the PREM model.");
        return;
    }

    for (const auto& upd : updates) {
        ChangeLayerProperties(
            upd.tableNumber,      // 1-based index, same as CSV row
            upd.densityContrast,  // %
            upd.chemicalContrast  // %
        );
    }
}




// Generates the interacting rates for a given final flavour and MH
TH3D*AsimovSimulation::GetIntRates() //Binning fixed to be in th
{  
    // Define paths
    std::string NuTomoPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder  = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string AziIndex   = NuTomoPath + IntFolder + "TrueIndexTable.txt"; 

    // Final neutrino flavours
    int nue  = 0;  // Electron neutrino  
    int numu = 1;  // Muon neutrino

    // Particle types
    int nu    =  1;  // Neutrino
    int nubar = -1;  // Antineutrino

    double R_ijk = 0;  // Rate: Poisson mean for bin (i, j, k)
    double N_ijk = 0;  // Events: Poisson mean for bin (i, j, k)

    // Define histogram (X=energy, Y=zenith, Z=azimuth)
    TH3D* IntHist = new TH3D("IntHist", "Interacting Rate Histogram", 
                             ibins, Emin, Emax,
                             jbins, thmin, thmax,
                             kbins, phimin, phimax); 

    // Initialize PMNS for oscillation probabilities
    OscProb::PMNS_Fast PMNS_H;
    PMNS_H.SetStdPars();  // Standard 3-flavor oscillation parameters

    // Load Honda flux data
    NuFlux HondaFlux;

    std::vector<std::vector<double>> FluxData = HondaFlux.SetFluxData(HondaTable);

    // Generate flux histograms by flavor
    TH2D* muflux  = HondaFlux.GetFluxHist(1, FluxData);  // Muon neutrino
    TH2D* mubflux = HondaFlux.GetFluxHist(2, FluxData);  // Muon antineutrino
    TH2D* eflux   = HondaFlux.GetFluxHist(3, FluxData);  // Electron neutrino
    TH2D* ebflux  = HondaFlux.GetFluxHist(4, FluxData);  // Electron antineutrino

    // Set up 3D Earth model
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable);
    Earth3D.SetDetector(pos);
    Earth3D.PileThickness = PileHeight;
    Earth3D.aWidth = aperture;
    Earth3D.SetPile(MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
    //Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);

    // If we have updates, register them in the Earth model
    if (!LayerUpdates_intsim.empty()) {
        Earth3D.AddUpdatesList(LayerUpdates_intsim);
    } 




    // Create default PREM Model
    OscProb::PremModel prem_oscprob(PremTable);

    double l, d, z; //Path length, Path density, Path Z/A

    // Loop over azimuth bins
    for (int k = 1; k <= IntHist->GetNbinsZ(); ++k) 
    {
        double phi  = IntHist->GetZaxis()->GetBinCenter(k); // Bin center azimuth
        //Wraping around the circle
        if( phi < 0.0 ) phi = phi + 2*TMath::Pi();
        if(phi > 2*TMath::Pi()) phi = phi - 2*TMath::Pi();
        
        double dphi = IntHist->GetZaxis()->GetBinWidth(k);  // Bin width azimuth
        double checkphi = phi;

        // Loop over zenith bins 
        for (int j = 1; j <= IntHist->GetNbinsY(); ++j) 
        {
            double th  = IntHist->GetYaxis()->GetBinCenter(j); // Bin center zenith
            double cth  = cos(th);
            double sth  = sin(th); 

            double dth = IntHist->GetYaxis()->GetBinWidth(j); //  Bin width zenith 

            if (cth < -1 || cth > 1) continue;  // Skip unphysical values

            // Currently, my 3D code cannot generate path when \zenith = \pi/2 \pm \eps
            // A fix for this is to use PREM from OSCPROB in that region; "Model Patch"

            if( th < (91.0)*(TMath::Pi()/180.0) && th > (89.0)*(TMath::Pi()/180.0))
            {
                prem_oscprob.SetDetPos(abs(pos[2])); // Detector Position
               
                prem_oscprob.FillPath(cth); // Set the nu-direction and calculate the total path

                PMNS_H.SetPath(prem_oscprob.GetNuPath()); // Set the total neutrino path

                std::cout << "Using Oscprob to calculate horizontal paths" << std::endl;
            }//Calculate path using OscProb models; End
            

            else //If neutrino direction is not in the PATCH region use my 3D model
            { 
                // Set direction in Earth model and compute 3D path
                Earth3D.SetDirection(cth, phi);
                std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

                // Initialize oscillation path with first segment
                l = EarthPath[0][0]; //Path lenght
                d = EarthPath[0][1]; //Path density
                z = EarthPath[0][2]; //Path Z/A
                PMNS_H.SetPath(l, d, z);

                // Add remaining segments
                for (int np = 1; np < EarthPath.size(); ++np) 
                { 
                    l = EarthPath[np][0];
                    d = EarthPath[np][1];
                    z = EarthPath[np][2];    
                    PMNS_H.AddPath(l, d, z);
                }


            }//Calculate path using 3D model; End

            // Compute numerical rates for neutrinos

            // Loop over energy bins
            for (int i = 1; i <= IntHist->GetNbinsX(); ++i)
            {
                double e     = IntHist->GetXaxis()->GetBinCenter(i); // Bin center energy
                double dE    = IntHist->GetXaxis()->GetBinWidth(i);  // Bin width energy
                double logEi = log10(e); //Log10(e); used during interpolation of Neutrino flux

                // Interpolate logarithmic flux values
                double logdPsiMu  = muflux->Interpolate(logEi, cth);
                double logdPsiMub = mubflux->Interpolate(logEi, cth);
                double logdPsiE   = eflux->Interpolate(logEi, cth);
                double logdPsiEb  = ebflux->Interpolate(logEi, cth);

                // Convert to linear flux; factor of "sth" is added to change variables from cos(\theta) to \theta
                double dPsiMudEdth    = sth*pow(10, logdPsiMu);
                double dPsiMubardEdth = sth*pow(10, logdPsiMub);
                double dPsiEdEdth     = sth*pow(10, logdPsiE);
                double dPsiEbardEdth  = sth*pow(10, logdPsiEb);

                // Neutrino contribution
                PMNS_H.SetIsNuBar(false);
                double Rijk_e  = XSec(e, nu)*(PMNS_H.Prob(nue, flvf, e)* dPsiEdEdth); //Electron neutrino contribution 
                double Rijk_mu = XSec(e, nu)*(PMNS_H.Prob(numu, flvf, e)* dPsiMudEdth);//Muon neutrino contribution 
                double Rijk_nu = Rijk_e + Rijk_mu; //Total neutrino contribution 

                // Antineutrino contribution
                PMNS_H.SetIsNuBar(true);
                double Rijk_eb  = XSec(e, nubar)*(PMNS_H.Prob(nue, flvf, e)*dPsiEbardEdth); //Electron neutrino contribution
                double Rijk_mub = XSec(e, nubar)*(PMNS_H.Prob(numu, flvf, e)*dPsiMubardEdth);//Muon neutrino contribution
                double Rijk_nubar = Rijk_eb + Rijk_mub; //Total neutrino contribution

                //Total interaction rate
                R_ijk =  (Rijk_nu + Rijk_nubar);
                //Total interaction rate
                //N_ijk =  (MT/mN)*(R_ijk)*(dE*dth*dphi);
                
                //std::cout <<"Energy:" << e <<" width "<< dE << " Bin " << i << ":"<< IntHist->GetNbinsZ() <<  std::endl;
                //std::cout <<"Zenith:" << th*(180.0)/TMath::Pi()<< "Width" << dth << " Bin" << j << ":"<< IntHist->GetNbinsY() <<  std::endl;
                //std::cout <<"Azimuth:" << phi*(180.0)/TMath::Pi() << " Width" << dphi << " Bin" << k << ":"<< IntHist->GetNbinsZ() <<  std::endl;
                //std::cout << "Rate: "<< R_ijk <<  " Events " << N_ijk << std::endl;

                // Fill histogram
                IntHist->SetBinContent(i, j, k, R_ijk);
            }  // End energy loop
        }  // End zenith loop
    }  // End azimuth loop

    //TH3D* intRates = (TH3D*) IntHist->Clone("Interacting Rates Cache");

    return IntHist;
}

// Generates the event distribution from a rate distribution:
TH3D* AsimovSimulation::IntEvents(TH3D* ratesHist) {
   
    if (!ratesHist) return nullptr;

    TH3D* eventsHist = (TH3D*) ratesHist->Clone("InteractingEvents");
    eventsHist->SetTitle("Interacting Events Histogram");

    for (int i = 1; i <= ratesHist->GetNbinsX(); ++i) {
        double dE = ratesHist->GetXaxis()->GetBinWidth(i);
        for (int j = 1; j <= ratesHist->GetNbinsY(); ++j) {
            double dth = ratesHist->GetYaxis()->GetBinWidth(j);
            for (int k = 1; k <= ratesHist->GetNbinsZ(); ++k) {
                double dphi = ratesHist->GetZaxis()->GetBinWidth(k);

                double Rijk = ratesHist->GetBinContent(i, j, k);
                double Nijk = (MT/mN)*Rijk * dE * dth * dphi;

                eventsHist->SetBinContent(i, j, k, Nijk);
            }
        }
    }

    return eventsHist;
}

// Acces event rates do:
TH3D* AsimovSimulation::GetIntEvents() {
    TH3D* rateHist = GetIntRates();    // Step 1: Call your rate-calculation method
    return IntEvents(rateHist);        // Step 2: Convert rates to events
}

//Avaraged Interacting distribution over Azimuth (Useful for Marginalized vMF)
TH2D* AsimovSimulation::GetAzimuthAvaraged( TH3D* Nint_hist)
{
    TH2D* Nint_histavg = new TH2D("TrueHist_avg", "Azimuth-Averaged Interacting Rate",
                        Nint_hist->GetNbinsX(), Nint_hist->GetXaxis()->GetXbins()->GetArray(),
                        Nint_hist->GetNbinsY(), Nint_hist->GetYaxis()->GetXbins()->GetArray());

    double Nbins = static_cast<double>(Nint_hist->GetNbinsZ()); // Preferred C++ way
 
    for (int i = 1; i <= Nint_hist->GetNbinsX(); ++i) 
    {
        for (int j = 1; j <= Nint_hist->GetNbinsY(); ++j) 
        {
            double sum = 0.0;

            for (int k = 1; k <= Nint_hist->GetNbinsZ(); ++k) 
            {
                sum += Nint_hist->GetBinContent(i, j, k);
            }

            double avg = sum / Nbins;
            Nint_histavg->SetBinContent(i, j, avg);
        }
    }
return Nint_histavg;
}



TH3D*AsimovSimulation::GetIntSignal() //Binning fixed to be in th
{  
    // Define paths
    std::string NuTomoPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder  = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string AziIndex   = NuTomoPath + IntFolder + "TrueIndexTable.txt";
    std::ofstream IndexAzi(AziIndex); 

    std::cout << "Not here 1" << std::endl;

    // Final neutrino flavours
    int nue  = 0;  // Electron neutrino  
    int numu = 1;  // Muon neutrino

    // Particle types
    int nu    =  1;  // Neutrino
    int nubar = -1;  // Antineutrino

    double N_ijk = 0;  // Poisson mean for bin (i, j, k)

    // Define histogram (energy, cos(theta), azimuth)
    TH3D* IntHist = new TH3D("IntHist", "Interacting Signal Histogram", 
                             ibins, Emin, Emax,
                             jbins, thmin, thmax,
                             kbins, phimin, phimax); 

    // Initialize PMNS for oscillation probabilities
    OscProb::PMNS_Fast PMNS_H;
    PMNS_H.SetStdPars();  // Standard 3-flavor oscillation parameters

    // Load Honda flux data
    NuFlux HondaFlux;
    std::vector<std::vector<double>> FluxData = HondaFlux.SetFluxData(HondaTable);

    std::cout << "Not here 3" << std::endl;

    // Generate flux histograms by flavor
    TH2D* muflux  = HondaFlux.GetFluxHist(1, FluxData);  // Muon neutrino
    TH2D* mubflux = HondaFlux.GetFluxHist(2, FluxData);  // Muon antineutrino
    TH2D* eflux   = HondaFlux.GetFluxHist(3, FluxData);  // Electron neutrino
    TH2D* ebflux  = HondaFlux.GetFluxHist(4, FluxData);  // Electron antineutrino

    // Set up 3D Earth model
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable);
    Earth3D.SetDetector(pos);
    Earth3D.PileThickness = PileHeight;
    Earth3D.aWidth = aperture;
    Earth3D.SetPile(MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
    Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);

    // Create default PREM Model
    OscProb::PremModel prem_oscprob(PremTable);

    double l, d, z;

    // Loop over azimuth bins
    for (int k = 1; k <= kbins; ++k) 
    {
        double phi  = IntHist->GetZaxis()->GetBinCenter(k);
        double dphi = IntHist->GetZaxis()->GetBinWidth(k);
        double checkphi = phi;

        //Wraping around the circle
        if( phi < 0.0 ) phi = phi + 2*TMath::Pi();

        if(phi > 2*TMath::Pi()) phi = phi - 2*TMath::Pi();

        //std::cout << "True phi" << checkphi*(180.0/TMath::Pi()) << " " << phi*(180.0/TMath::Pi()) << std::endl;


        // Loop over zenith bins (cos(theta))
        for (int j = 1; j <= jbins; ++j) 
        {
            double th  = IntHist->GetYaxis()->GetBinCenter(j);
            double cth  = cos(th);
            double sth  = sin(th); 
            double dth = IntHist->GetYaxis()->GetBinWidth(j);

            if (cth < -1 || cth > 1) continue;  // Skip unphysical values

            
            if( th < (91.0)*(TMath::Pi()/180.0) && th > (89.0)*(TMath::Pi()/180.0))
            {
                prem_oscprob.SetDetPos( abs(pos[2]));
                prem_oscprob.FillPath(cth);
                PMNS_H.SetPath(prem_oscprob.GetNuPath());
            }
            

            else{ // --------------------------------------------Set Direction Using 3D model
                // Set direction in Earth model and compute 3D path
                Earth3D.SetDirection(cth, phi);
                std::cout << "Direction of the path: phi = " << phi*180.0/TMath::Pi() << " th"<< th*180.0/TMath::Pi() << ", cth = " << cth << std::endl;
                std::cout << "bins phi = " << kbins << " th"<< jbins <<  std::endl;

                std::cout << "min th" << IntHist->GetYaxis()->GetXmin()*180.0/TMath::Pi()<< " Max th" << IntHist->GetYaxis()->GetXmax()*180.0/TMath::Pi() << std::endl; // Minimum X
                

                std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();
                std::cout << "Path is the problem" << std::endl;

                // Initialize oscillation path with first segment
                l = EarthPath[0][0];
                d = EarthPath[0][1];
                z = EarthPath[0][2];
                PMNS_H.SetPath(l, d, z);

                // Add remaining segments
                for (int np = 1; np < EarthPath.size(); ++np) 
                { 
                    l = EarthPath[np][0];
                    d = EarthPath[np][1];
                    z = EarthPath[np][2];    
                    PMNS_H.AddPath(l, d, z);
                }// --------------------------------------- Set direction using 3D model
            }               
            // Loop over energy bins
            for (int i = 1; i <= ibins; ++i)
            {
                double e     = IntHist->GetXaxis()->GetBinCenter(i);
                double dE    = IntHist->GetXaxis()->GetBinWidth(i);
                double logEi = log10(e);

                // Interpolate logarithmic flux values
                double logdPsiMu  = muflux->Interpolate(logEi, cth);
                double logdPsiMub = mubflux->Interpolate(logEi, cth);
                double logdPsiE   = eflux->Interpolate(logEi, cth);
                double logdPsiEb  = ebflux->Interpolate(logEi, cth);

                // Convert to linear flux
                double dPsiMudEdth    = sth*pow(10, logdPsiMu);
                double dPsiMubardEdth = sth*pow(10, logdPsiMub);
                double dPsiEdEdth     = sth*pow(10, logdPsiE);
                double dPsiEbardEdth  = sth*pow(10, logdPsiEb);

                // Neutrino contribution
                PMNS_H.SetIsNuBar(false);
                double Ri_e  = XSec(e, nu)    * (PMNS_H.Prob(nue, flvf, e)  * dPsiEdEdth);
                double Ri_mu = XSec(e, nu)    * (PMNS_H.Prob(numu, flvf, e) * dPsiMudEdth);
                double Ri_nu = Ri_e + Ri_mu;

                // Antineutrino contribution
                PMNS_H.SetIsNuBar(true);
                double Ri_eb  = XSec(e, nubar) * (PMNS_H.Prob(nue, flvf, e)  * dPsiEbardEdth);
                double Ri_mub = XSec(e, nubar) * (PMNS_H.Prob(numu, flvf, e) * dPsiMubardEdth);
                double Ri_nubar = Ri_eb + Ri_mub;

                // Final expected interaction rate
                N_ijk = (MT / mN) * (Ri_nu + Ri_nubar) * dE * dth * dphi;

                // Fill histogram
                IntHist->SetBinContent(i, j, k, N_ijk);
            }  // End energy loop
        }  // End zenith loop
    }  // End azimuth loop

    return IntHist;
}

//Simulation with binning in cos(theta)
std::vector< TH2D* > AsimovSimulation::GetTrueEvents3D()
{  
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string AziIndex = NuTomoPath+IntFolder+"TrueIndexTable.txt";
    std::ofstream IndexAzi(AziIndex); 

    std::cout << "Not here 1" << std::endl;

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    //Binning Scheme 
    double cthmin = cos(thmax); // Min value for cth (cth = -1) 
    double cthmax = cos(thmin); // Max value for cth (cth = 1)
    //ibins - Total number of zenith("theta") bins | cos(theta)
    //jbins - Total number of azimuth bins
    //kbins - Total number of Energy bins 
    double N_ijk = 0   ;        // Poisson mean for bin (i,j,k).

    // Histogram definition (3D: cos(theta), azimuth, energy)
    TH3D * EventHist3D = new TH3D("TrueHist","True Event Histrogram", 
                                  ibins,cthmin,cthmax,
                                  jbins,phimin,phimax,
                                  kbins,Emin,Emax); 

    std::cout << "Not here 2" << std::endl;

    //Neutrino Oscillation Probabilities calculation
    OscProb::PMNS_Fast PMNS_H; // PMNS object
    PMNS_H.SetStdPars();       // Set PDG 3-flavor parameters

    //Load Honda flux data
    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);

    std::cout << "Not here 3" << std::endl;


    //Generate flux histograms for different neutrino flavours
    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); // Muon neutrino
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData);// Muon antineutrino
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData);  // Electron neutrino
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); // Electron antineutrino


    std::cout << "Not here 4" << std::endl;

    // Set up the Earth's 3D model    
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable);
    Earth3D.SetDetector(pos);
    Earth3D.PileThickness = PileHeight;  
    Earth3D.aWidth = aperture;
    Earth3D.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
    Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);

    std::cout << "Not here 3" << PremTable <<  std::endl;


    
    // Initialize vector for 2D histograms and other variables
    std::vector< TH2D* > HistVec;
    TH2D * EventHist2D[jbins];
    double l,d,z;

    //Loop over azimuth angles
    for (int j = 1; j <= jbins; j++) 
    {

        std::cout << "Not here 4 " << j <<  std::endl;


        // Get azimuth center and width for the current bin
        double phi = EventHist3D->GetYaxis()->GetBinCenter(j); 
        double dphi = EventHist3D-> GetYaxis()->GetBinWidth(j); 

        IndexAzi << j << " " << phi*180/TMath::Pi() << " " << dphi*180/TMath::Pi() << std::endl; 

        // Create a new 2D histogram for each azimuth bin
        //char histchar[64];
        //strcat(histchar, "truehist%d" , j);
        //const char *histname = histchar;

        std::cout << "Not here 5 IndexTable is not a problem " << std::endl;

        TString histname = Form("OscHist%d",j);

        std::cout << "  histname is note the problem " << std::endl;

        EventHist2D[j] = new TH2D(histname,Form("Oscillogram%d",j), ibins,cthmin,cthmax,kbins,Emin,Emax); 

        std::cout << "  Vector is the problem  " << std::endl;


        // Loop over zenith bins (cos(theta))
        for(int i=1; i<= ibins ; i++) 
        {    
            // Get cos(theta) center and width for the current bin.
            // The distance traveled through the earth "Baseline" is L=-2 * R_earth * cos(theta)

            double cth = EventHist2D[j]->GetXaxis()->GetBinCenter(i);
            double dcth = EventHist2D[j] -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip unphysical cos(theta)

            // Set Earth model direction and calculate Earth path
            
            Earth3D.SetDirection(cth,phi);
            std::cout << " Direction of the path  " << phi << " " << cth  << std::endl;
 


            std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();
            std::cout << " Path is the problem  " << std::endl;



            // Set PMNS neutrino oscillation path
            l = EarthPath[0][0];
            d = EarthPath[0][1];
            z = EarthPath[0][2];
            
            PMNS_H.SetPath(l,d,z); // Initialize path
            
            for (int np = 1; np < EarthPath.size(); np++) 
            { 
                l = EarthPath[np][0];
                d = EarthPath[np][1];
                z = EarthPath[np][2];    
                PMNS_H.AddPath(l,d,z);
            } 
            
            // Loop over energy bins
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = EventHist2D[j]->GetYaxis()->GetBinCenter(k); 
                double dE = EventHist2D[j]->GetYaxis()->GetBinWidth(k);
                double logEi = log10(e);
                

                // Interpolate neutrino flux values
                double logdPsiMu = muflux->Interpolate(logEi,cth);
                double logdPsiMub = mubflux->Interpolate(logEi,cth);
                double logdPsiE = eflux->Interpolate(logEi,cth);
                double logdPsiEb = ebflux->Interpolate(logEi,cth);

                // Convert logarithmic flux values to linear
                double dPsiMudEdct = pow(10,logdPsiMu);     // Muon neutrino flux
                double dPsiMubardEdct = pow(10,logdPsiMub); // Muon anti-neutrino flux
                double dPsiEdEdct = pow(10,logdPsiE);       // Electron neutrino flux
                double dPsiEbardEdct = pow(10,logdPsiEb);   // Electron anti-neutrino flux

                PMNS_H.SetIsNuBar(false);// neutrino
                double Ri_e = XSec(e,nu)*( PMNS_H.Prob(nue,flvf,e)*dPsiEdEdct); 
                double  Ri_mu = XSec(e,nu)*(PMNS_H.Prob(numu, flvf, e)*dPsiMudEdct); 
                double Ri_nu = Ri_e + Ri_mu;
                
                PMNS_H.SetIsNuBar(true); //Antineutrino 
                double Ri_eb=XSec(e,nubar)*(PMNS_H.Prob(nue,flvf,e)*dPsiEbardEdct ); 
                double Ri_mub=XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardEdct );
                double Ri_nubar = Ri_eb + Ri_mub;

                //double  mN = 1.67492749804E-27; // Approximate mass of nucleon
                // Calculate Poisson mean for the current bins (N_ijk)
                // Note : MT/mN  refers to target nucleons per year
                double N_ijk = (MT/mN)*(Ri_nu + Ri_nubar)*dE*dcth*dphi; 

                // Set the bin content in the 2D histogram
                EventHist2D[j]->SetBinContent(i,k, N_ijk); 

            
            } // End energy loop 

        } // End zenith loop

        // Add the current 2D histogram to the vector

        std::cout << "All good for phi " << j <<  std::endl;

        HistVec.push_back(EventHist2D[j]);

    } //End  Azimuth Loop 


    // Close index file and return the vector of histograms      
    IndexAzi.close();
    return HistVec;
}

//Simulation with binning in theta
std::vector< TH2D* > AsimovSimulation::GetTrueEvents3Dth()
{  
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string AziIndex = NuTomoPath+IntFolder+"TrueIndexTable.txt";
    std::ofstream IndexAzi(AziIndex); 

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    //Binning Scheme 

    //double cthmin = cos(thmax); // Min value for cth (cth = -1) 
    //double cthmax = cos(thmin); // Max value for cth (cth = 1)
    //ibins - Total number of zenith("theta") bins | cos(theta)
    //jbins - Total number of azimuth bins
    //kbins - Total number of Energy bins 
    double N_ijk = 0   ;        // Poisson mean for bin (i,j,k).

    // Histogram definition (3D: cos(theta), azimuth, energy)
    TH3D * EventHist3D = new TH3D("TrueHist","True Event Histrogram", 
                                  ibins,thmin,thmax,
                                  jbins,phimin,phimax,
                                  kbins,Emin,Emax); 

    //Neutrino Oscillation Probabilities calculation
    OscProb::PMNS_Fast PMNS_H; // PMNS object
    PMNS_H.SetStdPars();       // Set PDG 3-flavor parameters

    //Load Honda flux data
    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);

    //Generate flux histograms for different neutrino flavours
    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); // Muon neutrino
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData);// Muon antineutrino
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData);  // Electron neutrino
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); // Electron antineutrino

    // Set up the Earth's 3D model    
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable);
    Earth3D.SetDetector(pos);
    Earth3D.PileThickness = PileHeight;  
    Earth3D.aWidth = aperture;
    Earth3D.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
    Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
    // Initialize vector for 2D histograms and other variables
    std::vector< TH2D* > HistVec;
    TH2D * EventHist2D[jbins];
    double l,d,z;

    //Loop over azimuth angles
    for (int j = 1; j <= jbins; j++) 
    {
        // Get azimuth center and width for the current bin
        double phi = EventHist3D->GetYaxis()->GetBinCenter(j); 
        double dphi = EventHist3D-> GetYaxis()->GetBinWidth(j); 

        IndexAzi << j << " " << phi*180/TMath::Pi() << " " << dphi*180/TMath::Pi() << std::endl; 

        // Create a new 2D histogram for each azimuth bin
        //char histchar[64];
        //strcat(histchar, "truehist%d" , j);
        //const char *histname = histchar;

        TString histname = Form("OscHist%d",j);
        EventHist2D[j] = new TH2D(histname,Form("Oscillogram%d",j), ibins,thmin,thmax,kbins,Emin,Emax); 

        // Loop over zenith bins (cos(theta))
        for(int i=1; i<= ibins ; i++) 
        {    
            // Get cos(theta) center and width for the current bin.
            // The distance traveled through the earth "Baseline" is L=-2 * R_earth * cos(theta)

            double th = EventHist2D[j]->GetXaxis()->GetBinCenter(i);
            double cth = cos(th);
            double sth = sin(th);
            double dth = EventHist2D[j] -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip unphysical cos(theta)

            // Set Earth model direction and calculate Earth path
            
            Earth3D.SetDirection(cth,phi); 
            std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

            // Set PMNS neutrino oscillation path
            l = EarthPath[0][0];
            d = EarthPath[0][1];
            z = EarthPath[0][2];
            
            PMNS_H.SetPath(l,d,z); // Initialize path
            
            for (int np = 1; np < EarthPath.size(); np++) 
            { 
                l = EarthPath[np][0];
                d = EarthPath[np][1];
                z = EarthPath[np][2];    
                PMNS_H.AddPath(l,d,z);
            } 
            
            // Loop over energy bins
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = EventHist2D[j]->GetYaxis()->GetBinCenter(k); 
                double dE = EventHist2D[j]->GetYaxis()->GetBinWidth(k);
                double logEi = log10(e);
                

                // Interpolate neutrino flux values
                double logdPsiMu = muflux->Interpolate(logEi,cth);
                double logdPsiMub = mubflux->Interpolate(logEi,cth);
                double logdPsiE = eflux->Interpolate(logEi,cth);
                double logdPsiEb = ebflux->Interpolate(logEi,cth);

                // Convert logarithmic flux values to linear
                double dPsiMudEdth    = sth*pow(10,logdPsiMu);     // Muon neutrino flux
                double dPsiMubardEdth = sth*pow(10,logdPsiMub); // Muon anti-neutrino flux
                double dPsiEdEdth     = sth*pow(10,logdPsiE);       // Electron neutrino flux
                double dPsiEbardEdth  = sth*pow(10,logdPsiEb);   // Electron anti-neutrino flux

                PMNS_H.SetIsNuBar(false);// neutrino
                double Ri_e = XSec(e,nu)*( PMNS_H.Prob(nue,flvf,e)*dPsiEdEdth); 
                double  Ri_mu = XSec(e,nu)*(PMNS_H.Prob(numu, flvf, e)*dPsiMudEdth); 
                double Ri_nu = Ri_e + Ri_mu;
                
                PMNS_H.SetIsNuBar(true); //Antineutrino 
                double Ri_eb=XSec(e,nubar)*(PMNS_H.Prob(nue,flvf,e)*dPsiEbardEdth ); 
                double Ri_mub=XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardEdth );
                double Ri_nubar = Ri_eb + Ri_mub;

                //double  mN = 1.67492749804E-27; // Approximate mass of nucleon
                // Calculate Poisson mean for the current bins (N_ijk)
                // Note : MT/mN  refers to target nucleons per year
                double N_ijk = (MT/mN)*(Ri_nu + Ri_nubar)*dE*dth*dphi; 

                // Set the bin content in the 2D histogram
                EventHist2D[j]->SetBinContent(i,k, N_ijk); 

            } // End energy loop 

        } // End zenith loop

        // Add the current 2D histogram to the vector

        HistVec.push_back(EventHist2D[j]);

    } //End  Azimuth Loop 


    // Close index file and return the vector of histograms      
    IndexAzi.close();
    return HistVec;
}

//Simulation with binning in theta
std::vector< TH2D* > AsimovSimulation::GetTrueEvents3Dthlog10E()
{  

    /*
    Note for the reader
    Energy bins: Even space bins in log10(E)
    Zenith bins: Even space bins in  \theta

    Neutrino flux is interpolated in the space log10(dPhi/dEdcth) and (log10E, cos(\theta)) 
    */
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string AziIndex = NuTomoPath+IntFolder+"TrueIndexTable.txt";
    std::ofstream IndexAzi(AziIndex); 

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    //Binning Scheme 

    //double cthmin = cos(thmax); // Min value for cth (cth = -1) 
    //double cthmax = cos(thmin); // Max value for cth (cth = 1)
    //ibins - Total number of zenith("theta") bins | cos(theta)
    //jbins - Total number of azimuth bins
    //kbins - Total number of Energy bins 
    double N_ijk = 0   ;        // Poisson mean for bin (i,j,k).

    // Histogram definition (3D: cos(theta), azimuth, energy)
    TH3D * EventHist3D = new TH3D("TrueHist","True Event Histrogram", 
                                  ibins,thmin,thmax,
                                  jbins,phimin,phimax,
                                  kbins,log10(Emin),log10(Emax)); 

    //Neutrino Oscillation Probabilities calculation
    OscProb::PMNS_Fast PMNS_H; // PMNS object
    PMNS_H.SetStdPars();       // Set PDG 3-flavor parameters

    //Load Honda flux data
    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);

    //Generate flux histograms for different neutrino flavours
    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); // Muon neutrino
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData);// Muon antineutrino
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData);  // Electron neutrino
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); // Electron antineutrino

    // Set up the Earth's 3D model    
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable);
    Earth3D.SetDetector(pos);
    Earth3D.PileThickness = PileHeight;  
    Earth3D.aWidth = aperture;
    Earth3D.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
    Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
    // Initialize vector for 2D histograms and other variables
    std::vector< TH2D* > HistVec;
    TH2D * EventHist2D[jbins];
    double l,d,z;

    //Loop over azimuth angles
    for (int j = 1; j <= jbins; j++) 
    {
        // Get azimuth center and width for the current bin
        double phi = EventHist3D->GetYaxis()->GetBinCenter(j); 
        double dphi = EventHist3D-> GetYaxis()->GetBinWidth(j); 

        IndexAzi << j << " " << phi*180/TMath::Pi() << " " << dphi*180/TMath::Pi() << std::endl; 

        // Create a new 2D histogram for each azimuth bin
        //char histchar[64];
        //strcat(histchar, "truehist%d" , j);
        //const char *histname = histchar;

        TString histname = Form("OscHist%d",j);
        EventHist2D[j] = new TH2D(histname,Form("Oscillogram%d",j), ibins,thmin,thmax,kbins,log10(Emin),log10(Emax) ); 

        // Loop over zenith bins (cos(theta))
        for(int i=1; i<= ibins ; i++) 
        {    
            // Get cos(theta) center and width for the current bin.
            // The distance traveled through the earth "Baseline" is L=-2 * R_earth * cos(theta)

            double th = EventHist2D[j]->GetXaxis()->GetBinCenter(i);
            double cth = cos(th);
            double sth = sin(th);
            double dth = EventHist2D[j] -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip unphysical cos(theta)

            // Set Earth model direction and calculate Earth path
            
            Earth3D.SetDirection(cth,phi); 
            std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

            // Set PMNS neutrino oscillation path
            l = EarthPath[0][0];
            d = EarthPath[0][1];
            z = EarthPath[0][2];
            
            PMNS_H.SetPath(l,d,z); // Initialize path
            
            for (int np = 1; np < EarthPath.size(); np++) 
            { 
                l = EarthPath[np][0];
                d = EarthPath[np][1];
                z = EarthPath[np][2];    
                PMNS_H.AddPath(l,d,z);
            } 
            
            // Loop over energy bins
            for (int k = 1; k <=kbins; ++k)
            { 
                double log10e = EventHist2D[j]->GetYaxis()->GetBinCenter(k); 
                double e = pow(10, log10e);
                double dlog10E = EventHist2D[j]->GetYaxis()->GetBinWidth(k);
                //double log10Ei = log10(e);
                

                // Interpolate neutrino flux values
                double logdPsiMu = muflux->Interpolate(log10e,cth);
                double logdPsiMub = mubflux->Interpolate(log10e,cth);
                double logdPsiE = eflux->Interpolate(log10e,cth);
                double logdPsiEb = ebflux->Interpolate(log10e,cth);

                // Convert logarithmic flux values to linear
                double dPsiMudlog10Edth    = log(10)*e*sth*pow(10,logdPsiMu);     // Muon neutrino flux
                double dPsiMubardlog10Edth = log(10)*e*sth*pow(10,logdPsiMub); // Muon anti-neutrino flux
                double dPsiEdlog10Edth     = log(10)*e*sth*pow(10,logdPsiE);       // Electron neutrino flux
                double dPsiEbardlog10Edth  = log(10)*e*sth*pow(10,logdPsiEb);   // Electron anti-neutrino flux

                PMNS_H.SetIsNuBar(false);// neutrino
                double Ri_e = XSec(e,nu)*( PMNS_H.Prob(nue,flvf,e)*dPsiEdlog10Edth); 
                double  Ri_mu = XSec(e,nu)*(PMNS_H.Prob(numu, flvf, e)*dPsiMudlog10Edth); 
                double Ri_nu = Ri_e + Ri_mu;
                
                PMNS_H.SetIsNuBar(true); //Antineutrino 
                double Ri_eb=XSec(e,nubar)*(PMNS_H.Prob(nue,flvf,e)*dPsiEbardlog10Edth ); 
                double Ri_mub=XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardlog10Edth );
                double Ri_nubar = Ri_eb + Ri_mub;

                //double  mN = 1.67492749804E-27; // Approximate mass of nucleon
                // Calculate Poisson mean for the current bins (N_ijk)
                // Note : MT/mN  refers to target nucleons per year
                double N_ijk = (MT/mN)*(Ri_nu + Ri_nubar)*dlog10E*dth*dphi; 

                // Set the bin content in the 2D histogram
                EventHist2D[j]->SetBinContent(i,k, N_ijk); 

            } // End energy loop 

        } // End zenith loop

        // Add the current 2D histogram to the vector

        HistVec.push_back(EventHist2D[j]);

    } //End  Azimuth Loop 


    // Close index file and return the vector of histograms      
    IndexAzi.close();
    return HistVec;
}


std::vector< TH2D* > AsimovSimulation::GetTrueEvents2D( ) //To be Deleted
{   
    //Define cos(theta) range and azimuthal angle range
    double cthmin = cos(thmax); // min cos(theta) = -1 
    double cthmax = cos(thmin); // max cos(theta) = 0
    double dphi = phimax - phimin; //Azimutham angle range

    //Bin labels
    // ibins - Bins in Zenith
    // jbins - Bins in Azimuth
    // kbins - Bins in Energy

    //Initialize histograms for neutrino events
    TH2D * EventHist = new TH2D("EventHist","True Event Histrogram", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtseflv = new TH2D("Evtseflv","#nu_{e} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtseflvb = new TH2D("Evtseflvb","#bar{#nu}_{e} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtsmuflv = new TH2D("Evtsmuflv","#nu_{#mu} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax);
    TH2D * Evtsmuflvb = new TH2D("Evtsmuflvb","#bar{#nu}_{#mu} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
                            
     std::vector < TH2D * >  Histflv; //Store Oscillograms for (muflv, muflvbar, eflv, eflvbar, flv+ flvbar)
     
     TH2D * Hist[5];
    
    //Neutrino properties

    // Neutrino initial flavour
    int nue        = 0;  // Electron neutrino  
    int numu       = 1;  // Muon neutrino
    // Particle type
    int nu         = 1;  // Neutrino
    int nubar      = -1; // Antineutrino

    // PMNS matrix setup
    OscProb::PMNS_Fast PMNS_H; 
    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters

    // Load Honda flux data and create histograms for different neutrino types
    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);
    //Generate flux histograms for different neutrino flavours
    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); // Muon neutrino
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData);// Muon antineutrino
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData);  // Electron neutrino
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); // Electron antineutrino

    OscProb::PremModel prem(PremTable);


    // Set Earth model and Detector    
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable); // Earth model
    Earth3D.SetDetector(pos);    // Detector position
    // Modify Specific layer properties
    Earth3D.SetLayerProp(PremTableNumber, 
                         DensityContrast, ChemicalContrast); 
    // Set heterogeneities in lower mantle (e.g., LLVPs)
    Earth3D.PileThickness = PileHeight;  // Height of LLVP
    Earth3D.aWidth = aperture;           // Aperture of LLVP
    Earth3D.SetPile( MantleAnomaly, 
                   AnomalyShape, PileDensityContrast, PileChemContrast);
        
        // Loop over zenith bins (cos(theta))
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    
            // Get the center and width of a bin in cos(theta)
            double cth  = EventHist->GetXaxis()->GetBinCenter(i);
            double dcth = EventHist -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip unphysical values 
            
            // Set the path through the Earth for direction (cos(theta),phi)
            prem.FillPath(cth); // Fill paths from PREM model
            prem.FillPath(cth); // Fill paths from PREM model
            PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
            //Earth3D.SetDirection(cth,0);  
            //std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

            // Distance through "ith" layer : l - EarthPath[i][0];
            // Densty of "ith" layer : rho - EarthPath[i][2];
            // z/a of "ith" layer : zoa - EarthPath[i][3];
            
            // Set initial Earth path in the PMNS model
            //PMNS_H.SetPath(EarthPath[0][0],EarthPath[0][1],EarthPath[0][2]);
            
            // Add subsequent segments of the path
            //for (int i = 1; i < EarthPath.size(); i++) 
            //{ 
              //PMNS_H.AddPath(EarthPath[i][0],EarthPath[i][1],EarthPath[i][2]);
            //} 

            // Loop through energy bins 
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = EventHist->GetYaxis()->GetBinCenter(k); 
                double dE = EventHist->GetYaxis()->GetBinWidth(k);
                double logEi = log10(e);
                
                //Interpolate neutrino fluxes
                 double logdPsiMu = muflux->Interpolate(logEi,cth);
                 double logdPsiMub = mubflux->Interpolate(logEi,cth);
                 double logdPsiE = eflux->Interpolate(logEi,cth);
                 double logdPsiEb = ebflux->Interpolate(logEi,cth);

                //Evaluate differential Flux
                 double dPsiMudEdct = pow(10,logdPsiMu);     //Muon neutrino flux
                 double dPsiMubardEdct = pow(10,logdPsiMub); //Muon anti-neutrino flux
                 double dPsiEdEdct = pow(10,logdPsiE);        //Electron neutrino flux
                 double dPsiEbardEdct = pow(10,logdPsiEb);    //Electron anti-neutrino flux


                //Neutrino Contribution;
                PMNS_H.SetIsNuBar(false); // Set neutrino probabilities 
                // iflv -> fflv oscillation probabilities
                double Pef = PMNS_H.Prob(nue,flvf,e);   // nue -> nufinal
                double Pmuf = PMNS_H.Prob(numu,flvf, e);// numu -> nufinal

                // Neutrino Interacting rate
                double R_ef = XSec(e,nu)*( Pef*dPsiEdEdct);
                double  R_muf = XSec(e,nu)*(Pmuf*dPsiMudEdct); 
                double R_f = R_ef + R_muf; // Events in the final-neutrino channel

                //Antineutrino contribution
                PMNS_H.SetIsNuBar(true); // Set antineutrino probabilities 
                // iflv -> fflv oscillation probabilities
                double Pefb = PMNS_H.Prob(nue,flvf,e);   // nuebar -> nufinalbar
                double Pmufb = PMNS_H.Prob(numu,flvf, e);// numubar -> nufinalbar
                
                // Antineutrino Interacting rate
                double R_efb=XSec(e,nubar)*(  Pefb*dPsiEbardEdct ); 
                double R_mufb=XSec(e,nubar)*( Pmufb *dPsiMubardEdct ); 
                double R_fbar = R_efb + R_mufb; //// Events in the final-antineutrino channel

                // Total Event rate
                double R = R_f+R_fbar;

                //double  mN = 1.67492749804E-27; // Approximate mass of nucleon
                //double  Meff = 14.6*pow(log10(e),1.8 )*MTon;
                double    Meff = 1*MTon;
                //double  N_A = 6.02214E23;
                double  T = 1*years2sec;
                //double  N_f = (2*TMath::Pi())*(Meff*T/mN)*R*dE*dcth;  // Number of events for bin (cth,e)

                  
                //Note : MT/mN  refers to target nucleons per year
                //double  N_f = (2*TMath::Pi())*(MT/mN)*R*dE*dcth;  // Number of events for bin (cth,e)
                double  N_f = (e*e)*(2*TMath::Pi())*(Meff*T/mN)*R*dcth;  // Number of events for bin (cth,e) 

                Evtseflv->SetBinContent(i,k,  ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_ef );
                Evtseflvb->SetBinContent(i,k, ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_efb );
                Evtsmuflv->SetBinContent(i,k, ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_muf );
                Evtsmuflvb->SetBinContent(i,k,((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_mufb );
                EventHist->SetBinContent(i,k, N_f ); //Create histogram for  kth Pseudo-Experimens

            } // End of loop in energy

            // Clear PMNS path for the next iteration
            PMNS_H.ClearPath();

        }   //   End of loop in zenith

    //Store all histrogram in a vector
    
    Histflv.push_back(Evtsmuflv); // Muon neutrino contribution to events
    Histflv.push_back(Evtsmuflvb);// Muon antineutrino contribution to events
    Histflv.push_back(Evtseflv);  // Electron neutrino contribution to events
    Histflv.push_back(Evtseflvb); // Electron antineutrino contribution to events
    Histflv.push_back(EventHist); // flvf neutrino events

    // Simulation summary output
    std::cout << "Detector location" << pos[2] << std::endl;
    std::cout << "PREM model: " << PremTable  <<std::endl;
    std::cout << "Flux Table: " << HondaTable  <<std::endl;
    std::cout << "Pile| h:" <<  PileHeight << " aWidth:" << aperture 
              << " Shape: " << AnomalyShape << " rhopct: " << PileDensityContrast << std::endl;
    std::cout <<" Simulation Summary| flvf:" << flvf << " IsNuBar:" <<  (nunubar <= 0) << " E: " 
              << Emin <<  " " <<  Emax << " cth:" << cthmin << " " << cthmax << " th:"<< thmin*180.0/TMath::Pi() 
              << " " << thmax*180.0/TMath::Pi() <<  std::endl;

    return Histflv;
}

// Interactive events marginilized over azimuth
std::vector< TH2D* > AsimovSimulation::GetTrueEventsMarginal( ) //To be Deleted
{   
    //Define cos(theta) range and azimuthal angle range
    double cthmin = cos(thmax); // min cos(theta) = -1 
    double cthmax = cos(thmin); // max cos(theta) = 0
    double dphi = phimax - phimin; //Azimutham angle range

    //Bin labels
    // ibins - Bins in Zenith
    // jbins - Bins in Azimuth
    // kbins - Bins in Energy

    //Initialize histograms for neutrino events
    TH2D * EventHist = new TH2D("EventHist","True Event Histrogram", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtseflv = new TH2D("Evtseflv","#nu_{e} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtseflvb = new TH2D("Evtseflvb","#bar{#nu}_{e} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtsmuflv = new TH2D("Evtsmuflv","#nu_{#mu} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax);
    TH2D * Evtsmuflvb = new TH2D("Evtsmuflvb","#bar{#nu}_{#mu} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
                            
     std::vector < TH2D * >  Histflv; //Store Oscillograms for (muflv, muflvbar, eflv, eflvbar, flv+ flvbar)
     
     TH2D * Hist[5];
    
    //Neutrino properties

    // Neutrino initial flavour
    int nue        = 0;  // Electron neutrino  
    int numu       = 1;  // Muon neutrino
    // Particle type
    int nu         = 1;  // Neutrino
    int nubar      = -1; // Antineutrino

    // PMNS matrix setup
    OscProb::PMNS_Fast PMNS_H; 
    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters

    // Load Honda flux data and create histograms for different neutrino types
    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);
    //Generate flux histograms for different neutrino flavours
    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); // Muon neutrino
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData);// Muon antineutrino
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData);  // Electron neutrino
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); // Electron antineutrino

    OscProb::PremModel prem(PremTable);


    // Set Earth model and Detector    
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable); // Earth model
    Earth3D.SetDetector(pos);    // Detector position
    // Modify Specific layer properties
    Earth3D.SetLayerProp(PremTableNumber, 
                         DensityContrast, ChemicalContrast); 
    // Set heterogeneities in lower mantle (e.g., LLVPs)
    Earth3D.PileThickness = PileHeight;  // Height of LLVP
    Earth3D.aWidth = aperture;           // Aperture of LLVP
    Earth3D.SetPile( MantleAnomaly, 
                   AnomalyShape, PileDensityContrast, PileChemContrast);
        
        // Loop over zenith bins (cos(theta))
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    
            // Get the center and width of a bin in cos(theta)
            double cth  = EventHist->GetXaxis()->GetBinCenter(i);
            double dcth = EventHist -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip unphysical values 
            
            // Set the path through the Earth for direction (cos(theta),phi)
            prem.FillPath(cth); // Fill paths from PREM model
            prem.FillPath(cth); // Fill paths from PREM model
            PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
            //Earth3D.SetDirection(cth,0);  
            //std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

            // Distance through "ith" layer : l - EarthPath[i][0];
            // Densty of "ith" layer : rho - EarthPath[i][2];
            // z/a of "ith" layer : zoa - EarthPath[i][3];
            
            // Set initial Earth path in the PMNS model
            //PMNS_H.SetPath(EarthPath[0][0],EarthPath[0][1],EarthPath[0][2]);
            
            // Add subsequent segments of the path
            //for (int i = 1; i < EarthPath.size(); i++) 
            //{ 
              //PMNS_H.AddPath(EarthPath[i][0],EarthPath[i][1],EarthPath[i][2]);
            //} 

            // Loop through energy bins 
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = EventHist->GetYaxis()->GetBinCenter(k); 
                double dE = EventHist->GetYaxis()->GetBinWidth(k);
                double logEi = log10(e);
                
                //Interpolate neutrino fluxes
                 double logdPsiMu = muflux->Interpolate(logEi,cth);
                 double logdPsiMub = mubflux->Interpolate(logEi,cth);
                 double logdPsiE = eflux->Interpolate(logEi,cth);
                 double logdPsiEb = ebflux->Interpolate(logEi,cth);

                //Evaluate differential Flux
                 double dPsiMudEdct = pow(10,logdPsiMu);     //Muon neutrino flux
                 double dPsiMubardEdct = pow(10,logdPsiMub); //Muon anti-neutrino flux
                 double dPsiEdEdct = pow(10,logdPsiE);        //Electron neutrino flux
                 double dPsiEbardEdct = pow(10,logdPsiEb);    //Electron anti-neutrino flux


                //Neutrino Contribution;
                PMNS_H.SetIsNuBar(false); // Set neutrino probabilities 
                // iflv -> fflv oscillation probabilities
                double Pef = PMNS_H.Prob(nue,flvf,e);   // nue -> nufinal
                double Pmuf = PMNS_H.Prob(numu,flvf, e);// numu -> nufinal

                // Neutrino Interacting rate
                double R_ef = XSec(e,nu)*( Pef*dPsiEdEdct);
                double  R_muf = XSec(e,nu)*(Pmuf*dPsiMudEdct); 
                double R_f = R_ef + R_muf; // Events in the final-neutrino channel

                //Antineutrino contribution
                PMNS_H.SetIsNuBar(true); // Set antineutrino probabilities 
                // iflv -> fflv oscillation probabilities
                double Pefb = PMNS_H.Prob(nue,flvf,e);   // nuebar -> nufinalbar
                double Pmufb = PMNS_H.Prob(numu,flvf, e);// numubar -> nufinalbar
                
                // Antineutrino Interacting rate
                double R_efb=XSec(e,nubar)*(  Pefb*dPsiEbardEdct ); 
                double R_mufb=XSec(e,nubar)*( Pmufb *dPsiMubardEdct ); 
                double R_fbar = R_efb + R_mufb; //// Events in the final-antineutrino channel

                // Total Event rate
                double R = R_f+R_fbar;

                //double  mN = 1.67492749804E-27; // Approximate mass of nucleon
                //double  Meff = 14.6*pow(log10(e),1.8 )*MTon;
                double    Meff = 1*MTon;
                //double  N_A = 6.02214E23;
                double  T = 1*years2sec;
                //double  N_f = (2*TMath::Pi())*(Meff*T/mN)*R*dE*dcth;  // Number of events for bin (cth,e)

                  
                //Note : MT/mN  refers to target nucleons per year
                //double  N_f = (2*TMath::Pi())*(MT/mN)*R*dE*dcth;  // Number of events for bin (cth,e)
                double  N_f = (e*e)*(2*TMath::Pi())*(Meff*T/mN)*R*dcth;  // Number of events for bin (cth,e) 

                Evtseflv->SetBinContent(i,k,  ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_ef );
                Evtseflvb->SetBinContent(i,k, ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_efb );
                Evtsmuflv->SetBinContent(i,k, ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_muf );
                Evtsmuflvb->SetBinContent(i,k,((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_mufb );
                EventHist->SetBinContent(i,k, N_f ); //Create histogram for  kth Pseudo-Experimens

            } // End of loop in energy

            // Clear PMNS path for the next iteration
            PMNS_H.ClearPath();

        }   //   End of loop in zenith

    //Store all histrogram in a vector
    
    Histflv.push_back(Evtsmuflv); // Muon neutrino contribution to events
    Histflv.push_back(Evtsmuflvb);// Muon antineutrino contribution to events
    Histflv.push_back(Evtseflv);  // Electron neutrino contribution to events
    Histflv.push_back(Evtseflvb); // Electron antineutrino contribution to events
    Histflv.push_back(EventHist); // flvf neutrino events

    // Simulation summary output
    std::cout << "Detector location" << pos[2] << std::endl;
    std::cout << "PREM model: " << PremTable  <<std::endl;
    std::cout << "Flux Table: " << HondaTable  <<std::endl;
    std::cout << "Pile| h:" <<  PileHeight << " aWidth:" << aperture 
              << " Shape: " << AnomalyShape << " rhopct: " << PileDensityContrast << std::endl;
    std::cout <<" Simulation Summary| flvf:" << flvf << " IsNuBar:" <<  (nunubar <= 0) << " E: " 
              << Emin <<  " " <<  Emax << " cth:" << cthmin << " " << cthmax << " th:"<< thmin*180.0/TMath::Pi() 
              << " " << thmax*180.0/TMath::Pi() <<  std::endl;

    return Histflv;
}


TGraph * AsimovSimulation::GetOscProb( int flvi, int flvf, bool nunubar, double cth ) //To be Deleted
{  
    
    //Energy Intervals
    //double Emin      = EnuMin;//Lower limit for Energy
    //double Emax      = EnuMax;//Upper limit for Energy
    int n = 1000;   //Bins in Energy
    
//Neutrino Oscillation Probabilities calculation--------------------------------------------------------------------
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects

    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters


///----------------

    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);

//Set earth model -------------------------------------------------------------------------------------------------------------
    
        

     std::cout << "3D" << std::endl;
    /*
     Earth3DModel Earth3D;
     Earth3D.SetModel(PremTable);
     Earth3D.SetDetector(pos);
     Earth3D.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
     Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
    double l,d,z,ly;
    double phi = 0.0; //< This will defined a constant L for different values of ct provided Dct is Small    
    Earth3D.SetDirection(cth,phi); 
    std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

    l = EarthPath[0][0];
    d = EarthPath[0][1];
    z = EarthPath[0][2];
    //ly =  PathMatrix[0][3]; 
        
    PMNS_H.SetPath(l,d,z);
    for (int i = 1; i < EarthPath.size(); i++) 
    { 
    
        l = EarthPath[i][0];
        d = EarthPath[i][1];
        z = EarthPath[i][2];
        PMNS_H.AddPath(l,d,z);
        
    } 
    
    */
    
    std::cout << "OscProb" << std::endl;
    OscProb::PremModel prem(PremTable);
    prem.FillPath(cth); // Fill paths from PREM model
    prem.FillPath(cth); // Fill paths from PREM model
    PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
    
    PMNS_H.SetIsNuBar(nunubar);
    
    //Oscillation Probabilities

    double dE = (Emax-Emin)/1000.0;

    //TH1D * EHist = new TH1D("EHist", "Energy Value", 1000 , Emin, Emax);

    std::cout <<  Emin << " " << Emax <<  std::endl; 

    TGraph * g = new TGraph(1001);

    for (int i = 0; i < 1001; ++i)
        {
            //double e=EHist->GetXaxis()->GetBinCenter(i);
            double E = Emin + i*dE;
            double pab = PMNS_H.Prob(flvi,flvf,E);
            g->SetPoint(i,E,pab);

          //  std::cout << E << " " << pab << std::endl;
            
        }    

   
    return g;
}


