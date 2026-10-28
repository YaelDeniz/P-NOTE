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
#include "InteractingEvents.h"

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
#include "Earth3DModel.h"



using namespace std;

//Simulation with binning in cos(theta)
std::vector< TH2D* > NintSimulation::GetNint_ijk()
{  
    /*
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string AziIndex = NuTomoPath+IntFolder+"TrueIndexTable.txt";
    std::ofstream IndexPhi(AziIndex); 

    std::cout << "Not here 1" << std::endl;
    */

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    //Binning Scheme 
    double cthmin = cos(thmax); // Min value for cth (cth = -1) 
    double cthmax = cos(thmin); // Max value for cth (cth = 1)

    //ibins - Total number of Energy bins
    //jbins - Total number of zenith("theta") bins | cos(theta)
    //kbins - Total number of azimuth bins
    double N_ijk = 0   ;        // Poisson mean for bin (ii,jj,kk).

    // Histogram definition (3D: cos(theta), azimuth, energy)
    // SetDirectory(0) detaches from gDirectory so repeated calls don't collide on name
    static int callCount = 0;
    ++callCount;
    TH3D * EventHist3D = new TH3D(Form("TrueHist_%d", callCount),"True Event Histrogram",
                                  jbins,cthmin,cthmax,
                                  kbins,phimin,phimax,
                                  ibins,Emin,Emax);
    EventHist3D->SetDirectory(0);


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
    Earth3D.SetPile( MantleAnomaly, AnomalyShape, Pilerho_pct, PileChemContrast);
    //Earth3D.SetLayerProp(PremLayer, rho_pct, zoa_pct);

    std::cout << "Not here 3" << PremTable <<  std::endl;


    // Initialize vector for 2D histograms and other variables
    std::vector< TH2D* > HistVec;

    TH2D * EventHist2D[kbins]; // kbins - Total number of azimuth bins

    double l,d,z; // path variables length, density, depth, resppectively

    // Cache file: azimuth bin index table (written only when OutputDir is set)
    std::ofstream IndexPhi;
    bool writeCacheFile = !OutputDir.empty();
    if (writeCacheFile) {
        IndexPhi.open(OutputDir / "Nint_PhiIndexTable.txt");
    }

    //Loop over azimuth angles
    for (int kk = 1; kk <= kbins; kk++)
    {

        std::cout << "Not here 4 " << kk <<  std::endl;


        // Get azimuth center and width for the current bin
        double phi = EventHist3D->GetYaxis()->GetBinCenter(kk);
        double dphi = EventHist3D-> GetYaxis()->GetBinWidth(kk);

        if (writeCacheFile && IndexPhi.is_open()) {
            IndexPhi << kk << " " << phi*180/TMath::Pi() << " " << dphi*180/TMath::Pi() << std::endl;
        }

        std::cout << "Not here 5 IndexTable is not a problem " << std::endl;

        TString histname = Form("OscHist_%d_%d", callCount, kk);

        std::cout << "  histname is note the problem " << std::endl;

        EventHist2D[kk] = new TH2D(histname,Form("Oscillogram%d",kk), jbins,cthmin,cthmax,ibins,Emin,Emax);
        EventHist2D[kk]->SetDirectory(0);

        std::cout << "  Vector is the problem  " << std::endl;


        // Loop over zenith bins (cos(theta))
        for(int jj=1; jj<= jbins ; jj++)
        {    
            // Get cos(theta) center and width for the current bin.
            // The distance traveled through the earth "Baseline" is L=-2 * R_earth * cos(theta)

            double cth = EventHist2D[kk]->GetXaxis()->GetBinCenter(jj);
            double dcth = EventHist2D[kk] -> GetXaxis()->GetBinWidth(jj);

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
            for (int ii = 1; ii <=ibins; ++ii)
            { 
                double e = EventHist2D[kk]->GetYaxis()->GetBinCenter(ii); 
                double dE = EventHist2D[kk]->GetYaxis()->GetBinWidth(ii);
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
                // Note : /mN  refers to target nucleons per year
                double N_ijk = (Exposure_det)*(Ri_nu + Ri_nubar)*dE*dcth*dphi; 

                // Set the bin content in the 2D histogram
                EventHist2D[kk]->SetBinContent(jj,ii, N_ijk); 

            
            } // End energy loop 

        } // End zenith loop

        // Add the current 2D histogram to the vector

        std::cout << "All good for phi " << kk <<  std::endl;

        HistVec.push_back(EventHist2D[kk]);

    } //End  Azimuth Loop 


    // Close index file and return the vector of histograms
    if (IndexPhi.is_open()) IndexPhi.close();
    return HistVec;
}

