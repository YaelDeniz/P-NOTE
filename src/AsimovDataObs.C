#include "GetObservedEvents.h"

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
#include <array>

//Cern Root
#include <math.h>
#include "TF1.h"
#include "TMath.h"
#include "TH2.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Interpolator.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/RootFinder.h"
//#include "TGraphTrueD.h"

//OSCPROB

//#include "PremModel.h"
//#include "PMNS_Fast.h"


// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"
#include "Earth3DModel.h"

//Layer Updates
#include "LayerUpdate.h"


using namespace std;


// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

# define R_earth 6368.0 //km
void AsimovObsSimulation::ChangeLayerProperties(int n,
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
    for (auto& upd : LayerUpdates_recosim) {
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

    LayerUpdates_recosim.push_back(NewUpdate);
}


void AsimovObsSimulation::AddUpdatesList(const std::vector<LayerUpdate>& updates)
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





//RESOLUTION FUNCITONS----------------------------------------------------------

//Energy Resolution Parameters
void AsimovObsSimulation::SetEnergyResolution(double A, double B){
    Ae = A;
    Be = B;}

//Energy Resolution function
double AsimovObsSimulation::PDFE(double Ereco , double Etrue ) {
    double SigmaE = Ae*Etrue +Be*sqrt(Etrue); //Detector Energy Resolution
    double pdfe = ROOT::Math::gaussian_pdf( Ereco, SigmaE, Etrue) ;
    return pdfe; } //Gaussian Distribution

//Angular Resolution Parameters (rad)
void AsimovObsSimulation::SetAngularResolution(double A1, double B1, double A2, double B2){
    Ath = A1;
    Bth = B1;
    Aphi = A2;
    Bphi = B2;}

//Angular(Zenith only) Resolution function
double AsimovObsSimulation::PDFNth(double threco, double Etrue, double thtrue){
    double Sigmath = Ath + Bth/(sqrt(Etrue)); //Detector Angular Resolution
    double pdfth = ROOT::Math::gaussian_pdf( threco, Sigmath, thtrue) ;
    return pdfth; } //Gaussian Distribution

//Angular(azimuth only) Resolution function
double AsimovObsSimulation::PDFNphi(double phireco, double Etrue, double phitrue){
    double Sigmaphi= Aphi + Bphi/(sqrt(Etrue)); //Detector Angular Resoluton
    double pdfphi = ROOT::Math::gaussian_pdf( phireco, Sigmaphi, phitrue) ;
    return pdfphi; } //Gaussian Distribution

double AsimovObsSimulation::PDFth(double threco, double Etrue, double thtrue){

    /// @param sigma - std. dev. in radians. Angular resolution
    /// kappa = 1/sigma^2 - VMF concentrarion parameter. 

    double sigmath = Ath + Bth/(sqrt(Etrue)); // Angular resolution in radians

    //double cos_z = cthreco;
    //double cos_z0 = cthtrue;

    double cos_z  = cos(threco);
    double cos_z0 = cos(thtrue);

    double sin_z = sqrt(1 - pow(cos_z, 2));
    double sin_z0 = sqrt(1 - pow(cos_z0, 2));

    double kappa = 1/(sigmath*sigmath);

    double exp_arg = kappa * cos_z * cos_z0;

    double norm = kappa;
    
    if (kappa < 50)
        norm /= 2 * sinh(kappa);
    else
        exp_arg -= kappa;

    double arg = kappa * sin_z * sin_z0;
    double bess = 1;
    if (arg < 50)
        bess = TMath::BesselI0(arg);
    else 
    {
        exp_arg += arg;
        bess *=
            (1 + 1. / (8 * arg) + 4.5 / pow(8 * arg, 2) + 37.5 / pow(8 * arg, 3));
        bess /= sqrt(2 * M_PI * arg);
    }

    double expo = exp(exp_arg);

    return norm * expo * bess * sin_z;} //Von Mises Fisher for Zenith

double AsimovObsSimulation::PDFangle(std::array<double, 2> recodir, std::array<double, 2> truedir, double Etrue){

/// Computes the angular probability using a von Mises–Fisher-like distribution
    /// @param recodir: pointer to reconstructed direction [theta_rec, phi_rec]
    /// @param truedir: pointer to true direction [theta_true, phi_true]
    /// @param Etrue: true neutrino energy (for resolution dependence)

    // Angular resolution sigma in radians (same for both theta and phi)

    double sigma = Ath + Bth/(sqrt(Etrue)); 

    //double cos_z = cthreco;
    //double cos_z0 = cthtrue;

    double cth  = cos(recodir[0]);
    double sth  = sqrt(1-cth*cth);
    double cth0 = cos(truedir[0]);
    double sth0 = sqrt(1-cth0*cth0);

    double phi = recodir[1];
    double phi0 = truedir[1];

    double sphi =  sin(recodir[1]);
    double sphi0 = sin(truedir[1]);
    
    // von Mises–Fisher concentration parameter
    double kappa = 1/(sigma*sigma);

    // Argument of the exponential in the vMF kernel on the sphere
    double exp_arg = kappa * cth * cth0;
    exp_arg += kappa * sth * sth0*cos(phi - phi0);

    // Normalization factor
    double norm = kappa;
    
    if (kappa < 50)
        norm /= 4* TMath::Pi()* sinh(kappa);
    else
    {
        exp_arg -= kappa;// Compensate for exponential overflow
    
        // For large kappa (kappa > 50), sinh(kappa) ≈ 0.5 * e^kappa
        norm /= (2.0 * TMath::Pi());
    }

    double expo = exp(exp_arg);

    return norm * expo * sth;} //Von Mises Fisher for

// Observed events in cth binning
std::vector < TH2D* >  AsimovObsSimulation::GetObsEvents3Dcth(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin    = EnuMin;
    double ErecoMax    = EnuMax;
    
    //Observed Angle
    double threcoMin   = ZenMin;
    double threcoMax   = ZenMax;

    double CthrecoMin  = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    double CthrecoMax  = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,ErecoMin);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoMax);
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double EtrueMin = finder_min.Root();
    double EtrueMax = finder_max.Root();

    //True Angle-
    double thtrueMin_a= (threcoMin*TMath::Pi()/180.0)-4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= (threcoMax*TMath::Pi()/180.0)+4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;


    //std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    //std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;


    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location



    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;




    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf=flvf;


    std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3D(); //cos(theta) binning



    std::vector<TH2D*> HistVec;
    TH2D* RecoHist[jbins]; 

    //Observed distribution:
    //TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    //Loop Over Histogram in Azimuth
    for (int nphi = 0 ; nphi < TrueHist.size(); ++nphi)
    {

        RecoHist[nphi] = new TH2D(Form("RecoHist%d",nphi),Form("ObsOscillogram%d",nphi),
                                  mbins,CthrecoMin,CthrecoMax,
                                  nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

        //Loop over Reconstructed cth
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

            double cthreco = RecoHist[nphi]->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
            double threco  =  acos(cthreco);
            double sthreco =  sin(threco); 
            double dcthreco= RecoHist[nphi]->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram

            //Loop over Reconstructed energy
            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                double Ereco = RecoHist[nphi]->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
                double dEreco = RecoHist[nphi]->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
           
                double Nsmear = 0;   // Reset total number of events.
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    // Get cos(theta_z) from bin center, It fix a value of L
                    double cthtrue = TrueHist[nphi]->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                    double thtrue = acos(cthtrue);
                    double dcthtrue = TrueHist[nphi]->GetXaxis()->GetBinWidth(i); 

                    for (int k = 1; k <= kbins; ++k)
                    { 
                        //NUMERICAL INTEGRATION


                        double Etrue = TrueHist[nphi]->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                        double dEtrue = TrueHist[nphi]->GetYaxis()->GetBinWidth(k); 

                        double NikTrue = TrueHist[nphi]->GetBinContent(i,k);

                        double dEdcthReco = dEreco*dcthreco;


                         ZenithBinSize = TrueHist[nphi]->GetXaxis()->GetXmax();
                         AzimuthBinSize =0;
                         EnergyBinSize = TrueHist[nphi]->GetYaxis()->GetXmax();
                       
                        Nsmear += (PDFth(threco,Etrue,thtrue)/sthreco)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdcthReco);
        
                    } // loop e

                } // Loop ct
                
                double Nmn = Nsmear;

                RecoHist[nphi]->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 
             
            } // End reconstructed E loop 

        } // End reconstructed cth loop

        HistVec.push_back(RecoHist[nphi]); 

    } // End Histogram loop
   
   return HistVec;}

// Observed events in th binning
std::vector < TH2D* >  AsimovObsSimulation::GetObsEvents3Dth(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle (rad)
    double threcoMin   = ZenMin*TMath::Pi()/180.0;
    double threcoMax   = ZenMax*TMath::Pi()/180.0;

    //double CthrecoMin = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    //double CthrecoMax = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,ErecoMin);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoMax);
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double EtrueMin = finder_min.Root();
    double EtrueMax = finder_max.Root();

    //True Angle-
    double thtrueMin_a= threcoMin - 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= threcoMax + 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;
    

    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location

    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;


    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf=flvf; 

    //Binning in zenith angle
    std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3Dth();


    std::vector<TH2D*> HistVec;
    TH2D* RecoHist[jbins]; 

    //Observed distribution:
    //TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    //Loop Over Histogram in Azimuth
    for (int nphi = 0 ; nphi < TrueHist.size(); ++nphi)
    {

        RecoHist[nphi] = new TH2D(Form("RecoHist%d",nphi),Form("ObsOscillogram%d",nphi),
                                  mbins,threcoMin,threcoMax,
                                  nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

        //Loop over Reconstructed cth
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

            double threco = RecoHist[nphi]->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
            double dthreco= RecoHist[nphi]->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram

            double cthreco  =  cos(threco);

            //Loop over Reconstructed energy
            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                double Ereco = RecoHist[nphi]->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
                double dEreco = RecoHist[nphi]->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
           
                double Nsmear = 0;   // Reset total number of events.
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    // Get cos(theta_z) from bin center, It fix a value of L
                    double thtrue = TrueHist[nphi]->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dthtrue = TrueHist[nphi]->GetXaxis()->GetBinWidth(i); 

                    double cthtrue = cos(thtrue);

                    for (int k = 1; k <= kbins; ++k)
                    { 
                        //NUMERICAL INTEGRATION

                        double Etrue = TrueHist[nphi]->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                        double dEtrue = TrueHist[nphi]->GetYaxis()->GetBinWidth(k); 

                        double NikTrue = TrueHist[nphi]->GetBinContent(i,k);

                        double dEdthReco = dEreco*dthreco;
                       
                        Nsmear += PDFth(threco,Etrue,thtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdthReco);

                    } // loop e

                } // Loop th
                
                double Nmn = Nsmear;

                RecoHist[nphi]->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 
             
            } // End reconstructed E loop 

        } // End reconstructed th loop

        HistVec.push_back(RecoHist[nphi]); 

    } // End Histogram loop

    std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   return HistVec;}

// Observed events in th binning
std::vector < TH2D* >  AsimovObsSimulation::GetObsEvents3Dthlog10E(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle (rad)
    double threcoMin   = ZenMin*TMath::Pi()/180.0;
    double threcoMax   = ZenMax*TMath::Pi()/180.0;

    //double CthrecoMin = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    //double CthrecoMax = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,ErecoMin);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoMax);
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double EtrueMin = finder_min.Root();
    double EtrueMax = finder_max.Root();

    //True Angle-
    double thtrueMin_a= threcoMin - 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= threcoMax + 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;
    

    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location

    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;


    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf=flvf; 

    //Binning in zenith angle
    std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3Dthlog10E();


    std::vector<TH2D*> HistVec;
    TH2D* RecoHist[jbins]; 

    //Observed distribution:
    //TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    //Loop Over Histogram in Azimuth
    for (int nphi = 0 ; nphi < TrueHist.size(); ++nphi)
    {

        RecoHist[nphi] = new TH2D(Form("RecoHist%d",nphi),Form("ObsOscillogram%d",nphi),
                                  mbins,threcoMin,threcoMax,
                                  nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

        //Loop over Reconstructed cth
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

            double threco = RecoHist[nphi]->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
            double dthreco= RecoHist[nphi]->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram

            double cthreco  =  cos(threco);

            //Loop over Reconstructed energy
            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                double log10Ereco = RecoHist[nphi]->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
                double Ereco = pow(10, log10Ereco);
                double dlog10Ereco = RecoHist[nphi]->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
           
                double Nsmear = 0;   // Reset total number of events.
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    // Get cos(theta_z) from bin center, It fix a value of L
                    double thtrue = TrueHist[nphi]->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dthtrue = TrueHist[nphi]->GetXaxis()->GetBinWidth(i); 

                    double cthtrue = cos(thtrue);

                    for (int k = 1; k <= kbins; ++k)
                    { 
                        //NUMERICAL INTEGRATION

                        double log10Etrue = TrueHist[nphi]->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                        double Etrue = pow(10, log10Etrue);
                        double dlog10Etrue = TrueHist[nphi]->GetYaxis()->GetBinWidth(k); 

                        double NikTrue = TrueHist[nphi]->GetBinContent(i,k);

                        double dlog10EdthReco = dlog10Ereco*dthreco;
                       
                        Nsmear += PDFth(threco,Etrue,thtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(log(10)*Ereco)*(dlog10EdthReco);

                    } // loop e

                } // Loop th
                
                double Nmn = Nsmear;

                RecoHist[nphi]->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 
             
            } // End reconstructed E loop 

        } // End reconstructed th loop

        HistVec.push_back(RecoHist[nphi]); 

    } // End Histogram loop

    std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   return HistVec;}

// Observed events in th binning - Variable binning
std::vector < TH2D* >  AsimovObsSimulation::GetObsEvents3Dthvb(){  

   
    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of Zenith bins of Observed event distribution
    //int nbins = recobinsAzi; //Number of Azimuth bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle (rad)
    double threcoMin   = ZenMin*TMath::Pi()/180.0;
    double threcoMax   = ZenMax*TMath::Pi()/180.0;

    //double CthrecoMin = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    //double CthrecoMax = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0

    //Create bins of variable size

    std::vector<double> ErecoEdges; 
    ErecoEdges.push_back(ErecoMin);
    double Erecoi = ErecoMin;

    while (ErecoEdges.back() < ErecoMax) 
    {
        //Find center
        TF1 *frecoC = new TF1("fmin","-1*x + [0] + 0.5*( [1]*x+[2]*sqrt(x) )", 0, 100);
        frecoC->SetParameter(0,Erecoi);
        frecoC->SetParameter(1,Ae);
        frecoC->SetParameter(2,Be);
        ROOT::Math::RootFinder finder_center;
        finder_center.Solve(*frecoC ,0,100);
        double ErecoC = finder_center.Root();
        
        //Find next edge
        double Erecoii = ErecoC + 0.5*(Ae*ErecoC + Be*sqrt(ErecoC));
        
        ErecoEdges.push_back(Erecoii);
        
        Erecoi = Erecoii;
        
    }

    double ErecobinEdges [ErecoEdges.size()];
        
    for (int i=0 ; i < ErecoEdges.size(); i++ )
    {
          ErecobinEdges[i] = ErecoEdges[i];
            
    }
    
    int nbins = sizeof(ErecobinEdges)/sizeof(ErecobinEdges[0]) - 1; // Number of bins

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,ErecoMin);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoEdges.back());
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double EtrueMin = finder_min.Root();
    double EtrueMax = finder_max.Root();

    //True Angle-
    double thtrueMin_a= threcoMin - 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= threcoMax + 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;
    

    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location

    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;


    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf=flvf; 

    //Binning in zenith angle
    std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3Dth();


    std::vector<TH2D*> HistVec;
    TH2D* RecoHist[jbins]; 

    //Observed distribution:
    //TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    //Loop Over Histogram in Azimuth
    for (int nphi = 0 ; nphi < TrueHist.size(); ++nphi)
    {

        RecoHist[nphi] = new TH2D(Form("RecoHist%d",nphi),Form("ObsOscillogram%d",nphi),
                                  mbins,threcoMin,threcoMax,
                                  nbins,ErecobinEdges); // xbins correspond to energy values and ybins to zenith angle cost

        //Loop over Reconstructed cth
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

            double threco = RecoHist[nphi]->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
            double dthreco= RecoHist[nphi]->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram

            double cthreco  =  cos(threco);

            //Loop over Reconstructed energy
            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                double Ereco = RecoHist[nphi]->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
                double dEreco = RecoHist[nphi]->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
           
                double Nsmear = 0;   // Reset total number of events.
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    // Get cos(theta_z) from bin center, It fix a value of L
                    double thtrue = TrueHist[nphi]->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dthtrue = TrueHist[nphi]->GetXaxis()->GetBinWidth(i); 

                    double cthtrue = cos(thtrue);

                    for (int j = 1; j <= kbins; ++j)
                    { 
                        //NUMERICAL INTEGRATION

                        double Etrue = TrueHist[nphi]->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small
                        double dEtrue = TrueHist[nphi]->GetYaxis()->GetBinWidth(j); 

                        double NikTrue = TrueHist[nphi]->GetBinContent(i,j); //Number of neutrin

                        double dEdthReco = dEreco*dthreco;
                       
                        Nsmear += PDFth(threco,Etrue,thtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdthReco);

                    } // loop e

                } // Loop th
                
                double Nmn = Nsmear;

                RecoHist[nphi]->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 
             
            } // End reconstructed E loop 

        } // End reconstructed th loop

        HistVec.push_back(RecoHist[nphi]); 

    } // End Histogram loop

    std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   return HistVec;}

// Observed events in th binning - Variable binning

TH3D*  AsimovObsSimulation::GetRecoN(){  

   
   // ------------------------------------------------------------------------------------------------------
   // 1. OBSERVED VARIABLES & BINNING SETUP - Binning scheme: Reco Energy(l), Reco Zenith(m), Reco Azimuth(n) 
   // ------------------------------------------------------------------------------------------------------
    
    //int lbins = recobinsE
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsAzi; //Number of Energy bins of Observed event distribution 
    
    //Reconstructed Energy range [GeV]
    double E_reco_min      = EnuMin;
    double E_reco_max      = EnuMax;
    
    //Reconstructed Zenith range (rad)
    double theta_reco_min = ZenMin*TMath::Pi()/180.0;
    double theta_reco_max = ZenMax*TMath::Pi()/180.0;

    //Reconstructed Azimuth range (rad)
    double phi_reco_min = AziMin*TMath::Pi()/180.0 ;
    double phi_reco_max = AziMax*TMath::Pi()/180.0 ;


    //Create bins of variable size

    std::vector<double> ErecoEdges; 
    ErecoEdges.push_back(E_reco_min);
    double Erecoi = E_reco_min;

    while (ErecoEdges.back() < E_reco_max) 
    {
        //Find center
        TF1 *frecoC = new TF1("fmin","-1*x + [0] + 0.5*( [1]*x+[2]*sqrt(x) )", 0, 100);
        frecoC->SetParameter(0,Erecoi);
        frecoC->SetParameter(1,Ae);
        frecoC->SetParameter(2,Be);
        ROOT::Math::RootFinder finder_center;
        finder_center.Solve(*frecoC ,0,100);
        double ErecoC = finder_center.Root();
        
        //Find next edge
        double Erecoii = ErecoC + 0.5*(Ae*ErecoC + Be*sqrt(ErecoC));
        
        ErecoEdges.push_back(Erecoii);
        
        Erecoi = Erecoii;
        
    }

    double ErecobinEdges [ErecoEdges.size()];
        
    for (int i=0 ; i < ErecoEdges.size(); i++ )
    {
          ErecobinEdges[i] = ErecoEdges[i];
            
    }
    
    int lbins = sizeof(ErecobinEdges)/sizeof(ErecobinEdges[0]) - 1; // Number of bins

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsE; 
    int jbins = truebinsZen;
    int kbins = truebinsAzi; 

    
    //True Energy range [GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,E_reco_min);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoEdges.back()); // "Back" indicate upper bound in Energy obtained in variable binning
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double E_true_min = finder_min.Root();
    double E_true_max = finder_max.Root();

    std::cout << " ****************************************** " << E_true_min << std::endl;
   
    // True Zenith Range [deg]
    double theta_true_min = max(0.0,theta_reco_min - 4*(Ath + Bth/sqrt(E_true_min))) *180.0/TMath::Pi(); 
    double theta_true_max = min(TMath::Pi(),theta_reco_max + 4*(Ath + Bth/sqrt(E_true_min))) *180.0/TMath::Pi(); 


    //True Azimuth range [deg]
    double phi_true_min = ( (phi_reco_min) - 4*(Aphi + Bphi/sqrt(E_true_min)) )*180.0/TMath::Pi(); 

     //std::cout << "phi  ****************************************** " <<phi_reco_min*180.0/TMath::Pi() << " " << phi_true_min << std::endl;

    //if(phi_true_min < 0.0 ) phi_true_min = phi_true_min + 360.0;

     //std::cout << "phi fix ****************************************** " <<phi_reco_min*180.0/TMath::Pi() << " " << phi_true_min << std::endl;



    double phi_true_max = ( (phi_reco_max) + 4*(Aphi + Bphi/sqrt(E_true_min)) )*180.0/TMath::Pi(); 

     //std::cout << "phi  ****************************************** " <<phi_reco_max*180.0/TMath::Pi() << " " << phi_true_max << std::endl;

    //if(phi_true_max > 360.0 ) phi_true_max = phi_true_max - 360.0;

    std::cout << "Theta  ****************************************** " <<theta_reco_min*180.0/TMath::Pi() << " " << theta_true_min << std::endl;
    std::cout << "Theta  ****************************************** " <<theta_reco_max*180.0/TMath::Pi() << " " << theta_true_max << std::endl;

    std::cout << "phi fix min ****************************************** " <<phi_reco_min*180.0/TMath::Pi() << " " << phi_true_min << std::endl;
    std::cout << "phi fix max****************************************** " <<phi_reco_max*180.0/TMath::Pi() << " " << phi_true_max << std::endl;

    //delete frecoC;
    //delete fmin;

    
    // --------------------------------------------------------------------------------------
    // 3. HISTOGRAM SETUP & RECONSTRUCTION SIMULATION
    // --------------------------------------------------------------------------------------


    //Observed distribution:
    TH3D* RecoHist = new TH3D("RecoHist","Reconstructed Signal",lbins,E_reco_min,E_reco_max,
                                                                mbins,theta_reco_min,theta_reco_max,
                                                                nbins,phi_reco_min,phi_reco_max); // xbins correspond to energy values and ybins to zenith angle cost


    double x_min = RecoHist->GetXaxis()->GetXmin();
    double x_max = RecoHist->GetXaxis()->GetXmax();
    double x_bins =RecoHist->GetXaxis()->GetNbins();

    double y_min = RecoHist->GetYaxis()->GetXmin()*180.0/TMath::Pi();
    double y_max = RecoHist->GetYaxis()->GetXmax()*180.0/TMath::Pi();
    double y_bins =RecoHist->GetYaxis()->GetNbins();

    double z_min = RecoHist->GetZaxis()->GetXmin()*180.0/TMath::Pi();
    double z_max = RecoHist->GetZaxis()->GetXmax()*180.0/TMath::Pi();
    double z_bins =RecoHist->GetZaxis()->GetNbins();

 

    // --------------------------------------------------------------------------------------
    // 3. HISTOGRAM SETUP & RECONSTRUCTION SIMULATION
    // --------------------------------------------------------------------------------------

    // Set up true simulation properties
    AsimovSimulation TrueSimulation;
    TrueSimulation.PremTable = ThePremTable;
    TrueSimulation.HondaTable = TheHondaTable;
    TrueSimulation.SetDetectorPosition(xyzTelescope);

    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape = ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight;
    TrueSimulation.aperture = ThePileAperture;
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;

    //TrueSimulation.ModifyLayer(PremLayer, DensityContrast, ChemicalContrats);
    if (!LayerUpdates_recosim.empty()) {
    TrueSimulation.AddUpdatesList(LayerUpdates_recosim);
    }

    TrueSimulation.SetIntervals(theta_true_min, theta_true_max, phi_true_min, phi_true_max, E_true_min, E_true_max);
    TrueSimulation.SetBinning(ibins,jbins,kbins);

    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf = flvf;

    // Retrieve the true histogram
    TH3D* TrueHist = TrueSimulation.GetIntRates();

    double xo_min = TrueHist->GetXaxis()->GetXmin();
    double xo_max = TrueHist->GetXaxis()->GetXmax();
    double xo_bins =TrueHist->GetXaxis()->GetNbins();

    double yo_min = TrueHist->GetYaxis()->GetXmin() * 180.0/TMath::Pi();
    double yo_max = TrueHist->GetYaxis()->GetXmax() * 180.0/TMath::Pi();
    double yo_bins =TrueHist->GetYaxis()->GetNbins() ;

    double zo_min = TrueHist->GetZaxis()->GetXmin() *180.0/TMath::Pi();
    double zo_max = TrueHist->GetZaxis()->GetXmax() *180.0/TMath::Pi();
    double zo_bins =TrueHist->GetZaxis()->GetNbins();

    std::cout << "TH3D RECO Axis Bounds:" << std::endl;
    std::cout << "  X axis: [" << x_min << ", " << x_max << ": " << x_bins << "]" << std::endl;
    std::cout << "  Y axis: [" << y_min << ", " << y_max << ": " << y_bins << "]" << std::endl;
    std::cout << "  Z axis: [" << z_min << ", " << z_max << ": " << z_bins << "]" << std::endl;

    std::cout << "TH3D INT Axis Bounds:" << std::endl;
    std::cout << "  X axis: [" << xo_min << ", " << xo_max << ": " << xo_bins << "]" << std::endl;
    std::cout << "  Y axis: [" << yo_min << ", " << yo_max << ": " << yo_bins << "]" << std::endl;
    std::cout << "  Z axis: [" << zo_min << ", " << zo_max << ": " << zo_bins << "]" << std::endl;

    double Ireco = 0; //Integrals


    //Binning in zenith angle
    //std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3Dth();
    //std::vector<TH2D*> HistVec;
    //TH2D* RecoHist[jbins]; ErecoMin,ErecoMax


    // --------------------------------------------------------------------------------------
    // 4. EVENT RECONSTRUCTION SIMULATION
    // --------------------------------------------------------------------------------------

    //double E_reco,theta_reco,phi_reco;
    //double dE_reco, dTheta_reco, dPhi_reco, dV_reco;
    //double Rlmn , Nlmn;
    
    //double E_true,theta_true,phi_true;
    //double dE_true, dth_true, dphi_true, dV_true;
    //double Rijk;



    for (int n = 1; n <= nbins; ++n) {
        double phi_reco = RecoHist->GetZaxis()->GetBinCenter(n);
        double dPhi_reco = RecoHist->GetZaxis()->GetBinWidth(n);

        for (int m = 1; m <= mbins; ++m) {
            double theta_reco = RecoHist->GetYaxis()->GetBinCenter(m);
            double dTheta_reco = RecoHist->GetYaxis()->GetBinWidth(m);

            for (int l = 1; l <= lbins; ++l) {
                double E_reco = RecoHist->GetXaxis()->GetBinCenter(l);
                double dE_reco = RecoHist->GetXaxis()->GetBinWidth(l);

                double Rlmn = 0;
                double Nlmn = 0;

                for (int k = 1; k <= kbins; ++k) {
                    
                    double phi_true = TrueHist->GetZaxis()->GetBinCenter(k);

                    if( phi_true < 0.0 ) phi_true = phi_true + 2*TMath::Pi();
                    if( phi_true > 2*TMath::Pi()) phi_true = phi_true - 2*TMath::Pi();


                    //std::cout << " phi: " << TrueHist->GetZaxis()->GetBinCenter(k)*180.0/TMath::Pi() << " phi wrap " << phi_true*180.0/TMath::Pi() << std::endl;

                    for (int j = 1; j <= jbins; ++j) {
                        
                        double theta_true = TrueHist->GetYaxis()->GetBinCenter(j);

                        for (int i = 1; i <= ibins; ++i) {

                            std::array<double, 2> recodir = {theta_reco,phi_reco};
                            std::array<double, 2> truedir = {theta_true,phi_true}; // [theta, phi] in radians

                            
                            double E_true = TrueHist->GetXaxis()->GetBinCenter(i);

                            //True Bin widths
                            double dE_true = TrueHist->GetXaxis()->GetBinWidth(i);
                            double dth_true = TrueHist->GetYaxis()->GetBinWidth(j);
                            double dphi_true = TrueHist->GetZaxis()->GetBinWidth(k);
  
                            //True rates
                            double Rijk = TrueHist->GetBinContent(i, j, k);


                            //double PDF_dir = PDFNphi(phi_reco,E_true,phi_true)* PDFth(theta_reco, E_true, theta_true);
                            double dV_true = dE_true * dth_true * dphi_true;
                            

                            //std::cout << theta_true_min << " " << theta_true_max <<" "<< phi_true_min<<" "<<  phi_true_max <<" "<< E_true_min<< " " <<E_true_max<<std::endl;
                            //std::cout << theta_reco_min << " " << theta_reco_max <<" "<< phi_reco_min<<" "<<  phi_reco_max <<" "<< E_reco_min<< " " <<E_reco_max<<std::endl;
                            //std::cout << theta_true << " " << theta_reco << " " << phi_true << " " << phi_reco << " " << E_true << " " << E_reco <<std::endl;

                            //std::cout << "True Indeixng [" << Rijk << "] " << i << ":" << ibins << " " << j << ":" << jbins<< " " << k << ":" << kbins <<  std::endl;
                            //std::cout << TrueHist->GetXaxis()->GetXmin() << " " << TrueHist->GetXaxis()->GetXmax() << " " << TrueHist->GetYaxis()->GetXmin() << " " << TrueHist->GetYaxis()->GetXmax() << " " << TrueHist->GetZaxis()->GetXmin() << " " << TrueHist->GetZaxis()->GetXmax() << std::endl;
                            //std::cout << dE_true << " " << dth_true << " " << dphi_true << std::endl;

                            //std::cout << "Observed Indeixng" << l << ":" << lbins << " " << m << ":" << mbins<< " " << n << ":" << nbins <<  std::endl;
                            //std::cout << RecoHist->GetXaxis()->GetXmin() << " " << RecoHist->GetXaxis()->GetXmax() << " " << RecoHist->GetYaxis()->GetXmin() << " " << RecoHist->GetYaxis()->GetXmax() << " " << RecoHist->GetZaxis()->GetXmin() << " " << RecoHist->GetZaxis()->GetXmax() << std::endl;
                            //std::cout << dE_reco << " " << dTheta_reco << " " << dPhi_reco << std::endl;

                            //std::cout << recodir[0] * 180.0/TMath::Pi() << " | " << truedir[0] * 180.0/TMath::Pi()<< " " < PDFth(theta_reco,E_true, theta_true) << std::endl;
                            //std::cout << i << " " << j << " " << k <<  " " << E_true <<" PDFx" << recodir[0] * 180.0/TMath::Pi()<< " " << recodir[1]* 180.0/TMath::Pi() << " | " << truedir[0] * 180.0/TMath::Pi()<< " " << truedir[1] * 180.0/TMath::Pi()<< " | " <<PDFangle(recodir,truedir,E_true) << std::endl; 
                            
                            std::cout << "Bin size" << std::endl;
                            std::cout << dE_true << " " << dth_true << " " << dphi_true << std::endl;
                            std::cout << dE_reco << " " << dTheta_reco << " " << dPhi_reco<< std::endl;



                            Rlmn += PDFangle(recodir,truedir,E_true) * PDFE(E_reco, E_true) * Rijk* dV_true;
                        }
                    }
                }


                double dV_reco = dE_reco * dTheta_reco * dPhi_reco;

                Nlmn = (MT/mN)*Rlmn*dV_reco; // Reconstructed signal

                //std::cout << "Observed N: " << Nlmn << std::endl;

                RecoHist->SetBinContent(l, m, n, Nlmn);
                
            }
        }
    }


    //std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    //std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   return RecoHist;
}



/*
TH2D*  AsimovObsSimulation::MarginalRecoRates(){  
   // ------------------------------------------------------------------------------------------------------
   // 1. OBSERVED VARIABLES & BINNING SETUP - Binning scheme: Reco Energy(l), Reco Zenith(m), Reco Azimuth(n) 
   // ------------------------------------------------------------------------------------------------------
    
    //int lbins = recobinsE
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsAzi; //Number of Energy bins of Observed event distribution 
    
    //Reconstructed Energy range [GeV]
    double E_reco_min      = EnuMin;
    double E_reco_max      = EnuMax;
    
    //Reconstructed Zenith range (rad)
    double theta_reco_min = ZenMin*TMath::Pi()/180.0;
    double theta_reco_max = ZenMax*TMath::Pi()/180.0;

    //Reconstructed Azimuth range (rad)
    double phi_reco_min = AziMin*TMath::Pi()/180.0 ;
    double phi_reco_max = AziMax*TMath::Pi()/180.0 ;


    //Create bins of variable size

    std::vector<double> ErecoEdges; 
    ErecoEdges.push_back(E_reco_min);
    double Erecoi = E_reco_min;

    while (ErecoEdges.back() < E_reco_max) 
    {
        //Find center
        TF1 *frecoC = new TF1("fmin","-1*x + [0] + 0.5*( [1]*x+[2]*sqrt(x) )", 0, 100);
        frecoC->SetParameter(0,Erecoi);
        frecoC->SetParameter(1,Ae);
        frecoC->SetParameter(2,Be);
        ROOT::Math::RootFinder finder_center;
        finder_center.Solve(*frecoC ,0,100);
        double ErecoC = finder_center.Root();
        
        //Find next edge
        double Erecoii = ErecoC + 0.5*(Ae*ErecoC + Be*sqrt(ErecoC));
        
        ErecoEdges.push_back(Erecoii);
        
        Erecoi = Erecoii;
        
    }

    double ErecobinEdges [ErecoEdges.size()];
        
    for (int i=0 ; i < ErecoEdges.size(); i++ )
    {
          ErecobinEdges[i] = ErecoEdges[i];
            
    }
    
    int lbins = sizeof(ErecobinEdges)/sizeof(ErecobinEdges[0]) - 1; // Number of bins

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsE; 
    int jbins = truebinsZen;
    int kbins = truebinsAzi; 
    
    //True Energy range [GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,E_reco_min);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoEdges.back()); // "Back" indicate upper bound in Energy obtained in variable binning
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double E_true_min = finder_min.Root();
    double E_true_max = finder_max.Root();
   
    // True Zenith Range [deg]
    double theta_true_min = max(0.0,theta_reco_min - 4*(Ath + Bth/sqrt(E_true_min))) *180.0/TMath::Pi(); 
    double theta_true_max = min(TMath::Pi(),theta_reco_max + 4*(Ath + Bth/sqrt(E_true_min))) *180.0/TMath::Pi(); 

    //True Azimuth range [deg]
    double phi_true_min = ( (phi_reco_min) - 4*(Aphi + Bphi/sqrt(E_true_min)) )*180.0/TMath::Pi(); 

    double phi_true_max = ( (phi_reco_max) + 4*(Aphi + Bphi/sqrt(E_true_min)) )*180.0/TMath::Pi(); 

    
    // --------------------------------------------------------------------------------------
    // 3. HISTOGRAM SETUP & RECONSTRUCTION SIMULATION
    // --------------------------------------------------------------------------------------


    //Reconstructed rate distribution:
    
    TH2D* RecoHistMar = new TH2D("RecoHistMar","Azimuth-Averaged Reconstructed rate",lbins,ErecobinEdges,
                                                                mbins,theta_reco_min,theta_reco_max); // xbins correspond to energy values and ybins to zenith angle cost

    // --------------------------------------------------------------------------------------
    // 3. HISTOGRAM SETUP & RECONSTRUCTION SIMULATION
    // --------------------------------------------------------------------------------------

    // Set up true simulation properties
    AsimovSimulation TrueSimulation;
    TrueSimulation.PremTable = ThePremTable;
    TrueSimulation.HondaTable = TheHondaTable;
    TrueSimulation.SetDetectorPosition(xyzTelescope);

    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape = ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight;
    TrueSimulation.aperture = ThePileAperture;
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;

    TrueSimulation.ModifyLayer(PremLayer, DensityContrast, ChemicalContrats);
    TrueSimulation.SetIntervals(theta_true_min, theta_true_max, phi_true_min, phi_true_max, E_true_min, E_true_max);
    TrueSimulation.SetBinning(ibins,jbins,kbins);

    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf = flvf;

    // Retrieve the true histogram
    TH3D* TrueHist = TrueSimulation.GetIntRates();

    //Avarage over azimuth 
    TH2D* TrueHist_avg = new TH2D("TrueHist_avg", "Azimuth-Averaged Interacting Rate",
                        ibins, TrueHist->GetXaxis()->GetXbins()->GetArray(),
                        jbins, TrueHist->GetYaxis()->GetXbins()->GetArray());

    double Nbins = static_cast<double>(TrueHist->GetNbinZ()); // Preferred C++ way
   
    for (int i = 1; i <= TrueHist->GetNbinsX(); ++i) 
    {
        for (int j = 1; j <= TrueHist->GetNbinsY(); ++j) 
        {
            double sum = 0.0;

            for (int k = 1; k <= TrueHist->GetNbinZ(); ++k) 
            {
                sum += TrueHist->GetBinContent(i, j, k);
            }

            double avg = sum / Nbins;
            TrueHist_avg->SetBinContent(i, j, avg);
        }
    }

    // --------------------------------------------------------------------------------------
    // 4. EVENT RECONSTRUCTION SIMULATION
    // --------------------------------------------------------------------------------------

        for (int m = 1; m <= RecoHistMar->GetYaxis()->GetNbins(); ++m) {

            double theta_reco = RecoHistMar->GetYaxis()->GetBinCenter(m);
            double dTheta_reco = RecoHistMar->GetYaxis()->GetBinWidth(m);

            for (int l = 1; l <= RecoHistMar->GetXaxis()->GetNbins(); ++l) {
                double E_reco = RecoHistMar->GetXaxis()->GetBinCenter(l);
                double dE_reco = RecoHistMar->GetXaxis()->GetBinWidth(l);

                double Rlm = 0;

                    for (int j = 1; j <= TrueHist_avg->GetYaxis()->GetNbins(); ++j) {
                        
                        double theta_true = TrueHist_avg->GetYaxis()->GetBinCenter(j);
                        double dtheta_true = TrueHist_avg->GetYaxis()->GetBinWidth(j);

                        for (int i = 1; i <=  TrueHist_avg->GetXaxis()->GetNbins(); ++i) {
                          
                            double E_true = TrueHist_avg->GetXaxis()->GetBinCenter(i);
                            double dE_true = TrueHist_avg->GetXaxis()->GetBinWidth(i);

                            double Rij = TrueHist_avg->GetBinContent(i, j); // Interacting rate

                            double dV_true = dE_true * dtheta_true;
   
                            Rlm += PDFth(theta_reco,E_true, theta_true) * PDFE(E_reco, E_true) * Rij * dV_true;
                        
                        }  // Int Energy
                    }  // Int Zenith
                
                RecoHistMar->SetBinContent(l, m, Rlm);
            } //Reco energy
        } // reco Zenith
   
   return RecoHistMar;}

// Generates the event distribution from a rate distribution:
TH2D* AsimovSimulation::MarginalRecoEvents(TH2D* rateHist) {
   
    if (!ratesHist) return nullptr;

    TH2D* eventsHist = (TH2D*) ratesHist->Clone("Marginalized Reconstructed Events");
    eventsHist->SetTitle("Marginalized Reconstructed Events Histogram");

    for (int i = 1; i <= ratesHist->GetNbinsX(); ++i) {
        double dE = ratesHist->GetXaxis()->GetBinWidth(i);
        for (int j = 1; j <= ratesHist->GetNbinsY(); ++j) {
            double dth = ratesHist->GetYaxis()->GetBinWidth(j);
          
                double Rij = ratesHist->GetBinContent(i, j);
                double Nij = (MT/mN)*Rij * dE * dth;

                eventsHist->SetBinContent(i, j, k, Nij);
        }
    }

    return eventsHist;
}

// Acces event rates do:
TH2D* AsimovSimulation::GetMarginalRecoEvents() {
    TH2D* rateHist = MarginalRecoRates();    // Step 1: Call your rate-calculation method
    return MarginalRecoEvents(rateHist);        // Step 2: Convert rates to events
}

*/

/*
TH2D*  AsimovObsSimulation::GetObsEvents2Dcth(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle
    double threcoMin   = ZenMin;
    double threcoMax   = ZenMax;

    double CthrecoMin = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    double CthrecoMax = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,ErecoMin);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoMax);
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double EtrueMin = finder_min.Root();
    double EtrueMax = finder_max.Root();

    //True Angle-
    double thtrueMin_a= (threcoMin*TMath::Pi()/180.0)-4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= (threcoMax*TMath::Pi()/180.0)+4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;


    std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;


    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location

    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;


    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(NnT);
    TrueSimulation.flvf=flvf;
    TH2D * TrueHist = TrueSimulation.GetTrueEvents2D();

    //Observed distribution:
    TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

   //Calculation of observed events------------------------------------------------------------------------------------- 
    for (int m = 1; m <= mbins; ++m) //Observed Angular bins
    {

        double cthreco = RecoHist->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
        double threco = acos(cthreco);
        double dcthreco = RecoHist->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram
       
        for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
        {
     
            double Ereco = RecoHist->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
            double dEreco = RecoHist->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
       
            double Nsmear = 0;   // Reset total number of events.
            
            for(int i=1; i<= ibins ; i++) //Loop in cosT
            {    
                // Get cos(theta_z) from bin center, It fix a value of L
                double cthtrue = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                double thtrue = acos(cthtrue);
                double dcthtrue = TrueHist->GetXaxis()->GetBinWidth(i); 

                for (int k = 1; k <= kbins; ++k)
                { 
                    //NUMERICAL INTEGRATION

                    double Etrue = TrueHist->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dEtrue = TrueHist->GetYaxis()->GetBinWidth(k); 

                    double NikTrue = TrueHist->GetBinContent(i,k);

                  //  std::cout << "True bin conent: " <<  NikTrue << std::endl;

                    double dEdcthReco = dEreco*dcthreco;
                   
                    //Nsmear += PDFcth(cthreco,Etrue,cthtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdcthReco);

                    Nsmear += PDFE(Ereco,Etrue)*(NikTrue);

    
                } // loop e

            } // Loop ct
            
            double Nmn = Nsmear;

            RecoHist->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 

         
        } // loop in Ereco

    } // loop in CthReco

   std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "] " << phiTrueMin << " " <<  phiTrueMax << std::endl;
   std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   delete TrueHist;
   
   return RecoHist;}
*/
