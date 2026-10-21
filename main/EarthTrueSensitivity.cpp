
/* Direction of neutrinos a re given in degrees*/
#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"


# define myPi 3.14159265358979323846  /* pi */

//Some important Earth Radius
# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)

//Some units
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds

using namespace std;

void ExportToCSV(TH2D* hist, std::string filename); //Export File to CSV format


int main(int argc, char **argv)
{
    //Setup the detector
    double R =  Rocean;                // Underwater detector
    double rdet[3] = {0.0,0.0, -1*R }; // Detector location
    double DetMass = 10.0*MTon;        // Mass in megaton units
    double T       = 10.0*years2sec;   // Detector Exposure time in sec: One Year
    double MT = DetMass*T;             // Exposure [Mton*years]

    //Set the neutrino Neutrino Flux
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "GRN_AziAveraged_solmin/grn-nu-20-01-000.d"; //Gran Sasso
    std::string FluxTable = FluxFolder + FluxFile; 

    //Set the Earth models
    std::string PremFolder ="/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_425layers";
    int MaxLayers = 425;
    std::string PremFile = PremName+".txt";
    std::string path2prem = PremFolder+PremFile;
    //Set Layer Density contrats
    int rhopct = 10;

    //std::string PremFolderAlt ="/home/dehy0499/OscProb/PremTables/Prem44pct10/";

    //Simulation Setup

    //Energy interval (in GeV)
    double Emin=1.0 ; 
    double Emax=40.0 ;

    int intE1= 1;
    int intE2= 40;

    //Zenith Angle Interval
    double thmin = 78; // Not Able to use 90 degrees
    double thmax = 180.0;
    double cthmin = cos(thmax*TMath::Pi()/180.0);
    double cthmax = cos(thmin*TMath::Pi()/180.0);

    //Azimuthal angle Interval
    double phimin = 0.0;
    double phimax = 360.0 ; //Whole Earth
    
    // Number of bins
    int cthbins=50; // # Bins in zenith/cos(zenith)
    int abins = 1;
    int ebins=  50; // bins in energy

    //TMultiGraph *Pnu = new TMultiGraph(); // Oscillations to muon neutrinos
    //TMultiGraph *Pnubar = new TMultiGraph(); // Oscillation to electron neutrinos

    // Signal from a Standard Earth model
    AsimovSimulation TrueEvents;
    TrueEvents.HondaTable = FluxTable;
    TrueEvents.SetIntervals(thmin,thmax,phimin,phimax,Emin,Emax);
    TrueEvents.SetBinning(cthbins,abins,ebins);
    TrueEvents.SetDetectorPosition(rdet);
    TrueEvents.SetExposure(MT);
    TrueEvents.PremTable = path2prem;

    //Muon neutrinio like events
    TrueEvents.flvf=1;
    std::vector< TH2D * > TrueStd_mu = TrueEvents.GetTrueEvents3D();

    // Electron netrino like events
    TrueEvents.flvf=0;
    std::vector< TH2D * > TrueStd_e = TrueEvents.GetTrueEvents3D();

    AsimovSimulation AlternativeEarth;
    AlternativeEarth.SetDetectorPosition(rdet);
    AlternativeEarth.HondaTable = FluxTable;
    AlternativeEarth.MantleAnomaly = false;
    AlternativeEarth.AnomalyShape="pancake";
    AlternativeEarth.SetIntervals(thmin,thmax,phimin,phimax,Emin,Emax);
    AlternativeEarth.SetBinning(cthbins,abins,ebins);
    AlternativeEarth.SetExposure(MT);
    

    TGraph * chi2plot = new TGraph(MaxLayers);

    for (int i = 0; i < MaxLayers; ++i)
    {
        
    int layer = i+1;


    AlternativeEarth.PremTable = path2prem;
    AlternativeEarth.ModifyLayer(layer,10,0.0);
    
    AlternativeEarth.flvf=1;    
    std::vector< TH2D * > TrueAlt_mu = AlternativeEarth.GetTrueEvents3D();
    double chi2mu = Get2DChi2(TrueStd_mu[0],TrueAlt_mu[0]);

    
    AlternativeEarth.flvf=0;
    std::vector< TH2D * > TrueAlt_e = AlternativeEarth.GetTrueEvents3D();
    double chi2e = Get2DChi2(TrueStd_e[0],TrueAlt_e[0]);
    
    
    double chi2total = chi2e+chi2mu;

    //chi2file << std::setprecision(10) << i+1 << " " << chi2e <<  " " << chi2mu << chi2total << std::endl;

    chi2plot->SetPoint(i,chi2total,i+1);


    }

    TApplication app("app", &argc, argv);
    TCanvas *c1 = new TCanvas();
    //TColor::InvertPalette();
    chi2plot->Draw("apl");
    gPad->SetLogx();
    c1->Modified(); 
    c1->Update();
    TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();

    app.Run();

    
    
    /*
    TH2D * hist = new TH2D("EventHist", "Charge Current events", ebins,Emin,Emax,cthbins,cthmin,cthmax);
    for (int i = 1; i <=cthbins; ++i)
    {   
        std::cout << "*-*"<< std::endl;
        for (int j = 1; j <= ebins; ++j)
        {
            double N = TrueStd_mu[4]->GetBinContent(i,j);
            hist->SetBinContent(j,i, N);
        }
    }
    
    TApplication app("app", &argc, argv);
    TCanvas *c1 = new TCanvas();
    gStyle->SetPalette( kInvertedDarkBodyRadiator);
    //TColor::InvertPalette();
    hist->Draw("COLZ");
    hist->SetStats(0);
    gPad->SetLogx();
    c1->Modified(); 
    c1->Update();
    TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();

    TCanvas *c2 = new TCanvas();

    c2->Divide(2,2);

    c2->cd(1);
    //gStyle->SetPalette(kSolar);
    //TColor::InvertPalette();
    TrueStd_e[0]->Draw("COLZ");
    TrueStd_e[0]->SetStats(0);
    c2->Modified(); 
    c2->Update();
    
    c2->cd(2);
    //gStyle->SetPalette(kSolar);
    //TColor::InvertPalette();
    TrueStd_e[1]->Draw("COLZ");
    TrueStd_e[1]->SetStats(0);
    c2->Modified(); 
    c2->Update();
    
    c2->cd(3);
    //gStyle->SetPalette(kSolar);
    //TColor::InvertPalette();
    TrueStd_e[2]->Draw("COLZ");
    TrueStd_e[2]->SetStats(0);
    c2->Modified(); 
    c2->Update();
    
    c2->cd(4);
    //gStyle->SetPalette(kSolar);
    //TColor::InvertPalette();
    TrueStd_e[3]->Draw("COLZ");
    TrueStd_e[3]->SetStats(0);
    c2->Modified(); 
    c2->Update();

    TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();

    app.Run();

    */


  return 0;

}

void ExportToCSV(TH2D* hist, std::string filename)
{
   
    std::ofstream outfile(filename);

    // Recorre los bins y guarda el contenido
    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= hist->GetNbinsY(); ++j)
        {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);
            double content = hist->GetBinContent(i, j);

            outfile << x << "," << y << "," << content << std::endl;
        }
    }

    outfile.close();
}