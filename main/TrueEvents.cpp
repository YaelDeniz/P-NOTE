
/* Direction of neutrinos a re given in degrees*/
#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>


//Cern ROOT
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "PlotTool.h"
#include "StatsAna.h"
#include "GetTrueEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need

# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds

using namespace std;

void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename);

int main(int argc, char **argv)
{

    // Constants
    double detectorMassMT = 10.0;              // Mass in megaton units
    double detectorExposureYears = 10.0;       // Exposure time in years
    double exposureMtonYears = detectorMassMT * MTon * detectorExposureYears * years2sec; // Exposure in [Mton*years]

    // Detector Position
    double detectorRadius = Rearth;            // Radius at South Pole
    double detectorPos[3] = {0.0, 0.0, -1.0 * detectorRadius};

    // Paths for input and output files
    std::string fluxPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string premPath = "/home/dehy0499/OscProb/PremTables/prem_44layers.txt";

    // LLVP Configuration
    double llvpHeight = 1000;                  // LLVP height in km
    double llvpRadius = Rcmb + llvpHeight;
    double depthMin = Rcmb - 2500;
    double depthMax = llvpRadius + 500;        // Max distance from Earth's center
    double llvpDensityContrast = 2.0;          // LLVP density contrast percentage
    std::string llvpShape = "pancake";

    // Simulation Configuration
    int zenithBins = 100;                      // Zenith bins
    int azimuthBins = 100;                     // Azimuth bins
    int energyBins = 100;                      // Energy bins

    // Energy [GeV]
    double enuMin = 1.0;
    double enuMax = 10.0;

    // Zenith angle [degrees]
    double zenithMin = 180 - TMath::ASin(depthMax / detectorRadius) * (180.0 / TMath::Pi());
    double zenithMax = 180 - TMath::ASin(depthMin / detectorRadius) * (180.0 / TMath::Pi());

    double cosZenithMin = cos(zenithMax * TMath::Pi() / 180.0);
    double cosZenithMax = cos(zenithMin * TMath::Pi() / 180.0);

    // Azimuthal angle [degrees]
    double azimuthStart = -50;
    double azimuthEnd = 50;
    double dAzimuth = (azimuthEnd - azimuthStart) / (2.0 * azimuthBins);
    double azimuthMin = azimuthStart - dAzimuth;
    double azimuthMax = azimuthEnd + dAzimuth;

    // Neutrino flavor (0: nue, 1: numu, 2: nutau)
    int nuflv = 0;

    // File and folder structure configuration
    int exposureYearsLabel = static_cast<int>(detectorMassMT * detectorExposureYears);
    std::string BinLabel = std::to_string(zenithBins)+"Zen"+ std::to_string(azimuthBins)+"Az"+ std::to_string(energyBins)+"Enu";
    std::string SimLabel = "Simulation_"+llvpShape+"_"+std::to_string(exposureYearsLabel)+"Mton"+"_"+ std::to_string(static_cast<int>(enuMin)) +"-"+ 
                                 std::to_string(static_cast<int>(enuMax))+"GeV_" + std::to_string(static_cast<int>(zenithMin)) +"-"+ 
                                 std::to_string(static_cast<int>(zenithMax))+"Zen_"  + std::to_string(static_cast<int>(azimuthMin)) +"-"+ 
                                 std::to_string(static_cast<int>(azimuthMax))+"Az";
    
    std::string pathToResults = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/PreliminaryResults/";
    std::string pathToEvents = pathToResults + "IntEvents/";
    std::string eventFolderPath = pathToEvents + "Events_" + SimLabel + "/";
    std::string pathToChiSquares = pathToResults + "IntChi2/";
    std::string chiSquareFolderPath = pathToChiSquares + "Chi2_" + SimLabel + "/";

    std::string chiSquareFileName = "Int_cth_Chi2LLVP" + std::to_string(exposureYearsLabel) + "myrs_" + llvpShape + 
                                          "nu" + std::to_string(nuflv) + "_" + BinLabel + ".csv";

    // Create output directories if they do not exist
    mkdir(eventFolderPath.c_str(), 0777);
    mkdir(chiSquareFolderPath.c_str(), 0777);

    // Initialize Simulation (Standard Earth Model)
    AsimovSimulation standardEarth;
    standardEarth.PremTable = premPath;
    standardEarth.HondaTable = fluxPath;
    standardEarth.SetDetectorPosition(detectorPos);
    standardEarth.MantleAnomaly = false;

    // Set intervals, binning, and exposure for the simulation
    standardEarth.SetIntervals(zenithMin, zenithMax, azimuthMin, azimuthMax, enuMin, enuMax);
    standardEarth.SetBinning(zenithBins, azimuthBins, energyBins);
    standardEarth.SetExposure(exposureMtonYears);
    standardEarth.flvf = nuflv;

    // Run simulation and generate events for Standard and Alternative

    std::vector<TH2D*> trueEventsStandard = standardEarth.GetTrueEvents3D();

    // Configure for Anomalous Earth Model (with LLVP)
    standardEarth.MantleAnomaly = true;
    standardEarth.PileHeight = llvpHeight;
    standardEarth.aperture = 45;
    standardEarth.AnomalyShape = llvpShape;
    standardEarth.PileDensityContrast = llvpDensityContrast;

    // Generate events for the Anomalous Earth model
    std::vector<TH2D*> trueEventsAnomalous = standardEarth.GetTrueEvents3D();

    // Exporting Events
    for (int nhist = 0; nhist < trueEventsAnomalous.size(); ++nhist)
    {
        std::string eventFileName = "IntcthLLVP_nu" + std::to_string(nuflv) + "_"+ BinLabel + 
                                          "_" + std::to_string(nhist) +".csv";
        std::string eventFilePath = eventFolderPath + eventFileName;

        TH2D* oscDifference = new TH2D(Form("OscDiff2D%d", nhist), Form("Oscillogram%d", nhist), zenithBins, cosZenithMin, cosZenithMax, energyBins, enuMin, enuMax);
        
        // Compute differences and export to CSV
        GetDiff2D(trueEventsStandard[nhist], trueEventsAnomalous[nhist], oscDifference);
        ExportToCSV(trueEventsStandard[nhist], trueEventsAnomalous[nhist], oscDifference, eventFilePath);
    }

    // Simulation completed
    std::cout << "Simulation completed and events generated." << std::endl;

    return 0;
}


void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename)
{
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";

    //std::string Evtfname  = EvtFolder.c_str()+chi2namecth;

    //std::string path2file  = NuTomoPath+IntFolder+filename;
   
    std::ofstream outfile(filename);

    // Recorre los bins y guarda el contenido


    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= hist->GetNbinsY(); ++j)
        {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);

            double nexp = stdhist->GetBinContent(i, j);
            double nobs = althist->GetBinContent(i, j); 
            double diff = hist->GetBinContent(i, j);

            outfile << x << "," << y << "," <<  nexp <<" , " << nobs << " , " << diff << std::endl;
        }
    }

    outfile.close();
}

/*
  // Visualization of Events

   TApplication app("app", &argc, argv);

   //cth lines
   //double cthbottom = -0.8376; 
   double cthbottom = TMath::Cos(TMath::Pi()-TMath::ASin( 3480.0/detectorRadius ));
   double cthmid_bottom = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+300)/detectorRadius ));;
   double cthmid_top = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+600)/detectorRadius ));;
   double cthtop = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/detectorRadius ));;

   std::cout << cthtop << std::endl;

   TLine * lbottom = new TLine(cthbottom,enuMin,cthbottom,enuMax);
   TLine * lmid_bottom = new TLine(cthmid_bottom,enuMin,cthmid_bottom,enuMax);
   TLine * lmid_top = new TLine(cthmid_top,enuMin,cthmid_top,enuMax);
   TLine * ltop = new TLine(cthtop,enuMin,cthtop,enuMax);

   lbottom -> SetLineColor(kRed);
   lmid_bottom -> SetLineColor(kRed);
   lmid_top -> SetLineColor(kRed);
   ltop -> SetLineColor(kRed);

   lbottom -> SetLineStyle(2);
   lmid_bottom -> SetLineStyle(2);
   lmid_top -> SetLineStyle(2);
   ltop -> SetLineStyle(2);

   //th lines 
   double thbottom     = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( 3480.0/detectorRadius ));
   double thmid_bottom = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 300)/detectorRadius ));
   double thmid_top    = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 600)/detectorRadius ));
   double thtop        = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/detectorRadius ));

   TLine * lbottomth = new TLine(thbottom,enuMin,thbottom,enuMax);
   TLine * lmid_bottomth = new TLine(thmid_bottom,enuMin,thmid_bottom,enuMax);
   TLine * lmid_topth = new TLine(thmid_top,enuMin,thmid_top,enuMax);
   TLine * ltopth = new TLine(thtop,enuMin,thtop,enuMax);

   lbottomth -> SetLineColor(kBlue);
   lmid_bottomth -> SetLineColor(kBlue);
   lmid_topth -> SetLineColor(kBlue);
   ltopth -> SetLineColor(kBlue);

 
   TCanvas *c1 = new TCanvas();
   gStyle->SetPalette(kBlueGreenYellow);
   
   TCanvas *c2 = new TCanvas();
   gStyle->SetPalette(kBlueGreenYellow);


   c1->Divide(5,5);
   c2->Divide(5,5);



c1->cd(i+1);
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    lbottom->Draw("same");
    lmid_bottom->Draw("same");
    lmid_top->Draw("same");
    ltop->Draw("same");
    

    c2->cd(i+1);
    TrueDiff2Dth->Draw("COLZ");
    TrueDiff2Dth->SetStats(0);
    lbottomth->Draw("same");
    lmid_bottomth->Draw("same");
    lmid_topth->Draw("same");
    ltopth->Draw("same");


   gPad->Update();
   c1->Modified(); c1->Update();
   TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
   rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

   gPad->Update();
   c2->Modified(); c2->Update();
   TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();
   rc2->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

   app.Run();

    //c1->Print("AzimuthHist.png");




*/