
/* Direction of neutrinos a re given in degrees*/
#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <ctime>
#include <chrono>


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

# define CMB 3480.0 // Core-Mantle-Boundary Radius (km)
# define OCEAN_BOTTOM 6368.0 // Earth radius up to oceanic crust in (km)
# define CONTIENTAL_CRUST 6371.0 // Earth radius in (km)
# define ATM_RADIUS 6386.0 // Radius of the Atmosphere (km)
# define NUCLEON_MASS   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MEGATON_KG  1E9  //Metric MegaTon in kg
# define SECONDS_PER_YEAR 3.154E7 // Years in Seconds

using namespace std;

void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename);
//std::string SimulationID();

int main(int argc, char **argv)
{   

    auto start = std::chrono::high_resolution_clock::now();

    // Get the current date and time; Format the date (YYYY-MM-DD)
    std::time_t now = std::time(0);
    std::tm* localTime_unit = std::localtime(&now);
    char datetimeBuffer[21];
    char dateBuffer[21];

    std::strftime(datetimeBuffer, sizeof(datetimeBuffer), "%Y-%m-%d_%H%M%S", localTime_unit);
    std::strftime(dateBuffer, sizeof(dateBuffer), "%Y-%m-%d", localTime_unit);


    // Detector parameters
    double Mass_unit = 10.0;              // Mass in megaton units
    double Time_unit = 10.0;       // Exposure time in years
    double Exposure = (Mass_unit * MEGATON_KG )*(Time_unit * SECONDS_PER_YEAR); // Exposure in [Mton*years]

    // Detector Position
    double Detector_location = CONTIENTAL_CRUST;            // Radius at South Pole
    double Coordinates[3] = {0.0, 0.0, -1.0 * Detector_location};

    // Paths for input and output files
    std::string fluxPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string premPath = "/home/dehy0499/OscProb/PremTables/prem_44layers.txt";

    // LLVP Configuration
    double h_llvp = 1000;                  // LLVP height in km
    double r_llvp = CMB + h_llvp;          // LLVP radii
    double depthMin = CMB - 2500;
    double depthMax = r_llvp + 500;        // Max distance from Earth's center
    //double llvpDensityContrast = 2.0;          // LLVP density contrast percentage
    std::string mod_llvp = "cake";

    // Simulation Configuration
    int zenithBins = 40;                      // Zenith bins
    int energyBins = 80;                      // Energy bins
    int azimuthBins =100;                     // Azimuth bins
    

    // Energy [GeV]
    double enuMin = 1.0;
    double enuMax = 10.0;

    // Zenith angle [degrees]
    double zenithMin = 180 - TMath::ASin(depthMax / Detector_location) * (180.0 / TMath::Pi());
    double zenithMax = 180 - TMath::ASin(depthMin / Detector_location) * (180.0 / TMath::Pi());

    double cosZenithMin = cos(zenithMax * TMath::Pi() / 180.0);
    double cosZenithMax = cos(zenithMin * TMath::Pi() / 180.0);

    // Azimuthal angle [degrees]
    double azimuthStart = 90-50;
    double azimuthEnd = 90+50;

    double azimuthMin = azimuthStart;
    double azimuthMax = azimuthEnd;

    //double dAzimuth = (azimuthEnd - azimuthStart) / (2.0 * azimuthBins);
    //double azimuthMin = azimuthStart - dAzimuth;
    //double azimuthMax = azimuthEnd + dAzimuth;

    //double rhosmain[3] = {1.0,2.0,3.0};

    //for (int k = 0; k < 1; ++k)
    //{
        
    //int k = 1;
    //double llvpDensityContrast = rhosmain[k]; // 2% density contrats for LLVP

    std::string testVariable = "PLMtest";
    std::string binVariable = "th";

    std::string BinLabel = binVariable+std::to_string(zenithBins)+"-"+ std::to_string(azimuthBins)+"-"+ std::to_string(energyBins);
    std::string SimResultsDir = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/Single-Detector-Senv/";
    std::string ParSenvDir = SimResultsDir + "ne_senv/";
    std::string ExpSenvDir = ParSenvDir+"IntSenv-"+testVariable+"_"+std::string(dateBuffer)+"/";
    mkdir(ExpSenvDir.c_str(), 0777);
    
    for (int nuflv = 0; nuflv < 2; ++nuflv)
    {
        
    
             // File and folder structure configuration
            //int exposureYearsLabel = static_cast<int>(Mass_unit * Time_unit);

            //std::string BinLabel = "Obs"+std::to_string(zenithBinsReco)+"Zen"+ std::to_string(azimuthBinsReco)+"Az"+ std::to_string(energyBinsReco)+"Enu_"
                                //+"Int"+std::to_string(zenithBins)+"Zen"+ std::to_string(azimuthBins)+"Az"+ std::to_string(energyBins)+"Enu";

    std::string ExpSenv = ExpSenvDir + "IntSenv-flv"+std::to_string(nuflv)+"-"+BinLabel+"_"+std::string(datetimeBuffer)+".csv";
        // Create output directories if they do not exist
        

    // Initialize Simulation (Standard Earth Model)
    AsimovSimulation experminet;
    experminet.PremTable = premPath;
    experminet.HondaTable = fluxPath;
    experminet.SetDetectorPosition(Coordinates);
    experminet.SetIntervals(zenithMin, zenithMax, azimuthMin, azimuthMax, enuMin, enuMax);
    experminet.SetBinning(energyBins,zenithBins, azimuthBins);
    experminet.SetExposure(Exposure);
    experminet.flvf = nuflv;

    // Run simulation and generate events for Standard and Alternative


    // Configure for Expected Earth Model 
    experminet.MantleAnomaly = false;
    experminet.PileHeight = h_llvp;
    experminet.aperture = 45;
    experminet.AnomalyShape = mod_llvp;
    experminet.PileDensityContrast = 2;
    //std::vector<TH2D*> Nint_expected = experminet.GetTrueEvents3D();

    //Write experiment details
    std::string DetailsFile = ExpSenvDir + "IntExp_details_"+datetimeBuffer+".txt";
    experminet.GetDetails(DetailsFile);

    TH3D*Nint_expected = experminet.GetIntEvents();
    //TH2D* Nint_ExpMar = experminet.GetAzimuthAvaraged(Nint_expected);

    std::ofstream SenvData(ExpSenv);
    // Run simulation and generate events for Standard and Alternative

    double  rho[9] = {-4.0,-3.0,-2.0,-1,0,1.0, 2.0, 3.0,4.0};
            
    double  h[8] = {100,200,300,400,500,600,700,800};

    double azi[5] ={15,25,35,45,50}; 

    for (int i = 0; i < 9; ++i)
    {
        

        // Configure for Observed Earth Model 

        experminet.PremTable = premPath;
        experminet.MantleAnomaly = true;

        //experminet.PileHeight = h[i]; //km
        experminet.PileDensityContrast = rho[i];
        //experminet.aperture = azi[i];

        // Generate events for the Anomalous Earth model
        //std::vector<TH2D*> Nint_Observed = experminet.GetTrueEvents3D();
        TH3D*Nint_Observed = experminet.GetIntEvents();
        //TH2D* Nint_ObsMar = experminet.GetAzimuthAvaraged(Nint_Observed);

        double chi2tot = 0;
        //for (int n = 0; n < Nint_Observed.size(); ++ n)
        //{   
        //chi2tot += IntChi2( Nint_expected, Nint_Observed);
        //}
        
        chi2tot = IntChi2(Nint_expected, Nint_Observed);
        //chi2tot = Get2DChi2(Nint_ExpMar, Nint_ObsMar);
        SenvData << i <<" , "<< experminet.aperture << " , "<< experminet.PileDensityContrast << " , "  <<    experminet.PileHeight << " , " << chi2tot  <<  std::endl; 

    }

    SenvData.close();
    //std::string DetailsFileName = SimulationID();


    }//Loop over flavors

   // Record the end time
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    std::chrono::duration<double> duration = end - start;

    // Print the execution time in seconds
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    std::cout << "rates test" << std::endl;



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
std::string SimulationID() {
    // Get current time
    std::time_t now = std::time(nullptr);
    std::tm localTime;
    localtime_r(&localTime, &now);  // Use localtime_r(&now, &localTime) on Linux/Mac

    // Format: "Simulation_YYYYMMDD_HHMMSS.txt"
    char filename[100];
    std::strftime(filename, sizeof(filename), "Nint_details_%Y%m%d_%H%M%S.txt", &localTime);

    return std::string(filename);
}
*/


/*
  // Visualization of Events

   TApplication app("app", &argc, argv);

   //cth lines
   //double cthbottom = -0.8376; 
   double cthbottom = TMath::Cos(TMath::Pi()-TMath::ASin( 3480.0/Detector_location ));
   double cthmid_bottom = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+300)/Detector_location ));;
   double cthmid_top = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+600)/Detector_location ));;
   double cthtop = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/Detector_location ));;

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
   double thbottom     = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( 3480.0/Detector_location ));
   double thmid_bottom = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 300)/Detector_location ));
   double thmid_top    = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 600)/Detector_location ));
   double thtop        = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/Detector_location ));

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