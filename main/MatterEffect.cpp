
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
    std::string premPath = "/home/dehy0499/OscProb/PremTables/";
    
    std::string earthModel = "prem_44layers";
    //std::string PremHeavy  = premPath+"Prem_heavy.txt";
    //std::string PremSimple = premPath+"Prem_simple.txt";

    std::string PremAlt  = premPath+earthModel+".txt";
    
    std::string PremVac    = premPath+"prem_vac.txt";

    // LLVP Configuration
   
    // Simulation Configuration
    int zenithBins = 1000;                      // Zenith bins
    int azimuthBins =0;                     // Azimuth bins
    int energyBins = 1000;                      // Energy bins

    // Energy [GeV]
    double enuMin = 1.0;
    double enuMax = 10.0;

    // Zenith angle [degrees]
    double zenithMin = 160;
    double zenithMax = 165;

    double cosZenithMin = cos(zenithMax * TMath::Pi() / 180.0);
    double cosZenithMax = cos(zenithMin * TMath::Pi() / 180.0);

    // Azimuthal angle [degrees]
    double azimuthStart = 0;
    double azimuthEnd   = 360;
    double dAzimuth = (azimuthEnd - azimuthStart) / (2.0 * azimuthBins);
    double azimuthMin = azimuthStart;
    double azimuthMax = azimuthEnd;

    // Neutrino flavor (0: nue, 1: numu, 2: nutau)
    int nuflv = 1;

    // File and folder structure configuration
    int exposureYearsLabel = static_cast<int>(detectorMassMT * detectorExposureYears);
    std::string BinLabel = std::to_string(zenithBins)+"Zen"+ std::to_string(azimuthBins)+"Az"+ std::to_string(energyBins)+"Enu";
    std::string SimLabel = "MSWExamples_"+std::to_string(exposureYearsLabel)+"Mton"+"_"+ std::to_string(static_cast<int>(enuMin)) +"-"+ 
                                 std::to_string(static_cast<int>(enuMax))+"GeV_" + std::to_string(static_cast<int>(zenithMin)) +"-"+ 
                                 std::to_string(static_cast<int>(zenithMax))+"Zen_"  + std::to_string(static_cast<int>(azimuthMin)) +"-"+ 
                                 std::to_string(static_cast<int>(azimuthMax))+"Az";
    
    std::string pathToResults = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/PreliminaryResults/";
    std::string pathToEvents = pathToResults + "MatterExample/";
    std::string eventFolderPath = pathToEvents + "Events_" + SimLabel + "/";


    // Create output directories if they do not exist
    mkdir(eventFolderPath.c_str(), 0777);


    // Initialize Simulation (Standard Earth Model)
    AsimovSimulation Simulation;
    Simulation.PremTable = PremVac;
    Simulation.HondaTable = fluxPath;
    Simulation.SetDetectorPosition(detectorPos);
    Simulation.MantleAnomaly = false;

    // Set intervals, binning, and exposure for the simulation
    Simulation.SetIntervals(zenithMin, zenithMax, azimuthMin, azimuthMax, enuMin, enuMax);
    Simulation.SetBinning(zenithBins, azimuthBins, energyBins);
    Simulation.SetExposure(exposureMtonYears);
    Simulation.flvf = nuflv;

    // Run simulation and generate events for Standard and Alternative

    std::cout << "Before events ********************************************" << std::endl;


    std::vector<TH2D*> trueEventsStandard = Simulation.GetTrueEvents3D();

   std::cout << "No problem with the Events ********************************************" << std::endl;


    // Configure for Anomalous Earth Model (with LLVP)
    Simulation.PremTable = PremAlt;
    // Generate events for the Anomalous Earth model
    std::vector<TH2D*> trueEventsAnomalous = Simulation.GetTrueEvents3D();

    std::cout << "No problem with the Events ********************************************" << std::endl;


    // Exporting Events
    for (int nhist = 0; nhist < trueEventsAnomalous.size(); ++nhist)
    {
        std::string eventFileName = "MSWExamples_"+earthModel+"_" + std::to_string(nuflv) + "_"+ BinLabel + 
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

