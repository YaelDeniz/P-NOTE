/*
Author: Yael Deniz
PhD Candidate, University of Idaho

Description: Code that calculates chi2 for a reconstructed signal of a neutrino detector

Notes:
3D code has a Break memory violation if theta is 90 ± 0.1 -> Fixed with a patch
that changes to oscprob prem for that range in zenith. A potential issue appears
if we consider the relative position of the LLVPs-detector
*/

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

// CERN ROOT
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

// My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetObservedEvents.h"

#define Rcmb     3480.0     // Core-Mantle-Boundary Radius (km)
#define Rocean   6368.0     // Earth radius (km)
#define Rearth   6371.0     // Earth radius (km)
#define Ratm     6386.0     // Radius of the Atmosphere (km)
#define mN       1.67492749804E-27  // Nucleon mass (kg)
#define MTon     1E9        // Metric MegaTon in kg
#define years2sec 3.154E7   // Years in seconds

using namespace std;

void ExportToCSV(TH2D* stdhist, TH2D* althist, std::string filename);

int main(int argc, char **argv)
{
    auto start = std::chrono::high_resolution_clock::now();

    // Get current date and time
    std::time_t now = std::time(0);
    std::tm* localTime_unit = std::localtime(&now);
    char datetimeBuffer[21];
    char dateBuffer[21];

    std::strftime(datetimeBuffer, sizeof(datetimeBuffer), "%Y-%m-%d_%H%M%S", localTime_unit);
    std::strftime(dateBuffer, sizeof(dateBuffer), "%Y-%m-%d", localTime_unit);

    // Detector settings
    double detectorMass = 10.0; // Mass in megatons
    double detectorExposureYears = 10.0;
    double exposureMtonYears = detectorMass * MTon * detectorExposureYears * years2sec;

    double detectorRadius = Rearth;
    double detectorPos[3] = {0.0, 0.0, -1.0 * detectorRadius};

    // File paths
    std::string fluxPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string premPath = "/home/dehy0499/OscProb/PremTables/prem_44layers.txt";

    // Pile setup
    double llvpHeight = 1000; // km
    double llvpRadius = Rcmb + llvpHeight;
    double depthMin = Rcmb - 2500;
    double depthMax = llvpRadius + 500;
    double Aperture = 45;
    double llvpDensityContrast = 0.0;
    std::string llvpShape = "cake";

    // Simulation configuration
    int energyBins = 120;
    int zenithBins = 80;
    int azimuthBins = 1;

    int energyBinsReco = 11;
    int zenithBinsReco = 10;
    int azimuthBinsReco = 1;

    double enuMin = 2.0;
    double enuMax = 10.0;

    double zenithMin = 180 - TMath::ASin(depthMax / Rocean) * (180.0 / TMath::Pi());
    double zenithMax = 180 - TMath::ASin(depthMin / Rocean) * (180.0 / TMath::Pi());

    double cosZenithMin = cos(zenithMax * TMath::Pi() / 180.0);
    double cosZenithMax = cos(zenithMin * TMath::Pi() / 180.0);

    //double azimuthMin = 90 - 50;
    //double azimuthMax = 90 + 50;
    double azimuthMin = 89;
    double azimuthMax = 90;
    

    std::string binVariable = "th";
    std::string testVariable = "PPtest-" + llvpShape + "-" + binVariable;

    std::string BinLabel = "ObsOsc" + std::to_string(zenithBinsReco) + "-" + 
                           std::to_string(azimuthBinsReco) + "-" + 
                           std::to_string(energyBinsReco) + "Int" +
                           std::to_string(zenithBins) + "-" + 
                           std::to_string(azimuthBins) + "-" + 
                           std::to_string(energyBins);

    std::string SimResultsDir = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/Single-Detector-Senv/";
    std::string ParSenvDir = SimResultsDir + "Oscillograms/";
    std::string ExpSenvDir = ParSenvDir + "ObsOsc_" + testVariable + "_" + std::string(dateBuffer) + "/";
    //mkdir(ExpSenvDir.c_str(), 0777);

    for (int flv = 0; flv < 2; ++flv)
    {
        int nuflv = flv;

        AsimovObsSimulation Simulation;
        Simulation.ThePremTable = premPath;
        Simulation.TheHondaTable = fluxPath;
        Simulation.SetDetectorXYZ(detectorPos);
        Simulation.SetIntervals(zenithMin, zenithMax, azimuthMin, azimuthMax, enuMin, enuMax);
        Simulation.SetTrueBinning(energyBins, zenithBins, azimuthBins);
        Simulation.SetRecoBinning(energyBinsReco, zenithBinsReco, azimuthBinsReco);
        Simulation.SetExposure(exposureMtonYears);
        Simulation.flvf = nuflv;

        Simulation.PileInModel = false;
        Simulation.ShapeOfPile = llvpShape;
        Simulation.ThePileHight = llvpHeight;
        Simulation.ThePileAperture = Aperture;
        Simulation.ThePileDensityContrast = llvpDensityContrast;
        //Simulation.ThePileChemicalContrast = 0.0;

        TH3D* ObservedEventsStandard = Simulation.GetRecoN();

        Simulation.ThePremTable = premPath;
        Simulation.PileInModel = true;
        Simulation.ShapeOfPile = llvpShape;
        Simulation.ThePileDensityContrast = 2.0;

        TH3D* ObservedEventsAlternative = Simulation.GetRecoN();

        for (int n = 1; n <= ObservedEventsStandard->GetZaxis()->GetNbins(); ++n)
        {
            double phi = ObservedEventsStandard->GetZaxis()->GetBinCenter(n) * 180.0 / TMath::Pi();

            int nE = ObservedEventsStandard->GetXaxis()->GetNbins();
            int nTheta = ObservedEventsStandard->GetYaxis()->GetNbins();

            double emin = ObservedEventsStandard->GetXaxis()->GetXmin();
            double emax = ObservedEventsStandard->GetXaxis()->GetXmax();
            double thmin = ObservedEventsStandard->GetYaxis()->GetXmin();
            double thmax = ObservedEventsStandard->GetYaxis()->GetXmax();

            std::string title_std = Form("OscGram2D_nstd_phiBin%d", n);
            TH2D* OscGram2D_nstd = new TH2D("OscGram2D_nstd", title_std.c_str(), nE, emin, emax, nTheta, thmin, thmax);

             std::string title_alt = Form("OscGram2D_nalt_phiBin%d", n);
            TH2D* OscGram2D_nalt = new TH2D("OscGram2D_nalt", title_alt.c_str(), nE, emin, emax, nTheta, thmin, thmax);

            for (int l = 1; l <= nE; ++l)
            {
                for (int m = 1; m <= nTheta; ++m)
                {
                    OscGram2D_nstd->SetBinContent(l, m, ObservedEventsStandard->GetBinContent(l, m, n));
                    OscGram2D_nalt->SetBinContent(l, m, ObservedEventsAlternative->GetBinContent(l, m, n));
                }
            }

            std::string DirPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/";
            std::string OscGramPath = "SimulationResults/Single-Detector-Senv/Oscillograms/";
            std::string OscGramFile = DirPath + OscGramPath + "ObsOsc_" + std::to_string(n) +"_"+std::string(dateBuffer)+ ".csv";

            ExportToCSV(OscGram2D_nstd, OscGram2D_nalt, OscGramFile);

            delete OscGram2D_nstd;
            delete OscGram2D_nalt;

        } // End of azimuth loop

        //std::string DetailsFile = ExpSenvDir + "ObsExp_details_" + datetimeBuffer + ".txt";
        //Simulation.GetDetails(DetailsFile);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Execution time: " << duration.count() << "Application finish running " << std::endl;
    return 0;
}

void ExportToCSV(TH2D* stdhist, TH2D* althist, std::string filename)
{
    std::ofstream outfile(filename);

    for (int i = 1; i <= stdhist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= stdhist->GetNbinsY(); ++j)
        {
            double x = stdhist->GetXaxis()->GetBinCenter(i);
            double y = stdhist->GetYaxis()->GetBinCenter(j);

            double nexp = stdhist->GetBinContent(i, j);
            double nobs = althist->GetBinContent(i, j);
            double ndiff = (nexp != 0) ? 100.0 * ((nobs - nexp) / nexp) : 0.0;


            outfile << x << "," << y << "," << nexp << "," << nobs << "," << ndiff << std::endl;
        }
    }

    outfile.close();
}
