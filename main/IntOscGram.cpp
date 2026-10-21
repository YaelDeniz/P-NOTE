/*
Author: Yael Deniz
PhD Candidate, University of Idaho

Description: Code that calculates chi2 for a reconstructed signal of a neutrino detector

Notes:
3D code has a Break memory violation if theta is 90 ± 0.1 -> Fixed with a patch
that changes to oscprob prem for that range in zenith. A potential issue appears
if we consider the relative position of the LLVPs-detector
*/

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
    std::string llvpShape = "pancake";

    // Simulation configuration
    int energyBins = 120;
    int zenithBins = 80;
    int azimuthBins = 1;



    double enuMin = 2.0;
    double enuMax = 10.0;

    double zenithMin = 180 - TMath::ASin(depthMax / Rocean) * (180.0 / TMath::Pi());
    double zenithMax = 180 - TMath::ASin(depthMin / Rocean) * (180.0 / TMath::Pi());

    double cosZenithMin = cos(zenithMax * TMath::Pi() / 180.0);
    double cosZenithMax = cos(zenithMin * TMath::Pi() / 180.0);

    double azimuthMin = 89;
    double azimuthMax = 90;
    
    for (int flv = 0; flv < 2; ++flv)
    {
        int nuflv = flv;

    // Initialize Simulation (Standard Earth Model)
        AsimovSimulation experminet;
        experminet.PremTable = premPath;
        experminet.HondaTable = fluxPath;
        experminet.SetDetectorPosition(detectorPos);
        experminet.SetIntervals(zenithMin, zenithMax, azimuthMin, azimuthMax, enuMin, enuMax);
        experminet.SetBinning(energyBins,zenithBins, azimuthBins);
        experminet.SetExposure(exposureMtonYears);
        experminet.flvf = nuflv;

        //Expected
        experminet.MantleAnomaly = false;
        experminet.PremTable = premPath;
        experminet.PileHeight = llvpHeight;
        experminet.aperture = 45;
        experminet.AnomalyShape = "pancake";
        experminet.PileDensityContrast = 0;
        TH3D*Nint_expected = experminet.GetIntEvents();

        //Obseved model
        experminet.PremTable = premPath;
        experminet.MantleAnomaly = true;
        experminet.PileDensityContrast = 2.0;
        TH3D*Nint_Observed = experminet.GetIntEvents();


        for (int n = 1; n <= Nint_Observed->GetZaxis()->GetNbins(); ++n)
        {
            double phi = Nint_Observed->GetZaxis()->GetBinCenter(n) * 180.0 / TMath::Pi();

            int nE = Nint_Observed->GetXaxis()->GetNbins();
            int nTheta = Nint_Observed->GetYaxis()->GetNbins();

            double emin = Nint_Observed->GetXaxis()->GetXmin();
            double emax = Nint_Observed->GetXaxis()->GetXmax();
            double thmin = Nint_Observed->GetYaxis()->GetXmin();
            double thmax = Nint_Observed->GetYaxis()->GetXmax();

            std::string title_std = Form("IntOscGram2D_nstd_phiBin%d", n);
            TH2D* OscGram2D_nstd = new TH2D("IntOscGram2D_nstd", title_std.c_str(), nE, emin, emax, nTheta, thmin, thmax);

             std::string title_alt = Form("IntOscGram2D_nalt_phiBin%d", n);
            TH2D* OscGram2D_nalt = new TH2D("IntOscGram2D_nalt", title_alt.c_str(), nE, emin, emax, nTheta, thmin, thmax);

            for (int l = 1; l <= nE; ++l)
            {
                for (int m = 1; m <= nTheta; ++m)
                {
                    OscGram2D_nstd->SetBinContent(l, m, Nint_expected->GetBinContent(l, m, n));
                    OscGram2D_nalt->SetBinContent(l, m, Nint_Observed->GetBinContent(l, m, n));
                }
            }

            std::string DirPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/";
            std::string OscGramPath = "SimulationResults/Single-Detector-Senv/Oscillograms/";
            std::string OscGramFile = DirPath + OscGramPath + "IntOscPan_" + std::to_string(n) +"_"+std::string(dateBuffer)+ ".csv";

            ExportToCSV(OscGram2D_nstd, OscGram2D_nalt, OscGramFile);

            delete OscGram2D_nstd;
            delete OscGram2D_nalt;

        } // End of azimuth loop

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
