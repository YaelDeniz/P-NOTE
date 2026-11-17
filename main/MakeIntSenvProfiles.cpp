
/* Direction of neutrinos a re given in degrees*/
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <filesystem>
#include <string>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>

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
#include "InteractingEvents.h"



using namespace std;
namespace fs = std::filesystem;

void RecordIntSenvSimInputs(
    const fs::path& PATH_fluxtables,
    const fs::path& PATH_premtables,
    const std::string& geom_mllvp,
    double M_det, double T_det,
    double R_det, const double Location_det[3],
    int nbins_th, int nbins_phi, int nbins_enu,
    double enu_min, double enu_max,
    double th_min, double th_max,
    double phi_lowedge, double phi_upedge,
    double h_mllvp, double R_mllvp,
    double rhopct_mllvp_min, double rhopct_mllvp_max, double rhopct_mllvp_step,
    fs::path& outputDir
);

int main(int argc, char **argv)
{
    // As discussed in section 3 of Deniz et al, a oscillation tomography simulation requires the following inputs:
    // 1. An Earth model, here constructed from PREM tables (Joao Coelho, 2020)
    // 2. A neutirno flux model, generated from here from differential neutrino flux tables (Honda et al. 2015)
    // 3. Detector parameters and exposure, in this specific case, a perfect detector with 10 Mton mass and 10 years of exposure.

    // Paths for tables defining the neutrino flux and the Earth's matter profile

    fs::path PATH_fluxtables =
    fs::current_path() /
    "fluxtables" /
    "SP_AziAveraged_solmin" / // Differential flux tables for the South Pole (SP), averaged over azimuthal angles(AziAveraged), and for solar minimum conditions (solmin).
    "spl-nu-20-01-000.d";

    std::cout << "Using flux tables from: " << PATH_fluxtables << std::endl;
    fs::path PATH_premtables =
    fs::current_path() /
    "premtables" /
    "prem_44layers.txt";
    std::cout << "Using PREM tables from: " << PATH_premtables << std::endl;

    // Setting up the perfect detector location and exposure
    // Detector location in the Earth model:
    //  Location is specified by the detector radius from the Earth's center, R_det , and
    //  the cartesian coordinates, Location_det , of the detector location in the Earth model.

    // Detector radius from the Earth's center
    double R_det = R_EARTH;            // detector located somewhere in the earth's continental crust

    // Detector cartesian coordinates in the Earth model
    double Location_det[3] = {0.0, 0.0, -1.0 * R_det}; // detector located at the model's South Pole

    double M_det = 10.0;       // detector mass in megaton units.
    double T_det = 10.0;       // detector time of data taking in years.
    double Exposure_det = (M_det *T_det); //*MTonYear ; // detector exposure in [Mton*years] units.


    // To construct the pertubed Earth model, which containts an m-LLVP volume representing a LLVP,
    // we need to specify the m-LLVP parameters, including its height,
    // radius, density contrast, and shape. Additionally, we need to set up the simulation parameters,
    // such as the number of bins for zenith angle, azimuthal angle, and energy.


    std::string geom_mllvp = "pancake";   // m-LLVP geometry: pancake (single shell) or cake (three-shell, narrowing toward top)
    double h_mllvp = 1000;                  // m-LLVP height in km
    double R_mllvp = R_CMB + h_mllvp;

    double R_min = R_CMB - 2500;
    double R_max = R_mllvp + 500;

    // Binning of the angular variables defining the neutrino direction

    // Zenith angle [degrees]
    double th_min = 180 - TMath::ASin(R_max / R_det) * (180.0 / TMath::Pi());
    double th_max = 180 - TMath::ASin(R_min / R_det) * (180.0 / TMath::Pi());
    // Cosine of zenith angle
    double costh_min = cos(th_max * TMath::Pi() / 180.0);
    double costh_max = cos(th_min * TMath::Pi() / 180.0);

    int nbins_th = 10;                      // number of zenith bins [90, 180] degrees

    // Azimuthal angle [degrees]
    double phi_min = 40;
    double phi_max = 150;
    int nbins_phi = 10;                     // number of azimuth bins [0, 360] degrees
    double dphi = (phi_max - phi_min) / (2.0 * nbins_phi);
    double phi_lowedge = phi_min - dphi; // Lower edge of the first azimuthal bin
    double phi_upedge = phi_max + dphi;  // Upper edge of the last azimuthal bin


    // Binning fo the neutrino energy range in GeV
    double enu_min = 1.0;
    double enu_max = 10.0;
    int nbins_enu = 10;                      // number of neutrino energy bins

    // Neutrino flavor labels for file naming
    std::string flvLabels[3] = {"nue", "numu", "nutau"};

    // rhopct scan range
    double rhopct_mllvp_min = -4.0;
    double rhopct_mllvp_max =  4.0;
    double rhopct_mllvp_step = 1.0;
    double rhopct_range[9] = {-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
    int n_rhopct = 9;

    // Create parent output directory
    std::time_t simTime = std::time(nullptr);
    std::tm localTime = *std::localtime(&simTime);
    std::ostringstream timestampStream;
    timestampStream << std::put_time(&localTime, "%Y-%m-%d_%H%M%S");
    fs::path outputDir = fs::current_path() / ("IntSenv_" + timestampStream.str());
    fs::create_directories(outputDir);

    // Initialize Simulation
    NintSimulation Nint_sim;
    Nint_sim.SetDetPosition(Location_det);
    Nint_sim.SetTrueIntervals(enu_min, enu_max, th_min, th_max, phi_lowedge, phi_upedge);
    Nint_sim.SetTrueBinning(nbins_enu, nbins_th, nbins_phi);
    Nint_sim.SetDetExposure(M_det * T_det);
    Nint_sim.SetOutputDir(outputDir);

    // Configure anomaly parameters once (applied only when MantleAnomaly = true)
    Nint_sim.PileHeight = h_mllvp;
    Nint_sim.aperture = 45;
    Nint_sim.AnomalyShape = geom_mllvp;

    for (int flv_nu : {0, 1, 2})
    {
        std::string flvLabel = flvLabels[flv_nu];
        Nint_sim.SetFlavor(flv_nu);

        // Generate reference event distribution for the standard Earth model
        Nint_sim.MantleAnomaly = false;
        std::vector<TH2D*> Nint_reference = Nint_sim.GetNint_ijk();

        // Open sensitivity output file for this flavor
        std::ofstream SenvData(outputDir / ("Senv_" + flvLabel + ".csv"));

        for (int ir = 0; ir < n_rhopct; ++ir)
        {
            double rhopct = rhopct_range[ir];

            // Generate events for the perturbed Earth model at this rhopct
            Nint_sim.MantleAnomaly = true;
            Nint_sim.Pilerho_pct = rhopct;
            std::vector<TH2D*> Nint_perturbed = Nint_sim.GetNint_ijk();

            // Accumulate chi2 over all phi bins
            double chi2_flv = 0.0;
            for (int nphi = 0; nphi < static_cast<int>(Nint_reference.size()); ++nphi)
            {
                chi2_flv += Get2DChi2(Nint_reference[nphi], Nint_perturbed[nphi]);
                delete Nint_perturbed[nphi];
            }

            SenvData << rhopct << " , " << chi2_flv << "\n";
        }

        SenvData.close();
        for (TH2D* h : Nint_reference) delete h;
    }

    std::cout << "Simulation completed and sensitivity generated." << std::endl;

    // Record simulation inputs to the already-created output directory
    RecordIntSenvSimInputs(
        PATH_fluxtables, PATH_premtables,
        geom_mllvp,
        M_det, T_det,
        R_det, Location_det,
        nbins_th, nbins_phi, nbins_enu,
        enu_min, enu_max,
        th_min, th_max,
        phi_lowedge, phi_upedge,
        h_mllvp, R_mllvp,
        rhopct_mllvp_min, rhopct_mllvp_max, rhopct_mllvp_step,
        outputDir
    );

    return 0;
}


void RecordIntSenvSimInputs(
    const fs::path& PATH_fluxtables,
    const fs::path& PATH_premtables,
    const std::string& geom_mllvp,
    double M_det,
    double T_det,
    double R_det,
    const double Location_det[3],
    int nbins_th,
    int nbins_phi,
    int nbins_enu,
    double enu_min,
    double enu_max,
    double th_min,
    double th_max,
    double phi_lowedge,
    double phi_upedge,
    double h_mllvp,
    double R_mllvp,
    double rhopct_mllvp_min,
    double rhopct_mllvp_max,
    double rhopct_mllvp_step,
    fs::path& outputDir
) {
    // outputDir is already created by the caller; write the inputs file into it
    fs::path inputPath = outputDir / "IntN_sim_input.txt";

    std::ofstream inputFile(inputPath);

    if (!inputFile) {
        throw std::runtime_error(
            "Could not create data file: " + inputPath.string()
        );
    }

    int label_expdet = static_cast<int>(M_det * T_det);

    inputFile << "\n--- Interacting events (Perfect Detector) Simulation ---\n";

    inputFile << "Detector: R_det = " << R_det
            << " | Location = ("
            << Location_det[0] << ", "
            << Location_det[1] << ", "
            << Location_det[2] << ")\n";

    inputFile << "Detector: M_det = " << M_det
            << " | T_det = " << T_det << '\n';
    inputFile << "Exposure: " << label_expdet << " Mton*years\n";


    inputFile << "\n--- Simulation Parameters ---\n";


    inputFile << "Flux tables: "
              << fs::absolute(PATH_fluxtables).lexically_normal()
              << '\n';
    inputFile << "PREM tables: "
              << fs::absolute(PATH_premtables).lexically_normal()
              << '\n';

    inputFile << "Energy: enu_min = " << enu_min
            << " | enu_max = " << enu_max
            << " | nbins_enu = " << nbins_enu << '\n';

    inputFile << "Angle: th_min = " << th_min
            << " | th_max = " << th_max
            << " | nbins_th = " << nbins_th << '\n';

    inputFile << "Angle: phi_lowedge = " << phi_lowedge
        << " | phi_upedge = " << phi_upedge
        << " | nbins_phi = " << nbins_phi << '\n';

    inputFile << "\n--- MLLVP Geometry ---\n";
    inputFile << "Dimensions: geom_mllvp = " << geom_mllvp
            << " | h_mllvp = " << h_mllvp
            << " | R_mllvp = " << R_mllvp << '\n';

    inputFile << "Density scan: rhopct_min = " << rhopct_mllvp_min
            << " | rhopct_max = " << rhopct_mllvp_max
            << " | rhopct_step = " << rhopct_mllvp_step << '\n';
}
