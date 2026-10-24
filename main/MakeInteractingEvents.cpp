
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
#include "StatsAna.h"
#include "GetTrueEvents.h"



// Boundaries of major layers of the Earth
# define R_CMB 3480.0 // Core-Mantle-Boundary Radius (km)
# define R_OCEAN 6368.0 // Earth radius in (km)
# define R_EARTH 6371.0 // Earth radius in (km)
# define R_ATM 6386.0 // Radius of the Atmosphere (km)

// Some relevant units
# define M_n   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTonYear  3.154E16  //Metric MegaTon-Year unit in kg-seconds
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds

using namespace std;
namespace fs = std::filesystem;

void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename);
void RecordIntSimInputs(
    const fs::path& PATH_fluxtables,
    const fs::path& PATH_premtables,
    int flv_nu,
    const std::string& geom_mllvp,
    double M_det, double T_det,
    double R_det, const double Location_det[3],
    int nbins_th, int nbins_phi, int nbins_enu,
    double enu_min, double enu_max,
    double th_min, double th_max,
    double phi_lowedge, double phi_upedge,
    double h_mllvp, double R_mllvp, double rhopct_mllvp,
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
    fs::current_path().parent_path() /
    "fluxtables" /
    "SP_AziAveraged_solmin" / // Differential flux tables for the South Pole (SP), averaged over azimuthal angles(AziAveraged), and for solar minimum conditions (solmin).
    "spl-nu-20-01-000.d";

    std::cout << "Using flux tables from: " << PATH_fluxtables << std::endl;
    fs::path PATH_premtables =
    fs::current_path().parent_path() /
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
    double Exposure_det = (M_det *T_det)*MTonYear ; // detector exposure in [Mton*years] units.


    // To construct the pertubed Earth model, which containts an m-LLVP volume representing a LLVP,
    // we need to specify the m-LLVP parameters, including its height,
    // radius, density contrast, and shape. Additionally, we need to set up the simulation parameters,
    // such as the number of bins for zenith angle, azimuthal angle, and energy.


    std::string geom_mllvp = "pancake";   // m-LLVP geometry: pancake (single shell) or cake (three-shell, narrowing toward top)
    double h_mllvp = 1000;                  // m-LLVP height in km
    double R_mllvp = R_CMB + h_mllvp;

    double R_min = R_CMB - 2500;
    double R_max = R_mllvp + 500;

    double rhopct_mllvp = 2.0;          // m-LLVP percentage density anomaly


    // Binning of the angular variables defining the neutrino direction

    // Zenith angle [degrees]
    double th_min = 180 - TMath::ASin(R_max / R_det) * (180.0 / TMath::Pi());
    double th_max = 180 - TMath::ASin(R_min / R_det) * (180.0 / TMath::Pi());
    // Cosine of zenith angle
    double costh_min = cos(th_max * TMath::Pi() / 180.0);
    double costh_max = cos(th_min * TMath::Pi() / 180.0);

    int nbins_th = 100;                      // number of zenith bins [90, 180] degrees

    // Azimuthal angle [degrees]
    double phi_min = -50;
    double phi_max = 50;
    int nbins_phi = 100;                     // number of azimuth bins [0, 360] degrees
    double dphi = (phi_max - phi_min) / (2.0 * nbins_phi);
    double phi_lowedge = phi_min - dphi; // Lower edge of the first azimuthal bin
    double phi_upedge = phi_max + dphi;  // Upper edge of the last azimuthal bin


    // Binning fo the neutrino energy range in GeV
    double enu_min = 1.0;
    double enu_max = 10.0;
    int nbins_enu = 100;                      // number of neutrino energy bins

    // flavor of the interacting neutrino  (0: nue, 1: numu, 2: nutau)
    int flv_nu = 0;

    // Initialize Simulation
    AsimovSimulation Nint_sim;
    Nint_sim.PremTable = PATH_premtables;
    Nint_sim.HondaTable = PATH_fluxtables;
    Nint_sim.SetDetectorPosition(Location_det);
    // Set intervals, binning, and exposure for the simulation
    Nint_sim.SetIntervals(th_min, th_max, phi_lowedge, phi_upedge, enu_min, enu_max);
    Nint_sim.SetBinning(nbins_th, nbins_phi, nbins_enu);
    Nint_sim.SetExposure(Exposure_det);



    Nint_sim.MantleAnomaly = false; // Standard Earth model without m-LLVP

    // Run simulation and generate events for Standard and Alternative
    Nint_sim.flvf = flv_nu;
    std::vector<TH2D*> Nint_reference = Nint_sim.GetTrueEvents3D();

    // Configure for Anomalous Earth Model (with LLVP)
    Nint_sim.MantleAnomaly = true;
    Nint_sim.PileHeight = h_mllvp;
    Nint_sim.aperture = 45;
    Nint_sim.AnomalyShape = geom_mllvp;
    Nint_sim.PileDensityContrast = rhopct_mllvp;

    // Generate events for the perturbed  Earth model
    Nint_sim.flvf = flv_nu;
    std::vector<TH2D*> Nint_perturbed = Nint_sim.GetTrueEvents3D();

    // Create output folder and record simulation inputs
    fs::path outputDir;
    RecordIntSimInputs(
        PATH_fluxtables, PATH_premtables,
        flv_nu, geom_mllvp,
        M_det, T_det,
        R_det, Location_det,
        nbins_th, nbins_phi, nbins_enu,
        enu_min, enu_max,
        th_min, th_max,
        phi_lowedge, phi_upedge,
        h_mllvp, R_mllvp, rhopct_mllvp,
        outputDir
    );

    // Neutrino flavor label for file naming
    std::string flvLabels[3] = {"nue", "numu", "nutau"};
    std::string flvLabel = flvLabels[flv_nu];

    // Exporting Events
    for (int nhist = 0; nhist < static_cast<int>(Nint_perturbed.size()); ++nhist)
    {
        fs::path csvPath = outputDir / ("Nint_" + flvLabel + "_" + std::to_string(nhist) + ".csv");

        TH2D* dNintpct_nphi = new TH2D(Form("OscDiff2D%d", nhist), Form("Oscillogram%d", nhist), nbins_th, costh_min, costh_max, nbins_enu, enu_min, enu_max);

        // Compute the percentage differences between distributions and export result to a CSV
        GetDiff2D(Nint_reference[nhist], Nint_perturbed[nhist], dNintpct_nphi);
        ExportToCSV(Nint_reference[nhist], Nint_perturbed[nhist], dNintpct_nphi, csvPath.string());
    }

    // Simulation completed
    std::cout << "Simulation completed and events generated." << std::endl;

    return 0;
}


void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename)
{
    std::ofstream outfile(filename);

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



void RecordIntSimInputs(
    const fs::path& PATH_fluxtables,
    const fs::path& PATH_premtables,
    int flv_nu,
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
    double rhopct_mllvp,
    fs::path& outputDir
) {
    std::string flvLabels[3] = {"nue", "numu", "nutau"};
    std::string flvLabel = flvLabels[flv_nu];

    // Generate timestamp
    std::time_t now = std::time(nullptr);
    std::tm localTime = *std::localtime(&now);

    std::ostringstream timestampStream;
    timestampStream << std::put_time(&localTime, "%Y-%m-%d_%H%M%S");

    // Create output directory in the project root (parent of main/)
    outputDir = fs::current_path().parent_path() /
        ("Nint_" + flvLabel + "_" + timestampStream.str());

    fs::create_directories(outputDir);

    fs::path inputPath = outputDir / "IntN_sim_input.txt";

    // Write simulation metadata
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

    inputFile << "Interacting flavor: " << flvLabel << '\n';

    inputFile << "\n--- MLLVP Geometry ---\n";
    inputFile << "Dimensions: geom_mllvp = " << geom_mllvp
            << " | h_mllvp = " << h_mllvp
            << " | R_mllvp = " << R_mllvp << '\n';

    inputFile << "Density: rhopct_mllvp = " << rhopct_mllvp << '\n';

}
