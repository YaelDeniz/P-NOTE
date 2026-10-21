
/* Direction of neutrinos a re given in degrees*/
#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <string>
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
#include "SimPaths.h" 
#include "MathTools.h"
#include "PhyTools.h"
#include "PlotTool.h"
#include "StatsAna.h"
#include "GetTrueEvents.h"
#include "ExportData.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need

# define CMB 3480.0 // Core-Mantle-Boundary Radius (km)
# define OCEAN_BOTTOM 6368.0 // Earth radius up to oceanic crust in (km)
# define CONTIENTAL_CRUST 6371.0 // Earth radius in (km)
# define ATM_RADIUS 6386.0 // Radius of the Atmosphere (km)
# define NUCLEON_MASS   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MEGATON_KG  1E9  //Metric MegaTon in kg
# define YEAR_SEC 3.154E7 // Years in Seconds

using namespace std;


int main(int argc, char** argv)
{
    using Clock = std::chrono::high_resolution_clock;

    // -----------------------------------------------------------------------------
    // 0. Start timing
    // -----------------------------------------------------------------------------
    auto start = Clock::now();

    // -----------------------------------------------------------------------------
    // 1. Time stamp and date strings (for file naming)
    // -----------------------------------------------------------------------------
    std::time_t now = std::time(nullptr);
    std::tm* local_time = std::localtime(&now);

    char datetime_buffer[21];  // "YYYY-MM-DD_HHMMSS" + '\0'
    char date_buffer[21];      // "YYYY-MM-DD" + '\0'

    std::strftime(datetime_buffer, sizeof(datetime_buffer), "%Y-%m-%d_%H%M%S", local_time);
    std::strftime(date_buffer,    sizeof(date_buffer),    "%Y-%m-%d",        local_time);

    // -----------------------------------------------------------------------------
    // 2. Detector parameters and exposure
    // -----------------------------------------------------------------------------
    double mton_unit = 10.0;        // Detector mass [Mton]
    double year_unit = 10.0;        // Exposure time [years]
    double exposure_unit = mton_unit*year_unit; // Exposure [Mton*years]

    // Total exposure in [kg * s]: (Mton -> kg) × (years -> seconds)
    double exposure_unit_si = exposure_unit*(MEGATON_KG* YEAR_SEC); //Exporure in SI units (Requiered units for the simulation)

    // Detector position (Cartesian coordinates, units: Earth's radius or km depending on convention)
    double radius_det = CONTIENTAL_CRUST;
    double coordinates[3] = { 0.0, 0.0, -1.0 * radius_det };

    // -----------------------------------------------------------------------------
    // 3. Input file paths (flux and Earth model)
    // -----------------------------------------------------------------------------
    std::string flux_file = BASE_DIR + "/NuFlux/SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string prem_file = OSCPROB_DIR + "/PremTables/prem_default.txt";

    // -----------------------------------------------------------------------------
    // 4. LAYER/ LLVP / depth configuration
    //    (radius_min / radius_max: radial distances from Earth's center)
    // -----------------------------------------------------------------------------
    double height_llvp = 1000;                  // LLVP height in km
    double radius_min = CMB - 500.0;
    double radius_max = (CMB + height_llvp) + 500.0;

    // -----------------------------------------------------------------------------
    // 5. Simulation configuration: binning and ranges
    // -----------------------------------------------------------------------------
    // Binning
    int energy_bins  = 80;
    int zenith_bins  = 40;
    int azimuth_bins = 101;

    // Energy [GeV]
    double energy_min =  1.0;
    double energy_max = 10.0;

    // Zenith angles [degrees], derived from depth and detector radius
    double zenith_min = 180.0 - TMath::ASin(radius_max / radius_det) * (180.0 / TMath::Pi());
    double zenith_max = 180.0 - TMath::ASin(radius_min / radius_det) * (180.0 / TMath::Pi());

    double cos_zenith_min = std::cos(zenith_max * TMath::Pi() / 180.0);
    double cos_zenith_max = std::cos(zenith_min * TMath::Pi() / 180.0);

    (void)cos_zenith_min; // currently unused
    (void)cos_zenith_max; // currently unused

    // Azimuthal angle [degrees]
    double azimuth_min = 0.0;
    double azimuth_max = 359.0;

    // -----------------------------------------------------------------------------
    // 6. Output directory and file naming
    // -----------------------------------------------------------------------------
    std::string experiment_test = "MTZ-ne";  // "MTZ-H", "LLVP-ne"
    std::string bin_variables  = "E-th-phi";

    // Label for binning configuration
    std::string bin_labels =
        bin_variables + "_"
        + std::to_string(zenith_bins) + "-"
        + std::to_string(azimuth_bins) + "-"
        + std::to_string(energy_bins);

    // Base directory for simulation results
    std::string RESULTS_DIR = BASE_DIR + "/PredictedSensitivity/SingleDetector/";

    std::string SENSITIVITY_DIR = RESULTS_DIR + "neSenv/";

    // Directory for this particular run (depends on test variable and date)
    std::string EXPERIMENT_DIR =
        SENSITIVITY_DIR + "IntSenv_" + experiment_test + "_" + std::string(datetime_buffer) + "/";

    // Create output directory (ignore errors, 0777 permissions)
    mkdir(EXPERIMENT_DIR.c_str(), 0777);

    // -----------------------------------------------------------------------------
    // 7. Loop over neutrino flavors
    // -----------------------------------------------------------------------------
    for (int flavor = 0; flavor < 2; ++flavor)
    {
        // Output CSV file for this flavor
        std::string exp_senv_file =
            EXPERIMENT_DIR
            + "IntSenv_flv_" + std::to_string(flavor)
            + "_" + bin_labels + ".csv";
        // Open CSV output for sensitivity scan
        std::ofstream senv_data(exp_senv_file);

        // File with run details (same for all flavors; overwritten or updated each time)
        std::string details_file =
            EXPERIMENT_DIR + "IntSenv_details_" + std::string(datetime_buffer) + ".txt";

        // -------------------------------------------------------------------------
        // 7.1 Initialize Asimov simulation with standard Earth model
        // -------------------------------------------------------------------------
        AsimovSimulation experiment;
        experiment.PremTable   = prem_file;
        experiment.HondaTable  = flux_file;
        experiment.SetDetectorPosition(coordinates);
        experiment.SetIntervals(zenith_min, zenith_max,
                                azimuth_min, azimuth_max,
                                energy_min,    energy_max);
        experiment.SetBinning(energy_bins, zenith_bins, azimuth_bins);
        experiment.SetExposure(exposure_unit_si);
        experiment.flvf = flavor;

        // -----------------------------------------------------------------------------
        // Configure the "expected" (reference) Earth model
        // This setup represents the standard Earth without any LLVP / mantle anomaly.
        // The resulting event distribution will be used as the reference ("Asimov") data.
        // -----------------------------------------------------------------------------

        experiment.MantleAnomaly      = false;      // Use the unmodified mantle density profile (no anomaly)
        experiment.PileHeight         = height_llvp;     // LLVP vertical extent [km] measured from the CMB
        experiment.aperture           = 45;         // Half-opening angle of the LLVP region [degrees]
        experiment.AnomalyShape = "cake";   // Parameterization of the LLVP geometry (e.g. "cake")
        experiment.PileDensityContrast = 2;         // Density contrast [%] of the LLVP relative to background mantle


        // Expected (reference) events for the unmodified Earth model
        TH3D* Nint_expected = experiment.GetIntEvents();

        // -------------------------------------------------------------------------
        // 7.2 Parameter scan over density contrast rho
        // -------------------------------------------------------------------------
        constexpr int kNumRho = 9;
        double rho[kNumRho] = { -4.0, -3.0, -2.0, -1.0,
                                 0.0,
                                 1.0,  2.0,  3.0,  4.0 };

        for (int i = 0; i < kNumRho; ++i)
        {
            // -----------------------------------------------------------------------------
            // Configure "observed" Earth model for this rho[i] step
            // -----------------------------------------------------------------------------

            // Reset to the baseline PREM table before applying modifications
            experiment.PremTable = prem_file;

            // Apply a density modification to layer 34:
            //   - 34      : index of the PREM layer to modify
            //   - rho[i]  : density contrast (e.g. in % or relative units, depending on implementation)
            //   - 0       : additional flag / option (here: no extra behavior)
            // experiment.ModifyLayer(34, rho[i], 0);

            // If instead of a single-layer modification you want a full mantle anomaly / LLVP model,
            // you would enable and configure it like this:
            //
            experiment.MantleAnomaly        = true;     // Turn on parameterized mantle anomaly
            experiment.PileDensityContrast  = rho[i];   // Density contrast of the anomaly
            //experiment.PileHeight           = h[i];     // LLVP height [km]
            //experiment.aperture             = azi[i];   // Half-opening angle [degrees]
            //experiment.ChangeLayerProperties()

            // -----------------------------------------------------------------------------
            // Generate events for the modified (anomalous) Earth model
            // -----------------------------------------------------------------------------

            // 3D binned interaction rates: N(energy, zenith, azimuth)
            TH3D* Nint_observed = experiment.GetIntEvents();

            // Compute chi^2 between expected and observed event distributions
            double chi2_tot = IntChi2(Nint_expected, Nint_observed);

            // Write: index, aperture, density contrast, pile height, chi^2
            senv_data << i << " , "
                      << experiment.aperture << " , "
                      << experiment.PileDensityContrast << " , "
                      << experiment.PileHeight << " , "
                      << chi2_tot << '\n';
        }

        senv_data.close();

        // -------------------------------------------------------------------------
        // 7.3 Write run details (geometry, configuration, etc.)
        // -------------------------------------------------------------------------
        experiment.GetDetails(details_file);
    } // end loop over flavors

    // -----------------------------------------------------------------------------
    // 8. End timing and print execution time
    // -----------------------------------------------------------------------------
    auto end = Clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Execution time: " << duration.count() << " seconds\n";
    std::cout << "rates test" << std::endl;

    return 0;
}







