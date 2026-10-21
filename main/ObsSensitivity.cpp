/*

Author: Yael Deniz
Phd Candidate, University of Idaho

Description: Code that calculates chi2 for a reconstructed signal of a neutrino
detector

Notes:

3D code has a Break memory violation if theta is 90 \pm 0.1 -> Fixed with a patch
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
#include "GetObservedEvents.h"
#include "ExportData.h"



# define CMB 3480.0 // Core-Mantle-Boundary Radius (km)
# define OCEAN_BOTTOM 6368.0 // Earth radius up to oceanic crust in (km)
# define CONTIENTAL_CRUST 6371.0 // Earth radius in (km)
# define ATM_RADIUS 6386.0 // Radius of the Atmosphere (km)
# define NUCLEON_MASS   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MEGATON_KG  1E9  //Metric MegaTon in kg
# define YEAR_SEC 3.154E7 // Years in Seconds


using namespace std;

void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename);

int main(int argc, char **argv)
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

    std::strftime(datetime_buffer, sizeof(datetime_buffer), 
                                    "%Y-%m-%d_%H%M%S", local_time);
    std::strftime(date_buffer,    sizeof(date_buffer),    
                                    "%Y-%m-%d",        local_time);

    // -----------------------------------------------------------------------------
    // 2. Detector parameters and exposure
    // -----------------------------------------------------------------------------

    double mton_unit = 10.0;        // Detector mass [Mton]
    double year_unit = 20.0;        // Exposure time [years]
    double exposure_unit = mton_unit*year_unit; // Exposure [Mton*years]

   // Total exposure in [kg * s]: (Mton -> kg) × (years -> seconds)
    double exposure_unit_si = exposure_unit*(MEGATON_KG* YEAR_SEC); 

    // Detector position (Cartesian coordinates, units: Earth's radius or km depending on convention)
    double radius_det = CONTIENTAL_CRUST;
    double coordinates[3] = { 0.0, 0.0, -1.0 * radius_det };

    // -----------------------------------------------------------------------------
    // 3. Input file paths (flux and Earth model)
    // -----------------------------------------------------------------------------
    std::string flux_file = BASE_DIR + "/NuFlux/SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string prem_file = OSCPROB_DIR + "/PremTables/prem_46layers_MTZ.txt";
    //std::string prem_file = OSCPROB_DIR + "/PremTables/prem_44layers.txt";

    // -----------------------------------------------------------------------------
    // 4. LAYER/ LLVP / depth configuration
    //    (radius_min / radius_max: radial distances from Earth's center)
    // -----------------------------------------------------------------------------

     //Pile Set Up---------------------------------------------------------------
    //double height_llvp = 1000; // Height of LLVP in km
    //double llvp_radius = CMB + height_llvp; //km
    
    double depth_max  =  6000;
    double depth_min =   300;
    double radius_min = 6371 - depth_max;
    double radius_max = 6371 - depth_min;

    
    //double aperture = 45;
    //double llvp_density_contrasts = 2.0; // 2% density contrats for LLVP
    std::string llvp_morphology = "cake";


    // -----------------------------------------------------------------------------
    // 5. Simulation configuration: binning and ranges
    // -----------------------------------------------------------------------------
    
    // Int binning
    int energy_bins  = 120;  // # 80Bins in energy
    int zenith_bins  = 80; // # 160Bins in zenith/cos(zenith)
    int azimuth_bins = 10; // # 100Bins in azimuth (optimal bins are 110 or 22)
    

    // Reco binning
    int energy_bins_reco  =25; // # 11Bins in energy
    int zenith_bins_reco  =20; // # 10Bins in zenith/cos(zenith)
    int azimuth_bins_reco =5; // # 10Bins in azimuth (optimal bins are 110 or 22)


    // Energy [GeV]    
    double energy_min = 2.0; 
    double energy_max = 10.0;

    //Reconstructed Directions

    // Zenith angles [degrees], derived from depth and detector radius
    double zenith_min = 180.0 - TMath::ASin(radius_max / radius_det) * (180.0 / TMath::Pi());
    double zenith_max = 180.0 - TMath::ASin(radius_min / radius_det) * (180.0 / TMath::Pi());

    double cos_zenith_min = std::cos(zenith_max * TMath::Pi() / 180.0);
    double cos_zenith_max = std::cos(zenith_min * TMath::Pi() / 180.0);

    // Azimuth angle [degrees]
    //double azimuth_min =  90-50;
    //double azimuth_max =  90+50;

     // Azimuthal angle [degrees]
    double azimuth_min = 0.0;
    double azimuth_max = 359.0;

    // -----------------------------------------------------------------------------
    // 6. Output directory and file naming
    // -----------------------------------------------------------------------------
    std::string experiment_test = "MTZ-ne";  // "MTZ-H", "LLVP-ne"
    std::string bin_variables  = "E-th-phi";


    std::string bin_labels = bin_variables+"_Obs-"+std::to_string(zenith_bins_reco)+"-"+ std::to_string(azimuth_bins_reco)+"-"+ std::to_string(energy_bins_reco)+"_Int-"
                            +std::to_string(zenith_bins)+"-"+ std::to_string(azimuth_bins)+"-"+ std::to_string(energy_bins);

    std::string RESULTS_DIR = BASE_DIR + "/PredictedSensitivity/SingleDetector/";

    std::string SENSITIVITY_DIR = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/PredictedSensitivity/SingleDetector/neSenv/";
    std::string EXPERIMENT_DIR = 
                SENSITIVITY_DIR+"ObsSenv_"+experiment_test+"_"+std::string(date_buffer)+"/";

    std::cout << EXPERIMENT_DIR << std::endl;

    mkdir(EXPERIMENT_DIR.c_str(), 0777);

    if (mkdir(EXPERIMENT_DIR.c_str(), 0777) != 0) {
    if (errno != EEXIST) { // ignore "already exists"
        std::cerr << "Failed to create directory " << EXPERIMENT_DIR
                  << ": " << std::strerror(errno) << "\n";
        std::exit(EXIT_FAILURE);  
        }
    }
 
    for (int flavor = 0; flavor < 2; ++flavor) // Neutrino flavor (0: nue, 1: numu, 2: nutau)
    {
    
        /*
         // File and folder structure configuration
        int exposureYearsLabel = static_cast<int>(mton_unit * year_unit);

        //std::string bin_labels = "Obs"+std::to_string(zenith_bins_reco)+"Zen"+ std::to_string(azimuth_bins_reco)+"Az"+ std::to_string(energy_bins_reco)+"Enu_"
                            //+"Int"+std::to_string(zenith_bins)+"Zen"+ std::to_string(azimuth_bins)+"Az"+ std::to_string(energy_bins)+"Enu";





        // Folder nomeclature Obs_[ "LLVP shape"][height_llvp]_[Exposure in Mton]_[Emin-Emax]_[thmin]-[thmax]_[phimin]-[phimax] 
        std::string folderName = "ObsChi2_"+testVariable+llvp_morphology+std::to_string(static_cast<int>(height_llvp))+"_"+std::to_string(exposureYearsLabel)+"_"+ std::to_string(static_cast<int>(energy_min)) +"-"+ 
                                     std::to_string(static_cast<int>(energy_max))+"_" + std::to_string(static_cast<int>(zenith_min)) +"-"+ 
                                     std::to_string(static_cast<int>(zenith_max))+"_"  + std::to_string(static_cast<int>(azimuth_min)) +"-"+ 
                                     std::to_string(static_cast<int>(azimuth_max))+"_"+std::string(dateBuffer)+"/";
        
        std::string pathToResults = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/PreliminaryResults/";
        std::string pathToChiSquares = pathToResults + "ObsChi2/";
        std::string chiSquareFolderPath = pathToChiSquares + folderName;

        std::string chiSquareFileName = chiSquareFolderPath + "ObsChi2_cth_"+testVariable+llvp_morphology+std::to_string(static_cast<int>(height_llvp))+"_"+
                                        std::to_string(exposureYearsLabel)+"nu" + std::to_string(flavor) + "_" + bin_labels + ".csv";

        // Create output directories if they do not exist
        mkdir(chiSquareFolderPath.c_str(), 0777);
        */

        // Output CSV file for this flavor
        std::string exp_senv_file =
            EXPERIMENT_DIR
            + "RecoSenv_flv_" + std::to_string(flavor)
            + "_" + bin_labels + ".csv";
        // Open CSV output for sensitivity scan
        std::ofstream senv_data(exp_senv_file);

        // File with run details (same for all flavors; overwritten or updated each time)
        std::string details_file =
            EXPERIMENT_DIR + "ObsSenv_details_" + std::string(datetime_buffer) + ".txt";

        // -------------------------------------------------------------------------
        // 7.1 Initialize Asimov simulation with standard Earth model
        // -------------------------------------------------------------------------
        AsimovObsSimulation experiment;

       experiment.ThePremTable = prem_file;
       experiment.TheHondaTable = flux_file;
       experiment.SetDetectorXYZ(coordinates);
       experiment.SetIntervals(zenith_min,zenith_max,azimuth_min,azimuth_max,energy_min,energy_max);
       experiment.SetTrueBinning(energy_bins,zenith_bins,azimuth_bins);
       experiment.SetRecoBinning(energy_bins_reco,zenith_bins_reco,azimuth_bins_reco);
       experiment.SetExposure(exposure_unit_si);
       experiment.flvf=flavor;

        // -----------------------------------------------------------------------------
        // Configure the "expected" (reference) Earth model
        // This setup represents the standard Earth without any LLVP / mantle anomaly.
        // The resulting event distribution will be used as the reference ("Asimov") data.
        // -----------------------------------------------------------------------------

       experiment.PileInModel = false;
       experiment.ShapeOfPile = llvp_morphology;
       experiment.ThePileHight = 0.0;
       experiment.ThePileAperture = 0.0;
       experiment.ThePileDensityContrast = 0.0;
       experiment.ThePileChemicalContrast = 0.0;


       //std::vector<TH2D*>  RecoEventsStandard = experiment.GetObsEvents3Dthvb();

       //TH2D*  ObservedRatesStandard = experiment.GetRecoRate_vMFMarginal();
       //TH2D*  RecoEventsStandard = experiment.GetNreco_vMFMarginal(ObservedRatesStandard);

       TH3D * N_reco_expected = experiment.GetRecoN();

        // -------------------------------------------------------------------------
        // 7.2 Parameter scan over density contrast rho
        // -------------------------------------------------------------------------
        
        constexpr int kNumRho = 9;

        double rho[kNumRho] = { -2.5, -2.0, -1.5,
                                -1.0, -0.5, 0.0,
                                 0.5,  1.5,  2.0 };
        
        double  h[kNumRho] = {100,200,300,
                              400,500,600,
                              700,800,900};

        double azi[kNumRho] ={15,15,15,
                              15,15,25,
                              35,45,50}; 

        for (int i = 0; i < kNumRho ; ++i)
        {
            // -----------------------------------------------------------------------------
            // Configure "observed" Earth model for this rho[i] step
            // -----------------------------------------------------------------------------

            // Reset to the baseline PREM table before applying modifications

            //experiment.ThePremTable = prem_file;
            //experiment.PileInModel = false;
            //experiment.ShapeOfPile = llvp_morphology;
            //experiment.ThePileDensityContrast = rho[i];
            experiment.ChangeLayerProperties(37, rho[i], 0);
            experiment.ChangeLayerProperties(38, rho[i], 0);
            experiment.ChangeLayerProperties(39, rho[i], 0);
            
            // Apply a density modification to layer 34:
            //   - 34      : index of the PREM layer to modify
            //   - rho[i]  : density contrast (e.g. in % or relative units, depending on implementation)
            //   - 0       : additional flag / option (here: no extra behavior)
            // experiment.ModifyLayer(34, rho[i], 0);

            // If instead of a single-layer modification you want a full mantle anomaly / LLVP model,
            // you would enable and configure it like this:

            //experiment.ThePileHight = h[i];
            //experiment.ThePileAperture = azi[i];


            //std::vector<TH2D*>  RecoEventsAlternative = experiment.GetObsEvents3Dthvb();
            //TH2D*  ObservedRatesAlternative = experiment.GetRecoRate_vMFMarginal();
            //TH2D*  RecoEventsAlternative = experiment.GetNreco_vMFMarginal(ObservedRatesAlternative);

            // -----------------------------------------------------------------------------
            // Generate events for the modified (anomalous) Earth model
            // -----------------------------------------------------------------------------

            TH3D * N_reco_observed = experiment.GetRecoN();

            //double chi2_tot  = IntChi2( N_reco_expected , N_reco_observed);
             double chi2_tot  = IntChi2(N_reco_expected, N_reco_observed);

            //for (int n = 0; n < RecoEventsAlternative.size(); ++n)
            //{
            //chi2tot += Get2DChi2( N_reco_expected[n] , N_reco_observed[n]);

            //chi2tot+= Get2DChi2( RecoEventsStandard[n] , RecoEventsAlternative[n]);
                //testchi2 = chi2tot;
            //}
            // Write: index, aperture, density contrast, pile height, chi^2
            senv_data << i << " , "
                      << rho[i] << " , "
                      << chi2_tot << " , "
                      << zenith_bins_reco << " , "
                      << azimuth_bins_reco << " , " 
                      << energy_bins_reco <<  '\n';

            
        }

        senv_data.close();

        // -------------------------------------------------------------------------
        // 7.3 Write run details (geometry, configuration, etc.)
        // -------------------------------------------------------------------------
        experiment.GetDetails(details_file);

    }   // end loop over flavors

    // -----------------------------------------------------------------------------
    // 8. End timing and print execution time
    // -----------------------------------------------------------------------------
    auto end = Clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Execution time: " << duration.count() << " seconds\n";
    std::cout << "rates test" << std::endl;
   
    return 0;
   


}


