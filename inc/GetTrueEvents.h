/*

Author: Yael Deniz
Phd Candidate, University of Idaho

Notes:

All angular inputs are in degrees but must be implemented as radians in codes

Neutrino flux is interpolated in the space log10(dPhi/dEdcth) and (log10E, cos(\theta))
    **Rember to adjust the result of the iinterpolation into dPhi/dEdth (dPhi/dEdcth; dPhi/dlog10Edth)

Bnning scheme Convention:

    ibins - Total number of Energy [GeV] bins
    jbins - Total number of zenith("theta") [rad] bins | cos(theta)
    kbins - Total number of azimuth [rad] bins 

*/

#ifndef TRUEEVENTGENERATOR_H
#define TRUEEVENTGENERATOR_H

//C++
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <math.h> 


//CERN ROOT
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"
#include "TGraph.h"

//Layer Updates
#include "LayerUpdate.h"



using namespace std;
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)

//TH2D*  AsimovTrueEvents(std::string modelname,bool LLVP,  std::vector<int> layers, int flvf, double Region[], int Bins[], double NnT);

class AsimovSimulation
{
    public:
    

    std::string label = "std";

    double MT; //Detector Exposure
  
    
    double pos[3] = {0.0, 0.0, -1.0*Rearth};
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string HondaTable = FluxFolder + FluxFile; //Class Assumes Sout Pole flux as default

    //Earth Settings

    //The Earth model
    std::string PremFolder = "/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_44layers.txt";
    std::string PremTable = PremFolder+PremName;
 
    //The LLVP
    bool MantleAnomaly = false;
    double PileHeight = 1000;
    double aperture = 45.0;
    std::string AnomalyShape ="cake";
    double PileDensityContrast = 2.0;
    double PileChemContrast = 0.0;
    

    //Modify specific Layers
    int PremTableNumber = 44 ;
    double DensityContrast = 0 ;
    double ChemicalContrast = 0 ;

    //Neutrino Settings
    int flvf;
    int nunubar = 1; //Neutrinos

    double Emin = 1.00; //GeV
    double Emax = 10.00; //GeV

    double thmin; // Min 90;
    double thmax; // Max 180;

    double phimin=0; // Min 0;
    double phimax=TMath::Pi(); //Max 360;

    //Simulation Settings
    int ibins; //Bins in Energy
    int jbins; //Bins in Zenith
    int kbins; //Bins in Azimuth

    void SetDetectorPosition( double x[3])
    {

        pos[0] = x[0];
        pos[1] = x[1];
        pos[2] = x[2];

    }

    //PATCH 12/3/2025 -- Multiply layer modification

    std::vector <LayerUpdate> LayerUpdates_intsim; 

    void ChangeLayerProperties(int n , double rhopct , double zoapct);

    void AddUpdatesList(const std::vector<LayerUpdate>& updates);


    void ModifyLayer (int n , double rhopct , double zoapct)
    {

        PremTableNumber = n; // Index of layer
        DensityContrast = rhopct; // Density percentage difference
        ChemicalContrast = zoapct; // zoa(compositional) perctengae differnce

    }

    void SetIntervals(double zLow, double zUpper, double azLow, double azUpper, double eLow, double eUpper)
    {
        thmin = zLow*TMath::Pi()/180;
        thmax = zUpper*TMath::Pi()/180;

        phimin = azLow*TMath::Pi()/180;
        phimax = azUpper*TMath::Pi()/180;
        
        Emin = eLow;
        Emax = eUpper;

    }

    void SetBinning(int i, int j, int k)
    {
        ibins = i;   //Energy bins
        jbins = j;   //Zenith bins  
        kbins = k+1; // Azimuth bins centered at values
    }

    void SetExposure(double MassTime) 
    {

        MT = MassTime; //MTon-Years

    }

    TH3D * GetIntRates(); // Creates interacting rates;  Def in AsimovDataTrue.C, line 102
    TH3D * IntEvents(TH3D* ratesHist);   // Integrates interacting rates to obtain int events;  Def in AsimovDataTrue.C, line 266
    TH3D * GetIntEvents(); // Access interacting event distribution; Def in  
    TH2D* GetAzimuthAvaraged( TH3D* Nint_hist);
    TH3D * GetIntSignal();

    std::vector< TH2D* > GetTrueEvents2D(); //Azimuth independet
    std::vector< TH2D* > GetTrueEventsMarginal(); // Marginalized over azimuth
    
    std::vector< TH2D* >  GetTrueEvents3D();
    std::vector< TH2D* >  GetTrueEvents3Dth();
    std::vector< TH2D* >  GetTrueEvents3Dthlog10E();

   //TH2D * GetOscProb2D(int flvi, int flvf, bool nunubar  ); //Azimuth independet
    
    TGraph * GetOscProb( int flvi, int flvf, bool nunubar, double cth ); //To be Deleted
    
   // TH2D * SensitivityTrueEvents2D( int n, double pct ); //To be Del

    //TH2D * TestTrueEvents2D(std::string model_std , std::string model_alt);

   // std::vector< std::vector<double> > GetPremMatrix( std::string path2table  );


   void GetDetails(const std::string& PATH)
   {

        std::ofstream details(PATH);
    
        if (!details) {
            std::cout << "Error: Unable to open file " << PATH << std::endl;
            return;
        }

        details << "--------------------------------------------------------------------------------\n";
        details << "                                Detector Details                               \n";
        details << "--------------------------------------------------------------------------------\n";
        details << " Detector Radius: " << pos[2] << " m\n";
        details << " Neutrino Flux: " << FluxFile << " \n";
        details << " Number of Detectors: 1\n";
        details << " Mass: " << 10.0 << " Mton\n";
        details << " Time: " << 10.0 << " years\n";
        details << " Exposure: " << MT << " Mton-Year\n";

        details << "\n--------------------------------------------------------------------------------\n";
        details << "                                Earth Details                                  \n";
        details << "--------------------------------------------------------------------------------\n";
        details << " Earth Model: " << PremName << "\n";
        details << "Index of modified layer (default is 0): " << PremTableNumber << "\n";
        details << "Density Contrast (drho/rho): " << DensityContrast<< "\n";
        details << " Chem Contrast (H%): " << ChemicalContrast << "\n";
        
        details << "\n--------------------------------------------------------------------------------\n";
        details << "                        Lower Mantle Anomalies (LLVP)                          \n";
        details << "--------------------------------------------------------------------------------\n";
        details << " LLVP: " << MantleAnomaly << "\n";
        details << " LLVP Model: " << AnomalyShape << "\n";
        details << " LLVP Height: " << PileHeight << " km\n";
        details << " Aperture: " << aperture << "\n";
        details << " Density Contrast (drho/rho): " << PileDensityContrast << "\n";
        details << " Chem Contrast (H%): " << PileChemContrast << "\n";

        details << "\n--------------------------------------------------------------------------------\n";
        details << "                             Simulation Details                                \n";
        details << "--------------------------------------------------------------------------------\n";
        details << " Neutrino Energy Range: [" << Emin << " - " << Emax << "] GeV | Bins: " << kbins << "\n";
        details << " Zenith Range: [" << thmin*(180.0/TMath::Pi()) << " (" << cos(thmin) << ") - " << thmax*(180.0/TMath::Pi()) << " (" << cos(thmax) << ")] | Bins: " << ibins << "\n";
        details << " Azimuth Range: [" << phimin*(180.0/TMath::Pi()) << " - " << phimax*(180.0/TMath::Pi()) << "] Bins: " << jbins<< "\n";

        details.close();

   }



};




TH2D*  OscProbEarth(std::string modelname, int flvf, double Region[], int Bins[], double NnT);

//TObjArray*  GetTrueEventsMC(std::string modelname, int flvf , double Region[], int Bins[], double NnT, int K);
//TH2D*  GetTrueEvents(std::string modelname, int flvf , double Region[], int Bins[], double NnT);
/*
List of Variables 

Region[]={Emin, Emax, Etamin, Etamax, dAz}

Bins[]={Angular bins, Energy bins}

NnT = Gigaton-year exposure

K = Number of Pseudo Experiments

*/


#endif
