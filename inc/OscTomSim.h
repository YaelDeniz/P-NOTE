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



using namespace std;
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)

//TH2D*  AsimovTrueEvents(std::string modelname,bool LLVP,  std::vector<int> layers, int flvf, double Region[], int Bins[], double NnT);

class OscTom
{
    public:

    //          Neutrino properties
    
    double bar;
    double flvi;
    double flvf;
    double energy;
    double zenith;
    double azimuth;

    //          Detector properties
    double pos[3] = {0.0, 0.0, -1.0*Rearth};
    double mass;
    double time;
    
    // Energy resolution 
    double Ae = 0.05;
    double Be = 0.1;

    // Zenith resolution
    double Ath = 0.04;  //0.04 + 0.18/(sqrt(e))
    double Bth = 0.18;

    // Azimuth resolution


    //Default constructor
    OscTom() : mass(10.0), time(10.0) {}

    //Earth Model

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

    // Neutrino Flux
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string HondaTable = FluxFolder + FluxFile; //Class Assumes Sout Pole flux as default
 

    

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
    int ibins; //Bins in Zenith
    int jbins; //Bins in Azimuth
    int kbins; //Bins in Energy

    void SetDetectorPosition( double x[3])
    {

        pos[0] = x[0];
        pos[1] = x[1];
        pos[2] = x[2];

    }


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
        ibins = i;
        jbins = j+1; // Bins centered at values
        kbins = k;
    }

    void SetExposure(double MassTime) 
    {

        MT = MassTime; //MTon-Years

    }

   // TH3D * GetTrueEvents3D();

    std::vector< TH2D* > GetTrueEvents2D(); //Azimuth independet
    
    std::vector< TH2D* >  GetTrueEvents3D();
    std::vector< TH2D* >  GetTrueEvents3Dth();

   // TH2D * GetOscProb2D(int flvi, int flvf, bool nunubar  ); //Azimuth independet
    
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