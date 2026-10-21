#ifndef NUASIMOVSIMULATION_H
#define NUASIMOVSIMULATION_H

//C++
#include <string.h>
#include <stdlib.h>
#include <vector>

//CERN ROOT
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"
#include "TGraph.h"

# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)

class NuAsimovSimulation
{
    public:

    //Detector Position
    double xyzTelescope[3] = {0.0, 0.0, -1.0*Rearth};

    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string TheHondaTable = FluxFolder + FluxFile; //Class Assumes Sout Pole flux as default

    //The Earth model
    std::string PremFolder = "/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_44layers.txt";
    std::string ThePremTable = PremFolder+PremName;

    //Modify specific Layers
    int PremLayer = 44 ;
    double LayerDensityContrast = 0;
    double LayerChemicalContrats = 0;

    //LLVP

    bool PileInModel = false; 
    std::string ShapeOfPile = "cake";
    double ThePileHight = 1000.0;
    double ThePileAperture=45.0;
    double ThePileDensityContrast = 10.0;
    double ThePileChemicalContrast = 10.0;

    //Neutrino Settings
    int flvf;
    
    double Emin = 1.00; //GeV
    double Emax = 10.00; //GeV

    double Zenmin; // Min 90;
    double Zenmax; // Max 180;

    double Azimin=0; // Min 0;
    double Azimax=360; //Max 360;

    //Simulation Settings

    //Original Bins
    int thbins0; //Bins in Zenith
    int phibins0; //Bins in Azimuth
    int ebins0; //Bins in Energy

    //Reconstructed Bins
    int thbins; //Bins in Zenith
    int phbins; //Bins in Azimuth
    int ebins; //Bins in Energy


    double NnT; //Detector Exposure //0.05*e +0.1*sqrt(e)
    

    void ChangeLayer (int n , double rhopct , double zoapct)
    {

        PremLayer = n; // Index of layer in .txt file 
        LayerDensityContrast = rhopct; // Density percentage difference
        LayerChemicalContrats= zoapct; // zoa(compositional) perctengae differnce

    }

    void SetIntervals(double th1, double th2, double azi1, double azi2, double e1, double e2)
    {
        Zenmin = z1;
        Zenmax = z2;

        Azimin = a1;
        Azimax = a2;
        
        Emin = e1;
        Emax = e2;

    }

    void SetTrueBinning(int ibins, int jbins, int kbins)
    {
        thbins0 = ibins;
        phibins0 = jbins;
        ebins0   = kbins;
    }

    void SetRecoBinning(int mbins, int nbins, int lbins) 
    {
        thbins = mbins;
        phibins = nbins;
        ebins  = lbins;
    }


    void SetExposure(double exposure) 
    {

        NnT = exposure;

    }

    void SetDetectorXYZ(double xyz[3])
    {
        xyzTelescope[0]=xyz[0];
        xyzTelescope[1]=xyz[1];
        xyzTelescope[2]=xyz[2];

    }

    TH2D * IntEventGen2Dcth();

    TH2D * RecoEvents2Dcth();

};



#endif