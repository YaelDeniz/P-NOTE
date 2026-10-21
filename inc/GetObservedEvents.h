#ifndef OBSEVENTGENERATOR_H
#define OBSEVENTGENERATOR_H

//C++
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>

//CERN ROOT
#include "TH2.h"
#include "TMath.h"
#include "TH3.h"
#include "TObjArray.h"
#include "TGraph.h"

//Layer Updates
#include "LayerUpdate.h"


# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass


class AsimovObsSimulation
{
    public:

    //Detector
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
    double DensityContrast = 0;
    double ChemicalContrats = 0;

    //LLVP

    bool PileInModel = false; 
    std::string ShapeOfPile = "cake";
    double ThePileHight = 1000.0;
    double ThePileAperture=45.0;
    double ThePileDensityContrast = 10.0;
    double ThePileChemicalContrast = 10.0;

    //Neutrino Settings
    int flvf;
    double EnuMin = 1.00; //GeV
    double EnuMax = 10.00; //GeV

    double ZenMin; // Min 90;
    double ZenMax; // Max 180;

    double AziMin=0; // Min 0;
    double AziMax=360; //Max 360;

    double ZenithBinSize;
    double AzimuthBinSize;
    double EnergyBinSize;

    //Simulation Settings
    int truebinsZen; //Bins in Zenith
    int truebinsAzi; //Bins in Azimuth
    int truebinsE; //Bins in Energy

    int recobinsZen; //Bins in Zenith
    int recobinsAzi; //Bins in Azimuth
    int recobinsE; //Bins in Energy


    double MT; //Detector Exposure //0.05*e +0.1*sqrt(e)
    
    //Detector Resolution-------------------------------------------------------

    double Ae = 0.05;
    double Be = 0.1;

    double Ath = 0.04;  //0.04 + 0.18/(sqrt(e))
    double Bth = 0.18;

    double Aphi = 0.04;  //0.04 + 0.18/(sqrt(e))
    double Bphi = 0.18;

    //Patch 12/3/2025
    std::vector<LayerUpdate> LayerUpdates_recosim;
    void ChangeLayerProperties(int n, double rhopct, double zoapct);
    void AddUpdatesList(const std::vector<LayerUpdate>& updates);


    void SetEnergyResolution(double A, double B);

	double PDFE(double Ereco , double Etrue );//Gaussian

    void SetAngularResolution(double A1, double B1, double A2, double B2);

    double PDFth(double threco, double Etrue, double thtrue);

    double PDFNth(double threco, double Etrue, double thtrue);

    double PDFNphi(double phireco, double Etrue, double phitrue);

    double PDFangle(std::array<double, 2> recodir, std::array<double, 2> truedir, double Etrue);

    
    void ModifyLayer (int LayerNumber , double rhopct , double zoapct)
    {

        PremLayer = LayerNumber; // Index of layer
        DensityContrast = rhopct; // Density percentage difference
        ChemicalContrats= zoapct; // zoa(compositional) perctengae differnce

    }

    //Beware Angles are in degrees here
    void SetIntervals(double zlow, double zup, double alow, double aup, double elow, double eup)
    {
        ZenMin = zlow;
        ZenMax = zup;

        AziMin = alow;
        AziMax = aup;
        
        EnuMin = elow;
        EnuMax = eup;

    }

    void SetTrueBinning(int ebins, int zbins, int abins)
    {
        truebinsE   = ebins;
        truebinsZen = zbins;
        truebinsAzi = abins;
        
    }

    void SetRecoBinning(int ebins, int zbins, int abins) 
    {
        recobinsE   = ebins;
        recobinsZen = zbins;
        recobinsAzi = abins+1;
        
    }


    void SetExposure(double exposure) 
    {

        MT = exposure; // MTon Units

    }

    void SetDetectorXYZ(double xyz[3])
    {
        xyzTelescope[0]=xyz[0];
        xyzTelescope[1]=xyz[1];
        xyzTelescope[2]=xyz[2];

    }

    TH3D * GetRecoN(); //Generates Neutrino rates using full vMF

    //TH3D * GetObsEvents3D();

    //TH2D * GetObsEvents2Dcth(); //Azimuth independet

    //TH2D * GetRecoRate_vMFMarginal();

    TH2D * GetNreco_vMFMarginal( TH2D * Rreco) 
    {
        //Avarage over azimuth 
        TH2D* Nrecohist = new TH2D("Nrecohist_avg", "Azimuth-Averaged Reconstructed Events",
                        Rreco->GetXaxis()->GetNbins(), Rreco->GetXaxis()->GetXbins()->GetArray(),
                        Rreco->GetXaxis()->GetNbins(), Rreco->GetYaxis()->GetXbins()->GetArray());


        for (int m = 1; m <= Rreco->GetYaxis()->GetNbins(); ++m) {

            double dTheta_reco = Rreco->GetYaxis()->GetBinWidth(m);

            for (int l = 1; l <= Rreco->GetXaxis()->GetNbins(); ++l) {
                
                double dE_reco = Rreco->GetXaxis()->GetBinWidth(l);

                double Rlm = Rreco->GetBinContent(l, m); // reconstructed rate

                double Nlm = (MT / mN) * Rlm * dTheta_reco * dE_reco; 
                
                Nrecohist->SetBinContent(l, m, Nlm);
            } //Reco energy
        } // reco Zenith
   
   return Nrecohist;}



    std::vector < TH2D* >  GetObsEvents3Dcth();
    std::vector < TH2D* >  GetObsEvents3Dth();
    std::vector < TH2D* >  GetObsEvents3Dthvb(); //Variable bin
    std::vector < TH2D* >  GetObsEvents3Dthlog10E();


    //TH2D * TestObsEvents2D(std::string model_std , std::string model_alt);

   void GetDetails(const std::string& PATH)
   {

        std::ofstream details(PATH);
    
        if (!details) {
            std::cout << "Error: Unable to open file " << PATH << std::endl;
            return;
        }

        details << "--------------------------------------------------------------------------------\n";
        details << "                               REALISTIC DETECTOR                               \n";
        details << "--------------------------------------------------------------------------------\n";
        details << "                                Detector Details                               \n";
        details << "--------------------------------------------------------------------------------\n";
        details << " Detector Radius: " << xyzTelescope[2] << " m\n";
        details << " Neutrino Flux: " << FluxFile << " \n";
        details << " Number of Detectors: 1\n";
        details << " Mass: " << 10.0 << " Mton\n";
        details << " Time: " << 10.0 << " years\n";
        details << " Exposure: " << MT << " Mton-Year\n";
        details << " Energy resolution Ae:" << Ae << " Be:" << Be << "\n";
        details << " Zenith resolution Ath:" << Ath << " Bth:" << Bth << "\n";


        details << "\n--------------------------------------------------------------------------------\n";
        details << "                                Earth Details                                  \n";
        details << "--------------------------------------------------------------------------------\n";
        details << " Earth Model: " << PremName << "\n";
        details << "Index of modified layer (default is 0): " << PremLayer << "\n";
        details << "Density Contrast (drho/rho): " << DensityContrast<< "\n";
        details << " Chem Contrast (H%): " << ChemicalContrats << "\n";
        
        details << "\n--------------------------------------------------------------------------------\n";
        details << "                        Lower Mantle Anomalies (LLVP)                          \n";
        details << "--------------------------------------------------------------------------------\n";
        details << " LLVP: " << PileInModel << "\n";
        details << " LLVP Model: " << ShapeOfPile << "\n";
        details << " LLVP Height: " << ThePileHight << " km\n";
        details << " Aperture: " << ThePileAperture << "\n";
        details << " Density Contrast (drho/rho): " << ThePileDensityContrast << "\n";
        details << " Chem Contrast (H%): " << ThePileChemicalContrast << "\n";

        details << "\n--------------------------------------------------------------------------------\n";
        details << "                             Simulation Details                                \n";
        details << "--------------------------------------------------------------------------------\n";
        
        details << " Neutrino Energy Range: [" << EnuMin << " - " << EnuMax << "] GeV | Reco Bins: " << recobinsE << " True bins:" << truebinsE <<  "\n";
         
        details << " Energy bin size:" << EnergyBinSize << "\n";
        
        details << " Zenith Range: [" << ZenMin << " (" << cos(ZenMin*TMath::Pi()/180.0) << ") - " << ZenMax << " (" << cos(ZenMax*TMath::Pi()/180.0) << ")] | Reco Bins: " << recobinsZen << " True bins:" << truebinsZen << "\n";
        
        details << " Zenith bin size:" << ZenithBinSize << "\n";

        details << " Azimuth Range: [" << AziMin << " - " << AziMax << "] Reco Bins: " << recobinsAzi << " True bins:" << truebinsAzi <<  "\n";
        
        details << " Azimuth bin size:" << AzimuthBinSize << "\n";
        
        details.close();

   }


};

//TH2D*  AsimovObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT);

//TH2D*  DetectorRes( double Region[], int Bins[],double Det_par[] );

//TH2D*  GetObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT);

//TObjArray*  SmearEvents( TObjArray* TrueEvents, int flvf, double Region[], int Bins[],double Det_par[], double NnT, int K );


//Generate energy dependent Oscillation Probabilities  (fixed Baseline) 
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif