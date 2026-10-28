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

#ifndef INTEVENTGENERATOR_H
#define INTEVENTGENERATOR_H

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

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>


//CERN ROOT
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "PhyTools.h"


using namespace std;


class NintSimulation
{
    private:

    fs::path PATH_fluxtables = fs::current_path() / "fluxtables";
     //The Earth model
    fs::path PATH_premtables = fs::current_path() / "premtables";

    public:
    

    std::string label = "std";

    double Exposure_det; //Detector Exposure
  
    
    double pos[3] = {0.0, 0.0, -1.0*R_EARTH}; //Detector position in the Earth model. Cartersian coordinates (x,y,z) in km




    std::string FluxFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string HondaTable = (PATH_fluxtables / FluxFile).string(); //Class Assumes South Pole flux as default


   
    std::string PremFile = "prem_44layers.txt";
    std::string PremTable = (PATH_premtables / PremFile).string();
 
    //The LLVP
    bool MantleAnomaly = false;
    double PileHeight = 1000;
    double aperture = 45.0;
    std::string AnomalyShape ="cake";
    double Pilerho_pct = 2.0;
    double PileChemContrast = 0.0;
    

    //Modify specific Layers
    int PremLayer = 44 ;
    double rho_pct = 0 ;
    double zoa_pct = 0 ;

    //Neutrino Settings
    int flvf = 0; // default: electron neutrino

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

    void SetDetPosition( double x[3]) //Detector position in the Earth model. Cartersian coordinates (x,y,z) in km
    {

        pos[0] = x[0];
        pos[1] = x[1];
        pos[2] = x[2];

    }

    void SetDetExposure(double MTonsYears_units) 
    {

        Exposure_det = (MTonsYears_units*MTonYear)/M_n; //MTon-Years

    }

    void ModifyLayer (int PremLayer_ , double rho_pct_ , double zoa_pct_)
    {

        PremLayer = PremLayer_; // Index of layer
        rho_pct = rho_pct_; // Density percentage difference
        zoa_pct = zoa_pct_; // zoa(compositional) perctengae differnce

    }

    void SetTrueIntervals(double Emin_, double Emax_, double thmin_, double thmax_, double phimin_, double phimax_)
    {
        thmin = thmin_*TMath::Pi()/180;
        thmax = thmax_*TMath::Pi()/180;

        phimin = phimin_*TMath::Pi()/180;
        phimax = phimax_*TMath::Pi()/180;
        
        Emin = Emin_;
        Emax = Emax_;

    }

    void SetTrueBinning(int i, int j, int k)
    {
        ibins = i;   //Energy bins
        jbins = j;   //Zenith bins
        kbins = k;   // Azimuth bins
    }

    void SetFlavor(int f) { flvf = f; }

    fs::path OutputDir;
    void SetOutputDir(const fs::path& dir) { OutputDir = dir; }

    std::vector< TH2D* >  GetNint_ijk();

};


#endif
