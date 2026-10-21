#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need
# define R_earth 6368.0 //km
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;


int main()
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;

    // Number of bins
    
    int Ebins=200 ; // # of Bins of Energy
    int Tbins=200 ; // # of Bins of cosEta

    int Bins[]={Tbins, Ebins};

    //Range in Theta and E for Event Oscillogram

    //Energy interval (in GeV):
    double Emin=1.0 ; 
    double Emax=10.0 ;

    //Zenith Angle Interval:
    //double Etamin = 33.0;
    //double Etamax = 38.0;

    // Alternative Earth Model:

    double h_llsvp = 700; //km
    double R_cmb = 3500.0;
    double R_llsvp = R_cmb + h_llsvp; //km

    double R_min = 3000.0;
    //double R_min = 3500.0; // Distance from the center of the Earth , 3500 Km CMB
    //double R_llsvp = R_cmb + h_llsvp; //Km
    double R_max = 4500; // Distance from the center of the Earth

    
    double Etamin = TMath::ASin( (R_min)/R_earth )*(180.0/TMath::Pi()) ;
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    double Etamax = TMath::ASin( R_max/R_earth )*(180.0/TMath::Pi()) ;

    //double Etamin = 0.0 ;
    //double Etamax = 180.0 ;
    

    //double Etamin = 10 ;
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    //double Etamax = 30 ;


    //MOdify prem to fit H_LLSVP

    //--------------------------------------------------------------------------------------------

    //double Etamin =10 ;
    //double Etamax =30 ;


    double Phim = 0.0;
    double PhiM = 80.0 ;
    //double PhiM = 360.0 ; 

     //DETECTOR PROPERTIES
    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year

    double NnT = Nn*T; 


    /*  
    double Emin = 4.0; 
    double Emax = 6.0;

    double Etamin = 10.0;
    double Etamax = 30.0;
    
    double Phim    = 0.0;
    double PhiM    = 360.0 ;
    */

    double dAz = PhiM-Phim;


    double Region[] = {Emin,Emax,Etamin,Etamax,dAz};

    //Density contrast
    double drho_dp = 2.0; // 2 % density contrats

    // neutrino state options nue (0), numu (1) or nutau (2)
    int flvf = 1;


    //

    std::string prem_llsvp, prem_default;
    
    std::ifstream PREM_DEFAULT;

    //Density profile PREM altered--------------------------------------------------------------------------

    PREM_DEFAULT.open("/home/dehy0499/OscProb/PremTables/prem_default.txt");
               
    std::ofstream PREM_LLSVP("/home/dehy0499/OscProb/PremTables/prem_llsvp.txt"); //Create a file with Alternative Earth density profile 
    
    // Variables for storing table rows
    float radius, density, zoa, layer;

    
    // Loop over table rows
    while(PREM_DEFAULT >> radius >> density >> zoa >> layer)
    {

        if( radius >=  R_cmb  && radius <= R_llsvp )
            {

                PREM_LLSVP << radius << " " << density*(1.0 + drho_dp/100) << " " << zoa << " " << layer << "\n";
                            

            }
                        
        else 
            {

                PREM_LLSVP << radius << " " << density << " " << zoa << " " << layer << "\n";

            }

    }

    PREM_DEFAULT.close();
    PREM_LLSVP.close();



    //------------------------------------------------------------------------------------------------------
    //std::ifstream StdPrem;
    
    
    //prem_test   = "/home/dehy0499/OscProb/PremTables/prem_test.txt"; //Specify PREM table from OscProb
    std::cout<< "PREM DATA--------------"<< std::endl;  
    
    prem_default   = "prem_default"; //Specify PREM table from OscProb
    TH2D* nullhist = OscProbEarth(prem_default, flvf , Region,  Bins,  NnT);

    std::cout<< "ALT DATA--------------"<< std::endl;  

    prem_llsvp   = "prem_llsvp"; //Specify PREM table from OscProb
    TH2D* althist = OscProbEarth(prem_llsvp, flvf , Region,  Bins,  NnT);


    

    return 0;

}