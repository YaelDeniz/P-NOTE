#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "../inc/MathTools.h"
#include "../inc/PhyTools.h"
#include "../inc/GetObservedEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds
# define R_earth 6368.0 //km



using namespace std;


int main()
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;
    std::cout << " Observed Events: track-like only. "<< std::endl; 

    int Ebins=30; // # of Bins of Energy
    int Etabins=30; // # of Bins of cosEta

    int Ebins_o=100; // # of Bins of Energy
    int Etabins_o=100; // # of Bins of cosEta

    int Bins[]={Etabins, Ebins, Etabins_o, Ebins_o};


    //Energy interval (in GeV):
    double Emin = 1.0; 
    double Emax = 20.0;

    //------------------


    double thetamin = 0*(180.0/TMath::Pi()) ;
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    double thetamax = 180*(180.0/TMath::Pi()) ;




    //------------------

    double Phim = 0.0;
    double PhiM = 80.0 ; 

    double Det_par[] = {0.2,0.25};
     

    double dAz = PhiM-Phim;

    double Region[] = {Emin,Emax,thetamin,thetamax,dAz};

    
    //-------------------------------------------------------------------------------------------------------------------

    TH2D* Detector_res = DetectorRes( Region, Bins , Det_par );



    TCanvas *c = new TCanvas();

    Detector_res->Draw("COLZ");

    
    c->Print("SimulationResults/ResolutionFunctions/test.png");

    return 0;

}


