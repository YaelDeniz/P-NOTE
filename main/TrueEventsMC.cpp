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
    int K = 1000; //Number of Pseudo experiments.
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

    double R_min = 3500.0; // Distance from the center of the Earth , 3500 Km CMB
    //double R_llsvp = R_cmb + h_llsvp; //Km
    double R_max = 4500; // Distance from the center of the Earth

    double Etamin = TMath::ASin( (R_min)/R_earth )*(180.0/TMath::Pi()) ;
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    double Etamax = TMath::ASin( R_max/R_earth )*(180.0/TMath::Pi()) ;

    //MOdify prem to fit H_LLSVP

    //--------------------------------------------------------------------------------------------

    //double Etamin =10 ;
    //double Etamax =30 ;


    double Phim = 0.0;
    double PhiM = 80.0 ;
    //double PhiM = 360.0 ; 

     //DETECTOR PROPERTIES
    double DetMass = 1.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 1.0*years2sec; //Detector Exposure time in sec: One Year

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

    //Simulation Information
      std::cout << "Monte Carlo simulation for neutrino mu-like true events " << std::endl;
      std::cout << "Energy Range in GeV: [" << Emin << " - " << Emax << "]" << "Angular Range(zenith): [" << Etamin << " - " << Etamax << "]" <<std::endl;
      std::cout << "LLSVP information-  density contrats(%): " << drho_dp << " Height(km): "<< h_llsvp << std::endl;
      std::cout << "Simulation set up- Angular bins: " << Tbins << " Energy bins: "  << Ebins<< " # of pseudo experiments: "<< K <<std::endl;
      std::cout << "PREM tables located in /OscProb/PremTables"<< std::endl;


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
    TObjArray* PREM_Events = GetTrueEvents(prem_default, flvf, Region, Bins, NnT, K);

    std::cout<< "ALT DATA--------------"<< std::endl;  

    prem_llsvp   = "prem_llsvp"; //Specify PREM table from OscProb
    TObjArray* Alt_Events = GetTrueEvents(prem_llsvp, flvf, Region, Bins, NnT, K);

    std::cout << "True data is stored inside './SimulationResults" << std::endl;

    /*
    TCanvas *c = new TCanvas();
    TH2D*e;
    e=(TH2D*)PREM_Events->At(10);
    e->Draw();
    c->Print("histtestLAST.png");
    */


    TH2D* nullhist;
    TH2D* althist;

    nullhist=(TH2D*)PREM_Events->At(2);
    althist=(TH2D*)Alt_Events->At(2);

    TH2D* diffhist = new TH2D("diffhist"," Events percentage difference; #eta ; E",Tbins,Etamin,Etamax,Ebins,Emin,Emax); 

    // Data visualization
   
   std::ofstream EventDiff("SimulationResults/TrueEventsResults/perdiff_MCTest200200.csv"); 
   double eta, e, nexp, nobs, dn;

    for(int i=1; i<= Tbins  ; i++) 
    {    
        eta = nullhist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small

                    
                    for (int j=1; j <= Ebins ; j++)
                    { 
                        e = nullhist->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small
                        
                        nexp = nullhist->GetBinContent(i,j); // expected
                        nobs = althist->GetBinContent(i,j); // obserbed

                        dn = 100.0*abs ((nexp-nobs))/nexp;

                        diffhist->SetBinContent(i,j, dn); //Create histogram for  kth Pseudo-Experiment

                        EventDiff<< eta << ", " << e << ", "<< nexp << ", "<< nobs<< ", " << dn << "\n";



                        
                    } // loop energy

    } // Loop eta

    EventDiff.close();


    TCanvas *c = new TCanvas();
    diffhist->Draw("COLZ");
    c->Print("SimulationResults/Histograms/perdiff_MCTestmean2.png");
    

    

    return 0;

}