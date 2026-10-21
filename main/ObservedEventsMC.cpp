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
#include "TObjArray.h"
#include "Math/Interpolator.h"
#include "Math/PdfFuncMathCore.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"
#include "GetObservedEvents.h"



# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need
# define R_earth 6368.0 //km
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;


double PDF_angle(double eta_o, double e,double eta, double a_a); // Angular resolution;

double PDF_energy(double  e_o,double e, double a_e);       // Energy resolution;



int main()
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;

    // Number of bins
    int K = 1; //Number of Pseudo experiments.
    int Ebins=10; // # of Bins of Energy
    int Etabins=10; // # of Bins of cosEta

    int Ebins_o=10; // # of Bins of Energy
    int Etabins_o=10; // # of Bins of cosEta

    int Bins[]={Etabins, Ebins};

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

    double R_min = 3000.0; // Distance from the center of the Earth , 3500 Km CMB
    //double R_llsvp = R_cmb + h_llsvp; //Km
    double R_max = 4500; // Distance from the center of the Earth

    double Etamin = TMath::ASin( (R_min)/R_earth )*(180.0/TMath::Pi()) ;
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    double Etamax = TMath::ASin( R_max/R_earth )*(180.0/TMath::Pi()) ;



    //Bins
    int ibins = Etabins; // Number of  angular bins of True event distribution
    int jbins = Ebins; // Number of  energy bins of True event distribution
    double deta = (Etamax - Etamin)/(ibins); //< Bin width
    double dE = (Emax - Emin)/(jbins); //< Bin width

    int mbins = Etabins_o; // Number of  energy bins of True event distribution
    int nbins = Ebins_o; // Number of  angular bins of True event distribution

    int Bins_o[]={Etabins, Ebins, Etabins_o, Ebins_o};
    //double deta = (etamax - etamin)/(2.0*ibins); //< Bin width/2
    //double dE = (Emax - Emin)/(2.0*jbins); //< Bin width/2


    //MOdify prem to fit H_LLSVP

    //--------------------------------------------------------------------------------------------

    //double Etamin =10 ;
    //double Etamax =30 ;


    double Phim = 0.0;
    double PhiM = 80.0 ;
    //double PhiM = 360.0 ; 

    // neutrino state options nue (0), numu (1) or nutau (2)
    int flvf = 1;

     //DETECTOR PROPERTIES
    double DetMass = 1.0*MTon; //Mass in megaton units
    
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    
    double T       = 1.0*years2sec; //Detector Exposure time in sec: One Year

    double NnT = Nn*T; 
    
    double a_eta=0.25;

    double a_e=0.2;

    double Det_par[]={a_eta,a_e};

    double dAz = PhiM-Phim;

    double Region[] = {Emin,Emax,Etamin,Etamax,dAz};

    //--- Generated Basic earth model---------------------------------------------------------------------------------------------------------------------

    //Density contrast
    double drho_dp = 2.0; // 2 % density contrats

    //Simulation Information
      std::cout << "Monte Carlo simulation for neutrino mu-like true events " << std::endl;
      std::cout << "Energy Range in GeV: [" << Emin << " - " << Emax << "]" << "Angular Range(zenith): [" << Etamin << " - " << Etamax << "]" <<std::endl;
      std::cout << "LLSVP information-  density contrats(%): " << drho_dp << " Height(km): "<< h_llsvp << std::endl;
      std::cout << "Simulation set up- Angular bins: " << Etabins << " Energy bins: "  << Ebins<< " # of seudo experiments: "<< K <<std::endl;
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



    //SMEAR EVENTS------------------------------------------------------------------------------------------------------
    //Smearing variables
    
    //prem_test   = "/home/dehy0499/OscProb/PremTables/prem_test.txt"; //Specify PREM table from OscProb

    std::cout<< "Smearing distribution for null hypothesis (Standard Earth:PREM)"<< std::endl;  
    
    prem_default   = "prem_default"; //Specify PREM table from OscProb

    TObjArray* PREM_Events = GetTrueEvents(prem_default, flvf, Region, Bins, NnT, K); //True event distribution

    TObjArray* PREM_Events_O = SmearEvents( PREM_Events, flvf, Region, Bins_o, Det_par, NnT, K ); //Add detector resolution to true events.

    /*

    std::cout<< "Smearing distribution for Alternative hypothesis (Alternative Earth:LLSVP)"<< std::endl;   

    prem_llsvp   = "prem_llsvp"; //Specify PREM table from OscProb
    
    TObjArray* Alt_Events = GetTrueEvents(prem_llsvp, flvf, Region, Bins, NnT, K); //Alternative event distribution

    TObjArray* Alt_Events_O = SmearEvents( PREM_Events, flvf, Region, Bins_o, Det_par, NnT, K ); //Add detector resolution to Alternative events.
    
 


    

    //------------------------------------------------------------------------------------------------------

    std::cout<< "ALT DATA--------------"<< std::endl;  

    prem_llsvp   = "prem_llsvp"; //Specify PREM table from OscProb
    
    TObjArray* Alt_Events = GetTrueEvents(prem_llsvp, flvf, Region, Bins, NnT, K); //True event distribution

    TObjArray* SmearExp_alt =  new TObjArray(K+1);

    TH2D* Alt_Events_O = new TH2D("Alt_Events_O","Observed Events (Alternative Earth Model):Kth pseudo-experiment; #eta ; E ",mbins,Etamin,Etamax,nbins,Emin,Emax); //Observed event distribution
    
    TH2D* Alt_Events_Ok = new TH2D("Alt_Events_Ok","Observed Events (Alternative Earth Model):Kth pseudo-experiment; #eta ; E ",mbins,Etamin,Etamax,nbins,Emin,Emax); //Observed event distribution

    TH2D* kthexp_a; //alt 

    //Smearing of True events-Default model

    for (int k = 0; k < K; ++k) //Loop in each sudo experiment.
    {
        // Observed Events 
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {   
            kthexp_a = (TH2D*)Alt_Events->At(k);

            eta_o = Alt_Events_Ok->GetXaxis()->GetBinCenter(m); // Select angular bin from the "Observed" histogram
           


            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                e_o = Alt_Events_Ok->GetYaxis()->GetBinCenter(n); // Select Energy bin from the "Observed" histogram
                
                //Smearing---------------------------------------------------------------------------------------------------------------------------
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    eta = kthexp_a->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small

                    for (int j = 1; j <= jbins; ++j)
                    {

                        e = kthexp_a->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small 

                        Nint  = kthexp_a->GetBinContent(i,j);
                        //Nint  = PREM_Events_k->GetBinContent(i,j);
                        
                        // bin width

                        n_mn  += Nint*PDF_angle(eta_o,e,eta, a_a)*PDF_energy(e_o,e,a_e)*deta*dE; 
                        
                    } // loop e

                } // Loop ct

                //----------------------------------------------------------------------------------------------------------------------------------

                Alt_Events_Ok->SetBinContent(m,n, n_mn); //Histogram for smeard psudo experiment.

                gbin = Alt_Events_O->GetBin(m,n); //Generalized mean

                Nmn = Alt_Events_O->GetBinContent(gbin)+ n_mn;

                Alt_Events_O->SetBinContent(m,n, Nmn);      //Expected Events

                n_mn = 0;

            } // loop in eo

        } // loop in eta_o

        SmearExp_alt->Add(PREM_Events_Ok);

    }

        SmearExp_alt->Add(PREM_Events_O);


    kthexp_o = (TH2D*)SmearExp_null->At(K); 

    kthexp_a = (TH2D*)SmearExp_alt->At(K);





    std::cout << "True data is stored inside './SimulationResults" << std::endl;


    // Data visualization

    TH2D* diffhist = new TH2D("diffhist"," Observed Events percentage difference; #eta ; E",mbins,Etamin,Etamax,nbins,Emin,Emax); 


   std::ofstream EventDiff("SimulationResults/ObsEventsResults/dNobs_test.csv"); 
   double nexp, nobs, dn, dnr;

    for(int m=1; m<= mbins  ; m++) 
    {    
        eta = kthexp_o->GetXaxis()->GetBinCenter(m); //< This will defined a constant L por different values of ct provided Dct is Small

                    
                    for (int n=1; n <= nbins ; n++)
                    { 
                        e = kthexp_o->GetYaxis()->GetBinCenter(n); //< This will defined a constant L por different values of ct provided Dct is Small
                        
                        nexp = kthexp_o->GetBinContent(m,n); // expected
                        nobs = kthexp_a->GetBinContent(m,n); // obserbed



                        dn = 100.0*((nexp-nobs)/nexp);
                        dnr = 100.0*((nexp-nobs)/nobs);

                        std::cout << nexp <<" "<< nobs <<" "<<dn << std::endl;

                        diffhist->SetBinContent(m,n, dn);

                        EventDiff<< eta << ", " << e << ", "<< nexp << ", "<< nobs<< ", " << dn << ", "<< dnr <<"\n";



                        
                    } // loop energy

    } // Loop eta

    EventDiff.close();
    */


    TH2D* hist_test;

    hist_test = (TH2D*)PREM_Events_O->At(0); 

    TCanvas *c = new TCanvas();

    hist_test->Draw("COLZ");
    
    c->Print("SimulationResults/Histograms/hist_test.png");


    return 0;

}
