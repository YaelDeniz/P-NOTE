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
#include "Math/ProbFuncMathCore.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
//#include "GetObservedEvents.h"
#include "GetTrueEvents.h"



# define myPi 3.14159265358979323846  /* pi */

# define R_earth 6368.0 //km

// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;


int main()
{
    
    std::cout << " Neutrino Oscillation tomography. " << std::endl;
    std::cout << " Generating Log-Likelihood Ratio Test (LLRT) / Chi^2 test. "<< std::endl; 

    
    // Alternative Earth Model:-------------------------------------------------------------

    double h_llsvp = 700; //km
    double R_cmb = 3500.0; //Km
    double R_llsvp = R_cmb + h_llsvp; //Km

    double rho_min = -2;
    double rho_max = 2;

    int m = 10; // # Number of points

    double rho_arr [m+1]; 

    double drho_arr = (rho_max-rho_min)/(m);

    for (int i = 0; i <= m; ++i)
    {

        rho_arr [i] = rho_min + i*drho_arr;

        std::cout<< rho_arr[i] << std::endl; 
        
    }




    //Detector-------------------------------------------------------------------------------

    //Energy interval (in GeV):
    double Emin = 3.0; //GeV
    double Emax = 7.0; //GeV

    //Zenith and Azimuthal Angle Intervals (deg):

    double Etamin = TMath::ASin( (R_cmb)/R_earth )*(180.0/TMath::Pi()) ;
    double Etamax = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;

    double Phim = 0.0;
    double PhiM = 80.0 ; 
    double dAz = PhiM-Phim;


    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    //double T       = 30.0*years2sec; //Detector Exposure time in sec: One Year


    double E[] = {Emin,Emax};

    double Eta[] = {Etamin,Etamax};


    //Detector Parameters

    double aE = 0.2;

    double aTh = 0.25;
    


    //STATS--------------------------------------------------------------------------------------------

    // Number of bins
    int Ebins =9; // # of Bins of Energy
    int Tbins =4; // # of Bins of cosThetaZ

    //Final neutrino Flavor nue(0) numu(1) nutau(2):
    int flvf = 1;

    // Get the model table paths ------------
    std::string Path2Original;
    std::ifstream StdPrem;
    ofstream ChiSqrtVals("SimulationResults/Chi2test.csv", std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file


    //Loops---------------------------------------------------------------------------------------------

    for (int t = 1; t <= 4 ; ++t)
    {   
        double T       = t*10.0*years2sec; //Detector Exposure time in sec: One Year
        double NnT=Nn*T; //DETECTOR EXPOSURE IN MTON*YEARS
        double Det_par[] = {aE,aTh,NnT};

        for (int i = 0; i <= m; ++i)
        {
                
                StdPrem.open("/home/dehy0499/OscProb/PremTables/prem_default.txt");
               
                std::ofstream AltPrem("/home/dehy0499/OscProb/PremTables/prem_alt.txt"); //Create a file with Alternative Earth density profile 
    

                double drho_dp = rho_arr[i];

                std::cout << "density contrats:" << drho_dp << std::endl;


                // Variables for storing table rows
                float radius, density, zoa, layer;

                // Keep track of previous radius
                double rprev = 0;


                // Loop over table rows
                while(StdPrem >> radius >> density >> zoa >> layer){


                    if( radius > R_cmb  && radius < R_llsvp )
                    {


                        AltPrem << radius << " " << density*(1.0 + drho_dp/100) << " " << zoa << " " << layer << "\n";
                        

                    }
                    
                    else 
                    {

                        AltPrem << radius << " " << density << " " << zoa << " " << layer << "\n";

                    }

                }

                AltPrem.close();
                StdPrem.close();

                //---------------------------------------

                

                // Null hypothesis: Standard Earth

                std::string prem_default;
                
                prem_default   = "/home/dehy0499/OscProb/PremTables/prem_default.txt"; //Specify PREM table from OscProb

                // Alternative hypothesis

                std::string prem_alt;
                
                prem_alt   = "/home/dehy0499/OscProb/PremTables/prem_alt.txt"; //Specify PREM table from OscProb

                // Generate Observed events for two different models/hypothesis

                //TH2D* EventsExp = GetObservedEvents(prem_default,flvf, E, Eta , dAz , Ebins, Tbins, Det_par ); //Expected data from Null

                //TH2D* EventsObs = GetObservedEvents(prem_alt,flvf, E, Eta , dAz , Ebins, Tbins, Det_par ); //Observed data from Alternative


                TH2D* EventsExp = GetTrueEvents(prem_default,flvf, E, Eta , dAz , Ebins, Tbins, NnT); //Expected data from Null

                TH2D* EventsObs = GetTrueEvents(prem_alt,flvf, E, Eta , dAz , Ebins, Tbins, NnT); //Observed data from Alternative
                   
                //Chi2 Analysis: Log-likelihood Analysis.

                double dChi2 = 0;

                for(int m=1; m <= Tbins  ; m++) 
                {    

                    
                    for (int n = 1; n <= Ebins ; ++n)
                    { 
                        double nexp = EventsExp->GetBinContent(m,n); // expected
                        double nobs = EventsObs->GetBinContent(m,n); // obserbed

                        dChi2 +=  2*( nexp - nobs + nobs*TMath::Log(nobs/nexp) ); // LLRT

                        //dChi2 +=  (nexp - nobs)*(nexp - nobs)/(nexp); // Pearson's Chi^2

                        //dChi2 +=  (nexp - nobs)*(nexp - nobs)/(nobs); // Neyman's Chi^2


                        
                    } // loop energy

                } // Loop eta

                 //P-Value calculation

                int r = Tbins*Ebins - 2; //Degrees of freedorm

                double pval_c = ROOT::Math::chisquared_cdf_c  (dChi2, r) ;
                double pval = ROOT::Math::chisquared_cdf   (dChi2, r);  

                std::cout << "exposure" << t*10 <<std::endl;

                std::cout << "Chi^2 Result = " << dChi2 << std::endl;       

                std::cout << "Pvalue 1 = " << pval_c << " " << "Pvalue 2= " << pval<< std::endl;

                std::cout << " " << std::endl;

                ChiSqrtVals << t*10 <<" " << drho_dp << " "<< dChi2 << " "<< Ebins << " "<< Tbins <<"\n";

                delete EventsExp;
                delete EventsObs;

        }//End Loop

    }

    ChiSqrtVals.close();


/*            
 // It work
            //Quick Test

            double ni [] = {6, 23, 29, 31, 27, 13, 8, 13}; //Observed
            double Ei [] = {5.54, 18.26, 30.12, 33.14, 27.35, 18.05, 9.93, 7.63}; //Expected

            int r_o = sizeof(ni)/sizeof(ni[0]); //length calculation

            double dChi2 = 0;

            for(int i=0; i < r_o  ; i++) 
            {    

                  dChi2 +=  (ni[i] - Ei[i])*(ni[i] - Ei[i])/(Ei[i]); // Chi^2
                std::cout << ni[i] << std::endl;

            } // Loop eta
            

            //P-Value calculation

            int r = r_o - 2; //Degrees of freedorm

            double pval_c = ROOT::Math::chisquared_cdf_c  (dChi2, r) ;
            double pval = ROOT::Math::chisquared_cdf   (dChi2, r);  

             std::cout << "Chi^2 Result = " << dChi2 << std::endl;       

            std::cout << "Pvalue 1 = " << pval_c << " " << "Pvalue 2= " << pval<< std::endl;
            */

    return 0;

}