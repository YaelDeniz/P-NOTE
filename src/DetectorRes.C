#include "GetObservedEvents.h"

//C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//Cern Root
#include <math.h>
#include "TF1.h"
#include "TMath.h"
#include "TH2.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Interpolator.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/RootFinder.h"
//#include "TGraphTrueD.h"

//OSCPROB
#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif

// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"

using namespace std;


// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

# define R_earth 6368.0 //km

//Resolution functions

double PDF_energy(double e_o , double e , double a_e);

double PDF_angle(double eta_o, double e, double eta, double a_eta);

//double w_E(double e_o , double e, double dE_o, double a_E, double E[] ); 

//double w_eta(double eta_o, double e, double eta, double deta_o,double a_eta ,double Eta[]);

// Make oscillogram for given final flavour and MH


TH2D*  DetectorRes( double Region[], int Bins[],double Det_par[] )
{   
    std::cout << "Detector Resolution function" << std::endl;

    std::string location = "SimulationResults/ResolutionFunctions/" ;
   
    std::string title = "Resolution"+std::to_string(Bins[2])+"_"+std::to_string(Bins[3])+"_"+std::to_string(Region[0])+"_"+std::to_string(Region[1])+".csv";
    std::string filename = location+title;
    
    ofstream DetRes(filename, std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file
    

    //Detector Response-------------------------------------------------------------------------------------------------
    
    double a_E =  Det_par[0] ;
    
    double a_theta = Det_par[1] ;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = Bins[2]; //Number of angular bins of Observed event distribution
    
    int nbins = Bins[3]; //Number of Energy bins of Observed event distribution
    
    //Energy[GeV]
    double EOmin      = Region[0];
    
    double EOmax      = Region[1];
    
    double dE_o = (EOmax - EOmin)/(nbins); //< Bin width
    
    double E_itv[] = {EOmin, EOmax};
    
    //Observed Angle
    double thetaOmax   = (180.0-Region[2])*TMath::Pi()/180.0;
    
    double thetaOmin   = (180.0-Region[3])*TMath::Pi()/180.0;

    double ctOmin   = cos(thetaOmax);
    
    double ctOmax   = cos(thetaOmin);
    
    double dtheta_o = (thetaOmax - thetaOmin)/(mbins); //< Bin width

    double Theta_itv[] = {thetaOmin, thetaOmax};
    
    //Azimuthal 
    double dAz = Region[4]*TMath::Pi()/180.0;

    //Observed distribution:
    TH2D* hObs = new TH2D("hObs","Observed Variable",mbins,ctOmin,ctOmax,nbins,EOmin,EOmax); // xbins correspond to energy values and ybins to zenith angle cost
    
    //TH2D* hObs = new TH2D("hObs","Observed Variable",mbins,thetaOmin,thetaOmax,nbins,EOmin,EOmax); // xbins correspond to energy values and ybins to zenith angle cost
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = Bins[0]; // Number of  angular bins of True event distribution

    int jbins = Bins[1]; // Number of  energy bins of True event distribution


    
      //True Energy[GeV]
    //double Emin=(1/(1+4*a_E))*(EOmin);

    //double Emax=(1/(1-4*a_E))*(EOmax);
    TF1 *fmin = new TF1("fmin","-1*x + [0]-4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,EOmin);
    fmin->SetParameter(1,0.05);
    fmin->SetParameter(2,0.1);

    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    double Emin = finder_min.Root();

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,EOmax);
    fmax->SetParameter(1,0.05);
    fmax->SetParameter(2,0.1);

    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double Emax = finder_max.Root();

    
    //NextGen Energy


    
    double dE = (Emax - Emin)/(jbins); //< Bin width

    //True Angle
    /*
    double etamin_a= (etaOmin)-4*(a_eta/sqrt(Emin));
    double etamin_b= 0;

    double etamax_a= (etaOmax)+4*(a_eta/sqrt(Emin));
    double etamax_b= TMath::Pi();
    */

    double tpar_a = 0.04;
    double tpar_b = 0.18;

    double thetamin_a= (thetaOmin)-4*(tpar_a + tpar_b/sqrt(Emin));
    double thetamin_b= 0;

    double thetamax_a= (thetaOmax)+4*(tpar_a + tpar_b/sqrt(Emin));
    double thetamax_b= TMath::Pi();


    double thetamin = max(thetamin_b,thetamin_a);

    double thetamax = min(thetamax_b,thetamax_a);

    double dtheta = (thetamax - thetamin)/(ibins); //< Bin width
    
    
    //True distribution        
    TH2D* hTrue = new TH2D("hTrue","True Variable",jbins,thetamin,thetamax,ibins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost


    std::cout << "Observed variable region ["<< EOmin << "-"<<EOmax<< "]*["<<thetaOmin<< "-"<<thetaOmax <<"]" << std::endl;

    //std::cout << "Extended variable region ["<< Emin << "-"<<Emax<< "]*["<<etamin_a<< "-"<<etamax_a <<"]" << std::endl;

    std::cout << "True variable region ["<< Emin << "-"<<Emax<< "]*["<<thetamin<< "-"<<thetamax <<"]" << std::endl;

    
    



   //Calculation of detector resolution------------------------------------------------------------------------------------- 
    for (int m = 1; m <= mbins; ++m) //Observed Angular bins
    {

        double ct_o = hObs->GetXaxis()->GetBinCenter(m); // Select angular bin from the "Observed" histogram

        double theta_o = TMath::ACos(ct_o);
       
        for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
        {
     
            double e_o = hObs->GetYaxis()->GetBinCenter(n); // Select Energy bin from the "Observed" histogram
       
            double Ntot = 0;   // Reset total number of events.
            
           // for(int i=1; i<= ibins ; i++) //Loop in cosT
            //{    

                // Get cos(theta_z) from bin center, It fix a value of L

              //  double theta = hTrue->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                
                // double cosT =cos(theta);

                //if(cosEta < -1 || cosEta > 1) break; // Skip if cosT is unphysical 
                
                

 

                //for (int j = 1; j <= jbins; ++j)
                //{ 
                    //NUMERICAL INTEGRATION

                  ///  double e = hTrue->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small

            		double theta = 134.0*TMath::Pi()/180.0;

            		double e = 5;

                    double f_mn = PDF_angle(theta_o,e,theta,a_theta)*PDF_energy(e_o,e,a_E);


                    
                   
                    
                //} // loop e

           // } // Loop ct
            

            DetRes << theta_o << ", " << e_o << ", "<< f_mn << "\n"; //Write im file.

            
            hObs->SetBinContent(m,n,f_mn); // Set events in Observed histogram. 

         
        } // loop in eo

    } // loop in eta_o


   DetRes.close();
       
   hObs->SetTitle(" detector resolution");
    
   hObs->SetStats(0);

   //Delete 
   delete hTrue;

   return hObs;
}


//double w_E(double a_E, double E , std::vector<double> eo, double E_GeV[] ) 

double PDF_energy(double e_o , double e, double a_e ) 
{

//double sigma_ene = a_e*e;

double sigma_ene = 0.05*e +0.1*sqrt(e);

double pdf_ene = ROOT::Math::gaussian_pdf( e_o, sigma_ene, e) ;
 
return pdf_ene; 



}

//double w_eta(double a_eta, double eta, double E, std::vector<double> eta_o, double Eta[])

double PDF_angle(double eta_o, double e, double eta, double a_eta)
{

//double sigma_angle = a_eta/(sqrt(e));

double sigma_angle = 0.04 + 0.18/(sqrt(e));

double pdf_angle = ROOT::Math::gaussian_pdf     ( eta_o, sigma_angle, eta) ;
 
return pdf_angle; 


}
