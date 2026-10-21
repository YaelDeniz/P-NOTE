#include "NeutrinoDetector.h"

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

//#include "PremModel.h"
//#include "PMNS_Fast.h"


// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"



using namespace std;

//RESOLUTION FUNCITONS----------------------------------------------------------

//Energy Resolution Parameters
void NeutrinoTelescope::SetEnergyResolution(double A, double B){
    Ae = A;
    Be = B;}

//Energy Resolution function
double NeutrinoTelescope::PDFE(double Ereco , double Etrue ) {
    double SigmaE = Ae*Etrue +Be*sqrt(Etrue); //Detector Energy Resolution
    double pdfe = ROOT::Math::gaussian_pdf( Ereco, SigmaE, Etrue) ;
    return pdfe; } //Gaussian Distribution

//Angular (Zenith only) Resolution Parameters
void NeutrinoTelescope::SetAngularResolution(double A, double B){
    Ath = A;
    Bth = B;}

//Angular(Zenith only) Resolution function
double NeutrinoTelescope::PDFth(double thReco, double Etrue, double thTrue){
    double Sigmath = Ath + Bth/(sqrt(Etrue)); //Detector Angular Resolution
    double pdfth = ROOT::Math::gaussian_pdf( thReco, Sigmath, thTrue) ;
    return pdfth; } //Gaussian Distribution

double NeutrinoTelescope::PDFcth(double cthreco, double Etrue, double cthtrue){

    /// @param sigma - std. dev. in radians. Angular resolution
    /// kappa = 1/sigma^2 - VMF concentrarion parameter. 

    double sigmath = Ath + Bth/(sqrt(Etrue)); // Angular resolution in radians

    double cos_z = cthreco;
    double cos_z0 = cthtrue;

    double sin_z = sqrt(1 - pow(cos_z, 2));
    double sin_z0 = sqrt(1 - pow(cos_z0, 2));

    double kappa = 1/(sigmath*sigmath);

    double exp_arg = kappa * cos_z * cos_z0;

    double norm = kappa;
    
    if (kappa < 50)
        norm /= 2 * sinh(kappa);
    else
        exp_arg -= kappa;

    double arg = kappa * sin_z * sin_z0;
    double bess = 1;
    if (arg < 50)
        bess = TMath::BesselI0(arg);
    else 
    {
        exp_arg += arg;
        bess *=
            (1 + 1. / (8 * arg) + 4.5 / pow(8 * arg, 2) + 37.5 / pow(8 * arg, 3));
        bess /= sqrt(2 * M_PI * arg);
    }

    double expo = exp(exp_arg);

    return norm * expo * bess * sin_z;} //Von Mises Fisher