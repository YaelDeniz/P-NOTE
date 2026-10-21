#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
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
#include "../inc/NuRateGenerator.h"


# define myPi 3.14159265358979323846  /* pi */

using namespace std;


int main()
{
    // Plot and Histogram2D settings

    

    // Number of bins
    int Ebins =1; // # of Bins of Energy
    int Tbins =1; // # of Bins of cosThetaZ

    //Range in Theta and E for Event Oscillogram

    //Zenith angle interval:
    //double ctmin = cos(147.0*myPi/180); 
    //double ctmax = cos(142*myPi/180);

    double ctmin = cos(170.0*myPi/180); 
    double ctmax = cos(150*myPi/180);
    
    double dcT = (ctmax - ctmin)/(2.0*Tbins); //< Bin width/2

    //Energy interval (in GeV):
    //double Emin = 1.0; 
    //double Emax = 20.0;


    //Energy interval (in GeV):
    double Emin = 4.0; 
    double Emax = 6.0;

    double Es[] = {Emin,Emax};
    double cts[] = {ctmin,ctmax};

    int flvf = 1;

    TH2D* EventOsc = ObservedEvents(flvf, Es, cts , Ebins, Tbins);

    //TH2D*  ObservedEvents(int flvf, double Energy[], double CosT[] ,int Ebins, int Tbins)

    TCanvas *c = new TCanvas();
    
    //EventOsc->Draw("colz");    
    EventOsc->Draw("SURF1 Z");
    
    gStyle->SetPalette(55);

    gPad->SetTheta(30.0); // default is 30
    gPad->SetPhi(330.0); // default is 30 --Azimuthal
    gPad->Update();
    
    //gPad->SetLogx();
    
    gPad->SetRightMargin(0.18);

    
/*
    int flvi = 0 ;
    int flvf = 1 ;

    double Emin = 1;
    double Emax = 100;

    double L =  1000; // Neutrino Baseline in Km

    double rho = 0; // Vacuum;  
 
    TGraph* Pem =  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

    TCanvas *c = new TCanvas();

    Pem->Draw();

    c->Print("Pem_test.png");
    */
  
    return 0;
    



}
