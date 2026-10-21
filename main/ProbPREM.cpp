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

//My Libraries
#include "GetProbPREM.h"


using namespace std;


int main()
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation Probability. " << std::endl;
    std::cout << " Generating Oscillation Probabilities. "<< std::endl; 

    // Number of bins
    int Ebins =101; // # of Bins of Energy
    int Tbins =101; // # of Bins of cosThetaZ

    //Range in Theta and E for Event Oscillogram

    //Energy interval (in GeV):
    double Emin = 1.0; 
    double Emax = 10.0;

    //Zenith Angle Interval:
    double Etamin = 33.0;
    double Etamax = 38.0;




    double E_GeV[] = {Emin,Emax};
    double Eta[] = {Etamin,Etamax};

     // neutrino state options nue (0), numu (1) or nutau (2)
    int flvi = 1;
    int flvf = 0;

    TH2D* Pab = GetProbPREM( flvi, flvf, E_GeV, Eta, Ebins, Tbins);

     std::cout << " Probability data is stored inside './SimulationResults'  " << std::endl;

    TCanvas *c = new TCanvas();
    


    Pab->Draw("colz");
    Pab->SetStats(0);
    gStyle->SetPalette(55);

    // Add space for the colz bar
    gPad->SetLogy();
    
   
    
    gPad->SetRightMargin(0.18);
    
    c->Print("test.png" );

    return 0;
    



}