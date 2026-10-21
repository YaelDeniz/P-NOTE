#ifndef NEUTRINOTELESCOPE_H
#define NEUTRINOTELESCOPE_H

//C++
#include <string.h>
#include <stdlib.h>
#include <vector>

//CERN ROOT
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"
#include "TGraph.h"

# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)

class NeutrinoTelescope
{ 
 
   //Detector Resolution-------------------------------------------------------

    double Ae = 0.05;
    double Be = 0.1;

    double Ath = 0.04;  //0.04 + 0.18/(sqrt(e))
    double Bth = 0.18;

    void SetEnergyResolution(double A, double B);

	double PDFE(double Ereco , double Etrue );//Gaussian

    void SetAngularResolution(double A, double B);

    double PDFth(double thReco, double Etrue, double thTrue); //Gaussian

    double PDFcth(double cthreco, double Etrue, double cthtrue); //Von Mises Fisher





};

#endif