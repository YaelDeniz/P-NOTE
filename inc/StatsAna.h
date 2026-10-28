
#ifndef  STATSANA_H
#define STATSANA_H


#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TMath.h"

# define myPi 3.14159265358979323846  /* pi */

# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg

using namespace std;
 
double Get2DChi2( TH2D * histstd, TH2D * histalt);

double Get3DChi2( TH3D * histstd, TH3D * histalt);

#endif