
#ifdef  STATSANA_H
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
// Some Constants i need

# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6368.0 // Radius of the Atmosphere (km)
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds

using namespace std;
 
double Get2DChi2( TH2D * histstd, TH2D * histalt);

double Get3DChi2( TH3D * histstd, TH3D * histalt);

#endif