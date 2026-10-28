#ifndef PHYTOOLS_H
#define PHYTOOLS_H

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

//ROOT

#include "TH2.h"

using namespace std;


// Earth structure radii (km)
# define R_CMB     3480.0               // Core-Mantle-Boundary radius
# define R_OCEAN   6368.0               // Ocean bottom / mean sea level radius
# define R_EARTH   6371.0               // Mean Earth radius
# define R_ATM     6386.0               // Atmosphere base radius

// Physical units
# define M_n        1.67492749804E-27   // Nucleon mass (kg, ~neutron mass)
# define MTon       1E9                 // 1 Megaton in kg
# define years2sec  3.154E7             // 1 year in seconds
# define MTonYear   3.154E16            // 1 Megaton-year in kg·s (= MTon * years2sec)


/*

Particle Physics tool for simulations

*/

//Neutrino Cross Sections
double XSec(double E, int ptype);

//Neutrino Flux Generation

class NuFlux
{

    public: 

    //Neutrino Energies
    double Enu_min= 1.0; //GeV
    double Enu_max= 100.0;

    //Detector Locations
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string HondaTable = FluxFolder + FluxFile; //Class Assumes Sout Pole flux as default

    bool DataHeader (const std::string &aline); // Get data into tabular form but ignore Headers from file

    std::vector<std::vector<double>> reshape(const std::vector<double>& vec, int a, int b); //Reshape colum vector into matrix

    std::vector< std::vector<double> > SetFluxData( std::string FluxTable  );
    
    TH2D*GetFluxHist(int flavor, std::vector<std::vector<double>> FluxData); // Create 2D Histogram

    
};




#endif
