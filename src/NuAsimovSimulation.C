/*
Notes:

All angular inputs are in degrees

All angular variables are in radias

*/
#include "NuAsimovSimulation.h"

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

//Random number generator
#include <random>


//Cern Root
#include <math.h>
#include "TF1.h"
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include <Math/Interpolator.h>
#include "TGraph.h"
//#include "TGraphEventsD.h"

//OSCPROB
#include "PremModel.h"
#include "PMNS_Fast.h"

// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"
#include "Earth3DModel.h"
#include "ObsEventGen.h"
#include "TrueEventGen.h"

using namespace std;

// Some Constants i need
//# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
//# define MTon  1E9  //Metric MegaTon
//# define years2sec 3.154E7 // Years in Seconds

#define MyPremTables "/home/dehy0499/OscProb/PremTables/"


TH2D* NuAsimovSimulation::IntEventsGen2D( ) //To be Deleted
{   
    AsimovSimulation TrueSimulation;
    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;
    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(ZenMin,ZenMax,PhiMin,PhiMax,EnuMin,EnuMax);
    TrueSimulation.SetBinning(truebinsZen,trueBinsAzi,truebinsE);
    TrueSimulation.SetExposure(NnT);
    TrueSimulation.flvf=flvf;
    TH2D * EventHist = TrueSimulation.GetTrueEvents2D();

    return EventHist;
}
