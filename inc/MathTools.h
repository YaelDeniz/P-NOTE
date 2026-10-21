
#ifndef MATHTOOLS_H
#define MATHTOOLS_H

/* 

This files contains declarations and provides Interface to other 
parts of a program 

*/


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <valarray>

//CERN ROOT libraries
#include "TH2.h"
#include "TH3.h"










void GetDiff2D( TH2D * histstd , TH2D * histalt, TH2D * diff); // Create a histogram containing the difference between oscillograms

double IntChi2( TH3D * histstd, TH3D * histalt); // Calculate chi2 values in a 3 variable binning scheme

double Get2DChi2( TH2D * histstd, TH2D * histalt); // Calculate chi2 values in a 2 variable binning scheme

//Factorial of a function
double factorial(const int& n);

//Generates random Samples from Poisson Distribution
double  Sample_Pois(double lambda, double u);

// Numericla Integration

double ncquad(std::vector<double> xdat, std::vector<double> ydat); //Composite Simpson's Rule

double simpson(std::vector<double> x, std::vector<double> f); //Simpsons Rule 1/3



//Logarithm Spaced points (Even spaced points in LogSpace)
std::vector<double> logspace(double a, double b, int k);

#endif
