#include "Earth3DModel.h"

//C++
#include <ctime>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>  

#include "TGraph.h"
#include "TMultiGraph.h"

#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2D.h"


//OSCPROB

#include "PremModel.h"

//Some import Earth Radius
# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)

//Some units
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds

void ExportToCSV(TH2D* hist, std::string filename);


int main(int argc, char **argv)
{


  // Simulation Configuration
  int zenithBins = 100;                      // Zenith bins
  int azimuthBins = 100;                     // Azimuth bins


  Earth3DModel Baselines;

  Baselines.PileThickness = 1000;

  Baselines.aWidth = 45;

  std::string llvpShape = "pancake";
  double llvpRadius = Rcmb + Baselines.PileThickness;

  double depthMin = Rcmb - 2500;
  double depthMax = llvpRadius + 500;        // Max distance from Earth's center


  // Zenith angle [rad]
  double zenithMin = 180 - TMath::ASin(depthMax / Rearth) * (180.0 / TMath::Pi());
  double zenithMax = 180 - TMath::ASin(depthMin / Rearth) * (180.0 / TMath::Pi());

  double cosZenithMin = cos(zenithMax * TMath::Pi() / 180.0);
  double cosZenithMax = cos(zenithMin * TMath::Pi() / 180.0);

  // Azimuthal angle [rad]
  double azimuthStart = -50;
  double azimuthEnd = 50;
  double dAzimuth = (azimuthEnd - azimuthStart) / (2.0 * azimuthBins);
  double azimuthMin = (azimuthStart - dAzimuth);
  double azimuthMax = (azimuthEnd + dAzimuth);

  TH2D* LLVPBaseLineHIst = new TH2D("LLVPBaseLine", "LLVP Base Line", zenithBins, cosZenithMin, cosZenithMax, azimuthBins, azimuthMin, azimuthMax);

  double cthi = LLVPBaseLineHIst->GetXaxis()->GetBinCenter(1);
  double phij = LLVPBaseLineHIst->GetYaxis()->GetBinCenter(1)*TMath::Pi()/180.0;

  //-----------------------------------------------------------------------------------------------------------------------------------------

    double PileHeight = 1000; // Height of LLVP in km
    double PileRadius = Rcmb + PileHeight; //km
    double DepthMin = Rcmb;
    double DepthMax = PileRadius; // Distance from the center of the Earth
    double PileDensityPct = 30.0; // 2% density contrats for LLVP

  double zenmin = 180-TMath::ASin( DepthMax/Rocean )*(180.0/TMath::Pi()) ; // min 90
  double zenmax = 180-TMath::ASin( (DepthMin)/Rocean )*(180.0/TMath::Pi()) ; // max 180

  std::cout << zenmin << " " << zenmax << std::endl;

  double eta = 35;
  double th = 180-eta;

  //double cth = cos(th*TMath::Pi()/180.0);

  //double cth = -0.8194; //Mid point lower center

  //double cth = -0.95;

  double cth = cos(th*TMath::Pi()/180.0);

  //double cth =  -0.75;

  double phi = 0*TMath::Pi()/180.0;

  // OSCPROB

  std::string PremFolder = "/home/dehy0499/OscProb/PremTables/";
  std::string PremName = "prem_44layers.txt";
  std::string path2prem = PremFolder+PremName;

  OscProb::PremModel prem(path2prem);

  prem.FillPath(cth); // Fill paths from PREM model

  std::vector<OscProb::NuPath> paths = prem.GetNuPath();




  //Mine


  Earth3DModel test;

  test.PileThickness = 1000;

  test.aWidth = 45;
  test.SetModel(path2prem);

  test.SetDirection(cth, phi);
  
  test.SetPile( true, "cake", 3.0, 0.0);

  test.SetLayerProp(27,0.0,0.0);

  std::vector<std::vector<double>> EarthPath = test.Create3DPath( );

  std::cout << "path Seems good" << std::endl;

  //test.SetPile( true, "cake", 0.0, 0.0);

  //std::vector<std::vector<double>> EarthPathPREM = test.Create3DPath( );

  //std::cout << "Mine" << std::endl;

  TApplication app("app", &argc, argv);

  std::cout << "Before Canvas" << std::endl;
  TCanvas * c = test.EarthCanvas;
  std::cout << "After Canvas" << std::endl;

  //c->Print("./TEST.png");

  c->Modified(); c->Update();
  
  std::cout << "Not here " << std::endl;

  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  
  std::cout << "Not here " << std::endl;



  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  std::cout << "Not here " << std::endl;
  
  app.Run();


/*
  //Caculate Baselines

  // Get the current date and time
  std::time_t now = std::time(0);
  std::tm* localTime = std::localtime(&now);

  // Format the date (YYYY-MM-DD)
  char dateBuffer[11];
  std::strftime(dateBuffer, sizeof(dateBuffer), "%Y-%m-%d", localTime);


  std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
  
  std::string DistanceFolder = "/SimulationResults/PreliminaryResults/BaselinesLLVPs/";

  std::string BaselineFile = "BaselineLLVP_"+llvpShape+"_"+std::to_string(zenithBins)+"zenithBins"+"_"+ std::to_string(static_cast<int>(zenithMin)) +"-"+ 
                                 std::to_string(static_cast<int>(zenithMax))+"Zen_"+ std::to_string(azimuthBins)+ std::to_string(static_cast<int>(azimuthMin)) +"-"+ 
                                 std::to_string(static_cast<int>(azimuthMax))+"Az_"+std::string(dateBuffer)+".csv";



  double Ltot_up = 0;
  double Lij = 0;

  for (int i = 1; i <= LLVPBaseLineHIst->GetNbinsX(); ++i)
  {
        for (int j = 1; j <= LLVPBaseLineHIst->GetNbinsY(); ++j)
        {
            double cthi = LLVPBaseLineHIst->GetXaxis()->GetBinCenter(i);
            double phij = LLVPBaseLineHIst->GetYaxis()->GetBinCenter(j)*TMath::Pi()/180.0;

            Baselines.SetModel(path2prem);

            Baselines.SetDirection(cthi, phij);
  
            Baselines.SetPile( true, llvpShape, 3.0, 0.0);

            std::vector<std::vector<double>> EarthPathAlt = Baselines.Create3DPath( );

            Baselines.SetPile( true, llvpShape, 0.0, 0.0);

            std::vector<std::vector<double>> EarthPathStd = Baselines.Create3DPath( );


            for (int i = 0; i < EarthPathStd.size() ; ++i)
            {

              Ltot_up += EarthPathStd[i][0]; 
              double Ltot_low = Ltot_up - EarthPathStd[i][0];  

              std::cout << "******Safe Check " <<EarthPathStd[i][0] << " " << Ltot_up << " " << Ltot_low << " " << (EarthPathAlt[i][1]-EarthPathStd[i][1]) << std::endl;


              if ( 100*abs(EarthPathAlt[i][1]-EarthPathStd[i][1])/EarthPathStd[i][1] > 0.01 )
              {
                Lij += (Ltot_up-Ltot_low);


              }

            }

          LLVPBaseLineHIst->SetBinContent(i, j, Lij);
          Lij=0;
          Ltot_up=0;

        }
    }

  ExportToCSV(LLVPBaseLineHIst, NuTomoPath+DistanceFolder+BaselineFile);
*/
 

return 0;

}


void ExportToCSV(TH2D* hist, std::string filename)
{


    std::ofstream outfile(filename);

    // Recorre los bins y guarda el contenido


    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= hist->GetNbinsY(); ++j)
        {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);

            double N = hist->GetBinContent(i, j);
     

            outfile << x << "," << y << "," << N << std::endl;
        }
    }

    outfile.close();
}