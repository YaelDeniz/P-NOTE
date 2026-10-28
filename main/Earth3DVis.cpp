#include "Earth3DModel.h"

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

#include "TGraph.h"
#include "TMultiGraph.h"

#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TF1.h"


//OSCPROB

#include "PremModel.h"
#include "LayerUpdate.h"
#include "PhyTools.h"

int main(int argc, char **argv)
{

    double PileHeight = 1000; // Height of LLVP in km
    double PileRadius = R_CMB + PileHeight; //km
    double DepthMin = R_CMB;
    double DepthMax = PileRadius; // Distance from the center of the Earth
    double PileDensityPct = 3.0; // 2% density contrats for LLVP

  double zenmin = 180-TMath::ASin( DepthMax/R_OCEAN )*(180.0/TMath::Pi()) ; // min 90
  double zenmax = 180-TMath::ASin( (DepthMin)/R_OCEAN )*(180.0/TMath::Pi()) ; // max 180

  std::cout << zenmin << " " << zenmax << std::endl;

  double th = 179.9;

  //double cth = cos(th*TMath::Pi()/180.0);

  //double cth = -0.8194; //Mid point lower center

  //double cth = -0.95;

  double cth = cos(th*TMath::Pi()/180.0);

  //double cth =  -0.75;

  double phi = (90)*TMath::Pi()/180.0;

  // OSCPROB

  const char* osc_env = getenv("OSCPROB_PATH");
  std::string PremFolder = std::string(osc_env ? osc_env : (std::string(getenv("HOME")) + "/OscProb")) + "/PremTables/";
  std::string PremName = "prem_default.txt";
  std::string path2prem = PremFolder+PremName;

  OscProb::PremModel prem(path2prem);
  prem.SetDetPos(6368.0);
  prem.FillPath(cth); // Fill paths from PREM model

  std::vector<OscProb::NuPath> paths = prem.GetNuPath();

  double sumL = 0;

  if (th < 90.2 && th > 89.8)
  {
    

    for (int i = 0; i < paths.size() ; ++i)
    {

    sumL += paths[i].length ; 
    std::cout << i << " " << sumL << " " << paths[i].length   << " " <<  paths[i].density  << " " << paths[i].zoa << std::endl;

    }
   
   std::cout << "It work" << std::endl;

  }




  //Mine

  TApplication app("app", &argc, argv);

  // =========================================================================
  // Earth3DModel — available commands reference
  // =========================================================================
  //
  // DISPLAY
  //   test.ModelVisibility = true/false;
  //       Show or hide the 3D OpenGL canvas.
  //
  // PREM TABLE
  //   test.SetModel(std::string path2table);
  //       Path to the PREM .txt file (e.g. prem_default.txt or prem_44layers.txt).
  //
  // NEUTRINO DIRECTION
  //   test.SetDirection(double cth, double phi);
  //       cth = cos(zenith angle), must be in [-1, 1].
  //       phi = azimuthal angle in radians.
  //
  // DETECTOR POSITION
  //   double pos[3] = {x, y, z};  test.SetDetector(pos);
  //       Cartesian position in km. Default = KM3NeT site (0, 0, -6368).
  //
  // LLVP (Large Low-Velocity Province) — direct member access + method
  //   test.PileThickness = 1000;      // LLVP height above CMB in km
  //   test.aWidth       = 25.0;       // angular half-width in degrees
  //   test.SetPile(bool on, std::string shape, double densityPct, double zoaPct);
  //       shape: "pancake" (single shell) or "cake" (three-tiered, narrowing toward top)
  //       densityPct: density contrast as percentage (e.g. 3.0 = +3%)
  //       zoaPct:     Z/A (chemical composition) contrast as percentage
  //
  // SINGLE-LAYER MODIFICATION (legacy — use batch below when modifying multiple layers)
  //   test.SetLayerProp(int row, double densityPct, double zoaPct);
  //       row is 1-based PREM table row index.
  //
  // BATCH LAYER MODIFICATION
  //   std::vector<LayerUpdate> updates;
  //   LayerUpdate mod; mod.tableNumber = N; mod.densityContrast = X; mod.chemicalContrast = Y;
  //   updates.push_back(mod);  // repeat for each layer
  //   test.AddUpdatesList(updates);
  //
  // PATH COMPUTATION
  //   std::vector<std::vector<double>> path = test.Create3DPath();
  //       Each inner vector: { Li (km), density (g/cm³), Z/A, Ri (km) }
  //       Call after all configuration is set. Also populates test.EarthCanvas.
  //
  // CANVAS
  //   TCanvas* c = test.EarthCanvas;
  //       Valid after Create3DPath(). Shows the 3D geometry with the neutrino path.
  //
  // =========================================================================

  Earth3DModel test;
  test.ModelVisibility = true;

  test.PileThickness = 1000;
  test.aWidth = 45;
  test.SetModel(path2prem);

  test.SetDirection(cth, phi);

  test.SetPile( true, "pancake", 3.0, 0.0);
  test.SetLayerProp(27, 0.0, 0.0);  // single-layer override: row 27, no density/chem contrast

  std::vector<LayerUpdate> Updates_;
  LayerUpdate mod;

  mod.tableNumber=44;      // 1-based PREM CSV row
  mod.densityContrast = 1000000 ;  // relative change in density (%)
  mod. chemicalContrast = 0 ; // relative change in Z/A (%)
  Updates_.push_back(mod);

  mod.tableNumber=14;      // 1-based PREM CSV row
  mod.densityContrast = 100 ;  // relative change in density (%)
  mod. chemicalContrast = 0 ; // relative change in Z/A (%)
  Updates_.push_back(mod);

   mod.tableNumber=16;      // 1-based PREM CSV row
  mod.densityContrast = 100 ;  // relative change in density (%)
  mod. chemicalContrast = 0 ; // relative change in Z/A (%)
  Updates_.push_back(mod);


   mod.tableNumber=18;      // 1-based PREM CSV row
  mod.densityContrast = 100 ;  // relative change in density (%)
  mod. chemicalContrast = 0 ; // relative change in Z/A (%)
  Updates_.push_back(mod);

  test.AddUpdatesList(Updates_);

  

  
  std::vector<std::vector<double>> EarthPath = test.Create3DPath( );


  std::cout << "Mine" << std::endl;


  sumL = 0;

  TMultiGraph * ThePaths = new TMultiGraph();

  TGraph * MyModel = new TGraph(EarthPath.size());

  TGraph * apcModel = new TGraph(paths.size());


 
  std::ofstream ModelPath("ModelPath.csv"); 

  for (int i = 0; i < EarthPath.size() ; ++i)
  {

    sumL += EarthPath[i][0]; 
    std::cout<< i << "  me " << sumL<< " " << std::setprecision(10) <<EarthPath[i][0] << " " <<EarthPath[i][1]  << " " <<EarthPath[i][2]<< std::endl;

    ModelPath << sumL << " " << EarthPath[i][0] << " " << EarthPath[i][1]<< std::endl;

    MyModel->SetPoint(i,sumL,EarthPath[i][1] );

  }

  std::cout << "     " << std::endl;
  std::cout << "OscProb" << std::endl;

  sumL = 0;


  std::ofstream OscProbPath("OscProbPath.csv"); 
  
  for (int i = 0; i < paths.size() ; ++i)
  {

    sumL += paths[i].length ; 
    std::cout << i << " " << sumL << " " << paths[i].length - EarthPath[i][0]  << " " <<  paths[i].density  << " " << paths[i].zoa << std::endl;

    OscProbPath << sumL << " " << paths[i].length << " " << paths[i].density << std::endl;

    apcModel->SetPoint(i,sumL,paths[i].density  );

  }

  OscProbPath.close();
  ModelPath.close();

  ThePaths->Add(MyModel);
  ThePaths->Add(apcModel);

 std::cout << zenmin << " " << zenmax <<  std::endl;

  
  TCanvas * c = test.EarthCanvas;

  //c->Print("./TEST.png");

  c->Modified(); c->Update();

  
  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");


  TCanvas * c1 = new TCanvas();

  ThePaths->Draw("apl");

  c1->Modified(); c1->Update();

  
  TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();
  rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");



  app.Run();
  

return 0;

}
