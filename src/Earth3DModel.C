/*
  
  Description: 
  Class generates a 3D Sphercial model of the Earth based of PREM density pro-
  files.Model Includes the presence of a LLVP-like object in the Lowermost ma-
  ntle.Class is used to generate the Neutrino Path at specific regions og the
  Earth.

  Author: Deniz, Yael.
  Earth and Spatial Science Deparment, University of Idaho.

*/

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




//Cern Root
#include <math.h>
#include "TMath.h"
#include <Math/SVector.h>
#include "TCanvas.h"
#include "TSystem.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TView.h"


//Geometry manager

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoVolume.h"
#include "TGLViewer.h"


using namespace std;

# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)


//Set the Prem Model

void Earth3DModel::SetModel(std::string path2table)
{
    premtable = path2table; //.txt files from PREM Tables
} 
  
//Set the Detector Location

void Earth3DModel::SetDetector(double Position[3]) //Detetctor in Cartesian coordinates
{
   Det[0] = Position[0]; // x coordinate
   Det[1] = Position[1]; // y coordinate
   Det[2] = Position[2]; // z coordinate
}



void Earth3DModel::SetDirection(double cth , double phi ) //Set Neutrino Direction, [phi] = rads
{
    
    if (cth < -1 || cth > 1) 
    {
        throw std::invalid_argument("Invalid value as neutrino direction, zenith direction must be in terms of cos(zen)");
    }

    else
    {

    zenith = TMath::Pi() - acos(cth); //Polar(zenith) angle in spherical coordinates

    azimuth = phi;   //azimut angle in spherical coordinates
    
    //std::cout <<" Neutrino direction: " <<  zenith*180.0/TMath::Pi() << " " << azimuth*180.0/TMath::Pi() << std::endl;

    }
} 

//Define if LLVP exist

void Earth3DModel::SetPile( bool value , std::string shape, double density, double zoa){ Pile = value; // Presence of Pile
  PileShape = shape; //Pile shape: "Slab" or "Segmented"
  PileDensityContrats = density;
  PileZoAContrats = zoa;} 


// Create a matrix from PremTables
std::vector< std::vector<double> > Earth3DModel::GetPremData( std::string path2table ){

   std::vector< std::vector<double> > PremMatrix; // Matrix data  form of Prem tables

   std::vector<double> PremRow(4);               // PremRow = [radius density z/a, layer ID]

   double radius, density, zoa, layer; // Variables for storing Prem table values

   std::ifstream PREM_DATA;

   PREM_DATA.open(path2table);

      // Loop over table rows
   while(PREM_DATA >> radius >> density >> zoa >> layer)
   {

      PremRow[0] = radius;
      PremRow[1] = density;
      PremRow[2] = zoa;
      PremRow[3] = layer;

      PremMatrix.push_back(PremRow);
    }

   return PremMatrix;}


// PATCH 12/2/25: Modification of multiple layers

void Earth3DModel::AddUpdatesList(const std::vector<LayerUpdate>& updates)
{
    // AddUpdatesList
    // --------------
    // Registers a batch of layer property updates for the PREM-based Earth model.
    //
    // Parameters
    // ----------
    //   updates : vector of LayerUpdate objects. Each entry corresponds to a
    //             PREM layer (1-based row index in the CSV) and contains:
    //               - tableNumber      : 1-based index of the PREM layer
    //               - densityContrast  : relative change in density (%)
    //               - chemicalContrast : relative change in Z/A (%)
    //
    // For each element, this method forwards the values to ChangeLayerProperties(),
    // which stores or updates the corresponding entry in UpdatesList.
    // All stored modifications are later applied in UpdateLayers().

    if (updates.empty()) {
        Error("AddUpdatesList",
              "No layer updates were provided. At least one LayerUpdate "
              "is required to modify the PREM model.");
        return;
    }

    for (const auto& upd : updates) {
        ChangeLayerProperties(
            upd.tableNumber,      // 1-based index, same as CSV row
            upd.densityContrast,  // %
            upd.chemicalContrast  // %
        );
    }
}


void Earth3DModel::ChangeLayerProperties(int TableNumber,
                                         double DensityContrast,
                                         double ChemicalContrast)
{
    // ChangeLayerProperties
    // ---------------------
    // Registers or updates a modification for a given PREM layer.
    //
    // Parameters:
    //   TableNumber      : 1-based index of the layer (row) in the PREM CSV file
    //   DensityContrast  : relative change in density (%) for this layer
    //   ChemicalContrast : relative change in Z/A (%) for this layer
    //
    // This function keeps TableNumber as 1-based in UpdatesList to match the CSV.
    // Conversion to 0-based happens only when indexing EarthMatrix in UpdateLayers().

    const int tableNumber = TableNumber;  // keep 1-based

    // If this layer already has a modification, update it instead of duplicating
    for (auto& mod : UpdatesList) {
        if (mod.tableNumber == tableNumber) {
            mod.densityContrast  = DensityContrast;
            mod.chemicalContrast = ChemicalContrast;
            return;
        }
    }

    // Otherwise add a new modification
    LayerUpdate newUpt;
    newUpt.tableNumber      = tableNumber;      // 1-based
    newUpt.densityContrast  = DensityContrast;
    newUpt.chemicalContrast = ChemicalContrast;

    UpdatesList.push_back(newUpt);
}


std::vector<std::vector<double>>
Earth3DModel::UpdateLayers(std::vector<std::vector<double>> EarthMatrix)
{
    // UpdateLayers
    // ------------
    // Applies all user-defined layer modifications to the PREM-based Earth model.
    //
    // Each modification is stored in the LayerUpdate struct and specifies:
    //   - tableNumber      : 1-based index of the PREM layer (CSV row)
    //   - densityContrast  : relative change in density (%)
    //   - chemicalContrast : relative change in Z/A (%)
    //
    // For each valid row, updates:
    //   EarthMatrix[row][1] ← density scaled by (1 + densityContrast  / 100)
    //   EarthMatrix[row][2] ← Z/A    scaled by (1 + chemicalContrast / 100)
    //
    // It is intended to be called inside Earth3DPath, after loading the PREM table
    // and after all desired calls to ChangeLayerProperties/AddUpdatesList.

    for (const auto& upd : UpdatesList) {
        // Convert 1-based table number to 0-based row index
        int row = upd.tableNumber - 1;

        // Safety check: skip if out of range
        if (row < 0 || row >= static_cast<int>(EarthMatrix.size()))
            continue;

        // Assuming:
        // EarthMatrix[row][0] = radius
        // EarthMatrix[row][1] = density
        // EarthMatrix[row][2] = Z/A

        // Apply density contrast (%)
        EarthMatrix[row][1] *= (1.0 + upd.densityContrast / 100.0);

        // Apply Z/A (chemical) contrast (%)
        EarthMatrix[row][2] *= (1.0 + upd.chemicalContrast / 100.0);
    }

    return EarthMatrix;
}



//

//Modify Specific Layers of the Earth

void Earth3DModel::SetLayerProp(int TableNumber, double DensityContrats, double ChemicalContrats) { 
   //Table number is the row that correspond to the desire layer in the .txt files


  PremRow     = TableNumber-1;
  LayerDensityContrats   = DensityContrats;  
  LayerZoAContrats   = ChemicalContrats;}

 std::vector< std::vector<double> > Earth3DModel::ChangeLayerProp(std::vector< std::vector<double> > EarthMatrix){
 
  EarthMatrix[PremRow][1] = (1.0 + LayerDensityContrats/100.0) * EarthMatrix[PremRow][1];

  double test = (1.0 + LayerDensityContrats/100.0) * EarthMatrix[PremRow][1];

  EarthMatrix[PremRow][2] = (1.0 + LayerZoAContrats/100.0) * EarthMatrix[PremRow][2];

  //std::cout << " ******************************************************************************+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+* " <<  100.0*(test  - EarthMatrix[PremRow][1])/EarthMatrix[PremRow][1] << std::endl;

  return EarthMatrix;}

//CONSTRUCT LLVP

void Earth3DModel::CreatePanCake(std::vector<TGeoVolume*> LAYER, std::vector< std::vector<double> > PremMatrix ){

  //Baker----------------------------------------------------------------------------------------------------------------------
  
  std::vector<TGeoMaterial*> LLVPMat; //Vector to define material for each LLVP segment ( A X% anomaly in local Density or Z/A)

  std::vector<TGeoMedium*> LLVPMed;   //Vector to define medium for each LLVP segment

  std::vector<TGeoVolume*> LLVPLAYER; //Vector to define Volume for each LLVP segment
  
  TGeoRotation   *rot1 = new TGeoRotation("rot1", 0.0, 270.0, 0.0);// Some Geometrical trasformation that move LLVPs to the right place 

  //TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);// Some Geometrical trasformation that move LLVPs to the right place 




  //double PileThickness = 1000;

  double Bottom,Top;

  double Localdensity;

  double Localchem;

  Bottom = Rcmb;  //Bottom of LLVP is the cmb
  
  Top = Rcmb + PileThickness; //Top of the LLVP 

   // Segment 1 test
   double slabLength =0;

   double InnerR = 0;
   double OuterR = 0;

   int PileId = 0;

   //std::cout << " " << std::endl;
   //std::cout << "Pile radius: " << Bottom << " - " << Top << std::endl;

   double lowerSegment = Bottom; //Lower Radius of a segment in Pile section

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInPile = PremMatrix[i][0]; //Lower boundary of current layer

      //std::cout << PremMatrix.size() << " Current Layer Properties " << PremMatrix[i][0] << " " << PremMatrix[i][1] << " " << PremMatrix[i][2] << " " << PremMatrix[i][3] << " " << Top <<std::endl;

      if (RadiusInPile > Bottom )
      {

          char integer_string[32];
            
          sprintf(integer_string, "%d", i);//Layer Index

          char llvpmaterial_string[64]="LLVPMaterial";
          strcat(llvpmaterial_string, integer_string);
          const char *llvpMatName = llvpmaterial_string;
          

          char llvpmedium_string[64]="LLVPMedium";
          strcat(llvpmedium_string, integer_string);
          const char *llvpMedName = llvpmedium_string;
            
          char llvplayer_string[64]="LLVPLayer";
          strcat(llvplayer_string, integer_string);
          const char *llvpLayerName= llvplayer_string;
        
         slabLength = slabLength + (RadiusInPile - lowerSegment ); //Calculate segement lenght 

         if (slabLength > PileThickness)
         { 

            double excess = PileThickness - slabLength;
            InnerR= lowerSegment;
            OuterR= RadiusInPile + excess;
   

            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            ////std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            ////std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
            Localdensity = (1.0 + PileDensityContrats/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + PileZoAContrats/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(llvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(llvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            //std::cout<< " Last " << InnerR << " " << OuterR << " " << Localdensity << " " << PremMatrix[i][1] << std::endl; 


            LLVPLAYER.push_back( gGeoManager -> MakeSphere(llvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, aWidth ,0 ,360 ) ); // DEFINE LLVP Segment

            //LLVPLAYER[PileId]->SetLineColor(kBlue);
            LLVPLAYER[PileId]->SetLineColor(46);

            LAYER[i]->AddNodeOverlap(LLVPLAYER[PileId],1,rot1);  

            LLVPMat.clear();
            
            LLVPMed.clear();
            
            LLVPLAYER.clear();

            break;
         }
    
         else
         {
            InnerR = lowerSegment;
            OuterR = RadiusInPile;
            
            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            ////std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            ////std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
          
            ////std::cout << "Radius of lower section's segements: "<< InnerR << " " << OuterR << std::endl;

            Localdensity = (1.0 + PileDensityContrats/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + PileZoAContrats/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(llvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(llvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            //std::cout << "Still a layer" << InnerR << " " << OuterR << " " << Localdensity << " " << PremMatrix[i][1] << std::endl; 

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(llvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, aWidth ,0 ,360 ) ); // DEFINE LLVP Segment

            //LLVPLAYER[PileId]->SetLineColor(kBlue);
            LLVPLAYER[PileId]->SetLineColor(46);

            LAYER[i]->AddNodeOverlap(LLVPLAYER[PileId],1,rot1);  

            ++PileId;

         }

         lowerSegment = RadiusInPile;
      }
   }}

void Earth3DModel::CreateCake(std::vector<TGeoVolume*> LAYER, std::vector< std::vector<double> > PremMatrix ){

  //Baker----------------------------------------------------------------------------------------------------------------------
  
  std::vector<TGeoMaterial*> LLVPMat; //Vector to define material for each LLVP segment ( A X% anomaly in local Density or Z/A)

  std::vector<TGeoMedium*> LLVPMed;   //Vector to define medium for each LLVP segment

  std::vector<TGeoVolume*> LLVPLAYER; //Vector to define Volume for each LLVP segment
  
  TGeoRotation   *rot1 = new TGeoRotation("rot1", 0.0, 270.0, 0.0);// Some Geometrical trasformation that move LLVPs to the right place

   //double PileThickness = 1000;

   double thicknessOfSection = PileThickness/3.0; // Salabas have same thickness

   double lowerBottom, lowerTop, midBottom, midTop, upperBottom, upperTop;

   double Localdensity;

   double Localchem;

   lowerBottom = Rcmb;

   lowerTop = Rcmb + thicknessOfSection;

   midTop = Rcmb + 2*thicknessOfSection;

   upperTop = Rcmb + 3*thicknessOfSection;

   // Segment 1 test
   double slabLength =0;

   double InnerR = 0;
   double OuterR = 0;

   int PileId = 0;

   //std::cout << " " << std::endl;
   //std::cout << " Segement 1 radius: " << lowerBottom << " - " << lowerTop << std::endl;

   double lowerSegment = lowerBottom; //Lower Radius of a segment in Pile section

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInBottom = PremMatrix[i][0]; //Lower boundary of current layer

      if (RadiusInBottom > lowerBottom )
      {

          char lowinteger_string[32];
            
          sprintf(lowinteger_string, "%d", i);//Layer Index

          char lowllvpmaterial_string[64]="lowLLVPMaterial";
          strcat(lowllvpmaterial_string, lowinteger_string);
          const char *lowllvpMatName = lowllvpmaterial_string;
          

          char lowllvpmedium_string[64]="lowLLVPMedium";
          strcat(lowllvpmedium_string, lowinteger_string);
          const char *lowllvpMedName = lowllvpmedium_string;
            
          char lowllvplayer_string[64]="lowLLVPLayer";
          strcat(lowllvplayer_string, lowinteger_string);
          const char *lowllvpLayerName= lowllvplayer_string;
        
         slabLength = slabLength + (RadiusInBottom - lowerSegment ); //Calculate segement lenght 

         if (slabLength > thicknessOfSection)
         { 

            double excess = thicknessOfSection - slabLength;
            InnerR= lowerSegment;
            OuterR= RadiusInBottom + excess;
  
          ////std::cout << "top of lower segment: " << InnerR << " " << OuterR << std::endl;
 
            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            ////std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            ////std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
            
            Localdensity = (1.0 + PileDensityContrats/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + PileZoAContrats/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(lowllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(lowllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(lowllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, aWidth ,0 ,360 ) ); // DEFINE LLVP Segment

            //LLVPLAYER[PileId]->SetLineColor(kBlue);
            LLVPLAYER[PileId]->SetLineColor(46);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            LLVPMat.clear();
            
            LLVPMed.clear();
            
            LLVPLAYER.clear();

            midBottom = OuterR;
            
            break;
         }
    
         else
         {
            InnerR = lowerSegment;
            OuterR = RadiusInBottom;
            
            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            ////std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            ////std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
          
            ////std::cout << "Radius of lower section's segements: "<< InnerR << " " << OuterR << std::endl;

            Localdensity = (1.0 + PileDensityContrats/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + PileZoAContrats/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(lowllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(lowllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(lowllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, aWidth ,0 ,360 ) ); // DEFINE LLVP Segment

            //LLVPLAYER[PileId]->SetLineColor(kBlue);
            LLVPLAYER[PileId]->SetLineColor(46);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            ++PileId;

         }

         lowerSegment = RadiusInBottom;
      }
   }

  
   //Segment 2------------------------------------------------------------------------------------------------------------------------
   PileId=0;

   slabLength =0;
   InnerR = 0;
   OuterR = 0;
   
   //std::cout << " " << std::endl;
   //std::cout << " Segement 2 radius: " << midBottom << " - " << midTop << std::endl;
   
   double midSegment = midBottom; //Lower Radius of a segment in Pile section

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInMiddle = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (RadiusInMiddle > midBottom )
      {
        char midinteger_string[32];
            
        sprintf(midinteger_string, "%d", i);//Layer Index

        char midllvpmaterial_string[64]="midLLVPMaterial";
        strcat(midllvpmaterial_string, midinteger_string);
        const char *midllvpMatName = midllvpmaterial_string;
          

        char midllvpmedium_string[64]="midLLVPMedium";
        strcat(midllvpmedium_string, midinteger_string);
        const char *midllvpMedName = midllvpmedium_string;
            
        char midllvplayer_string[64]="midLLVPLayer";
        strcat(midllvplayer_string, midinteger_string);
        const char *midllvpLayerName= midllvplayer_string;

         //std::cout << " Radius of layer above middle bottom  " << RadiusInMiddle << std::endl;


         slabLength = slabLength + (RadiusInMiddle - midSegment ); //Calculate segment lenght 

         if (slabLength > thicknessOfSection)
         {
            double midExcess = thicknessOfSection-slabLength;
            
            InnerR= midSegment;

            OuterR= RadiusInMiddle + midExcess;
            
            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();

            Localdensity = (1.0 + PileDensityContrats/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + PileZoAContrats/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(midllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(midllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            double midwidth= (2.0/3.0)*aWidth;

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(midllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, midwidth ,0 ,360 ) ); // DEFINE LLVP Segment

            //LLVPLAYER[PileId]->SetLineColor(kGreen);
            LLVPLAYER[PileId]->SetLineColor(46);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            LLVPMat.clear();
            
            LLVPMed.clear();
            
            LLVPLAYER.clear();

            upperBottom = OuterR;

            break;
         }

         else
         {

          InnerR = midSegment;

          OuterR = RadiusInMiddle;

          double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();

          double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
       
          Localdensity = (1.0 + PileDensityContrats/100.0)*PremMatrix[i][1];

          Localchem = (1.0 + PileZoAContrats/100.0 )*PremMatrix[i][2];
  
          LLVPMat.push_back( new TGeoMaterial(midllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material

          LLVPMed.push_back( new TGeoMedium(midllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

          double midwidth= (2.0/3.0)*aWidth;

          LLVPLAYER.push_back( gGeoManager -> MakeSphere(midllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, midwidth ,0 ,360 ) ); // DEFINE LLVP Segment

          //LLVPLAYER[PileId]->SetLineColor(kGreen);
          LLVPLAYER[PileId]->SetLineColor(46);

          LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1); 

          ++PileId; 

         }
         midSegment = RadiusInMiddle;
      }

   }
   //Segment 3------------------------------------------------------------------------------------------------------------------------
   PileId=0;
   slabLength =0;
   InnerR = 0;
   OuterR = 0;
   
   //std::cout << " " << std::endl;
   
   //std::cout << " Segement 3 radius: " << upperBottom << " - " << upperTop << std::endl;

   double upperSegment = upperBottom;

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInTop = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (RadiusInTop > upperBottom )
      {


        char upinteger_string[32];
            
        sprintf(upinteger_string, "%d", i);//Layer Index

        char upllvpmaterial_string[64]="upLLVPMaterial";
        strcat(upllvpmaterial_string, upinteger_string);
        const char *upllvpMatName = upllvpmaterial_string;
          

        char upllvpmedium_string[64]="upLLVPMedium";
        strcat(upllvpmedium_string, upinteger_string);
        const char *upllvpMedName = upllvpmedium_string;
            
        char upllvplayer_string[64]="upLLVPLayer";
        strcat(upllvplayer_string, upinteger_string);
        const char *upllvpLayerName= upllvplayer_string;

         //std::cout << " Radius of layer above uppper bottom  " << RadiusInTop << std::endl;

         slabLength = slabLength + (RadiusInTop - upperSegment );

         //r.push_back(dh-mslabL);

         if (slabLength > thicknessOfSection)
         {
            double topExcess = thicknessOfSection-slabLength;

            InnerR= upperSegment;

            OuterR= RadiusInTop + topExcess;

            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();

            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            
            Localdensity = (1.0 + PileDensityContrats/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + PileZoAContrats/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(upllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(upllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            double upwidth= (1.0/3.0)*aWidth;

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(upllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, upwidth ,0 ,360 ) ); // DEFINE LLVP Segment

            //LLVPLAYER[PileId]->SetLineColor(kRed);
            LLVPLAYER[PileId]->SetLineColor(46);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            LLVPMat.clear();
            
            LLVPMed.clear();
            
            LLVPLAYER.clear();

            break;
         }

         else
         {
            InnerR = upperSegment;
           
            OuterR = RadiusInTop;

            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
           
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
           
            Localdensity = (1.0 + PileDensityContrats/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + PileZoAContrats/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(upllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(upllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            double upwidth= (1.0/3.0)*aWidth;

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(upllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, upwidth ,0 ,360 ) ); // DEFINE LLVP Segment

            //LLVPLAYER[PileId]->SetLineColor(kRed);
            LLVPLAYER[PileId]->SetLineColor(46);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            ++PileId;

         }
         upperSegment = RadiusInTop;
      }

   }}

// CONSTRUCT THE 3D MODEL

std::vector<std::vector<double>> Earth3DModel::Earth3DPath( double th, double phi, std::string path2table ){

   gSystem->Load("libGeom");

   std::vector< std::vector<double> > PremMatrix = GetPremData(path2table); // Sort PREM model into a readable matrix

   //PremMatrix = ChangeLayerProp(PremMatrix); //Single layer modifiction
   PremMatrix = UpdateLayers(PremMatrix);  //Multiple layer modification

   //std::cout << " " << std::endl;

   // Section to check Prem Matrix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

       //std::cout<< i << " " << PremMatrix[i][0] << " " << PremMatrix[i][1] << " " <<PremMatrix[i][2] << " " <<PremMatrix[i][3] << std::endl;

   }
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  std::vector< std::vector<double> > LLVPMatrix = GetPremData(path2table); // A Copy of PREM model matrix to construct LLVPs

  //double Ratm = PremMatrix.back()[0]; // Outermost Layer in the PREM model

  std::vector< std::vector<double> > EarthPath;  // Store Paths inside the Earth for oscillation calculations

  //std::vector< int > EarthPathId;  // Store ID layer for each Paths inside the Earth for oscillation calculations

  EarthCanvas = new TCanvas("3D Earth", "3D Earth",0,0,600,600);
 
  new TGeoManager("EarthModel3D", "Spherical Earth ");

  gGeoManager->SetVerboseLevel(0);

  //Defining Vacuum Medium
  TGeoMaterial *VAC = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);
  TGeoMedium *vac = new TGeoMedium("VACUUM",99,VAC);
  TGeoVolume *top = gGeoManager->MakeBox("TOP",vac,Ratm+5000,Ratm+5000,Ratm+5000); // Container Volume

  gGeoManager->SetTopVolume(top);

  //Create Earth

  // Each Layer is represented of a Spherical Shell defined by a Inner Radius (rmin) and Outer Radius(rmax).
  // Each Layer Medium is defined by a Material. Material is created from a specific density value and Z[proton]/A[neutron + proton] (zoa)

  double rmin, rmax, density, zoa, layer; 

  std::vector<TGeoMaterial*> MAT; //Vector to define material for each layer en in PREM (Density & Z/A)
  std::vector<TGeoMedium*>  MED;  //Vector to define medium for each layer en in PREM 
  std::vector<TGeoVolume*> LAYER; //Vector to define a  volume for each layer en in PREM
  std::vector<std::vector<double>> check;

  for (int i = 0; i < PremMatrix.size(); ++i)
   {
      if (i==0) { rmin = 0; }

      else { rmin = PremMatrix[i-1][0]; }

      rmax = PremMatrix[i][0]; //Layer Outer radius from PREM 

      density = PremMatrix[i][1];  //Layer density from PREM
      
      zoa = PremMatrix[i][2]; //Layer Z/A from PREM
      
      layer = PremMatrix[i][3]; //Layer ID number from PREM

      //Generate Names for Material, Medium and Layer
      char integer_string[32];

      sprintf(integer_string, "%d", i+1);

      char material_string[64]="Material";
      strcat(material_string, integer_string);
      const char *MatName = material_string;

      char medium_string[64]="Medium";
      strcat(medium_string, integer_string);
      const char *MedName = medium_string;

      char layer_string[64]="Layer";
      strcat(layer_string, integer_string);
      const char *LayerName= layer_string;

      //Definiton of world "Daugther" volumes (i.e., Earth Layers)

      MAT.push_back( new TGeoMaterial( MatName, 1.0 , zoa , density) ); //Define Material for ith Layer 
       
      MED.push_back( new TGeoMedium( MedName , 1 , MAT[i]) ); //Define Medium for ith Layer 

      check.push_back({rmin,rmax});

      LAYER.push_back( gGeoManager->MakeSphere(LayerName,MED[i], rmin,rmax,0,180,0,360) ); // Define Volume for the ith Layer

      if (rmax == 3480.0 ) { 
      LAYER[i]->SetVisibility(kTRUE);
      LAYER[i]->SetTransparency(0);
      LAYER[i]->SetLineColor(28);


      } //Outer Core is visible

      else { LAYER[i]->SetVisibility(kFALSE); }
  }

      LAYER[LAYER.size()-1]-> SetVisibility(kTRUE); // Crust is vissible
      LAYER[LAYER.size()-1]-> SetTransparency(95); // Crust is vissible

  
  /*
  //SECTION TO BE DELETED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double R_cmb = 3480.0;
  double thtest = TMath::ASin( (R_cmb)/6368.0 ) ; // max 180

  bool TestBoolLLVP = Pile;

  //std::cout << " Test of angles: " << th*180.0/TMath::Pi() << " " << thtest*180.0/TMath::Pi() << " " << Pile << " " << TestBoolLLVP <<std::endl;
  if (th  <  thtest )
  {
    TestBoolLLVP = false;
    //std::cout << "LLVP IS OFF " << th*180.0/TMath::Pi() << " " << thtest*180.0/TMath::Pi() << " " << Pile << " " << TestBoolLLVP <<std::endl;
  }
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    */

  //Create LLVP segments 
  if (Pile /*After test, variable should be Pile*/)
  {

      if (PileShape == "cake")
      {
         CreateCake(LAYER, LLVPMatrix); // Adjust LLVP Matrix for PremMatrix if modified a particlar layer and it contains LLVpS
      }

      else if (PileShape == "pancake")
      {
         CreatePanCake(LAYER, LLVPMatrix);
      }

  }

  //And add Earth Layers to TOP Volume
  for (int i = 0; i < PremMatrix.size(); ++i)
  {
        top->AddNode(LAYER[i],1);
  }

  gGeoManager->CloseGeometry(); // Finish Geometry
  //gGeoManager->SetTopVisible(); // the TOP is invisible
  

  
 
  // Tracking----------------------------------------------------------------------------------------------------------------------------------------

  //GetInitialPosition(zen,azi)
  //Calculate Paths inside the Earth

  //Direction of neutrino in spherical coordiates: https://mathworld.wolfram.com/SphericalCoordinates.html
  //double zen = 180.0; // zenith angle respect to the detector location in degrees (90,180]

 // double azi = 0.0;  //Azimuthal angle respect to the detector location in degrees [0,360]

  
  // double th = TMath::Pi()*( 1-(zen/180.0) ); //Polar angle in spherical coordinates in rads [pi/2, pi] 

  //double phi = TMath::Pi()*(azi)/180.0; //Azimuthal angle in spehericla coordinates in rads [0, 2*pi]
   


  double r = Ratm; // Max Radius
  
  //Incoming neutrino direction Cosines (respect to detector location)
  double nx = sin(th)*cos(phi);
  double ny = sin(th)*sin(phi);
  double nz = cos(th);


  ROOT::Math::SVector<double, 3> u(nx, ny, nz); // Direction of incoming neutrino
  ROOT::Math::SVector<double, 3> o(Det[0], Det[1], Det[2]); // Origin of the line (Set to be Detector location)
  ROOT::Math::SVector<double, 3> cc(0,0,0); // Origin of the Sphere 
 
  /* To calculate the initial position of the neutrino we need to compute 
  the intersection betweem 3D line (Neutrino direction) and Sphere(Earth) */

  double a = ROOT::Math::Dot( u, u);
  double b = 2*ROOT::Math::Dot( u, o - cc);
  double c = ROOT::Math::Dot( o-cc, o-cc) - r*r ;

  //Initial position of the neutrino... It should be located somewhere in the Atmosphere (Outer Layer of the Earth)
  double d1 = (-b+sqrt(b*b -4*a*c))/(2*a); // Distance 1 from the origin of the Line (Detector Location)
  double d2 = (-b-sqrt(b*b -4*a*c))/(2*a); // Distance 2 from the origin of the Line (Detector Location)

  
  //Coordinates of the Neutrino Positionin the Earth

  double xo = o[0]+d1*u[0];
  double yo = o[1]+d1*u[1];
  double zo = o[2]+d1*u[2];

  // Neutrino Direction

  double ni = -u[0];
  double nj = -u[1];
  double nk = -u[2];
  
  gGeoManager->InitTrack(xo, yo, zo, ni, nj, nk); // d*u is the point in the sphere

  TPolyMarker3D *l1 = new TPolyMarker3D(2,2); //Markers indicating DetLoc and Inical Neutrino
  TPolyMarker3D *LLVPpoints = new TPolyMarker3D(50,4);
  l1->SetPoint( 0 , o[0],o[1],o[2]); //Detector
  l1->SetPoint( 1 ,xo,yo,zo); //Incoming Neutrino


  TPolyLine3D *l2 = new TPolyLine3D(); // Lines that represent neutrino Paths.

 int s = 0;

  int i = 0;

  double Dx2, Dy2, Dz2, Dnorm;

  double sumL = 0; //Track Total Baseile

  //Neutrino Propagation inside the Earth

  //double v = 0;

  double l_LLVP = 0; // Distance traveled through LLVP
  while (!gGeoManager->IsOutside ())
  {  
         //++v;

         ////std::cout << v << std::endl;
      
         int nodeID= gGeoManager->GetCurrentNode()->GetIndex();

         const Double_t *cpoint = gGeoManager->GetCurrentPoint(); // Current nuetrino Position

         Dx2 = (cpoint[0]-Det[0])*(cpoint[0]-Det[0]);
         Dy2 = (cpoint[1]-Det[1])*(cpoint[1]-Det[1]);
         Dy2 = (cpoint[2]-Det[2])*(cpoint[2]-Det[2]);

         Dnorm = round(sqrt( Dx2 + Dy2 + Dz2 ));


         if (Dnorm <= 0.0 )  //Passing detector
         {
            //std::cout<< "Neutrino has reach detectio point" << Dnorm << std::endl; 
            break;
         
         }

         

         const char *path = gGeoManager->GetPath();

         TGeoVolume *cvol = gGeoManager->GetCurrentVolume();

         TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial(); //Material of the current Boundary

         Double_t safety = gGeoManager->GetSafeDistance(); //Distance to the next boundary/Neutrino Baseline

         gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.
          
         Double_t Li = gGeoManager->GetStep(); //Baseline Segement


         
         std::cout << "*** Distance from current location to " << gGeoManager->GetCurrentNode()->GetName() << " is: " << Li << " | "<< th*180.0/TMath::Pi() << " " << phi*180.0/TMath::Pi()  << std::endl;

         
            if (strstr(gGeoManager->GetCurrentNode()->GetName(), "LLVPLayer") != nullptr) { // strstr returns nullptr if not found
               l_LLVP += Li; // Distance traveled through LLVP
               //std::cout << "Selecting only LLVP info************************* "<< l_LLVP  << std::endl;

               LLVPpoints->SetPoint( s , cpoint[0],cpoint[1],cpoint[2]); //Detector
               s += 1;
            }
            else if(strstr(gGeoManager->GetCurrentNode()->GetName(), "LLVPLayer") == nullptr && l_LLVP > 0)
            {
               l_LLVP += Li; // Distance traveled through LLVP
               //std::cout << " All layer crossed "<< l_LLVP  << std::endl;

               LLVPpoints->SetPoint( s , cpoint[0],cpoint[1],cpoint[2]); //DetectoR

               Baseline_LLVP = l_LLVP;

               l_LLVP = 0;

            }
         
         // 

         double R_i = sqrt(cpoint[0]*cpoint[0]+ cpoint[1]*cpoint[1] + cpoint[2]*cpoint[2]); //Current Radius

        if (cvol->GetMedium()->GetId() == 99)
        {

          continue;

        }

         EarthPath.push_back({Li, cmat->GetDensity(), cmat->GetZ(),R_i}); //Store Path information [Baseline segment, Density, Z/A, R_i]
        
         sumL = sumL + Li; //Track Total Segment (Max = -2*R_earth*cos(zen))

         l2->SetPoint( i , cpoint[0],cpoint[1],cpoint[2]);

         i += 1;       
  }

   l1->SetMarkerColor(6);
   l1->SetMarkerSize(3);
   l2->SetLineColor(14);
   l2->SetLineWidth(3);


   //Set model visibility
   if (ModelVisibility)
   {
         /* code */
   

      top->Draw("apl");
      l1->Draw("same");
      l2->Draw("same");
      TView *view = gPad->GetView();
      view->ShowAxis();

      //LLVPpoints->Draw("same");

      EarthCanvas->Modified();
      EarthCanvas->Update();
      
   }
   //std::cout << "Canvas is good" << std::endl;

   

   //std::cout << " PileThickness " << PileThickness << " " << aWidth << std::endl;
   //std::cout << " Direction of incoming Neutrino: " << xo << " " << yo << " " << zo << std::endl;
   //std::cout << " Direction of coordinattes: " << o[0] << " " << o[1] << " " << o[2] << std::endl;
   //std::cout << " Neutrino Source location "<< sqrt(xo*xo + yo*yo + zo*zo) << " Detector location: " << sqrt( o[0]*o[0] + o[1]*o[1] + o[2]*o[2] ) <<std::endl;

   return EarthPath;}

// Create LLVP

 std::vector<std::vector<double>> Earth3DModel::Create3DPath () { return Earth3DPath ( zenith , azimuth , premtable); }