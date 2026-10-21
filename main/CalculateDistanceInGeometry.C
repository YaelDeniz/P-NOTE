#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoMedium.h"

void CalculateDistanceInGeometry() {
    // Create geometry manager
    TGeoManager *geom = new TGeoManager("Geometry", "Example Geometry");

    // Define materials and media
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMedium *vacuum = new TGeoMedium("Vacuum", 1, matVacuum);

    // Create the world volume (a big box that contains all geometry)
    TGeoVolume *world = geom->MakeBox("World", vacuum, 210, 210, 210);
    geom->SetTopVolume(world);

    // Add a target volume (e.g., a smaller box within the world)
    TGeoVolume *targetBox = geom->MakeBox("TargetBox", vacuum, 100, 100, 100);
    world->AddNode(targetBox, 1, new TGeoTranslation(0, 0, 0));

    // Close the geometry
    geom->CloseGeometry();

    // Define starting point and direction
    Double_t startPoint[3] = { -110.0, -110.0, 0.0 }; // Starting outside the target
    Double_t direction[3] = { 0.70710678118, 0.70710678118, 0.0 };      // Moving along the x-axis

    // Initialize navigation
    geom->InitTrack(startPoint, direction);

    // Track distance
    Double_t distance = 0.0;

    TPolyLine3D *l = new TPolyLine3D(); // Lines that represent neutrino Paths.

    TPolyMarker3D *ledges = new TPolyMarker3D(10,2); //Markers indicating DetLoc and Inical Neutrino
    ledges->SetPoint( 0 , startPoint[0], startPoint[1], startPoint[2]); //Detector

    int i = 0;

    bool insideTarget = false;
    double totalDistanceInTarget = 0;

while (true) {
    // Get the current volume and step distance to the next boundary
    TGeoNode *currentNode = geom->GetCurrentNode();
    const Double_t *cpoint = geom->GetCurrentPoint(); // Current position
    Double_t step = geom->GetStep();

    std::cout << "Distance from current location to " << currentNode->GetName() << " is: " 
              << step << " | Step index: " << i 
              << " | Next boundary: " << geom->FindNextBoundary() << std::endl;

    // Store the current point for plotting
    l->SetPoint(i, cpoint[0], cpoint[1], cpoint[2]);
    ledges->SetPoint(i + 1, cpoint[0], cpoint[1], cpoint[2]); // Detector

    // Check if the particle has entered the target volume
    if (currentNode->GetName() == TString("TargetBox_1")) {
        
        std::cout << "if #1: Evaluating if particle is entering the box, current state: " << insideTarget << std::endl;

        if (!insideTarget) {

            std::cout << "Entering TargetBox_1" << std::endl;
            insideTarget = true; // Mark that we are inside the target
            std::cout << "inner if: particle is entering the box, changing current state: " << insideTarget << std::endl;
        }

        // Move the particle to the next boundary
        geom->FindNextBoundaryAndStep();

        // Increment the step index
        i += 1;

        continue;
        
        //totalDistanceInTarget += step; // Accumulate distance in target volume
    }
    // Check if the particle has exited the target volume
    else if (insideTarget) {
        // We were inside TargetBox_1 but now are not
        totalDistanceInTarget += step; // Accumulate distance in target volume

        std::cout << " Particle is leaving TargetBox_1, Total distance traveled in TargetBox_1: " 
                  << totalDistanceInTarget << " units" << std::endl;
        insideTarget = false;
        
        totalDistanceInTarget = 0.0; // Reset for the next time we enter
    }

    // Calculate the distance to the next boundary and check if we're still inside
    if (geom->IsOutside()) {
        std::cout << "Outside the geometry" << std::endl;
        break; // Exit the loop if the particle is outside the entire geometry
    }

    // Move the particle to the next boundary
    geom->FindNextBoundaryAndStep();

    // Increment the step index
    i += 1;
}


  TCanvas * c =   new TCanvas("3D Earth", "3D Earth",0,0,600,600);

  gGeoManager->CloseGeometry(); // Finish Geometry
  gGeoManager->SetTopVisible(); // the TOP is invisible
  world->Draw();
  
  l->SetLineColor(14);
  l->SetLineWidth(3);

  ledges->SetMarkerColor(6);
  ledges->SetMarkerSize(3);

  l->Draw("same");
  ledges->Draw("same");


  TView *view = gPad->GetView();
  view->ShowAxis();

  c->Modified();
  c->Update();

  // Output the total distance traveled in the target volume
  printf("Total distance traveled in TargetBox: %f cm\n", totalDistanceInTarget);




}



/*

// Add step to total distance
        if (currentNode && currentNode->GetName() == TString("TargetBox")) {
            
                    std::cout<< "Inside" << std::endl;

            
            insideTarget = true;

            distance += step;
        }else if (insideTarget) {
            // Exit if we were inside and now left the target volume
            ledges->SetPoint( 1 ,cpoint[0],cpoint[1],cpoint[2]); //Incoming Neutrino

            break;
        }

        // Stop if the track exits the world volume
        if (!geom->FindNextBoundaryAndStep() || geom->IsOutside())
        {


            break;

        }

        // Move to the next point
       
        l->SetPoint( i , cpoint[0],cpoint[1],cpoint[2]);
        i += 1; 

while (true) {
    // Get current volume and next boundary
    TGeoNode *currentNode = geom->GetCurrentNode();
    const Double_t *cpoint = geom->GetCurrentPoint(); // Current neutrino position
    Double_t step = geom->GetStep();

    std::cout << "Distance from current location to " <<currentNode->GetName() << " is: " 
              << step << " | Step index: " << i << "Next boundary " << geom -> FindNextBoundary() <<std::endl;

    // Store the current point for plotting
    l->SetPoint(i, cpoint[0], cpoint[1], cpoint[2]);
    ledges->SetPoint(i + 1, cpoint[0], cpoint[1], cpoint[2]); // Detector

    if (currentNode->GetName() == TString("TargetBox_1"))
    {
        std::cout << "Entering Target" << std::endl;
        
    } else if (step > 0 && currentNode->GetName() != TString("TargetBox_1") )

    {

        if  (geom->IsOutside()) {
        std::cout << " Also Outside target the geometry" << std::endl;
        }else{
            std::cout << "Distance to Target " << step << std::endl;

        }


    }

    // Calculate the distance to the next boundary and check if we're still inside
    if  (geom->IsOutside()) {
        // Exit the loop if the particle has no more steps or is outside the geometry
        std::cout << " Outside the geometry" << std::endl;
        break;
    }

    // Move the particle to the next boundary
    geom->FindNextBoundaryAndStep();

    // Increment the index
    i += 1;
}
*/