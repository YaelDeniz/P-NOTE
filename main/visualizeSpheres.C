#include "TGeoManager.h"
#include "TGeoSphere.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"

void visualizeFinerSpheres() {
    // Create a TGeoManager instance
    TGeoManager *geom = new TGeoManager("spheres", "Visualization of fine concentric spheres");

    // Define materials and media
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *matSolid = new TGeoMaterial("SolidMaterial", 1, 1, 1);
    TGeoMedium *vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
    TGeoMedium *solidMed = new TGeoMedium("SolidMaterial", 2, matSolid);

    // Create the world volume (outermost container)
    TGeoVolume *top = geom->MakeBox("World", vacuum, 100, 100, 100); // Large enough to contain the spheres
    geom->SetTopVolume(top);

    // Create spheres with finer segmentation
    TGeoSphere *outerSphereShape = new TGeoSphere(50, 60, 0, 180, 0, 360); // Outer sphere
    outerSphereShape->SetNseg(100); // Finer segmentation for latitude/longitude

    TGeoSphere *middleSphereShape = new TGeoSphere(30, 50, 0, 180, 0, 360); // Middle sphere
    middleSphereShape->SetNseg(100); // Finer segmentation for latitude/longitude

    TGeoSphere *innerSphereShape = new TGeoSphere(0, 30, 0, 180, 0, 360); // Inner sphere
    innerSphereShape->SetNseg(100); // Finer segmentation for latitude/longitude

    // Create volumes from the custom sphere shapes
    TGeoVolume *outerSphere = new TGeoVolume("OuterSphere", outerSphereShape, vacuum);
    TGeoVolume *middleSphere = new TGeoVolume("MiddleSphere", middleSphereShape, vacuum);
    TGeoVolume *innerSphere = new TGeoVolume("InnerSphere", innerSphereShape, solidMed);

    // Add spheres to the world volume
    top->AddNode(outerSphere, 1);
    top->AddNode(middleSphere, 1);
    top->AddNode(innerSphere, 1);

    // Set visualization styles
    outerSphere->SetLineColor(kBlue);       // Outer sphere color
    outerSphere->SetTransparency(90);       // Transparent
    middleSphere->SetLineColor(kGreen);     // Middle sphere color
    middleSphere->SetTransparency(50);      // Semi-transparent
    innerSphere->SetLineColor(kRed);        // Inner sphere color
    innerSphere->SetTransparency(0);        // Opaque (solid)

    // Draw the geometry
    geom->CloseGeometry();
    top->Draw("ogl"); // Use OpenGL viewer for better visualization
}