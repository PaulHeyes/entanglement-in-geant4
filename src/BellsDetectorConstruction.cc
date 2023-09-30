//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

//Paul Heyes, 16.09.2015
//Brief Implementation of the BellsDetectorConstruction class

/*
Construct() - method
The geometry is similar to that of the bells experiment, consisting of two scatterers and two detectors (placed in a world). The material used for all these objects is sodium iodide. In order to detect all possible angles at once the detectors are designed as cylindrical rings around the scatterers (detector one is blue, detector 2 is red). For closer resemblance to actual PET systems, crystals may be placed around these rings and act as detectors. In the extent of this simulation things are kept simpler though.

ConstructSDandField() - method
Here the before constructed detectors are assigned the attribute of actually detecting, i.e. G4VSensitiveDetector objects are assigned to the corresponding volumes. 
The scatterers also become sensitive detectors, this however is mainly implemented in order to check the program is working properly (i.e. check the energy of the photons and count the hits).
*/

#include "BellsDetectorConstruction.hh"
#include "QERunManager.hh"
#include "G4NistManager.hh"
#include "ScatterSD.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4VUserRegionInformation.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4ios.hh"

#include "G4TransportationManager.hh"
//#include "G4GlobalMagFieldMessenger.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Constructor
BellsDetectorConstruction::BellsDetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true),
  flogicScatterDetector1(0),
  flogicScatterDetector2(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Destructor
BellsDetectorConstruction::~BellsDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4VPhysicalVolume* BellsDetectorConstruction::Construct()
{  
  //
  // ***** Materials *****
  //
  G4NistManager* nist = G4NistManager::Instance();

  // **Retrieve Nist Materials** 
  //G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* vacuum_mat = nist->FindOrBuildMaterial("G4_Galactic");
  //G4Material* detector_mat = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  G4Material* scatter_mat = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  //G4Material* cryst_mat = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");



  //     
  // ***** World *****
  //
  G4double world_sizeXY = 4.0*m;
  G4double world_sizeZ  = 4.0*m;

  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size
  
  // World Logical Volume definition   
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        vacuum_mat,         //its material
                        "World");            //its name
                      
  // World Physical Volume Placement at (0,0,0)              
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps 
          

  //
  // *** ScatterDetectors ***
  //
  G4double scatter_radius = 30*cm;
  G4double scatter_dZ = 50*cm;
  G4double dScatter_dZ = 50*cm;
    
  G4Tubs* solidScatter =
    new G4Tubs("Scatterer", 0., scatter_radius, 0.5*scatter_dZ, 0., twopi);
      
  flogicScatterDetector1 =
    new G4LogicalVolume(solidScatter, scatter_mat, "ScatterLV");   

  G4Colour blue(0.0, 1.0, 0.0); //really blue is (0,0,1)... but it's a joke :P
  flogicScatterDetector1->SetVisAttributes(new G4VisAttributes(blue));
                                 
    new G4PVPlacement(0, G4ThreeVector(0., 0., dScatter_dZ), flogicScatterDetector1, "ScatterDet1", logicWorld, false, 0, fCheckOverlaps);

  // Repeat for Second Scatterer

  G4Tubs* solidScatter2 =
    new G4Tubs("Scatterer2", 0., scatter_radius, 0.5*scatter_dZ, 0., twopi);
      
  flogicScatterDetector2 =
    new G4LogicalVolume(solidScatter2, scatter_mat, "Scatter2LV"); 

  G4Colour red(1.0, 0.0, 0.0);
  flogicScatterDetector2->SetVisAttributes(new G4VisAttributes(red));
                                  
    new G4PVPlacement(0, G4ThreeVector(0., 0., -dScatter_dZ), flogicScatterDetector2, "ScatterDet2", logicWorld, false, 0, fCheckOverlaps);


  // Create a region 
  G4String regName = "Matter";
  G4Region* matter = new G4Region(regName); 
  // Attach a logical volume to the region 
  matter->AddRootLogicalVolume(flogicScatterDetector1); 
  matter->AddRootLogicalVolume(flogicScatterDetector2); 


  /*

  // Here is a rough guideline presented how to incorporate crystals into the
  // detector(s), as it is done in the tutorial. This code only creates one 
  // detector ring and isn't compatible with the rest of this class without 
  // some modification! I have mainly left it here since it might be useful 
  // in the future.

  //
  // ***** Crystal-detector *****
  //

  //Parameters

  G4double cryst_dX = 6*cm, cryst_dY = 10*cm, cryst_dZ = 3*cm;

  G4int nb_cryst = 16;
  G4int nb_rings = 5;
 
  G4double dPhi = twopi/nb_cryst, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  
  G4double ring_R1 = 0.5*cryst_dY/tandPhi;
  G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;

  G4double detector_dZ = nb_rings*cryst_dX;
 
  //
  // *** Ring ***
  //
  // define one ring as an "envelope" made of air. This will be filled 
  // by the crystals (daughter volumes). The logical volume of the ring 
  // (which contains all daughters) will be then placed many times
  //
  G4Tubs* solidRing =
    new G4Tubs("Ring",      //name
	       ring_R1,      //inner radius
	       ring_R2,      //outer radius
	       0.5*cryst_dX, //height
	       0.,           //start angle
	       twopi);       //spanning angle
      
  G4LogicalVolume* logicRing =                         
    new G4LogicalVolume(solidRing,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name

  //     
  // *** Crystal ***
  //
  G4double gap = 0.5*mm;        //a gap for wrapping
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
  G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);
                     
  G4LogicalVolume* logicCryst = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name

  G4Colour red(1.0, 0.0, 0.0);
  logicCryst->SetVisAttributes(new G4VisAttributes(red));
               

  // Place crystals within a ring (loop over crystals) 
  //
  for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {

    G4double phi = icrys*dPhi;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg); 
    rotm.rotateZ(phi);

    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);
                      
    // Place the crystal with the appropriate transformation
    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "crystal",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }

  // 
  // place rings within detector 
  //
  G4double OG = -0.5*(detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,OG), //position
                      logicRing,             //its logical volume
                      "ring",                //its name
                      flogicDetector1,         //its mother  volume
                      false,                 //no boolean operation
                      iring,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }

 */


                                        
  //Return the physical World
  return physWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void BellsDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname;

  G4VSensitiveDetector* ScatterDet1
    = new ScatterSD(SDname="/ScatterDet1");
   SDman->AddNewDetector(ScatterDet1);
  flogicScatterDetector1->SetSensitiveDetector(ScatterDet1);

  G4VSensitiveDetector* ScatterDet2
    = new ScatterSD(SDname="/ScatterDet2");
   SDman->AddNewDetector(ScatterDet2);
  flogicScatterDetector2->SetSensitiveDetector(ScatterDet2);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

