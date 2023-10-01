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
//

//Paul Heyes, 2019
//Brief Implementation of the Hit class

/*
The SDHit class records values of interest (as defined in the ScatterSD class too):
the scatterer ID, the particle time and the particle local and global positions, energy,
momentum direction, polarization (for most of these parameters the values both before and
after the hit, i.e., preStep and postStep points).

In the Print() method many values are written to the console. Useful for testing, but should
definitely be switched off for 'larger' simulations.
*/

#include "SDHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<SDHit>* HitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDHit::SDHit()
: G4VHit(), fId(-1), fTime(0.), fLocalPos(0), fWorldPos(0), fLocalPos2(0), fWorldPos2(0), fMomDir(0), fMomDir2(), fPol(), fE(), fE2()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDHit::SDHit(G4double t)
: G4VHit(), fId(-1), fTime(t), fLocalPos(0), fWorldPos(0), fLocalPos2(0), fWorldPos2(0), fMomDir(0), fMomDir2(), fPol(), fE(), fE2()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

SDHit::~SDHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDHit::SDHit(const SDHit &right)
: G4VHit() {
    fId = right.fId;
    fTime = right.fTime;
    fLocalPos = right.fLocalPos;
    fWorldPos = right.fWorldPos;
    fLocalPos2 = right.fLocalPos2;
    fWorldPos2 = right.fWorldPos2;
    fMomDir = right.fMomDir;
    fMomDir2 = right.fMomDir2;
    fPol = right.fPol;
    fE = right.fE;
    fE2 = right.fE2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const SDHit& SDHit::operator=(const SDHit &right)
{
    fId = right.fId;
    fTime = right.fTime;
    fLocalPos = right.fLocalPos;
    fWorldPos = right.fWorldPos; 
    fLocalPos2 = right.fLocalPos2;
    fWorldPos2 = right.fWorldPos2;
    fMomDir = right.fMomDir;
    fMomDir2 = right.fMomDir2;
    fPol = right.fPol;
    fE = right.fE;
    fE2 = right.fE2;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int SDHit::operator==(const SDHit &/*right*/) const
{
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDHit::Draw()  //is this really necessary?
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (pVVisManager)
    {
        G4Circle circle(fWorldPos);
        circle.SetScreenSize(2);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1.,1.,0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* SDHit::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("Hit",isNew);

    if (isNew) {
        (*store)["HitType"] 
          = G4AttDef("HitType","Hit Type","Physics","","G4String");
        
        (*store)["ID"] 
          = G4AttDef("ID","ID","Physics","","G4int");
        
        (*store)["Time"] 
          = G4AttDef("Time","Time","Physics","G4BestUnit","G4double");
        
        (*store)["Pos"] 
          = G4AttDef("Pos", "Position", "Physics","G4BestUnit","G4ThreeVector");

        (*store)["Pol"]
          = G4AttDef("Pol", "Polarization", "Physics","G4BestUnit","G4ThreeVector");
    }
    return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* SDHit::CreateAttValues() const
{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
    
    values
      ->push_back(G4AttValue("HitType","Hit",""));
    values
      ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fId),""));
    values
      ->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
    values
      ->push_back(G4AttValue("Pos",G4BestUnit(fWorldPos,"Length"),""));
    values
      ->push_back(G4AttValue("Pol",G4BestUnit(fPol,"Length"),""));
    
    return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDHit::Print()
{
  //     G4cout << "  Scatterer[" << fId << "] : time " << fTime/ns
  //  << " (nsec) --- local (x,y,z) " << fLocalPos.x()
  //  << ", " << fLocalPos.y() << ", " << fLocalPos.z() << G4endl;

  //     G4cout << "  Scatterer[" << fId << "] : time " << fTime/ns
  //  << " (nsec) " << G4endl;


    G4cout << " ** PreStepPoint ** " << G4endl;

    /*
      // Global Pos
      G4cout << " Global (x,y,z) " << fWorldPos.x()
    << ", " << fWorldPos.y() << ", " << fWorldPos.z() << G4endl;

      // Local Pos
     G4cout << " Local (x,y,z) " << fLocalPos.x()
    << ", " << fLocalPos.y() << ", " << fLocalPos.z() << G4endl;
*/

    // Momentum Direction
       G4cout << "  ScatterDet: Direction (x,y,z) " << fMomDir.x()
    << ", " << fMomDir.y() << ", " << fMomDir.z() << G4endl;

    // Energy
    G4cout << "  ScatterDet: Energy " << fE << G4endl;


    G4cout << " ** PostStepPoint ** " << G4endl;

    /*
    // Global Pos
    G4cout << " Global (x,y,z) " << fWorldPos2.x()
  << ", " << fWorldPos2.y() << ", " << fWorldPos2.z() << G4endl;

    // Local Pos
   G4cout << " Local (x,y,z) " << fLocalPos2.x()
  << ", " << fLocalPos2.y() << ", " << fLocalPos2.z() << G4endl;
*/

   // Momentum Direction
       G4cout << "  ScatterDet: Direction2 (x,y,z) " << fMomDir2.x()
    << ", " << fMomDir2.y() << ", " << fMomDir2.z() << G4endl;

   // Energy
       G4cout << "  ScatterDet: Energy2 " << fE2 << G4endl;

   // Polarization
       G4cout << "  ScatterDet: Polarization (x,y,z) " << fPol.x()
    << ", " << fPol.y() << ", " << fPol.z() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
