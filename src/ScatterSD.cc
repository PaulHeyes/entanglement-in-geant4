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
//Brief Implementation of the ScatterSD class

/*
Initialize()
This method allocates a hits collection object to the sensitive detector (scatterer).

ProcessHits()
This method is invoked by G4SteppingManager at some point during the processing or analysis of an event (unfortunately I wasn't able to find the exact point at which the function is called). It creates an object 'hit' as an instance of the Hit class, passing on variables such as time and position (local & global).
*/

#include "ScatterSD.hh"
#include "SDHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ParticleChange.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Constructor
ScatterSD::ScatterSD(G4String name) :
  G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
  G4String HCname = "scattererColl";
  collectionName.insert(HCname);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Destructor
ScatterSD::~ScatterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScatterSD::Initialize(G4HCofThisEvent *HCE)
{
  fHitsCollection = new SDHitsCollection(SensitiveDetectorName,collectionName[0]);
 
  if (fHCID<0)
    { fHCID=G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }

  HCE->AddHitsCollection(fHCID,fHitsCollection);
  
  // return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ScatterSD::ProcessHits (G4Step *step, G4TouchableHistory*)
{
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  G4Track *aTrack = step->GetTrack() ;

  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());

  G4double edep = step->GetTotalEnergyDeposit();
  // if(edep==0.) return false;    

  //G4double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
  //if (charge==0.) return true;
 
  G4VPhysicalVolume* motherPhysical = touchable->GetVolume(1); // mother
  G4int copyNo = motherPhysical->GetCopyNo();

  G4ThreeVector worldPos1 = preStepPoint->GetPosition();
  G4ThreeVector localPos1
    = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos1);

  G4ThreeVector worldPos2 = postStepPoint->GetPosition();
  G4ThreeVector localPos2
     = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos2);

  G4ThreeVector momentumDirection1 = preStepPoint->GetMomentumDirection();
  G4ThreeVector momentumDirection2 = postStepPoint->GetMomentumDirection();

  G4double Energy1 = preStepPoint->GetTotalEnergy();
  G4double Energy2 = postStepPoint->GetTotalEnergy();


  //G4StepStatus status = preStepPoint->GetStepStatus();
  //G4VProcess* process = preStepPoint->GetProcessDefinedStep();

  //PartChange = new G4ParticleChange;
  //const G4ThreeVector* Direction = PartChange->GetMomentumDirection();


  SDHit* hit = new SDHit(copyNo);
  hit->SetWorldPos(worldPos1);
  hit->SetWorldPos2(worldPos2);
  hit->SetLocalPos(localPos1);
  hit->SetLocalPos2(localPos2);
  hit->SetTime(preStepPoint->GetGlobalTime());
  hit->SetMomentumDirection(momentumDirection1);
  hit->SetMomentumDirection2(momentumDirection2);
  hit->SetEnergy(Energy1);
  hit->SetEnergy2(Energy2);
  hit->SetPolarization(postStepPoint->GetPolarization());

    
  fHitsCollection->insert(hit);

  // Stop and kill the particle in order to prevent any further events it could have!
  // aTrack->SetTrackStatus(fStopAndKill);
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScatterSD::EndOfEvent(G4HCofThisEvent* )
{}

