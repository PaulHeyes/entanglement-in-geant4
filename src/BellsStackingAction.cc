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

//Paul Heyes, 18.10.2016
//Implementation of the BellsStackingAction class

/*
Description
*/

#include "BellsStackingAction.hh"

#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "QERunManager.hh"
#include "G4Track.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor
BellsStackingAction::BellsStackingAction()
    : G4UserStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Destructor
BellsStackingAction::~BellsStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
BellsStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
    G4ClassificationOfNewTrack classification = fWaiting;

    G4int parent_ID = aTrack->GetParentID();
    G4ParticleDefinition* particleType = aTrack->GetDefinition();
    const G4String& particleName = particleType->GetParticleName();

    //G4cout << "Particle = " << particleName << G4endl;
    //G4cout << "Ekin = " << aTrack->GetKineticEnergy() << G4endl;
    //G4cout << "Etot = " << aTrack->GetTotalEnergy() << G4endl;

    /*
    const G4ThreeVector partDirection = aTrack->GetMomentumDirection();
    G4cout << "MomentumDirection (x,y,z): " << "(" << partDirection.x() << ", "
           << partDirection.y() << ", " << partDirection.z() << ")" << G4endl;
    */

    if(parent_ID == 0)
    {
        //G4cout << "-----------" << " Primary particle " << "-----------" << G4endl;
        classification = fUrgent;

        //G4ParticleDefinition* particleType = aTrack->GetDefinition();
        //if(particleType == G4Photon::G4PhotonDefinition) {}
        //    else G4cout << "This is not a Photon primary..." << G4endl;
    }
    
    if(particleName == "QEgamma" || particleName == "gamma")
	{
        classification = fUrgent;
    }    

    if(particleName != "gamma" && particleName != "QEgamma")
    {
        //G4cout << "-----------" << " Particle killed " << "-----------" << G4endl;
        classification = fKill;
    }

    if(parent_ID < 0)
    {
        G4cout << "Classification problem: Something is waiting..." << G4endl;
    }

    return classification;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BellsStackingAction::NewStage()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BellsStackingAction::PrepareNewEvent()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
