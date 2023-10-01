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

// Paul Heyes, 2019
// Brief implementation of the EntangledGeneratorAction class

/*
Instead of generating actual annihilations in the Geant4 framework, an easier route is taken here,
namely, two photons with opposite momentum direction and respective energy of 511 keV are created. 

In the constructor, the values for the second gamma particle of the same energy, opposite momentum direction 
(toward the second scatterer) and perpendicular polarisation, compared to the first particle created in 
'BellsPrimaryGeneratorAction' are set and the particle origin placed in the centre of the world. 

Again, the 'physical' event is created in the GeneratePrimaries() function. In case the 'first' events have 
been created with random direction and polarisation, these values must be passed on and the 
opposite/perpendicular values set for the second photons (in GeneratePrimaries()).
*/

#include "EntangledGeneratorAction.hh"

#include "QERunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

//#include "G4Gamma2.hh"

#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "BellsAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor
EntangledGeneratorAction::EntangledGeneratorAction()
  : BellsPrimaryGeneratorAction(),
    fParticleGun(0)
{
 
  G4int n_particle = 1;

  fParticleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* particle
    = particleTable->FindParticle("QEgamma");

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(511*keV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  
  extern G4double polPhi;  
  if (G4UniformRand() < 0.5) { polPhi = -pi/2.; }
  else{ polPhi = pi/2.; }
  
  fParticleGun->SetParticlePolarization(G4ThreeVector(std::sin(polPhi),0.,0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Destructor
EntangledGeneratorAction::~EntangledGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EntangledGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  // Create the vertex of the particle passing it to the G4Event pointer
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
