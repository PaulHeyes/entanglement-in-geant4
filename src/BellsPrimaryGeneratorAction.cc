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

// Paul Heyes, 2019
// This class file describes how to generate the primary events

/*
In the constructor of BellsPrimaryGeneratorAction the number of particles per event is set to one and the type 
of particle is specified (i.e. the novel QEGamma particle is selected from the particleTable). Its point of origin 
(the centre of the world), energy of 511keV, momentum direction (toward the first scatterer) and polarisation are set. 
Note: At this point however no particle is 'physically' created.

While all the values are initially set in the constructor, the actual event is created in the GeneratePrimaries() function. 
Here a random momentum direction and polarisation may be implemented in the future. 
(Note: This must happen here, since only the GeneratePrimaries() function is called every event. 
If it is done in the constructor, every particle during one run will have the same identical values!)

Note: The annihilation is replaced by two particle guns (particles being generated in BellsPrimaryGeneratorAction and 
EntangledGeneratorAction), with the physical parameters being hard-coded. This is done for the sake of simplicity - 
randomizing the direction (as described briefly above) or even implementing the model for 'actual' annihilations taking
place during a simulation run (creating QEGammas and passing relevant values) should be possible. 
*/

#include "BellsPrimaryGeneratorAction.hh"

#include "QERunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

//#include "G4Gamma2.hh"

#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Constructor
BellsPrimaryGeneratorAction::BellsPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
 
  G4int n_particle = 1; 			//Number of particles

  fParticleGun  = new G4ParticleGun(n_particle);
    
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* particle
    = particleTable->FindParticle("QEgamma");		//Particle

//Particle kinematics
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(511*keV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticlePolarization(G4ThreeVector(0.,1.,0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Destructor
BellsPrimaryGeneratorAction::~BellsPrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void BellsPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //extern G4ThreeVector QEPolarization;

  /* 
  // Randomized Polarization (Momentum direction may be done the same way)
  //
  G4double phi = G4UniformRand()*twopi *rad;
  G4double theta = G4UniformRand()*pi *rad;
  G4double x0  = std::sin(theta)*std::cos(phi);
  G4double y0  = std::sin(theta)*std::sin(phi);
  G4double z0  = std::cos(theta); 
  // x0 += dx0*(G4UniformRand()-0.5);
  //  y0 += dy0*(G4UniformRand()-0.5);
  // z0 += dz0*(G4UniformRand()-0.5);
  fParticleGun->SetParticlePolarization(G4ThreeVector(x0,y0,z0));
  fParticleGun2->SetParticlePolarization(G4ThreeVector(-x0,-y0,-z0));
  */

    /*
    G4double r = G4UniformRand()*twopi;       //r in [0,2*pi]
    G4double x0 = std::sin(r);
    G4double y0 = std::cos(r);
    G4double z0 = 0;
    //x0 += x0*(G4UniformRand()-0.5);
    //y0 += y0*(G4UniformRand()-0.5);

    QEPolarization = G4ThreeVector(x0,y0,z0);

    fParticleGun->SetParticlePolarization(G4ThreeVector(x0,y0,z0));
      */

  //Creation of the vertex of the primary particle passing it to the G4Event pointer
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
