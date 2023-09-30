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
//Brief implementation of the EntangledGeneratorAction class

/*
In the constructor, the values for the second gamma particle of the same energy, opposite momentum direction (toward the second scatterer) and perpendicular polarisation, compared to the first particle created in 'BellsPrimaryGeneratorAction' are set and the particle origin placed in the centre of the world. 
Again, the 'physical' event is created in the GeneratePrimaries() function. In case the 'first' events have been created with random direction and polarisation, these values must be passed on (as arguments maybe? - then the function must be also changed where called in the QERunManager) and the opposite/perpendicular values set for the second photons.
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

//Constructor
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

//Destructor
EntangledGeneratorAction::~EntangledGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EntangledGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //G4cerr << "here" << G4endl;
  /* 
  // Randomized Polarization
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
  // Perpendicular polarization to 'first correlated photon'
  // This could be done a lot simpler, given that we are essentially dealing with a two-component vector ^^
  extern G4ThreeVector QEPolarization;

  G4double x0 = QEPolarization.x();
  G4double y0 = QEPolarization.y();

  G4cout << "Pol(x) = " << QEPolarization.x() << ", Pol(y) = " << QEPolarization.y() << G4endl;

  G4double r1 = G4UniformRand()*twopi;
  G4double x1;
  G4double y1;

  if(x0 == 0){
      y1 = 0;
      if(std::sin(r1) < 0) x1 = -1;
      else x1 = 1;
  }
  else if(y0 == 0){
      x1 = 0;
      if(std::sin(r1) < 0) y1 = -1;
      else y1 = 1;
  }

  x1 = std::sin(r1);
  y1 = (-1)*x0*x1/y0;

  //G4cout << "x1: " << x1 << " y1: " << y1 << G4endl;

  G4ThreeVector NewPolarization = G4ThreeVector(x1,y1,0).unit();
  x1 = NewPolarization.x();
  y1 = NewPolarization.y();

  // Test
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH2(0, x0, y0);
  analysisManager->FillH2(1, x1, y1);

  //G4double angle = NewPolarization.angle(QEPolarization);
  //G4double scalarProd = NewPolarization.dot(QEPolarization);

  //G4cout << "NewPol(x) = " << x1 << ", NewPol(y) = " << y1 << G4endl;
  //G4cout << "Angle = " << angle << ", Scalar Product = " << scalarProd << G4endl;

  fParticleGun->SetParticlePolarization(NewPolarization);
*/
  //Create the vertex of the particle passing it to the G4Event pointer
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
