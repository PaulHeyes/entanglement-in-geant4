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
//Brief Implementation of the BellsPhysicsList class

/*
The used physics is registered. Here the idea is to implement a new physics class, probably using a standard one as a base class (e.g. EmStandardPhysics or EmLivermorePolarizedPhysics) containing the 'new' physics for the second photon. I.e. the information about the scattering angle of the first photon must be passed on, modifying the scattering angle probability distribution of the second photon (see paper by Snyder, Pasternack et. al.).

Open: 
How is the angle going to be passed on?
How is the physics going to be 'activated' just for the second particle?

- A new particle class is created (gamma2), for which the physics is switched on. The particle class must be registered in the particleTable (how is this done?). After that, switching on the physics for this new particle class is simple (i.e. if (particle='gamma2') QEPhysics). The particle class could contain the scattering angle of the incident photon.

- The scattering angle is set to zero or deleted by default after every 2 events. The new physics can be switched on if a scattering angle (not equal 0) exists. In this case no new particle class is needed, a standard gamma may be used. This case may prove simpler!
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*

#include "BellsPhysicsList.hh"
//#include "BellsPhysicsListMessenger.hh"

#include "PhysListEmStandard.hh"
#include "PhysListEmLivermore.hh"
//#include "PhysListEmPenelope.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BellsPhysicsList::BellsPhysicsList()
: G4VModularPhysicsList(),
  fEmPhysicsList(0)
  //fMessenger(0)
{
  G4LossTableManager::Instance();
  SetDefaultCutValue(1.0*mm);

  //fMessenger = new BellsPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // EM physics
  fEmName = G4String("livermore");
  fEmPhysicsList = new PhysListEmLivermore(fEmName);

  //add new units for cross sections
  //
  new G4UnitDefinition( "mm2/g", "mm2/g","Surface/Mass", mm2/g);
  new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BellsPhysicsList::~BellsPhysicsList()
{
  //delete fMessenger;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BellsPhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "standard") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmStandard(name);

  } else if (name == "livermore") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmLivermore(name);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
*/


#include "BellsPhysicsList.hh"

//#include "PhysListEmStandard.hh"
//#include "PhysListEmLivermore.hh"
#include "EmLivermorePolarizedPhysics.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4PhysicsConstructorFactory.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4QEGamma.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4VUserRegionInformation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Constructor
BellsPhysicsList::BellsPhysicsList()
: G4VModularPhysicsList(),
  fLivPolPhysicsList(0)
{

  // Create a modular physics list and register only a 
  // few modules for it: EM interactions (decay of 
  // particles and radioactive decay) and QEPhysics (quantum 
  // entanglement physics).

  SetVerboseLevel(1);

  fLivPolPhysicsList = new EmLivermorePolarizedPhysics();

  // Default Decay Physics
  //RegisterPhysics(new G4DecayPhysics());

  // Default Radioactive Decay Physics
  //RegisterPhysics(new G4RadioactiveDecayPhysics());

  // Standard EM Physics
  // RegisterPhysics(new PhysListEmStandard());

  // Extended EM Physics
  // RegisterPhysics(new G4EmLivermorePolarizedPhysics());
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Destructor
BellsPhysicsList::~BellsPhysicsList()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BellsPhysicsList::ConstructParticle()
{
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
    
    G4QEGamma::QEGammaDefinition();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmProcessOptions.hh"

void BellsPhysicsList::ConstructProcess()
{
  // Transportation
  //
  AddTransportation();

  // Electromagnetic physics list
  //
  fLivPolPhysicsList->ConstructProcess();

  // Em options
  //
  G4EmProcessOptions emOptions;
  emOptions.SetIntegral(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void BellsPhysicsList::AddPhysicsList()
{
    fLivPolPhysicsList = new EmLivermorePolarizedPhysics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BellsPhysicsList::SetCuts()
{
  // Default SetCuts() provided by the base class
  G4VUserPhysicsList::SetCuts();
  
    // Production thresholds for detector regions
  G4Region* region;
  G4String regName;
  G4ProductionCuts* cuts;
  
  regName = "Matter";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(0.1*um,G4ProductionCuts::GetIndex("gamma"));
  //cuts->SetLowEdge(1*eV);
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);
  region->SetProductionCuts(cuts);
}  


