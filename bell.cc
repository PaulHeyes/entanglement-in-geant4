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

/*
Main method:
  - Difference to standard Geant4: specific run manager QERunManager constructed
  - Common user action classes set

In addition all external variables necessary are defined. They are later set/used as 
follows:
  - QEPhi_ref is set in the EventAction class (only used for data collection in histogram)
  - polPhi is set in EntangledGeneratorAction class
  - QEPhi, dQEPhi, QEeps1 & QEeps2 are all sampled in the QEComptonModel class and used
    for the calculation of the differential cross section of the entangled Compton scattering
Some of these variables are (re)set to zero in the event loop of the QERunManager (necessary 
for correct calculation of the entangled Compton scattering).
*/

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "QERunManager.hh"
#endif

#include "QERunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "BellsDetectorConstruction.hh"
#include "BellsPhysicsList.hh"
#include "BellsActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

// External Variables
G4double QEPhi_ref;
G4double QEPhi;
G4double polPhi;
G4double dQEPhi;
G4double QEeps1;
G4double QEeps2;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  //
  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
     
  // Construct the default run manager. Pick the proper run 
  // manager depending if the multi-threading option is 
  // active or not.
  //
   #ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  #else
  QERunManager* runManager = new QERunManager;
  #endif  

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new BellsDetectorConstruction);
  //
  runManager->SetUserInitialization(new BellsPhysicsList);

  // Set user action initialization
  //
  runManager->SetUserInitialization(new BellsActionInitialization());  
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // if an argument is given after the name of the executable 
  // (i.e. argc > 1), then take the argument as a Geant4 macro 
  // and execute it
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  // otherwise (only the executable is given), start a user 
  // interface session. An initialization macro is executed 
  // by default. The macro which is executed depends on the 
  // activation (or not) of the visualization
  else
    {  // interactive mode : define UI session
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
      UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
      // start the session here: make the Geant4 prompt Idle>
      // available to the user
      ui->SessionStart();
      delete ui;
#endif
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
