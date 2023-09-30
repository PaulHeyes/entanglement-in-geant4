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
//Brief Implementation of the QERunManager class

/*
This version is derived from Geant4's G4RunManager class. It is overridden in order to create a quantum entangled gamma event.
In the event loop a second event is created, processed, analysed and terminated after every first event, which is implemented in the ProcessQEEvent(), GenerateQEEvent and TerminateQEEvent() methods.

ProcessQEEvent()
The 'GeneratePrimaries()' method of the EntangledGeneratorAction is called. I am unsure at what point all the values are set, since the constructor wasn't explicitly called (but I might be wrong, and since it works, I must be wrong...).

TerminateQEEvent()
Works in complete analogy to TerminateOneEvent(), i.e. current event is terminated.
*/

#include "QERunManager.hh"

#include "G4UImanager.hh"
#include "G4ScoringManager.hh"
#include "G4Timer.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "EntangledGeneratorAction.hh"
#include "G4Run.hh"
#include "EventAction.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Constructor (If the constructor is called a second time, the
//constructor of base will generate an exception)
QERunManager::QERunManager():
  G4RunManager(),
  currentRun(0),
  userPrimaryGeneratorAction(0),
  //QEGenerator(0),
  currentQEEvent(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Destructor
QERunManager::~QERunManager(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Implementation of second (quantum entangled) event 
//in the event loop

 void QERunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
 {
  InitializeEventLoop(n_event,macroFile,n_select);
 
 //Event loop
  for(G4int i_event=0; i_event<2*n_event; i_event+=2 )
  {
    
  extern G4double dQEPhi;
  dQEPhi = 0.;
  extern G4double QEPhi;
  QEPhi = 0.;
  extern G4double polPhi;
  polPhi = 0.;

  // 'First event'
  ProcessOneEvent(i_event);
  TerminateOneEvent();

  // 'Quantum entangled event'
  ProcessQEEvent(i_event+1);
  TerminateQEEvent();

  if(runAborted) break;  
  }
  TerminateEventLoop();
 }
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 void QERunManager::ProcessQEEvent(G4int i_event)
 {
   //G4cerr << "ProcessQE started" << G4endl;
  currentQEEvent = GenerateQEEvent(i_event);
  //G4cerr << "QEEvent created" << G4endl;
  eventManager->ProcessOneEvent(currentQEEvent);
  //G4cerr << "PreAnalysis" << G4endl;
  AnalyzeEvent(currentQEEvent);
  UpdateScoring();
  if(i_event<n_select_msg) G4UImanager::GetUIpointer()->ApplyCommand(msgText);
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 G4Event* QERunManager::GenerateQEEvent(G4int i_event)
 {
  EntangledGeneratorAction* QEGenerator = new EntangledGeneratorAction;

  if(!QEGenerator)
  {
  G4Exception("QERunManager::GenerateQEEvent()", "Run0032", FatalException,
  "QEGenerator is not defined!");
  return 0;
  }

  G4Event* anEvent = new G4Event(i_event);
 
  if(storeRandomNumberStatusToG4Event==1 || storeRandomNumberStatusToG4Event==3)
  {
  std::ostringstream oss;
  G4Random::saveFullState(oss);
  randomNumberStatusForThisEvent = oss.str();
  anEvent->SetRandomNumberStatus(randomNumberStatusForThisEvent);
  }
 
  if(storeRandomNumberStatus) {
  G4String fileN = "currentEvent";
  if ( rngStatusEventsFlag ) {
  std::ostringstream os;
  os << "run" << currentRun->GetRunID() << "evt" << anEvent->GetEventID();
  fileN = os.str();
  }
  StoreRNGStatus(fileN);
  } 

  if(printModulo > 0 && anEvent->GetEventID()%printModulo == 0 )
  { G4cout << "--> Event " << anEvent->GetEventID() << " starts." << G4endl; }

  QEGenerator->GeneratePrimaries(anEvent);
 
  return anEvent;
 }
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 void QERunManager::TerminateQEEvent()
 {
  StackPreviousEvent(currentQEEvent);
  currentQEEvent = 0;
  numberOfEventProcessed++;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
