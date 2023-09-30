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
//Brief Definition of teh QERunManager class


#ifndef QERunManager_h
#define QERunManager_h 1

// Override the G4RunManager class so that the quantum entangled // event can be implemented in the event loop. 

//     Notes:
//     1) G4RunManager is an oddly implemented singleton.  Those who access methods
//        via the instance pointer in the base will never want to access any of the
//        methods in this class.  So we do not need an instance pointer to this class.
//

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "globals.hh"
#include "EntangledGeneratorAction.hh"
#include "EventAction.hh"

class G4VUserPrimaryGeneratorAction;
class G4Run;

class QERunManager : public G4RunManager{

public:
  QERunManager();
  virtual ~QERunManager();

public:

  virtual void DoEventLoop(G4int n_event,const char* macroFile=0,G4int n_select=-1);
  virtual void ProcessQEEvent(G4int i_event);
  virtual G4Event* GenerateQEEvent(G4int i_event);
  virtual void TerminateQEEvent();

  G4Run* currentRun;
  G4VUserPrimaryGeneratorAction* userPrimaryGeneratorAction;
  //EntangledGeneratorAction* QEGenerator;
  G4Event* currentQEEvent;

  //  inline const G4VUserPrimaryGeneratorAction* GetUserPrimaryGeneratorAction() const
  // { return userPrimaryGeneratorAction; }
 
};

#endif
