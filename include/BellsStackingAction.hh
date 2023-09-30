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

// Paul Heyes, 18.10.2016
// Definition of the BellsStackingAction class

#ifndef BellsStackingAction_h
#define BellStackingAction_h 1

class G4StackManager;
class G4Track;

#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"

class BellsStackingAction : public G4UserStackingAction
{
  public:
      BellsStackingAction();
      virtual ~BellsStackingAction();

  //protected:
  //    G4StackManager * stackManager;

  public:
     // inline void SetStackManager(G4StackManager * value)
     // { stackManager = value; }

      virtual G4ClassificationOfNewTrack
        ClassifyNewTrack(const G4Track* aTrack);

//    enum G4ClassificationOfNewTrack
//    {
//      fUrgent,    // put into the urgent stack
//      fWaiting,   // put into the waiting stack
//      fPostpone,  // postpone to the next event
//      fKill       // kill without stacking
//    };
//
//    The parent_ID of the track indicates the origin of it.

   // G4int parent_ID = aTrack->get_parentID();

//      parent_ID = 0 : primary particle
//                > 0 : secondary particle
//                < 0 : postponed from the previous event

  //  if(parent_ID > 0) fKill;        // Kill secondaries (because they are not needed and get in the way)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

      virtual void NewStage();

//    This method is called by G4StackManager when the urgentStack
//  becomes empty and contents in the waitingStack are transtered
//  to the urgentStack.
//    Note that this method is not called at the begining of each
//  event, but "PrepareNewEvent" is called.
//
//    In case re-classification of the stacked tracks is needed,
//  use the following method to request to G4StackManager.
//
//    stackManager->ReClassify();
//
//  All of the stacked tracks in the waitingStack will be re-classified
//  by "ClassifyNewTrack" method.
//    To abort current event, use the following method.
//
//    stackManager->clear();
//
//  Note that this way is valid and safe only for the case it is called
//  from this user class. The more global way of event abortion is
//
//    G4UImanager * UImanager = G4UImanager::GetUIpointer();
//    UImanager->ApplyCommand("/event/abort");

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

      virtual void PrepareNewEvent();

//    This method is called by G4StackManager at the begining of
//  each event.
//    Be careful that the urgentStack and the waitingStack of
//  G4StackManager are empty at this moment, because this method
//  is called before accepting primary particles. Also, note that
//  the postponeStack of G4StackManager may have some postponed
//  tracks.


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

};

#endif
