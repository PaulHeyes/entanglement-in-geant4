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

//Paul Heyes, 2019
//Brief Implementation of the EventAction class

/*
Aside from the constructor and destructor, the two methods in this class are called during event-processing.

The 'BeginOfEventAction()' method sets the ID's of the individual sensitive detectors in order to collect 
the hits recorded in each of them. 

In the 'EndOfEventAction()' method the hits collections from the aforementioned sensitive detectors are 
retrieved. Using the G4AnalysisManager the ntuples and histograms are filled with the collected data. Here 
1d histograms are filled with the number of hits in each of the detectors. (Note: There are multiple hits per 
one event in some cases because the photons may scatter multiple times in the each detector and each scatterer! 
This needs to be put in consideration when working out the scattering angle: The scatterers may have to be of 
smaller dimensions and then the first hit in time in the detectors used to compute the angle.) The 2d histograms
are filled with the local (with respect to the detector) x and y coordinates of the hits. The ntuples are filled 
with the hits and their times in all detectors and scatterers.
Some of the data is also currently printed directly to the console.
*/

#include "EventAction.hh"
#include "SDHit.hh"
#include "BellsAnalysis.hh"
#include "EntangledGeneratorAction.hh"

#include "G4Event.hh"
#include "QERunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"

#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Constructor
EventAction::EventAction()
: G4UserEventAction(), 
  fSHC1ID(-1),
  fSHC2ID(-1)
  //psi()
{
  // set printing per each event
  QERunManager::GetRunManager()->SetPrintProgress(1);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Destructor
EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    if (fSHC1ID==-1) {
      fSHC1ID 
        = G4SDManager::GetSDMpointer()->GetCollectionID("ScatterDet1/scattererColl");
    }


    if (fSHC2ID==-1) {
      fSHC2ID 
        = G4SDManager::GetSDMpointer()->GetCollectionID("ScatterDet2/scattererColl");
    }
    

}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) 
    {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found.\n"; 
        G4Exception("EventAction::EndOfEventAction()",
                    "Code001", JustWarning, msg);
        return;
    }   


    // Get hits collections
    SDHitsCollection* sHC1
      = static_cast<SDHitsCollection*>(hce->GetHC(fSHC1ID));
      
    SDHitsCollection* sHC2
      = static_cast<SDHitsCollection*>(hce->GetHC(fSHC2ID));

  
    //    if ( (!sHC1) || (!sHC2) )
    // {
    //    G4ExceptionDescription msg;
    //    msg << "Some of hits collections of this event not found.\n"; 
    //    G4Exception("EventAction::EndOfEventAction()",
    //                "Code001", JustWarning, msg);
    //    return;
    // }   
    

    // Fill histograms & ntuple
    
    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
    // Fill histograms

    G4int ID = event->GetEventID();

    G4int n_hit1 = sHC1->entries();
    G4int n_hit2 = sHC2->entries();

    extern G4double QEPhi_ref;
    // extern G4double QEPhi;
    // extern G4double QETheta;


    if(ID%2==0) // All even events (incl. 0)
    {

        //G4cout << "||||||||| " << event->GetEventID() << " ||||||||||" << G4endl;

        //Detector1
        //G4int n_hit1 = dHC1->entries();                 //Number of hits
        //analysisManager->FillH1(0, n_hit1);

        if(n_hit1 > 0)
        {
           G4cout << "**** ScatterDetector1 Hit **** #" << n_hit1 << G4endl;

        // Here we are solely interested in the first ScatterDetector hit (that of argument '0') for now
        // indeed, this is true when only Compton scattering is simulated!
            
           SDHit* hit1 = (*sHC1)[1];
           G4ThreeVector momentumDirection1 = hit1->GetMomentumDirection();
           G4ThreeVector pos1 = hit1->GetWorldPos();
           G4ThreeVector pos2 = hit1->GetWorldPos2();
           
           // G4ThreeVector Pol1 = hit1->GetPolarization();

           G4double x1 = pos1.x();
           G4double y1 = pos1.y();
           G4double z1 = pos1.z();
           
           G4double x2 = pos2.x();
           G4double y2 = pos2.y();
           G4double z2 = pos2.z();
           
           // G4cout << "x1,y1,z1 = " << x1 << ", " << y1 << ", " << z1 << G4endl;
           // G4cout << "x2,y2,z2 = " << x2 << ", " << y2 << ", " << z2 << G4endl;
           
           G4double x = x2-x1;
           G4double y = y2-y1;
           G4double z = std::abs(z2-z1);           

           // G4double x = momentumDirection1.x();
           // G4double y = momentumDirection1.y();
           // G4double z = momentumDirection1.z();

           //analysisManager->FillH2(0, x, y);          //Fill 2D-histogram with ScatterDetector1-hits (used for polarization instead...)

		   ///
		   G4double r = std::sqrt(x*x+y*y+z*z);
		   //G4double phi = std::acos(x/r);
		   //G4double validPhi = std::asin(y/r);
		   
		   G4double phi = std::atan2(y/r,x/r); 	
		   ///

           G4double theta = std::acos(z/r);

           G4cout << "x,y,z = " << x << ", " << y << ", " << z << "  phi = " << phi << "  theta = " << theta << G4endl;
        
  
           if(phi != 0) analysisManager->FillH1(0, phi);
            analysisManager->FillH1(1, theta);

           QEPhi_ref = phi;

        }

     }

     else
     {

        //Detector2
        //G4int n_hit2 = dHC2->entries();
        //analysisManager->FillH1(1, n_hit2);

        if(n_hit2 > 0)
        {
           G4cout << "**** ScatterDetector2 Hit **** #" << n_hit2 << G4endl;
           
           // I think there should be a loop over the hit-IDs until scattering angle is no longer 0 (only if more physics processes are enabled)
        
           SDHit* hit2 = (*sHC2)[1];
           G4ThreeVector momentumDirection2 = hit2->GetMomentumDirection();
           G4ThreeVector pos1 = hit2->GetWorldPos();
           G4ThreeVector pos2 = hit2->GetWorldPos2();
           
           G4double x1 = pos1.x();
           G4double y1 = pos1.y();
           G4double z1 = pos1.z();  

		   G4double x2 = pos2.x();
           G4double y2 = pos2.y();
           G4double z2 = pos2.z();
           
           // G4cout << "x1,y1,z1 = " << x1 << ", " << y1 << ", " << z1 << G4endl;
           // G4cout << "x2,y2,z2 = " << x2 << ", " << y2 << ", " << z2 << G4endl;
           
           G4double x = x2-x1;
           G4double y = y2-y1;
           G4double z = std::abs(z2-z1);     
           
           // G4double x = momentumDirection2.x();
           // G4double y = momentumDirection2.y();
           // G4double z = momentumDirection2.z();

           //analysisManager->FillH2(1, x, y);      //Fill 2D-histogram with ScatterDetector2-hits (used for polarization instead...)


		   ///
		   G4double r = std::sqrt(x*x+y*y+z*z);
		   //G4double phi2 = std::acos(x/r);
		   //G4double validPhi2 = std::asin(y/r);
		   
		   G4double phi2 = std::atan2(y/r,x/r); 	
		   ///

           // z = -z;
           G4double theta2 = std::acos(z/r);

           //G4cout << "phi2 = " << phi2 << G4endl;
           G4cout << "x,y,z = " << x << ", " << y <<  ", " << z <<"  phi2 = " << phi2 << "  theta2 = " << theta2 << G4endl;
        
           if(phi2 != 0) analysisManager->FillH1(2, phi2);
           analysisManager->FillH1(3, theta2);

           G4double delta_phi = std::abs(QEPhi_ref-phi2);

           /// delta_phi handling: must be between 0 and pi
           if(delta_phi > pi) delta_phi = std::abs(delta_phi - 2*pi);
           
           G4cout << "delta_phi = " << delta_phi << G4endl;
           

           // if(QEPhi != 0){
            analysisManager->FillH1(4,delta_phi);
           // }          

        }

    }

  
    // Fill ntuple

    // ScatterDetector1Hits
    analysisManager->FillNtupleIColumn(2, sHC1->entries());

    // ScatterDetector2Hits
    analysisManager->FillNtupleIColumn(3, sHC2->entries());
    
     
    // Time 2
    for (G4int i=0;i<sHC2->entries();i++)
    {
      analysisManager->FillNtupleDColumn(4,(*sHC2)[i]->GetTime());
    }
 
  
    analysisManager->AddNtupleRow();  
    
    //
    // Print diagnostics: UI command /run/printProgress can be used
    // to set frequency of how often info should be dumpled
    // 
    
    G4int printModulo = QERunManager::GetRunManager()->GetPrintProgress();
    if ( printModulo==0 || event->GetEventID() % printModulo != 0) return;
    
    G4PrimaryParticle* primary = event->GetPrimaryVertex(0)->GetPrimary(0);
    G4cout << G4endl
           << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
           << primary->GetG4code()->GetParticleName()
           << " " << primary->GetMomentum() << G4endl;
    
/*

    //ScatterDetector 1
    //G4int n_hit = sHC1->entries();
    G4cout << "Scatterer 1 has " << n_hit1 << " hits." << G4endl;
    for (G4int i=0;i<n_hit1;i++)
    {
        SDHit* hit1 = (*sHC1)[i];
        hit1->Print();
    }


    //ScatterDetector 2
    //n_hit = sHC2->entries();
    G4cout << "Scatterer 2 has " << n_hit2 << " hits." << G4endl;
    for (G4int i=0;i<n_hit2;i++)
    {
        SDHit* hit2 = (*sHC2)[i];
        hit2->Print();
    }

*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
