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
//Brief Implementation of the BellsRunAction class

/*
In the function 'GenerateRun()' an instance of the BellsRun class is created (i.e. new run is generated).

The 'BeginOfRunAction()' method prepares Ntuples and Histograms and opens a file ('BellsData') for saving 
the data during the run. Ntuples are created for saving the Detector and Scatterer hits, as well as the 
time of hits in Scatterer 2. 1d histograms are created to display the number of hits per run (this is no 
longer very useful data for our purpose). 2d histograms are created to show the local positions of the hits 
in both detectors respectively.

Both of the aforementioned methods are called by the RunManager during the RunInitialization (i.e. before 
the run).

The 'EndOfRunAction()' method is called at the end of a run. Output of information at the end of the run is 
possible here, mainly however the collected data is written to the file, which is subsequently closed.
*/

#include "BellsRunAction.hh"
#include "BellsPrimaryGeneratorAction.hh"
#include "BellsRun.hh"
#include "BellsAnalysis.hh"

#include "G4PhysicalConstants.hh"
#include "G4Run.hh"
#include "QERunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Constructor
BellsRunAction::BellsRunAction()
 : G4UserRunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Destructor
BellsRunAction::~BellsRunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* BellsRunAction::GenerateRun()
{ return new BellsRun; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BellsRunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetVerboseLevel(1);
  //man->SetFirstNtupleId(1); 
  //man->SetFirstHistoId(1);  

  //Create an n-column ntuple
  man->CreateNtuple("Bell", "Hits");
  man->CreateNtupleDColumn("NPhiNinety");
  man->CreateNtupleDColumn("NPhiZero");
  man->CreateNtupleIColumn("ScatterDet1Hits");		//Integer Ntuple
  man->CreateNtupleIColumn("ScatterDet2Hits");
  man->CreateNtupleDColumn("Scatter2Time");             //Double Ntuple
  man->FinishNtuple();

  //Create a 1d Histo
  man->CreateH1("Phi1", "Angle", 180, -pi, pi);         //G4int CreateH1(const G4String& name, const G4String& title, G4int nbins, G4double xmin, G4double xmax)
  man->CreateH1("Theta1", "Angle", 36, 0, pi);
  man->CreateH1("Phi2", "Angle", 180, -pi, pi);
  man->CreateH1("Theta2", "Angle", 36, 0, pi);
  man->CreateH1("delta_phi", "N", 180, 0, pi);
  man->CreateH1("N(delta_phi=90)", "N", 36, 0, pi);
  man->CreateH1("delta_phi=0", "N", 72, 0, 2*pi);
  man->CreateH1("N(delta_phi=0)", "N", 36, 0, pi);
  man->CreateH1("Phi_Distribution", "Angle", 180, 0, 2*pi);
  man->CreateH1("QEPhi_Distribution", "Angle", 180, 0, 2*pi);
  man->CreateH1("delta_phi_ref", "Angle", 180, 0, pi);

  //Create a 2d Histo
  man->CreateH2("Initial Polarization XY", "X vs Y", 100, -2, 2, 100, -2, 2);
  man->CreateH2("QE Polarization XY", "X vs Y", 100, -2, 2, 100, -2, 2);

  //Create a new output file
  man->OpenFile("BellsData");
 
  //Inform the runManager to save random number seed
  QERunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BellsRunAction::EndOfRunAction(const G4Run* run)
{
  //retrieve the number of events produced in the run
  G4int nofEvents = run->GetNumberOfEvent();

  //do nothing, if no events were processed
  if (nofEvents == 0) return;
  
  const BellsPrimaryGeneratorAction* generatorAction = static_cast<const BellsPrimaryGeneratorAction*>(
        QERunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 

  G4String partName;
  if (generatorAction) 
  { G4ParticleDefinition* particle = generatorAction -> GetParticleGun() -> GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
  
  //results
  //
  //const BellsRun* bellsRun = static_cast<const BellsRun*>(run);
  //G4int nbGoodEvents = bellsRun->GetNbGoodEvents();
  //G4double sumDose   = bellsRun->GetSumDose();
        
  //print
  //
  if (IsMaster())
  {
    G4cout
     << "\n--------------------End of Global Run-----------------------"
     << " \n The run was " << nofEvents << " events ";
  }
  else
  {
    G4cout
     << "\n--------------------End of Local Run------------------------"
     << " \n The run was " << nofEvents << " "<< partName;
  }


  //save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
