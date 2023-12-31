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
//
// $Id: G4QEGamma.cc 67971 2019-10-24 gcosmo $
//
// 
// ---------------------------------------------------------------------
//  Paul Heyes, 2019
// ---------------------------------------------------------------------

 /*
QEGamma class is identical to the G4Gamma class. We need to create a 'special'
particle class in order to apply the correct physics -> in this case the 'standard'
Compton scattering or the quantum entangled Compton scattering. 
 */	

#include "G4QEGamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                          QEGAMMA                               ###
// ######################################################################
G4QEGamma* G4QEGamma::theQEGamma = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4QEGamma::G4QEGamma(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable,  G4bool              shortlived,
       const G4String&     subType,      G4int               antiEncoding)
 : G4ParticleDefinition( aName, mass, width, charge, iSpin, iParity,
           iConjugation, iIsospin, iIsospin3, gParity, pType,
           lepton, baryon, encoding, stable, lifetime, decaytable,
           shortlived, subType, antiEncoding)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4QEGamma::~G4QEGamma()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//    Arguments for constructor are as follows 
//               name             mass          width         charge
//             2*spin           parity  C-conjugation
//          2*Isospin       2*Isospin3       G-parity
//               type    lepton number  baryon number   PDG encoding
//             stable         lifetime    decay table 
//             shortlived      subType    anti_encoding
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4QEGamma* G4QEGamma::Definition() 
{  
  if(!theQEGamma) {
    theQEGamma = new G4QEGamma(
	        "QEgamma",    0.0*MeV,       0.0*MeV,         0.0, 
		            2,         -1,            -1,          
		            0,          0,             0,             
	          "boson",          0,             0,          22,
	             true,        0.0,          NULL,
                false,   "photon",        22);
  }
  return theQEGamma;
}

G4QEGamma*  G4QEGamma::QEGammaDefinition() 
{
  return Definition();
}

G4QEGamma*  G4QEGamma::QEGamma() 
{
  return Definition();
}
