//
// ********************************************************************
// * License and Disclaimer *
// * *
// * The Geant4 software is copyright of the Copyright Holders of *
// * the Geant4 Collaboration. It is provided under the terms and *
// * conditions of the Geant4 Software License, included in the file *
// * LICENSE and available at http://cern.ch/geant4/license . These *
// * include a list of copyright holders. *
 // * *
 // * Neither the authors of this software system, nor their employing *
 // * institutes,nor the agencies providing financial support for this *
 // * work make any representation or warranty, express or implied, *
 // * regarding this software system or assume any liability for its *
 // * use. Please see the license in the file LICENSE and URL above *
 // * for the full disclaimer and the limitation of liability. *
 // * *
 // * This code implementation is the result of the scientific and *
 // * technical work of the GEANT4 collaboration. *
 // * By using, copying, modifying or distributing the software (or *
 // * any work based on the software) you agree to acknowledge its *
 // * use in resulting scientific publications, and indicate your *
 // * acceptance of all terms of the Geant4 Software license. *
 // ********************************************************************
 //
 // $Id: EmLivermorePolarizedPhysics.hh 66704 2013-01-10 18:20:17Z gunter $

 #ifndef EmLivermorePolarizedPhysics_h
 #define EmLivermorePolarizedPhysics_h 1

 #include "G4VPhysicsConstructor.hh"
 #include "globals.hh"

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 class EmLivermorePolarizedPhysics : public G4VPhysicsConstructor
 {
 public:
  EmLivermorePolarizedPhysics(G4int ver = 1);

  // obsolete
  // EmLivermorePolarizedPhysics(G4int ver, const G4String&);

  virtual ~EmLivermorePolarizedPhysics();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

 private:
  G4int verbose;
 };

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 #endif

/*

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
/// \file electromagnetic/TestEm14/include/PhysListEmLivermore.hh
/// \brief Definition of the PhysListEmLivermore class
//
// $Id: PhysListEmLivermore.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysListEmLivermore_h
#define PhysListEmLivermore_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysListEmLivermore : public G4VPhysicsConstructor
{
  public: 
    PhysListEmLivermore(const G4String& name = "livermore");
   ~PhysListEmLivermore();

  public: 
    // This method is dummy for physics
    virtual void ConstructParticle() {};
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


*/




