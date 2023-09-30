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
//Brief Definition of the Hit class

#ifndef SDHit_h
#define SDHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// SDHit
///
/// It records:
/// - the scatterer ID
/// - the particle time
/// - the particle local and global position
/// - and more

class SDHit : public G4VHit
{
public:
  SDHit();
    SDHit(G4double t);
    SDHit(const SDHit &right);
    virtual ~SDHit();

    const SDHit& operator=(const SDHit &right);
    int operator==(const SDHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void*aHit);
    
    void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    void Print();
    
    G4int GetID() const { return fId; }

    void SetTime(G4double val) { fTime = val; }
    G4double GetTime() const { return fTime; }

    void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    G4ThreeVector GetLocalPos() const { return fLocalPos; }

    void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    G4ThreeVector GetWorldPos() const { return fWorldPos; }

    void SetLocalPos2(G4ThreeVector xyz) { fLocalPos2 = xyz; }
    G4ThreeVector GetLocalPos2() const { return fLocalPos2; }

    void SetWorldPos2(G4ThreeVector xyz) { fWorldPos2 = xyz; }
    G4ThreeVector GetWorldPos2() const { return fWorldPos2; }

    void SetMomentumDirection(G4ThreeVector xyz) { fMomDir = xyz; }
    G4ThreeVector GetMomentumDirection() const { return fMomDir; }

    void SetMomentumDirection2(G4ThreeVector xyz) { fMomDir2 = xyz; }
    G4ThreeVector GetMomentumDirection2() const { return fMomDir2; }

    void SetPolarization(G4ThreeVector xyz) { fPol = xyz; }
    G4ThreeVector GetPolarization() const { return fPol; }

    void SetEnergy(G4double ene) { fE = ene; }
    G4double GetEnergy() const { return fE; }

    void SetEnergy2(G4double ene) { fE2 = ene; }
    G4double GetEnergy2() const { return fE2; }

private:
    G4int fId;
    G4double fTime;
    G4ThreeVector fLocalPos;
    G4ThreeVector fWorldPos;
    G4ThreeVector fLocalPos2;
    G4ThreeVector fWorldPos2;
    G4ThreeVector fMomDir;
    G4ThreeVector fMomDir2;
    G4ThreeVector fPol;
    G4double fE;
    G4double fE2;

};

typedef G4THitsCollection<SDHit> SDHitsCollection;

extern G4ThreadLocal G4Allocator<SDHit>* HitAllocator;

inline void* SDHit::operator new(size_t)
{
    if (!HitAllocator) HitAllocator = new G4Allocator<SDHit>;
    return (void*)HitAllocator->MallocSingle();
}

inline void SDHit::operator delete(void*aHit)
{
    HitAllocator->FreeSingle((SDHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
