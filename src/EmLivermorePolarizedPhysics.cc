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
// $Id: G4EmLivermorePolarizedPhysics.cc 68750 2013-04-05 10:19:04Z gcosmo $

/*
Only minor adaptions are made to the G4EmLivermorePolarizedPhysics class:

- Header files QEComptonModel.hh and BoringComptonModel.hh, as well as GEgamma.hh are included.
- In the method ConstructParticle(), the QEgamma particle is added.
- In the method ConstructProcess(), the relevant process models are added depending on the particle type.
  For a QEgamma particle, an instance of the the QEComptonModel is registered as the applicale model within
  the G4ComptonScattering; for a standard gamma, the BoringComptonModel is registered.

To make things simpler, all other possible EM processes have been commented out, but can of course be switched
back on again.
*/

#include "EmLivermorePolarizedPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

// *** Processes and models

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePolarizedPhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermorePolarizedComptonModel.hh"
#include "QEComptonModel.hh"
#include "BoringComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermorePolarizedGammaConversionModel.hh"

#include "G4RayleighScattering.hh"
#include "G4LivermorePolarizedRayleighModel.hh"

// e+-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+
#include "G4eplusAnnihilation.hh"

// mu+-
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hPairProductionModel.hh"

// hadrons
#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"

#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4alphaIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

// msc models
#include "G4UrbanMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

// interfaces
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"
#include "G4UAtomicDeexcitation.hh"

// particles
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4QEGamma.hh"

//
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4VPhysicsConstructor.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(EmLivermorePolarizedPhysics);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmLivermorePolarizedPhysics::EmLivermorePolarizedPhysics(G4int ver)
 : G4VPhysicsConstructor("EmLivermorePolarizedPhysics"), verbose(ver)
{
 G4LossTableManager::Instance();
 SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
EmLivermorePolarizedPhysics::EmLivermorePolarizedPhysics(G4int ver, const String&)
 : G4VPhysicsConstructor("EmLivermorePolarizedPhysics"), verbose(ver)
{
 G4LossTableManager::Instance();
 SetPhysicsType(bElectromagnetic);
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmLivermorePolarizedPhysics::~EmLivermorePolarizedPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmLivermorePolarizedPhysics::ConstructParticle()
{
// gamma
 G4Gamma::Gamma();
 G4QEGamma::QEGamma();

// leptons
 G4Electron::Electron();
 G4Positron::Positron();
 G4MuonPlus::MuonPlus();
 G4MuonMinus::MuonMinus();

// mesons
 G4PionPlus::PionPlusDefinition();
 G4PionMinus::PionMinusDefinition();
 G4KaonPlus::KaonPlusDefinition();
 G4KaonMinus::KaonMinusDefinition();

// baryons
 G4Proton::Proton();
 G4AntiProton::AntiProton();

// ions
 G4Deuteron::Deuteron();
 G4Triton::Triton();
 G4He3::He3();
 G4Alpha::Alpha();
 G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmLivermorePolarizedPhysics::ConstructProcess()
{
 if(verbose > 1) {
 G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
 }
 G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

 /*
 // muon & hadron bremsstrahlung and pair production
 G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
 G4MuPairProduction* mup = new G4MuPairProduction();
 G4hBremsstrahlung* pib = new G4hBremsstrahlung();
 G4hPairProduction* pip = new G4hPairProduction();
 G4hBremsstrahlung* kb = new G4hBremsstrahlung();
 G4hPairProduction* kp = new G4hPairProduction();
 G4hBremsstrahlung* pb = new G4hBremsstrahlung();
 G4hPairProduction* pp = new G4hPairProduction();

 // muon & hadron multiple scattering
 G4MuMultipleScattering* mumsc = new G4MuMultipleScattering();
 mumsc->AddEmModel(0, new G4WentzelVIModel());
 G4MuMultipleScattering* pimsc = new G4MuMultipleScattering();
 pimsc->AddEmModel(0, new G4WentzelVIModel());
 G4MuMultipleScattering* kmsc = new G4MuMultipleScattering();
 kmsc->AddEmModel(0, new G4WentzelVIModel());
 G4MuMultipleScattering* pmsc = new G4MuMultipleScattering();
 pmsc->AddEmModel(0, new G4WentzelVIModel());
 G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");

 // high energy limit for e+- scattering models
 G4double highEnergyLimit = 100*MeV;

 // nuclear stopping
 G4NuclearStopping* ionnuc = new G4NuclearStopping();
 G4NuclearStopping* pnuc = new G4NuclearStopping();
*/

 // Add Livermore EM Processes
 auto theParticleIterator = GetParticleIterator(); 
 theParticleIterator->reset();

 while( (*theParticleIterator)() ){

 G4ParticleDefinition* particle = theParticleIterator->value();
 G4String particleName = particle->GetParticleName();
 
  // G4cout << "Particle iterated: " << particleName << G4endl;

 //Applicability range for Livermore models
 //for higher energies, the Standard models are used
 G4double LivermoreHighEnergyLimit = GeV;

 if (particleName == "QEgamma") {

 //G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
 //G4LivermorePolarizedPhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePolarizedPhotoElectricModel();
 //theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
 //ph->RegisterProcess(thePhotoElectricEffect, particle);

 G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
 
 // ********** Switch between models here *********** //

 //G4LivermorePolarizedComptonModel* theLivermoreComptonModel = new G4LivermorePolarizedComptonModel();
 //theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
 
 //BoringComptonModel* theBoringComptonModel = new BoringComptonModel();
 //theBoringComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //theComptonScattering->AddEmModel(0, theBoringComptonModel);

 QEComptonModel* theQECorrelatedComptonModel = new QEComptonModel();            // Line added
 theQECorrelatedComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);     // Line added
 theComptonScattering->AddEmModel(0, theQECorrelatedComptonModel);              // Line added

 // G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;

 ph->RegisterProcess(theComptonScattering, particle);
 
  // ********** Switch between models here *********** //

 //G4GammaConversion* theGammaConversion = new G4GammaConversion();
 //G4LivermorePolarizedGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermorePolarizedGammaConversionModel();
 //theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
 //ph->RegisterProcess(theGammaConversion, particle);

 //G4RayleighScattering* theRayleigh = new G4RayleighScattering();
 //G4LivermorePolarizedRayleighModel* theRayleighModel = new G4LivermorePolarizedRayleighModel();
 //theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //theRayleigh->AddEmModel(0, theRayleighModel);
 //ph->RegisterProcess(theRayleigh, particle);

 }
 
  else if (particleName == "gamma") {

 //G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
 //G4LivermorePolarizedPhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePolarizedPhotoElectricModel();
 //theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
 //ph->RegisterProcess(thePhotoElectricEffect, particle);

 G4ComptonScattering* theNormalComptonScattering = new G4ComptonScattering();
 
 // ********** Switch between models here *********** //

 //G4LivermorePolarizedComptonModel* theLivermoreComptonModel = new G4LivermorePolarizedComptonModel();
 //theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //theNormalComptonScattering->AddEmModel(0, theLivermoreComptonModel);
 
 BoringComptonModel* theBoringComptonModel = new BoringComptonModel();
 theBoringComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 theNormalComptonScattering->AddEmModel(0, theBoringComptonModel);

 //QEComptonModel* theQECorrelatedComptonModel = new QEComptonModel();            // Line added
 //theQECorrelatedComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);     // Line added
 //theComptonScattering->AddEmModel(0, theQECorrelatedComptonModel);              // Line added

 ph->RegisterProcess(theNormalComptonScattering, particle);
 
  // ********** Switch between models here *********** //

 //G4GammaConversion* theGammaConversion = new G4GammaConversion();
 //G4LivermorePolarizedGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermorePolarizedGammaConversionModel();
 //theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
 //ph->RegisterProcess(theGammaConversion, particle);

 //G4RayleighScattering* theRayleigh = new G4RayleighScattering();
 //G4LivermorePolarizedRayleighModel* theRayleighModel = new G4LivermorePolarizedRayleighModel();
 //theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
 //theRayleigh->AddEmModel(0, theRayleighModel);
 //ph->RegisterProcess(theRayleigh, particle);
 
}

 /*
 else if (particleName == "e-") {

 // multiple scattering
 G4eMultipleScattering* msc = new G4eMultipleScattering;
 msc->SetStepLimitType(fUseDistanceToBoundary);
 G4UrbanMscModel* msc1 = new G4UrbanMscModel();
 G4WentzelVIModel* msc2 = new G4WentzelVIModel();
 msc1->SetHighEnergyLimit(highEnergyLimit);
 msc2->SetLowEnergyLimit(highEnergyLimit);
 msc->AddEmModel(0, msc1);
 msc->AddEmModel(0, msc2);

 G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
 G4CoulombScattering* ss = new G4CoulombScattering();
 ss->SetEmModel(ssm, 1);
 ss->SetMinKinEnergy(highEnergyLimit);
 ssm->SetLowEnergyLimit(highEnergyLimit);
 ssm->SetActivationLowEnergyLimit(highEnergyLimit);
 ph->RegisterProcess(msc, particle);
 ph->RegisterProcess(ss, particle);

 // Ionisation
 G4eIonisation* eIoni = new G4eIonisation();
 G4LivermoreIonisationModel* theIoniLivermore = new
 G4LivermoreIonisationModel();
 theIoniLivermore->SetHighEnergyLimit(0.1*MeV);
 eIoni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );
 eIoni->SetStepFunction(0.2, 100*um); //
 ph->RegisterProcess(eIoni, particle);

 // Bremsstrahlung from standard
 G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
 ph->RegisterProcess(eBrem, particle);

 } else if (particleName == "e+") {

 // multiple scattering
 G4eMultipleScattering* msc = new G4eMultipleScattering;
 msc->SetStepLimitType(fUseDistanceToBoundary);
 G4UrbanMscModel* msc1 = new G4UrbanMscModel();
 G4WentzelVIModel* msc2 = new G4WentzelVIModel();
 msc1->SetHighEnergyLimit(highEnergyLimit);
 msc2->SetLowEnergyLimit(highEnergyLimit);
 msc->AddEmModel(0, msc1);
 msc->AddEmModel(0, msc2);

 G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
 G4CoulombScattering* ss = new G4CoulombScattering();
 ss->SetEmModel(ssm, 1);
 ss->SetMinKinEnergy(highEnergyLimit);
 ssm->SetLowEnergyLimit(highEnergyLimit);
 ssm->SetActivationLowEnergyLimit(highEnergyLimit);

 // Ionisation
 G4eIonisation* eIoni = new G4eIonisation();
 eIoni->SetStepFunction(0.2, 100*um);

 ph->RegisterProcess(msc, particle);
 ph->RegisterProcess(eIoni, particle);
 ph->RegisterProcess(new G4eBremsstrahlung(), particle);
 ph->RegisterProcess(new G4eplusAnnihilation(), particle);
 ph->RegisterProcess(ss, particle);

 } else if (particleName == "mu+" ||
 particleName == "mu-" ) {

 G4MuIonisation* muIoni = new G4MuIonisation();
 muIoni->SetStepFunction(0.2, 50*um);

 ph->RegisterProcess(mumsc, particle);
 ph->RegisterProcess(muIoni, particle);
 ph->RegisterProcess(mub, particle);
 ph->RegisterProcess(mup, particle);
 ph->RegisterProcess(new G4CoulombScattering(), particle);

 } else if (particleName == "alpha" ||
 particleName == "He3" ) {

 // Identical to G4EmStandardPhysics_option3

 G4hMultipleScattering* msc = new G4hMultipleScattering();
 G4ionIonisation* ionIoni = new G4ionIonisation();
 ionIoni->SetStepFunction(0.1, 10*um);

 ph->RegisterProcess(msc, particle);
 ph->RegisterProcess(ionIoni, particle);
 ph->RegisterProcess(ionnuc, particle);

 } else if (particleName == "GenericIon") {

 // Identical to G4EmStandardPhysics_option3

 G4ionIonisation* ionIoni = new G4ionIonisation();
 ionIoni->SetEmModel(new G4IonParametrisedLossModel());
 ionIoni->SetStepFunction(0.1, 1*um);

 ph->RegisterProcess(hmsc, particle);
 ph->RegisterProcess(ionIoni, particle);
 ph->RegisterProcess(ionnuc, particle);

 } else if (particleName == "pi+" ||
 particleName == "pi-" ) {

 //G4hMultipleScattering* pimsc = new G4hMultipleScattering();
 G4hIonisation* hIoni = new G4hIonisation();
 hIoni->SetStepFunction(0.2, 50*um);

 ph->RegisterProcess(pimsc, particle);
 ph->RegisterProcess(hIoni, particle);
 ph->RegisterProcess(pib, particle);
 ph->RegisterProcess(pip, particle);

 } else if (particleName == "kaon+" ||
 particleName == "kaon-" ) {

 //G4hMultipleScattering* kmsc = new G4hMultipleScattering();
 G4hIonisation* hIoni = new G4hIonisation();
 hIoni->SetStepFunction(0.2, 50*um);

 ph->RegisterProcess(kmsc, particle);
 ph->RegisterProcess(hIoni, particle);
 ph->RegisterProcess(kb, particle);
 ph->RegisterProcess(kp, particle);

 } else if (particleName == "proton" ||
 particleName == "anti_proton") {

 //G4hMultipleScattering* pmsc = new G4hMultipleScattering();
 G4hIonisation* hIoni = new G4hIonisation();
 hIoni->SetStepFunction(0.2, 50*um);

 ph->RegisterProcess(pmsc, particle);
 ph->RegisterProcess(hIoni, particle);
 ph->RegisterProcess(pb, particle);
 ph->RegisterProcess(pp, particle);
 ph->RegisterProcess(pnuc, particle);

 } else if (particleName == "B+" ||
 particleName == "B-" ||
 particleName == "D+" ||
 particleName == "D-" ||
 particleName == "Ds+" ||
 particleName == "Ds-" ||
 particleName == "anti_He3" ||
 particleName == "anti_alpha" ||
 particleName == "anti_deuteron" ||
 particleName == "anti_lambda_c+" ||
 particleName == "anti_omega-" ||
 particleName == "anti_sigma_c+" ||
 particleName == "anti_sigma_c++" ||
 particleName == "anti_sigma+" ||
 particleName == "anti_sigma-" ||
 particleName == "anti_triton" ||
 particleName == "anti_xi_c+" ||
 particleName == "anti_xi-" ||
 particleName == "deuteron" ||
 particleName == "lambda_c+" ||
 particleName == "omega-" ||
 particleName == "sigma_c+" ||
 particleName == "sigma_c++" ||
 particleName == "sigma+" ||
 particleName == "sigma-" ||
 particleName == "tau+" ||
 particleName == "tau-" ||
 particleName == "triton" ||
 particleName == "xi_c+" ||
 particleName == "xi-" ) {

 // Identical to G4EmStandardPhysics_option3

 ph->RegisterProcess(hmsc, particle);
 ph->RegisterProcess(new G4hIonisation(), particle);
 ph->RegisterProcess(pnuc, particle);
 }
 */
 }

 // Em options
 //
 G4EmProcessOptions opt;
 opt.SetVerbose(verbose);

 // Multiple Coulomb scattering
 //
 opt.SetPolarAngleLimit(CLHEP::pi);

 // Physics tables
 //

 opt.SetMinEnergy(10*eV);
 opt.SetMaxEnergy(10*TeV);
 opt.SetDEDXBinning(220);
 opt.SetLambdaBinning(220);

 // Nuclear stopping
 //pnuc->SetMaxKinEnergy(MeV);

 // Ionization
 //
 //opt.SetSubCutoff(true);

 // Deexcitation
 //
 //G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
 //G4LossTableManager::Instance()->SetAtomDeexcitation(de);
 //de->SetFluo(true);
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
