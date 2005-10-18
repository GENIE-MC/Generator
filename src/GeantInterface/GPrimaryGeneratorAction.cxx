//____________________________________________________________________________
/*!

\class   genie::geant::GPrimaryGeneratorAction

\brief   A GENIE/GEANT4 PrimaryGeneratorAction for driving the GENIE/GEANT4
         G4VPrimaryGenerator feeding GENIE events into GEANT4

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created Sepember 12, 2005

*/
//____________________________________________________________________________

#include <cassert>

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>

#include "GeantInterface/GPrimaryGeneratorAction.h"
#include "GeantInterface/GNtpPrimaryGenerator.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::geant;

//____________________________________________________________________________
GPrimaryGeneratorAction::GPrimaryGeneratorAction(GPrimaryGeneratorType_t pgt):
fPrimaryGeneratorType(pgt)
{
  this->Init();

  if(!fPrimaryGenerator) {
     LOG("GEANTi", pFATAL) << "No G4VPrimaryGenerator was instantiated";
  }
  assert(fPrimaryGenerator);
}
//____________________________________________________________________________
GPrimaryGeneratorAction::~GPrimaryGeneratorAction()
{

}
//____________________________________________________________________________
void GPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // ask the G4VPrimaryGenerator to generate (or simply read) an event
  fPrimaryGenerator->GeneratePrimaryVertex(event);
}
//____________________________________________________________________________
void GPrimaryGeneratorAction::Init(void)
{
  fPrimaryGenerator = 0;

  LOG("GEANTi", pINFO) 
            << "You requested: " 
            << GPrimaryGeneratorType::AsString(fPrimaryGeneratorType);

  switch (fPrimaryGeneratorType) {
     case kPGNtpRd:     
          fPrimaryGenerator = new GNtpPrimaryGenerator;
          break;

     case kPGEVGDrv:  
     case kPGMCJDrv:
          LOG("GEANTi", pERROR) 
             << "The requested G4VPrimaryGenerator is not supported yet!";
          break;

     case kPGUndefined: 
     default:   
          LOG("GEANTi", pERROR) 
                     << "You didn't request a valid G4VPrimaryGenerator!";
          break;
  }
}
//____________________________________________________________________________

