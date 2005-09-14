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

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>

#include "GeantInterface/GPrimaryGeneratorAction.h"

using namespace genie;
using namespace genie::geant;

//____________________________________________________________________________
GPrimaryGeneratorAction::GPrimaryGeneratorAction() :
fPrimaryGeneratorType(pgt)
{

}
//____________________________________________________________________________
GPrimaryGeneratorAction::~GPrimaryGeneratorAction()
{

}
//____________________________________________________________________________
void GPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{

}
//____________________________________________________________________________
