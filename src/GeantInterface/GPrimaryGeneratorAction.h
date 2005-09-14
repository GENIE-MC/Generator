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

#ifndef _G_PRIMARY_GENERATOR_ACTION_
#define _G_PRIMARY_GENERATOR_ACTION_

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VPrimaryGenerator.hh>
#include <globals.hh>

#include "GeantInterface/GPrimaryGeneratorType.h"

class G4Event;

namespace genie {
namespace geant {

class GPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  GPrimaryGeneratorAction(GPrimaryGeneratorType_t pgt);
  virtual ~GPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event * event);

private:
  G4VPrimaryGenerator *   fPrimaryGenerator;
  GPrimaryGeneratorType_t fPrimaryGeneratorType;
};

} // geant namespace
} // genie namespace

#endif // _G_PRIMARY_GENERATOR_ACTION_


