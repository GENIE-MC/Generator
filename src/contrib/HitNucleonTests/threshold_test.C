#include <iostream>
#include <memory>

#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "Framework/Interaction/Interaction.h"

#include "Math/Vector4D.h"
#include "Math/Boost.h" 


using namespace genie;
using namespace std;

void threshold_test() {


  unique_ptr<Interaction> i_p {Interaction::QELCC( 1000180400, kPdgNeutron, kPdgNuMu, 2. )} ;

  auto th = i_p -> PhaseSpace().Threshold();

  cout << th << endl;

  TLorentzVector on_shell_target;
  on_shell_target.SetXYZM( 0.1, -0.2, 0, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() );
    
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( on_shell_target ) ;

  i_p -> InitState().Tgt().HitNucP4().Print();

  th = i_p -> PhaseSpace().Threshold();

  cout << th << endl;

  

}
