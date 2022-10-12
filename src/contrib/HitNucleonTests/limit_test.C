
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "Math/Vector4D.h"


void limit_test() {

  using namespace genie;

  // neutrino case on shell
  ROOT::Math::PxPyPzEVector probe( 0., 0., 2., 2. );

  double proton_mass = PDGLibrary::Instance()->Find(kPdgProton)->Mass();

  ROOT::Math::PxPyPzMVector nucleon( 0.1, 0.2, 0.1, proton_mass );
  
  std::cout << probe ;
  std::cout << nucleon ;
  

}
