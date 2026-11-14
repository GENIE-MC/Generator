//____________________________________________________________________________
/*
Dark Matter model parameters calculated in the same form as
the Feynmann-Kislinger-Ravndall (FKR) baryon excitation model parameters.

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/BoostedDarkMatter/XSection/FKRDM.h"

using std::endl;

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const FKRDM & parameters)
  {
     parameters.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
FKRDM::FKRDM()
{
  this->Reset();
}
//____________________________________________________________________________
FKRDM::~FKRDM()
{

}
//____________________________________________________________________________
void FKRDM::Print(ostream & stream) const
{
  stream << endl;
  stream << " lamda = " << Lamda   << endl;
  stream << " Tv    = " << Tv      << endl;
  stream << " Rv    = " << Rv      << endl;
  stream << " S     = " << S       << endl;
  stream << " Ta    = " << Ta      << endl;
  stream << " Ra    = " << Ra      << endl;
  stream << " Bs    = " << Bs      << endl;
  stream << " Cs    = " << Cs      << endl;
  stream << " Bz    = " << Bz      << endl;
  stream << " Cz    = " << Cz      << endl;
}
//____________________________________________________________________________
void FKRDM::Reset(void)
{
  Lamda        =  0.0;
  Tv           =  0.0;
  Rv           =  0.0;
  S            =  0.0;
  Ta           =  0.0;
  Ra           =  0.0;
  Bs           =  0.0;
  Cs           =  0.0;
  Bz           =  0.0;
  Cz           =  0.0;
}
//____________________________________________________________________________
