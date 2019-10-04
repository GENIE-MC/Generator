//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/FKR.h"

using std::endl;

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const FKR & parameters)
  {
     parameters.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
FKR::FKR()
{
  this->Reset();
}
//____________________________________________________________________________
FKR::~FKR()
{

}
//____________________________________________________________________________
void FKR::Print(ostream & stream) const
{
  stream << endl;
  stream << " lamda = " << Lamda   << endl;
  stream << " Tv    = " << Tv      << endl;
  stream << " Rv    = " << Rv      << endl;
  stream << " S     = " << S       << endl;
  stream << " Ta    = " << Ta      << endl;
  stream << " Ra    = " << Ra      << endl;
  stream << " B     = " << B       << endl;
  stream << " C     = " << C       << endl;
  stream << " R     = " << R       << endl;
  stream << " T     = " << T       << endl;
  stream << " T+    = " << Tplus   << endl;
  stream << " T-    = " << Tminus  << endl;
  stream << " R+    = " << Rplus   << endl;
  stream << " R-    = " << Rminus  << endl;
}
//____________________________________________________________________________
void FKR::Reset(void)
{
  Lamda        =  0.0;
  Tv           =  0.0;
  Rv           =  0.0;
  S            =  0.0;
  Ta           =  0.0;
  Ra           =  0.0;
  B            =  0.0;
  C            =  0.0;
  R            =  0.0;
  T            =  0.0;
  Tplus        =  0.0;
  Tminus       =  0.0;
  Rplus        =  0.0;
  Rminus       =  0.0;
}
//____________________________________________________________________________
