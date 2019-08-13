//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Igor Kakorin <kakorin@jinr.ru>
         Joint Institute for Nuclear Research  Nov 12, 2019
         

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 November 12, 2019 - IK
 Extended for MK model

*/
//____________________________________________________________________________
#include "Physics/Resonance/XSection/FKR_MK.h"

using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const FKR_MK & parameters)
  {
     parameters.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
FKR_MK::FKR_MK() : FKR()
{
  this->Reset();
}
//____________________________________________________________________________
FKR_MK::~FKR_MK()
{

}
//____________________________________________________________________________
void FKR_MK::Print(ostream & stream) const
{
  FKR::Print(stream);
  stream << " S+    = " << Splus   << endl;
  stream << " B+    = " << Bplus   << endl;
  stream << " C+    = " << Cplus   << endl;
  stream << " S-    = " << Sminus  << endl;
  stream << " B-    = " << Bminus  << endl;
  stream << " C-    = " << Cminus  << endl;

}
//____________________________________________________________________________
void FKR_MK::Reset(void)
{
  FKR::Reset();
  Splus        =  0.0;
  Bplus        =  0.0;
  Cplus        =  0.0;
  Sminus       =  0.0;
  Bminus       =  0.0;
  Cminus       =  0.0;

}
//____________________________________________________________________________
