//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResParams.h"
#include "BaryonResonance/BaryonResDataSetI.h"

using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator<<(ostream & stream, const BaryonResParams & res_params)
  {
     res_params.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
BaryonResParams::BaryonResParams()
{
  fDataSet = 0;
  InitParams();
}
//____________________________________________________________________________
BaryonResParams::BaryonResParams(const BaryonResParams & res_params)
{
  fDataSet           = res_params.fDataSet;
  
  fResonanceIndex    = res_params.fResonanceIndex;
  fOrbitalAngularMom = res_params.fOrbitalAngularMom;
  fIsDeltaResonance  = res_params.fIsDeltaResonance;
  fIsNResonance      = res_params.fIsNResonance;
  fMass              = res_params.fMass;
  fWidth             = res_params.fWidth;
  fBreitWignerNorm   = res_params.fBreitWignerNorm;
}
//____________________________________________________________________________
BaryonResParams::~BaryonResParams()
{

}
//____________________________________________________________________________
void BaryonResParams::SetDataSet(const BaryonResDataSetI * data_set)
{
  InitParams();

  fDataSet = data_set;  
}
//____________________________________________________________________________
void BaryonResParams::RetrieveData(Resonance_t resonance)
{
  fResonanceIndex    = fDataSet -> ResonanceIndex    (resonance);
  fOrbitalAngularMom = fDataSet -> OrbitalAngularMom (resonance);
  fIsDeltaResonance  = fDataSet -> IsDeltaResonance  (resonance);
  fIsNResonance      = fDataSet -> IsNResonance      (resonance);
  fMass              = fDataSet -> Mass              (resonance);
  fWidth             = fDataSet -> Width             (resonance);
  fBreitWignerNorm   = fDataSet -> BreitWignerNorm   (resonance);
}
//____________________________________________________________________________
void BaryonResParams::Print(ostream & stream) const
{
  stream << "ResonanceIndex    = " << fResonanceIndex     << endl;
  stream << "OrbitalAngularMom = " << fOrbitalAngularMom  << endl;
  stream << "IsDeltaResonance  = " << fIsDeltaResonance   << endl;
  stream << "IsNResonance      = " << fIsNResonance       << endl;
  stream << "Mass              = " << fMass               << endl;
  stream << "Width             = " << fWidth              << endl;
  stream << "BreitWignerNorm   = " << fBreitWignerNorm    << endl;
}
//____________________________________________________________________________
void BaryonResParams::InitParams(void)
{
  fResonanceIndex    = 0;
  fOrbitalAngularMom = 0;
  fIsDeltaResonance  = false;
  fIsNResonance      = false;
  fMass              = 0.0;
  fWidth             = 0.0;
  fBreitWignerNorm   = 0.0;
}
//____________________________________________________________________________

