//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <sstream>

#include "Framework/Interaction/ProcessInfo.h"

using std::ostringstream;
using std::endl;

using namespace genie;

ClassImp(ProcessInfo)

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const ProcessInfo & proc)
 {
   proc.Print(stream);
   return stream;
 }
}
//____________________________________________________________________________
ProcessInfo::ProcessInfo() :
TObject()
{
  this->Reset();
}
//____________________________________________________________________________
ProcessInfo::ProcessInfo(
              ScatteringType_t sc_type, InteractionType_t  int_type) :
TObject(),
fScatteringType(sc_type),
fInteractionType(int_type)
{

}
//____________________________________________________________________________
ProcessInfo::ProcessInfo(const ProcessInfo & proc) :
TObject()
{
  this->Copy(proc);
}
//____________________________________________________________________________
ProcessInfo::~ProcessInfo()
{

}
//____________________________________________________________________________
void ProcessInfo::Reset(void)
{
  fScatteringType  = kScNull;
  fInteractionType = kIntNull;
}
//____________________________________________________________________________
bool ProcessInfo::IsQuasiElastic(void) const
{
  return (fScatteringType == kScQuasiElastic);
}
//____________________________________________________________________________
bool ProcessInfo::IsDarkMatterElastic(void) const
{
  return (fScatteringType == kScDarkMatterElastic);
}
//____________________________________________________________________________
bool ProcessInfo::IsSingleKaon(void) const
{
  return (fScatteringType == kScSingleKaon);
}
//____________________________________________________________________________
bool ProcessInfo::IsDeepInelastic(void) const
{
  return (fScatteringType == kScDeepInelastic);
}
//____________________________________________________________________________
bool ProcessInfo::IsDarkMatterDeepInelastic(void) const
{
  return (fScatteringType == kScDarkMatterDeepInelastic);
}
//____________________________________________________________________________
bool ProcessInfo::IsResonant(void) const
{
  return (fScatteringType == kScResonant);
}
//____________________________________________________________________________
bool ProcessInfo::IsCoherent(void) const
{
  return (fScatteringType == kScCoherent);
}
//____________________________________________________________________________
bool ProcessInfo::IsCoherentElas(void) const
{
  return (fScatteringType == kScCoherentElas);
}
//____________________________________________________________________________
bool ProcessInfo::IsElectronScattering(void) const
{
  return (fScatteringType == kScNuElectronElastic ||
          fScatteringType == kScInverseMuDecay ||
          fScatteringType == kScIMDAnnihilation);
}
//____________________________________________________________________________
bool ProcessInfo::IsNuElectronElastic(void) const
{
  return (fScatteringType == kScNuElectronElastic);
}
//____________________________________________________________________________
bool ProcessInfo::IsInverseMuDecay(void) const
{
  return (fScatteringType == kScInverseMuDecay);
}
//____________________________________________________________________________
bool ProcessInfo::IsIMDAnnihilation(void) const
{
  return (fScatteringType == kScIMDAnnihilation);
}
//____________________________________________________________________________
bool ProcessInfo::IsInverseBetaDecay(void) const
{
  return (fScatteringType == kScInverseBetaDecay);
}
//____________________________________________________________________________
bool ProcessInfo::IsGlashowResonance(void) const
{
  return (fScatteringType == kScGlashowResonance);
}
//____________________________________________________________________________
bool ProcessInfo::IsAMNuGamma(void) const
{
  return (fScatteringType == kScAMNuGamma);
}
//____________________________________________________________________________
bool ProcessInfo::IsMEC(void) const
{
  return (fScatteringType == kScMEC);
}
//____________________________________________________________________________
bool ProcessInfo::IsDiffractive(void) const
{
  return (fScatteringType == kScDiffractive);
}
//____________________________________________________________________________
bool ProcessInfo::IsEM(void) const
{
  return (fInteractionType == kIntEM);
}
//____________________________________________________________________________
bool ProcessInfo::IsWeak(void) const
{
  return ( this->IsWeakCC() || this->IsWeakNC() || this->IsWeakMix());
}
//____________________________________________________________________________
bool ProcessInfo::IsWeakCC(void) const
{
  return (fInteractionType == kIntWeakCC);
}
//____________________________________________________________________________
bool ProcessInfo::IsWeakNC(void) const
{
  return (fInteractionType == kIntWeakNC);
}
//____________________________________________________________________________
bool ProcessInfo::IsDarkMatter(void) const
{
  return (fInteractionType == kIntDarkMatter);
}
//____________________________________________________________________________
bool ProcessInfo::IsWeakMix(void) const
{
  return (fInteractionType == kIntWeakMix);
}
//____________________________________________________________________________
InteractionType_t ProcessInfo::InteractionTypeId(void) const
{
  return fInteractionType;
}
//____________________________________________________________________________
ScatteringType_t ProcessInfo::ScatteringTypeId(void) const
{
  return fScatteringType;
}
//____________________________________________________________________________
string ProcessInfo::AsString(void) const
{
  ostringstream stream;

  stream << "<"
         << this->ScatteringTypeAsString()
         << " - "
         << this->InteractionTypeAsString()
         << ">";

  return stream.str();
}
//____________________________________________________________________________
string ProcessInfo::ScatteringTypeAsString(void) const
{
  string scattering_type = ScatteringType::AsString(fScatteringType);
  return scattering_type;
}
//____________________________________________________________________________
string ProcessInfo::InteractionTypeAsString(void) const
{
  string interaction_type = InteractionType::AsString(fInteractionType);
  return interaction_type;
}
//____________________________________________________________________________
void ProcessInfo::Set(ScatteringType_t sc_type, InteractionType_t  int_type)
{
  fScatteringType  = sc_type;
  fInteractionType = int_type;
}
//____________________________________________________________________________
bool ProcessInfo::Compare(const ProcessInfo & proc) const
{
  return (
       fScatteringType  == proc.fScatteringType &&
       fInteractionType == proc.fInteractionType
  );
}
//____________________________________________________________________________
void ProcessInfo::Copy(const ProcessInfo & proc)
{
  fScatteringType  = proc.fScatteringType;
  fInteractionType = proc.fInteractionType;
}
//____________________________________________________________________________
void ProcessInfo::Print(ostream & stream) const
{
  stream << "[-] [Process-Info]  " << endl
         << " |--> Interaction : " << this->InteractionTypeAsString() << endl
         << " |--> Scattering  : " << this->ScatteringTypeAsString()  << endl;
}
//____________________________________________________________________________
bool ProcessInfo::operator == (const ProcessInfo & proc) const
{
  return this->Compare(proc);
}
//___________________________________________________________________________
ProcessInfo & ProcessInfo::operator = (const ProcessInfo & proc)
{
  this->Copy(proc);
  return (*this);
}
//____________________________________________________________________________

