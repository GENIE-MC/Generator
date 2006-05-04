//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include "Interaction/ProcessInfo.h"

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
ProcessInfo::ProcessInfo()
{
  this->Reset();
}
//____________________________________________________________________________
ProcessInfo::ProcessInfo(
                      ScatteringType_t sc_type, InteractionType_t  int_type) :
fScatteringType  (sc_type),
fInteractionType (int_type)
{

}
//____________________________________________________________________________
ProcessInfo::ProcessInfo(const ProcessInfo & proc)
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
bool ProcessInfo::IsNuElectronElastic(void) const
{
  return (fScatteringType == kScNuElectronElastic);
}
//____________________________________________________________________________
bool ProcessInfo::IsQuasiElastic(void) const
{
  return (fScatteringType == kScQuasiElastic);
}
//____________________________________________________________________________
bool ProcessInfo::IsDeepInelastic(void) const
{
  return (fScatteringType == kScDeepInelastic);
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
bool ProcessInfo::IsInverseMuDecay(void) const
{
  return (fScatteringType == kScInverseMuDecay);
}
//____________________________________________________________________________
bool ProcessInfo::IsEM(void) const
{
  return (fInteractionType == kIntEM);
}
//____________________________________________________________________________
bool ProcessInfo::IsWeak(void) const
{
  return ( IsWeakCC() || IsWeakNC() );
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

