//____________________________________________________________________________
/*!

\class    genie::ProcessInfo

\brief    A class encapsulating an enumeration of interaction types (EM,
          Weak-CC, Weak-NC) and scattering types (Elastic, Quasi Elastic,
          Deep Inelastic, Resonant Single Pion Production, Coherent Pion
          Production).

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004
 
*/
//____________________________________________________________________________

#include <sstream>
#include <string>

#include "Interaction/ProcessInfo.h"

using std::ostringstream;
using std::string;
using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const ProcessInfo & proc)
 {
   proc.Print(stream);

   return stream;
 }
}
//____________________________________________________________________________
ProcessInfo::ProcessInfo() :
fScatteringType  (kScNull ),
fInteractionType (kIntNull)
{

}
//____________________________________________________________________________
ProcessInfo::ProcessInfo(ScatteringType_t sc_type, InteractionType_t  int_type) :
fScatteringType  (sc_type),
fInteractionType (int_type)
{

}
//____________________________________________________________________________
ProcessInfo::ProcessInfo(const ProcessInfo & proc)
{
  fScatteringType  = proc.fScatteringType;
  fInteractionType = proc.fInteractionType;
}
//____________________________________________________________________________
ProcessInfo::~ProcessInfo()
{

}
//____________________________________________________________________________
bool ProcessInfo::IsElastic(void) const
{
  return (fScatteringType == kScElastic);
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
const char * ProcessInfo::AsString(void) const
{
  ostringstream stream;
  
  stream << "<" << ScatteringTypeAsString() << " - "
                                          << InteractionTypeAsString() << ">";

  return stream.str().c_str();
}
//____________________________________________________________________________
const char * ProcessInfo::ScatteringTypeAsString(void) const
{  
  string scattering_type = ScatteringType::AsString(fScatteringType);

  return scattering_type.c_str();
}
//____________________________________________________________________________
const char * ProcessInfo::InteractionTypeAsString(void) const
{
  string interaction_type = InteractionType::AsString(fInteractionType);

  return interaction_type.c_str();
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
void ProcessInfo::Print(ostream & stream) const
{
  stream << "[-] [Process-Info]  " << endl
         << " |--> Interaction : " << InteractionTypeAsString() << endl
         << " |--> Scattering  : " << ScatteringTypeAsString()  << endl;
}
//____________________________________________________________________________

