//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - December 08, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include "BaryonResonance/BaryonResUtils.h"
#include "Interaction/XclsTag.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

using std::endl;
using std::ostringstream;

using namespace genie;
using namespace genie::utils;

ClassImp(XclsTag)

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const XclsTag & xcls)
 {
   xcls.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
XclsTag::XclsTag() :
TObject()
{
  this->Reset();
}
//___________________________________________________________________________
XclsTag::XclsTag(const XclsTag & xcls) :
TObject()
{
  this->Reset();
  this->Copy(xcls);
}
//___________________________________________________________________________
XclsTag::~XclsTag()
{

}
//___________________________________________________________________________
bool XclsTag::IsInclusiveCharm(void) const
{
  return ( this->IsCharmEvent() && (this->CharmHadronPdg() == 0) );
}
//___________________________________________________________________________
void XclsTag::SetCharm(int charm_pdgc)
{
  fIsCharmEvent     = true;
  fCharmedHadronPdg = charm_pdgc; // leave as 0 (default) for inclusive charm
}
//___________________________________________________________________________
void XclsTag::UnsetCharm(void)
{
  fIsCharmEvent     = false;
  fCharmedHadronPdg = 0;
}
//___________________________________________________________________________
void XclsTag::SetNPions(int npi_plus, int npi_0, int npi_minus)
{
  fNPiPlus  = npi_plus;
  fNPi0     = npi_0;
  fNPiMinus = npi_minus;
}
//___________________________________________________________________________
void XclsTag::SetNNucleons(int np, int nn)
{
  fNProtons  = np;
  fNNeutrons = nn;
}
//___________________________________________________________________________
void XclsTag::ResetNPions(void)
{
  fNPi0     = 0;
  fNPiPlus  = 0;
  fNPiMinus = 0;
}
//___________________________________________________________________________
void XclsTag::ResetNNucleons(void)
{
  fNProtons  = 0;
  fNNeutrons = 0;
}
//___________________________________________________________________________
void XclsTag::SetResonance(Resonance_t res)
{
  fResonance = res;
}
//___________________________________________________________________________
void XclsTag::Reset(void)
{
  fIsCharmEvent     = false;
  fCharmedHadronPdg = 0;
  fNProtons         = 0;
  fNNeutrons        = 0;
  fNPi0             = 0;
  fNPiPlus          = 0;
  fNPiMinus         = 0;
  fResonance        = kNoResonance;
}
//___________________________________________________________________________
void XclsTag::Copy(const XclsTag & xcls)
{
  fIsCharmEvent     = xcls.fIsCharmEvent;
  fCharmedHadronPdg = xcls.fCharmedHadronPdg;
  fNProtons         = xcls.fNProtons;
  fNNeutrons        = xcls.fNNeutrons;
  fNPi0             = xcls.fNPi0;
  fNPiPlus          = xcls.fNPiPlus;
  fNPiMinus         = xcls.fNPiMinus;
  fResonance        = xcls.fResonance;
}
//___________________________________________________________________________
string XclsTag::AsString(void) const
{
// codifies XclsTag state into a compact string

  ostringstream tag;

  bool need_separator = false;

  if(fIsCharmEvent) {
    tag << "charm:";
    if(fCharmedHadronPdg) tag << fCharmedHadronPdg;
    else tag << "incl";
    need_separator = true;
  }

  bool multset = 
       fNProtons>0 || fNNeutrons>0 || fNPiPlus>0 || fNPiMinus>0 || fNPi0>0;
  if(multset) {
    if(need_separator) tag << ";";
    tag << "hmult:"
        << "(p=" << fNProtons << ",n=" << fNNeutrons
        << ",pi+=" << fNPiPlus << ",pi-=" << fNPiMinus << ",pi0=" << fNPi0 
        << ")";
  }

  if(this->KnownResonance()) {
    if(need_separator) tag << ";";
    tag << "res:" << fResonance;
  }

  return tag.str();
}
//___________________________________________________________________________
void XclsTag::Print(ostream & stream) const
{
  stream << "[-] [Exclusive Process Info] " << endl;

  stream << " |--> charm prod.  : "
         << utils::print::BoolAsString(fIsCharmEvent);
  if(fIsCharmEvent) {
     if(!fCharmedHadronPdg) stream << " [inclusive]";
     else  {
       stream << " - Charm hadron PDG-code = " << fCharmedHadronPdg;

       TParticlePDG * chadr = PDGLibrary::Instance()->Find( fCharmedHadronPdg );
       if(chadr)
           stream << " (" << chadr->GetName() << ")";
     }
  }
  stream << endl;

  stream << " |--> f/s nucleons :"
         << " N(p) = "    << fNProtons
         << " N(n) = "    << fNNeutrons
         << endl;

  stream << " |--> f/s pions    :"
         << " N(pi^0) = "    << fNPi0
         << " N(pi^+) = "    << fNPiPlus
         << " N(pi^-) = "    << fNPiMinus
         << endl;

  stream << " |--> resonance    : ";
  if(this->KnownResonance()) {
     stream << res::AsString(fResonance);
  } else {
     stream << "[not set]";
  }
  stream << endl;
}
//___________________________________________________________________________
XclsTag & XclsTag::operator = (const XclsTag & xcls)
{
  this->Copy(xcls);
  return (*this);
}
//___________________________________________________________________________


