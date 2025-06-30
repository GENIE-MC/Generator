//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <sstream>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Interaction/XclsTag.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

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
bool XclsTag::IsInclusiveStrange(void) const
{
  return ( this->IsStrangeEvent() && (this->StrangeHadronPdg() == 0) );
}
//___________________________________________________________________________
void XclsTag::SetStrange(int strange_pdgc)
{
  fIsStrangeEvent     = true;
  fStrangeHadronPdg   = strange_pdgc; // leave as 0 (default) for inclusive strange
}
//___________________________________________________________________________
void XclsTag::UnsetStrange(void)
{
  fIsStrangeEvent     = false;
  fStrangeHadronPdg   = 0;
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
void XclsTag::SetNRhos(int nrho_plus, int nrho_0, int nrho_minus)
{
  fNRhoPlus  = nrho_plus;
  fNRho0     = nrho_0;
  fNRhoMinus = nrho_minus;
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
void XclsTag::ResetNRhos(void)
{
  fNRho0     = 0;
  fNRhoPlus  = 0;
  fNRhoMinus = 0;
}
//___________________________________________________________________________
void XclsTag::SetResonance(Resonance_t res)
{
  fResonance = res;
}
//___________________________________________________________________________
void XclsTag::SetDecayMode(int decay_mode)
{
  fDecayMode = decay_mode;
}
//___________________________________________________________________________
void XclsTag::SetFinalQuark(int finalquark_pdgc)
{
  fIsFinalQuarkEvent     = true;
  fFinalQuarkPdg = finalquark_pdgc; // leave as 0 (default) for inclusive charm
}
//___________________________________________________________________________
void XclsTag::SetFinalLepton(int finallepton_pdgc)
{
  fIsFinalLeptonEvent     = true;
  fFinalLeptonPdg = finallepton_pdgc; // leave as 0 (default) for inclusive charm
}
//___________________________________________________________________________
void XclsTag::Reset(void)
{
  fIsStrangeEvent   = false ;
  fIsCharmEvent     = false ;
  fStrangeHadronPdg = 0 ;
  fCharmedHadronPdg = 0 ;
  fNProtons         = 0 ;
  fNNeutrons        = 0 ;
  fNPi0             = 0 ;
  fNPiPlus          = 0 ;
  fNPiMinus         = 0 ;
  fNSingleGammas    = 0 ;
  fNRho0            = 0 ;
  fNRhoPlus         = 0 ;
  fNRhoMinus        = 0 ;
  fResonance        = kNoResonance ;
  fDecayMode        = -1 ;
  fIsFinalQuarkEvent  = false;
  fFinalQuarkPdg      = 0;
  fIsFinalLeptonEvent = false;
  fFinalLeptonPdg     = 0;
}
//___________________________________________________________________________
void XclsTag::Copy(const XclsTag & xcls)
{
  fIsCharmEvent     = xcls.fIsCharmEvent;
  fCharmedHadronPdg = xcls.fCharmedHadronPdg;
  fIsStrangeEvent   = xcls.fIsStrangeEvent;
  fStrangeHadronPdg = xcls.fStrangeHadronPdg;
  fNProtons         = xcls.fNProtons;
  fNNeutrons        = xcls.fNNeutrons;
  fNPi0             = xcls.fNPi0;
  fNPiPlus          = xcls.fNPiPlus;
  fNPiMinus         = xcls.fNPiMinus;
  fNSingleGammas    = xcls.fNSingleGammas ;
  fNRho0             = xcls.fNRho0;
  fNRhoPlus          = xcls.fNRhoPlus;
  fNRhoMinus         = xcls.fNRhoMinus;
  fResonance        = xcls.fResonance;
  fDecayMode        = xcls.fDecayMode;
  fIsFinalQuarkEvent  = xcls.fIsFinalQuarkEvent;
  fFinalQuarkPdg      = xcls.fFinalQuarkPdg;
  fIsFinalLeptonEvent = xcls.fIsFinalLeptonEvent;
  fFinalLeptonPdg     = xcls.fFinalLeptonPdg;
}
//___________________________________________________________________________
/*
bool XclsTag::Compare(const XclsTag & xcls) const
{
  return (
     fIsCharmEvent     == xcls.fIsCharmEvent      &&
     fCharmedHadronPdg == xcls.fCharmedHadronPdg  &&
     fNProtons         == xcls.fNProtons          &&
     fNNeutrons        == xcls.fNNeutrons         &&
     fNPi0             == xcls.fNPi0              &&
     fNPiPlus          == xcls.fNPiPlus           &&
     fNPiMinus         == xcls.fNPiMinus          &&
     fResonance        == xcls.fResonance         &&
     fResonance        == xcls.fResonance         &&
  );
}*/
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

  if(fIsStrangeEvent) {
    tag << "strange:";
    if(fStrangeHadronPdg) tag << fStrangeHadronPdg;
    else tag << "incl";
    need_separator = true;
  }

  bool multset =
       fNProtons>0 || fNNeutrons>0 ||
       fNPiPlus>0 || fNPiMinus>0 || fNPi0>0 ||
       fNSingleGammas>0 ||
       fNRho0>0 || fNRhoPlus>0 || fNRhoMinus>0 ;
  if(multset) {
    if(need_separator) tag << ";";
    tag << "hmult:"
        << "(p=" << fNProtons << ",n=" << fNNeutrons
        << ",pi+=" << fNPiPlus << ",pi-=" << fNPiMinus << ",pi0=" << fNPi0
        << ",gamma=" << fNSingleGammas
        << ",rho+=" << fNRhoPlus << ",rho-=" << fNRhoMinus << ",rho0=" << fNRho0
        << ")";
  }

  if(this->KnownResonance()) {
    if(need_separator) tag << ";";
    tag << "res:" << fResonance;
  }

  if(fDecayMode != -1) {
    tag << "dec:" << fDecayMode;
  }

  if(fIsFinalQuarkEvent) {
    tag << "finalquark:" << fFinalQuarkPdg;
  }

  if(fIsFinalLeptonEvent) {
    tag << "finallepton:" << fFinalLeptonPdg;
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

  stream << " |--> strange prod.  : "
         << utils::print::BoolAsString(fIsStrangeEvent);
  if(fIsStrangeEvent) {
     if(!fStrangeHadronPdg) stream << " [inclusive]";
     else  {
       stream << " - Strange hadron PDG-code = " << fStrangeHadronPdg;

       TParticlePDG * chadr = PDGLibrary::Instance()->Find( fStrangeHadronPdg );
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

  stream << " |--> f/s Other    :"
         << " N(gamma) = "    << fNSingleGammas
         << " N(Rho^0) = "    << fNRho0
         << " N(Rho^+) = "    << fNRhoPlus
         << " N(Rho^-) = "    << fNRhoMinus
         << endl;

  stream << " |--> resonance    : ";
  if(this->KnownResonance()) {
     stream << res::AsString(fResonance);
  } else {
     stream << "[not set]";
  }

  stream << endl;

  stream << " |--> final quark prod.  : "
         << utils::print::BoolAsString(fIsFinalQuarkEvent);
  if(fIsFinalQuarkEvent) {
    stream << " - Final Quark PDG-code = " << fFinalQuarkPdg;
    TParticlePDG * chadr = PDGLibrary::Instance()->Find( fFinalQuarkPdg );
    if(chadr) stream << " (" << chadr->GetName() << ")";
  }

  stream << endl;

  stream << " |--> final lepton prod.  : "
         << utils::print::BoolAsString(fIsFinalLeptonEvent);
  if(fIsFinalLeptonEvent) {
    stream << " - Final Lepton PDG-code = " << fFinalLeptonPdg;
    TParticlePDG * chadr = PDGLibrary::Instance()->Find( fFinalLeptonPdg );
    if(chadr) stream << " (" << chadr->GetName() << ")";
  }

  stream << endl;
}
//___________________________________________________________________________
//bool XclsTag::operator == (const XclsTag & xcls) const
//{
//  return this->Compare(xcls);
//}
//___________________________________________________________________________
XclsTag & XclsTag::operator = (const XclsTag & xcls)
{
  this->Copy(xcls);
  return (*this);
}
//___________________________________________________________________________
