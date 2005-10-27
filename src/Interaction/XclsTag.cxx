//____________________________________________________________________________
/*!

\class    genie::XclsTag

\brief    Contains minimal information for tagging exclusive processes.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  December 08, 2004

*/
//___________________________________________________________________________

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
XclsTag::XclsTag()
{
  this->Initialize();
}
//___________________________________________________________________________
XclsTag::XclsTag(const XclsTag & xcls)
{
  this->Initialize();

  this->Copy(xcls);
}
//___________________________________________________________________________
XclsTag::~XclsTag()
{

}
//___________________________________________________________________________
bool XclsTag::IsInclusiveCharm(void) const
{
  return ( this->IsCharmEvent() && (this->CharmHadronPDGCode() == 0) );
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
void XclsTag::Initialize(void)
{
  fIsCharmEvent     = false;
  fCharmedHadronPdg = 0;

  fNProtons     = 0;
  fNNeutrons    = 0;
  fNPi0         = 0;
  fNPiPlus      = 0;
  fNPiMinus     = 0;

  fResonance    = kNoResonance;
}
//___________________________________________________________________________
void XclsTag::Copy(const XclsTag & xcls)
{
  fIsCharmEvent     = xcls.fIsCharmEvent;
  fCharmedHadronPdg = xcls.fCharmedHadronPdg;

  fNProtons     = xcls.fNProtons;
  fNNeutrons    = xcls.fNNeutrons;
  fNPi0         = xcls.fNPi0;
  fNPiPlus      = xcls.fNPiPlus;
  fNPiMinus     = xcls.fNPiMinus;

  fResonance    = xcls.fResonance;
}
//___________________________________________________________________________
string XclsTag::AsString(void) const
{
// codifies XclsTag state into a compact string:
// c=is-charm,charm-pdgc;nucl(p,n)=np,nn;pi(+,-,0)=npi+,npi-,npi0;res=respdg

  ostringstream tag;

  tag << "c=" << fIsCharmEvent << "," << fCharmedHadronPdg << ";";
  tag << "nucl(p,n)=" << fNProtons << "," << fNNeutrons << ";";
  tag << "pi(+,-,0)=" << fNPiPlus << "," << fNPiMinus << "," << fNPi0 << ";";
  tag << "res=" << fResonance;

  return tag.str();
}
//___________________________________________________________________________
void XclsTag::Print(ostream & stream) const
{
  stream << "[-] [Exclusive Process Info] " << endl;

  stream << " |--> charm        : "
         << utils::print::BoolAsString(fIsCharmEvent)
         << " - Charm hadron PDG-code = ";

  if(!fCharmedHadronPdg) stream << "[inclusive]";
  else  {
     stream << fCharmedHadronPdg;

     TParticlePDG * chadr = PDGLibrary::Instance()->Find( fCharmedHadronPdg );
     if(chadr)
        stream << " (" << chadr->GetName() << ")";
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

  if(this->KnownResonance()) {
     stream << " |--> resonance    : " << res::AsString(fResonance) << endl;
  }
}
//___________________________________________________________________________


