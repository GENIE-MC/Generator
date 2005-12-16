//____________________________________________________________________________
/*!

\class    genie::RSSingleRESPXSec

\brief    Computes the (+/- breit-wigner weighted) differential RES xsec
          d^2 xsec/ dQ^2 dW for a *single* Baryon Resonance.

          The actual Rein-Seghal model is coded at ReinSeghalPartialXSec.

          This algorithm merely decides which of the ReinSeghalPartialXSec
          algorithm instances to run depending on the initial state:
          1) CC (n+p), 2) NC (n), 3) NC (p).

          Concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "ReinSeghal/RSSingleRESPXSec.h"

using namespace genie;

//____________________________________________________________________________
RSSingleRESPXSec::RSSingleRESPXSec() :
XSecAlgorithmI("genie::RSSingleRESPXSec")
{

}
//____________________________________________________________________________
RSSingleRESPXSec::RSSingleRESPXSec(string config) :
XSecAlgorithmI("genie::RSSingleRESPXSec", config)
{

}
//____________________________________________________________________________
RSSingleRESPXSec::~RSSingleRESPXSec()
{

}
//____________________________________________________________________________
double RSSingleRESPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghalRes", pDEBUG) << *fConfig;

  //-- Find out the weak current & the struck nucleon pdg code
  const InitialState & init_state = interaction->GetInitialState();
  int nucleon_pdgc = init_state.GetTarget().StruckNucleonPDGCode();

  bool is_CC = interaction->GetProcessInfo().IsWeakCC();
  bool is_NC = interaction->GetProcessInfo().IsWeakNC();
  bool is_p  = pdg::IsProton(nucleon_pdgc);
  bool is_n  = pdg::IsNeutron(nucleon_pdgc);

  //-- Complain if you can not handle it & return 0
  if(!is_CC && !is_NC) {
     LOG("ReinSeghalRes", pERROR) << "*** Not a CC or a NC interaction";
     return 0;
  }
  if(!is_p && !is_n ) {
     LOG("ReinSeghalRes", pERROR) << "*** Undefined struck nucleon";
     return 0;
  }

  //-- Figure out which algorithm to use
  const XSecAlgorithmI * xsec_alg = 0;
  if(is_CC) xsec_alg = fRSCC;
  else {
      if (is_p) xsec_alg = fRSNCp;
      else      xsec_alg = fRSNCn;
  }
  assert(xsec_alg);

  //-- Ask the selected algorithm to compute d^2 xsec/ dQ^2 dW
  double xsec = xsec_alg->XSec(interaction);
  return xsec;
}
//____________________________________________________________________________
void RSSingleRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
}
//____________________________________________________________________________
void RSSingleRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
}
//____________________________________________________________________________
void RSSingleRESPXSec::LoadSubAlg(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed
  fRSCC  = 0;
  fRSNCp = 0;
  fRSNCn = 0;

  fRSCC  = dynamic_cast<const XSecAlgorithmI *>
                             (this->SubAlg("CC-alg-name","CC-param-set"));
  fRSNCp = dynamic_cast<const XSecAlgorithmI *>
                             (this->SubAlg("CC-alg-name","NC-p-param-set"));
  fRSNCn = dynamic_cast<const XSecAlgorithmI *>
                             (this->SubAlg("CC-alg-name","NC-n-param-set"));
  assert(fRSCC);
  assert(fRSNCp);
  assert(fRSNCn);
}
//____________________________________________________________________________

