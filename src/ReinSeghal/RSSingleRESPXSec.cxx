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

#include "AlgFactory/AlgFactory.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "ReinSeghal/RSSingleRESPXSec.h"

using namespace genie;

//____________________________________________________________________________
RSSingleRESPXSec::RSSingleRESPXSec() :
XSecAlgorithmI()
{
  fName = "genie::RSSingleRESPXSec";
}
//____________________________________________________________________________
RSSingleRESPXSec::RSSingleRESPXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::RSSingleRESPXSec";

  FindConfig();
}
//____________________________________________________________________________
RSSingleRESPXSec::~RSSingleRESPXSec()
{

}
//____________________________________________________________________________
double RSSingleRESPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghal", pDEBUG) << *fConfig;

  //-- Find out the weak current & the struck nucleon pdg code
  
  const InitialState &     init_state = interaction -> GetInitialState();

  bool is_CC = interaction->GetProcessInfo().IsWeakCC();
  bool is_NC = interaction->GetProcessInfo().IsWeakNC();

  int nucleon_pdgc = init_state.GetTarget().StruckNucleonPDGCode();

  //-- Complain if you can not handle it & return 0

  if( !is_CC && !is_NC) {
    
     LOG("ReinSeghal", pERROR) << "*** Not a CC or a NC interaction";
     return 0;
  }
    
  if(! pdg::IsNeutronOrProton(nucleon_pdgc) ) {

     LOG("ReinSeghal", pERROR) << "*** Undefined struck nucleon";
     return 0;
  }

  //-- Create the correct keys for requesting the right Rein-Seghal
  //   partial cross section algorithm for the interaction at hand

  string alg_key, config_key;

  if(is_CC)
  {
      alg_key    = "CC-alg-name";
      config_key = "CC-param-set";
  }
  else {
      if ( pdg::IsProton(nucleon_pdgc) )
      {
         alg_key    = "CC-alg-name";
         config_key = "CC-param-set";      
      }
      else
      {
         alg_key    = "CC-alg-name";
         config_key = "CC-param-set";      
      }      
  }  

  LOG("ReinSeghal", pDEBUG)
          << "\n alg-key = " << alg_key << " / config-key = " << config_key;

  //-- Look-up the configuration registry to get the values that
  //   correspond to these keys & retrieve the requested algorithm
  //   from the AlgFactory

  assert( fConfig->Exists(alg_key) && fConfig->Exists(config_key) );

  string alg_name  = fConfig->GetString( alg_key    );
  string param_set = fConfig->GetString( config_key );
  
  LOG("ReinSeghal", pDEBUG)
          << "\n alg-name = " << alg_name << " / param-set = " << param_set;

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);

  const XSecAlgorithmI * xsec_alg =
                             dynamic_cast<const XSecAlgorithmI *> (algbase);

  assert(xsec_alg);

  
  //-- Ask the selected algorithm to compute d^2 xsec/ dQ^2 dW
  
  double xsec = xsec_alg->XSec(interaction);

  Resonance_t res = res_utils::FromInteraction(interaction);
  
  LOG("ReinSeghal", pDEBUG)
         << "d^2xsec/dQ^2 dW [" << res_utils::AsString(res) << "] = " << xsec;
         
  return xsec;
}
//____________________________________________________________________________
