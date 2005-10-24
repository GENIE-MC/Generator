//____________________________________________________________________________
/*!

\class    genie::RSListRESPXSec

\brief    Computes the total differential RES cross section RES for a list
          of baryon resonances.

          The cross section is the double differential d^2 xsec/ dQ^2 dW \n
          where \n
          \li Q^2 : momentum transfer ^ 2
          \li W   : invariant mass of the final state hadronic system

          The cross section is computed for an input list of resonances
          as the sum of the Rein-Seghal single resonance cross sections
          weighted with the values of their Breit-Wigner distributions at
          the given W,Q^2 (The user needs to make sure that he does not
          run the single resonance cross section code with a configuration
          that inhibits weighting).

          Concrete implementation of the XSecAlgorithmI interface

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#include "Algorithm/AlgFactory.h"
#include "BaryonResonance/BaryonResList.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Messenger/Messenger.h"
#include "ReinSeghal/RSListRESPXSec.h"

using namespace genie;

//____________________________________________________________________________
RSListRESPXSec::RSListRESPXSec() :
XSecAlgorithmI("genie::RSListRESPXSec")
{

}
//____________________________________________________________________________
RSListRESPXSec::RSListRESPXSec(string config) :
XSecAlgorithmI("genie::RSListRESPXSec", config)
{

}
//____________________________________________________________________________
RSListRESPXSec::~RSListRESPXSec()
{

}
//____________________________________________________________________________
double RSListRESPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghalRes", pDEBUG) << *fConfig;

  //-- Get the requested single resonance cross section model

  const XSecAlgorithmI * res_xsec_model = this->SingleResXSecModel();

  //-- Create a BaryonResList by decoding the resonance list from
  //   the XML input
  //   The list of resonances can be specified as a string with
  //   comma separated resonance names (eg "P33(1233),S11(1535),D13(1520)")
  //   The BaryonResList can also decode lists of pdg-codes or
  //   resonance-ids (Resonance_t enumerations).
  //   Support for this will be added here as well.
  
  assert( fConfig->Exists("resonance-name-list") );
    
  string resonanes = fConfig->GetString("resonance-name-list");
  
  BaryonResList res_list;

  res_list.DecodeFromNameList(resonanes);
  
  //-- Loop over the specified list of baryon resonances and compute
  //   the total weighted cross section

  unsigned int nres = res_list.NResonances();
  
  LOG("ReinSeghalRes", pDEBUG)
      << "Computing a weighted cross section for = " << nres << " resonances";
    
  double xsec = 0;
  
  for(unsigned int ires = 0; ires < nres; ires++) {

     //-- Current resonance
     Resonance_t res = res_list.ResonanceId(ires);
  
     //-- Set current resonance to interaction object
     interaction->GetScatParamsPtr()->Set("resonance-id", (int) res);
     
     //-- Get the Breit-Wigner weighted xsec for the current resonance     
     double rxsec = res_xsec_model->XSec(interaction);

     xsec += rxsec;
  }    

  LOG("ReinSeghalRes", pDEBUG) << "Res List: d^2 xsec/ dQ^2 dW = " << xsec;
  
  return xsec;
}
//____________________________________________________________________________
const XSecAlgorithmI * RSListRESPXSec::SingleResXSecModel(void) const
{
// Retrieves the specifed cross section model for a single resonance

  assert(
       fConfig->Exists("single-res-xsec-alg-name")   &&
       fConfig->Exists("single-res-xsec-param-set")
  );

  string alg_name  = fConfig->GetString("single-res-xsec-alg-name");
  string param_set = fConfig->GetString("single-res-xsec-param-set");
  

  AlgFactory * algf = AlgFactory::Instance();
  
  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);

  const XSecAlgorithmI * xs = dynamic_cast<const XSecAlgorithmI *> (algbase);

  assert(xs);

  return xs;
}
//____________________________________________________________________________
