//____________________________________________________________________________
/*!

\class    genie::RSExclusiveRESPXSec

\brief    Computes the differential resonance cross section  for an exclusive
          resonance reaction.

          The computed xsec is the double differential d^2 xsec/ dQ^2 dW \n
          where \n
          \li Q^2 : momentum transfer ^ 2
          \li W   : invariant mass of the final state hadronic system

          The cross section is computed for an input list of resonances
          as the sum of the Rein-Seghal single resonance cross sections
          weighted:

          \li With the value of their Breit-Wigner distributions at the given
              W,Q^2 (The code for BW weighting is included in the single
              resonance cross section algorithm. The user needs to make sure
              that he does not run the single resonance cross section code with
              a configuration that inhibits weighting).

          \li With the isospin Glebsch-Gordon coefficient determining the
              contribution of each resonance to the exclusive final state.

          \li With the BR for the produced resonance to decay into the given
              exclusive final state.

          In this algorithm we follow the non-coherent approach: we sum
          the weighted resonance production cross sections rather than the
          resonance production amplitudes.

          Is a concrete implementation of the XSecAlgorithmI initerface.

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#include "Algorithm/AlgFactory.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Interaction/SppChannel.h"
#include "Messenger/Messenger.h"
#include "ReinSeghal/RSExclusiveRESPXSec.h"

using namespace genie;

//____________________________________________________________________________
RSExclusiveRESPXSec::RSExclusiveRESPXSec() :
XSecAlgorithmI("genie::RSExclusiveRESPXSec")
{

}
//____________________________________________________________________________
RSExclusiveRESPXSec::RSExclusiveRESPXSec(string config) :
XSecAlgorithmI("genie::RSExclusiveRESPXSec", config)
{

}
//____________________________________________________________________________
RSExclusiveRESPXSec::~RSExclusiveRESPXSec()
{

}
//____________________________________________________________________________
double RSExclusiveRESPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghalRes", pINFO) << *fConfig;

  //-- Get the requested SPP channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  if( spp_channel == kSppNull ) {
    LOG("ReinSeghalRes", pERROR)
            << "\n *** Insufficient SPP exclusive final state information!";
    return 0;
  } else {
     LOG("ReinSeghalRes", pINFO)
                       << "Reaction: " << SppChannel::AsString(spp_channel);
  }

  //-- Loop over the specified list of baryon resonances and compute
  //   the total weighted cross section

  unsigned int nres = fResList.NResonances();
  LOG("ReinSeghalRes", pINFO)
      << "Computing a weighted cross section for " << *interaction
      << "using " << nres << " resonances";

  double xsec = 0;

  for(unsigned int ires = 0; ires < nres; ires++) {

     //-- Get next resonance from the resonance list
     Resonance_t res = fResList.ResonanceId(ires);

     //-- Find out the charge that the produced resonance should have to
     //   yield the requested exclusive final state.
     //   Eg if the resonance is kP33_1232, find out which of the Delta-,
     //      Delta0, Delta+ or Delta++ is relevant here.
     //int QRes = SppChannel::ResonanceCharge(spp_channel);

     //-- Set current resonance to interaction object
     interaction->GetExclusiveTagPtr()->SetResonance(res);

     //-- Get the Breit-Wigner weighted xsec for the current resonance
     double rxsec = fSingleResXSecModel->XSec(interaction);

     //-- Get the BR for the (resonance) -> (exclusive final state)
     double br = SppChannel::BranchingRatio(spp_channel, res);

     //-- Get the Isospin Glebsch-Gordon coefficient for the given resonance
     //   and exclusive final state
     double igg = SppChannel::IsospinWeight(spp_channel, res);

     //-- Compute the weighted xsec
     //  (total weight = Breit-Wigner * BR * isospin Glebsch-Gordon)
     double res_xsec_contrib = rxsec*br*igg;

     SLOG("ReinSeghalRes", pINFO)
         << "Contrib. from [" << utils::res::AsString(res) << "] = "
         << "<Glebsch-Gordon = " << igg
         << "> * <BR(->1pi) = " << br
         << "> * <Breit-Wigner * d^2xsec/dQ^2dW = " << rxsec
         << "> = " << res_xsec_contrib;

     //-- Add contribution of this resonance to the cross section
     xsec += res_xsec_contrib;
  }

  LOG("ReinSeghalRes", pINFO) << "d^2 xsec/ dQ^2 dW = " << xsec;

  return xsec;
}
//____________________________________________________________________________
void RSExclusiveRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
  this->GetResonanceList();
}
//____________________________________________________________________________
void RSExclusiveRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
  this->GetResonanceList();
}
//____________________________________________________________________________
void RSExclusiveRESPXSec::LoadSubAlg(void)
{
// load the single resonance cross section algorithm specified in the config.

  fSingleResXSecModel = 0;

  fSingleResXSecModel =
       dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                    "single-res-xsec-alg-name", "single-res-xsec-param-set"));
  assert(fSingleResXSecModel);
}
//____________________________________________________________________________
void RSExclusiveRESPXSec::GetResonanceList(void)
{
// create the baryon resonance list specified in the config.

  fResList.Clear();

  assert( fConfig->Exists("resonance-name-list") );
  string resonanes = fConfig->GetString("resonance-name-list");

  //-- Create a BaryonResList by decoding the resonance list from
  //   the XML input
  //   The list of resonances can be specified as a string with
  //   comma separated resonance names (eg "P33(1233),S11(1535),D13(1520)")
  //   The BaryonResList can also decode lists of pdg-codes or
  //   resonance-ids (Resonance_t enumerations).
  //   Support for this will be added here as well.

  fResList.DecodeFromNameList(resonanes);
}
//____________________________________________________________________________

