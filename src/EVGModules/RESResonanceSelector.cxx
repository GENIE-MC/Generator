//____________________________________________________________________________
/*!

\class   genie::RESResonanceSelector

\brief   Generates a baryon resonance for (v+N->Resonance->pi+X) events &
         adds it to the event record.

         Is a concrete implementation of the VtxGeneratorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 18, 2004

*/
//____________________________________________________________________________

#include <vector>
#include <sstream>

#include "BaryonResonance/BaryonResUtils.h"
#include "Base/XSecAlgorithmI.h"
#include "EVGModules/RESResonanceSelector.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepOrder.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"

using std::vector;
using std::ostringstream;

using namespace genie;

//___________________________________________________________________________
RESResonanceSelector::RESResonanceSelector() :
EventRecordVisitorI()
{
  fName = "genie::RESResonanceSelector";
}
//___________________________________________________________________________
RESResonanceSelector::RESResonanceSelector(const char * param_set) :
EventRecordVisitorI(param_set)
{
  fName = "genie::RESResonanceSelector";

  this->FindConfig();
}
//___________________________________________________________________________
RESResonanceSelector::~RESResonanceSelector()
{

}
//___________________________________________________________________________
void RESResonanceSelector::ProcessEventRecord(GHepRecord * evrec) const
{
  //-- select a baryon resonance
  Resonance_t res = this->SelectResonance(evrec);

  //-- add the resonance at the event summary
  Interaction * interaction = evrec->GetInteraction();
  interaction->GetScatParamsPtr()->Set("resonance-id", res );

  //-- add an entry at the GHep event record & the event summary
  this->AddResonance(evrec);
}
//___________________________________________________________________________
Resonance_t RESResonanceSelector::SelectResonance(GHepRecord * evrec) const
{
// Select a baryon resonance from a list of considered resonances based on
// their differential d^2xsec/dWdQ^2 cross sections

  LOG("RESSelector", pINFO) << "Selecting a baryon resonance";

  //-- Get a list of all baryon resonances considered
  BaryonResList res_list;
  this->FillResonanceList(res_list);

  //-- Get a d^2xsec/dWdQ^2|RES cross section calculator (with B/W wgt ON)

  LOG("RESSelector", pDEBUG) << "Getting the differential xsec algorithm";
  const XSecAlgorithmI * xsec_alg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                 "single-res-xsec-alg-name", "single-res-xsec-param-set"));

  //-- Figure out what the resonance charge should be.
  Interaction * interaction = evrec->GetInteraction();
  int q_res = this->ResQ(interaction);

  //-- Loop over all considered baryon resonances and compute the double
  //   differential cross section for the selected kinematical variables

  double xsec_sum  = 0;
  unsigned int nres = res_list.NResonances();
  vector<double> xsec_vec(nres);

  for(unsigned int ires = 0; ires < nres; ires++) {

     //-- Current resonance
     Resonance_t res = res_list.ResonanceId(ires);

     //-- Pass the current resonance as a 'scattering parameter'
     //   compute the differential cross section d^2xsec/dWdQ^2
     //   (do it only for resonances that can conserve charge)
     interaction->GetScatParamsPtr()->Set("resonance-id", (int) res );
     double xsec = 0;
     bool   skip = (q_res==2 && !res_utils::IsDelta(res));

     if(!skip) xsec = xsec_alg->XSec(interaction);
     else {
       SLOG("RESSelector", pINFO)
                 << "RES: " << res_utils::AsString(res)
                         << " would not conserve charge -- skipping it";
     }
     //-- For the ith resonance store the sum of (xsec) * (breit-wigner)
     //   for the resonances in the range [0,i]
     xsec_sum      += xsec;
     xsec_vec[ires] = xsec_sum;

     SLOG("RESSelector", pINFO)
          << "Resonances (0->" << ires << "): "
           << "Sum{ BW(W) * d^2xsec(E,W,Q^2)/dWd*Q^2 } = "  << xsec_sum;
  } // res

  //-- Use the computed differential cross sections to select a resonance

  RandomGen * rnd = RandomGen::Instance();
  double R = xsec_sum * rnd->Random2().Rndm();

  SLOG("RESSelector", pINFO) << "R = " << R;

  for(unsigned int ires = 0; ires < nres; ires++) {
     SLOG("RESSelector", pDEBUG)
                    << "SUM-XSEC(0->" << ires <<") = " << xsec_vec[ires];

     if( R < xsec_vec[ires] ) {
        Resonance_t sres = res_list.ResonanceId(ires); // selected RES.
        LOG("RESSelector", pINFO)
                   << "Selected RES = " << res_utils::AsString(sres);
        return sres;
     }
  }
  LOG("RESSelector", pERROR) << "** Failed to select a resonance";
  return kNoResonance;
}
//___________________________________________________________________________
void RESResonanceSelector::FillResonanceList(BaryonResList& res_list) const
{
// Creates a list with all the resonances that the user wants me to consider.
// The list of resonances is found in this algorithm's config file.

  LOG("RESSelector", pDEBUG) << "Getting the baryon resonance list";

  assert( fConfig->Exists("resonance-name-list") );
  string resonances = fConfig->GetString("resonance-name-list");

  SLOG("RESSelector", pDEBUG) << "Resonance list: " << resonances;

  res_list.DecodeFromNameList(resonances);

  LOG("RESSelector", pINFO) << res_list;
}
//___________________________________________________________________________
void RESResonanceSelector::AddResonance(GHepRecord * evrec) const
{
  //-- Get the interaction & initial state objects

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  Resonance_t res = res_utils::FromInteraction(interaction);

  //-- Get all initial & final state particles 4-momenta (in the LAB frame)

  //incoming v:
  TLorentzVector * nu_p4 = init_state.GetProbeP4(kRfLab);
  assert(nu_p4);

  //struck nucleon:
  TLorentzVector * nucl_p4 = init_state.GetTarget().StruckNucleonP4();
  assert(nucl_p4);

  //final state primary lepton:
  int fsl_pdgc = interaction->GetFSPrimaryLepton()->PdgCode();
  GHepParticle * fsl = evrec->FindParticle(
                                      fsl_pdgc, kIStStableFinalState, 0);
  assert(fsl);

  //-- Compute the resonance 4-momentum
  // Pv(Ev,pxv,pyv,pzv) + Pnucl(En,pxn,pyn,pzn) = Pl(El,pxl,pyl,pzl) + P

  double E  = nu_p4->Energy() + nucl_p4->Energy() - fsl->E();
  double px = nu_p4->Px()     + nucl_p4->Px()     - fsl->Px();
  double py = nu_p4->Py()     + nucl_p4->Py()     - fsl->Py();
  double pz = nu_p4->Pz()     + nucl_p4->Pz()     - fsl->Pz();

  delete nu_p4;

  //-- Determine the RES pdg code (from the selected Resonance_t & charge)
  int q_res    = this->ResQ(interaction);
  int res_pdgc = res_utils::PdgCode(res,q_res);
  LOG("RESSelector", pINFO)
           << "Adding RES with PDGC = " << res_pdgc << ", Q = " << q_res;

  //-- Add the resonance at the EventRecord

  GHepStatus_t ist = kIstPreDecayResonantState;
  int mom = GHepOrder::StruckNucleonPosition(interaction);

  evrec->AddParticle(res_pdgc, ist, mom,-1,-1,-1, px,py,pz,E, 0,0,0,0);
}
//___________________________________________________________________________
int RESResonanceSelector::ResQ(const Interaction * interaction) const
{
// Figure out what th e resonance charge should be to conserve charge

  const InitialState & init_state = interaction->GetInitialState();

  int nuc_pdgc = init_state.GetTarget().StruckNucleonPDGCode();
  int fsl_pdgc = interaction->GetFSPrimaryLepton()->PdgCode();

  int q_nuc    = int( PDGLibrary::Instance()->Find(nuc_pdgc)->Charge() );
  int q_fsl    = int( PDGLibrary::Instance()->Find(fsl_pdgc)->Charge() );
  int q_res    = (q_nuc - q_fsl) /3;

  return q_res;
}
//___________________________________________________________________________
