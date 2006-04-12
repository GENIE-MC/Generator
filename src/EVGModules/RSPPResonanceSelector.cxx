//____________________________________________________________________________
/*!

\class   genie::RSPPResonanceSelector

\brief   Generates an intermediate baryon resonance for exclusive interactions
         proceeding through resonance productions and adds it to the event
         record. The resonance is selected based on its contribution to the
         selected exclusive reaction cross section.
         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 18, 2004

*/
//____________________________________________________________________________

#include <vector>
#include <sstream>

#include "BaryonResonance/BaryonResUtils.h"
#include "Base/XSecAlgorithmI.h"
#include "EVGModules/RSPPResonanceSelector.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"

using std::vector;
using std::ostringstream;

using namespace genie;

//___________________________________________________________________________
RSPPResonanceSelector::RSPPResonanceSelector() :
HadronicSystemGenerator("genie::RSPPResonanceSelector")
{

}
//___________________________________________________________________________
RSPPResonanceSelector::RSPPResonanceSelector(string config) :
HadronicSystemGenerator("genie::RSPPResonanceSelector", config)
{

}
//___________________________________________________________________________
RSPPResonanceSelector::~RSPPResonanceSelector()
{

}
//___________________________________________________________________________
void RSPPResonanceSelector::ProcessEventRecord(GHepRecord * evrec) const
{
  //-- select a baryon resonance
  Resonance_t res = this->SelectResonance(evrec);
  assert(res != kNoResonance);

  //-- add the resonance at the event summary
  Interaction * interaction = evrec->GetInteraction();
  interaction->GetExclusiveTagPtr()->SetResonance(res);

  //-- add an entry at the GHep event record & the event summary
  this->AddResonance(evrec);
}
//___________________________________________________________________________
Resonance_t RSPPResonanceSelector::SelectResonance(GHepRecord * evrec) const
{
// Select a baryon resonance from a list of considered resonances based on
// their differential d^2xsec/dWdQ^2 cross sections

  LOG("RESSelector", pNOTICE) << "Selecting a baryon resonance";

  //-- Figure out what the resonance charge should be.
  Interaction * interaction = evrec->GetInteraction();
  int q_res = utils::res::ResonanceCharge(interaction);

  //-- Use selected kinematics
  interaction->GetKinematicsPtr()->UseSelectedKinematics();

  //-- Trust kinematics and process type already set.
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Loop over all considered baryon resonances and compute the double
  //   differential cross section for the selected kinematical variables

  double xsec_sum  = 0;
  unsigned int nres = fResList.NResonances();
  vector<double> xsec_vec(nres);

  for(unsigned int ires = 0; ires < nres; ires++) {

     //-- Current resonance
     Resonance_t res = fResList.ResonanceId(ires);

     //-- Set the current resonance at the interaction summary
     //   compute the differential cross section d^2xsec/dWdQ^2
     //   (do it only for resonances that can conserve charge)
     interaction->GetExclusiveTagPtr()->SetResonance(res);

     double xsec = 0;
     bool   skip = (q_res==2 && !utils::res::IsDelta(res));

     if(!skip) xsec = fXSecAlg->XSec(interaction);
     else {
       SLOG("RESSelector", pNOTICE)
                 << "RES: " << utils::res::AsString(res)
                         << " would not conserve charge -- skipping it";
     }
     //-- For the ith resonance store the sum of (xsec) * (breit-wigner)
     //   for the resonances in the range [0,i]
     xsec_sum      += xsec;
     xsec_vec[ires] = xsec_sum;

     SLOG("RESSelector", pNOTICE)
          << "Resonances (0->" << ires << "): "
           << "Sum{ BW(W) * d^2xsec(E,W,Q^2)/dWd*Q^2 } = "  << xsec_sum;
  } // res

  //-- Reset 'trust' bits
  interaction->ResetBit(kISkipProcessChk);
  interaction->ResetBit(kISkipKinematicChk);

  //-- Reset running kinematics
  interaction->GetKinematicsPtr()->ClearRunningValues();

  //-- Use the computed differential cross sections to select a resonance
  RandomGen * rnd = RandomGen::Instance();
  double R = xsec_sum * rnd->Random1().Rndm();

  SLOG("RESSelector", pDEBUG) << "R = " << R;

  for(unsigned int ires = 0; ires < nres; ires++) {
     SLOG("RESSelector", pDEBUG)
                    << "SUM-XSEC(0->" << ires <<") = " << xsec_vec[ires];

     if(R < xsec_vec[ires]) {
        Resonance_t sres = fResList.ResonanceId(ires); // selected RES.
        LOG("RESSelector", pNOTICE)
                   << "Selected RES = " << utils::res::AsString(sres);
        return sres;
     }
  }
  LOG("RESSelector", pERROR) << "** Failed to select a resonance";
  return kNoResonance;
}
//___________________________________________________________________________
void RSPPResonanceSelector::AddResonance(GHepRecord * evrec) const
{
  // compute RES p4 = p4(neutrino) + p4(hit nucleon) - p4(primary lepton)
  TLorentzVector p4 = this->Hadronic4pLAB(evrec);

  //-- Determine the RES pdg code (from the selected Resonance_t & charge)
  Interaction * interaction = evrec->GetInteraction();
  Resonance_t res = interaction->GetExclusiveTag().Resonance();
  int charge = utils::res::ResonanceCharge(interaction);
  int pdgc   = utils::res::PdgCode(res,charge);

  LOG("RESSelector", pNOTICE)
               << "Adding RES with PDGC = " << pdgc << ", Q = " << charge;

  //-- Add the resonance at the EventRecord
  GHepStatus_t ist = kIStPreDecayResonantState;
  int mom = evrec->StruckNucleonPosition();

  evrec->AddParticle(
        pdgc, ist, mom,-1,-1,-1, p4.Px(),p4.Py(),p4.Pz(),p4.E(), 0,0,0,0);
}
//___________________________________________________________________________
void RSPPResonanceSelector::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
  this->LoadConfigData();
}
//____________________________________________________________________________
void RSPPResonanceSelector::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadSubAlg();
  this->LoadConfigData();
}
//____________________________________________________________________________
void RSPPResonanceSelector::LoadConfigData(void)
{
  // Create the list with all the baryon resonances that the user wants me to
  // consider (from this algorithm's config file).

  LOG("RESSelector", pDEBUG) << "Getting the baryon resonance list";

  fResList.Clear();
  assert( fConfig->Exists("resonance-name-list") );
  string resonances = fConfig->GetString("resonance-name-list");
  SLOG("RESSelector", pDEBUG) << "Resonance list: " << resonances;

  fResList.DecodeFromNameList(resonances);
  LOG("RESSelector", pINFO) << fResList;
}
//____________________________________________________________________________
void RSPPResonanceSelector::LoadSubAlg(void)
{
  fXSecAlg = 0;

  //-- Get a d^2xsec/dWdQ^2|SPP cross section calculator 
  //   The algorithm should be able to compute xsec for SPP channels and also
  //   allow one to compute the contribution of specific resonances
  LOG("RESSelector", pDEBUG) << "Getting the differential xsec algorithm";
  fXSecAlg = dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                              "spp-xsec-alg-name", "spp-xsec-param-set"));

  assert(fXSecAlg);
}
//____________________________________________________________________________

