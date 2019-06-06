//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new RES package from its previous location (EVGModules).
 @ Jul 23, 2010 - CA
   Use ResonanceCharge() from base class. Function removed from utils::res.

*/
//____________________________________________________________________________

#include <vector>
#include <sstream>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen//RunningThreadInfo.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Physics/Resonance/EventGen/RSPPResonanceSelector.h"

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
  Interaction * interaction = evrec->Summary();
  interaction->ExclTagPtr()->SetResonance(res);

  //-- add an entry at the GHep event record & the event summary
  this->AddResonance(evrec);
}
//___________________________________________________________________________
Resonance_t RSPPResonanceSelector::SelectResonance(GHepRecord * evrec) const
{
// Select a baryon resonance from a list of considered resonances based on
// their differential d^2xsec/dWdQ^2 cross sections

  LOG("RESSelector", pNOTICE) << "Selecting a baryon resonance";

  Interaction * interaction = evrec->Summary();

  //-- Figure out what the resonance charge should be.
  int q_res = this->ResonanceCharge(evrec);

  //-- Use selected kinematics
  interaction->KinePtr()->UseSelectedKinematics();

  //-- Trust kinematics and process type already set.
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Access cross section algorithm for running thread
  //   The algorithm must be able to compute the RES contribution to
  //   the specified RES/SPP channel
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  const XSecAlgorithmI * xsecalg = evg->CrossSectionAlg();

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
     interaction->ExclTagPtr()->SetResonance(res);

     double xsec = 0;
     bool   skip = (q_res==2 && !utils::res::IsDelta(res));

     if(!skip) xsec = xsecalg->XSec(interaction,kPSWQ2fE);
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
  interaction->KinePtr()->ClearRunningValues();

  //-- Use the computed differential cross sections to select a resonance
  RandomGen * rnd = RandomGen::Instance();
  double R = xsec_sum * rnd->RndGen().Rndm();

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
  Interaction * interaction = evrec->Summary();
  Resonance_t res = interaction->ExclTag().Resonance();
  int charge = this->ResonanceCharge(evrec);
  int pdgc   = utils::res::PdgCode(res,charge);

  LOG("RESSelector", pNOTICE)
               << "Adding RES with PDGC = " << pdgc << ", Q = " << charge;

  //-- Add the resonance at the EventRecord
  GHepStatus_t ist = kIStPreDecayResonantState;
  int mom = evrec->HitNucleonPosition();

  evrec->AddParticle(
        pdgc, ist, mom,-1,-1,-1, p4.Px(),p4.Py(),p4.Pz(),p4.E(), 0,0,0,0);
}
//___________________________________________________________________________
void RSPPResonanceSelector::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSPPResonanceSelector::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSPPResonanceSelector::LoadConfig(void)
{
  fResList.Clear();
  string resonances = "";
  this->GetParam("ResonanceNameList", resonances);
  SLOG("RESSelector", pDEBUG) << "Resonance list: " << resonances;

  fResList.DecodeFromNameList(resonances);
  LOG("RESSelector", pINFO) << fResList;
}
//____________________________________________________________________________
