//____________________________________________________________________________
/*
 See the corresponding header file for the class documentation.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - March 09, 2006

 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "ReinSeghal/ReinSeghalRESXSec.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Cache.h"
#include "Utils/CacheBranchFx.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::units;

//____________________________________________________________________________
ReinSeghalRESXSec::ReinSeghalRESXSec() :
ReinSeghalRESXSecWithCache("genie::ReinSeghalRESXSec")
{

}
//____________________________________________________________________________
ReinSeghalRESXSec::ReinSeghalRESXSec(string config) :
ReinSeghalRESXSecWithCache("genie::ReinSeghalRESXSec", config)
{

}
//____________________________________________________________________________
ReinSeghalRESXSec::~ReinSeghalRESXSec()
{

}
//____________________________________________________________________________
double ReinSeghalRESXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  assert(kps==kPSfE);

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //-- Get cache
  Cache * cache = Cache::Instance();

  //-- Get init state and process information
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  InteractionType_t it = proc_info.InteractionTypeId();
  int nucleon_pdgc = target.HitNucPdg();
  int nu_pdgc      = init_state.ProbePdg();

  //-- Get neutrino energy in the struck nucleon rest frame
  double Ev = init_state.ProbeE(kRfHitNucRest);

  //-- Get the requested resonance
  Resonance_t res = interaction->ExclTag().Resonance();

  //-- Build a unique name for the cache branch
  string key = this->CacheBranchName(res, it, nu_pdgc, nucleon_pdgc);

  LOG("ReinSeghalResT", pINFO) << "Finding cache branch with key: " << key;

  CacheBranchFx * cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));

  if(!cache_branch) {
     LOG("ReinSeghalResT", pWARN)  
         << "No cached RES v-production data for input neutrino"
         << " (pdgc: " << nu_pdgc << ")";
     LOG("ReinSeghalResT", pWARN)  
         << "Wait while computing/caching RES production xsec first...";

     this->CacheResExcitationXSec(interaction); 

     LOG("ReinSeghalResT", pINFO) << "Done caching resonance xsec data";
     LOG("ReinSeghalResT", pINFO) 
               << "Finding newly created cache branch with key: " << key;
     cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
     assert(cache_branch);
  }
  const CacheBranchFx & cbranch = (*cache_branch);
  
  //-- Get cached resonance neutrinoproduction xsec
  //   (If E>Emax, assume xsec = xsec(Emax) - but do not evaluate the
  //    cross section spline at the end of its energy range-)
  double rxsec = (Ev<fEMax-1) ? cbranch(Ev) : cbranch(fEMax-1);

  SLOG("ReinSeghalResT", pNOTICE)  
    << "XSec[RES/" << utils::res::AsString(res)<< "/free] (Ev = " 
               << Ev << " GeV) = " << rxsec/(1E-38 *cm2)<< " x 1E-38 cm^2";

  //-- If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return rxsec;

  //-- number of scattering centers in the target
  int NNucl = (pdg::IsProton(nucleon_pdgc)) ? target.Z() : target.N();

  rxsec*=NNucl; // nuclear xsec 

  return rxsec;
}
//____________________________________________________________________________
bool ReinSeghalRESXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return fSingleResXSecModel->ValidProcess(interaction);
}
//____________________________________________________________________________
bool ReinSeghalRESXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> InitState();
  double Ev  = init_state.ProbeE(kRfHitNucRest);

  double EvThr = interaction->EnergyThreshold();
  if(Ev <= EvThr) return false;

  return true;
}
//____________________________________________________________________________
void ReinSeghalRESXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSeghalRESXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSeghalRESXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fSingleResXSecModel = 0;
  fIntegrator = 0;

  //-- get the requested d^2xsec/dxdy xsec algorithm to use
  fSingleResXSecModel =
       dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                    "single-res-xsec-alg-name", "single-res-xsec-param-set"));

  fIntegrator = dynamic_cast<const IntegratorI *> (
                 this->SubAlg("integrator-alg-name", "integrator-param-set"));

  assert (fSingleResXSecModel);
  assert (fIntegrator);

  // user cuts in W,Q2
  fWmin  = fConfig->GetDoubleDef("Wmin", - 1.0);
  fWmax  = fConfig->GetDoubleDef("Wmax",   1e9);
  fQ2min = fConfig->GetDoubleDef("Q2min", -1.0);
  fQ2max = fConfig->GetDoubleDef("Q2max",  1e9);

  // get upper E limit on res xsec spline (=f(E)) before assuming xsec=const
  fEMax = fConfig->GetDoubleDef("ESplineMax", 40);
  fEMax = TMath::Max(fEMax,10.); // don't accept user Emax if less than 10 GeV

  // create the baryon resonance list specified in the config.
  fResList.Clear();
  string resonances = fConfig->GetStringDef(
                   "resonance-name-list", gc->GetString("ResonanceNameList"));
  fResList.DecodeFromNameList(resonances);

  //-- Use algorithm within a DIS/RES join scheme. If yes get Wcut
  fUsingDisResJoin = fConfig->GetBoolDef("use-dis-res-joining-scheme", false);
  fWcut = 999999;
  if(fUsingDisResJoin) {
    fWcut = fConfig->GetDoubleDef("Wcut",gc->GetDouble("Wcut"));
  }
}
//____________________________________________________________________________

