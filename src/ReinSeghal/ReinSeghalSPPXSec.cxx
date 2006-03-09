//____________________________________________________________________________
/*!

\class    genie::ReinSeghalSPPXSec

\brief    Computes the cross section for an exclusive 1pi reaction through
          resonance neutrinoproduction according to the Rein-Seghal model.

          This algorithm produces in principle what you could also get from 
          the genie::RESXSec algorithm (RES cross section integrator) by 
          specifying the genie::ReinSeghalSPPPXSec as the differential 
          (d2xsec/fQ2dW) cross section model. However, ReinSeghalSPPXSec
          offers a faster alternative. Before computing any SPP cross section
          this algorithm computes and caches splines for resonance neutrino-
          production cross sections. This improves the speed since it is 
          reducing the number of calculations (the generic algorithm needs to
          recompute resonance production xsec for every exclusive channel).

          In this algorithm we follow the non-coherent approach: we sum
          the weighted resonance production cross sections rather than the
          resonance production amplitudes.

          Is a concrete implementation of the XSecAlgorithmI interface.\n

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  March 09, 2006

*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "CrossSections/GXSecFunc.h"
#include "Interaction/SppChannel.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "ReinSeghal/ReinSeghalSPPXSec.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Cache.h"
#include "Utils/CacheBranchFx.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSeghalSPPXSec::ReinSeghalSPPXSec() :
XSecAlgorithmI("genie::ReinSeghalSPPXSec")
{

}
//____________________________________________________________________________
ReinSeghalSPPXSec::ReinSeghalSPPXSec(string config) :
XSecAlgorithmI("genie::ReinSeghalSPPXSec", config)
{

}
//____________________________________________________________________________
ReinSeghalSPPXSec::~ReinSeghalSPPXSec()
{

}
//____________________________________________________________________________
double ReinSeghalSPPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghalSpp", pDEBUG) << *fConfig;

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //-- Get 1pi exclusive channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  assert(1);

  Cache * cache = Cache::Instance();

  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();
  const Target &       target     = init_state.GetTarget();

  // Get neutrino energy in the struck nucleon rest frame
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);

  double xsec = 0;

  unsigned int nres = fResList.NResonances();
  for(unsigned int ires = 0; ires < nres; ires++) {

     //-- Get next resonance from the resonance list
     Resonance_t res = fResList.ResonanceId(ires);

     InteractionType_t it = proc_info.InteractionTypeId();
     int               iN = target.StruckNucleonPDGCode();

     //-- Build a unique name for the cache branch

     string res_name = utils::res::AsString(res);
     int    nu_code  = init_state.GetProbePDGCode();
     string it_name  = InteractionType::AsString(it);
     string nc_nuc   = "";
     if(it == kIntWeakNC) { nc_nuc = ((iN==kPdgProton) ? "p" : "n"); }

     ostringstream intk;
     intk << "ResExcitationXSec/R:" << res_name << ";nu:"  << nu_code
           << ";int:" << it_name << nc_nuc;
      
     string algkey = this->Id().Key();
     string ikey   = intk.str();
     string key    = cache->CacheBranchKey(algkey, ikey);

     LOG("Kinematics", pINFO) << "Finding cache branch with key: " << key;
     CacheBranchFx * cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));

     if(!cache_branch) {
       LOG("ReinSeghalSpp", pWARN)  
         << "No cached RES v-production data for input neutrino"
         << " (pdgc: " << nu_code << ")";
       LOG("ReinSeghalSpp", pWARN)  
         << "Wait while computing/caching RES production xsec first...";

       this->CacheResExcitationXSec(interaction); 

       LOG("ReinSeghalSpp", pINFO) << "Done caching resonance xsec data";
       LOG("ReinSeghalSpp", pINFO) 
               << "Finding newly created cache branch with key: " << key;
       cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
       assert(cache_branch);
     }

     //-- Get cached resonance neutrinoproduction xsec
     //   (If E>Emax, assume xsec = xsec(Emax) - but do not evaluate the
     //    cross section spline at the end of its energy range-)
     double rxsec = (Ev<fEMax-1) ? (*cache_branch)(Ev) : (*cache_branch)(fEMax-1);

     //-- Get the BR for the (resonance) -> (exclusive final state)
     double br = SppChannel::BranchingRatio(spp_channel, res);

     //-- Get the Isospin Glebsch-Gordon coefficient for the given resonance
     //   and exclusive final state
     double igg = SppChannel::IsospinWeight(spp_channel, res);

     //-- Compute the weighted xsec
     //  (total weight = Breit-Wigner * BR * isospin Glebsch-Gordon)
     double res_xsec_contrib = rxsec*br*igg;

     SLOG("ReinSeghalSpp", pINFO)
       << "Contrib. from [" << utils::res::AsString(res) << "] = "
       << "<Glebsch-Gordon = " << igg
       << "> * <BR(->1pi) = " << br
       << "> * <Breit-Wigner * d^2xsec/dQ^2dW = " << rxsec
       << "> = " << res_xsec_contrib;
   
     //-- Add contribution of this resonance to the cross section
     xsec += res_xsec_contrib;

  }//res

  SLOG("ReinSeghalSpp", pNOTICE)  
                         << "XSec[RES] (Ev = " << Ev << " GeV) = " << xsec;
  return xsec;
}
//____________________________________________________________________________
void ReinSeghalSPPXSec::CacheResExcitationXSec(const Interaction * in) const
{
  Cache * cache = Cache::Instance();

  Interaction * interaction = new Interaction(*in);

  const int kNSt = 3;
  const InteractionType_t it[kNSt] = {kIntWeakCC, kIntWeakNC, kIntWeakNC };
  const int iN[kNSt] = {kPdgProton, kPdgProton, kPdgNeutron};

  const double kEmin         = 0.120; // xsec_res(Emin) = 0
  const double kLogEmin      = TMath::Log(kEmin);
  const double kLogEmax      = TMath::Log(fEMax);
  const int    kMinNKnots    = (int) (10*(kLogEmax-kLogEmin)); // at least 10 knots per decade...
  const int    kNSplineKnots = TMath::Max(40, kMinNKnots); // but at least 40 knots
  const double kdLogE        = (kLogEmax-kLogEmin)/(kNSplineKnots-1);

  unsigned int nres = fResList.NResonances();
  for(int i=0; i<kNSt; i++) {

    interaction->GetInitialStatePtr()->GetTargetPtr()->SetStruckNucleonPDGCode(iN[i]);
    interaction->GetProcessInfoPtr()->Set(kScResonant, it[i]);

    for(unsigned int ires = 0; ires < nres; ires++) {

      //-- Get next resonance from the resonance list
      Resonance_t res = fResList.ResonanceId(ires);

      interaction->GetExclusiveTagPtr()->SetResonance(res);

      //-- Build a unique name for the cache branch

      string res_name = utils::res::AsString(res);
      int    nu_code  = interaction->GetInitialState().GetProbePDGCode();
      string it_name  = InteractionType::AsString(it[i]);
      string nc_nuc   = "";
      if(it[i] == kIntWeakNC) { nc_nuc = ((iN[i]==kPdgProton) ? "p" : "n"); }

      ostringstream intk;
      intk << "ResExcitationXSec/R:" << res_name << ";nu:"  << nu_code
           << ";int:" << it_name << nc_nuc;
      
      string algkey = this->Id().Key();
      string ikey   = intk.str();
      string key    = cache->CacheBranchKey(algkey, ikey);

      //-- Make sure the cache branch does not already exists
      CacheBranchFx * cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
      assert(!cache_branch);

      //-- Create the new cache branch
      LOG("ReinSeghalSpp", pNOTICE) 
                        << "\n ** Creating cache branch - key = " << key;
      cache_branch = new CacheBranchFx("RES Excitation XSec");
      cache->AddCacheBranch(key, cache_branch);
      assert(cache_branch);

      TLorentzVector p4(0,0,0,0);

      for(int ie=0; ie<kNSplineKnots; ie++) {

        double Ev = TMath::Exp(kLogEmin + ie*kdLogE);
        p4.SetPxPyPzE(0,0,Ev,Ev);
        interaction->GetInitialStatePtr()->SetProbeP4(p4);

        // Get W integration range
        Range1D_t rW = this->WRange(interaction);
        // Get the wider possible Q2 range for the input W range
        Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction, rW);

        LOG("ReinSeghalSpp", pINFO) 
	  << "*** Integrating d^2 XSec/dWdQ^2 for R: " 
     	                 << utils::res::AsString(res) << " at Ev = " << Ev;
        LOG("ReinSeghalSpp", pINFO) << "{W}   = " << rW.min  << ", " << rW.max;
	LOG("ReinSeghalSpp", pINFO) << "{Q^2} = " << rQ2.min << ", " << rQ2.max;

        double xsec = 0;

        if(rW.max<rW.min || rQ2.max<rQ2.min || rW.min<0 || rQ2.min<0) {
	  LOG("ReinSeghalSpp", pINFO) << "** Not allowed kinematically, xsec=0";
        } else {
          GXSecFunc * func = new Integrand_D2XSec_DWDQ2_E(
                                           fSingleResXSecModel, interaction);
          func->SetParam(0,"W",  rW);
          func->SetParam(1,"Q2", rQ2);
          xsec = fIntegrator->Integrate(*func);
          delete func;
        }
        cache_branch->AddValues(Ev,xsec);
        LOG("ReinSeghalSpp", pNOTICE) 
               << "RES XSec (R:" << res_name << ", E=" << Ev << ")= " << xsec;
      }//ie

      cache_branch->CreateSpline();

    }//ires
  }//i

  delete interaction;
}
//____________________________________________________________________________
bool ReinSeghalSPPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();

  if(!proc_info.IsResonant()) return false;
  if(!proc_info.IsWeak())     return false;

  int  nuc = init_state.GetTarget().StruckNucleonPDGCode();
  int  nu  = init_state.GetProbePDGCode();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
bool ReinSeghalSPPXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);

  double EvThr = utils::kinematics::EnergyThreshold(interaction);
  if(Ev <= EvThr) return false;

  return true;
}
//____________________________________________________________________________
void ReinSeghalSPPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSeghalSPPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSeghalSPPXSec::LoadConfig(void)
{
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
  fWminCut  = fConfig->GetDoubleDef("Wmin", -1.0);
  fWmaxCut  = fConfig->GetDoubleDef("Wmax",  1e9);
  fQ2minCut = fConfig->GetDoubleDef("Q2min", -1.0);
  fQ2maxCut = fConfig->GetDoubleDef("Q2max",  1e9);

  // get upper E limit on res xsec spline (=f(E)) before assuming xsec=const
  fEMax = fConfig->GetDoubleDef("ESplineMax", 40);
  fEMax = TMath::Max(fEMax,10.); // don't accept user Emax if less than 10 GeV

  // get list of all the resonances that the user wants to consider
  this->GetResonanceList();
}
//____________________________________________________________________________
void ReinSeghalSPPXSec::GetResonanceList(void)
{
// create the baryon resonance list specified in the config.

  fResList.Clear();

  assert( fConfig->Exists("resonance-name-list") );
  string resonances = fConfig->GetString("resonance-name-list");
  fResList.DecodeFromNameList(resonances);
}
//____________________________________________________________________________
Range1D_t ReinSeghalSPPXSec::WRange(const Interaction * interaction) const
{
  //-- Get the physically allowed W range for this interaction and allow the
  //   user inputs (if any) to narrow it

  Range1D_t rW = utils::kinematics::WRange(interaction); // physical range

  LOG("ReinSeghalSpp", pDEBUG)
       << "Physical W range: " << "[" << rW.min << ", " << rW.max << "] GeV";

  // apply user cuts
  utils::kinematics::ApplyCutsToKineLimits(rW, fWminCut,  fWmaxCut );

  LOG("ReinSeghalSpp", pDEBUG)
       << "Physical & User W range: "
                               << "[" << rW.min << ", " << rW.max << "] GeV";
  return rW;
}
//___________________________________________________________________________
Range1D_t ReinSeghalSPPXSec::Q2Range(const Interaction * interaction) const
{
  //-- Get the physically allowed Q2 range for this interaction and allow the
  //   user inputs (if any) to narrow it

  Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction); // physical range

  LOG("ReinSeghalSpp", pDEBUG) << "Physical Q2 range: "
                         << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";

  // apply user cuts
  utils::kinematics::ApplyCutsToKineLimits(rQ2, fQ2minCut, fQ2maxCut);

  LOG("ReinSeghalSpp", pDEBUG)
       << "Physical & User Q2 range: "
                         << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
  return rQ2;
}
//___________________________________________________________________________


