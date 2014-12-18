//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 09, 2009 - CA
   Modified to handle charged lepton - nucleon(nucleus) scattering.
   Renamed QPMDISPXSec from DISPartonModelPXSec following code reorganization.
 @ Oct 11, 2009 - CA
   Implemented ValidProcess()
 @ Jan 29, 2013 - CA
   Don't look-up depreciated $GDISABLECACHING environmental variable.
   Use the RunOpt singleton instead.

*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>
#include <TH1D.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Units.h"
#include "Base/DISStructureFuncModelI.h"
#include "Base/XSecIntegratorI.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "Fragmentation/HadronizationModelI.h"
#include "Messenger/Messenger.h"
#include "PartonModel/QPMDISPXSec.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/RunOpt.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"
#include "Utils/KineUtils.h"
#include "Utils/Cache.h"
#include "Utils/CacheBranchFx.h"

using std::ostringstream;

using namespace genie;
using namespace genie::constants;
using namespace genie::units;

//____________________________________________________________________________
QPMDISPXSec::QPMDISPXSec() :
XSecAlgorithmI("genie::QPMDISPXSec")
{
  fInInitPhase = true;
}
//____________________________________________________________________________
QPMDISPXSec::QPMDISPXSec(string config) :
XSecAlgorithmI("genie::QPMDISPXSec", config)
{
  fInInitPhase = true;
}
//____________________________________________________________________________
QPMDISPXSec::~QPMDISPXSec()
{

}
//____________________________________________________________________________
double QPMDISPXSec::XSec(
     const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Get kinematical & init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();

  double E     = init_state.ProbeE(kRfHitNucRest);
  double ml    = interaction->FSPrimLepton()->Mass();
  double Mnuc  = init_state.Tgt().HitNucMass();
  double x     = kinematics.x();
  double y     = kinematics.y();

  double E2    = E    * E;
  double ml2   = ml   * ml;
  double ml4   = ml2  * ml2;
  double Mnuc2 = Mnuc * Mnuc;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pDEBUG)  
   << "Computing d2xsec/dxdy @ E = " << E << ", x = " << x << ", y = " << y;
#endif

  // One of the xsec terms changes sign for antineutrinos @ DIS/CC

  bool is_nubar_cc = pdg::IsAntiNeutrino(init_state.ProbePdg()) && 
                     proc_info.IsWeakCC();
  int sign = (is_nubar_cc) ? -1 : 1;

  // Calculate the DIS structure functions
  fDISSF.Calculate(interaction); 

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pDEBUG) << fDISSF;
#endif

  //
  // Compute the differential cross section
  //

  double g2 = kGF2;
  // For EM interaction replace  G_{Fermi} with :
  // a_{em} * pi / ( sqrt(2) * sin^2(theta_weinberg) * Mass_{W}^2 }
  // See C.Quigg, Gauge Theories of the Strong, Weak and E/M Interactions,
  // ISBN 0-8053-6021-2, p.112 (6.3.57)
  // Also, take int account that the photon propagator is 1/p^2 but the
  // W propagator is 1/(p^2-Mass_{W}^2), so weight the EM case with
  // Mass_{W}^4 / q^4
  // So, overall:
  // G_{Fermi}^2 --> a_{em}^2 * pi^2 / (2 * sin^4(theta_weinberg) * q^{4})
  //
  double Q2 = utils::kinematics::XYtoQ2(E,Mnuc,x,y);
  double Q4 = Q2*Q2;
  if(proc_info.IsEM()) {
    g2 = kAem2 * kPi2 / (2.0 * fSin48w * Q4); 
  }
#ifdef __GENIE_VHE_ENABLED__  
  if (proc_info.IsWeakCC()) {
    g2 = kGF2 * kMw2 * kMw2 / TMath::Power((Q2 + kMw2), 2);
  } else if (proc_info.IsWeakNC()) {
    g2 = kGF2 * kMz2 * kMz2 / TMath::Power((Q2 + kMz2), 2);
  }
#endif
  double front_factor = (g2*Mnuc*E) / kPi;

  // Build all dxsec/dxdy terms
  double term1 = y * ( x*y + ml2/(2*E*Mnuc) );
  double term2 = 1 - y - Mnuc*x*y/(2*E) - ml2/(4*E2);
  double term3 = sign * (x*y*(1-y/2) - y*ml2/(4*Mnuc*E));
  double term4 = x*y*ml2/(2*Mnuc*E) + ml4/(4*Mnuc2*E2);
  double term5 = -1.*ml2/(2*Mnuc*E);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pDEBUG)  
    << "\nd2xsec/dxdy ~ (" << term1 << ")*F1+(" << term2 << ")*F2+(" 
                  << term3 << ")*F3+(" << term4 << ")*F4+(" << term5 << ")*F5";
#endif

  term1 *= fDISSF.F1();
  term2 *= fDISSF.F2();
  term3 *= fDISSF.F3();
  term4 *= fDISSF.F4();
  term5 *= fDISSF.F5();

  double xsec = front_factor * (term1 + term2 + term3 + term4 + term5);
  xsec = TMath::Max(xsec,0.);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pINFO)
        << "d2xsec/dxdy[FreeN] (E= " << E 
                      << ", x= " << x << ", y= " << y << ") = " << xsec;
#endif

  // If the DIS/RES joining scheme is enabled, modify the xsec accordingly
  if(fUsingDisResJoin) {
     double R = this->DISRESJoinSuppressionFactor(interaction);
     xsec*=R;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DISPXSec", pINFO) << "D/R Join scheme - suppression factor R = " << R;;
     LOG("DISPXSec", pINFO) << "d2xsec/dxdy[FreeN, D/R Join] " << xsec;
#endif
  }

  // The algorithm computes d^2xsec/dxdy
  // Check whether variable tranformation is needed
  if(kps!=kPSxyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSxyfE,kps);
    xsec *= J;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Compute nuclear cross section (simple scaling here, corrections must
  // have been included in the structure functions)
  const Target & target = init_state.Tgt();
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 
  xsec *= NNucl; 

  // Apply scaling / if required to reach well known asymmptotic value
  xsec *= fScale;

  // Subtract the inclusive charm production cross section
  interaction->ExclTagPtr()->SetCharm();
  double xsec_charm = fCharmProdModel->XSec(interaction,kps);
  interaction->ExclTagPtr()->UnsetCharm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pINFO) 
       << "Subtracting charm piece: " << xsec_charm << " / out of " << xsec;
#endif
  xsec = TMath::Max(0., xsec-xsec_charm);
  return xsec;
}
//____________________________________________________________________________
double QPMDISPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool QPMDISPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsDeepInelastic()) return false;

  const InitialState & init_state = interaction -> InitState();
  int probe_pdg = init_state.ProbePdg();
  if(!pdg::IsLepton(probe_pdg)) return false;

  if(! init_state.Tgt().HitNucIsSet()) return false;

  int hitnuc_pdg = init_state.Tgt().HitNucPdg();
  if(!pdg::IsNeutronOrProton(hitnuc_pdg)) return false;

  return true;
}
//____________________________________________________________________________
double QPMDISPXSec::DISRESJoinSuppressionFactor(
   const Interaction * in) const
{
// Computes suppression factors for the DIS xsec under the used DIS/RES join
// scheme. Since this is a 'low-level' algorithm that is being called many
// times per generated event or computed cross section spline, the suppression
// factors would be cached to avoid calling the hadronization model too often.
//
  double R=0, Ro=0;

  const double Wmin = kNeutronMass + kPionMass + 1E-3;

  const InitialState & ist = in->InitState();
  const ProcessInfo &  pi  = in->ProcInfo();

  double E    = ist.ProbeE(kRfHitNucRest);
  double Mnuc = ist.Tgt().HitNucMass();
  double x    = in->Kine().x(); 
  double y    = in->Kine().y();
  double Wo   = utils::kinematics::XYtoW(E,Mnuc,x,y);

  TH1D * mprob = 0;

  if(!fUseCache) {
    // ** Compute the reduction factor at each call - no caching 
    //
    mprob = fHadronizationModel->MultiplicityProb(in,"+LowMultSuppr");
    R = 1;
    if(mprob) {
       R = mprob->Integral("width");
       delete mprob;
    }
  }
  else {

    // ** Precompute/cache the reduction factors and then use the 
    // ** cache to evaluate these factors

    // Access the cache branch. The branch key is formed as:
    // algid/DIS-RES-Join/nu-pdg:N;hit-nuc-pdg:N/inttype
    Cache * cache = Cache::Instance();
    string algkey = this->Id().Key() + "/DIS-RES-Join";

    ostringstream ikey;
    ikey << "nu-pdgc:" << ist.ProbePdg() 
         << ";hit-nuc-pdg:"<< ist.Tgt().HitNucPdg() << "/"
         << pi.InteractionTypeAsString();

    string key = cache->CacheBranchKey(algkey, ikey.str());

    CacheBranchFx * cbr =
          dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));

    // If it does't exist then create a new one 
    // and cache DIS xsec suppression factors
    bool non_zero=false;
    if(!cbr) {
      LOG("DISXSec", pNOTICE) 
                        << "\n ** Creating cache branch - key = " << key;

      cbr = new CacheBranchFx("DIS Suppr. Factors in DIS/RES Join Scheme");
      Interaction interaction(*in);

      const int kN   = 300;
      double WminSpl = Wmin;
      double WmaxSpl = fWcut + 0.1; // well into the area where scaling factor = 1
      double dW      = (WmaxSpl-WminSpl)/(kN-1);

      for(int i=0; i<kN; i++) {
        double W = WminSpl+i*dW;
        interaction.KinePtr()->SetW(W);
        mprob = fHadronizationModel->MultiplicityProb(&interaction,"+LowMultSuppr");
        R = 1;
        if(mprob) {
           R = mprob->Integral("width");
           delete mprob;
        }
        // make sure that it takes enough samples where it is non-zero:
        // modify the step and the sample counter once I've hit the first
        // non-zero value
        if(!non_zero && R>0) {
	  non_zero=true;
          WminSpl=W;
          i = 0;
          dW = (WmaxSpl-WminSpl)/(kN-1);
        }
        LOG("DISXSec", pNOTICE) 
	    << "Cached DIS XSec Suppr. factor (@ W=" << W << ") = " << R;

        cbr->AddValues(W,R);
      }
      cbr->CreateSpline();

      cache->AddCacheBranch(key, cbr);
      assert(cbr);
    } // cache data

    // get the reduction factor from the cache branch
    if(Wo > Wmin && Wo < fWcut-1E-2) {
       const CacheBranchFx & cache_branch = (*cbr);
       R = cache_branch(Wo);
    }
  }

  // Now return the suppression factor
  if      (Wo > Wmin && Wo < fWcut-1E-2) Ro = R;
  else if (Wo <= Wmin)                   Ro = 0.0;
  else                                   Ro = 1.0;

  LOG("DISXSec", pDEBUG) 
      << "DIS/RES Join: DIS xsec suppr. (W=" << Wo << ") = " << Ro;

  return Ro;
}
//____________________________________________________________________________
void QPMDISPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDISPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDISPXSec::LoadConfig(void)
{
  // Access global defaults to use in case of missing parameters
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fDISSFModel = 0;
  fDISSFModel = 
     dynamic_cast<const DISStructureFuncModelI *> (this->SubAlg("SFAlg"));
  assert(fDISSFModel);

  fDISSF.SetModel(fDISSFModel); // <-- attach algorithm

  fUsingDisResJoin = fConfig->GetBoolDef(
                       "UseDRJoinScheme", gc->GetBool("UseDRJoinScheme"));
  fHadronizationModel = 0;
  fWcut = 0.;

  if(fUsingDisResJoin) {
     fHadronizationModel = 
       dynamic_cast<const HadronizationModelI *> (this->SubAlg("Hadronizer"));
     assert(fHadronizationModel);

     // Load Wcut determining the phase space area where the multiplicity prob.
     // scaling factors would be applied -if requested-
     fWcut = fConfig->GetDoubleDef("Wcut", gc->GetDouble("Wcut"));
  }

  // Cross section scaling factor
  fScale = fConfig->GetDoubleDef("Scale", gc->GetDouble("DIS-XSecScale"));

  // sin^4(theta_weinberg)
  double thw = fConfig->GetDoubleDef(
         "weinberg-angle", gc->GetDouble("WeinbergAngle"));
     
  fSin48w = TMath::Power( TMath::Sin(thw), 4 );

  // Caching the reduction factors used in the DIS/RES joing scheme?
  // In normal event generation (1 config -> many calls) it is worth caching
  // these suppression factors.
  // Depending on the way this algorithm is used during event reweighting,
  // precomputing (for all W's) & caching these factors might not be efficient.
  // Here we provide the option to turn the caching off at run-time (default: on)

  bool cache_enabled = RunOpt::Instance()->BareXSecPreCalc();

  fUseCache = fConfig->GetBoolDef("UseCache", true);
  fUseCache = fUseCache && cache_enabled;

  // Since this method would be called every time the current algorithm is 
  // reconfigured at run-time, remove all the data cached by this algorithm
  // since they depend on the previous configuration

  if(!fInInitPhase) {
     Cache * cache = Cache::Instance();
     string keysubstr = this->Id().Key() + "/DIS-RES-Join";
     cache->RmMatchedCacheBranches(keysubstr);
  }
  fInInitPhase = false;

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  // Load the charm production cross section model
  RgKey xkey    = "CharmXSecModel";
  RgKey xdefkey = "XSecModel@genie::EventGenerator/DIS-CC-CHARM";
  RgAlg xalg    = fConfig->GetAlgDef(xkey, gc->GetAlg(xdefkey));
  LOG("DISXSec", pDEBUG)
     << "Loading the cross section model: " << xalg;
  fCharmProdModel = dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(xkey));
  assert(fCharmProdModel);
}
//____________________________________________________________________________

