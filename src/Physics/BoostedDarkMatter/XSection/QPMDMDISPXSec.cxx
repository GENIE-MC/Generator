//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joshua Berger <jberger \at physics.wisc.edu>
         University of Wisconsin-Madison

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 
*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>
#include <TH1D.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Units.h"
#include "Physics/DeepInelastic/XSection/DISStructureFuncModelI.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Hadronization/HadronizationModelI.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/BoostedDarkMatter/XSection/QPMDMDISPXSec.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"

using std::ostringstream;

using namespace genie;
using namespace genie::constants;
//using namespace genie::units;

//____________________________________________________________________________
QPMDMDISPXSec::QPMDMDISPXSec() :
XSecAlgorithmI("genie::QPMDMDISPXSec")
{
  fInInitPhase = true;
}
//____________________________________________________________________________
QPMDMDISPXSec::QPMDMDISPXSec(string config) :
XSecAlgorithmI("genie::QPMDMDISPXSec", config)
{
  fInInitPhase = true;
}
//____________________________________________________________________________
QPMDMDISPXSec::~QPMDMDISPXSec()
{

}
//____________________________________________________________________________
double QPMDMDISPXSec::XSec(
     const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Get kinematical & init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  // const ProcessInfo &  proc_info  = interaction -> ProcInfo(); // comment-out unused variable to eliminate warnings

  LOG("DMDISPXSec", pDEBUG) << "Using v^" << fVelMode << " dependence";
  
  double E     = init_state.ProbeE(kRfHitNucRest);
  double ml    = interaction->FSPrimLepton()->Mass();
  double Mnuc  = init_state.Tgt().HitNucMass();
  double x     = kinematics.x();
  double y     = kinematics.y();

  double E2    = E    * E;
  double ml2   = ml   * ml;
  // double ml4   = ml2  * ml2; // comment-out unused variable to eliminate warnings
  // double Mnuc2 = Mnuc * Mnuc; // comment-out unused variable to eliminate warnings

  #ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMDISPXSec", pDEBUG)  
   << "Computing d2xsec/dxdy @ E = " << E << ", x = " << x << ", y = " << y;
  #endif

  // One of the xsec terms changes sign for antineutrinos @ DMDIS/CC

  // bool is_nubar_cc = pdg::IsAntiNeutrino(init_state.ProbePdg()) && 
  //                    proc_info.IsWeakCC(); // // comment-out unused variable to eliminate warnings
  // int sign = (is_nubar_cc) ? -1 : 1; // comment-out unused variable to eliminate warnings

  // Calculate the DMDIS structure functions
  fDISSF.Calculate(interaction); 

  #ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMDISPXSec", pDEBUG) << fDISSF;
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
  // double Q4 = Q2*Q2; // comment-out unused variable to eliminate warnings
  // temp: set the Z' mass to MZ and g' = 1 for now
  LOG("DMDISPXSec", pDEBUG)
    << "Using a mediator mass " << fMedMass;
  double Mzp2 = TMath::Power(fMedMass,2);
  // double gzp = RunOpt::Instance()->ZpCoupling();
  double gzp = fgzp;
  double gzp4 = TMath::Power(gzp,4);
  g2 = gzp4 / TMath::Power((Q2 + Mzp2), 2);
  double p2 = TMath::Max(E2 - ml2,0.);
  double front_factor = (g2*Mnuc*E) / kPi * (E2 / p2);
  
  // Build all dxsec/dxdy terms
  double term1 = 0.;
  double term2 = 0.;
  double term3 = 0.;
  double term4 = 0.;
  double term5 = 0.;
  // The cross-check of these expressions is that they should
  // give the elastic cross-section in the limit x -> 1, PDF -> 1,
  // and absent nuclear effects
  switch (fVelMode) {
  case 0:
    // Second lines contain longitudinal Z' coupling
    // If the mediator is relatively light, these terms are important
    // and can't be neglected like they are in the SM
    term1  = 0.125 * y * ( x*y + 3.*ml2/(E*Mnuc) );
    term1 += ml2*x*y*y/(2*Mzp2) *  (1. + Mnuc*E*x*y / Mzp2);
    term2  = 0.125 * (1 - y - Mnuc*x*y/(2*E) - ml2/E2);
    term2 += ml2*y*y/(4*Mzp2) * (1. + Mnuc*E*x*y / Mzp2);
    term4  = x*y*ml2/(4*Mnuc*E);
    term4 += ml2*x*x*y*y/Mzp2 * (1. + Mnuc*E*x*y / Mzp2);
    term5  = -y*ml2/(8*Mnuc*E);
    term5 += -ml2*x*y*y/(2*Mzp2) * (1. + Mnuc*E*x*y / Mzp2);
    break;
  case 2:
    // Scalar case has no longitudinal Z' coupling
    term1 = -0.125 * y * (0.5 * x * y + ml2/(E*Mnuc));
    term2 = 0.125 * (1. - y + 0.25 * y * y);
  }  
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMDISPXSec", pDEBUG)  
    << "\nd2xsec/dxdy ~ (" << term1 << ")*F1+(" << term2 << ")*F2+(" 
                  << term3 << ")*F3+(" << term4 << ")*F4+(" << term5 << ")*F5";
#endif

  term1 *= fDISSF.F1();
  term2 *= fDISSF.F2();
  term3 *= fDISSF.F3();
  term4 *= fDISSF.F4();
  term5 *= fDISSF.F5();

  LOG("DMDISPXSec", pDEBUG)  
    << "\nd2xsec/dxdy ~ (" << term1 << ")+(" << term2 << ")+(" 
                  << term3 << ")+(" << term4 << ")+(" << term5 << ")";

  
  double xsec = front_factor * (term1 + term2 + term3 + term4 + term5);
  xsec = TMath::Max(xsec,0.);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMDISPXSec", pINFO)
        << "d2xsec/dxdy[FreeN] (E= " << E 
                      << ", x= " << x << ", y= " << y << ") = " << xsec;
#endif

  // If the DMDIS/RES joining scheme is enabled, modify the xsec accordingly
  if(fUsingDisResJoin) {
     double R = this->DMDISRESJoinSuppressionFactor(interaction);
     xsec*=R;
     #ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DMDISPXSec", pINFO) << "D/R Join scheme - suppression factor R = " << R;;
     LOG("DMDISPXSec", pINFO) << "d2xsec/dxdy[FreeN, D/R Join] " << xsec;
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
  LOG("DMDISPXSec", pINFO) 
       << "Subtracting charm piece: " << xsec_charm << " / out of " << xsec;
#endif
  xsec = TMath::Max(0., xsec-xsec_charm);
  return xsec;
}
//____________________________________________________________________________
double QPMDMDISPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool QPMDMDISPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsDarkMatterDeepInelastic()) return false;

  const InitialState & init_state = interaction -> InitState();
  int probe_pdg = init_state.ProbePdg();
  if(!pdg::IsDarkMatter(probe_pdg)) return false;

  if(! init_state.Tgt().HitNucIsSet()) return false;

  int hitnuc_pdg = init_state.Tgt().HitNucPdg();
  if(!pdg::IsNeutronOrProton(hitnuc_pdg)) return false;

  return true;
}
//____________________________________________________________________________
double QPMDMDISPXSec::DMDISRESJoinSuppressionFactor(
   const Interaction * in) const
{
// Computes suppression factors for the DMDIS xsec under the used DMDIS/RES join
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
    // algid/DMDIS-RES-Join/nu-pdg:N;hit-nuc-pdg:N/inttype
    Cache * cache = Cache::Instance();
    string algkey = this->Id().Key() + "/DMDIS-RES-Join";

    ostringstream ikey;
    ikey << "nu-pdgc:" << ist.ProbePdg() 
         << ";hit-nuc-pdg:"<< ist.Tgt().HitNucPdg() << "/"
         << pi.InteractionTypeAsString();

    string key = cache->CacheBranchKey(algkey, ikey.str());

    CacheBranchFx * cbr =
          dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));

    // If it does't exist then create a new one 
    // and cache DMDIS xsec suppression factors
    bool non_zero=false;
    if(!cbr) {
      LOG("DMDISXSec", pNOTICE) 
                        << "\n ** Creating cache branch - key = " << key;

      cbr = new CacheBranchFx("DMDIS Suppr. Factors in DMDIS/RES Join Scheme");
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
        LOG("DMDISXSec", pNOTICE) 
	    << "Cached DMDIS XSec Suppr. factor (@ W=" << W << ") = " << R;

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

  LOG("DMDISXSec", pDEBUG) 
      << "DMDIS/RES Join: DMDIS xsec suppr. (W=" << Wo << ") = " << Ro;

  return Ro;
}
//____________________________________________________________________________
void QPMDMDISPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDMDISPXSec::Configure(string config)
{
  Algorithm::Configure(config);

  Registry r( "QPMDMDISPXSec_specific", false ) ;

  std::string  key = "XSecModel@genie::EventGenerator/DIS-CC-CHARM" ;
  RgKey local_key = "CharmXSec" ;
  r.Set( local_key , AlgConfigPool::Instance() -> GlobalParameterList() -> GetAlg(key) ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDMDISPXSec::LoadConfig(void)
{
  // Access global defaults to use in case of missing parameters

  fDISSFModel = 0;
  fDISSFModel = 
     dynamic_cast<const DISStructureFuncModelI *> (this->SubAlg("SFAlg"));
  assert(fDISSFModel);

  fDISSF.SetModel(fDISSFModel); // <-- attach algorithm

  this->GetParam( "UseDRJoinScheme", fUsingDisResJoin ) ;

  fHadronizationModel = 0;
  fWcut = 0.;

  if(fUsingDisResJoin) {
     fHadronizationModel = 
       dynamic_cast<const HadronizationModelI *> (this->SubAlg("Hadronizer"));
     assert(fHadronizationModel);

     // Load Wcut determining the phase space area where the multiplicity prob.
     // scaling factors would be applied -if requested-
     this->GetParam( "Wcut", fWcut ) ;
  }

  // Cross section scaling factor
  this->GetParam( "DIS-XSecScale", fScale ) ;

  // sin^4(theta_weinberg)
  double thw  ;
  this->GetParam( "WeinbergAngle", thw ) ;
  fSin48w = TMath::Power( TMath::Sin(thw), 4 );

  // Caching the reduction factors used in the DMDIS/RES joing scheme?
  // In normal event generation (1 config -> many calls) it is worth caching
  // these suppression factors.
  // Depending on the way this algorithm is used during event reweighting,
  // precomputing (for all W's) & caching these factors might not be efficient.
  // Here we provide the option to turn the caching off at run-time (default: on)

  bool cache_enabled = RunOpt::Instance()->BareXSecPreCalc();

  this->GetParamDef( "UseCache", fUseCache, true ) ;
  fUseCache = fUseCache && cache_enabled;

  // Since this method would be called every time the current algorithm is 
  // reconfigured at run-time, remove all the data cached by this algorithm
  // since they depend on the previous configuration

  if(!fInInitPhase) {
     Cache * cache = Cache::Instance();
     string keysubstr = this->Id().Key() + "/DMDIS-RES-Join";
     cache->RmMatchedCacheBranches(keysubstr);
  }
  fInInitPhase = false;

  // velocity dependence of the interaction
  this->GetParamDef("velocity-mode", fVelMode, 0);

  // mediator coupling
  this->GetParam("ZpCoupling", fgzp);

  // mediator mass ratio and mediator mass
  fMedMass = PDGLibrary::Instance()->Find(kPdgMediator)->Mass();
  
  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  RgKey xdefkey = "CharmXSec" ;
  RgAlg xalg ;
  GetParam( xdefkey, xalg) ;
  LOG("DMDISXSec", pDEBUG)
     << "Loading the cross section model: " << xalg;
  fCharmProdModel = dynamic_cast<const XSecAlgorithmI *> ( this -> SubAlg(xdefkey) ) ;
  assert(fCharmProdModel);
}
//____________________________________________________________________________

