//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

         Marco Roda <mroda \at liverpool.ac.uk>
         University of Liverpool

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
 @ This revision branches the old QPMDISPXSec in order to disentangle the pure DIS
   from Pachos et al. adn the internal GENIE tuning

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
#include "Framework/Messenger/Messenger.h"
#include "Physics/DeepInelastic/XSection/KNOTunedQPMDISPXSec.h"
#include "Framework/ParticleData/PDGCodes.h"
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
KNOTunedQPMDISPXSec::KNOTunedQPMDISPXSec() :
XSecAlgorithmI("genie::KNOTunedQPMDISPXSec") { ; }
//____________________________________________________________________________
KNOTunedQPMDISPXSec::KNOTunedQPMDISPXSec(string config) :
XSecAlgorithmI("genie::KNOTunedQPMDISPXSec", config) { ; }
//____________________________________________________________________________
KNOTunedQPMDISPXSec::~KNOTunedQPMDISPXSec()
{

}
//____________________________________________________________________________
double KNOTunedQPMDISPXSec::XSec( const Interaction * interaction,
		                          KinePhaseSpace_t kps) const
{

  double xsec = fDISModel -> XSec(interaction, kps) ;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pINFO)
        << "d2xsec/dxdy[FreeN] (E= " << E 
                      << ", x= " << x << ", y= " << y << ") = " << xsec;
#endif

  double R = this->DISRESJoinSuppressionFactor(interaction);

  xsec*=R;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DISPXSec", pINFO) << "D/R Join scheme - suppression factor R = " << R;;
     LOG("DISPXSec", pINFO) << "d2xsec/dxdy[FreeN, D/R Join] " << xsec;
#endif

  xsec = TMath::Max(0., xsec ) ;

  return xsec;
}
//____________________________________________________________________________
double KNOTunedQPMDISPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool KNOTunedQPMDISPXSec::ValidProcess(const Interaction * interaction) const
{
	return fDISModel -> ValidProcess(interaction) ;
}
//____________________________________________________________________________
double KNOTunedQPMDISPXSec::DISRESJoinSuppressionFactor(
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
void KNOTunedQPMDISPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KNOTunedQPMDISPXSec::Configure(string config)
{
  Algorithm::Configure(config);

  this->LoadConfig();
}
//____________________________________________________________________________
void KNOTunedQPMDISPXSec::LoadConfig(void)
{

  fHadronizationModel = nullptr ;

  fHadronizationModel =
    dynamic_cast<const KNOHadronization *> (this->SubAlg("Hadronizer"));
  assert(fHadronizationModel);

  GetParam( "Wcut", fWcut ) ;

  if ( fWcut <= 0. ) {

	LOG("KNOTunedQPMDISPXSec", pFATAL)
        << "Input configuration value for Wcut is not physical: Exiting" ;

      // From the FreeBSD Library Functions Manual
      //
      // EX_CONFIG (78)   Something was found in an unconfigured or miscon-
      //                  figured state.

    exit( 78 ) ;

  }

  // load the pure E.A.Paschos and J.Y.Yu model
  fDISModel = nullptr ;
  fDISModel = dynamic_cast<const QPMDISPXSec *>(this->SubAlg("DISModel") ) ;
  assert( fDISModel ) ;

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  // Caching the reduction factors used in the DIS/RES joing scheme?
  // In normal event generation (1 config -> many calls) it is worth caching
  // these suppression factors.
  // Depending on the way this algorithm is used during event reweighting,
  // precomputing (for all W's) & caching these factors might not be efficient.
  // Here we provide the option to turn the caching off at run-time (default: on)

  bool cache_enabled = RunOpt::Instance()->BareXSecPreCalc();

  GetParamDef( "UseCache", fUseCache, true ) ;
  fUseCache = fUseCache && cache_enabled;


  }
//____________________________________________________________________________

