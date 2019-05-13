//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Corey Reed <cjreed \at nikhef.nl>
         Nikhef - January 27, 2010

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TSpline.h>
#include <TGraph.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/InverseBetaDecay/XSection/Constants.h"
#include "Physics/InverseBetaDecay/XSection/KLVOxygenIBDPXSec.h"

using namespace genie;
using namespace genie::constants;

const double KLVOxygenIBDPXSec::kO16NubarThr = 0.0114; // GeV
const double KLVOxygenIBDPXSec::kO16NuMinE   = 0.0150; // GeV
const double KLVOxygenIBDPXSec::kMaxE        = 0.1000; // GeV

ClassImp(KLVOxygenIBDPXSec)

//____________________________________________________________________________
KLVOxygenIBDPXSec::KLVOxygenIBDPXSec() :
   XSecAlgorithmI("genie::KLVOxygenIBDPXSec"),
   fXsplNue(0),
   fXsplNuebar(0)
{

}
//____________________________________________________________________________
KLVOxygenIBDPXSec::KLVOxygenIBDPXSec(string config) :
   XSecAlgorithmI("genie::KLVOxygenIBDPXSec", config),
   fXsplNue(0),
   fXsplNuebar(0)
{

}
//____________________________________________________________________________
KLVOxygenIBDPXSec::~KLVOxygenIBDPXSec()
{
   delete fXsplNue;
   delete fXsplNuebar;
}
//____________________________________________________________________________
double KLVOxygenIBDPXSec::XSec(const Interaction * interaction,
                               KinePhaseSpace_t /* kps */) const
{
   // compute the differential cross section ds/dt
   // currently not implemented (only total)

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  LOG("KLVOxygen", pWARN)
     << "*** No differential cross section calculation is implemented yet";

  return 1;
}
//____________________________________________________________________________
void KLVOxygenIBDPXSec::MakeAntiNuESpline(void)
{
   // make a new xsec spline from the KLV paper's calculation
   // for the 16O + nu_e_bar reaction
   // remove any spline that might already exist

   delete fXsplNuebar;
   
   static const Int_t npts_nuebar = 21;
   static const Double_t Evnuebar[npts_nuebar] = {
      kO16NubarThr,
      15.0e-3, 17.5e-3, 20.0e-3, 22.5e-3,  25.0e-3,
      27.5e-3, 30.0e-3, 32.5e-3, 35.0e-3,  37.5e-3,
      40.0e-3, 45.0e-3, 50.0e-3, 55.0e-3,  60.0e-3,
      65.0e-3, 70.0e-3, 80.0e-3, 90.0e-3, 100.0e-3
   };
   static const Double_t xsunit = 1e-42*units::cm2;
   static const Double_t Onuebar[npts_nuebar] = {
      0,
      2.53e-2*xsunit, 7.27e-2*xsunit, 1.81e-1*xsunit, 4.21e-1*xsunit, 8.90e-1*xsunit,
      1.69*xsunit,    2.94*xsunit,    4.76*xsunit,    7.26*xsunit,    1.06e1*xsunit,
      1.48e1*xsunit,  2.64e1*xsunit,  4.29e1*xsunit,  6.46e1*xsunit,  9.17e1*xsunit,
      1.25e2*xsunit,  1.63e2*xsunit,  2.57e2*xsunit,  3.77e2*xsunit,  5.18e2*xsunit
   };
   // make spline via dummy TGraph because TSpline3's ctor isn't const correct
   const TGraph dummy(npts_nuebar,Evnuebar,Onuebar);
   fXsplNuebar = new TSpline3("16O_nu_e_bar_xsec",&dummy);
   fXsplNuebar->SetNpx(500);   
}
//____________________________________________________________________________
void KLVOxygenIBDPXSec::MakeNuESpline()
{
   // make a new xsec spline from the KLV paper's calculation
   // for the 16O + nu_e reaction
   // remove any spline that might already exist

   delete fXsplNue;

   static const Int_t npts_nue = 20;
   static const Double_t Evnue[npts_nue] = {
      15.0e-3, 17.5e-3, 20.0e-3, 22.5e-3,  25.0e-3,
      27.5e-3, 30.0e-3, 32.5e-3, 35.0e-3,  37.5e-3,
      40.0e-3, 45.0e-3, 50.0e-3, 55.0e-3,  60.0e-3,
      65.0e-3, 70.0e-3, 80.0e-3, 90.0e-3, 100.0e-3
   };
   static const Double_t xsunit = 1e-42*units::cm2;
   static const Double_t Onue[npts_nue] = {
      1.56e-6*xsunit, 8.42e-4*xsunit, 7.26e-3*xsunit, 3.99e-2*xsunit, 1.77e-1*xsunit,
      5.23e-1*xsunit, 1.25*xsunit,    2.58*xsunit,    4.76*xsunit,    8.05*xsunit,
      1.28e1*xsunit,  2.76e1*xsunit,  5.21e1*xsunit,  8.89e1*xsunit,  1.41e2*xsunit,
      2.12e2*xsunit,  3.02e2*xsunit,  5.52e2*xsunit,  8.92e2*xsunit,  1.32e3*xsunit
   };
   const TGraph dummy(npts_nue,Evnue,Onue);
   fXsplNue = new TSpline3("16O_nu_e_xsec",&dummy);
   fXsplNue->SetNpx(500);
}
//____________________________________________________________________________
double KLVOxygenIBDPXSec::Integral(const Interaction * interaction) const
{
   // Compute the total cross section for a free nucleon target

   assert(interaction!=0);
   if(! this -> ValidProcess    (interaction) ) return 0.;
   if(! this -> ValidKinematics (interaction) ) return 0.;

   const InitialState & init_state = interaction -> InitState();
   const double         Ev         = init_state.ProbeE(kRfHitNucRest);
   const int            prbpdg     = init_state.ProbePdg();
   
   double xsec = 0;
   
   if (pdg::IsNuE(prbpdg)) {
      assert(fXsplNue!=0);
      xsec = fXsplNue->Eval(Ev);
   } else if (pdg::IsAntiNuE(prbpdg)) {
      assert(fXsplNuebar!=0);
      xsec = fXsplNuebar->Eval(Ev);
   } else {
      LOG("KLVOxygen", pERROR) << "*** <Integral> Probe has invalid pdg ["
			       << init_state.ProbePdg() << "]";
   }
   
   return xsec;
}
//____________________________________________________________________________
bool KLVOxygenIBDPXSec::ValidProcess(const Interaction * interaction) const
{
   if(interaction->TestBit(kISkipProcessChk)) return true;
   
   // should be IBD and either nu_e + O16 or anu_e + O16
   if (interaction->ProcInfo().IsInverseBetaDecay()) {
      
      const InitialState & init_state = interaction -> InitState();
      if (init_state.TgtPdg() == kPdgTgtO16) {
	 
	 if ( (pdg::IsNuE(init_state.ProbePdg())) ||
	      (pdg::IsAntiNuE(init_state.ProbePdg())) ) {
	    
	    return true;
	    
	 } else {
	    LOG("KLVOxygen", pERROR) << "*** Probe has invalid pdg ["
				     << init_state.ProbePdg() << "]";
	 }
	 
      } else {
	 LOG("KLVOxygen", pERROR) << "*** Target has pdg ["
				  << init_state.TgtPdg()
				  << "], not 16O ("
				  << kPdgTgtO16 << ")!";
      }
      
   }
   
   return false;
}
//____________________________________________________________________________
bool KLVOxygenIBDPXSec::ValidKinematics(const Interaction* interaction) const
{
   // check energy range
   
   if(interaction->TestBit(kISkipKinematicChk)) return true;

   const InitialState & init_state = interaction -> InitState();
   const double Ev = init_state.ProbeE(kRfHitNucRest);
   if ( pdg::IsNuE(init_state.ProbePdg()) ) {
      if ( (Ev < kO16NuMinE) || (Ev > kMaxE) ) {
	 LOG("KLVOxygen", pERROR) << "*** Ev=" << Ev
				  << " outside range ("
				  << kO16NuMinE << ", " << kMaxE << ")!";
	 return false;
      }
   } else if ( pdg::IsAntiNuE(init_state.ProbePdg()) ) {
      if ( (Ev < kO16NubarThr) || (Ev > kMaxE) ) {
	 LOG("KLVOxygen", pERROR) << "*** Ev=" << Ev
				  << " outside range ("
				  << kO16NubarThr << ", " << kMaxE << ")!";
	 return false;
      }
   }
   
   const KPhaseSpace & kps = interaction->PhaseSpace();
   return kps.IsAboveThreshold();
}
//____________________________________________________________________________
void KLVOxygenIBDPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KLVOxygenIBDPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KLVOxygenIBDPXSec::LoadConfig(void)
{
   // make splines
   MakeAntiNuESpline();
   MakeNuESpline();
}
