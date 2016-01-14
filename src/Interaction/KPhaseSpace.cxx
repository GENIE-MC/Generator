//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 21, 2007 - CA
   Handle the introduction of a new type of coherent interactions (coherent 
   elastic). The kinematical limits for the old 'coherent' type now apply
   only to 'coherent pi production'. Computing kinematical limits for
   coherent elastic is not included at this moment.
 @ Jan 18, 2008 - CA
   Add protection against non-positive energy thresholds
 @ Feb 25, 2008 - CA
   Added threshold for anomaly-mediated single gamma interactions
 @ Mar 03, 2009 - CA
   Add threshold and kinematical limits for diffractive scattering.
   Adapt to naming changes made to the coherent generator for including 
   coherent vector meson production.
 @ Sep 19, 2009 - CR
   Add threshold and kinematical limits for inverse beta decay.
 @ Nov 28, 2011 - CA
   Add threshold for MEC.
 @ Feb 12, 2013 - CA (code from Rosen Matev)
  Handle the IMD annihilation channel. Add mass_electron^2 term in calculation
  of kinematic limits for neutrino-electron elastic scattering.
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Interaction/KPhaseSpace.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

ClassImp(KPhaseSpace)

//____________________________________________________________________________
KPhaseSpace::KPhaseSpace(void) :
TObject()
{
  this->UseInteraction(0);
}
//___________________________________________________________________________
KPhaseSpace::KPhaseSpace(const Interaction * in) :
TObject()
{
  this->UseInteraction(in);
}
//___________________________________________________________________________
KPhaseSpace::~KPhaseSpace(void)
{

}
//___________________________________________________________________________
void KPhaseSpace::UseInteraction(const Interaction * in) 
{
  fInteraction = in;
}
//___________________________________________________________________________
double KPhaseSpace::Threshold(void) const
{
  const ProcessInfo &  pi         = fInteraction->ProcInfo();
  const InitialState & init_state = fInteraction->InitState();
  const XclsTag &      xcls       = fInteraction->ExclTag();
  const Target &       tgt        = init_state.Tgt();

  double ml = fInteraction->FSPrimLepton()->Mass();

  if (pi.IsSingleKaon()) {
    int kaon_pdgc = xcls.StrangeHadronPdg();
    double Mi   = tgt.HitNucP4Ptr()->M(); // initial nucleon mass
    // Final nucleon can be different for K0 interaction
    double Mf = (xcls.NProtons()==1) ? kProtonMass : kNeutronMass;  
    double mk   = PDGLibrary::Instance()->Find(kaon_pdgc)->Mass();
    double ml   = PDGLibrary::Instance()->Find(fInteraction->FSPrimLeptonPdg())->Mass();
    double mtot = ml + mk + Mf; // total mass of FS particles
    double Ethresh = (mtot*mtot - Mi*Mi)/(2. * Mf);
    return Ethresh;
  }

  if (pi.IsCoherent()) {
    int tgtpdgc = tgt.Pdg(); // nuclear target PDG code (10LZZZAAAI)
    double MA   = PDGLibrary::Instance()->Find(tgtpdgc)->Mass();
    double m    = ml + kPionMass;
    double m2   = TMath::Power(m,2);
    double Ethr = m + 0.5*m2/MA;
    return TMath::Max(0.,Ethr);
  }

  if(pi.IsQuasiElastic()     || 
     pi.IsInverseBetaDecay() ||
     pi.IsResonant()         || 
     pi.IsDeepInelastic()    || 
     pi.IsDiffractive()) 
  {
    assert(tgt.HitNucIsSet());
    double Mn   = tgt.HitNucP4Ptr()->M();
    double Mn2  = TMath::Power(Mn,2);
    double Wmin = (pi.IsQuasiElastic() || pi.IsInverseBetaDecay()) ? 
                  kNucleonMass : kNucleonMass+kPionMass;

    if(xcls.IsCharmEvent()) {
       if(xcls.IsInclusiveCharm()) {
          Wmin = kNucleonMass+kLightestChmHad;
       } else {
          int cpdg = xcls.CharmHadronPdg();
          double mchm = PDGLibrary::Instance()->Find(cpdg)->Mass();
          if(pi.IsQuasiElastic() || pi.IsInverseBetaDecay()) { 
            Wmin = mchm + controls::kASmallNum; 
          } 
          else { 
            Wmin = kNeutronMass + mchm + controls::kASmallNum; 
          }
       }//incl.?
    }//charm?

    double smin = TMath::Power(Wmin+ml,2.);
    double Ethr = 0.5*(smin-Mn2)/Mn;

    return TMath::Max(0.,Ethr);
  }

  if(pi.IsInverseMuDecay() || pi.IsIMDAnnihilation()) {
    double Ethr = 0.5 * (kMuonMass2-kElectronMass2)/kElectronMass;
    return TMath::Max(0.,Ethr);
  }

  if(pi.IsNuElectronElastic()) {
    return 0;
  }
  if(pi.IsAMNuGamma()) {
    return 0;
  }
  if(pi.IsMEC()) {
    assert(tgt.HitNucIsSet());
    double Mn   = tgt.HitNucP4Ptr()->M();
    double Mn2  = TMath::Power(Mn,2);
    double Wmin = fInteraction->RecoilNucleon()->Mass(); // mass of the recoil nucleon cluster 
    double smin = TMath::Power(Wmin+ml,2.);
    double Ethr = 0.5*(smin-Mn2)/Mn;
    return TMath::Max(0.,Ethr);
  }

  SLOG("KPhaseSpace", pERROR) 
         << "Can't compute threshold for \n" << *fInteraction;
  exit(1);

  return 99999999;
}
//___________________________________________________________________________
Range1D_t KPhaseSpace::Limits(KineVar_t kvar) const
{
// Compute limits for the input kinematic variable irrespective of any other
// relevant kinematical variable
//
  assert(fInteraction);

  switch(kvar) {
  case(kKVW)  : return this->WLim();  break;
  case(kKVQ2) : return this->Q2Lim(); break;
  case(kKVq2) : return this->q2Lim(); break;
  case(kKVx)  : return this->XLim();  break;
  case(kKVy)  : return this->YLim();  break;
  default:
   LOG("KPhaseSpace", pERROR) 
        << "Couldn't compute limits for " << KineVar::AsString(kvar);
   Range1D_t R(-1.,-1);
   return R;
  }
}
//____________________________________________________________________________
double KPhaseSpace::Minimum(KineVar_t kvar) const
{
  Range1D_t lim = this->Limits(kvar);
  return lim.min;
}
//___________________________________________________________________________
double KPhaseSpace::Maximum(KineVar_t kvar) const
{
  Range1D_t lim = this->Limits(kvar);
  return lim.max;
}
//___________________________________________________________________________
bool KPhaseSpace::IsAboveThreshold(void) const
{
  double E    = 0.;
  double Ethr = this->Threshold();

  const ProcessInfo &  pi         = fInteraction->ProcInfo();
  const InitialState & init_state = fInteraction->InitState();

  if (pi.IsCoherent()       || 
      pi.IsInverseMuDecay() || 
      pi.IsIMDAnnihilation() || 
      pi.IsNuElectronElastic()) 
  {
      E = init_state.ProbeE(kRfLab);
  }

  if(pi.IsQuasiElastic()     || 
     pi.IsInverseBetaDecay() ||
     pi.IsResonant()         || 
     pi.IsDeepInelastic()    || 
     pi.IsDiffractive()      || 
     pi.IsSingleKaon()       ||
     pi.IsAMNuGamma()) 
  {
      E = init_state.ProbeE(kRfHitNucRest);
  }

  LOG("KPhaseSpace", pDEBUG) << "E = " << E << ", Ethr = " << Ethr; 
  return (E>Ethr);
}
//___________________________________________________________________________
bool KPhaseSpace::IsAllowed(void) const
{
  const ProcessInfo & pi   = fInteraction->ProcInfo();
  const Kinematics &  kine = fInteraction->Kine();

  // ASK single kaon:
  // XSec code returns zero when kinematics are not allowed
  // Here just let kinematics always be allowed
  if(pi.IsSingleKaon()) {
    return true;
  }

  // QEL: 
  //  Check the running Q2 vs the Q2 limits
  if(pi.IsQuasiElastic() || pi.IsInverseBetaDecay()) {
    Range1D_t Q2l = this->Q2Lim();
    double    Q2  = kine.Q2();
    bool in_phys = math::IsWithinLimits(Q2, Q2l);
    bool allowed = in_phys;
    return allowed;
  }

  // RES
  //   Check the running W vs the W limits 
  //   & the running Q2 vs Q2 limits for the given W
  if(pi.IsResonant()) {
    Range1D_t Wl  = this->WLim();
    Range1D_t Q2l = this->Q2Lim_W();
    double    W   = kine.W();
    double    Q2  = kine.Q2();
    bool in_phys = (math::IsWithinLimits(Q2, Q2l) && math::IsWithinLimits(W, Wl));
    bool allowed = in_phys;
    return allowed;
  }

  // DIS
  if(pi.IsDeepInelastic()) {
    Range1D_t Wl  = this->WLim();
    Range1D_t Q2l = this->Q2Lim_W();
    double    W   = kine.W();
    double    Q2  = kine.Q2();
    bool in_phys = (math::IsWithinLimits(Q2, Q2l) && math::IsWithinLimits(W, Wl));
    bool allowed = in_phys;
    return allowed;
  }

  //IMD
  if(pi.IsInverseMuDecay() || pi.IsIMDAnnihilation() || pi.IsNuElectronElastic()) {
    Range1D_t yl = this->YLim();
    double    y  = kine.y();
    bool in_phys = math::IsWithinLimits(y, yl);
    bool allowed = in_phys;
    return allowed;
  }

  //COH
  if (pi.IsCoherent()) {
    Range1D_t xl = this->XLim();
    Range1D_t yl = this->YLim();
    double    x  = kine.x();
    double    y  = kine.y();
    bool in_phys = (math::IsWithinLimits(x, xl) && math::IsWithinLimits(y, yl));
    bool allowed = in_phys;
    return allowed;
  }

  return false;
}
//___________________________________________________________________________
Range1D_t KPhaseSpace::WLim(void) const
{
// Computes hadronic invariant mass limits. 
// For QEL the range reduces to the recoil nucleon mass. 
// For DIS & RES the calculation proceeds as in kinematics::InelWLim(). 
// It is not computed for other interactions
//
  Range1D_t Wl;
  Wl.min = -1;
  Wl.max = -1;

  const ProcessInfo & pi = fInteraction->ProcInfo();
  bool is_qel  = pi.IsQuasiElastic()  || pi.IsInverseBetaDecay();
  bool is_inel = pi.IsDeepInelastic() || pi.IsResonant();

  if(is_qel) {
    double MR = fInteraction->RecoilNucleon()->Mass();
    Wl.min = MR;
    Wl.max = MR;
    return Wl;
  }
  if(is_inel) {
    const InitialState & init_state = fInteraction->InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double M  = init_state.Tgt().HitNucP4Ptr()->M(); //can be off m/shell
    double ml = fInteraction->FSPrimLepton()->Mass();
    Wl = kinematics::InelWLim(Ev,M,ml);  
    if(fInteraction->ExclTag().IsCharmEvent()) {
      //Wl.min = TMath::Max(Wl.min, kNeutronMass+kPionMass+kLightestChmHad);
      Wl.min = TMath::Max(Wl.min, kNeutronMass+kLightestChmHad);
      if(Wl.min>Wl.max) {Wl.min=-1; Wl.max=-1;}
    }
    return Wl;
  }
  return Wl;
}
//____________________________________________________________________________
Range1D_t KPhaseSpace::Q2Lim_W(void) const
{
// Computes momentum transfer (Q2>0) limits @ the input invariant mass
// The calculation proceeds as in kinematics::InelQ2Lim_W(). 
// For QEL, W is set to the recoil nucleon mass

  Range1D_t Q2l;
  Q2l.min = -1;
  Q2l.max = -1;

  const ProcessInfo & pi = fInteraction->ProcInfo();
  bool is_qel  = pi.IsQuasiElastic()  || pi.IsInverseBetaDecay();
  bool is_inel = pi.IsDeepInelastic() || pi.IsResonant();

  if(!is_qel && !is_inel) return Q2l;

  const InitialState & init_state = fInteraction->InitState();
  double Ev  = init_state.ProbeE(kRfHitNucRest);
  double M   = init_state.Tgt().HitNucP4Ptr()->M(); // can be off m/shell
  double ml  = fInteraction->FSPrimLepton()->Mass();

  double W = 0;
  if(is_qel) W = fInteraction->RecoilNucleon()->Mass();
  else       W = kinematics::W(fInteraction);

  if (pi.IsInverseBetaDecay()) {
     Q2l = kinematics::InelQ2Lim_W(Ev,M,ml,W, controls::kMinQ2Limit_VLE);
  } else {
     Q2l = kinematics::InelQ2Lim_W(Ev,M,ml,W);
  }
  return Q2l;
}
//____________________________________________________________________________
Range1D_t KPhaseSpace::q2Lim_W(void) const
{
// As Q2Lim_W(void) but with reversed sign (Q2 -> q2)
//
  Range1D_t Q2 = this->Q2Lim_W();
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t KPhaseSpace::Q2Lim(void) const
{
// Computes momentum transfer (Q2>0) limits irrespective of the invariant mass
// For QEL this is identical to Q2Lim_W (since W is fixed)
// For RES & DIS, the calculation proceeds as in kinematics::InelQ2Lim(). 
//
  Range1D_t Q2l;
  Q2l.min = -1;
  Q2l.max = -1;

  const ProcessInfo & pi = fInteraction->ProcInfo();
  bool is_qel  = pi.IsQuasiElastic()  || pi.IsInverseBetaDecay();
  bool is_inel = pi.IsDeepInelastic() || pi.IsResonant();

  if(!is_qel && !is_inel) return Q2l;

  const InitialState & init_state = fInteraction->InitState();
  double Ev  = init_state.ProbeE(kRfHitNucRest);
  double M   = init_state.Tgt().HitNucP4Ptr()->M(); // can be off m/shell
  double ml  = fInteraction->FSPrimLepton()->Mass();

  const XclsTag & xcls = fInteraction->ExclTag();

  // quasi-elastic
  if(is_qel) {
    double W = fInteraction->RecoilNucleon()->Mass();
    if(xcls.IsCharmEvent()) { 
      int charm_pdgc = xcls.CharmHadronPdg();           
      W = PDGLibrary::Instance()->Find(charm_pdgc)->Mass();
    }
    if (pi.IsInverseBetaDecay()) {
       Q2l = kinematics::InelQ2Lim_W(Ev,M,ml,W,controls::kMinQ2Limit_VLE);
    } else {
       Q2l = kinematics::InelQ2Lim_W(Ev,M,ml,W);
    }
    return Q2l;
  }

  // inelastic
  Q2l = kinematics::InelQ2Lim(Ev,M,ml);  
  return Q2l;
}
//____________________________________________________________________________
Range1D_t KPhaseSpace::q2Lim(void) const
{
// As Q2Lim(void) but with reversed sign (Q2 -> q2)
//
  Range1D_t Q2 = this->Q2Lim();
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t KPhaseSpace::XLim(void) const
{
// Computes x-limits;

  Range1D_t xl;
  xl.min = -1;
  xl.max = -1;

  const ProcessInfo & pi = fInteraction->ProcInfo();

  //RES+DIS
  bool is_inel = pi.IsDeepInelastic() || pi.IsResonant();
  if(is_inel) {
    const InitialState & init_state  = fInteraction->InitState();
    double Ev  = init_state.ProbeE(kRfHitNucRest);
    double M   = init_state.Tgt().HitNucP4Ptr()->M(); // can be off m/shell
    double ml  = fInteraction->FSPrimLepton()->Mass();
    xl = kinematics::InelXLim(Ev,M,ml);
    return xl;
  }
  //COH
  bool is_coh = pi.IsCoherent();
  if(is_coh) {
    xl = kinematics::CohXLim();
    return xl;
  }
  //QEL
  bool is_qel = pi.IsQuasiElastic() || pi.IsInverseBetaDecay();
  if(is_qel) {
    xl.min = 1;
    xl.max = 1;
    return xl;
  }
  bool is_dfr = pi.IsDiffractive();
  if(is_dfr) {
    xl.min =      controls::kASmallNum;
    xl.max = 1. - controls::kASmallNum;
    return xl;
  }

  return xl;
}
//____________________________________________________________________________
Range1D_t KPhaseSpace::YLim(void) const
{
  Range1D_t yl;
  yl.min = -1;
  yl.max = -1;

  const ProcessInfo & pi = fInteraction->ProcInfo();

  //RES+DIS
  bool is_inel = pi.IsDeepInelastic() || pi.IsResonant();
  if(is_inel) {
    const InitialState & init_state = fInteraction->InitState();
    double Ev  = init_state.ProbeE(kRfHitNucRest);
    double M   = init_state.Tgt().HitNucP4Ptr()->M(); // can be off m/shell
    double ml  = fInteraction->FSPrimLepton()->Mass();
    yl = kinematics::InelYLim(Ev,M,ml);
    return yl;
  }
  //COH
  bool is_coh = pi.IsCoherent();
  if(is_coh) {  
    const InitialState & init_state = fInteraction->InitState();
    double EvL = init_state.ProbeE(kRfLab);
    double ml  = fInteraction->FSPrimLepton()->Mass();
    yl = kinematics::CohYLim(EvL,ml);
    return yl;
  }
  // IMD
  if(pi.IsInverseMuDecay() || pi.IsIMDAnnihilation() || pi.IsNuElectronElastic()) {
    const InitialState & init_state = fInteraction->InitState();
    double Ev = init_state.ProbeE(kRfLab);
    double ml = fInteraction->FSPrimLepton()->Mass();
    double me = kElectronMass;
    yl.min = controls::kASmallNum;
    yl.max = 1 - (ml*ml + me*me)/(2*me*Ev) - controls::kASmallNum;
    return yl;
  }
  bool is_dfr = pi.IsDiffractive();
  if(is_dfr) {
    const InitialState & init_state = fInteraction -> InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest); 
    double ml = fInteraction->FSPrimLepton()->Mass();
    yl.min = kPionMass/Ev + controls::kASmallNum;
    yl.max = 1. -ml/Ev - controls::kASmallNum;
    return yl;
  }
  return yl;
}
//____________________________________________________________________________
Range1D_t KPhaseSpace::YLim_X(void) const
{
// Computes kinematical limits for y @ the input x

  Range1D_t yl;
  yl.min = -1;
  yl.max = -1;

  const ProcessInfo & pi = fInteraction->ProcInfo();

  //RES+DIS
  bool is_inel = pi.IsDeepInelastic() || pi.IsResonant();
  if(is_inel) {
    const InitialState & init_state = fInteraction->InitState();
    double Ev  = init_state.ProbeE(kRfHitNucRest);
    double M   = init_state.Tgt().HitNucP4Ptr()->M(); // can be off m/shell
    double ml  = fInteraction->FSPrimLepton()->Mass();
    double x   = fInteraction->Kine().x();
    yl = kinematics::InelYLim_X(Ev,M,ml,x);
    return yl;
  }
  //COH
  bool is_coh = pi.IsCoherent();
  if(is_coh) {  
    const InitialState & init_state = fInteraction->InitState();
    double EvL = init_state.ProbeE(kRfLab);
    double ml  = fInteraction->FSPrimLepton()->Mass();
    yl = kinematics::CohYLim(EvL,ml);
    return yl;
  }
  return yl;
}
//____________________________________________________________________________
