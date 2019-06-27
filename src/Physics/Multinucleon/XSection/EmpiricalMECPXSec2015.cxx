//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Steve Dytman <dytman+ \at pitt.edu>
         Pittsburgh University

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 28, 2011 - CA
   Integrated cross section for CCMEC is taken to be a fraction of the
   CCQE cross section for the given neutrino energy and nucleus.
 @ Dec 9, 2011 - SD
   Using a simple model now - 2N mass chosen with a Gaussian.
   Strength is tuned to get agreement with MiniBoone and NOMAD tot xs.
   Parameters tuned to get best agreement with MiniBoone muon angle/energy dist.
 @ May 14, 2012 - SD
   Changed empirical model to have correct Q2 behavior.  Shape is
   now Q2*Gaussian*dipole form factor.
   Add (e,e') MEC, start development to make them dependent on same model.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Units.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/XSection/EmpiricalMECPXSec2015.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
EmpiricalMECPXSec2015::EmpiricalMECPXSec2015() :
XSecAlgorithmI("genie::EmpiricalMECPXSec2015")
{

}
//____________________________________________________________________________
EmpiricalMECPXSec2015::EmpiricalMECPXSec2015(string config) :
XSecAlgorithmI("genie::EmpiricalMECPXSec2015", config)
{

}
//____________________________________________________________________________
EmpiricalMECPXSec2015::~EmpiricalMECPXSec2015()
{

}
//____________________________________________________________________________
double EmpiricalMECPXSec2015::XSec(
                      const Interaction * interaction, KinePhaseSpace_t kps) const
{

// meson exchange current contribution depends a lot on QE model.
// This is an empirical model in development, not used in default event generation.

//  if(! this -> ValidProcess    (interaction) ) return 0.;
//  if(! this -> ValidKinematics (interaction) ) return 0.;
  bool iscc = interaction->ProcInfo().IsWeakCC();
  bool isnc = interaction->ProcInfo().IsWeakNC();
  bool isem = interaction->ProcInfo().IsEM();

  const Kinematics &   kinematics = interaction -> Kine();
  double W  = kinematics.W();
  double Q2 = kinematics.Q2();
  //LOG("MEC", pINFO) << "W, Q2 trial= " << W << "  " << Q2 ;

  //
  // Do a check whether W,Q2 is allowed. Return 0 otherwise.
  //
  double Ev = interaction->InitState().ProbeE(kRfHitNucRest);  // kRfLab
  int nucleon_cluster_pdg = interaction->InitState().Tgt().HitNucPdg();
  double M2n = PDGLibrary::Instance()->Find(nucleon_cluster_pdg)-> Mass(); // nucleon cluster mass
  double M2n2 = M2n*M2n;
  double ml  = interaction->FSPrimLepton()->Mass();
  Range1D_t Wlim = isem ? genie::utils::kinematics::electromagnetic::InelWLim(Ev, ml, M2n) : genie::utils::kinematics::InelWLim(Ev, M2n, ml); 

  //LOG("MEC", pINFO) << "Ev, ml, M2n = " << Ev << "  " << ml << "  " << M2n;
  //LOG("MEC", pINFO) << "Wlim= " << Wlim.min << "  " <<Wlim.max ;
  if(W < Wlim.min || W > Wlim.max)
    {double xsec = 0.;
      return xsec;
    }
  //use proper Q2 limit from Controls.h
  Range1D_t Q2lim = isem ? genie::utils::kinematics::electromagnetic::InelQ2Lim_W(Ev, ml, M2n, W) : genie::utils::kinematics::InelQ2Lim_W (Ev, M2n, ml, W, kMinQ2Limit); 

  //LOG("MEC", pINFO) << "Q2lim= " << Q2lim.min << "  " <<Q2lim.max ;
  if(Q2 < Q2lim.min || Q2 > Q2lim.max)
    {double xsec = 0.;
      return xsec;
    }

  //get x and y
  double x = 0.;
  double y = 0.;
  genie::utils::kinematics::WQ2toXY(Ev,M2n,W,Q2,x,y);
  //  LOG("MEC", pINFO) << "x = " << x << ", y = " << y;
  // double Tmu = (1.-y)*Ev;  // UNUSED - comment to quiet compiler warnings

  // Calculate d^2xsec/dWdQ2 - first form factor which is common to both
  double Wdep  = TMath::Gaus(W, fMass, fWidth);
  double Q2dep = Q2*TMath::Power((1+Q2/fMq2d),-8.);
  //  double nudep = TMath::Power(Tmu,2.5);
  //  LOG("MEC", pINFO) << "Tmu = " << Tmu << ", nudep = " << nudep;
  double FF2  = Wdep * Q2dep;// * nudep;
  //LOG("MEC", pINFO) << "form factor = " << FF2 << ", Q2 = " << Q2 << ", W = " << W;

// using formulas in Bodek and Budd for (e,e') inclusive cross section
  double xsec = 1.;
  if(isem)  {
      // Calculate scattering angle
  //
  // Q^2 = 4 * E^2 * sin^2 (theta/2) / ( 1 + 2 * (E/M) * sin^2(theta/2) ) =>
  // sin^2 (theta/2) = MQ^2 / (4ME^2 - 2EQ^2)

    double E = Ev;
    double E2 = E*E;
    double sin2_halftheta = M2n*Q2 / (4*M2n*E2 - 2*E*Q2);
    //    double sin4_halftheta = TMath::Power(sin2_halftheta, 2.);
    double cos2_halftheta = 1.-sin2_halftheta;
    //    double cos_halftheta  = TMath::Sqrt(cos2_halftheta);
    double tan2_halftheta = sin2_halftheta/cos2_halftheta;
    double Q4 = Q2*Q2;

  // Calculate tau and the virtual photon polarization (epsilon)
  double tau     = Q2/(4*M2n2);
  //  double epsilon = 1. / (1. + 2.*(tau/x))*tan2_halftheta); //different than RosenbluthPXSec.cxx

  // Calculate the scattered lepton energy
  double Ep  = E / (1. + 2.*(E/M2n)*sin2_halftheta);
  double Ep2 = Ep*Ep;

  //calculate cross section - d2sig/dOmega dE for purely transverse process
  xsec = 4*kAem2*Ep2*cos2_halftheta/Q4 * FF2 * (tau/(1+tau) +2*tau*tan2_halftheta);
    }
  // use BB formula which seems to be same as Llewlyn-Smith
  // note B term is only difference between nu and antinu, so both same here
  else if(isnc||iscc){
    double tau     = Q2/(4*M2n2);
    double tau2 = tau*tau;
    double smufac = 4*M2n*Ev - Q2 - ml*ml;
    double A = (ml*ml+Q2)/M2n2 * (tau*(1+tau) - tau2*(1-tau)+4*tau2)/TMath::Power(1+tau,2.) * FF2;
    double C = tau/4/(1+tau) * FF2;
    xsec = A + smufac*smufac*C;   // CC or NC case - Llewelyn-Smith for transverse vector process.
  }
  // Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("MEC", pDEBUG)
     << "Jacobian for transformation to: "
                  << KinePhaseSpace::AsString(kps) << ", J = " << J;
#endif
    xsec *= J;
  }

  return xsec;
}
//____________________________________________________________________________
double EmpiricalMECPXSec2015::Integral(const Interaction * interaction) const
{
// Calculate the CCMEC cross section as a fraction of the CCQE cross section
// for the given nuclear target at the given energy.
// Alternative strategy is to calculate the MEC cross section as the difference
// of CCQE cross section for two different M_A values (eg ~1.3 GeV and ~1.0 GeV)
// Include hit-object combinatorial factor? Would yield different A-dependence
// for MEC and QE.
//

  bool iscc = interaction->ProcInfo().IsWeakCC();
  bool isnc = interaction->ProcInfo().IsWeakNC();
  bool isem = interaction->ProcInfo().IsEM();

  int    nupdg  = interaction->InitState().ProbePdg();
  int    tgtpdg = interaction->InitState().Tgt().Pdg();
  double E      = interaction->InitState().ProbeE(kRfLab);
  int nucleon_cluster_pdg = interaction->InitState().Tgt().HitNucPdg();
  double Z=interaction->InitState().Tgt().Z();
  double A=interaction->InitState().Tgt().A();
  double N=A-Z;

  if(iscc) {

     int nucpdg = 0;
     // neutrino CC: calculate the CCQE cross section resetting the
     // hit nucleon cluster to neutron
     if(pdg::IsNeutrino(nupdg)) {
         nucpdg = kPdgNeutron;
     }
     // anti-neutrino CC: calculate the CCQE cross section resetting the
     // hit nucleon cluster to proton
     else
     if(pdg::IsAntiNeutrino(nupdg)) {
         nucpdg = kPdgProton;
     }
     else {
         exit(1);
     }

     // Create a tmp QE process
     Interaction * in = Interaction::QELCC(tgtpdg,nucpdg,nupdg,E);

     // Calculate cross section for the QE process
     double xsec = fXSecAlgCCQE->Integral(in);

    // Add A dependence which is not known from theory
     double fFracADep = 1.;
     if(A>=12) fFracADep = TMath::Power((N/6.),fMECAPower-1.);

     // Use tunable fraction
     // FFracCCQE is fraction of QE going to MEC
     // fFracCCQE_cluster is fraction of MEC going to each NN pair
     //     double fFracCCQE = fFracCCQElo;

     double fFracCCQE_cluster=0.;
     if(pdg::IsNeutrino(nupdg) && nucleon_cluster_pdg==2000000201) fFracCCQE_cluster= fFracPN_CC;  //n+p
     if(pdg::IsNeutrino(nupdg) && nucleon_cluster_pdg==2000000200) fFracCCQE_cluster= 1.0-fFracPN_CC;  //n+n
     if(pdg::IsAntiNeutrino(nupdg) && nucleon_cluster_pdg==2000000201) fFracCCQE_cluster= fFracPN_CC;   //n+p
     if(pdg::IsAntiNeutrino(nupdg) && nucleon_cluster_pdg==2000000202) fFracCCQE_cluster= 1.0-fFracPN_CC;   //p+p


     xsec *= fFracCCQE*fFracCCQE_cluster*fFracADep;

     // Use gross combinatorial factor (number of 2-nucleon targets over number
     // of 1-nucleon targets) : (A-1)/2
     //     double combfact = (in->InitState().Tgt().A()-1)/2.;
     //     xsec *= combfact;

     delete in;
     return xsec;
  }

  else if(isnc) {
    int nucpdg = kPdgProton;
    // Create a tmp QE process
    Interaction * inp = Interaction::QELNC(tgtpdg,nucpdg,nupdg,E);
    nucpdg = kPdgNeutron;
    // Create a tmp QE process
    Interaction * inn = Interaction::QELNC(tgtpdg,nucpdg,nupdg,E);

    // Calculate cross section for the QE process - avg of p and n - best for isoscalar nuclei
    double xsec = (Z*fXSecAlgNCQE->Integral(inp) + N*fXSecAlgNCQE->Integral(inn))/A;

    // Add A dependence which is not known from theory
    double fFracADep = 1.;
    if(A>=12) fFracADep = TMath::Power((A/12.),fMECAPower-1.);

    // Use tunable fraction
    // FFracNCQE is fraction of QE going to MEC
    // fFracNCQE_cluster is fraction of MEC going to each NN pair
    double fFracNCQE_cluster=0.;
    if(nucleon_cluster_pdg==2000000200) fFracNCQE_cluster= 0.5*(1-fFracPN_NC);  //n+n
    if(nucleon_cluster_pdg==2000000201) fFracNCQE_cluster= fFracPN_NC;  //n+p
    if(nucleon_cluster_pdg==2000000202) fFracNCQE_cluster= 0.5*(1-fFracPN_NC);  //p+p
    xsec *= fFracNCQE*fFracNCQE_cluster*fFracADep;
    delete inn;
    delete inp;
    return xsec;
  }

  else if(isem) {
    int nucpdg = kPdgProton;
    // Create a tmp QE process
    Interaction * inp = Interaction::QELEM(tgtpdg,nucpdg,nupdg,E);
    nucpdg = kPdgNeutron;
    // Create a tmp QE process
    Interaction * inn = Interaction::QELEM(tgtpdg,nucpdg,nupdg,E);

    // Calculate cross section for the QE process - avg of p and n - best for isoscalar nuclei
    double xsec = (Z*fXSecAlgEMQE->Integral(inp) + N*fXSecAlgEMQE->Integral(inn))/A;

     // Add A dependence which is not known from theory, data wants high A suppression
    double fFracADep = 1.;
    if(A>=12) fFracADep = TMath::Power((A/12.),fMECAPower-1.);

    // Use tunable fraction
    // FFracEMQE is fraction of QE going to MEC
    // fFracEMQE_cluster is fraction of MEC going to each NN pair
    double fFracEMQE_cluster=0.;
    if(nucleon_cluster_pdg==2000000200) fFracEMQE_cluster= 0.5*(1-fFracPN_EM);  //n+n
    if(nucleon_cluster_pdg==2000000201) fFracEMQE_cluster= fFracPN_EM;  //n+p
    if(nucleon_cluster_pdg==2000000202) fFracEMQE_cluster= 0.5*(1-fFracPN_EM);  //p+p
    xsec *= fFracEMQE*fFracEMQE_cluster*fFracADep;
    delete inn;
    delete inp;
    return xsec;
  }

  return 0;
}
//____________________________________________________________________________
bool EmpiricalMECPXSec2015::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  if(!proc_info.IsMEC()) return false;

  return true;
}
//____________________________________________________________________________
void EmpiricalMECPXSec2015::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void EmpiricalMECPXSec2015::Configure(string config)
{
  Algorithm::Configure(config);

  Registry* algos = AlgConfigPool::Instance() -> GlobalParameterList() ;
  string global_key_head = "XSecModel@genie::EventGenerator/" ;
  string local_key_head = "XSecModel-" ;

  Registry r( "EmpiricalMECPXSec2015_specific", false ) ;
  r.Set( local_key_head + "QEL-NC", algos -> GetAlg( global_key_head + "QEL-NC") ) ;
  r.Set( local_key_head + "QEL-CC", algos -> GetAlg( global_key_head + "QEL-CC") ) ;
  r.Set( local_key_head + "QEL-EM", algos -> GetAlg( global_key_head + "QEL-EM") ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void EmpiricalMECPXSec2015::LoadConfig(void)
{
  fXSecAlgCCQE = 0;
  fXSecAlgNCQE = 0;
  fXSecAlgEMQE = 0;

  GetParam( "EmpiricalMEC-Mq2d", fMq2d ) ;
  GetParam( "EmpiricalMEC-Mass", fMass ) ;
  GetParam( "EmpiricalMEC-Width", fWidth ) ;
  GetParam( "EmpiricalMEC-APower", fMECAPower ) ;

  GetParam( "EmpiricalMEC-FracPN_NC", fFracPN_NC ) ;
  GetParam( "EmpiricalMEC-FracPN_CC", fFracPN_CC ) ;
  GetParam( "EmpiricalMEC-FracPN_EM", fFracPN_EM ) ;

  GetParam( "EmpiricalMEC-FracCCQE", fFracCCQE ) ;
  GetParam( "EmpiricalMEC-FracNCQE", fFracNCQE ) ;
  GetParam( "EmpiricalMEC-FracEMQE", fFracEMQE ) ;

  string key_head = "XSecModel-" ;

  fXSecAlgNCQE =
     dynamic_cast<const XSecAlgorithmI *> ( this -> SubAlg( key_head + "QEL-NC" ) ) ;
  assert(fXSecAlgNCQE);

  fXSecAlgCCQE =
     dynamic_cast<const XSecAlgorithmI *> ( this -> SubAlg( key_head + "QEL-CC" ) ) ;
  assert(fXSecAlgCCQE);

  fXSecAlgEMQE =
     dynamic_cast<const XSecAlgorithmI *> ( this -> SubAlg( key_head + "QEL-EM" ) ) ;
  assert(fXSecAlgEMQE);
}
//____________________________________________________________________________
