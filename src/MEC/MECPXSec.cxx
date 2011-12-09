//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

         Steve Dytman <dytman+ \at pitt.edu>
         Pittsburgh University

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 28, 2010 - CA
   Integrated cross section for CCMEC is taken to be a fraction of the 
   CCQE cross section for the given neutrino energy and nucleus.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "MEC/MECPXSec.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"

using namespace genie;

//____________________________________________________________________________
MECPXSec::MECPXSec() :
XSecAlgorithmI("genie::MECPXSec")
{

}
//____________________________________________________________________________
MECPXSec::MECPXSec(string config) :
XSecAlgorithmI("genie::MECPXSec", config)
{

}
//____________________________________________________________________________
MECPXSec::~MECPXSec()
{

}
//____________________________________________________________________________
double MECPXSec::XSec(
		      const Interaction * interaction, KinePhaseSpace_t kps) const
{

// We have no clue what the meson exchange current contribution is.
// This is a toy model and is not used in default event generation.

//  if(! this -> ValidProcess    (interaction) ) return 0.;
//  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  double W  = kinematics.W();
  double Q2 = kinematics.Q2();
  LOG("MEC", pINFO) << "W, Q2 trial= " << W << "  " << Q2 ;

  //
  // Do a check whether W,Q2 is allowed. Return 0 otherwise.
  // 
  double Ev = interaction->InitState().ProbeE(kRfHitNucRest);  // kRfLab
  int nucleon_cluster_pdg = interaction->InitState().Tgt().HitNucPdg();
  double M2n = PDGLibrary::Instance()->Find(nucleon_cluster_pdg)-> Mass(); // nucleon cluster mass  
  //  double M2n = 1.88;
  double ml  = interaction->FSPrimLepton()->Mass();
  Range1D_t Wlim = genie::utils::kinematics::InelWLim(Ev, M2n, ml);
  LOG("MEC", pINFO) << "Ev, ml, M2n = " << Ev << "  " << ml << "  " << M2n;
  //  LOG("MEC", pINFO) << "Wlim= " << Wlim.min << "  " <<Wlim.max ;
  if(W < Wlim.min || W > Wlim.max)
    {double xsec = 0.;
      return xsec;
    } 
  Range1D_t Q2lim = genie::utils::kinematics::InelQ2Lim_W (Ev, M2n, ml, W, 0.);
  //  LOG("MEC", pINFO) << "Q2lim= " << Q2lim.min << "  " <<Q2lim.max ;
  if(Q2 < Q2lim.min || Q2 > Q2lim.max)
    {double xsec = 0.;
      return xsec;
    }

  //get x and y
  double x = 0.;
  double y = 0.;
  genie::utils::kinematics::WQ2toXY(Ev,M2n,W,Q2,x,y);
  LOG("MEC", pINFO) << "x = " << x << ", y = " << y;
  double Tmu = (1.-y)*Ev;

  // Calculate d^2xsec/dWdQ2
  double Wdep  = TMath::Gaus(W, fMass, fWidth);
  double Q2dep = TMath::Power(1+Q2/fMq2d, -2.);
  double nudep = TMath::Power(Tmu,2.5);
  //  LOG("MEC", pINFO) << "Tmu = " << Tmu << ", nudep = " << nudep;
  double xsec  = Wdep * Q2dep;// * nudep;

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
double MECPXSec::Integral(const Interaction * interaction) const
{
// Calculate the CCMEC cross section as a fraction of the CCQE cross section
// for the given nuclear target at the given energy.
// Alternative strategy is to calculate the MEC cross section as the difference
// of CCQE cross section for two different M_A values (eg ~1.3 GeV and ~1.0 GeV)
// Include hit-object combinatorial factor? Would yield different A-dependence
// for MEC and QE.
//

  bool   iscc   = interaction->ProcInfo().IsWeakCC();

  int    nupdg  = interaction->InitState().ProbePdg();
  int    tgtpdg = interaction->InitState().Tgt().Pdg();
  double E      = interaction->InitState().ProbeE(kRfLab);
  int nucleon_cluster_pdg = interaction->InitState().Tgt().HitNucPdg();

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

     // Use tunable fraction 
     // FFracCCQE is fraction of QE going to MEC
     // fFracCCQE_cluster is fraction of MEC going to each NN pair
     double fFracCCQE=0.;
     double Elo=1.;
     double Ehi=5.;

     if (E<Elo) fFracCCQE=fFracCCQElo;
     if (E>Elo&&E<Ehi) fFracCCQE=fFracCCQElo*(1.-(E-Elo)/(Ehi-Elo));
     if (E>Ehi) fFracCCQE=0.;

     double fFracCCQE_cluster=0.;
     if(pdg::IsNeutrino(nupdg) && nucleon_cluster_pdg==2000000200) fFracCCQE_cluster= .8;  //n+n
     if(pdg::IsNeutrino(nupdg) && nucleon_cluster_pdg==2000000201) fFracCCQE_cluster= .2;  //n+p
     if(pdg::IsAntiNeutrino(nupdg) && nucleon_cluster_pdg==2000000201) fFracCCQE_cluster= .8;   //n+p
     if(pdg::IsAntiNeutrino(nupdg) && nucleon_cluster_pdg==2000000202) fFracCCQE_cluster= .2;   //p+p


     xsec *= fFracCCQE*fFracCCQE_cluster;

     // Use gross combinatorial factor (number of 2-nucleon targets over number
     // of 1-nucleon targets) : (A-1)/2
     //     double combfact = (in->InitState().Tgt().A()-1)/2.;
     //     xsec *= combfact;

     delete in;
     return xsec;
  }

  return 0;
}
//____________________________________________________________________________
bool MECPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  if(!proc_info.IsMEC()) return false;

  return true;
}
//____________________________________________________________________________
void MECPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MECPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MECPXSec::LoadConfig(void)
{
  fXSecAlgCCQE = 0;

  fMq2d   = 0.8; // GeV
  fMass   = 1.9; // GeV
  fWidth  = 0.4; // GeV
  fFracCCQElo = 0.45; //fraction of CCQE xsec at Miniboone energies to CCMEC xsec
                      //  at first, this is energy independent

  // Get the specified CCQE cross section model
  fXSecAlgCCQE = 
     dynamic_cast<const XSecAlgorithmI *> (this->SubAlg("CCQEXSecModel"));
  assert(fXSecAlgCCQE);
}
//____________________________________________________________________________

