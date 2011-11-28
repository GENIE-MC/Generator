//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

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
#include "Messenger/Messenger.h"
#include "MEC/MECPXSec.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
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

  double Wdep  = TMath::Gaus(W, fMass, fWidth);
  double Q2dep = TMath::Power(1+Q2/fMq2d, -1.5);
  double xsec  = Wdep * Q2dep;

  // The algorithm computes d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("QMEC", pDEBUG)
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

  Interaction in(*interaction);

  bool iscc  = in.ProcInfo().IsWeakCC();
  int  nupdg = in.InitState().ProbePdg();

  if(iscc) {

     // gross combinatorial factor (number of 2-nucleon targets over number
     // of 1-nucleon targets) : (A-1)/2
     double combfact = (in.InitState().Tgt().A()-1)/2.;

     // neutrino CC: calculate the CCQE cross section resetting the
     // hit nucleon cluster to neutron
     double xsec = 0;
     if(pdg::IsNeutrino(nupdg)) {
         in.InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgNeutron);
         xsec = fXSecAlgCCQE->Integral(&in);
         xsec *= fFracCCQE;
     }
     // anti-neutrino CC: calculate the CCQE cross section resetting the
     // hit nucleon cluster to proton
     else
     if(pdg::IsAntiNeutrino(nupdg)) {
         in.InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgProton);
         xsec = fXSecAlgCCQE->Integral(&in);
     }

     xsec *= (combfact*fFracCCQE);
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

  fMq2d   = 0.5; // GeV
  fMass   = 2.1; // GeV
  fWidth  = 0.3; // GeV
  fEc     = 0.4; // GeV
  fFracCCQE = 0.1;

  // Get the specified CCQE cross section model
  fXSecAlgCCQE = 
     dynamic_cast<const XSecAlgorithmI *> (this->SubAlg("CCQEXSecModel"));
  assert(fXSecAlgCCQE);
}
//____________________________________________________________________________

