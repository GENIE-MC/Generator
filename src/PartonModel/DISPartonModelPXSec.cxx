//____________________________________________________________________________
/*!

\class    genie::DISPartonModelPXSec

\brief    DIS differential (d^2xsec/dxdy) cross section

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.033002

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Base/DISStructureFuncModelI.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISPartonModelPXSec.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISPartonModelPXSec::DISPartonModelPXSec() :
XSecAlgorithmI("genie::DISPartonModelPXSec")
{

}
//____________________________________________________________________________
DISPartonModelPXSec::DISPartonModelPXSec(string config) :
XSecAlgorithmI("genie::DISPartonModelPXSec", config)
{

}
//____________________________________________________________________________
DISPartonModelPXSec::~DISPartonModelPXSec()
{

}
//____________________________________________________________________________
double DISPartonModelPXSec::XSec(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- Get kinematical & init-state parameters
  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  double E     = init_state.GetProbeE(kRfStruckNucAtRest);
  double ml    = interaction->GetFSPrimaryLepton()->Mass();
  double Mnuc  = init_state.GetTarget().StruckNucleonMass();
  double Mnuc2 = TMath::Power(Mnuc, 2);
  double x     = kinematics.x();
  double y     = kinematics.y();

  //----- One of the xsec terms changes sign for antineutrinos
  int sign = 1;
  if( pdg::IsAntiNeutrino(init_state.GetProbePDGCode()) ) sign = -1;

  //----- Calculate the DIS structure functions
  fDISSF.Calculate(interaction); 

  LOG("DISXSec", pDEBUG)  << "\n" << fDISSF;

  //-- calculate auxiliary parameters
  double ml2  = ml    * ml;
  double ml4  = ml2   * ml2;
  double E2   = E     * E;
  double Gfac = (kGF*kGF*Mnuc*E) / kPi;

  //----- Build all dxsec/dxdy terms
  double term1 = y * ( x*y + ml2/(2*E*Mnuc) );
  double term2 = 1 - y - Mnuc*x*y/(2*E) - ml2/(4*E2);
  double term3 = sign*x*y*(1-y/2) - y*ml2/(4*Mnuc*E);
  double term4 = x*y*ml2/(2*Mnuc*E) + ml4/(4*Mnuc2*E2);
  double term5 = -1.*ml2/(2*Mnuc*E);

  LOG("DISXSec", pDEBUG)  
    << "\nd^2xsec/dxdy ~ (" << term1 << ")*F1 + (" << term2 
    << ")*F2 +(" << term3 << ")*F3 + (" << term4 << ")*F4 + ("
    << term5 << ")*F5";

  //----- Compute the differential cross section
  term1 *= fDISSF.F1();
  term2 *= fDISSF.F2();
  term3 *= fDISSF.F3();
  term4 *= fDISSF.F4();
  term5 *= fDISSF.F5();

  double xsec = Gfac*(term1 + term2 + term3 + term4 + term5);
  xsec = TMath::Max(xsec,0.);

  LOG("DISXSec", pDEBUG)
      << "d^2xsec/dxdy (E = " << E << ", x = " << x << ", y = " << y << ") = "
      << xsec;

  return xsec;
}
//____________________________________________________________________________
bool DISPartonModelPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  return true;
}
//____________________________________________________________________________
bool DISPartonModelPXSec::ValidKinematics(
                                       const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  double E     = init_state.GetProbeE(kRfStruckNucAtRest);
  double Mnuc  = init_state.GetTarget().StruckNucleonMass();
  double Mnuc2 = TMath::Power(Mnuc, 2);
  double x     = kinematics.x();
  double y     = kinematics.y();
  double W2    = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W     = TMath::Sqrt(W2);
  double Q2    = 2*Mnuc*E*x*y;

  //----- Get the physical W and Q2 range and check whether the current W,Q2
  //      pair is allowed
  Range1D_t rW  = utils::kinematics::KineRange (interaction, kKVW);
  Range1D_t rQ2 = utils::kinematics::KineRange (interaction, kKVQ2);

  bool in_range = utils::math::IsWithinLimits(Q2, rQ2)
                                      && utils::math::IsWithinLimits(W, rW);
  if(!in_range) {
       LOG("DISXSec", pDEBUG)
             << "\n *** point (W = " << W
                           << ", Q2 = " << Q2 << " is not in physical range";
       LOG("DISXSec", pDEBUG)
             << "\n Physical W range: "
                               << "[" << rW.min << ", " << rW.max << "] GeV";
       LOG("DISXSec", pDEBUG)
             << "\n Physical Q2 range: "
                           << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
       return false;
   }

  return true;
}
//____________________________________________________________________________
void DISPartonModelPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISPartonModelPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISPartonModelPXSec::LoadConfig(void)
{
  fDISSFModel = 0;
  fDISSFModel = dynamic_cast<const DISStructureFuncModelI *> (
                                this->SubAlg("sf-alg-name", "sf-param-set"));
  assert(fDISSFModel);

  fDISSF.SetModel(fDISSFModel); // <-- attach algorithm
}
//____________________________________________________________________________

