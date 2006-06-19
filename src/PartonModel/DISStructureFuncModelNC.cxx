//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModelNC.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISStructureFuncModelNC::DISStructureFuncModelNC() :
DISStructureFuncModel("genie::DISStructureFuncModelNC")
{

}
//____________________________________________________________________________
DISStructureFuncModelNC::DISStructureFuncModelNC(string config):
DISStructureFuncModel("genie::DISStructureFuncModelNC", config)
{

}
//____________________________________________________________________________
DISStructureFuncModelNC::~DISStructureFuncModelNC()
{

}
//____________________________________________________________________________
void DISStructureFuncModelNC::Calculate(const Interaction * interaction) const
{
  // Reset mutable members
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;

  const Kinematics & kine  = interaction->GetKinematics();
  double x = kine.x();
  if(x<=0. || x>1) {
     LOG("DISSF", pERROR)
                 << "scaling variable x = " << x << "! Can not compute SFs";
     return;
  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Here all corrections to computing the slow rescaling variable and the
  // K factors are applied

  this->CalcPDFs(interaction);

  double u    = fPDF->UpValence() + fPDF->UpSea();
  double ubar = fPDF->UpSea();
  double d    = fPDF->DownValence() + fPDF->DownSea();
  double dbar = fPDF->DownSea();
  double s    = fPDF->Strange();
  double sbar = fPDF->Strange();
  double c    = fPDF->Charm();
  double cbar = fPDF->Charm();

  const InitialState & init_state = interaction->GetInitialState();

  bool isP = pdg::IsProton ( init_state.GetTarget().StruckNucleonPDGCode() );
  bool isN = pdg::IsNeutron( init_state.GetTarget().StruckNucleonPDGCode() );

  bool isNu    = pdg::IsNeutrino    ( init_state.GetProbePDGCode() );
  bool isNuBar = pdg::IsAntiNeutrino( init_state.GetProbePDGCode() );

  if(!isNu && !isNuBar) {
     LOG("DISSF", pWARN) << "v type is not handled" << *interaction;
     return;
  }

  double GL   = (isNu) ? ( 0.5 - (2./3.)*fSin2thw ) : (     - (2./3.)*fSin2thw);
  double GR   = (isNu) ? (     - (2./3.)*fSin2thw ) : ( 0.5 - (2./3.)*fSin2thw);
  double GLp  = (isNu) ? (-0.5 + (1./3.)*fSin2thw ) : (       (1./3.)*fSin2thw);
  double GRp  = (isNu) ? (       (1./3.)*fSin2thw ) : (-0.5 + (1./3.)*fSin2thw);

  double GL2  =  TMath::Power(GL,  2.);
  double GR2  =  TMath::Power(GR,  2.);
  double GLp2 =  TMath::Power(GLp, 2.);
  double GRp2 =  TMath::Power(GRp, 2.);

  double F2  = 0.;
  double xF3 = 0.;

  if(isP) {
      F2  = 2*( (GL2 + GR2) * (u + c + ubar + cbar) +
                            (GLp2 + GRp2) * (d + s + dbar + sbar) );
      xF3 = 2*( (GL2 - GR2) * (u + c - ubar - cbar) +
                            (GLp2 - GRp2) * (d + s - dbar - sbar) );

  } else if (isN) {
      F2  = 2*( (GL2 + GR2) * (d + c + dbar + cbar) +
                            (GLp2 + GRp2) * (u + s + ubar + sbar) );
      xF3 = 2*( (GL2 - GR2) * (d + c - dbar - cbar) +
                            (GLp2 - GRp2) * (u + s - ubar - sbar) );
  } else {
     LOG("DISSF", pWARN) << "N type is not handled" << *interaction;
     return;
  }

  // compute nuclear modification factor 
  double f  = this->NuclMod(interaction);

  LOG("DISSF", pDEBUG) << "Nucl. Factor = " << f;

  // compute the structure functions F1-F6
  fF3 = f * xF3/x;
  fF2 = f * F2;
  fF1 = 0.5 * fF2/x; 
}
//____________________________________________________________________________
void DISStructureFuncModelNC::LoadConfig(void) 
{
  DISStructureFuncModel::LoadConfig();

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  double thw = fConfig->GetDoubleDef(
                          "weinberg-angle", gc->GetDouble("WeinbergAngle"));
  fSin2thw = TMath::Power(TMath::Sin(thw), 2);
}
//____________________________________________________________________________

