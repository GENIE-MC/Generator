//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TLorentzVector.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Conventions/Units.h"

using namespace genie;
using namespace genie::constants;

//______________________________________________________________________
double genie::utils::mec::GetTmuCostFromq0q3(
    double dq0, double dq3, double Enu, double lmass,
    double &tmu, double &cost, double &area)
{
 tmu = Enu - dq0 - lmass;
  if(tmu < 0.0){
    cost = -999;
    tmu = -999;
    area = 0.0;
    return -999;
  }

  double thisE = tmu + lmass;
  double thisp = sqrt( thisE * thisE - lmass * lmass);
  double numerator =  Enu * Enu + thisp * thisp - dq3 * dq3;
  double denominator = 2.0 * thisp * Enu;
  if(denominator <= 0.0 ){
    cost = 0.0;
    if(denominator < 0.0){
      return -999;
    }
  }
  else cost = numerator / denominator;

  if(TMath::Abs(numerator) > TMath::Abs(denominator)){
    cost = -999;
    tmu = -999;
    area = 0.0;
    return -999;
  }

  // xCrossSect is not yet in the right units for this particular case.
  // need areaElement to go dsigma/dTmudcost to dsigma/dq0dq3
  // Recompute the area element jacobian

  // dT/dq0 x dC/dq3 - dT/dq3 x dC/dq0
  double areaElement = 0.0;
  //double veryCloseToZero = 0.000001;  // in GeV, this should be safe.
  double insqrt = 0.0;
  numerator = 0.0;
  insqrt = Enu * Enu - 2.0 * Enu * dq0 + dq0 * dq0 - lmass * lmass;
  numerator = dq3 / Enu;
  if(insqrt < 0.0)areaElement=0.0;
  else areaElement = numerator / TMath::Sqrt(insqrt);
  area = areaElement;

  return 0;
}
//______________________________________________________________________
bool genie::utils::mec::GetTlCostlFromq0q3(
 double q0, double q3, double Enu, double ml, double& Tl, double& costl)
{
  Tl = Enu - q0 - ml;
  if(Tl < 0.) {
    costl = -999;
    Tl    = -999;
    return false;
  }

  double El = Tl + ml;
  double plsq = El * El - ml * ml;
  if(plsq < 0.) {
    costl = -999;
    Tl    = -999;
    return false;
  }
  double pl = TMath::Sqrt(plsq);
  double numerator =  Enu * Enu + pl * pl - q3 * q3;
  double denominator = 2.0 * pl * Enu;
  if(denominator <= 0.0) {
    costl = 0.0;
    if(denominator < 0.0) {
      return false;
    }
  }
  else {
    costl = numerator / denominator;
  }

  if(TMath::Abs(numerator) > TMath::Abs(denominator)){
    costl = -999;
    Tl    = -999;
    return false;
  }

  return true;
}
//______________________________________________________________________
bool genie::utils::mec::Getq0q3FromTlCostl(
 double Tl, double costl, double Enu, double ml, double& q0, double& q3)
{
  q0 = -999;
  q3 = -999;

  double El = Tl + ml;
  q0 = Enu - El;
  if(q0 < 0) {
     q0 = -999;
     q3 = -999;
     return false;
  }

  double pl = TMath::Sqrt(El * El - ml * ml);
  double q3sq = pl * pl + Enu * Enu - 2.0 * pl * Enu * costl;
  if(q3sq < 0) {
     q0 = -999;
     q3 = -999;
     return false;
  }
  q3 = TMath::Sqrt(q3sq);

  return true;
}
//______________________________________________________________________
double genie::utils::mec::J(
  double dq0, double dq3, double Enu, double ml)
{
  // dT/dq0 x dC/dq3 - dT/dq3 x dC/dq0
  double area = 0.0;
  double insqrt = 0.0;
  insqrt = Enu * Enu - 2.0 * Enu * dq0 + dq0 * dq0 - ml * ml;
  double numerator = dq3 / Enu;
  if(insqrt < 0.0) {
     area=0.0;
  }
  else {
     area = numerator / TMath::Sqrt(insqrt);
  }
  return area;
}
//______________________________________________________________________
double genie::utils::mec::Qvalue(int targetpdg, int nupdg)
{
// The had tensor calculation requires parameters that describe the
// nuclear density. For the Valencia model they are taken by
// C. Garcia-Recio, Nieves, Oset Nucl.Phys.A 547 (1992) 473-487 Table I
// and simplifies that the p and n density are not necessarily the same?
// Standard tables such as deVries et al are similar ~5%
// This is the only nuclear input on the hadron side of the Hadron Tensor.
// and is what makes PseudoPb and PseudoFe different than Ni56 and Rf208
//

  int nearpdg = targetpdg;

  if(nupdg > 0){
    // get Z+1 for CC interaction e.g. 1000210400 from 1000200400
    nearpdg = targetpdg + 10000;
  } else {
    // get Z-1 for CC interaction e.g. 1000210400 from 1000200400
    nearpdg = targetpdg - 10000;
  }

  TParticlePDG *partf = PDGLibrary::Instance()->Find(nearpdg);
  TParticlePDG *parti = PDGLibrary::Instance()->Find(targetpdg);

  if(NULL == partf || NULL == parti){
    // maybe not every valid nucleus in the table has Z+1 or Z-1
    // for example, Ca40 did not.
    LOG("MECUtils", pFATAL)
      << "Can't get qvalue, nucleus " << targetpdg << " or " << nearpdg
      << " is not in the table of nuclei in /data/evgen/pdg ";
    exit(-1);
  }

  double massf = partf->Mass();
  double massi = parti->Mass();

  // keep this here as a reminder of what not to do
  // +/- emass;  // subtract electron from Z+1, add it to Z-1
  // the lookup table gives the nucleus mass, not the atomic
  // mass, so don't need this.

  return massf - massi;  // in GeV.
}
