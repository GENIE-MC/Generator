//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
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
//______________________________________________________________________
double genie::utils::mec::TensorContraction(
  const Interaction * interaction,  
  int tensor_pdg,
  MECHadronTensor::MECHadronTensorType_t tensor_type)
{
  int target_pdg = interaction->InitState().Tgt().Pdg();
  int nu_pdg     = interaction->InitState().ProbePdg();      
  double Enu     = interaction->InitState().ProbeE(kRfLab);
  double Ml      = interaction->FSPrimLepton()->Mass();
  double Tl      = interaction->Kine().GetKV(kKVTl);
  double costhl  = interaction->Kine().GetKV(kKVctl);
      
  return genie::utils::mec::TensorContraction(
    nu_pdg, target_pdg, Enu, Ml, Tl, costhl, tensor_pdg, tensor_type);
}
//______________________________________________________________________
double genie::utils::mec::TensorContraction(
  int nupdg, int targetpdg,
  double Enu, double Ml, double Tl, double costhl,
  int tensorpdg,
  MECHadronTensor::MECHadronTensorType_t tensor_type)
{
  TLorentzVector v4lep;
  TLorentzVector v4Nu(0,0,Enu,Enu); // assuming traveling along z:
  TLorentzVector v4q;
  double q0nucleus;  
  double facconv = 0.0389391289e15; // std::pow(0.19733,2)*1e15;
  
  double myQvalue = genie::utils::mec::Qvalue(targetpdg, nupdg);
  
  // Angles
  double sinthl = 1. - costhl * costhl;
  if(sinthl < 0.0) sinthl = 0.0;
  else sinthl = TMath::Sqrt(sinthl);
  
  double Cosh = TMath::Cos(TMath::ACos(costhl)/2.);
  double Sinh = TMath::Sin(TMath::ACos(costhl)/2.);
  
  // Lepton
  v4lep.SetE( Tl + Ml );
  // energy transfer from the lepton
  double q0 = v4Nu.E() - v4lep.E(); 
  // energy transfer that actually gets to the nucleons
  q0nucleus = q0 - myQvalue;
  
  // Define some calculation placeholders
  double part1, part2;
  double modkprime ;  
  double w1, w2, w3, w4, w5;
  double wtotd[5];

  double xsec = 0;
  if (q0nucleus <= 0){
    //nothing transfered to nucleus thus no 2p2h
    xsec = 0.;
  } else {
    // energy was transfered to the nucleus
    modkprime = std::sqrt(TMath::Power(v4lep.E(),2)-TMath::Power(Ml,2));
    v4lep.SetX(modkprime*sinthl);
    v4lep.SetY(0);
    v4lep.SetZ(modkprime*costhl);
    
    //q: v4q = v4Nu - v4lep; 
    v4q.SetE(q0nucleus);
    v4q.SetX(v4Nu.X() - v4lep.X());
    v4q.SetY(v4Nu.Y() - v4lep.Y());
    v4q.SetZ(v4Nu.Z() - v4lep.Z());
    
    MECHadronTensor * hadtensor = MECHadronTensor::Instance();
    const vector <genie::BLI2DNonUnifGrid *> &
         tensor_table = hadtensor->TensorTable(tensorpdg, tensor_type);
    
    for (int i=0 ; i < 5; i++){
      wtotd[i] = tensor_table[i]->Evaluate(v4q.Vect().Mag(),v4q.E());
    }
    
    // calculate hadron tensor components
    // these are footnote 2 of Nieves PRC 70 055503
    double W00 = wtotd[0];
    double W0Z = wtotd[1];
    double WXX = wtotd[2];
    double WXY = wtotd[3];
    double WZZ = wtotd[4];
         
    w1=WXX/2.;
    w2=(W00+WXX+(q0*q0/(v4q.Vect().Mag()*v4q.Vect().Mag())
                           *(WZZ-WXX))
          -(2.*q0/v4q.Vect().Mag()*W0Z))/2.;
    w3=WXY/v4q.Vect().Mag();
    w4=(WZZ-WXX)/(2.*v4q.Vect().Mag()*v4q.Vect().Mag());
    w5=(W0Z-(q0/v4q.Vect().Mag()*(WZZ-WXX)))/v4q.Vect().Mag();
    //w6 we have no need for w6, noted at the end of IIA.
  
    // adjust for anti neutrinos

    if (nupdg < 0) w3 = -1. * w3;
    
    // calculate cross section, in parts
    double xw1 = w1*costhl;
    double xw2 = w2/2.*costhl;
    double xw3 = w3/2.*(v4lep.E()+modkprime-(v4lep.E()+v4Nu.E())*costhl);
    double xw4 = w4/2.*(Ml*Ml*costhl+2.*v4lep.E()*(v4lep.E()+modkprime)*Sinh*Sinh);
    double xw5 = w5*(v4lep.E()+modkprime)/2.;
    part1 = xw1 - xw2 + xw3 + xw4 - xw5;
    
    double yw1 = 2.*w1*Sinh*Sinh;
    double yw2 = w2*Cosh*Cosh;
    double yw3 = w3*(v4lep.E()+v4Nu.E())*Sinh*Sinh;
    double yw4 = Ml*Ml*part1/(v4lep.E()*(v4lep.E()+modkprime));
    part2 = yw1 + yw2 - yw3 + yw4;
    
    xsec = modkprime*v4lep.E()*kGF2*2./kPi*part2*facconv;

    if( ! (xsec >= 0.0) ){
      // Traps negative numbers and nan's.
      // Sometimes costhl is just over 1.0 due to fluke or numerical
      // precision, so sinthl would be undefined.
      LOG("MECUtils", pDEBUG)
         << "Got xsec = " << xsec
         << "\n " << part1 << " " << part2
         << "\n w[] : " << w1 << ", " << w2 << ", " << w3 << ", " << w4 << ", " << w5
         << "\n wtotd[] : " << wtotd[0] << ", " << wtotd[1] << ", " << wtotd[2] << ", " << wtotd[3] << ", " << wtotd[4]
         << "\n vec " << v4q.Vect().Mag() << ", " << q0  << ", " << v4q.Px() << ", " << v4q.Py() << ", " << v4q.Pz()
         << "\n v4qX " << v4Nu.X() << ", " << v4lep.X() << ", " << costhl << ", " << sinthl << ", " << modkprime;
    
      xsec = 0;
    }
  }

//  LOG("MECUtils", pDEBUG)
//     << "xsec(Enu = " << Enu << " GeV, Ml = " << Ml << " GeV; "
//     << "Tl = " << Tl << " GeV, costhl = " << costhl << ") = " 
//     << xsec << " x 1E-41 cm^2";

  return (xsec * (1.0E-41 * units::cm2));
}
//______________________________________________________________________



