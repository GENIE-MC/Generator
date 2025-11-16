//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


 For documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TLorentzVector.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/HadronTensors/HadronTensorModelI.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Conventions/Units.h"
#include "Physics/HadronTensors/LabFrameHadronTensorI.h"
#include "Physics/HadronTensors/HadronTensorI.h"

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
double genie::utils::mec::OldTensorContraction(
  int nupdg, int targetpdg,
  double Enu, double Ml, double Tl, double costhl,
  int tensorpdg,
  HadronTensorType_t tensor_type,
  char* tensor_model)
{
  TLorentzVector v4lep;
  TLorentzVector v4Nu(0,0,Enu,Enu); // assuming traveling along z:
  TLorentzVector v4q;
  double q0nucleus;
  double facconv = 0.0389391289e15; // std::pow(0.19733,2)*1e15;

  double myQvalue = 0.0;

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

  //Get q3
  double pl = TMath::Sqrt(v4lep.E() * v4lep.E() - Ml * Ml);
  double q3sq = pl * pl + Enu * Enu - 2.0 * pl * Enu * costhl;
  double q3 = sqrt(q3sq);

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


    // Get the appropriate hadron tensor model object
    const genie::HadronTensorModelI* ht_model
      = dynamic_cast<const genie::HadronTensorModelI*>(
      genie::AlgFactory::Instance()->GetAlgorithm( tensor_model, "Default" ));

    const LabFrameHadronTensorI* tensor_table
      = dynamic_cast<const LabFrameHadronTensorI*>( ht_model->GetTensor(targetpdg,
      tensor_type) );

    double W00 = (tensor_table->tt(q0nucleus,q3)).real();
    double W0Z = (tensor_table->tz(q0nucleus,q3)).real();
    double WXX = (tensor_table->xx(q0nucleus,q3)).real();
    double WXY = -1.0*(tensor_table->xy(q0nucleus,q3)).imag();
    double WZZ = (tensor_table->zz(q0nucleus,q3)).real();

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
      xsec = 0;
    }
  }
  if((Enu>1.15 && Enu<1.25) && (v4q.Vect().Mag()>0.615 && v4q.Vect().Mag()<0.625) &&  (q0<0.500 && q0>0.490))
    LOG("MECUtils", pFATAL)//SB
         << "Got xsec = " << xsec
         << "\n " << part1 << " " << part2
         << "\n w[] : " << w1 << ", " << w2 << ", " << w3 << ", " << w4 << ", " << w5
         << "\n wtotd[] : " << wtotd[0] << ", " << wtotd[1] << ", " << wtotd[2] << ", " << wtotd[3] << ", " << wtotd[4]
         << "\n vec " << v4q.Vect().Mag() << ", " << q0  << ", " << v4q.Px() << ", " << v4q.Py() << ", " << v4q.Pz()
         << "\n v4qX " << v4Nu.X() << ", " << v4lep.X() << ", " << costhl << ", " << sinthl << ", " << modkprime
   << "\n input " << Enu << ", " << Ml << ", " << Tl << ", " <<costhl << ", " << v4q.Vect().Mag() << ", " << q0;

  //LOG("MECUtils", pFATAL) << "xsec(Enu = " << Enu << " GeV, Ml = " << Ml << " GeV; " << "Tl = " << Tl << " GeV, costhl = " << costhl << ") = " << xsec << " x 1E-41 cm^2";

  //return (xsec * (1.0E-39 * units::cm2));
  return (xsec);
}
//___________________________________________________________________________
double genie::utils::mec::GetMaxXSecTlctl( const XSecAlgorithmI& xsec_model,
  const Interaction& inter, const double tolerance, const double safety_factor,
  const int max_n_layers )
{
  // Clone the input interaction so that we can modify the Tl and ctl values
  // without messing anything up
  Interaction* interaction = new Interaction( inter );

  // Choose the appropriate minimum Q^2 value based on the interaction
  // mode (this is important for EM interactions since the differential
  // cross section blows up as Q^2 --> 0)
  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  if ( interaction->ProcInfo().IsEM() ) Q2min = genie::utils::kinematics
    ::electromagnetic::kMinQ2Limit; // EM limit

  const double Enu = interaction->InitState().ProbeE( kRfLab );
  const double ProbeMass = interaction->InitState().Probe()->Mass();
  const double LepMass = interaction->FSPrimLepton()->Mass();

  // Determine the bounds for the current layer being scanned.

  // Maximum energy transfer to consider
  const double Q0Max = std::min( utils::mec::Q0LimitMaxXSec, Enu );

  // Maximum momentum transfer to consider
  const double QMagMax = std::min( utils::mec::QMagLimitMaxXSec, 2.*Enu );

  // Translate these bounds into limits for the lepton
  // kinetic energy and Q^2
  const double TMax = std::max( 0., Enu - LepMass );
  double CurTMin = std::max( 0., TMax - Q0Max );
  double CurTMax = TMax;
  double CurQ2Min = Q2min;
  // Overshoots the true maximum Q^2, but not so severely that it's expected to
  // be a problem.
  double CurQ2Max = QMagMax * QMagMax;

  // This is based on the technique used to compute the maximum differential
  // cross section for rejection sampling in QELEventGenerator. A brute-force
  // scan over the two-dimensional kPSTlctl phase space is used to find the
  // maximum cross section. Multiple layers are used to "zoom in" on the
  // vicinity of the maximum. Tighter and tighter layers are made until
  // the maximum number of iterations is reached or the xsec stabilizes
  // to within a given tolerance.
  const int N_T = 10; // number of kinetic energy steps per layer
  const int N_Q2 = 10; // number of scattering cosine steps per layer
  double T_at_xsec_max = CurTMax;
  double Q2_at_xsec_max = CurQ2Max;
  double cur_max_xsec = -1.;

  for ( int ilayer = 0; ilayer < max_n_layers; ++ilayer ) {

    double last_layer_xsec_max = cur_max_xsec;

    double T_increment = ( CurTMax - CurTMin ) / ( N_T - 1 );
    double Q2_increment   = ( CurQ2Max - CurQ2Min ) / ( N_Q2 - 1 );

    // Scan through the 2D phase space using the bounds of the current layer
    for ( int iT = 0; iT < N_T; ++iT ) {
      double T = CurTMin + iT * T_increment;
      for ( int iQ2 = 0; iQ2 < N_Q2; ++iQ2 ) {
        double Q2 = CurQ2Min + iQ2 * Q2_increment;

        // Compute the scattering cosine at fixed Tl given a value of Q^2.
        double Elep = T + LepMass;
        double pnu = std::sqrt( std::max(0., Enu*Enu - ProbeMass*ProbeMass) );
        double plep = std::sqrt( std::max(0., std::pow(T + LepMass, 2)
          - LepMass*LepMass) );

        double Costh = plep==0?0:( 2.*Enu*Elep - ProbeMass*ProbeMass - LepMass*LepMass
          - Q2 ) / ( 2. * pnu * plep );
        // Respect the bounds of the cosine function.
        Costh = std::min( std::max(-1., Costh), 1. );

        // Update the values of Tl and ctl in the interaction, then
        // compute the differential cross section
        interaction->KinePtr()->SetKV( kKVTl, T );
        interaction->KinePtr()->SetKV( kKVctl, Costh );
        double xsec = xsec_model.XSec( interaction, kPSTlctl );

        // If the current differential cross section exceeds the previous
        // maximum value, update the stored value and information about where it
        // lies in the kPSTlctl phase space.
        if ( xsec > cur_max_xsec ) {
          T_at_xsec_max = T;
          Q2_at_xsec_max = Q2;
          cur_max_xsec = xsec;
        }

      } // Done with cos(theta) scan
    } // Done with lepton kinetic energy scan

    LOG("MEC", pDEBUG) << "Layer " << ilayer << " has max xsec = "
      << cur_max_xsec << " for T = " << T_at_xsec_max << ", Q2 = "
      << Q2_at_xsec_max;

    // ** Calculate the range for the next layer **

    // Don't let the minimum kinetic energy fall below zero
    CurTMin = std::max( 0., T_at_xsec_max - T_increment );

    // We know a priori that the maximum cross section should occur near the
    // maximum lepton kinetic energy, so keep the old maximum value when
    // recalculating the new range. It was originally set to something near the
    // maximum allowed value, so this should zero in on the appropriate region
    // in future scans. The updated maximum from the old code is shown below,
    // but due to problems with round-off in the kinematic limits, it can
    // sometimes miss the edge of the phase space if the T_increment value is
    // coarse enough. - S. Gardiner, 1 July 2020
    CurTMax = TMax; // T_at_xsec_max + T_increment;

    // We use the same idea here. A priori, we know that the maximum cross section
    // should lie near the minimum Q^2 value. Round-off gives us the same problems
    // in finding the low Q^2 edge of phase space, so just keep the minimum equal
    // to the cut for the next scan. The old assignment (which mostly worked
    // but had problems for a coarse Q2_increment value) is commented out below.
    // - S. Gardiner, 1 July 2020
    CurQ2Min = Q2min; // std::max( Q2min, Q2_at_xsec_max - Q2_increment );
    CurQ2Max = Q2_at_xsec_max + Q2_increment;

    // If the xsec has stabilized within the requested tolerance, then
    // stop adding more layers
    double improvement_factor = cur_max_xsec / last_layer_xsec_max;
    if ( ilayer && std::abs(improvement_factor - 1.) < tolerance ) break;

  } // Done with iterations over layers


  // We're done. Delete the cloned interaction object before returning the
  // estimate of the maximum differential cross section.
  delete interaction;

  double XSecMax = cur_max_xsec * safety_factor;
  return XSecMax;
}
//___________________________________________________________________________
genie::utils::mec::gsl::d2Xsec_dTCosth::d2Xsec_dTCosth(const XSecAlgorithmI * m, const Interaction & i, 
						       const double Enu, const double LepMass, const double Factor ) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i),
fEnu(Enu),
fLepMass(LepMass),
fFactor(Factor) 
{

}
//____________________________________________________________________________
genie::utils::mec::gsl::d2Xsec_dTCosth::~d2Xsec_dTCosth()
{

}
//____________________________________________________________________________
unsigned int genie::utils::mec::gsl::d2Xsec_dTCosth::NDim(void) const
{
  return 2;
}
//____________________________________________________________________________
double genie::utils::mec::gsl::d2Xsec_dTCosth::DoEval(const double * xin) const
{
// inputs:
//    T [GeV]
//    cos(theta)
// outputs:
//   differential cross section (hbar=c=1 units)
//

  double T     = xin[0];
  double costh = xin[1];

  Kinematics * kinematics = fInteraction.KinePtr();
  kinematics->SetKV(kKVTl, T);
  kinematics->SetKV(kKVctl, costh);

  double Q0 = 0 ;
  double Q3 = 0 ; 
  genie::utils::mec::Getq0q3FromTlCostl(T, costh, fEnu, fLepMass, Q0, Q3);
  
  kinematics ->SetKV(kKVQ0, Q0) ; 
  kinematics ->SetKV(kKVQ3, Q3) ; 

  double xsec = fModel->XSec( &fInteraction, kPSTlctl);
  return fFactor * xsec;
  
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim *
genie::utils::mec::gsl::d2Xsec_dTCosth::Clone() const
{
  return
    new genie::utils::mec::gsl::d2Xsec_dTCosth(fModel,fInteraction, fEnu, fLepMass, fFactor );
}
//____________________________________________________________________________

