//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - August 21, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TSystem.h>

#include "Conventions/Constants.h"
#include "Fragmentation/SchmitzMultiplicityModel.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
SchmitzMultiplicityModel::SchmitzMultiplicityModel() :
MultiplicityProbModelI("genie::SchmitzMultiplicityModel")
{                
  fKNO=0;
}
//____________________________________________________________________________
SchmitzMultiplicityModel::SchmitzMultiplicityModel(string config) :
MultiplicityProbModelI("genie::SchmitzMultiplicityModel", config)
{                 
  fKNO=0;
}
//____________________________________________________________________________
SchmitzMultiplicityModel::~SchmitzMultiplicityModel()
{
  if(fKNO) delete fKNO;
}
//____________________________________________________________________________
TH1D * SchmitzMultiplicityModel::ProbabilityDistribution(
                                        const Interaction * interaction) const
{
  // Compute the average multiplicity for the given interaction:
  // <n> = a + b * ln(W^2)
  
  double alpha = this->SelectOffset(interaction);
  double W     = interaction->GetKinematics().W();

  assert(W>kNeutronMass+kPionMass);
  
  // calculate average charged hadron multiplcity 
  double avn = alpha + fB * 2*TMath::Log(W);
 
  // calculate avergage hadron multiplicity (1.5 x charged)
  avn *= 1.5;

  if(avn <1) {
      LOG("Schmitz", pWARN) << "Average multiplicity too small: " << avn;
      return 0;
  }

  // Create a multiplicity probability distribution

  // set the maximum multiplicity as W = Mneutron + (maxmult-1)*Mpion
  double maxmult = (fForceNeuGenLimit) ?
                    10 : TMath::Floor(1 + (W-kNeutronMass)/kPionMass);
  double minmult = 2;
  int    nbins   = TMath::Nint(maxmult-minmult+1);

  TH1D * prob = new TH1D("", "", nbins, minmult-0.5, maxmult+0.5);

  for(int i = 1; i <= nbins; i++) {

     // KNO distribution is <n>*P(n) vs n/<n>
     double n       = prob->GetBinCenter(i); // bin centre
     double n_avn   = n/avn; // n/<n>
     bool   inrange = n_avn > fKNO->XMin() && n_avn < fKNO->XMax();
     double avnP    = (inrange) ? fKNO->Evaluate(n_avn) : 0; // <n>*P(n)
     double P       = avnP / avn; // P(n)

     SLOG("Schmitz", pDEBUG)
          << "W = " << W << ", <n> = " << avn << ", n/<n> = " << n_avn
          << ", <n>*P = " << avnP << ", P = " << P;

     prob->Fill(n,P);
  }
  //----- Normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  this->ApplyRijk(interaction, prob, true);

  for(int i = 1; i <= nbins; i++) {
     double n = prob->GetBinCenter(i); 
     double P = prob->GetBinContent(i);
     SLOG("Schmitz", pDEBUG)  << "n = " << n << ", P = " << P;
  }

  return prob; // Note: The calling function adopts the object
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::ApplyRijk(
                const Interaction * interaction, TH1D * prob, bool norm) const
{
// Apply the NEUGEN multiplicity probability scaling factors
//
  const InitialState & init_state = interaction->GetInitialState();
  int nu_pdg  = init_state.GetProbePDGCode();
  int nuc_pdg = init_state.GetTarget().StruckNucleonPDGCode();

  double R1=1., R2=1.;

  const ProcessInfo & proc_info = interaction->GetProcessInfo();
  bool isCC = proc_info.IsWeakCC();

  if(pdg::IsNeutrino(nu_pdg) && pdg::IsProton(nuc_pdg))  {
    R1 = (isCC) ? fRvpCCm1 : fRvpNCm1;
    R2 = (isCC) ? fRvpCCm2 : fRvpNCm2;
  } else 
  if(pdg::IsNeutrino(nu_pdg) && pdg::IsNeutron(nuc_pdg)) {
    R1 = (isCC) ? fRvnCCm1 : fRvnNCm1;
    R2 = (isCC) ? fRvnCCm2 : fRvnNCm2;
  } else 
  if(pdg::IsAntiNeutrino(nu_pdg) && pdg::IsProton(nuc_pdg))  {
    R1 = (isCC) ? fRvbpCCm1 :   fRvbpNCm1;
    R2 = (isCC) ? fRvbpCCm2 :   fRvbpNCm2;
  } else 
  if(pdg::IsAntiNeutrino(nu_pdg) && pdg::IsNeutron(nuc_pdg)) {
    R1 = (isCC) ? fRvbnCCm1 : fRvbnNCm1;
    R2 = (isCC) ? fRvbnCCm2 : fRvbnNCm2;
  } else {
    LOG("Schmitz", pERROR) << "Invalid initial state: " << init_state;
  }

  int nbins = prob->GetNbinsX();
  for(int i = 1; i <= nbins; i++) {
     int n = TMath::Nint( prob->GetBinCenter(i) ); 
     if(n==1) prob->SetBinContent(i, R1 * prob->GetBinContent(i));
     if(n==2) prob->SetBinContent(i, R2 * prob->GetBinContent(i));
  }

  // renormalize the histogram?
  if(norm) {
     prob->Scale( 1.0 / prob->Integral("width") );
  }
}
//____________________________________________________________________________
double SchmitzMultiplicityModel::SelectOffset(
                                        const Interaction * interaction) const
{
  const InitialState & init_state = interaction->GetInitialState();
  int nu_pdg  = init_state.GetProbePDGCode();
  int nuc_pdg = init_state.GetTarget().StruckNucleonPDGCode();

  if( pdg::IsNeutrino( nu_pdg ) ) {
      if ( pdg::IsProton(nuc_pdg)  )  return fAvp;
      if ( pdg::IsNeutron(nuc_pdg) )  return fAvn;
      else {
         LOG("Schmitz", pERROR)
                          << "PDG-Code = " << nuc_pdg << " is not a nucleon!";
      }
  } else  if (  pdg::IsAntiNeutrino(nu_pdg) ) {
      if ( pdg::IsProton(nuc_pdg)  )  return fAvbp;
      if ( pdg::IsNeutron(nuc_pdg) )  return fAvbn;
      else {
         LOG("Schmitz", pERROR)
                          << "PDG-Code = " << nuc_pdg << " is not a nucleon!";
      }
  } else {
    LOG("Schmitz", pERROR)<< "PDG-Code = " << nu_pdg << " is not a neutrino!";
  }
  return 0;        
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::LoadConfig(void)
{
// Load config parameters

  // Get the KNO Distribution
  if(fKNO) delete fKNO;

  assert(gSystem->Getenv("GENIE"));

  string basedir    = gSystem->Getenv("GENIE");
  string defknodata = basedir + "/data/kno/KNO.dat";
  string knodata    = fConfig->GetStringDef("kno-data",defknodata);

  LOG("Schmitz", pNOTICE) << "Loading KNO data from: " << knodata;

  fKNO = new Spline(knodata);
  // Load parameters determining the average multiplicity
  fAvp  = fConfig->GetDouble("alpha-vp");
  fAvn  = fConfig->GetDouble("alpha-vn");
  fAvbp = fConfig->GetDouble("alpha-vbp");
  fAvbn = fConfig->GetDouble("alpha-vbn");
  fB    = fConfig->GetDouble("beta");

  // Force NEUGEN upper limit in hadronic multiplicity (to be used only
  // NEUGEN/GENIE comparisons)
  fForceNeuGenLimit = fConfig->GetBoolDef("force-neugen-mult-limit", false);

  // Load NEUGEN multiplicity probability scaling parameters Rijk
  fRvpCCm1  = fConfig->GetDoubleDef("R-vp-CC-m1", 1.0);
  fRvpCCm2  = fConfig->GetDoubleDef("R-vp-CC-m2", 1.0);
  fRvpNCm1  = fConfig->GetDoubleDef("R-vp-NC-m1", 1.0);
  fRvpNCm2  = fConfig->GetDoubleDef("R-vp-NC-m2", 1.0);
  fRvnCCm1  = fConfig->GetDoubleDef("R-vn-CC-m1", 1.0);
  fRvnCCm2  = fConfig->GetDoubleDef("R-vn-CC-m2", 1.0);
  fRvnNCm1  = fConfig->GetDoubleDef("R-vn-NC-m1", 1.0);
  fRvnNCm2  = fConfig->GetDoubleDef("R-vn-NC-m2", 1.0);
  fRvbpCCm1 = fConfig->GetDoubleDef("R-vbp-CC-m1",1.0);
  fRvbpCCm2 = fConfig->GetDoubleDef("R-vbp-CC-m2",1.0);
  fRvbpNCm1 = fConfig->GetDoubleDef("R-vbp-NC-m1",1.0);
  fRvbpNCm2 = fConfig->GetDoubleDef("R-vbp-NC-m2",1.0);
  fRvbnCCm1 = fConfig->GetDoubleDef("R-vbn-CC-m1",1.0);
  fRvbnCCm2 = fConfig->GetDoubleDef("R-vbn-CC-m2",1.0);
  fRvbnNCm1 = fConfig->GetDoubleDef("R-vbn-NC-m1",1.0);
  fRvbnNCm2 = fConfig->GetDoubleDef("R-vbn-NC-m2",1.0);
}
//____________________________________________________________________________

