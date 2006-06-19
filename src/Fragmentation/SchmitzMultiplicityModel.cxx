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

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Controls.h"
#include "Conventions/Constants.h"
#include "Fragmentation/SchmitzMultiplicityModel.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;

//____________________________________________________________________________
SchmitzMultiplicityModel::SchmitzMultiplicityModel() :
MultiplicityProbModelI("genie::SchmitzMultiplicityModel")
{                
  fKNO      = 0;
  fMultProb = 0;
}
//____________________________________________________________________________
SchmitzMultiplicityModel::SchmitzMultiplicityModel(string config) :
MultiplicityProbModelI("genie::SchmitzMultiplicityModel", config)
{ 
  fKNO      = 0;
  fMultProb = 0;  
}
//____________________________________________________________________________
SchmitzMultiplicityModel::~SchmitzMultiplicityModel()
{
  if(fKNO)      delete fKNO;
  if(fMultProb) delete fMultProb;
}
//____________________________________________________________________________
const TH1D & SchmitzMultiplicityModel::ProbabilityDistribution(
                                        const Interaction * interaction) const
{
  // Reset previously computed multiplicity distribution
  SLOG("Schmitz", pDEBUG)<< "Resetting multiplicity probability distribution";
  fMultProb->Reset();

  // Compute the average charged hadron multiplicity as: <n> = a + b*ln(W^2)
  // Calculate avergage hadron multiplicity (= 1.5 x charged hadron mult.)

  double W = utils::kinematics::CalcW(interaction);
  if(W<kNeutronMass+kPionMass) {
    SLOG("Schmitz", pWARN) 
        << "Low mass, W=" << W 
               << "! Returning empty multiplicity probability distribution";
    return *fMultProb;
  }

  double alpha = this->SelectOffset(interaction);
  double avn   = alpha + fB * 2*TMath::Log(W);
  avn *= 1.5;

  SLOG("Schmitz", pINFO) 
             << "Average hadronic multiplicity (W=" << W << ") = " << avn;

  // Find the maximum multiplicity as W = Mneutron + (maxmult-1)*Mpion
  double maxmult = TMath::Floor(1 + (W-kNeutronMass)/kPionMass);

  // If required force the NeuGEN maximum multiplicity limit
  // Note: use for NEUGEN/GENIE comparisons, not physics MC production
  if(fForceNeuGenLimit & maxmult>10) maxmult=10;

  SLOG("Schmitz", pDEBUG) << "Computed maximum multiplicity = " << maxmult;

  // If it exceeds the upper end at the existing TH1D then recreate it
  if(maxmult + 0.5 > fMultProb->GetXaxis()->GetXmax()) {
    SLOG("Schmitz", pNOTICE)  
      << "Increasing the maximum  multiplicity from "
             <<  fMultProb->GetXaxis()->GetXmax() << " to " << maxmult;
    this->CreateProbHist(maxmult);
  }

  // Compute the multiplicity probabilities values up to the bin corresponding 
  // to the computed maximum multiplicity

  if(maxmult>2) {
    int nbins = fMultProb->FindBin(maxmult);

    for(int i = 1; i <= nbins; i++) {
       // KNO distribution is <n>*P(n) vs n/<n>
       double n       = fMultProb->GetBinCenter(i); // bin centre
       double n_avn   = n/avn; // n/<n>
       bool   inrange = n_avn > fKNO->XMin() && n_avn < fKNO->XMax();
       double avnP    = (inrange) ? fKNO->Evaluate(n_avn) : 0; // <n>*P(n)
       double P       = avnP / avn; // P(n)

       SLOG("Schmitz", pDEBUG)
          << "n = " << n << " (n/<n> = " << n_avn
                            << ", <n>*P = " << avnP << ") => P = " << P;
       fMultProb->Fill(n,P);
    }
  } else {
       fMultProb->Fill(2,1.);
  }
  double integral = fMultProb->Integral("width");

  if(integral>0) {
    // Normalize the probability distribution
    fMultProb->Scale( 1.0 / fMultProb->Integral("width") );

    // Apply the NeuGEN probability scaling factors -if requested-
    if(fApplyRijk) {
      SLOG("Schmitz", pDEBUG) << "Applying NeuGEN scaling factors";
      // Only do so for W<Wcut
      if(W<fWcut) {
        this->ApplyRijk(interaction, fRenormalize);
      } else {
        SLOG("Schmitz", pDEBUG)  
              << "W = " << W << " < Wcut = " << fWcut 
                                << " - Will not apply scaling factors";
      }//<wcut?
    }//apply?
  } else {
     SLOG("Schmitz", pDEBUG) << "probability distribution integral = 0";
  }

  return *fMultProb;
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::ApplyRijk(
                            const Interaction * interaction, bool norm) const
{
// Apply the NEUGEN multiplicity probability scaling factors
//
  const InitialState & init_state = interaction->GetInitialState();
  int nu_pdg  = init_state.GetProbePDGCode();
  int nuc_pdg = init_state.GetTarget().StruckNucleonPDGCode();

  double R2=1., R3=1.;

  const ProcessInfo & proc_info = interaction->GetProcessInfo();
  bool isCC = proc_info.IsWeakCC();

  if(pdg::IsNeutrino(nu_pdg) && pdg::IsProton(nuc_pdg))  {
    R2 = (isCC) ? fRvpCCm2 : fRvpNCm2;
    R3 = (isCC) ? fRvpCCm3 : fRvpNCm3;
  } else 
  if(pdg::IsNeutrino(nu_pdg) && pdg::IsNeutron(nuc_pdg)) {
    R2 = (isCC) ? fRvnCCm2 : fRvnNCm2;
    R3 = (isCC) ? fRvnCCm3 : fRvnNCm3;
  } else 
  if(pdg::IsAntiNeutrino(nu_pdg) && pdg::IsProton(nuc_pdg))  {
    R2 = (isCC) ? fRvbpCCm2 :   fRvbpNCm2;
    R3 = (isCC) ? fRvbpCCm3 :   fRvbpNCm3;
  } else 
  if(pdg::IsAntiNeutrino(nu_pdg) && pdg::IsNeutron(nuc_pdg)) {
    R2 = (isCC) ? fRvbnCCm2 : fRvbnNCm2;
    R3 = (isCC) ? fRvbnCCm3 : fRvbnNCm3;
  } else {
    LOG("Schmitz", pERROR) << "Invalid initial state: " << init_state;
  }

  int nbins = fMultProb->GetNbinsX();
  for(int i = 1; i <= nbins; i++) {
     int n = TMath::Nint( fMultProb->GetBinCenter(i) ); 

     double R=1;
     if      (n==2) R=R2;
     else if (n==3) R=R3;

     if(n==2 || n==3) {
        double P   = fMultProb->GetBinContent(i);
        double Psc = R*P;
        LOG("Schmitz", pDEBUG) 
          << "n=" << n << "/ Scaling factor R = " 
                              << R << "/ P " << P << " --> " << Psc;
        fMultProb->SetBinContent(i, Psc);
     }
     if(n>3) break;
  }

  // renormalize the histogram?
  if(norm) {
     fMultProb->Scale( 1.0 / fMultProb->Integral("width") );
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
void SchmitzMultiplicityModel::CreateProbHist(double maxmult) const
{
  if(fMultProb) delete fMultProb;

  double minmult = 2;
  int    nbins   = TMath::Nint(maxmult-minmult+1);

  fMultProb = new TH1D("multprob", 
      "hadronic multiplicity distribution", nbins, minmult-0.5, maxmult+0.5);
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
  this->CreateProbHist(kMaxMultiplicity);
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
  this->CreateProbHist(kMaxMultiplicity);
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::LoadConfig(void)
{
// Load config parameters

  // Access global defaults to use in case of missing parameters
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Delete the KNO spline from previous configuration of this instance
  if(fKNO) delete fKNO;

  assert(gSystem->Getenv("GENIE"));

  string basedir    = gSystem->Getenv("GENIE");
  string defknodata = basedir + "/data/kno/KNO.dat";
  string knodata    = fConfig->GetStringDef("kno-data",defknodata);

  LOG("Schmitz", pNOTICE) << "Loading KNO data from: " << knodata;

  fKNO = new Spline(knodata);

  // Load parameters determining the average multiplicity
  fAvp  = fConfig->GetDoubleDef("alpha-vp",  gc->GetDouble("KNO-Alpha-vp") );
  fAvn  = fConfig->GetDoubleDef("alpha-vn",  gc->GetDouble("KNO-Alpha-vn") ); 
  fAvbp = fConfig->GetDoubleDef("alpha-vbp", gc->GetDouble("KNO-Alpha-vbp")); 
  fAvbn = fConfig->GetDoubleDef("alpha-vbn", gc->GetDouble("KNO-Alpha-vbn"));
  fB    = fConfig->GetDoubleDef("beta",      gc->GetDouble("KNO-Beta")     );

  // Force NEUGEN upper limit in hadronic multiplicity (to be used only
  // NEUGEN/GENIE comparisons)
  fForceNeuGenLimit = fConfig->GetBoolDef("force-neugen-mult-limit", false);

  // Check whether to apply NEUGEN multiplicity probability scaling params
  // Rijk and whether to re-normalize the probability distribution 
  // (note that when called from a DIS cross section algorithm under a DIS
  // RES joining scheme the reduction in the integral of the multiplicity
  // probability distribution would be interpreted as a DIS reduction factor
  // so the distribution should not be re-normalized)
  fApplyRijk   = fConfig->GetBoolDef("apply-neugen-scaling-factors", true);
  fRenormalize = fConfig->GetBoolDef("renormalize-after-scaling",    true);

  // Load Wcut determining the phase space area where the multiplicity prob.
  // scaling factors would be applied -if requested-
  fWcut = fConfig->GetDoubleDef("Wcut",gc->GetDouble("Wcut"));

  // Load NEUGEN multiplicity probability scaling parameters Rijk
  fRvpCCm2  = fConfig->GetDoubleDef(
                      "R-vp-CC-m2", gc->GetDouble("DIS-HMultWgt-vp-CC-m2"));
  fRvpCCm3  = fConfig->GetDoubleDef(
                      "R-vp-CC-m3", gc->GetDouble("DIS-HMultWgt-vp-CC-m3"));
  fRvpNCm2  = fConfig->GetDoubleDef(
                      "R-vp-NC-m2", gc->GetDouble("DIS-HMultWgt-vp-NC-m2"));
  fRvpNCm3  = fConfig->GetDoubleDef(
                      "R-vp-NC-m3", gc->GetDouble("DIS-HMultWgt-vp-NC-m3"));
  fRvnCCm2  = fConfig->GetDoubleDef(
                      "R-vn-CC-m2", gc->GetDouble("DIS-HMultWgt-vn-CC-m2"));
  fRvnCCm3  = fConfig->GetDoubleDef(
                      "R-vn-CC-m3", gc->GetDouble("DIS-HMultWgt-vn-CC-m3"));
  fRvnNCm2  = fConfig->GetDoubleDef(
                      "R-vn-NC-m2", gc->GetDouble("DIS-HMultWgt-vn-NC-m2"));
  fRvnNCm3  = fConfig->GetDoubleDef(
                      "R-vn-NC-m3", gc->GetDouble("DIS-HMultWgt-vn-NC-m3"));
  fRvbpCCm2 = fConfig->GetDoubleDef(
                     "R-vbp-CC-m2",gc->GetDouble("DIS-HMultWgt-vbp-CC-m2"));
  fRvbpCCm3 = fConfig->GetDoubleDef(
                     "R-vbp-CC-m3",gc->GetDouble("DIS-HMultWgt-vbp-CC-m3"));
  fRvbpNCm2 = fConfig->GetDoubleDef(
                     "R-vbp-NC-m2",gc->GetDouble("DIS-HMultWgt-vbp-NC-m2"));
  fRvbpNCm3 = fConfig->GetDoubleDef(
                     "R-vbp-NC-m3",gc->GetDouble("DIS-HMultWgt-vbp-NC-m3"));
  fRvbnCCm2 = fConfig->GetDoubleDef(
                     "R-vbn-CC-m2",gc->GetDouble("DIS-HMultWgt-vbn-CC-m2"));
  fRvbnCCm3 = fConfig->GetDoubleDef(
                     "R-vbn-CC-m3",gc->GetDouble("DIS-HMultWgt-vbn-CC-m3"));
  fRvbnNCm2 = fConfig->GetDoubleDef(
                     "R-vbn-NC-m2",gc->GetDouble("DIS-HMultWgt-vbn-NC-m2"));
  fRvbnNCm3 = fConfig->GetDoubleDef(
                     "R-vbn-NC-m3",gc->GetDouble("DIS-HMultWgt-vbn-NC-m3"));
}
//____________________________________________________________________________
