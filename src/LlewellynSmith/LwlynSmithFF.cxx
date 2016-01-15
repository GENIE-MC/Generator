//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   Renamed LlewellynSmithModel -> LwlynSmithFF

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "ElFF/ELFormFactors.h"
#include "ElFF/ELFormFactorsModelI.h"
#include "ElFF/TransverseEnhancementFFModel.h"
#include "Conventions/Constants.h"
#include "LlewellynSmith/LwlynSmithFF.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LwlynSmithFF::LwlynSmithFF() :
QELFormFactorsModelI()
{

}
//____________________________________________________________________________
LwlynSmithFF::LwlynSmithFF(string name) :
QELFormFactorsModelI(name)
{

}
//____________________________________________________________________________
LwlynSmithFF::LwlynSmithFF(string name, string config) :
QELFormFactorsModelI(name, config)
{

}
//____________________________________________________________________________
LwlynSmithFF::~LwlynSmithFF()
{
  if (fCleanUpfElFFModel) {
    delete fElFFModel;
  }
}
//____________________________________________________________________________
double LwlynSmithFF::StrangeF1V(const Interaction * interaction) const
{
  double f1p = this->F1P(interaction); 
  double f1n = this->F1N(interaction);
  double value = 0.;

  const XclsTag &      xcls       = interaction->ExclTag();
  int pdgc = xcls.StrangeHadronPdg();

  if (pdgc == kPdgSigmaM)        value = (f1p + 2 * f1n);
  else if (pdgc == kPdgLambda)   value = -kSqrt3 / kSqrt2 * f1p;
  else if (pdgc == kPdgSigma0)   value = kSqrt2 / 2 * (f1p + 2 * f1n);

  return value;
}
//____________________________________________________________________________
double LwlynSmithFF::StrangexiF2V(const Interaction * interaction) const
{
  const XclsTag &      xcls       = interaction->ExclTag();
  int pdgc = xcls.StrangeHadronPdg();
  
  double f2p = this->F2P(interaction);
  double f2n = this->F2N(interaction);
  double MH = PDGLibrary::Instance()->Find(pdgc)->Mass();  
  const InitialState & init_state = interaction->InitState();
  double MN   = init_state.Tgt().HitNucMass();
  double value = 0.;

  if (pdgc == kPdgSigmaM)      
    value = (fMuP * f2p +  fMuN * f2n)* 2 * MN / (MN + MH);
  else if (pdgc == kPdgLambda) 
    value = (-kSqrt3 / kSqrt2 * fMuP * f2p) * 2 * MN / (MN + MH);
  else if (pdgc == kPdgSigma0) 
    value = kSqrt2 / 2 * (fMuP * f2p + fMuN * f2n) * 2 * MN / (MN + MH);

  return value;
}

//____________________________________________________________________________
double LwlynSmithFF::StrangeFA(const Interaction * interaction) const
{

  const XclsTag &      xcls       = interaction->ExclTag();
  int pdgc = xcls.StrangeHadronPdg();
  double value = 0.;

  if (pdgc == kPdgSigmaM)       value =  -1 * (1 - 2 * fFDratio);
  else if (pdgc == kPdgLambda)  value =  -1 / kSqrt6 * (1 + 2 * fFDratio);
  else if (pdgc == kPdgSigma0)  value =  -1 * kSqrt2 / 2 * (1 - 2 * fFDratio);

  return value;
}
//____________________________________________________________________________
double LwlynSmithFF::F1P(const Interaction * interaction) const
{ 
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gep() - t * fELFF.Gmp());
}
//____________________________________________________________________________
double LwlynSmithFF::F2P(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gmp() - fELFF.Gep());
}
//____________________________________________________________________________
double LwlynSmithFF::F1N(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gen() - t * fELFF.Gmn());
}
//____________________________________________________________________________
double LwlynSmithFF::F2N(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gmn() - fELFF.Gen());
}
//____________________________________________________________________________
double LwlynSmithFF::F1V(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gve = this->GVE(interaction);
  double gvm = this->GVM(interaction);

  double F1V = (gve - t*gvm) / (1-t);
  return F1V;
}
//____________________________________________________________________________
double LwlynSmithFF::xiF2V(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gve = this->GVE(interaction);
  double gvm = this->GVM(interaction);

  double xiF2V = (gvm-gve) / (1-t);
  return xiF2V;
}
//____________________________________________________________________________
double LwlynSmithFF::FA(const Interaction * interaction) const
{
  // get scattering parameters
  const Kinematics & kine = interaction->Kine();
  double q2 = kine.q2();

  // calculate FA(q2)
  double dn = TMath::Power(1.-q2/fMa2, 2);
  double FA = fFA0/dn;
  return FA;
}
//____________________________________________________________________________
double LwlynSmithFF::Fp(const Interaction * interaction) const
{
  // get momentum transfer
  const Kinematics & kine = interaction->Kine();
  double q2 = kine.q2();

  // get struck nucleon mass & set pion mass
  const InitialState & init_state = interaction->InitState();
  double MN   = init_state.Tgt().HitNucMass();
  double MN2  = TMath::Power(MN, 2);
  double Mpi  = kPionMass;
  double Mpi2 = TMath::Power(Mpi, 2);

  // calculate FA
  double fa = this->FA(interaction);

  // calculate Fp
  double Fp = 2. * MN2 * fa/(Mpi2-q2);
  return Fp;
}
//____________________________________________________________________________
void LwlynSmithFF::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LwlynSmithFF::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LwlynSmithFF::LoadConfig(void)
{
// Load configuration data from its configuration Registry (or global defaults)
// to private data members
	
  fElFFModel = 0;

  AlgConfigPool * confp = AlgConfigPool::Instance();
  Registry * gc = confp->GlobalParameterList();

  // load elastic form factors model
  RgAlg form_factors_model = fConfig->GetAlgDef(
                 "ElasticFormFactorsModel", gc->GetAlg("ElasticFormFactorsModel"));
  fElFFModel =  dynamic_cast<const ELFormFactorsModelI *> (
                                          this->SubAlg("ElasticFormFactorsModel"));
  assert(fElFFModel);
  fCleanUpfElFFModel = false;
  if(gc->GetBoolDef("UseElFFTransverseEnhancement", false)) {
    const ELFormFactorsModelI* sub_alg = fElFFModel;
    RgAlg transverse_enhancement = gc->GetAlg("TransverseEnhancement");
    fElFFModel = dynamic_cast<const ELFormFactorsModelI *> (
        AlgFactory::Instance()->AdoptAlgorithm(
            transverse_enhancement.name, transverse_enhancement.config));
    dynamic_cast<const TransverseEnhancementFFModel*>(fElFFModel)->SetElFFBaseModel(
        sub_alg);
    fCleanUpfElFFModel = true;
  }

  fELFF.SetModel(fElFFModel);  

  // axial mass and Fa(q2=0)
  fMa  = fConfig->GetDoubleDef("Ma",  gc->GetDouble("QEL-Ma"));  // Axial mass
  fFA0 = fConfig->GetDoubleDef("FA0", gc->GetDouble("QEL-FA0")); // FA(q2=0)

  fMa2 = TMath::Power(fMa,2);

  // anomalous magnetic moments
  fMuP = fConfig->GetDoubleDef("MuP", gc->GetDouble("AnomMagnMoment-P"));
  fMuN = fConfig->GetDoubleDef("MuN", gc->GetDouble("AnomMagnMoment-N"));

  // weinberg angle
  double thw = fConfig->GetDoubleDef(
                          "WeinbergAngle", gc->GetDouble("WeinbergAngle"));
  fSin28w  = TMath::Power(TMath::Sin(thw), 2);
  fFDratio = 0.58;  
}
//____________________________________________________________________________
double LwlynSmithFF::tau(const Interaction * interaction) const
{
// computes q^2 / (4 * MNucl^2)

  //-- get kinematics & initial state parameters
  const Kinematics &   kinematics = interaction->Kine();
  const InitialState & init_state = interaction->InitState();
  double q2     = kinematics.q2();
  double Mnucl  = init_state.Tgt().HitNucMass();
  double Mnucl2 = TMath::Power(Mnucl, 2);

  //-- calculate q^2 / (4*Mnuc^2)
  return q2/(4*Mnucl2);
}
//____________________________________________________________________________
double LwlynSmithFF::GVE(const Interaction * interaction) const
{
  //-- compute GVE using CVC

  fELFF.Calculate(interaction);
  double gve = fELFF.Gep() - fELFF.Gen();
  return gve;
}
//____________________________________________________________________________
double LwlynSmithFF::GVM(const Interaction * interaction) const
{
  //-- compute GVM using CVC

  fELFF.Calculate(interaction);
  double gvm = fELFF.Gmp() - fELFF.Gmn();
  return gvm;
}
//____________________________________________________________________________

