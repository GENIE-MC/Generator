//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   Renamed QELPXSec -> LwlynSmithQELCCPXSec
 @ Mar 18, 2016 - JJ (SD)
   Moved code to average over initial nucleons from QELXSec to the Integral()
   method here. For each nucleon, generate a struck nucleon position, then a
   momentum, then integrate.
 @ 2015 - AF
   Added FullDifferentialXSec method to work with QELEventGenerator
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "Base/QELFormFactors.h"
#include "Base/QELFormFactorsModelI.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "EVGModules/VertexGenerator.h"
#include "GHEP/GHepParticle.h"
#include "LlewellynSmith/LwlynSmithQELCCPXSec.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
LwlynSmithQELCCPXSec::LwlynSmithQELCCPXSec() :
XSecAlgorithmI("genie::LwlynSmithQELCCPXSec")
{

}
//____________________________________________________________________________
LwlynSmithQELCCPXSec::LwlynSmithQELCCPXSec(string config) :
XSecAlgorithmI("genie::LwlynSmithQELCCPXSec", config)
{

}
//____________________________________________________________________________
LwlynSmithQELCCPXSec::~LwlynSmithQELCCPXSec()
{

}
//____________________________________________________________________________
double LwlynSmithQELCCPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) {LOG("LwlynSmith",pWARN) << "not a valid process"; return 0.;}
  if(! this -> ValidKinematics (interaction) ) {LOG("LwlynSmith",pWARN) << "not valid kinematics"; return 0.;}

  if (kps == kPSTnctnBnctl){
    return this->FullDifferentialXSec(interaction);
  }

  // Get kinematics & init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  double E  = init_state.ProbeE(kRfHitNucRest);
  double E2 = TMath::Power(E,2);
  double ml = interaction->FSPrimLepton()->Mass();
  double M  = target.HitNucMass();
  double q2 = kinematics.q2();

  // One of the xsec terms changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = (is_neutrino) ? -1 : 1;

  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);    

  double F1V   = fFormFactors.F1V();
  double xiF2V = fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
  double Fp    = fFormFactors.Fp();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("LwlynSmith", pDEBUG) << "\n" << fFormFactors;
#endif

  // Calculate auxiliary parameters
  double ml2     = TMath::Power(ml,    2);
  double M2      = TMath::Power(M,     2);
  double M4      = TMath::Power(M2,    2);
  double FA2     = TMath::Power(FA,    2);
  double Fp2     = TMath::Power(Fp,    2);
  double F1V2    = TMath::Power(F1V,   2);
  double xiF2V2  = TMath::Power(xiF2V, 2);
  double Gfactor = M2*kGF2*fCos8c2 / (8*kPi*E2);
  double s_u     = 4*E*M + q2 - ml2;
  double q2_M2   = q2/M2;

  // Compute free nucleon differential cross section
  double A = (0.25*(ml2-q2)/M2) * (
	      (4-q2_M2)*FA2 - (4+q2_M2)*F1V2 - q2_M2*xiF2V2*(1+0.25*q2_M2)
              -4*q2_M2*F1V*xiF2V - (ml2/M2)*( 
               (F1V2+xiF2V2+2*F1V*xiF2V)+(FA2+4*Fp2+4*FA*Fp)+(q2_M2-4)*Fp2));
  double B = -1 * q2_M2 * FA*(F1V+xiF2V);
  double C = 0.25*(FA2 + F1V2 - 0.25*q2_M2*xiF2V2);

  double xsec = Gfactor * (A + sign*B*s_u/M2 + C*s_u*s_u/M4);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("LwlynSmith", pDEBUG)
     << "dXSec[QEL]/dQ2 [FreeN](E = "<< E << ", Q2 = "<< -q2 << ") = "<< xsec;
  LOG("LwlynSmith", pDEBUG) 
                 << "A(Q2) = " << A << ", B(Q2) = " << B << ", C(Q2) = " << C;
#endif

  //----- The algorithm computes dxsec/dQ2
  //      Check whether variable tranformation is needed
  if(kps!=kPSQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("LwlynSmith", pDEBUG)
     << "Jacobian for transformation to: " 
                  << KinePhaseSpace::AsString(kps) << ", J = " << J;
#endif
    xsec *= J;
  }

  //----- if requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- compute nuclear suppression factor
  //      (R(Q2) is adapted from NeuGEN - see comments therein)
  double R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);

  //----- number of scattering centers in the target
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("LwlynSmith", pDEBUG) 
       << "Nuclear suppression factor R(Q2) = " << R << ", NNucl = " << NNucl;
#endif

  xsec *= (R*NNucl); // nuclear xsec

  return xsec;
}
//____________________________________________________________________________
double LwlynSmithQELCCPXSec::FullDifferentialXSec(const Interaction *  interaction)const{

  // First we need access to all of the particles in the interaction
  // The particles were stored in the lab frame
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  const TLorentzVector leptonMom = kinematics.FSLeptonP4();
  const TLorentzVector outNucleonMom = kinematics.HadSystP4();

  TLorentzVector * neutrinoMom = init_state.GetProbeP4(kRfLab);
  TLorentzVector * inNucleonMom = init_state.TgtPtr()->HitNucP4Ptr();

  // Now we calculate q and qTilde
  //TLorentzVector qP4(0,0,0,0);
  TLorentzVector qTildeP4(0,0,0,0);
  //qP4 = *neutrinoMom - leptonMom;
  //qTildeP4 = outNucleonMom - *inNucleonMom;

  qTildeP4 = *neutrinoMom- leptonMom; // TESTING: Use q rather than qtilde
  
  double Q2tilde = -1 * qTildeP4.Mag2(); 
  interaction->KinePtr()->SetQ2(Q2tilde);

//  LOG("LwlynSmith",pDEBUG) << "Q2tilde = " << Q2tilde;
//  LOG("LwlynSmith",pDEBUG) << "Q2 (not tilde)= " << -1 * qP4.Mag2();
//  LOG("LwlynSmith",pDEBUG) << "Q2 difference (tilde - not) = " << Q2tilde + qP4.Mag2();

  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);

  double F1V   = fFormFactors.F1V();
  double xiF2V = fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
  double Fp    = fFormFactors.Fp();

  double Gfactor = kGF2*fCos8c2 / (8*kPi*kPi*inNucleonMom->E()*neutrinoMom->E()*outNucleonMom.E()*leptonMom.E());

  // Now, we can calculate the cross section
  double tau = Q2tilde / (4 * inNucleonMom->Mag2());
  double h1 = FA*FA*(1 + tau) + tau*(F1V + xiF2V)*(F1V + xiF2V);
  double h2 = FA*FA + F1V*F1V + tau*xiF2V*xiF2V;
  double h3 = 2.0 * FA * (F1V + xiF2V);
  double h4 = 0.25 * xiF2V*xiF2V *(1-tau) + 0.5*F1V*xiF2V + FA*Fp - tau*Fp*Fp;
  
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = (is_neutrino) ? -1 : 1;
  double l1 = 2*neutrinoMom->Dot(leptonMom)*(inNucleonMom->Mag2());
  double l2 = 2*(neutrinoMom->Dot(*inNucleonMom)) * (inNucleonMom->Dot(leptonMom)) - neutrinoMom->Dot(leptonMom)*inNucleonMom->Mag2();
  double l3 = (neutrinoMom->Dot(*inNucleonMom) * qTildeP4.Dot(leptonMom)) - (neutrinoMom->Dot(qTildeP4) * leptonMom.Dot(*inNucleonMom));
  l3 *= sign;
  double l4 = neutrinoMom->Dot(leptonMom) * qTildeP4.Dot(qTildeP4) - 2*neutrinoMom->Dot(qTildeP4)*leptonMom.Dot(qTildeP4);
  double l5 = neutrinoMom->Dot(*inNucleonMom) * leptonMom.Dot(qTildeP4) + leptonMom.Dot(*inNucleonMom)*neutrinoMom->Dot(qTildeP4) - neutrinoMom->Dot(leptonMom)*inNucleonMom->Dot(qTildeP4);

  double LH = 2 *(l1*h1 + l2*h2 + l3*h3 + l4*h4 + l5*h2);

  delete neutrinoMom;
  
  double xsec = Gfactor * LH;

  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- compute nuclear suppression factor
  //      (R(Q2) is adapted from NeuGEN - see comments therein)
  double R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);
  // LOG("LwlynSmith",pINFO)  << "Nuclear Suppression Factor = " << R;

  //----- number of scattering centers in the target
  const Target & target = init_state.Tgt();
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("LwlynSmith", pDEBUG) 
    << "Nuclear suppression factor R(Q2) = " << R << ", NNucl = " << NNucl;
#endif

  xsec *= (R*NNucl); // nuclear xsec

  return xsec;

}
//____________________________________________________________________________
double LwlynSmithQELCCPXSec::Integral(const Interaction * in) const
{
  bool nuclear_target = in->InitState().Tgt().IsNucleus();
  if(!nuclear_target || !fDoAvgOverNucleonMomentum) {
    return fXSecIntegrator->Integrate(this,in);
  }

  double E = in->InitState().ProbeE(kRfHitNucRest);
  if(fLFG || E < fEnergyCutOff) {
    // clone the input interaction so as to tweak the
    // hit nucleon 4-momentum in the averaging loop
    Interaction in_curr(*in);

    // hit target
    Target * tgt = in_curr.InitState().TgtPtr();

    // get nuclear masses (init & final state nucleus)
    int nucleon_pdgc = tgt->HitNucPdg();
    bool is_p = pdg::IsProton(nucleon_pdgc);
    int Zi = tgt->Z();
    int Ai = tgt->A();
    int Zf = (is_p) ? Zi-1 : Zi;
    int Af = Ai-1;
    PDGLibrary * pdglib = PDGLibrary::Instance();
    TParticlePDG * nucl_i = pdglib->Find( pdg::IonPdgCode(Ai, Zi) );
    TParticlePDG * nucl_f = pdglib->Find( pdg::IonPdgCode(Af, Zf) );
    if(!nucl_f) {
      LOG("LwlynSmith", pFATAL)
	<< "Unknwown nuclear target! No target with code: "
	<< pdg::IonPdgCode(Af, Zf) << " in PDGLibrary!";
      exit(1);
    }
    double Mi  = nucl_i -> Mass(); // initial nucleus mass
    double Mf  = nucl_f -> Mass(); // remnant nucleus mass

    // throw nucleons with fermi momenta and binding energies 
    // generated according to the current nuclear model for the
    // input target and average the cross section
    double xsec_sum = 0.;
    const int nnuc = 2000;
    // VertexGenerator for generating a position before generating
    // each nucleon
    VertexGenerator * vg = new VertexGenerator();
    vg->Configure("Default");
    for(int inuc=0;inuc<nnuc;inuc++){
      // Generate a position in the nucleus
      TVector3 nucpos = vg->GenerateVertex(&in_curr,tgt->A());
      tgt->SetHitNucPosition(nucpos.Mag());

      // Generate a nucleon
      fNuclModel->GenerateNucleon(*tgt, nucpos.Mag());
      TVector3 p3N = fNuclModel->Momentum3();
      double   EN  = Mi - TMath::Sqrt(p3N.Mag2() + Mf*Mf);
      TLorentzVector* p4N = tgt->HitNucP4Ptr();
      p4N->SetPx (p3N.Px());
      p4N->SetPy (p3N.Py());
      p4N->SetPz (p3N.Pz());
      p4N->SetE  (EN);

      double xsec = fXSecIntegrator->Integrate(this,&in_curr);
      xsec_sum += xsec;
    }
    double xsec_avg = xsec_sum / nnuc;
    delete vg;
    return xsec_avg;
  }else{
    return fXSecIntegrator->Integrate(this,in);
  }
}
//____________________________________________________________________________
bool LwlynSmithQELCCPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsQuasiElastic()) return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  bool prcok = proc_info.IsWeakCC() && ((isP&&isnub) || (isN&&isnu));
  if(!prcok) return false;

  return true;
}
//____________________________________________________________________________
void LwlynSmithQELCCPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LwlynSmithQELCCPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LwlynSmithQELCCPXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
  
  double thc = fConfig->GetDoubleDef(
                              "CabbiboAngle", gc->GetDouble("CabbiboAngle"));
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);

   // load QEL form factors model
  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
                                             this->SubAlg("FormFactorsAlg"));
  assert(fFormFactorsModel);
  fFormFactors.SetModel(fFormFactorsModel); // <-- attach algorithm

   // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  // Get nuclear model for use in Integral()
  RgKey nuclkey = "IntegralNuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);

  fLFG = fNuclModel->ModelType(Target()) == kNucmLocalFermiGas;

  // Always average over initial nucleons if the nuclear model is LFG
  fDoAvgOverNucleonMomentum =
    fLFG || fConfig->GetBoolDef("IntegralAverageOverNucleonMomentum", false);

  fEnergyCutOff = 0.;

  if(fDoAvgOverNucleonMomentum) {
    // Get averaging cutoff energy
    fEnergyCutOff = 
      fConfig->GetDoubleDef("IntegralNuclearInfluenceCutoffEnergy", 2.0);
  }
}
//____________________________________________________________________________
