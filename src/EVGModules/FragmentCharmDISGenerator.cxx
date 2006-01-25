//____________________________________________________________________________
/*!

\class   genie::FragmentCharmDISGenerator

\brief   Generates the final state hadronic system in v DIS charm production
         interactions. The charm hadron is generated according to the input
         charm fractions, its longitudinal momentum is generated from the
         input fragmentation function, its transverse momentum is generated
         from an exponential function for the input pT^2 scale. The hadronic
         remnants are generated so as to conserve the hadronic shower charge
         and are distributed uniformly in the remaining phase space.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created June 22, 2005

*/
//____________________________________________________________________________

#include <TGenPhaseSpace.h>
#include <TVector3.h>
#include <TF1.h>

#include "Conventions/Constants.h"
#include "EVGModules/FragmentCharmDISGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepOrder.h"
#include "Fragmentation/FragmentationFunctionI.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::print;

//___________________________________________________________________________
FragmentCharmDISGenerator::FragmentCharmDISGenerator() :
HadronicSystemGenerator("genie::FragmentCharmDISGenerator")
{

}
//___________________________________________________________________________
FragmentCharmDISGenerator::FragmentCharmDISGenerator(string config):
HadronicSystemGenerator("genie::FragmentCharmDISGenerator", config)
{

}
//___________________________________________________________________________
FragmentCharmDISGenerator::~FragmentCharmDISGenerator()
{

}
//___________________________________________________________________________
void FragmentCharmDISGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- If the struck nucleon was within a nucleus, then add the final state
  //   nucleus at the EventRecord
  this->AddTargetNucleusRemnant(evrec);

  //-- Add the charm hadron & the hadronic remnants

  bool charm_only = fConfig->Exists("model-charm-only") ?
                               fConfig->GetBool("model-charm-only") : false;

  if(charm_only) this->GenerateCharmHadronOnly(evrec, true);
  else           this->GenerateHadronicSystem(evrec);
}
//___________________________________________________________________________
bool FragmentCharmDISGenerator::GenerateCharmHadronOnly(
                              GHepRecord * evrec, bool ignore_remnants) const
{
// Generates the charmed hadron as described in the class documentation.
// The remnant hadronic system is added as a single 'void' entry (rootino)

  TLorentzVector P4Had = this->HadronicSystemP4(evrec);
  LOG("FragmentCharm", pINFO)
                     << "P4Had [LAB] = " << utils::print::P4AsString(&P4Had);

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction -> GetInitialState();

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double EHad = P4Had.Energy();
  double MHad = P4Had.M();

  //-- Get a fragmentation function
  const FragmentationFunctionI * fragmf =
     dynamic_cast<const FragmentationFunctionI *> (this->SubAlg(
        "fragmentation-func-alg-name", "fragmentation-func-param-set"));

  //-- Generate a charmed hadron PDG code
  int itry = 0, pdgc = 0;
  double m = 0, z = 0, EC = 0, m2 = 0, EC2 = 0, p2 = 0;

  while(itry<1000 && pdgc==0) {

    pdgc = this->CharmedHadronPdgCode(E); // generate hadron type
    z    = fragmf->GenerateZ(); // generate hadron fract. energy
    m    = PDGLibrary::Instance()->Find(pdgc)->Mass();
    EC   = z*EHad; // charm hadron energy
    m2   = TMath::Power(m, 2);
    EC2  = TMath::Power(EC,2);
    p2   = EC2-m2;

    if(m>MHad || p2<0) pdgc=0;
    itry++;
  }
  assert(pdgc!=0);

  //-- set charm hadron pdg code in the interaction summary
  XclsTag xcls;
  xcls.SetCharm(pdgc);
  interaction->SetExclusiveTag(xcls);

  LOG("FragmentCharm", pINFO)
    << "Generated: charm hadron pdg = " << pdgc << " (m = " << m << ")";
  LOG("FragmentCharm", pINFO)
       << "Generated: z = " << z << ", charm hadron E = " << EC
                                          << " / hadronic E = " << EHad;
  //-- Generate a charmed hadron pT
  double pT2 = this->GeneratePT2(p2);
  double pT  = TMath::Sqrt(pT2);

  LOG("FragmentCharm", pINFO) << "Maximum p2 = " << p2;
  LOG("FragmentCharm", pINFO) << "Generated: charm hadron pT = " << pT;

  //-- Compute charm hadron px, py, pz
  double t   = 2.*kPi*RandomGen::Instance()->Random1().Rndm();
  double pxC = pT * TMath::Sin(t);
  double pyC = pT * TMath::Cos(t);
  double pzC = TMath::Sqrt(EC2-pT2-m2);
  //double p   = TMath::Sqrt(p2);

  LOG("FragmentCharm", pINFO)
          << "Generated: charm hadron (px,py,pz) = ("
                           << pxC << ", " << pyC << ", " << pzC << ")";

  //-- Rotate charm hadron p along the direction of the hadronic shower
  //
  //??
  TVector3 P3Had(P4Had.Px(), P4Had.Py(), P4Had.Pz());
  TVector3 unitvec = P3Had.Unit(); // unit vector along the hadronic p
  TVector3 P3C(pxC, pyC, pzC);
  P3C.RotateUz(unitvec);
  pxC = P3C.Px();
  pyC = P3C.Py();
  pzC = P3C.Pz();

  SLOG("FragmentCharm", pINFO)
          << "Rotated: charm hadron (px,py,pz) = ("
                           << pxC << ", " << pyC << ", " << pzC << ")";

  //-- 4-p of the remaining hadronic system (remnants)
  double pRx = P4Had.Px()     - pxC;
  double pRy = P4Had.Py()     - pyC;
  double pRz = P4Had.Pz()     - pzC;
  double ER  = P4Had.Energy() - EC;

  //-- if selected not to ignore remnants at least check that there is
  //   sufficient mass to hadronize enough pions to conserve charge
  if(!ignore_remnants) {
    int qhs = this->HadronShowerCharge(interaction);
    int qch   = int (PDGLibrary::Instance()->Find(pdgc)->Charge() / 3.);
    int qremn = qhs - qch;

    const unsigned int nmult = TMath::Max(2,qremn);

    double mneed   = nmult * kPionMass;
    double mneed2  = TMath::Power(mneed, 2);
    double mavail2 = ER*ER - pRx*pRx - pRy*pRy - pRz*pRz;

    if(mavail2 < mneed2) {
       LOG("FragmentCharm", pINFO)
            << "(Available mass)^2 = " << mavail2
               << " < (mass needed)^2 = " << mneed2 << " - Retrying";
       return false;
    }
  }

  //-- add the entries at the event record

  int mom = GHepOrder::StruckNucleonPosition(interaction);

  evrec->AddParticle(
         pdgc, kIStStableFinalState, mom,-1,-1,-1, pxC,pyC,pzC,EC, 0,0,0,0);
  evrec->AddParticle(
            0, kIStStableFinalState, mom,-1,-1,-1, pRx,pRy,pRz,ER, 0,0,0,0);

  return true;
}
//___________________________________________________________________________
void FragmentCharmDISGenerator::GenerateHadronicSystem(
                                                   GHepRecord * evrec) const
{
// Do not use just yet - not tested.
//
  Interaction * interaction = evrec->GetInteraction();

  // generate the charm hadron
  int itry = 0;
  while (! this->GenerateCharmHadronOnly(evrec, false) ) {
    itry++;
    assert(itry<1000);
  }

  GHepParticle * ch = evrec->GetParticle(5);
  int pdgc = ch->PdgCode();

  // replace the 'void' remnant hadronic system with something realistic

  //-- compute the charge of the remnant system so that charge is conserved
  int qhs = this->HadronShowerCharge(interaction);
  int qch   = int (PDGLibrary::Instance()->Find(pdgc)->Charge() / 3.);
  int qremn = qhs - qch;

  LOG("FragmentCharm", pINFO)
         << "Hadron Charge (Shower, Charm, Remnants) = ("
                          << qhs << ", " << qch << ", " << qremn << ")";

  //-- find the hadronic remnants 'fake' particle and get its 4-p
  GHepParticle * remnants = evrec->FindParticle(0, kIStStableFinalState, 0);
  TLorentzVector * pR4 = remnants->GetP4();
  LOG("FragmentCharm", pINFO)
      << "P4(" << remnants->Name()
               << ") [Remnants/LAB] = " << utils::print::P4AsString(pR4);

  //-- Hadronic remnants multiplicity (= all - charm hadron)
/*
  double a = 0;
  double b = 1.3;
  double W2 = TMath::Power(remnants->Mass(),2);
  int    hmult = int( a+b*log(W2) );

  const unsigned int nmult = TMath::Max(
                                TMath::Max(2,TMath::Abs(qremn)), hmult);
*/
  const unsigned int nmult = TMath::Max(2,qremn);

  LOG("FragmentCharm", pINFO) << "Remnant multiplicity = " << nmult;

  //-- remnant hadronic system spectrum [Here I should be using a scheme
  //   for 'recombination' of the remaining quark 'lines']

  int rpdgc[nmult];
  for(unsigned int i=0; i<nmult; i++) {
    if      (qremn<0) { rpdgc[i] = kPdgPiMinus; qremn++; }
    else if (qremn>0) { rpdgc[i] = kPdgPiPlus;  qremn--; }
    else              { rpdgc[i] = kPdgPi0;              }
  }
  double mass[nmult];
  for(unsigned int i=0; i<nmult; i++) {
    mass[i] = PDGLibrary::Instance()->Find(rpdgc[i])->Mass();
  }

  //-- try to generate momenta for the remnant hadrons

  LOG("FragmentCharm", pINFO) << "Generating phase space";

  TGenPhaseSpace phase_space_generator;
  bool permitted = phase_space_generator.SetDecay(*pR4, nmult, mass);

  if(!permitted) {
     LOG("FragmentCharm", pERROR)
                           << "*** Decay forbidden by kinematics! ***";
     assert(false);
  } else {
     //-- generate kinematics in the Center-of-Mass (CM) frame
     phase_space_generator.Generate();

     int mom = GHepOrder::StruckNucleonPosition(interaction);

     //-- insert final state products into a TClonesArray of TMCParticles
     //   and return it
     for(unsigned int i = 0; i < nmult; i++) {

       TLorentzVector * p4 = phase_space_generator.GetDecay(i);

       double phx = p4->Px();
       double phy = p4->Py();
       double phz = p4->Pz();
       double Eh  = p4->Energy();

       LOG("FragmentCharm", pINFO)
               << "Adding final state particle PDGC = " << rpdgc[i]
               << " with 4-P = " << utils::print::P4AsString(p4);
       evrec->AddParticle(rpdgc[i], kIStStableFinalState,
                                   mom,-1,-1,-1, phx,phy,phz,Eh, 0,0,0,0);
    }//add remnants

    //-- change the 'rootino' status from kIStStableFinalState to
    //   kIstDISPreFragmHadronicState
    remnants->SetStatus(kIstDISPreFragmHadronicState);

  } //permitted decay

  delete pR4;
}
//___________________________________________________________________________
TLorentzVector FragmentCharmDISGenerator::HadronicSystemP4(
                                                    GHepRecord * evrec) const
{
// Hadronic System 4P in LAB

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  //incoming v:
  TLorentzVector * nu_p4 = init_state.GetProbeP4(kRfLab);
  assert(nu_p4);

  //struck nucleon:
  TLorentzVector * nucl_p4 = init_state.GetTarget().StruckNucleonP4();
  assert(nucl_p4);

  //final state primary lepton:
  int fsl_pdgc = interaction->GetFSPrimaryLepton()->PdgCode();
  GHepParticle * fsl = evrec->FindParticle(fsl_pdgc,kIStStableFinalState,0);
  assert(fsl);

  TLorentzVector * fsl_p4 = fsl->GetP4();

  LOG("HadronicVertex", pINFO)
                 << "\n v [LAB]: " << P4AsString( nu_p4   )
                 << "\n N [LAB]: " << P4AsString( nucl_p4 )
                 << "\n l [LAB]: " << P4AsString( fsl_p4  );

  //-- Compute the velocity of the LAB frame in the Final State Hadronic
  //   CM Frame (Pv + Pnuc = Pfsl + Sum{Phad_{i}})

  double PX = nu_p4->Px()     + nucl_p4->Px()     - fsl_p4->Px();
  double PY = nu_p4->Py()     + nucl_p4->Py()     - fsl_p4->Py();
  double PZ = nu_p4->Pz()     + nucl_p4->Pz()     - fsl_p4->Pz();
  double E  = nu_p4->Energy() + nucl_p4->Energy() - fsl_p4->Energy();

  delete fsl_p4;
  delete nu_p4;

  assert(E>0);
  TLorentzVector P4( PX, PY, PZ, E );

  return P4;
}
//___________________________________________________________________________
int FragmentCharmDISGenerator::CharmedHadronPdgCode(double E) const
{
  // generate a charmed hadron pdg code using a charm fraction table

  //------ eventually the charm fraction table should be an external object
  //       that gets loaded from its XML configuration file -----------------/

  RandomGen * rnd = RandomGen::Instance();

  double rndm = rnd->Random1().Rndm();

  if(E <= 20) {
     if      (rndm <= 0.32)                  return  421; // D^0
     else if (rndm >  0.32 && rndm <= 0.37)  return  411; // D^+
     else if (rndm >  0.37 && rndm <= 0.55)  return  431; // Ds^+
     else                                    return 4122; // Lamda_c^+

  } else if (E > 20 && E <=40)  {
     if      (rndm <= 0.50)                  return  421;
     else if (rndm >  0.50 && rndm <= 0.60)  return  411;
     else if (rndm >  0.60 && rndm <= 0.82)  return  431;
     else                                    return 4122;

  } else {
     if      (rndm <= 0.64)                  return  421;
     else if (rndm >  0.64 && rndm <= 0.86)  return  411;
     else if (rndm >  0.86 && rndm <= 0.95)  return  431;
     else                                    return 4122;
  }
  return 421;
}
//____________________________________________________________________________
double FragmentCharmDISGenerator::GeneratePT2(double pT2max) const
{
// Generate a charmed hadron pT

  double pT2scale = fConfig->Exists("pT2scale") ?
                                         fConfig->GetDouble("pT2scale") : 0.6;
  TF1 expo("expo","TMath::Exp(-x/[0])",0,TMath::Sqrt(TMath::Abs(pT2max)));
  expo.SetParameter(0,pT2scale);

  double pT2 = 9999999;

  while(pT2 > pT2max) pT2 = expo.GetRandom();

  return pT2;
}
//____________________________________________________________________________
int FragmentCharmDISGenerator::HadronShowerCharge(
                                        const Interaction * interaction) const
{
// Returns the hadron shower charge in units of +e
// eg in v n -> l- X the hadron shower charge is +1

  int HadronShowerCharge = 0;

  const InitialState & init_state = interaction->GetInitialState();

  int hit_nucleon = init_state.GetTarget().StruckNucleonPDGCode();

  assert( pdg::IsProton(hit_nucleon) || pdg::IsNeutron(hit_nucleon) );

  double qfsl  = interaction->GetFSPrimaryLepton()->Charge() / 3.;
  double qinit = PDGLibrary::Instance()->Find(hit_nucleon)->Charge() / 3.;

  HadronShowerCharge = (int) (qinit - qfsl);

  return HadronShowerCharge;
}
//____________________________________________________________________________


