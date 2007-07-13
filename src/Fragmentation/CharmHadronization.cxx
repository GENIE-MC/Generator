//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - August 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TPythia6.h>
#include <TVector3.h>
#include <TF1.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Fragmentation/CharmHadronization.h"
#include "Fragmentation/FragmentationFunctionI.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/FragmRecUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;

extern "C" void py1ent_(int *,  int *, double *, double *, double *);
extern "C" void py2ent_(int *,  int *, int *, double *);

//____________________________________________________________________________
CharmHadronization::CharmHadronization() :
HadronizationModelI("genie::CharmHadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
CharmHadronization::CharmHadronization(string config) :
HadronizationModelI("genie::CharmHadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
CharmHadronization::~CharmHadronization()
{

}
//____________________________________________________________________________
void CharmHadronization::Initialize(void) const
{
  fPythia = TPythia6::Instance();
}
//____________________________________________________________________________
TClonesArray * CharmHadronization::Hadronize(
                                        const Interaction * interaction) const
{
  LOG("CharmHad", pNOTICE) << "Running CHARM hadronizer";

  PDGLibrary * pdglib = PDGLibrary::Instance();
  RandomGen *  rnd    = RandomGen::Instance();
  
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  const TLorentzVector & p4Had = kinematics.HadSystP4();

  double Ev = init_state.ProbeE(kRfLab);
  double W  = kinematics.W(true);
  double Eh = p4Had.Energy();

  LOG("CharmHad", pNOTICE) << "Ehad (LAB) = " << Eh << ", W = " << W;

  //-- Generate a charmed hadron PDG code and generate its energy
  //   based on the input fragmentation function

  int nupdg = init_state.ProbePdg();
  int cpdgc = 0;

  double mC=0., mC2=0., EC=0., EC2=0., z=0., pC2=0.;

  unsigned int itry = 0;
  bool found = false;

  while(itry++<1000 && !found) {

    cpdgc = this->GenerateCharmHadron(nupdg,Ev); // generate hadron
    z     = fFragmFunc->GenerateZ();             // generate z(=Eh/Ev)
    mC    = pdglib->Find(cpdgc)->Mass();         // lookup mass
    EC    = z*Eh;                                // charm hadron energy
    mC2   = TMath::Power(mC,2);
    EC2   = TMath::Power(EC,2);
    pC2   = EC2-mC2;

    found = (mC<W && pC2>0);
  }
  
  if(!found){
     LOG("CharmHad", pWARN) << "Couldn't generate charm hadron id!";
     return 0;
  }

  LOG("CharmHad", pNOTICE) 
         << "Generated charm hadron = " << cpdgc << "(m = " << mC << ")";
  LOG("CharmHad", pNOTICE) 
         << "Generated charm hadron z = " << z << ", E = " << EC;

  //-- Generate the charm hadron pT^2 and pL^2 (with respect to the
  //   hadronic system direction @ the LAB)

  double plC2=0., ptC2=0.;
  found = false;
  while(!found) {
    ptC2 = fCharmPT2pdf->GetRandom();   
    plC2 = EC2 - ptC2 - mC2;
    found = (plC2>0);
  }

  if(!found){
     LOG("CharmHad", pWARN) << "Couldn't generate charm hadron 4p!";
     return 0;
  }

  //-- Generate the charm hadron momentum components (with respect to
  //   the hadronic system direction @ the LAB)
  double ptC = TMath::Sqrt(ptC2);
  double plC = TMath::Sqrt(plC2);   
  double phi = (2*kPi) * rnd->RndHadro().Rndm();
  double pxC = ptC * TMath::Cos(phi);
  double pyC = ptC * TMath::Sin(phi);
  double pzC = plC; 

  LOG("CharmHad", pNOTICE) 
           << "Generated charm hadron pT (tranv to pHad) = " << ptC;
  LOG("CharmHad", pNOTICE) 
           << "Generated charm hadron pL (along to pHad) = " << plC;

  //-- Rotate charm hadron 3-momentum from the system with z' along
  //   the hadronic momentum to the LAB

  TVector3 unitvq = p4Had.Vect().Unit();

  TVector3 p3C(pxC, pyC, pzC);
  p3C.RotateUz(unitvq);

  //-- Boost charm hadron 4-momentum from the LAB to the HCM frame

  TLorentzVector p4C(p3C,EC);

  LOG("CharmHad", pNOTICE) 
    << "Charm hadron p4 (LAB) = " << utils::print::P4AsString(&p4C);

  TVector3 beta = -1 * p4Had.BoostVector();
  p4C.Boost(beta);

  LOG("CharmHad", pNOTICE) 
    << "Charm hadron p4 (HCM) = " << utils::print::P4AsString(&p4C);

  //-- Hadronic non-charm remnant 4p at HCM
  TLorentzVector p4H(0,0,0,W);
  TLorentzVector p4R = p4H - p4C;

  double WR = p4R.M();

  LOG("CharmHad", pNOTICE) 
             << "Hadronic-blob (remnant) invariant mass = " << WR;

  //-- If the remnant mass is too low then return a null record so
  //   as to re-try.
  if(WR < kNucleonMass + kPionMass) {
     LOG("CharmHad", pWARN) << "Too small hadronic remnant mass!";
     return 0;
  }

  //-- Check whether I was only asked to generate the charm hadron and the
  //   hadronic blob and not to hadronize the blob as well
  if(fCharmOnly) {
    // Create particle list (fragmentation record)
    TClonesArray * particle_list = new TClonesArray("TMCParticle", 2);
    particle_list->SetOwner(true);

    // insert the generated charm hadron & the hadronic (non-charm) blob

    new ((*particle_list)[0]) TMCParticle (1,cpdgc,
	     -1,-1,-1, p4C.Px(),p4C.Py(),p4C.Pz(),p4C.E(),mC, 0,0,0,0,0);
    new ((*particle_list)[1]) TMCParticle (1,kPdgHadronicBlob,
             -1,-1,-1, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(),WR, 0,0,0,0,0);

    return particle_list;
  }
  
  //  *****  Hadronize non-charm hadronic blob *****

  //-- Figure out the quark systems to input to PYTHIA based on simple
  //   quark model arguments

  int  nuc  = target.HitNucPdg();
  //bool qpdg = target.HitQrkPdg();
  bool sea  = target.HitSeaQrk();
  bool isP  = pdg::IsProton (nuc);

  int  qrkSyst1 = 0;
  int  qrkSyst2 = 0;

  // scattering off valence quark
  if(!sea) {
     // Scattering off a valence d quark. The remaining target diquark is:
     // - uu for protons 
     // - ud for neutrons (50% in singlet, 50% in triplet state)
     if(isP) qrkSyst2 = kPdgUUDiquarkS1;
     else {
       qrkSyst2 = kPdgUDDiquarkS1;
       double Rqq = rnd->RndHadro().Rndm();
       if(Rqq<0.5) qrkSyst2 = kPdgUDDiquarkS0;
     }

     // The hit d quark is now a c quark within the charm hadron. 
     // The charm hadron was formed by picking up quarks from qqbar pairs.
     // Check the charm meson quark contents to figure out which quark
     // was left outside the meson so that it can go on hadronizing...
     
     if      (cpdgc == kPdgD0      ) qrkSyst1 = kPdgUQuark; // D0 (c ubar),      
     else if (cpdgc == kPdgDP      ) qrkSyst1 = kPdgDQuark; // D^+ (c dbar)
     else if (cpdgc == kPdgDPs     ) qrkSyst1 = kPdgSQuark; // Ds^+ (c sbar)
     else if (cpdgc == kPdgLambdaPc) qrkSyst1 = -2101;      // Lamda_c^+ (c ud)
  }
  else {
    LOG("CharmHad", pWARN) 
           << "I can't handle hadronization for scattering off sea quarks";
    return 0;
  }

  //-- Run PYTHIA.
  //   The hadronization takes place at the *remnant hadrons* centre of mass
  //   frame (not the hadronic centre of mass frame). Remember to boost.
  int ip     = 0;
  int mstj21 = fPythia->GetMSTJ(21);

  fPythia->SetMSTJ(21,0);                  // inhibit decays at this stage
  py2ent_(&ip, &qrkSyst1, &qrkSyst2, &WR); // hadronize
  fPythia->SetMSTJ(21,mstj21);             // restore

  //-- Get PYTHIA's LUJETS event record
  TClonesArray * remnants = 0;
  fPythia->GetPrimaries();
  remnants = dynamic_cast<TClonesArray *>(fPythia->ImportParticles("All"));
  if(!remnants) {
     LOG("CharmHad", pWARN) << "Couldn't hadronize (non-charm) remnants!";
     return 0;
  }

  //-- Create particle list (fragmentation record)

  unsigned int nremnants = remnants->GetEntries();
  unsigned int nhadrons  = 2 + nremnants;

  TClonesArray * particle_list = new TClonesArray("TMCParticle", nhadrons);
  particle_list->SetOwner(true);

  //-- Insert the generated charm hadron & the hadronic (non-charm) blob.
  //   In this case the hadronic blob is entered as a pre-fragm. state.

  new ((*particle_list)[0]) TMCParticle (1,cpdgc,
                  -1,-1,-1,  p4C.Px(),p4C.Py(),p4C.Pz(),p4C.E(),mC, 0,0,0,0,0);
  new ((*particle_list)[1]) TMCParticle (11,kPdgHadronicBlob,
       -1,2,3, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(),WR, 0,0,0,0,0);

  //-- Velocity for boosting fragments from the remnant hadrons CM frame
  //   to the hadronic -all hadrons- CM frame
  TVector3 rmnbeta = -1 * p4R.BoostVector();

  //-- Boost all hadronic blob fragments to the hadronic CM, fix their
  //   mother/daughter assignments and add them to the fragmentation record.

  TMCParticle * remnant = 0; 
  TIter remnant_iter(remnants);

  unsigned int rpos=2; // offset in event record

  while( (remnant = (TMCParticle *) remnant_iter.Next()) ) {

    // insert and get a pointer to inserted object for mods
    TMCParticle * boosted_remnant = 
      new ((*particle_list)[rpos++]) TMCParticle (*remnant);

    // boost 
    double px = remnant->GetPx();
    double py = remnant->GetPy();
    double pz = remnant->GetPz();
    double E  = remnant->GetEnergy();

    TLorentzVector p4(px,py,pz,E);
    p4.Boost(rmnbeta);

    boosted_remnant -> SetPx     (p4.Px());
    boosted_remnant -> SetPy     (p4.Py());
    boosted_remnant -> SetPz     (p4.Pz());
    boosted_remnant -> SetEnergy (p4.E() );

    // handle insertion of charmed hadron 
    int ip  = boosted_remnant->GetParent();
    int ifc = boosted_remnant->GetFirstChild();
    int ilc = boosted_remnant->GetLastChild();

    boosted_remnant -> SetParent     ( (ip  == 0 ?  1 : ip +1) );
    boosted_remnant -> SetFirstChild ( (ifc == 0 ? -1 : ifc+1) );
    boosted_remnant -> SetLastChild  ( (ilc == 0 ? -1 : ilc+1) );
  }

  //-- Print & return the fragmentation record
  //utils::fragmrec::Print(particle_list);
  return particle_list;
}
//____________________________________________________________________________
double CharmHadronization::Weight(void) const
{
  return 1.; // does not generate weighted events
}
//____________________________________________________________________________
PDGCodeList * CharmHadronization::SelectParticles(
                               const Interaction * /*interaction*/) const
{
  return 0;
}
//____________________________________________________________________________
TH1D * CharmHadronization::MultiplicityProb(
           const Interaction * /*interaction*/, Option_t * /*opt*/)  const
{
  return 0;
}
//____________________________________________________________________________
int CharmHadronization::GenerateCharmHadron(int nupdg, double EvLab) const
{
  // generate a charmed hadron pdg code using a charm fraction table

  RandomGen * rnd = RandomGen::Instance();
  double r  = rnd->RndHadro().Rndm();

  if(pdg::IsNeutrino(nupdg)) {
     double tf = 0;
     if      (r < (tf+=fD0FracSpl->Evaluate(EvLab)))  return kPdgD0;       // D^0
     else if (r < (tf+=fDpFracSpl->Evaluate(EvLab)))  return kPdgDP;       // D^+
     else if (r < (tf+=fDsFracSpl->Evaluate(EvLab)))  return kPdgDPs;      // Ds^+
     else                                             return kPdgLambdaPc; // Lamda_c^+

  } else if(pdg::IsAntiNeutrino(nupdg)) {
     if      (r < fD0BarFrac)          return kPdgAntiD0;
     else if (r < fD0BarFrac+fDmFrac)  return kPdgDM;
     else                              return kPdgDMs;
  }

  LOG("CharmHad", pERROR) << "Could not generate a charm hadron!";
  return 0;
}
//____________________________________________________________________________
void CharmHadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CharmHadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CharmHadronization::LoadConfig(void)
{
  //AlgConfigPool * confp = AlgConfigPool::Instance();
  //const Registry * gc = confp->GlobalParameterList();

  fCharmOnly = ! fConfig->GetBoolDef("HadronizeRemnants", true);

  //-- Get a fragmentation function
  fFragmFunc = dynamic_cast<const FragmentationFunctionI *> (
    this->SubAlg("FragmentationFunc"));
  assert(fFragmFunc);

  fCharmPT2pdf = new TF1("fCharmPT2pdf", "exp(-0.213362-6.62464*x)",0,0.6);

  // neutrino charm fractions: D^0, D^+, Ds^+ (remainder: Lamda_c^+)
  //
  const int nc = 15;
  double ec[nc] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,80.,90.,100.};

  double d0frac[nc] = { .000, .320, .460,  .500,  .520,  .530, .540, .540,
                        .540, .550,  .550,  .560,  .570, .580, .600  };
  double dpfrac[nc] = { .000, .120, .180,  .200,  .200,  .210, .210, .210,
                        .210, .210,  .220,  .220,  .220, .230, .230  };
  double dsfrac[nc] = { .000, .054, .078,  .130,  .130,  .140, .140, .140,
                        .140, .140,  .140,  .140,  .140, .150, .150  };
 
  fD0FracSpl = new Spline(nc, ec, d0frac);
  fDpFracSpl = new Spline(nc, ec, dpfrac);
  fDsFracSpl = new Spline(nc, ec, dsfrac);

  // anti-neutrino charm fractions: bar(D^0), D^-, (remainder: Ds^-)
  //
  fD0BarFrac = 0.667;
  fDmFrac    = 0.222;
}
//____________________________________________________________________________
