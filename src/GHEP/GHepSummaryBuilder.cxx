//____________________________________________________________________________
/*!

\class   genie::GHepSummaryBuilder

\brief   An object that knows how to look at the GHEP event record and extract
         summary information.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 10, 2005

*/
//____________________________________________________________________________

#include <cassert>

#include <TLorentzVector.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepSummaryBuilder.h"
#include "GHEP/GHepOrder.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GHepSummaryBuilder::GHepSummaryBuilder()
{
  this->Init();
}
//____________________________________________________________________________
GHepSummaryBuilder::~GHepSummaryBuilder()
{

}
//____________________________________________________________________________
void GHepSummaryBuilder::AnalyzeEventRecord(const GHepRecord & evrec)
{
  LOG("GHepSummaryBuilder", pINFO) << "Analyzing input GHEP event record";

  this->CleanUp();

  GHepParticle * probep = evrec.GetParticle(GHepOrder::ProbePosition());
  assert(probep);
  GHepParticle * targetp = evrec.GetParticle(1);
  assert(targetp);
  GHepParticle * fslp = evrec.GetParticle(probep->FirstDaughter());
  assert(fslp);

  fProbePdgC = probep  -> PdgCode();
  fFslPdgC   = fslp    -> PdgCode();
  fTgtPdgC   = targetp -> PdgCode();

  LOG("GHepSummaryBuilder", pINFO)
           << "PDG Codes: probe = " << fProbePdgC
                     << ", target = " << fTgtPdgC << ", fsl = " << fFslPdgC;

  bool tgt_is_nucleus = targetp->IsNucleus();
  bool tgt_is_nucleon = pdg::IsProton(fTgtPdgC) || pdg::IsNeutron(fTgtPdgC);

  GHepParticle * hitnuclp = 0;
  if(tgt_is_nucleus) {
   hitnuclp = evrec.GetParticle(2);
   if(!hitnuclp->Status() == kIstNucleonTarget) hitnuclp = 0;
  }
  if(tgt_is_nucleon) hitnuclp = targetp;

  fNuclPdgC = (hitnuclp) ? hitnuclp->PdgCode() : 0;

  LOG("GHepSummaryBuilder", pINFO)
                              << "PDG Codes: hit nucleon = " << fNuclPdgC;
  fIQrkPdgC  = 0;
  fFQrkPdgC  = 0;
  fResPdgC   = 0;
  fChHadPdgC = 0;

  if(tgt_is_nucleus || tgt_is_nucleon) {
    fTgtZ = pdg::IonPdgCodeToZ(fTgtPdgC);
    fTgtA = pdg::IonPdgCodeToA(fTgtPdgC);;
    fTgtN = fTgtA - fTgtZ;
  }

  fNProton  = evrec.NEntries(kPdgProton,  kIStStableFinalState);
  fNNeutron = evrec.NEntries(kPdgNeutron, kIStStableFinalState);
  fNPi0     = evrec.NEntries(kPdgPi0,     kIStStableFinalState);
  fNPiPlus  = evrec.NEntries(kPdgPiPlus,  kIStStableFinalState);
  fNPiMinus = evrec.NEntries(kPdgPiMinus, kIStStableFinalState);
  fNK0      = evrec.NEntries(0,           kIStStableFinalState);
  fNKPlus   = evrec.NEntries(0,           kIStStableFinalState);
  fNKMinus  = evrec.NEntries(0,           kIStStableFinalState);

  LOG("GHepSummaryBuilder", pINFO)
              << "N (p, n) = (" << fNProton << ", " << fNNeutron << ")";
  LOG("GHepSummaryBuilder", pINFO)
        << "N (pi+, pi0, pi-) = (" << fNPiPlus << ", "
                                   << fNPi0 << ", " << fNPiMinus << ")";
  LOG("GHepSummaryBuilder", pINFO)
        << "N (K+, K0, K-) = (" << fNKPlus << ", "
                                     << fNK0 << ", " << fNKMinus << ")";

  // Find the scattering type (QEL,DIS,RES,COH) using the following 'logic':
  // - has nuclear target && no struck nucleon                   : COH event
  // - has nuclear target && no struck nucleon + e in init state : IMD event
  // - has a baryon res. with status = kIstPreDecayResonantState : RES event
  // - final state (primary) hadronic system has mult = 1        : QEL event
  // - final state (primary) hadronic system has mult >= 2       : DIS event

  //COH:
  if(tgt_is_nucleus) { if(!hitnuclp) fScatType = kScCoherent; }

  //IMD
  GHepParticle * elec = evrec.GetParticle(2);
  if(elec->Status()==kIStInitialState &&
              elec->PdgCode()==kPdgElectron) fScatType = kScInverseMuDecay;

  //QEL or DIS
  if(fScatType != kScCoherent && fScatType != kScInverseMuDecay) {
    assert(hitnuclp);
    int dau1 = hitnuclp->FirstDaughter();
    int dau2 = hitnuclp->LastDaughter();
    assert( dau1 != -1 && dau2 != -1 && dau2 >= dau1);
    int multiplicity = 1 + dau2 - dau1;
    if(multiplicity == 1) fScatType = kScQuasiElastic;
    else                  fScatType = kScDeepInelastic;
  }
  //RES:
  GHepParticle * p = 0;
  TIter piter(&evrec);
  while ( (p = (GHepParticle *) piter.Next()) ) {
     bool is_res_pdg = utils::res::IsBaryonResonance(p->PdgCode());
     bool is_res_Ist = p->Status() == kIstPreDecayResonantState;
     if( is_res_Ist && is_res_pdg ) fScatType = kScResonant;
  }
  LOG("GHepSummaryBuilder", pINFO)
         << "Scattering Type = " << ScatteringType::AsString(fScatType);
  assert(fScatType != kScNull);

  // find the interaction type (CC,NC,...) by checking for charge
  // change (neutrino - final state primary lepton)
  int dQ = int(probep->Charge() - fslp->Charge());
  if (dQ != 0) fProcType = kIntWeakCC;
  else {
    // no charge change (either a weak NC or E/M)
    int  pdgc = probep->PdgCode();
    bool isnu = pdg::IsNeutrino(pdgc) || pdg::IsAntiNeutrino(pdgc);
    if (isnu) fProcType = kIntWeakNC;
    else      fProcType = kIntEM;
  }
  LOG("GHepSummaryBuilder", pINFO)
      << "Interaction Type = " << InteractionType::AsString(fProcType);
  assert(fProcType != kIntNull);


  fVtx      = probep->GetV4();
  fProbe4P  = probep->GetP4();
  fNucl4P   = (hitnuclp) ? hitnuclp->GetP4() : new TLorentzVector(0,0,0,0);
  fFsl4P    = fslp->GetP4();
  f2Fsl4P   = new TLorentzVector(0,0,0,0);
  fHadShw4P = new TLorentzVector(0,0,0,0);

  // q 4P

  TLorentzVector * k1 = probep->GetP4();
  TLorentzVector * k2 = fslp->GetP4();
  (*k1) -= (*k2);
  delete k2;
  fq4p = k1;

  // nu
  fNu = (*fq4p)*(*fNucl4P) / kNucleonMass;

  // Q2
  fQ2 = -1. * fq4p->M2();

  // x
  fX = 0.5*fQ2/(kNucleonMass*fNu);

  // Inelasticity, y = qP/kP
  float qP = (*fq4p)*(*fNucl4P);
  float kP = (*fProbe4P)*(*fNucl4P);
  fY = (kP>0) ? qP/kP : -1;

  LOG("GHepSummaryBuilder", pINFO)
              << "Kine.Vars: v = " << fNu << ", x = " << fX << ", y = " << fY;
  LOG("GHepSummaryBuilder", pINFO)
              << "Kine.Vars: Q2 = " << fQ2 << ", W = " << fW;

}
//____________________________________________________________________________
void GHepSummaryBuilder::CleanUp(void)
{
  LOG("GHepSummaryBuilder", pINFO) << "Cleaning up last GHEP Analyzer results";

  if (fVtx     ) delete fVtx;
  if (fProbe4P ) delete fProbe4P;
  if (fNucl4P  ) delete fNucl4P;
  if (fFsl4P   ) delete fFsl4P;
  if (f2Fsl4P  ) delete f2Fsl4P;
  if (fHadShw4P) delete fHadShw4P;
  if (fq4p     ) delete fq4p;

  this->Init();
}
//____________________________________________________________________________
void GHepSummaryBuilder::Init(void)
{
  LOG("GHepSummaryBuilder", pINFO) << "Initializing GHEP Analyzer results";

  fProbePdgC = 0;
  fFslPdgC   = 0;
  fTgtPdgC   = 0;
  fNuclPdgC  = 0;
  fIQrkPdgC  = 0;
  fFQrkPdgC  = 0;
  fResPdgC   = 0;
  fChHadPdgC = 0;
  fTgtZ      = 0;
  fTgtA      = 0;
  fTgtN      = 0;
  fScatType  = kScNull;
  fProcType  = kIntNull;
  fNu        = 0;
  fX         = 0;
  fY         = 0;
  fZ         = 0;
  fQ2        = 0;
  fW         = 0;
  fEmFrac    = 0;
  fVtx       = 0;
  fProbe4P   = 0;
  fNucl4P    = 0;
  fFsl4P     = 0;
  f2Fsl4P    = 0;
  fHadShw4P  = 0;
  fq4p       = 0;
  fNProton   = 0;
  fNNeutron  = 0;
  fNPi0      = 0;
  fNPiPlus   = 0;
  fNPiMinus  = 0;
  fNK0       = 0;
  fNKPlus    = 0;
  fNKMinus   = 0;
}
//____________________________________________________________________________
