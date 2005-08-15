//____________________________________________________________________________
/*!

\class   genie::NtpGHepAnalyzer

\brief   An object that knows how to look at the GHEP event record and extract
         summary information for the generated event so as to fill the output
         ntuple's NtpMCSummary.

         This information also exists at the Interaction summary attached to
         the event record, but the attached Interaction is NULL for the GHEP
         records not generated within the GENIE framework but extracted &
         translated from external generators.
         Being able to analyze the GHEP record allows us to properly generate
         the GENIE ntuple for other generators too facilitating cross-generator
         comparisons.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 10, 2005

*/
//____________________________________________________________________________

#include <cassert>

#include <TLorentzVector.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "EventGeneration/EventRecord.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpGHepAnalyzer.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
NtpGHepAnalyzer::NtpGHepAnalyzer()
{
  this->Init();
}
//____________________________________________________________________________
NtpGHepAnalyzer::~NtpGHepAnalyzer()
{

}
//____________________________________________________________________________
void NtpGHepAnalyzer::AnalyzeEventRecord(const EventRecord & evrec)
{
  LOG("NtpGHepAnalyzer", pINFO) << "Analyzing input GHEP event record";

  this->CleanUp();

  fEventRec = &evrec;

  GHepParticle * probep = fEventRec->GetParticle(0);
  assert(probep);
  GHepParticle * targetp = fEventRec->GetParticle(1);
  assert(targetp);
  GHepParticle * fslp = fEventRec->GetParticle(probep->FirstDaughter());
  assert(fslp);

  fProbePdgC = probep  -> PdgCode();
  fFslPdgC   = fslp    -> PdgCode();
  fTgtPdgC   = targetp -> PdgCode();

  LOG("NtpGHepAnalyzer", pINFO)
           << "PDG Codes: probe = " << fProbePdgC
                     << ", target = " << fTgtPdgC << ", fsl = " << fFslPdgC;

  bool tgt_is_nucleus = targetp->IsNucleus();
  bool tgt_is_nucleon = pdg::IsProton(fTgtPdgC) || pdg::IsNeutron(fTgtPdgC);

  GHepParticle * hitnuclp = 0;
  if(tgt_is_nucleus) {
   hitnuclp = fEventRec->GetParticle(2);
   if(!hitnuclp->Status() == kIstNucleonTarget) hitnuclp = 0;
  }
  if(tgt_is_nucleon) hitnuclp = targetp;

  fNuclPdgC = (hitnuclp) ? hitnuclp->PdgCode() : 0;

  LOG("NtpGHepAnalyzer", pINFO)
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

  fNProton  = this->NEntries(kPdgProton);
  fNNeutron = this->NEntries(kPdgNeutron);
  fNPi0     = this->NEntries(kPdgPi0);
  fNPiPlus  = this->NEntries(kPdgPiPlus);
  fNPiMinus = this->NEntries(kPdgPiMinus);
  fNK0      = this->NEntries(0);
  fNKPlus   = this->NEntries(0);
  fNKMinus  = this->NEntries(0);

  LOG("NtpGHepAnalyzer", pINFO)
              << "N (p, n) = (" << fNProton << ", " << fNNeutron << ")";
  LOG("NtpGHepAnalyzer", pINFO)
        << "N (pi+, pi0, pi-) = (" << fNPiPlus << ", "
                                   << fNPi0 << ", " << fNPiMinus << ")";
  LOG("NtpGHepAnalyzer", pINFO)
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
  TIter piter(fEventRec);
  while ( (p = (GHepParticle *) piter.Next()) ) {
     bool is_res_pdg = res_utils::IsBaryonResonance(p->PdgCode());
     bool is_res_Ist = p->Status() == kIstPreDecayResonantState;
     if( is_res_Ist && is_res_pdg ) fScatType = kScResonant;
  }
  LOG("NtpGHepAnalyzer", pINFO)
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
  LOG("NtpGHepAnalyzer", pINFO)
      << "Interaction Type = " << InteractionType::AsString(fProcType);
  assert(fProcType != kIntNull);

  // Cross section for this interaction at the selected energy and
  // Cross section for this interaction for the selected kinematics
  // note: this information is only stoted at the summary
  Interaction * interaction = fEventRec->GetInteraction();
  if(interaction) {
     fXSec  = interaction->XSec();
     fdXSec = interaction->DiffXSec();
  }
  LOG("NtpGHepAnalyzer", pINFO) << "XSec =      " << fXSec;
  LOG("NtpGHepAnalyzer", pINFO) << "Diff.XSec = " << fdXSec;

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

  LOG("NtpGHepAnalyzer", pINFO)
              << "Kine.Vars: v = " << fNu << ", x = " << fX << ", y = " << fY;
  LOG("NtpGHepAnalyzer", pINFO)
              << "Kine.Vars: Q2 = " << fQ2 << ", W = " << fW;

}
//____________________________________________________________________________
void NtpGHepAnalyzer::CleanUp(void)
{
  LOG("NtpGHepAnalyzer", pINFO) << "Cleaning up last GHEP Analyzer results";

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
void NtpGHepAnalyzer::Init(void)
{
  LOG("NtpGHepAnalyzer", pINFO) << "Initializing GHEP Analyzer results";

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
  fXSec      = 0;
  fdXSec     = 0;
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
unsigned int NtpGHepAnalyzer::NEntries(int pdgc) const
{
  unsigned int nentries = 0;

  TIter piter(fEventRec);
  GHepParticle * p = 0;

  while( (p = (GHepParticle *) piter.Next()) ) {
    if(p->PdgCode()==pdgc && p->Status()==kIStStableFinalState) nentries++;
  }
  return nentries;
}
//____________________________________________________________________________
