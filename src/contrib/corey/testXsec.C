// shell% genie
// genie[0] .x test_xsec.C
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"

#include "Algorithm/AlgConfigPool.h"
#include "EVGCore/InteractionList.h"
#include "VLE/StrumiaVissaniIBDPXSec.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGUtils.h"
#include "Messenger/Messenger.h"

#endif

#include "PDG/PDGCodes.h"
#include "Conventions/Units.h"
#include "Interaction/InteractionType.h"
#include "Interaction/ScatteringType.h"

const Int_t nknots = 44;
const Double_t E[nknots] = {
   /*1.806,*/ 2.01, 2.25, 2.51, 2.80,
   3.12,  3.48, 3.89, 4.33, 4.84,
   5.40,  6.02, 6.72, 7.49, 8.36,
   8.83,  9.85, 11.0, 12.3, 13.7,
   15.3,  17.0, 19.0, 21.2, 23.6,
   26.4,  29.4, 32.8, 36.6, 40.9,
   43.2,  48.2, 53.7, 59.9, 66.9,
   74.6,  83.2, 92.9, 104,  116,
   129,   144,  160,  179,  200
};
const Double_t s[nknots] = {
   /*0,*/  0.00351, 0.00735, 0.0127, 0.0202,
   0.0304, 0.0440,  0.06190, 0.0854, 0.116,
   0.155,  0.205,   0.269,   0.349,  0.451,
   0.511,  0.654,   0.832,   1.05,   1.33,
   1.67,   2.09,    2.61,    3.24,   4.01,
   4.95,   6.08,    7.44,    9.08,   11.0,
   12.1,   14.7,    17.6,    21.0,   25.0,
   29.6,   34.8,    40.7,    47.3,   54.6,
   62.7,   71.5,    81.0,    91.3,   102.
};
Double_t xsec[nknots], xsec_n[nknots];

using namespace genie;

TFile* outf;
TNtuple* nt;

void testXsec(const Char_t* outfn="VLExsecNT.root") {
   gROOT->Macro("loadlibs.C");
   
/*   
   PDGCodeList * targets = new PDGCodeList;
   targets->push_back(kPdgTgtFreeP);
   
   PDGCodeList * neutrinos = new PDGCodeList;
   neutrinos->push_back(kPdgAntiNuE);
*/
   
   // load VLE xsec parameters
   AlgConfigPool * pool = AlgConfigPool::Instance();
   const Registry * config = 
      pool->FindRegistry("genie::StrumiaVissaniIBDPXSec/Default");
   
   // our xsec calculator
   XSecAlgorithmI * alg = new StrumiaVissaniIBDPXSec;
   alg->Configure(*config);
   
   // prepare to save results
   outf = TFile::Open(outfn,"recreate");
   nt   = new TNtuple("nt","nt","Ev:xsec:xsnun:xspaper");
   nt->SetDirectory(outf);
   
   
   // nu_e_bar + p --> n + e+
   InitialState init_state(kPdgTgtFreeP, kPdgAntiNuE);
   ProcessInfo  proc_info(genie::kScInverseBetaDecay, genie::kIntWeakCC);
   genie::Interaction * interaction = 
      new genie::Interaction(init_state, proc_info);
   Target * target  = interaction->InitStatePtr()->TgtPtr();
   target->SetHitNucPdg(kPdgProton);
   Int_t i=0;
   for (i = 0; i < nknots; i++) {
      TLorentzVector p4(0,0,E[i]*1e-3,E[i]*1e-3);
      interaction->InitStatePtr()->SetProbeP4(p4);
      xsec[i] = alg->Integral(interaction);
   }
   
   
   // nu_e + n --> p + e-
   InitialState init_state_n(kPdgTgtFreeN, kPdgNuE);
   ProcessInfo  proc_info_n(genie::kScInverseBetaDecay, genie::kIntWeakCC);
   genie::Interaction * interaction_n = 
      new genie::Interaction(init_state_n, proc_info_n);
   Target * target_n  = interaction_n->InitStatePtr()->TgtPtr();
   target_n->SetHitNucPdg(kPdgNeutron);
   for (i = 0; i < nknots; i++) {
      TLorentzVector n4(0,0,E[i]*1e-3,E[i]*1e-3);
      interaction_n->InitStatePtr()->SetProbeP4(n4);
      xsec_n[i] = alg->Integral(interaction_n);
   }
   
   for (i = 0; i < nknots; i++) {
      Printf("xsec(E=%g MeV) = %g \t xsec_n = %g", E[i],
         (1E+41/units::cm2)*xsec[i],
         (1E+41/units::cm2)*xsec_n[i]);
      nt->Fill(E[i],
              (1E+41/units::cm2)*xsec[i],
              (1E+41/units::cm2)*xsec_n[i],
              s[i]);
   }
      
   outf->Write();
   
}
