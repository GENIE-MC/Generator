//____________________________________________________________________________
/*!

\program gtestQENuclSuppression

\brief   

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 20, 2004

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>

#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/NuclearUtils.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
  TNtuple * nt = new TNtuple("nt","","Z:A:nucl:Q2:R");

  Interaction * ipD2   = Interaction::QELNC(kPdgTgtDeuterium, kPdgProton,  kPdgNuMu, 0.);
  Interaction * inD2   = Interaction::QELNC(kPdgTgtDeuterium, kPdgNeutron, kPdgNuMu, 0.);
  Interaction * ipC12  = Interaction::QELNC(kPdgTgtC12,       kPdgProton,  kPdgNuMu, 0.);
  Interaction * inC12  = Interaction::QELNC(kPdgTgtC12,       kPdgNeutron, kPdgNuMu, 0.);
  Interaction * ipFe56 = Interaction::QELNC(kPdgTgtFe56,      kPdgProton,  kPdgNuMu, 0.);
  Interaction * inFe56 = Interaction::QELNC(kPdgTgtFe56,      kPdgNeutron, kPdgNuMu, 0.);

  const int nin = 6;
  Interaction * in[nin] = { ipD2, inD2, ipC12, inC12, ipFe56, inFe56 };

  const int    N        = 3001;
  const double Q2min    = 0.000001;
  const double Q2max    = 10;
  const double logQ2min = TMath::Log(Q2min);
  const double logQ2max = TMath::Log(Q2max);
  const double dlogQ2   = (logQ2max-logQ2min)/(N-1);

  for(int j=0; j<nin; j++) {
    Interaction * interaction = in[j];

    for(int i=0; i<N; i++) {
       double Q2 = TMath::Exp(logQ2min + i*dlogQ2); 
       interaction->KinePtr()->SetQ2(Q2);         
       double R = utils::nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);
       LOG("test", pNOTICE) << interaction->InitState().Tgt().AsString() << " ==> Q2 = " << Q2 << " GeV, R = " << R;

       nt->Fill(
        interaction->InitState().Tgt().Z(),
        interaction->InitState().Tgt().A(),
        interaction->InitState().Tgt().HitNucPdg(),
        Q2, R);
    }
  }

  TFile f("./nuclsupp.root","recreate");
  nt->Write();
  f.Close();

  return 0;
}
