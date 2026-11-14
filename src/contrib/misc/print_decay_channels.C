//
// Print-out the decay channel information for the input particle
//
// eg. To print-out the tau- (pdg code = 15) decay channels,
// type:
// root[0] .L print_decay_channels.C
// root[1] print_decay_channels(15)
//
// Costas Andreopoulos <c.andreopoulos \at cern.ch>
// University of Liverpool
//
#include <iostream>
#include <sstream>
using namespace std;

// root includes
#include "TSystem.h"
#include "TParticlePDG.h"
#include "TDecayChannel.h"
#include "TDatabasePDG.h"

void print_decay_channels(int pdgc)
{
  gSystem->Load("libEG");

  //cout << " *** Printing-out decay channels for: " << pdgc << endl;

  ostringstream pdg_table;
  pdg_table << gSystem->Getenv("GENIE")
            << "/data/evgen/catalogues/pdg/genie_pdg_table.txt";

  TDatabasePDG * pdglib = TDatabasePDG::Instance();
  pdglib->ReadPDGTable(pdg_table.str().c_str());

  TParticlePDG * p = pdglib->GetParticle(pdgc);

  cout << " *** Printing-out decay channels for: " << p->GetName() << endl;

  double brtot=0;
  for(int j=0; j<p->NDecayChannels(); j++) {
    cout << "\t - decay channel id = " << j << ", channel = "
         << p->GetName() << " --> ";
    TDecayChannel * dch = p->DecayChannel(j);
    for(int k=0; k<dch->NDaughters(); k++) {
      cout << pdglib->GetParticle(dch->DaughterPdgCode(k))->GetName();
      if ( k < dch->NDaughters() - 1 ) cout << " + ";
    } // k = daughters
    cout << ", BR = " << dch->BranchingRatio() << endl;
    brtot += dch->BranchingRatio();
  } // j = channels

  cout << "Sum{BR} = " << brtot << endl;
}
