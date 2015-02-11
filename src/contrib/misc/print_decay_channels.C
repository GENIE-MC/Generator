//
// Print-out the decay channel information for the input particle
//
// eg. To print-out the tau- (pdg code = 15) decay channels, 
// type:
// root[0] .L print_decay_channels.C
// root[1] print_decay_channels(15)
//
// Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
// University of Liverpool & STFC Rutherford Appleton Lab
//

gSystem->Load("libEG");

void print_decay_channels(int pdgc)
{
  ostringstream pdg_table;
  pdg_table << gSystem->Getenv("GENIE") << "/data/pdg/genie_pdg_table.txt";

  TDatabasePDG * pdglib = TDatabasePDG::Instance();
  pdglib->ReadPDGTable(pdg_table.str().c_str());

  TParticlePDG * p = pdglib->GetParticle(pdgc);

  cout << " *** Printing-out decay channels for: " << p->GetName() << endl;

  double brtot=0;
  for(int j=0; j<p->NDecayChannels(); j++) {
        cout << "\t - decay channel id = " << j << ", channel = " << p->GetName() << " --> ";
	TDecayChannel * dch = p->DecayChannel(j);
        for(int k=0; k<dch->NDaughters(); k++) {
	   cout << pdglib->GetParticle(dch->DaughterPdgCode(k))->GetName();
	   if(k < dch->NDaughters() - 1) cout << " + ";
	}//k
	cout << ", BR = " << dch->BranchingRatio() << endl;
        brtot += dch->BranchingRatio();
  }//j

  cout << "Sum{BR} = " << brtot << endl;
}
