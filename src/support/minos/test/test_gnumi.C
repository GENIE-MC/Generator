#include <string>

void test_gnumi()
{
  // need to have called load_genie_into_root.C

  // currently TFile not TChain, so better resolve to a single file
  string fluxpatt = "$FLUXPATH/v19/fluka05_le010z185i/root/fluka05_le010z185i_1.root";
  string fluxfilename(gSystem->ExpandPathName(fluxpatt.c_str()));

  genie::flux::GNuMIFlux* gnumi = new genie::flux::GNuMIFlux();
  gnumi->LoadBeamSimData(fluxfilename);

  const genie::PDGCodeList& pdgList = gnumi->FluxParticles();
  int nfluxtypes = pdgList.size();
  cout << "GNuMIFlux supplies " << nfluxtypes << " PDG particles of type: ";
  // don't try iterators .. apparently they aren't usable in CINT
  for (int indx=0; indx < nfluxtypes ; ++indx) 
    cout << pdgList[indx] << " ";
  cout << endl;

  cout << "MaxEnergy " << gnumi->MaxEnergy() << endl;

  cout << "scan for max weight at NearDet" << endl;
  gnumi->UseFluxAtNearDetCenter();
  gnumi->ScanForMaxWeight();

  bool genok = gnumi->GenerateNext();
  bool atend = gnumi->End();
  cout << "gnumi->GenerateNext() returned " << (genok?"true":"false") 
       << " End() returns " << (atend?"true":"false")
       << endl;

  cout << "Current entry: "
       << " PdgCode=" << gnumi->PdgCode()
       << " Weight=" << gnumi->Weight()
       << endl;
  gnumi->PrintCurrent();


}
