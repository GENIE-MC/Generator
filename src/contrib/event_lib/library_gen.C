#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include <iostream>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"


using namespace genie ;

// this macro is supposed to be an example on how to create an event library starting from GENIE files
// by any means it is supposed to be a complete and refined product

int library_gen( TString in_file_name, int tgt_pdg, int probe_pdg, bool CC, TString out_file_name = "" ) {


  if ( out_file_name == "" ) {
    out_file_name = in_file_name ;
    out_file_name.ReplaceAll( ".ghep.root",
                              ".evtlib.root" ) ;
  }


  // retrieving the genie tree 
  TFile in_file( in_file_name ) ;

  NtpMCTreeHeader * header = dynamic_cast<NtpMCTreeHeader*> (in_file.Get("header"));

  // Get the GENIE GHEP tree and set its branch address
  TTree * in_tree = dynamic_cast<TTree*> (in_file.Get("gtree"));
  NtpMCEventRecord * mcrec = 0;
  in_tree->SetBranchAddress("gmcrec", & mcrec);

  //create the library tree
  TFile out_file(out_file_name, "UPDATE");

  // create the proper structure
  TParticlePDG * pdg_tgt   = PDGLibrary::Instance() -> Find( tgt_pdg ) ; 
  TParticlePDG * pdg_probe = PDGLibrary::Instance() -> Find( probe_pdg ) ; 
  
  std::string current_name = CC ? "cc" : "nc" ;
  std::string probe_name = pdg_probe -> GetTitle() ;
  if ( ! CC ) {
    probe_name = probe_pdg > 0 ? "nu" : "nu_bar" ; 
  }
  
  TDirectory* out_dir = fout.mkdir( pdg_tgt -> GetTitle() ) 
    ->mkdir( current_name.c_str() )
    ->mkdir( probe_name.c_str() ) ;
  out_dir -> cd();

  //create out tree
  TTree trout("records", "records");
  float Enu_out, weight_out;
  trout.Branch("Enu",    &Enu_out);
  trout.Branch("weight", &weight_out);
  int prod_id;
  trout.Branch("prod_id", &prod_id);
  int nparts_out = 0;
  int pdg_out[1000];
  float E_out[1000];
  float px_out[1000];
  float py_out[1000];
  float pz_out[1000];
  trout.Branch("nparts", &nparts_out);
  trout.Branch("pdg", &pdg_out, "pdg[nparts]/I");
  trout.Branch("E", &E_out, "E[nparts]/F");
  trout.Branch("px", &px_out, "px[nparts]/F");
  trout.Branch("py", &py_out, "py[nparts]/F");
  trout.Branch("pz", &pz_out, "pz[nparts]/F");
  trin->SetBranchAddress("evType", &prod_id);
    
  
  // Event loop
  for(Long64_t i=0; i < in_tree->GetEntries(); i++) {

    in_tree->GetEntry(i);

    EventRecord & event = *(mcrec->event);

    const Interaction & inter = *( event.Summary() ) ;
    const InitState   & init_state = inter.InitState() ;
    const ProcessInfo & proc_info = inter.ProcInfo() ;

    // check the target
    if ( init_state.TgtPdg() == tgt_pdg ) {

      // check cc or nc
      if ( ( CC   && ProcessInfo.IsWeakCC() ) ||
	   ( !CC  && ProcessInfo.IsWeakNC() ) ) {
	
	// cehck the neutrino flavour
	if ( ( CC && InitState.ProbePdg() == probe_pdg ) ||
	     ( !CC && 
	       pdg::IsNeutrino( InitState.ProbePdg() ) && 
	       InitState.ProbePdg() * probe_pdg > 0 ) ) {

	  
	  // particle loop
	  const GHepParticle * part = nullptr ;
	  TObjArrayIter piter(event);
	  while( (part = (const GHepParticle *) piter.Next()) ) {
	    
	    if ( part -> Status() != kIStStableFinalState ) continue ;
	    // grep only stable particles
	    
	    // store them

	  } // particle loop
	  
	}  //good probe
      }  // proper current
    }  // right target

    mcrec->Clear() ;
  } // event loop

}



void usage()
{
  std::cout << "Usage: convert_gibuu_root_to_evt_lib NUC PROBE CURR IN.root OUT.root" << std::endl;
  std::cout << "  NUC - eg C12" << std::endl;
  std::cout << "  PROBE - eg nu_mu_bar" << std::endl;
  std::cout << "  CURR - cc or nc" << std::endl;
  std::cout << std::endl;
  std::cout << "For nc valid values of PROBE are 'nu' or 'nu_bar'" << std::endl;
  exit(2);
}
int main(int argc, char** argv)
{
  if(argc != 6) usage();
  
  const std::string nuc = argv[1];
  const std::string probe = argv[2];
  const std::string current = argv[3];
  const std::string inname = argv[4];
  const std::string outname = argv[5];
  TFile fin(inname.c_str());
  if(fin.IsZombie()) abort();
  TTree* trin = (TTree*)fin.Get("RootTuple");
  assert(trin);
  if(trin->GetEntries() == 0){
    std::cout << "RootTuple in " << inname << " is empty. Skipping." << std::endl;
    return 1;
  }
  if(current != "nc" && current != "cc") usage();
  if(current == "cc" &&
     (probe != "nu_e"   && probe != "nu_e_bar"  &&
      probe != "nu_mu"  && probe != "nu_mu_bar" &&
      probe != "nu_tau" && probe != "nu_tau_bar")) usage();
  if(current == "nc" &&
     (probe != "nu"     && probe != "nu_bar")) usage();
  // Take the nucleus on faith
  // NC
  int lep_pdg = 0;
  if(probe == "nu_e"      ) lep_pdg = +11;
  if(probe == "nu_e_bar"  ) lep_pdg = -11;
  if(probe == "nu_mu"     ) lep_pdg = +13;
  if(probe == "nu_mu_bar" ) lep_pdg = -13;
  if(probe == "nu_tau"    ) lep_pdg = +15;
  if(probe == "nu_tau_bar") lep_pdg = -15;
  // NC - use nu_mu by convention
  if(probe == "nu"    ) lep_pdg = +14;
  if(probe == "nu_bar") lep_pdg = -14;
  assert(lep_pdg != 0);
  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dirout = fout.mkdir(nuc.c_str())->mkdir(current.c_str())->mkdir(probe.c_str());
  dirout->cd();
  TTree trout("records", "records");
  float Enu_out, weight_out;
  trout.Branch("Enu", &Enu_out);
  trout.Branch("weight", &weight_out);
  int prod_id;
  trout.Branch("prod_id", &prod_id);
  int nparts_out = 0;
  int pdg_out[1000];
  float E_out[1000];
  float px_out[1000];
  float py_out[1000];
  float pz_out[1000];
  trout.Branch("nparts", &nparts_out);
  trout.Branch("pdg", &pdg_out, "pdg[nparts]/I");
  trout.Branch("E", &E_out, "E[nparts]/F");
  trout.Branch("px", &px_out, "px[nparts]/F");
  trout.Branch("py", &py_out, "py[nparts]/F");
  trout.Branch("pz", &pz_out, "pz[nparts]/F");
  trin->SetBranchAddress("evType", &prod_id);
  double Enu_in;
  trin->SetBranchAddress("lepIn_E", &Enu_in);
  double weight_in;
  trin->SetBranchAddress("weight", &weight_in);
  double nu_px, nu_py;
  trin->SetBranchAddress("lepIn_Px", &nu_px);
  trin->SetBranchAddress("lepIn_Py", &nu_py);
  double lepE_in, lep_px_in, lep_py_in, lep_pz_in;
  trin->SetBranchAddress("lepOut_E", &lepE_in);
  trin->SetBranchAddress("lepOut_Px", &lep_px_in);
  trin->SetBranchAddress("lepOut_Py", &lep_py_in);
  trin->SetBranchAddress("lepOut_Pz", &lep_pz_in);
  std::vector<double> *E_in = 0, *px_in = 0, *py_in = 0, *pz_in = 0;
  trin->SetBranchAddress("E", &E_in);
  trin->SetBranchAddress("Px", &px_in);
  trin->SetBranchAddress("Py", &py_in);
  trin->SetBranchAddress("Pz", &pz_in);
  std::vector<int>* pdg_in = 0;
  trin->SetBranchAddress("barcode", &pdg_in);
  // There is a bug in GiBUU generating NC antineutrino events which we can
  // work around here
  const bool ncbar_fix = (current == "nc" && probe == "nu_bar");
  for(long i = 0; i < trin->GetEntries(); ++i){
    trin->GetEntry(i);
    // These two modes should not be produced for NC - skip those events
    if(ncbar_fix && (prod_id == 32 || prod_id == 33)) continue;
    assert(nu_px == 0 && nu_py == 0); // Otherwise we'll need some rotations
    weight_out = weight_in;
    Enu_out = Enu_in;
    nparts_out = E_in->size()+1;
    assert(nparts_out <= 1000);
    E_out[0] = lepE_in;
    px_out[0] = lep_px_in;
    py_out[0] = lep_py_in;
    pz_out[0] = lep_pz_in;
    pdg_out[0] = lep_pdg;
    for(int j = 0; j < nparts_out; ++j){
      E_out[j+1] = (*E_in)[j];
      px_out[j+1] = (*px_in)[j];
      py_out[j+1] = (*py_in)[j];
      pz_out[j+1] = (*pz_in)[j];
      pdg_out[j+1] = (*pdg_in)[j];
    }
    trout.Fill();
  }
  dirout->cd();
  trout.Write();
  TVectorD n(1);
  n[0] = 1;
  n.Write("nfiles");
  std::cout << "Wrote " << trout.GetEntries() << " events to " << outname << std::endl;
  return 0;
}
