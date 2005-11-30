//____________________________________________________________________________
/*!

\program gmctest

\brief   A simple program that reads a generated event sample and plots some
         basic MC truth quantities [does not do much yet]

         Syntax :
           gmctest -f filename

         Options:
           -f specifies the GENIE/ROOT file with the generated event sample

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 02, 2005
*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TPostScript.h>
#include <TCanvas.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCSummary.h"
#include "Messenger/Messenger.h"

using std::string;
using namespace genie;

string getRootFilename    (int argc, char ** argv);
bool   checkRootFilename  (string filename);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- scan the command line arguments and get the ROOT filename
  string filename = getRootFilename(argc,argv);
  if ( !checkRootFilename(filename) ) return 1;

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile inpfile(filename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( inpfile.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( inpfile.Get("header") );

  LOG("Main", pINFO) << "Input tree header: " << *thdr;

  //-- figure out the TTree format (GENIE supports multiple formats)
  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord);

  //-- set the event record branch ptr
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- define a "summary" tree
  TTree * st = new TTree("st","");

  int   st_iev=0;
  float st_vtxx=0, st_vtxy=0, st_vtxz=0, st_vtxt=0;
  float st_pxnu=0, st_pynu=0, st_pznu=0, st_enu=0;

  st->Branch("iev",  &st_iev,  "iev/I" );
  st->Branch("vtxx", &st_vtxx, "vtxx/F");
  st->Branch("vtxy", &st_vtxy, "vtxy/F");
  st->Branch("vtxz", &st_vtxz, "vtxz/F");
  st->Branch("vtxt", &st_vtxt, "vtxt/F");
  st->Branch("px",   &st_pxnu, "px/F"  );
  st->Branch("py",   &st_pynu, "py/F"  );
  st->Branch("pz",   &st_pznu, "pz/F"  );
  st->Branch("E",    &st_enu,  "E/F"   );

  //-- loop over TTree NtpMC records & fill in the summary tree
  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("Main", pINFO) << rec_header;
    LOG("Main", pINFO) << event;

    TLorentzVector * vtx = event.GetVertex();

    st_iev  = i;
    st_vtxx = vtx->X();
    st_vtxy = vtx->Y();
    st_vtxz = vtx->Z();
    st_vtxt = vtx->T();
    st_pxnu = event.GetParticle(0)->P4()->Px();
    st_pynu = event.GetParticle(0)->P4()->Py();
    st_pznu = event.GetParticle(0)->P4()->Pz();
    st_enu  = event.GetParticle(0)->P4()->E();

    st->Fill();

    //delete vtx;
  }

  TCanvas * c = 0;

  TPostScript * ps = new TPostScript("gmctest.ps", 111);

  ps->NewPage();
  c = new TCanvas("c","",20,20,500,500);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  c->Draw();
  c->Divide(2,2);
  c->cd(1);
  st->Draw("vtxx");
  c->cd(2);
  st->Draw("vtxy");
  c->cd(3);
  st->Draw("vtxz");
  c->cd(4);
  st->Draw("vtxx:vtxy:vtxz");
  c->Update();
  delete c;

  ps->NewPage();
  c = new TCanvas("c","",20,20,500,500);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  c->Draw();
  c->Divide(2,2);
  c->cd(1);
  st->Draw("px");
  c->cd(2);
  st->Draw("py");
  c->cd(3);
  st->Draw("pz");
  c->cd(4);
  st->Draw("E");
  c->Update();
  delete c;

  ps->Close();
  delete ps;

  TFile outfile("gmctest.root","recreate");
  st->Write();
  outfile.Close();

  delete st;

  inpfile.Close();
  tree = 0;
  thdr = 0;

  LOG("Main", pINFO)  << "Done!";
  return 0;
}
//___________________________________________________________________
string getRootFilename(int argc, char ** argv)
{
  //-- Scan for filename from the command line argument (following -f)
  string filename="";
  for(int iarg = 0; iarg < argc-1; iarg++) {
   string argument(argv[iarg]);
   if( argument.compare("-f") == 0 ) filename = string(argv[++iarg]);
  }
  return filename;
}
//___________________________________________________________________
bool checkRootFilename(string filename)
{
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("Main", pERROR)
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//___________________________________________________________________
