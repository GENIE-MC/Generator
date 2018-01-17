//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec  4 16:19:50 2008 by ROOT version 5.21/07
// from TTree g3numi/neutrino
// found on file: generic_g3numi.root
//////////////////////////////////////////////////////////

#ifndef g3numi_h
#define g3numi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class g3numi {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           evtno;
   Float_t         Ndxdz;
   Float_t         Ndydz;
   Float_t         Npz;
   Float_t         Nenergy;
   Float_t         Ndxdznea;
   Float_t         Ndydznea;
   Float_t         Nenergyn;
   Float_t         Nwtnear;
   Float_t         Ndxdzfar;
   Float_t         Ndydzfar;
   Float_t         Nenergyf;
   Float_t         Nwtfar;
   Int_t           Norig;
   Int_t           Ndecay;
   Int_t           Ntype;
   Float_t         Vx;
   Float_t         Vy;
   Float_t         Vz;
   Float_t         pdpx;
   Float_t         pdpy;
   Float_t         pdpz;
   Float_t         ppdxdz;
   Float_t         ppdydz;
   Float_t         pppz;
   Float_t         ppenergy;
   Int_t           ppmedium;
   Int_t           ptype;
   Float_t         ppvx;
   Float_t         ppvy;
   Float_t         ppvz;
   Float_t         muparpx;
   Float_t         muparpy;
   Float_t         muparpz;
   Float_t         mupare;
   Float_t         Necm;
   Float_t         Nimpwt;
   Float_t         xpoint;
   Float_t         ypoint;
   Float_t         zpoint;
   Float_t         tvx;
   Float_t         tvy;
   Float_t         tvz;
   Float_t         tpx;
   Float_t         tpy;
   Float_t         tpz;
   Int_t           tptype;
   Int_t           tgen;
   Int_t           tgptype;
   Float_t         tgppx;
   Float_t         tgppy;
   Float_t         tgppz;
   Float_t         tprivx;
   Float_t         tprivy;
   Float_t         tprivz;
   Float_t         beamx;
   Float_t         beamy;
   Float_t         beamz;
   Float_t         beampx;
   Float_t         beampy;
   Float_t         beampz;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evtno;   //!
   TBranch        *b_Ndxdz;   //!
   TBranch        *b_Ndydz;   //!
   TBranch        *b_Npz;   //!
   TBranch        *b_Nenergy;   //!
   TBranch        *b_Ndxdznea;   //!
   TBranch        *b_Ndydznea;   //!
   TBranch        *b_Nenergyn;   //!
   TBranch        *b_Nwtnear;   //!
   TBranch        *b_Ndxdzfar;   //!
   TBranch        *b_Ndydzfar;   //!
   TBranch        *b_Nenergyf;   //!
   TBranch        *b_Nwtfar;   //!
   TBranch        *b_Norig;   //!
   TBranch        *b_Ndecay;   //!
   TBranch        *b_Ntype;   //!
   TBranch        *b_Vx;   //!
   TBranch        *b_Vy;   //!
   TBranch        *b_Vz;   //!
   TBranch        *b_pdpx;   //!
   TBranch        *b_pdpy;   //!
   TBranch        *b_pdpz;   //!
   TBranch        *b_ppdxdz;   //!
   TBranch        *b_ppdydz;   //!
   TBranch        *b_pppz;   //!
   TBranch        *b_ppenergy;   //!
   TBranch        *b_ppmedium;   //!
   TBranch        *b_ptype;   //!
   TBranch        *b_ppvx;   //!
   TBranch        *b_ppvy;   //!
   TBranch        *b_ppvz;   //!
   TBranch        *b_muparpx;   //!
   TBranch        *b_muparpy;   //!
   TBranch        *b_muparpz;   //!
   TBranch        *b_mupare;   //!
   TBranch        *b_Necm;   //!
   TBranch        *b_Nimpwt;   //!
   TBranch        *b_xpoint;   //!
   TBranch        *b_ypoint;   //!
   TBranch        *b_zpoint;   //!
   TBranch        *b_tvx;   //!
   TBranch        *b_tvy;   //!
   TBranch        *b_tvz;   //!
   TBranch        *b_tpx;   //!
   TBranch        *b_tpy;   //!
   TBranch        *b_tpz;   //!
   TBranch        *b_tptype;   //!
   TBranch        *b_tgen;   //!
   TBranch        *b_tgptype;   //!
   TBranch        *b_tgppx;   //!
   TBranch        *b_tgppy;   //!
   TBranch        *b_tgppz;   //!
   TBranch        *b_tprivx;   //!
   TBranch        *b_tprivy;   //!
   TBranch        *b_tprivz;   //!
   TBranch        *b_beamx;   //!
   TBranch        *b_beamy;   //!
   TBranch        *b_beamz;   //!
   TBranch        *b_beampx;   //!
   TBranch        *b_beampy;   //!
   TBranch        *b_beampz;   //!

   g3numi(TTree *tree=0);
   virtual ~g3numi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef g3numi_cxx
g3numi::g3numi(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("generic_g3numi.root");
      if (!f) {
         f = new TFile("generic_g3numi.root");
      }
      tree = (TTree*)gDirectory->Get("g3numi");

   }
   Init(tree);
}

g3numi::~g3numi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t g3numi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t g3numi::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void g3numi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evtno", &evtno, &b_evtno);
   fChain->SetBranchAddress("Ndxdz", &Ndxdz, &b_Ndxdz);
   fChain->SetBranchAddress("Ndydz", &Ndydz, &b_Ndydz);
   fChain->SetBranchAddress("Npz", &Npz, &b_Npz);
   fChain->SetBranchAddress("Nenergy", &Nenergy, &b_Nenergy);
   fChain->SetBranchAddress("Ndxdznea", &Ndxdznea, &b_Ndxdznea);
   fChain->SetBranchAddress("Ndydznea", &Ndydznea, &b_Ndydznea);
   fChain->SetBranchAddress("Nenergyn", &Nenergyn, &b_Nenergyn);
   fChain->SetBranchAddress("Nwtnear", &Nwtnear, &b_Nwtnear);
   fChain->SetBranchAddress("Ndxdzfar", &Ndxdzfar, &b_Ndxdzfar);
   fChain->SetBranchAddress("Ndydzfar", &Ndydzfar, &b_Ndydzfar);
   fChain->SetBranchAddress("Nenergyf", &Nenergyf, &b_Nenergyf);
   fChain->SetBranchAddress("Nwtfar", &Nwtfar, &b_Nwtfar);
   fChain->SetBranchAddress("Norig", &Norig, &b_Norig);
   fChain->SetBranchAddress("Ndecay", &Ndecay, &b_Ndecay);
   fChain->SetBranchAddress("Ntype", &Ntype, &b_Ntype);
   fChain->SetBranchAddress("Vx", &Vx, &b_Vx);
   fChain->SetBranchAddress("Vy", &Vy, &b_Vy);
   fChain->SetBranchAddress("Vz", &Vz, &b_Vz);
   fChain->SetBranchAddress("pdpx", &pdpx, &b_pdpx);
   fChain->SetBranchAddress("pdpy", &pdpy, &b_pdpy);
   fChain->SetBranchAddress("pdpz", &pdpz, &b_pdpz);
   fChain->SetBranchAddress("ppdxdz", &ppdxdz, &b_ppdxdz);
   fChain->SetBranchAddress("ppdydz", &ppdydz, &b_ppdydz);
   fChain->SetBranchAddress("pppz", &pppz, &b_pppz);
   fChain->SetBranchAddress("ppenergy", &ppenergy, &b_ppenergy);
   fChain->SetBranchAddress("ppmedium", &ppmedium, &b_ppmedium);
   fChain->SetBranchAddress("ptype", &ptype, &b_ptype);
   fChain->SetBranchAddress("ppvx", &ppvx, &b_ppvx);
   fChain->SetBranchAddress("ppvy", &ppvy, &b_ppvy);
   fChain->SetBranchAddress("ppvz", &ppvz, &b_ppvz);
   fChain->SetBranchAddress("muparpx", &muparpx, &b_muparpx);
   fChain->SetBranchAddress("muparpy", &muparpy, &b_muparpy);
   fChain->SetBranchAddress("muparpz", &muparpz, &b_muparpz);
   fChain->SetBranchAddress("mupare", &mupare, &b_mupare);
   fChain->SetBranchAddress("Necm", &Necm, &b_Necm);
   fChain->SetBranchAddress("Nimpwt", &Nimpwt, &b_Nimpwt);
   fChain->SetBranchAddress("xpoint", &xpoint, &b_xpoint);
   fChain->SetBranchAddress("ypoint", &ypoint, &b_ypoint);
   fChain->SetBranchAddress("zpoint", &zpoint, &b_zpoint);
   fChain->SetBranchAddress("tvx", &tvx, &b_tvx);
   fChain->SetBranchAddress("tvy", &tvy, &b_tvy);
   fChain->SetBranchAddress("tvz", &tvz, &b_tvz);
   fChain->SetBranchAddress("tpx", &tpx, &b_tpx);
   fChain->SetBranchAddress("tpy", &tpy, &b_tpy);
   fChain->SetBranchAddress("tpz", &tpz, &b_tpz);
   fChain->SetBranchAddress("tptype", &tptype, &b_tptype);
   fChain->SetBranchAddress("tgen", &tgen, &b_tgen);
   fChain->SetBranchAddress("tgptype", &tgptype, &b_tgptype);
   fChain->SetBranchAddress("tgppx", &tgppx, &b_tgppx);
   fChain->SetBranchAddress("tgppy", &tgppy, &b_tgppy);
   fChain->SetBranchAddress("tgppz", &tgppz, &b_tgppz);
   fChain->SetBranchAddress("tprivx", &tprivx, &b_tprivx);
   fChain->SetBranchAddress("tprivy", &tprivy, &b_tprivy);
   fChain->SetBranchAddress("tprivz", &tprivz, &b_tprivz);
   fChain->SetBranchAddress("beamx", &beamx, &b_beamx);
   fChain->SetBranchAddress("beamy", &beamy, &b_beamy);
   fChain->SetBranchAddress("beamz", &beamz, &b_beamz);
   fChain->SetBranchAddress("beampx", &beampx, &b_beampx);
   fChain->SetBranchAddress("beampy", &beampy, &b_beampy);
   fChain->SetBranchAddress("beampz", &beampz, &b_beampz);
   Notify();
}

Bool_t g3numi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void g3numi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t g3numi::Cut(Long64_t /* entry */)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef g3numi_cxx
