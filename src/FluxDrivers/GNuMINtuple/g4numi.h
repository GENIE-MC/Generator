//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec  4 16:19:52 2008 by ROOT version 5.21/07
// from TTree g4numi/g4numi Neutrino ntuple
// found on file: generic_g4numi.root
//////////////////////////////////////////////////////////

#ifndef g4numi_h
#define g4numi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>

#include <iostream>

class g4numi {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //data_t          *data;
   Int_t           run;
   Int_t           evtno;
   Double_t        beamHWidth;
   Double_t        beamVWidth;
   Double_t        beamX;
   Double_t        beamY;
   Double_t        protonX;
   Double_t        protonY;
   Double_t        protonZ;
   Double_t        protonPx;
   Double_t        protonPy;
   Double_t        protonPz;
   Double_t        nuTarZ;
   Double_t        hornCurrent;
   Double_t        Ndxdz;
   Double_t        Ndydz;
   Double_t        Npz;
   Double_t        Nenergy;
   Double_t        NdxdzNear[11];
   Double_t        NdydzNear[11];
   Double_t        NenergyN[11];
   Double_t        NWtNear[11];
   Double_t        NdxdzFar[2];
   Double_t        NdydzFar[2];
   Double_t        NenergyF[2];
   Double_t        NWtFar[2];
   Int_t           Norig;
   Int_t           Ndecay;
   Int_t           Ntype;
   Double_t        Vx;
   Double_t        Vy;
   Double_t        Vz;
   Double_t        pdPx;
   Double_t        pdPy;
   Double_t        pdPz;
   Double_t        ppdxdz;
   Double_t        ppdydz;
   Double_t        pppz;
   Double_t        ppenergy;
   Double_t        ppmedium;
   Int_t           ptype;
   Double_t        ppvx;
   Double_t        ppvy;
   Double_t        ppvz;
   Double_t        muparpx;
   Double_t        muparpy;
   Double_t        muparpz;
   Double_t        mupare;
   Double_t        Necm;
   Double_t        Nimpwt;
   Double_t        xpoint;
   Double_t        ypoint;
   Double_t        zpoint;
   Double_t        tvx;
   Double_t        tvy;
   Double_t        tvz;
   Double_t        tpx;
   Double_t        tpy;
   Double_t        tpz;
   Int_t           tptype;
   Int_t           tgen;

#ifndef SKIP_MINERVA_MODS
   //=========================================
   // The following was inserted by MINERvA
   //=========================================
   Int_t           ntrajectory;
   Bool_t          overflow;
   Int_t           pdg[10];
   Int_t           trackId[10];
   Int_t           parentId[10];
   Double_t        startx[10];
   Double_t        starty[10];
   Double_t        startz[10];
   Double_t        stopx[10];
   Double_t        stopy[10];
   Double_t        stopz[10];
   Double_t        startpx[10];
   Double_t        startpy[10];
   Double_t        startpz[10];
   Double_t        stoppx[10];
   Double_t        stoppy[10];
   Double_t        stoppz[10];
   Double_t        pprodpx[10];
   Double_t        pprodpy[10];
   Double_t        pprodpz[10];
   TString         proc[10];
   TString         ivol[10];
   TString         fvol[10];   
   //END of minerva additions
#endif

   Double_t        trkx[10];
   Double_t        trky[10];
   Double_t        trkz[10];
   Double_t        trkpx[10];
   Double_t        trkpy[10];
   Double_t        trkpz[10];

   // List of branches
   TBranch        *b_data_run;   //!
   TBranch        *b_data_evtno;   //!
   TBranch        *b_data_beamHWidth;   //!
   TBranch        *b_data_beamVWidth;   //!
   TBranch        *b_data_beamX;   //!
   TBranch        *b_data_beamY;   //!
   TBranch        *b_data_protonX;   //!
   TBranch        *b_data_protonY;   //!
   TBranch        *b_data_protonZ;   //!
   TBranch        *b_data_protonPx;   //!
   TBranch        *b_data_protonPy;   //!
   TBranch        *b_data_protonPz;   //!
   TBranch        *b_data_nuTarZ;   //!
   TBranch        *b_data_hornCurrent;   //!
   TBranch        *b_data_Ndxdz;   //!
   TBranch        *b_data_Ndydz;   //!
   TBranch        *b_data_Npz;   //!
   TBranch        *b_data_Nenergy;   //!
   TBranch        *b_data_NdxdzNear;   //!
   TBranch        *b_data_NdydzNear;   //!
   TBranch        *b_data_NenergyN;   //!
   TBranch        *b_data_NWtNear;   //!
   TBranch        *b_data_NdxdzFar;   //!
   TBranch        *b_data_NdydzFar;   //!
   TBranch        *b_data_NenergyF;   //!
   TBranch        *b_data_NWtFar;   //!
   TBranch        *b_data_Norig;   //!
   TBranch        *b_data_Ndecay;   //!
   TBranch        *b_data_Ntype;   //!
   TBranch        *b_data_Vx;   //!
   TBranch        *b_data_Vy;   //!
   TBranch        *b_data_Vz;   //!
   TBranch        *b_data_pdPx;   //!
   TBranch        *b_data_pdPy;   //!
   TBranch        *b_data_pdPz;   //!
   TBranch        *b_data_ppdxdz;   //!
   TBranch        *b_data_ppdydz;   //!
   TBranch        *b_data_pppz;   //!
   TBranch        *b_data_ppenergy;   //!
   TBranch        *b_data_ppmedium;   //!
   TBranch        *b_data_ptype;   //!
   TBranch        *b_data_ppvx;   //!
   TBranch        *b_data_ppvy;   //!
   TBranch        *b_data_ppvz;   //!
   TBranch        *b_data_muparpx;   //!
   TBranch        *b_data_muparpy;   //!
   TBranch        *b_data_muparpz;   //!
   TBranch        *b_data_mupare;   //!
   TBranch        *b_data_Necm;   //!
   TBranch        *b_data_Nimpwt;   //!
   TBranch        *b_data_xpoint;   //!
   TBranch        *b_data_ypoint;   //!
   TBranch        *b_data_zpoint;   //!
   TBranch        *b_data_tvx;   //!
   TBranch        *b_data_tvy;   //!
   TBranch        *b_data_tvz;   //!
   TBranch        *b_data_tpx;   //!
   TBranch        *b_data_tpy;   //!
   TBranch        *b_data_tpz;   //!
   TBranch        *b_data_tptype;   //!
   TBranch        *b_data_tgen;   //!

#ifndef SKIP_MINERVA_MODS
   //=========================================
   // The following was inserted by MINERvA
   //=========================================
   TBranch        *b_data_ntrajectory;   //!
   TBranch        *b_data_overflow;   //!
   TBranch        *b_data_pdg;   //!
   TBranch        *b_data_trackId;   //!
   TBranch        *b_data_parentId;   //!
   TBranch        *b_data_startx;   //!
   TBranch        *b_data_starty;   //!
   TBranch        *b_data_startz;   //!
   TBranch        *b_data_stopx;   //!
   TBranch        *b_data_stopy;   //!
   TBranch        *b_data_stopz;   //!
   TBranch        *b_data_startpx;   //!
   TBranch        *b_data_startpy;   //!
   TBranch        *b_data_startpz;   //!
   TBranch        *b_data_stoppx;   //!
   TBranch        *b_data_stoppy;   //!
   TBranch        *b_data_stoppz;   //!
   TBranch        *b_data_pprodpx;   //!
   TBranch        *b_data_pprodpy;   //!
   TBranch        *b_data_pprodpz;   //!
   TBranch        *b_data_proc;   //!
   TBranch        *b_data_ivol;   //!
   TBranch        *b_data_fvol;   //!
   //END of minerva additions
#endif

   TBranch        *b_data_trkx;   //!
   TBranch        *b_data_trky;   //!
   TBranch        *b_data_trkz;   //!
   TBranch        *b_data_trkpx;   //!
   TBranch        *b_data_trkpy;   //!
   TBranch        *b_data_trkpz;   //!

   g4numi(TTree *tree=0);
   virtual ~g4numi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef g4numi_cxx
g4numi::g4numi(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("generic_g4numi.root");
      if (!f) {
         f = new TFile("generic_g4numi.root");
      }
      tree = (TTree*)gDirectory->Get("g4numi");

   }
   Init(tree);
}

g4numi::~g4numi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t g4numi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t g4numi::LoadTree(Long64_t entry)
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

void g4numi::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_data_run);
   fChain->SetBranchAddress("evtno", &evtno, &b_data_evtno);
   fChain->SetBranchAddress("beamHWidth", &beamHWidth, &b_data_beamHWidth);
   fChain->SetBranchAddress("beamVWidth", &beamVWidth, &b_data_beamVWidth);
   fChain->SetBranchAddress("beamX", &beamX, &b_data_beamX);
   fChain->SetBranchAddress("beamY", &beamY, &b_data_beamY);
   fChain->SetBranchAddress("protonX", &protonX, &b_data_protonX);
   fChain->SetBranchAddress("protonY", &protonY, &b_data_protonY);
   fChain->SetBranchAddress("protonZ", &protonZ, &b_data_protonZ);
   fChain->SetBranchAddress("protonPx", &protonPx, &b_data_protonPx);
   fChain->SetBranchAddress("protonPy", &protonPy, &b_data_protonPy);
   fChain->SetBranchAddress("protonPz", &protonPz, &b_data_protonPz);
   fChain->SetBranchAddress("nuTarZ", &nuTarZ, &b_data_nuTarZ);
   fChain->SetBranchAddress("hornCurrent", &hornCurrent, &b_data_hornCurrent);
   fChain->SetBranchAddress("Ndxdz", &Ndxdz, &b_data_Ndxdz);
   fChain->SetBranchAddress("Ndydz", &Ndydz, &b_data_Ndydz);
   fChain->SetBranchAddress("Npz", &Npz, &b_data_Npz);
   fChain->SetBranchAddress("Nenergy", &Nenergy, &b_data_Nenergy);
   fChain->SetBranchAddress("NdxdzNear[11]", NdxdzNear, &b_data_NdxdzNear);
   fChain->SetBranchAddress("NdydzNear[11]", NdydzNear, &b_data_NdydzNear);
   fChain->SetBranchAddress("NenergyN[11]", NenergyN, &b_data_NenergyN);
   fChain->SetBranchAddress("NWtNear[11]", NWtNear, &b_data_NWtNear);
   fChain->SetBranchAddress("NdxdzFar[2]", NdxdzFar, &b_data_NdxdzFar);
   fChain->SetBranchAddress("NdydzFar[2]", NdydzFar, &b_data_NdydzFar);
   fChain->SetBranchAddress("NenergyF[2]", NenergyF, &b_data_NenergyF);
   fChain->SetBranchAddress("NWtFar[2]", NWtFar, &b_data_NWtFar);
   fChain->SetBranchAddress("Norig", &Norig, &b_data_Norig);
   fChain->SetBranchAddress("Ndecay", &Ndecay, &b_data_Ndecay);
   fChain->SetBranchAddress("Ntype", &Ntype, &b_data_Ntype);
   fChain->SetBranchAddress("Vx", &Vx, &b_data_Vx);
   fChain->SetBranchAddress("Vy", &Vy, &b_data_Vy);
   fChain->SetBranchAddress("Vz", &Vz, &b_data_Vz);
   fChain->SetBranchAddress("pdPx", &pdPx, &b_data_pdPx);
   fChain->SetBranchAddress("pdPy", &pdPy, &b_data_pdPy);
   fChain->SetBranchAddress("pdPz", &pdPz, &b_data_pdPz);
   fChain->SetBranchAddress("ppdxdz", &ppdxdz, &b_data_ppdxdz);
   fChain->SetBranchAddress("ppdydz", &ppdydz, &b_data_ppdydz);
   fChain->SetBranchAddress("pppz", &pppz, &b_data_pppz);
   fChain->SetBranchAddress("ppenergy", &ppenergy, &b_data_ppenergy);
   fChain->SetBranchAddress("ppmedium", &ppmedium, &b_data_ppmedium);
   fChain->SetBranchAddress("ptype", &ptype, &b_data_ptype);
   fChain->SetBranchAddress("ppvx", &ppvx, &b_data_ppvx);
   fChain->SetBranchAddress("ppvy", &ppvy, &b_data_ppvy);
   fChain->SetBranchAddress("ppvz", &ppvz, &b_data_ppvz);
   fChain->SetBranchAddress("muparpx", &muparpx, &b_data_muparpx);
   fChain->SetBranchAddress("muparpy", &muparpy, &b_data_muparpy);
   fChain->SetBranchAddress("muparpz", &muparpz, &b_data_muparpz);
   fChain->SetBranchAddress("mupare", &mupare, &b_data_mupare);
   fChain->SetBranchAddress("Necm", &Necm, &b_data_Necm);
   fChain->SetBranchAddress("Nimpwt", &Nimpwt, &b_data_Nimpwt);
   fChain->SetBranchAddress("xpoint", &xpoint, &b_data_xpoint);
   fChain->SetBranchAddress("ypoint", &ypoint, &b_data_ypoint);
   fChain->SetBranchAddress("zpoint", &zpoint, &b_data_zpoint);
   fChain->SetBranchAddress("tvx", &tvx, &b_data_tvx);
   fChain->SetBranchAddress("tvy", &tvy, &b_data_tvy);
   fChain->SetBranchAddress("tvz", &tvz, &b_data_tvz);
   fChain->SetBranchAddress("tpx", &tpx, &b_data_tpx);
   fChain->SetBranchAddress("tpy", &tpy, &b_data_tpy);
   fChain->SetBranchAddress("tpz", &tpz, &b_data_tpz);
   fChain->SetBranchAddress("tptype", &tptype, &b_data_tptype);
   fChain->SetBranchAddress("tgen", &tgen, &b_data_tgen);

#ifndef SKIP_MINERVA_MODS
   // RWH 2014-09-17 don't attach branches to addresses if they don't exist
   bool has_ancestry = ( fChain->GetBranch("ntrajectory") != 0 );
   std::cout << "g4numi ntuple " 
             << ((has_ancestry)?"has":"lacks") << " ancestry branches" 
             << std::endl;
   if ( has_ancestry ) {
   //=========================================
   // The following was inserted by MINERvA
   //=========================================
   fChain->SetBranchAddress("ntrajectory", &ntrajectory, &b_data_ntrajectory);
   fChain->SetBranchAddress("overflow", &overflow, &b_data_overflow);
   fChain->SetBranchAddress("pdg[10]", pdg, &b_data_pdg);
   fChain->SetBranchAddress("trackId[10]", trackId, &b_data_trackId);
   fChain->SetBranchAddress("parentId[10]", parentId, &b_data_parentId);
   fChain->SetBranchAddress("startx[10]", startx, &b_data_startx);
   fChain->SetBranchAddress("starty[10]", starty, &b_data_starty);
   fChain->SetBranchAddress("startz[10]", startz, &b_data_startz);
   fChain->SetBranchAddress("stopx[10]", stopx, &b_data_stopx);
   fChain->SetBranchAddress("stopy[10]", stopy, &b_data_stopy);
   fChain->SetBranchAddress("stopz[10]", stopz, &b_data_stopz);
   fChain->SetBranchAddress("startpx[10]", startpx, &b_data_startpx);
   fChain->SetBranchAddress("startpy[10]", startpy, &b_data_startpy);
   fChain->SetBranchAddress("startpz[10]", startpz, &b_data_startpz);
   fChain->SetBranchAddress("stoppx[10]", stoppx, &b_data_stoppx);
   fChain->SetBranchAddress("stoppy[10]", stoppy, &b_data_stoppy);
   fChain->SetBranchAddress("stoppz[10]", stoppz, &b_data_stoppz);
   fChain->SetBranchAddress("pprodpx[10]", pprodpx, &b_data_pprodpx);
   fChain->SetBranchAddress("pprodpy[10]", pprodpy, &b_data_pprodpy);
   fChain->SetBranchAddress("pprodpz[10]", pprodpz, &b_data_pprodpz);
   fChain->SetBranchAddress("proc[10]", proc, &b_data_proc);
   fChain->SetBranchAddress("ivol[10]", ivol, &b_data_ivol);
   fChain->SetBranchAddress("fvol[10]", fvol, &b_data_fvol);
   //END of minerva additions
   }
#endif

   fChain->SetBranchAddress("trkx[10]", trkx, &b_data_trkx);
   fChain->SetBranchAddress("trky[10]", trky, &b_data_trky);
   fChain->SetBranchAddress("trkz[10]", trkz, &b_data_trkz);
   fChain->SetBranchAddress("trkpx[10]", trkpx, &b_data_trkpx);
   fChain->SetBranchAddress("trkpy[10]", trkpy, &b_data_trkpy);
   fChain->SetBranchAddress("trkpz[10]", trkpz, &b_data_trkpz);
   Notify();
}

Bool_t g4numi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void g4numi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t g4numi::Cut(Long64_t /* entry */)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef g4numi_cxx
