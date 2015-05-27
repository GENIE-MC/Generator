#ifndef GST_H
#define GST_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class gst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           iev;
   Int_t           neu;
   Int_t           tgt;
   Int_t           Z;
   Int_t           A;
   Int_t           hitnuc;
   Int_t           hitqrk;
   Int_t           resid;
   Bool_t          sea;
   Bool_t          qel;
   Bool_t          res;
   Bool_t          dis;
   Bool_t          coh;
   Bool_t          dfr;
   Bool_t          imd;
   Bool_t          nuel;
   Bool_t          em;
   Bool_t          cc;
   Bool_t          nc;
   Bool_t          charm;
   Int_t           neut_code;
   Int_t           nuance_code;
   Double_t        wght;
   Double_t        xs;
   Double_t        ys;
   Double_t        ts;
   Double_t        Q2s;
   Double_t        Ws;
   Double_t        x;
   Double_t        y;
   Double_t        t;
   Double_t        Q2;
   Double_t        W;
   Double_t        Ev;
   Double_t        pxv;
   Double_t        pyv;
   Double_t        pzv;
   Double_t        En;
   Double_t        pxn;
   Double_t        pyn;
   Double_t        pzn;
   Double_t        El;
   Double_t        pxl;
   Double_t        pyl;
   Double_t        pzl;
   Int_t           nfp;
   Int_t           nfn;
   Int_t           nfpip;
   Int_t           nfpim;
   Int_t           nfpi0;
   Int_t           nfkp;
   Int_t           nfkm;
   Int_t           nfk0;
   Int_t           nfem;
   Int_t           nfother;
   Int_t           nip;
   Int_t           nin;
   Int_t           nipip;
   Int_t           nipim;
   Int_t           nipi0;
   Int_t           nikp;
   Int_t           nikm;
   Int_t           nik0;
   Int_t           niem;
   Int_t           niother;
   Int_t           ni;
   Int_t           pdgi[60];   //[ni]
   Int_t           resc[60];   //[ni]
   Double_t        Ei[60];   //[ni]
   Double_t        pxi[60];   //[ni]
   Double_t        pyi[60];   //[ni]
   Double_t        pzi[60];   //[ni]
   Int_t           nf;
   Int_t           pdgf[60];   //[nf]
   Double_t        Ef[60];   //[nf]
   Double_t        pxf[60];   //[nf]
   Double_t        pyf[60];   //[nf]
   Double_t        pzf[60];   //[nf]
   Double_t        vtxx;
   Double_t        vtxy;
   Double_t        vtxz;
   Double_t        vtxt;
   Double_t        calresp0;

   // List of branches
   TBranch        *b_iev;   //!
   TBranch        *b_neu;   //!
   TBranch        *b_tgt;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_A;   //!
   TBranch        *b_hitnuc;   //!
   TBranch        *b_hitqrk;   //!
   TBranch        *b_resid;   //!
   TBranch        *b_sea;   //!
   TBranch        *b_qel;   //!
   TBranch        *b_res;   //!
   TBranch        *b_dis;   //!
   TBranch        *b_coh;   //!
   TBranch        *b_dfr;   //!
   TBranch        *b_imd;   //!
   TBranch        *b_nuel;   //!
   TBranch        *b_em;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_nc;   //!
   TBranch        *b_charm;   //!
   TBranch        *b_neut_code;   //!
   TBranch        *b_nuance_code;   //!
   TBranch        *b_wght;   //!
   TBranch        *b_xs;   //!
   TBranch        *b_ys;   //!
   TBranch        *b_ts;   //!
   TBranch        *b_Q2s;   //!
   TBranch        *b_Ws;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_t;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_Ev;   //!
   TBranch        *b_pxv;   //!
   TBranch        *b_pyv;   //!
   TBranch        *b_pzv;   //!
   TBranch        *b_En;   //!
   TBranch        *b_pxn;   //!
   TBranch        *b_pyn;   //!
   TBranch        *b_pzn;   //!
   TBranch        *b_El;   //!
   TBranch        *b_pxl;   //!
   TBranch        *b_pyl;   //!
   TBranch        *b_pzl;   //!
   TBranch        *b_nfp;   //!
   TBranch        *b_nfn;   //!
   TBranch        *b_nfpip;   //!
   TBranch        *b_nfpim;   //!
   TBranch        *b_nfpi0;   //!
   TBranch        *b_nfkp;   //!
   TBranch        *b_nfkm;   //!
   TBranch        *b_nfk0;   //!
   TBranch        *b_nfem;   //!
   TBranch        *b_nfother;   //!
   TBranch        *b_np;   //!
   TBranch        *b_nn;   //!
   TBranch        *b_npip;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_nkp;   //!
   TBranch        *b_nkm;   //!
   TBranch        *b_nk0;   //!
   TBranch        *b_niem;   //!
   TBranch        *b_niother;   //!
   TBranch        *b_ni;   //!
   TBranch        *b_pdgi;   //!
   TBranch        *b_resc;   //!
   TBranch        *b_Ei;   //!
   TBranch        *b_pxi;   //!
   TBranch        *b_pyi;   //!
   TBranch        *b_pzi;   //!
   TBranch        *b_nf;   //!
   TBranch        *b_pdgf;   //!
   TBranch        *b_Ef;   //!
   TBranch        *b_pxf;   //!
   TBranch        *b_pyf;   //!
   TBranch        *b_pzf;   //!
   TBranch        *b_vtxx;   //!
   TBranch        *b_vtxy;   //!
   TBranch        *b_vtxz;   //!
   TBranch        *b_vtxt;   //!
   TBranch        *b_calresp0;   //!

   gst(char* fname);
   virtual ~gst();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
};

gst::gst(char* fname)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
   if (!f || !f->IsOpen()) {
      f = new TFile(fname);
   }
   f->GetObject("gst",fChain);
   Init(fChain);
}
gst::~gst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
void gst::Init(TTree *tree)
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

   fChain->SetBranchAddress("iev", &iev, &b_iev);
   fChain->SetBranchAddress("neu", &neu, &b_neu);
   fChain->SetBranchAddress("tgt", &tgt, &b_tgt);
   fChain->SetBranchAddress("Z", &Z, &b_Z);
   fChain->SetBranchAddress("A", &A, &b_A);
   fChain->SetBranchAddress("hitnuc", &hitnuc, &b_hitnuc);
   fChain->SetBranchAddress("hitqrk", &hitqrk, &b_hitqrk);
   fChain->SetBranchAddress("resid", &resid, &b_resid);
   fChain->SetBranchAddress("sea", &sea, &b_sea);
   fChain->SetBranchAddress("qel", &qel, &b_qel);
   fChain->SetBranchAddress("res", &res, &b_res);
   fChain->SetBranchAddress("dis", &dis, &b_dis);
   fChain->SetBranchAddress("coh", &coh, &b_coh);
   fChain->SetBranchAddress("dfr", &dfr, &b_dfr);
   fChain->SetBranchAddress("imd", &imd, &b_imd);
   fChain->SetBranchAddress("nuel", &nuel, &b_nuel);
   fChain->SetBranchAddress("em", &em, &b_em);
   fChain->SetBranchAddress("cc", &cc, &b_cc);
   fChain->SetBranchAddress("nc", &nc, &b_nc);
   fChain->SetBranchAddress("charm", &charm, &b_charm);
   fChain->SetBranchAddress("neut_code", &neut_code, &b_neut_code);
   fChain->SetBranchAddress("nuance_code", &nuance_code, &b_nuance_code);
   fChain->SetBranchAddress("wght", &wght, &b_wght);
   fChain->SetBranchAddress("xs", &xs, &b_xs);
   fChain->SetBranchAddress("ys", &ys, &b_ys);
   fChain->SetBranchAddress("ts", &ts, &b_ts);
   fChain->SetBranchAddress("Q2s", &Q2s, &b_Q2s);
   fChain->SetBranchAddress("Ws", &Ws, &b_Ws);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("W", &W, &b_W);
   fChain->SetBranchAddress("Ev", &Ev, &b_Ev);
   fChain->SetBranchAddress("pxv", &pxv, &b_pxv);
   fChain->SetBranchAddress("pyv", &pyv, &b_pyv);
   fChain->SetBranchAddress("pzv", &pzv, &b_pzv);
   fChain->SetBranchAddress("En", &En, &b_En);
   fChain->SetBranchAddress("pxn", &pxn, &b_pxn);
   fChain->SetBranchAddress("pyn", &pyn, &b_pyn);
   fChain->SetBranchAddress("pzn", &pzn, &b_pzn);
   fChain->SetBranchAddress("El", &El, &b_El);
   fChain->SetBranchAddress("pxl", &pxl, &b_pxl);
   fChain->SetBranchAddress("pyl", &pyl, &b_pyl);
   fChain->SetBranchAddress("pzl", &pzl, &b_pzl);
   fChain->SetBranchAddress("nfp", &nfp, &b_nfp);
   fChain->SetBranchAddress("nfn", &nfn, &b_nfn);
   fChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
   fChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
   fChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
   fChain->SetBranchAddress("nfkp", &nfkp, &b_nfkp);
   fChain->SetBranchAddress("nfkm", &nfkm, &b_nfkm);
   fChain->SetBranchAddress("nfk0", &nfk0, &b_nfk0);
   fChain->SetBranchAddress("nfem", &nfem, &b_nfem);
   fChain->SetBranchAddress("nfother", &nfother, &b_nfother);
   fChain->SetBranchAddress("nip", &nip, &b_np);
   fChain->SetBranchAddress("nin", &nin, &b_nn);
   fChain->SetBranchAddress("nipip", &nipip, &b_npip);
   fChain->SetBranchAddress("nipim", &nipim, &b_npim);
   fChain->SetBranchAddress("nipi0", &nipi0, &b_npi0);
   fChain->SetBranchAddress("nikp", &nikp, &b_nkp);
   fChain->SetBranchAddress("nikm", &nikm, &b_nkm);
   fChain->SetBranchAddress("nik0", &nik0, &b_nk0);
   fChain->SetBranchAddress("niem", &niem, &b_niem);
   fChain->SetBranchAddress("niother", &niother, &b_niother);
   fChain->SetBranchAddress("ni", &ni, &b_ni);
   fChain->SetBranchAddress("pdgi", pdgi, &b_pdgi);
   fChain->SetBranchAddress("resc", resc, &b_resc);
   fChain->SetBranchAddress("Ei", Ei, &b_Ei);
   fChain->SetBranchAddress("pxi", pxi, &b_pxi);
   fChain->SetBranchAddress("pyi", pyi, &b_pyi);
   fChain->SetBranchAddress("pzi", pzi, &b_pzi);
   fChain->SetBranchAddress("nf", &nf, &b_nf);
   fChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
   fChain->SetBranchAddress("Ef", Ef, &b_Ef);
   fChain->SetBranchAddress("pxf", pxf, &b_pxf);
   fChain->SetBranchAddress("pyf", pyf, &b_pyf);
   fChain->SetBranchAddress("pzf", pzf, &b_pzf);
   fChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
   fChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
   fChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
   fChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);
   fChain->SetBranchAddress("calresp0", &calresp0, &b_calresp0);
}


#endif  // GST_H
