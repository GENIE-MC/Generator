//____________________________________________________________________________
/*!

\class   HadPlots

\brief   Class to analyze mc events

\author  Tingjun Yang (Stanford Univ)

\created Feb 28, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HADPLOTS_H_
#define _HADPLOTS_H_

#include <vector>
#include <string>

class TGraph;
class TGraphErrors;
class TH1D;
class TProfile;

using namespace std;

namespace genie {
namespace vld_hadronization {

class HadPlots {
      
public:
  
  HadPlots() {};

  HadPlots(string name);

  virtual ~HadPlots() { }
  
  void SetName(string name);
  void LoadData(string mcfile);
  void Analyze();

  TGraphErrors *nch_w[4];         //charged multiplicity vs w2
  TGraph *D_nneg[4];              //dispersion vs N-
  TGraph *D_W2[4];                //dispersion vs w2
  TGraphErrors *npi0_w[4];        //pi0 multiplicity vs w2
  TGraph *nchpi_w[4];             //charged pion multiplicity vs w2
  TGraph *D_npi0[4];              //dispersion vs n for pi0
  TGraph *pnch_w2[13][4];         //topological crosssection
  TH1D *kno[5][4];                //KNO distributions
  TProfile *npi0_nneg[5][4];      //Npi0 vs N-
  TProfile *npi0_nch[3][4];       //Npi0 vs Nch
  TProfile *npi0_nm[4];
  TProfile *npi0_nm_lo[4];
  TProfile *npi0_nm_hi[4];
  TGraph *nch_w_f[4];             //F/B multiplicity
  TGraph *nch_w_b[4];
  TGraph *npos_w_f[4];
  TGraph *nneg_w_f[4];
  TGraph *npos_w_b[4];
  TGraph *nneg_w_b[4];
  TGraph *npos_w2_f[4];
  TGraph *nneg_w2_f[4];
  TGraph *npos_w2_b[4];
  TGraph *nneg_w2_b[4];
  TH1D *xf_pip[4];                //xF distributions
  TH1D *xf_pim[4];
  TH1D *fxf_pip[4];
  TH1D *fxf_pim[4];
  TH1D *Fxf_pos1[4];
  TH1D *Fxf_pro1[4];
  TH1D *Fxf_pos_kno1[4];
  TH1D *Fxf_pos2[4];
  TH1D *Fxf_neg1[4];
  TH1D *Fxf_neg_kno1[4];
  TH1D *Fxf_neg2[4];
  TH1D *Fxf_pos1_hi[4];
  TH1D *Fxf_pro1_hi[4];
  TH1D *Fxf_pos2_hi[4];
  TH1D *Fxf_neg1_hi[4];
  TH1D *Fxf_neg2_hi[4];  
  TH1D *z_pos[4];                //z distributions
  TH1D *z_neg[4];
  TH1D *z_E1[4];
  TH1D *z_E2[4];
  TH1D *z_E3[4];
  TH1D *z1_pos[4];               //nubar
  TH1D *z1_pro[4];
  TH1D *z1_neg[4];
  TH1D *z2_pos[4];               //nubar
  TH1D *z2_pro[4];
  TH1D *z2_neg[4];
  TH1D *pt2_xf1[4];              //Pt distributions
  TH1D *pt2_xf2[4];
  TH1D *pt2_xf3[4];
  TH1D *pt2_xf4[4];
  TH1D *pt2_pip[4];
  TH1D *pt2_pim[4];
  TProfile *pt2_W2_F[4];
  TProfile *pt2_W2_B[4];
  TProfile *pt_W_F[4];
  TProfile *pt_W_B[4];
  TProfile *pt_W[4];
  TProfile *pt2_xf_loW[4];       //seagull plots
  TProfile *pt2_xf_hiW[4];

  string modelName;

private:

  vector<string> mcFiles;

  void BookHists();

};

} //vld_hadronization
} //genie

#endif
