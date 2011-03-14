//____________________________________________________________________________
/*!

\class   HadPlots

\brief   Class to analyze mc events

\author  Tingjun Yang (Stanford Univ)

\created Feb 28, 2009

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
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

const int kMaxFiles=40; 

class HadPlots {
      
public:  
  HadPlots();
  HadPlots(string name);
 ~HadPlots();
  
  void SetName  (string name);
  void LoadData (string mcfile);
  void Analyze  (void);


  string         modelName;
  TGraphErrors * nch_w[kMaxFiles];             ///< charged multiplicity vs w2
  TGraphErrors * nk0_w[kMaxFiles];             ///< neutral kaons vs W2
  TGraphErrors * nlam_w[kMaxFiles];            ///< lambdas vs W2
  TGraph *       D_nneg[kMaxFiles];            ///< dispersion vs N-
  TGraph *       D_W2[kMaxFiles];              ///< dispersion vs w2
  TGraphErrors * npi0_w[kMaxFiles];            ///< pi0 multiplicity vs w2
  TGraph *       nchpi_w[kMaxFiles];           ///< charged pion multiplicity vs w2
  TGraph *       D_npi0[kMaxFiles];            ///< dispersion vs n for pi0
  TGraph *       pnch_w2[13][kMaxFiles];       ///< normalized topological cross section
  TH1D *         kno[5][kMaxFiles];            ///< KNO distributions
  TProfile *     npi0_nneg[5][kMaxFiles];      ///< Npi0 vs N-
  TProfile *     npi0_nch[3][kMaxFiles];       ///< Npi0 vs Nch
  TProfile *     npi0_nm[kMaxFiles];
  TProfile *     npi0_nm_lo[kMaxFiles];
  TProfile *     npi0_nm_hi[kMaxFiles];
  TGraph *       nch_w_f[kMaxFiles];           ///< F/B multiplicity
  TGraph *       nch_w_b[kMaxFiles];
  TGraph *       npos_w_f[kMaxFiles];
  TGraph *       nneg_w_f[kMaxFiles];
  TGraph *       npos_w_b[kMaxFiles];
  TGraph *       nneg_w_b[kMaxFiles];
  TGraph *       npos_w2_f[kMaxFiles];
  TGraph *       nneg_w2_f[kMaxFiles];
  TGraph *       npos_w2_b[kMaxFiles];
  TGraph *       nneg_w2_b[kMaxFiles];
  TH1D *         xf_pip[kMaxFiles];            ///< xF distributions
  TH1D *         xf_pim[kMaxFiles];
  TH1D *         fxf_pip[kMaxFiles];
  TH1D *         fxf_pim[kMaxFiles];
  TH1D *         Fxf_pos1[kMaxFiles];
  TH1D *         Fxf_pro1[kMaxFiles];
  TH1D *         Fxf_pos_kno1[kMaxFiles];
  TH1D *         Fxf_pos2[kMaxFiles];
  TH1D *         Fxf_neg1[kMaxFiles];
  TH1D *         Fxf_neg_kno1[kMaxFiles];
  TH1D *         Fxf_neg2[kMaxFiles];
  TH1D *         Fxf_pos1_hi[kMaxFiles];
  TH1D *         Fxf_pro1_hi[kMaxFiles];
  TH1D *         Fxf_pos2_hi[kMaxFiles];
  TH1D *         Fxf_neg1_hi[kMaxFiles];
  TH1D *         Fxf_neg2_hi[kMaxFiles];  
  TH1D *         z_pos[kMaxFiles];              ///< z distributions
  TH1D *         z_neg[kMaxFiles];
  TH1D *         z_E1[kMaxFiles];
  TH1D *         z_E2[kMaxFiles];
  TH1D *         z_E3[kMaxFiles];
  TH1D *         z1_pos[kMaxFiles];             ///< nubar
  TH1D *         z1_pro[kMaxFiles];
  TH1D *         z1_neg[kMaxFiles];
  TH1D *         z2_pos[kMaxFiles];             ///< nubar
  TH1D *         z2_pro[kMaxFiles];
  TH1D *         z2_neg[kMaxFiles];
  TH1D *         pt2_xf1[kMaxFiles];            ///< Pt distributions
  TH1D *         pt2_xf2[kMaxFiles];
  TH1D *         pt2_xf3[kMaxFiles];
  TH1D *         pt2_xf4[kMaxFiles];
  TH1D *         pt2_pip[kMaxFiles];
  TH1D *         pt2_pim[kMaxFiles];
  TProfile *     pt2_W2_F[kMaxFiles];
  TProfile *     pt2_W2_B[kMaxFiles];
  TProfile *     pt_W_F[kMaxFiles];
  TProfile *     pt_W_B[kMaxFiles];
  TProfile *     pt_W[kMaxFiles];
  TProfile *     pt2_xf_loW[kMaxFiles];         ///< seagull plots
  TProfile *     pt2_xf_hiW[kMaxFiles];

private:

  vector<string> mcFiles;

  void BookHists();

};

} //vld_hadronization
} //genie

#endif
