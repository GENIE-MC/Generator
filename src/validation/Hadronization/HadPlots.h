//____________________________________________________________________________
/*!

\class   HadPlots

\brief   Class to analyze mc events

\author  Tingjun Yang (Stanford Univ)

\created Feb 28, 2009

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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
namespace mc_vs_data {

const int kMaxFiles=40; 

class HadPlots {
      
public:  
  HadPlots();
  HadPlots(string name);
 ~HadPlots();
  
  void SetName  (string name);
  void LoadData (string mcfile);
  void Analyze  (void);

  string         modelName;                    ///< name for the current model shown on plot legends

  TGraphErrors * nch_w[kMaxFiles];             ///< average charged hadon multiplicity, <nch>, vs W^2 (p. 1, 13)
  TGraph *       D_nneg[kMaxFiles];            ///< negative hadron dispersion D- vs average multiplicity <n-> (p. 2, 14)
  TGraph *       D_W2[kMaxFiles];              ///< charged hadron dispersion over average multiplicity, D/<nch>, vs W^2 (p. 2)
  TH1D *         kno[5][kMaxFiles];            ///< KNO distributions (p. 3)
  TGraphErrors * npi0_w[kMaxFiles];            ///< average pi0 multiplicity <npi0> vs W^2 (p. 4, 16)
  TGraph *       D_npi0[kMaxFiles];            ///< pi0 dispersion, Dpi0, vs average multiplicty, <npi0> (p. 4)
  TProfile *     npi0_nneg[5][kMaxFiles];      ///< <npi0> vs n-, for various W ranges (p. 5)
  TGraph *       pnch_w2[13][kMaxFiles];       ///< normalized topological cross section, for various multiplicities (p. 6, 15)
  TGraph *       nch_w_f[kMaxFiles];           ///< charged hadron multiplicity <nch> vs W^2, xF > 0 (p. 7, 17)
  TGraph *       nch_w_b[kMaxFiles];           ///< charged hadron multiplicity <nch> vs W^2, xF < 0 (p. 7, 17)
  TH1D *         xf_pip[kMaxFiles];            ///< xF distribution, pi+ (p. 8)
  TH1D *         xf_pim[kMaxFiles];            ///< xF distribution, pi- (p. 8)
  TH1D *         z_pos[kMaxFiles];             ///< z distribution, h+ (p. 9)
  TH1D *         z_neg[kMaxFiles];             ///< z distribution, h- (p. 9)
  TProfile *     pt2_W2_F[kMaxFiles];          ///< <pT^2> vs W2, xF > 0.3 (p. 10)
  TProfile *     pt2_W2_B[kMaxFiles];          ///< <pT^2> vs W2, xF < 0.3 (p. 10)
  TProfile *     pt2_xf_loW[kMaxFiles];        ///< seagull plot, 9 < W^2 < 25 GeV^2 - Q2 > 1 GeV^2, x < 0.95, El > 4 GeV (p. 11)
  TProfile *     neta_W[kMaxFiles];            ///< eta multiplicity vs W (p. 12)
  TProfile *     neta_W_F[kMaxFiles];          ///< eta multiplicity vs W for xF > 0 (p. 12)
  TH1D *         Fxf_pos1[kMaxFiles];          ///< xF distribution, positive hadrons, 2 < W <  4 GeV (p. 18)
  TH1D *         Fxf_neg1[kMaxFiles];          ///< xF distribution, negative hadrons, 2 < W <  4 GeV (p. 18)
  TH1D *         Fxf_pos1_hi[kMaxFiles];       ///< xF distribution, positive hadrons, 4 < W < 10 GeV (p. 18)
  TH1D *         Fxf_neg1_hi[kMaxFiles];       ///< xF distribution, negative hadrons, 4 < W < 10 GeV (p. 18)
  TH1D *         fxf_pip[kMaxFiles];           ///< xF distribution, pi+ (p. 19)
  TH1D *         fxf_pim[kMaxFiles];           ///< xF distribution, pi- (p. 19)
  TH1D *         z_E1[kMaxFiles];              ///< z distribution, Ev < 15 GeV (p. 20)
  TH1D *         z_E2[kMaxFiles];              ///< z distribution, 15 < Ev < 30 GeV (p. 20)
  TH1D *         z_E3[kMaxFiles];              ///< z distribution, Ev > 30 GeV (p. 20)
  TProfile *     pt_W_F[kMaxFiles];            ///< <pT> vs W, xF > 0 (p. 21)
  TProfile *     pt_W_B[kMaxFiles];            ///< <pT> vs W, xF < 0 (p. 21)
  TProfile *     pt_W[kMaxFiles];              ///< <pT> vs W, all xF (p. 21)

  //
  // plots filled, but not shown
  //
  TProfile *     pt2_xf_hiW[kMaxFiles];        ///< seagull plot, W^2 > 40 GeV^2 (relevant plot at p. 11)
  TGraph *       nchpi_w[kMaxFiles];           ///< charged pion multiplicity vs W^2
  TGraph *       npos_w_f[kMaxFiles];
  TGraph *       nneg_w_f[kMaxFiles];
  TGraph *       npos_w_b[kMaxFiles];
  TGraph *       nneg_w_b[kMaxFiles];
  TGraph *       npos_w2_f[kMaxFiles];
  TGraph *       nneg_w2_f[kMaxFiles];
  TGraph *       npos_w2_b[kMaxFiles];
  TGraph *       nneg_w2_b[kMaxFiles];
  TProfile *     npi0_nch[3][kMaxFiles];       ///< npi0 vs nch
  TProfile *     npi0_nm[kMaxFiles];           ///<
  TProfile *     npi0_nm_lo[kMaxFiles];        ///<
  TProfile *     npi0_nm_hi[kMaxFiles];        ///<
  TH1D *         Fxf_pro1[kMaxFiles];
  TH1D *         Fxf_pos_kno1[kMaxFiles];
  TH1D *         Fxf_pos2[kMaxFiles];
  TH1D *         Fxf_neg_kno1[kMaxFiles];
  TH1D *         Fxf_neg2[kMaxFiles];
  TH1D *         Fxf_pro1_hi[kMaxFiles];
  TH1D *         Fxf_pos2_hi[kMaxFiles];
  TH1D *         Fxf_neg2_hi[kMaxFiles];  
  TH1D *         z1_pos[kMaxFiles];    
  TH1D *         z1_pro[kMaxFiles];
  TH1D *         z1_neg[kMaxFiles];
  TH1D *         z2_pos[kMaxFiles];
  TH1D *         z2_pro[kMaxFiles];
  TH1D *         z2_neg[kMaxFiles];
  TH1D *         pt2_xf1[kMaxFiles];  
  TH1D *         pt2_xf2[kMaxFiles];
  TH1D *         pt2_xf3[kMaxFiles];
  TH1D *         pt2_xf4[kMaxFiles];
  TH1D *         pt2_pip[kMaxFiles];
  TH1D *         pt2_pim[kMaxFiles];

/*
  //
  // unused plots: not filled and not shown
  //
  TGraphErrors * nk0_w[kMaxFiles];             ///< neutral kaons vs W^2
  TGraphErrors * nlam_w[kMaxFiles];            ///< lambdas vs W^2
*/

private:

  vector<string> mcFiles;

  void BookHists();

};

} //mc_vs_data
} //genie

#endif
