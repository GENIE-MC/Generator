//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenCommon

\brief    Encapsulation of NeuGEN's Physics Variables common block

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_COMMON_H_
#define _NEUGEN_COMMON_H_

namespace genie   {
namespace nuvld   {
namespace facades {

struct physconst_ {

  double gv_fermi;
  double gv_cabibbo;
  double gv_mv;
  double gv_ma_qel;
  double gv_ma_res;
  double gv_ma_coh;
  double gv_fa0;
  double gv_mup;
  double gv_mun;
  double gv_s2thetaw;
  double gv_eta;
  double gv_omega;
  double gv_rsz;
  double gv_fpi;
  double gv_r0;
  double gv_reim;
  double gv_hbarcsq;
  double gv_pi;
  int iok_physconst;
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif
