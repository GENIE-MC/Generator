//____________________________________________________________________________
/*!

\namespace genie::utils::mec

\brief     MEC utilities

\author    

\created   

\cpright   Copyright (c) 2003-2019, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEC_UTILS_H_
#define _MEC_UTILS_H_

#include "Physics/Multinucleon/XSection/MECHadronTensor.h"

namespace genie {

class Interaction;

namespace utils {
namespace mec   {

  // ---------------------- this should be removed (is replaced by code below)
  //
  // Kinematic calculation: Give q0q3 (and Enu and lepton mass) and 
  // return muon KE and costheta, and jacobian area.
  // Contributed by R.Gran.
  double GetTmuCostFromq0q3(
    double dq0, double dq3, double Enu, double lmass, double &tmu, double &cost, double &area);
  // -----------------------


  //----------------------- once in trunk, this could be copied in KineUtils
  // Kinematic calculations:
  // Get lepton KE and costheta from q0, q3 (and Enu and lepton mass) 
  bool GetTlCostlFromq0q3(
    double q0, double q3, double Enu, double ml, double & Tl, double & costl);
  // Get q0, q3 from lepton KE and costheta (and Enu and lepton mass) 
  bool Getq0q3FromTlCostl(
    double Tl, double costl, double Enu, double ml, double & q0, double & q3);
  //----------------------- once in trunk, this could be embedded in KineUtils::Jacobian()
  // Jacobian for tranformation of d2sigma / dT dcostl to d2sigma / dq0 dq3
  double J(double q0, double q3, double Enu, double ml);
  //----------------------


  // Utility that encodes the Qvalues for the kinematic calculation.
  // this is used in the code that contracts the hadron tensor with the lepton tensor
  // Contributed by R.Gran.
  double Qvalue(int targetpdg, int nupdg);

  // This method implements the lepton tensor contraction with a hadron
  // tensor provided in tabular form.
  // The lepton tensor is expressed formally in Nieves PRC 70 (2004) 055503.
  // Returns  d2sigma/(dTmu dcos_mu) in 10^{41} cm^2 / GeV 
  // Contributed by R.Gran.
  double TensorContraction(
     const Interaction * interaction, 
     int tensor_pdg, 
     MECHadronTensor::MECHadronTensorType_t tensor_type);
  double TensorContraction(
     int nu_pdg, int target_pdg, 
     double Enu,  // neutrino energy (GeV)
     double M_l, double T_l, double costh_l, // f/s lepton mass (GeV), kinetic energy (GeV) and cos(angle)
     int tensor_pdg, 
     MECHadronTensor::MECHadronTensorType_t tensor_type);

} // mec   namespace
} // utils namespace
} // genie namespace

#endif // _MEC_UTILS_H_
