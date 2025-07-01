//____________________________________________________________________________
/*!

\namespace genie::utils::mec

\brief     MEC utilities

\author

\created

\cpright   Copyright (c) 2003-2025, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MEC_UTILS_H_
#define _MEC_UTILS_H_

#include "Physics/HadronTensors/HadronTensorI.h"
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Framework/Numerical/MathUtils.h"

namespace genie {

class Interaction;
class XSecAlgorithmI;

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


  //Version of the tesor contraction in GENIE 2.12.X (modified to use new hadron tensor pool) for debugging purposes
  //double OldTensorContraction(int nupdg, int targetpdg, double Enu, double Ml, double Tl, double costhl, int tensorpdg, genie::HadronTensorType_t tensor_type, char* tensor_model);
  double OldTensorContraction(int nupdg, int targetpdg, double Enu, double Ml, double Tl, double costhl, int tensorpdg, genie::HadronTensorType_t tensor_type, char* tensor_model );

  // Performs a brute-force scan of the kPSTlctl phase space to compute the
  // maximum value of the differential cross section within a specified
  // tolerance. An optional safety factor can be applied to the final result.
  // This function is used by MECGenerator::SelectSuSALeptonKinematics() during
  // rejection sampling. -- S. Gardiner, 16 March 2020
  double GetMaxXSecTlctl( const XSecAlgorithmI& xsec_model,
    const Interaction& inter, const double tolerance = 0.01,
    const double safety_factor = 1.2, const int max_n_layers = 100 );

  // Hard-coded bounds for the scan over phase space in GetMaxXSecTlctl.
  // These are based on the upper limits in the current set of SuSAv2 hadron
  // tensor tables. In the future, it would be good to refactor GetMaxXSecTlctl
  // to retrieve these bounds from the relevant hadron tensor model itself.
  // - S. Gardiner

  // Maximum tabulated energy transfer (GeV)
  const double Q0LimitMaxXSec = 2.;

  // Maximum tabulated momentum transfer (GeV)
  const double QMagLimitMaxXSec = 2.;

  namespace gsl {
    //_____________________________________________________________________________________
    //
    // GSL wrappers
    //
    //_____________________________________________________________________________________
    
    class d2Xsec_dTCosth: public ROOT::Math::IBaseFunctionMultiDim
    {
    public:
      d2Xsec_dTCosth(const XSecAlgorithmI * m, const Interaction & i, const double Enu, 
		     const double LepMass, const double Factor = 1. );
      ~d2Xsec_dTCosth();
      // ROOT::Math::IBaseFunctionMultiDim interface
      unsigned int                        NDim   (void)               const;
      double                              DoEval (const double * xin) const;
      ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
      
    private:
      const XSecAlgorithmI * fModel;
      const Interaction fInteraction;
      const double fEnu ; 
      const double fLepMass ;
      const double fFactor;
    };
    
  }
} // mec   namespace
} // utils namespace
} // genie namespace

#endif // _MEC_UTILS_H_
