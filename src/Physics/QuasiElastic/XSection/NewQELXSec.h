//____________________________________________________________________________
/*!

\class    genie::NewQELXSec

\brief    Computes the Quasi Elastic (QEL) total cross section. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  February 26, 2019

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NEW_QEL_XSEC_H_
#define _NEW_QEL_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

#include "TMath.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

namespace genie {

class NuclearModelI;
class VertexGenerator;

namespace utils {
  namespace gsl   {

    class FullQELdXSec : public ROOT::Math::IBaseFunctionMultiDim
    {
     public:
       FullQELdXSec(const XSecAlgorithmI* xsec_model, const Interaction* interaction,
         QELEvGen_BindingMode_t binding_mode, double min_angle_EM);
       virtual ~FullQELdXSec();

       // ROOT::Math::IBaseFunctionMultiDim interface
       unsigned int NDim(void) const;
       double DoEval(const double* xin) const;
       ROOT::Math::IBaseFunctionMultiDim* Clone(void) const;

       Interaction* GetInteractionPtr();
       const Interaction& GetInteraction() const;

     private:
       const XSecAlgorithmI* fXSecModel;
       const NuclearModelI* fNuclModel;
       Interaction* fInteraction;
       QELEvGen_BindingMode_t fHitNucleonBindingMode;
       double fMinAngleEM;
    };

  } // gsl   namespace
} // utils namespace

class NewQELXSec : public XSecIntegratorI {

public:

  NewQELXSec(void);
  NewQELXSec(std::string config);

  /// XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI* model, const Interaction* i) const;

  /// Overload the Algorithm::Configure() methods to load private data
  /// members from configuration options
  void Configure(const Registry& config);
  void Configure(std::string config);

private:

  void LoadConfig (void);

  // Configuration obtained from cross section model
  //QELEvGen_BindingMode_t fBindingMode;

  // XML configuration parameters
  std::string fGSLIntgType;
  double fGSLRelTol;
  unsigned int fGSLMaxEval;
  AlgId fVertexGenID;
  int fNumNucleonThrows;
  double fMinAngleEM;

  // If false, the total cross section will be computed by integrating over
  // lepton scattering angles while preserving the momentum and removal energy
  // of the initial hit nucleon (specified in the input Interaction object).
  // This approach is needed for "CCQE shape" tweak dials in Reweight. If true
  // (default), then an MC integration is also performed over the distribution
  // of initial nucleons. This approach is needed to create total cross section
  // splines.
  bool fAverageOverNucleons;
};


} // genie namespace

#endif  // _NEW_QEL_XSEC_H_
