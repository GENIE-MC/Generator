//____________________________________________________________________________
/*!

\class    genie::mueloss::KokoulinPetrukhinModel

\brief    Kokoulin-Petrukhin model for the energy loss of muons due to direct
          e+e- pair production.
          Concrete implementation of the MuELossI interface.

\ref      W.Lohmann, R.Kopp and R.Voss,
          Energy Loss of Muons in the Energy Range 1-10000 GeV, CERN 85-03

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  December 10, 2003

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KOKOULIN_PETRUKHIN_MODEL_H_
#define _KOKOULIN_PETRUKHIN_MODEL_H_

#include <Math/IFunction.h>

#include "Physics/MuonEnergyLoss/MuELossI.h"
//#include "Numerical/GSFunc.h"

namespace genie {

//class IntegratorI;

namespace mueloss {

class KokoulinPetrukhinModel : public MuELossI
{
public:

  KokoulinPetrukhinModel();
  KokoulinPetrukhinModel(string config);
  virtual ~KokoulinPetrukhinModel();

  //! Implement the MuELossI interface
  double        dE_dx   (double E, MuELMaterial_t material) const;
  MuELProcess_t Process (void) const { return eMupPairProduction; }

//  //! overload the Algorithm::Configure() methods to load private data
//  //! members from configuration options
//  void Configure(const Registry & config);
//  void Configure(string config);
//
//private:
//  void LoadConfig (void);
// // const IntegratorI * fIntegrator;
};

} // mueloss namespace
} // genie   namespace

//____________________________________________________________________________
/*!
\class    genie::mueloss::KokoulinPetrukhinIntegrand

\brief    Auxiliary scalar function for the internal integration in Kokulin
          Petrukhin model

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  December 10, 2003
*/
//____________________________________________________________________________

namespace genie {
 namespace mueloss {
   namespace gsl {

    class KokoulinPetrukhinIntegrand : public ROOT::Math::IBaseFunctionMultiDim
    {
     public:
       KokoulinPetrukhinIntegrand(double E, double Z);
      ~KokoulinPetrukhinIntegrand();
       // ROOT::Math::IBaseFunctionMultiDim interface
       unsigned int                        NDim   (void)               const;
       double                              DoEval (const double * xin) const;
       ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
     private:
       double fE;
       double fZ;
     };

  }  // gsl namespace  
 }  // mueloss namespace
}  // genie   namespace

#endif  // _KOKOULIN_PETRUKHIN_MODEL_
