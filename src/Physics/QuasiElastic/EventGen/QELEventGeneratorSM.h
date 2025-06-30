//____________________________________________________________________________
/*!

\class    genie::QELEventGeneratorSM

\brief    Generates values for the kinematic variables describing QEL neutrino
          interaction events for Smith-Moniz model.
          Is a concrete implementation of the EventRecordVisitorI interface.

\ref      [1] R.A.Smith and E.J.Moniz, Nuclear Physics  B43, (1972) 605-622 \n
          [2] K.S. Kuzmin, V.V. Lyubushkin, V.A.Naumov Eur. Phys. J. C54, (2008) 517-538

\author   Igor Kakorin <kakorin@jinr.ru> \n
          Joint Institute for Nuclear Research \n

          adapted from  fortran code provided by:

          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>, \n
          Joint Institute for Nuclear Research,
          Institute for Theoretical and Experimental Physics \n

          Vadim Naumov <vnaumov@theor.jinr.ru>, \n
          Joint Institute for Nuclear Research  \n

          based on code of:
          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 05, 2017

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _QEL_EVENT_GENERATORSM_H_
#define _QEL_EVENT_GENERATORSM_H_

#include <Math/IntegratorMultiDim.h>

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"

namespace genie {

class QELEventGeneratorSM: public KineGeneratorWithCache {

public :
  QELEventGeneratorSM();
  QELEventGeneratorSM(string config);
 ~QELEventGeneratorSM();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  mutable SmithMonizUtils * sm_utils;

  void   LoadConfig     (void);
  double ComputeMaxXSec(const Interaction * in) const;
  double ComputeMaxXSec (const Interaction * in, const int nkey) const;
  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant

  

  mutable KinePhaseSpace_t fkps;


  bool fGenerateNucleonInNucleus;           ///< generate struck nucleon in nucleus
  double fQ2Min;                            ///< Q2-threshold for seeking the second maximum


}; // class definition

class XSecAlgorithmI;
class Interaction;

namespace utils {
namespace gsl   {
//.....................................................................................
//
// genie::utils::gsl::d3XSecSM_dQ2dvdkF_E
// A 3-D cross section function: d3XSecSM_dQ2dvdkF_E = f(Q2, v, kF=fixed)|(fixed E)
//
class d3XSecSM_dQ2dvdkF_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d3XSecSM_dQ2dvdkF_E(const XSecAlgorithmI *, const Interaction *, double pF);
 ~d3XSecSM_dQ2dvdkF_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double *)     const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction    * fInteraction;
  const double fpF;
};
//
// genie::utils::gsl::d1XSecSM_dQ2_E
// A 1-D cross section function: d1XSecSM_dQ2_E = f(Q2)|(fixed E)
//
class d1XSecSM_dQ2_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d1XSecSM_dQ2_E(const XSecAlgorithmI *, const Interaction *);
 ~d1XSecSM_dQ2_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double *)     const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction    * fInteraction;
};
//
// genie::utils::gsl::dv_dQ2_E=f(Q2)|(fixed E)
// A 1-D dependence of allowable \nu-range from Q2
//
class dv_dQ2_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  dv_dQ2_E(const Interaction *);
 ~dv_dQ2_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double *)     const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const Interaction    * fInteraction;
  mutable SmithMonizUtils * sm_utils;
};
} // gsl   namespace
} // utils namespace


} // genie namespace

#endif // _QEL_EVENT_GENERATORSM_H_
