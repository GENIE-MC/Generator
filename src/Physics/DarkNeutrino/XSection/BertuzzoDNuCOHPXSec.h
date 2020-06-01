//____________________________________________________________________________
/*!
  \class    genie::DNuPXSec

  \brief    Differential cross section for v+As coherent elastic scattering.Coherent DNu. \n
            Is a concrete implementation of the XSecAlgorithmI interface. \n

  \ref      E.Bertuzzo, S.Jana, P.A.N.Machado, R. Zukanovich Funcal
            PhysRevLett.121.241801 (2018)

  \author   Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
            University of Sussex

            Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
            University of Liverpool & STFC Rutherford Appleton Laboratory

  \created  September 4, 2017

  \cpright  Copyright (c) 2003-2020, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _BERTUZZO_DNu_COH_CROSS_SECTION_H_
#define _BERTUZZO_DNu_COH_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class BertuzzoDNuCOHPXSec : public XSecAlgorithmI {

public:
  BertuzzoDNuCOHPXSec();
  BertuzzoDNuCOHPXSec(string config);
  virtual ~BertuzzoDNuCOHPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const override;
  double Integral        (const Interaction * i) const override;
  bool   ValidProcess    (const Interaction * i) const override;
  bool   ValidKinematics (const Interaction * i,
                          const double DNu_energy,
                          const double DNu_mass2,
                          const double E,
                          const double M) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config) override;
  void Configure (string param_set) override;

private:

  void LoadConfig(void);

  // // Calculate nuclear density moments
  // double NuclearDensityMoment(int A, int k) const;

  const XSecIntegratorI * fXSecIntegrator;  ///< cross section integrator
  // double fSin2thw;                          ///< sin^2(weinberg angle)

  // // Parameters used for the numerical integration yielding nuclear density moments
  // double fNuclDensMomentCalc_UpperIntegrationLimit; // upper integration limit (in units of nuclear radii R0*A^1/3)
  // double fNuclDensMomentCalc_RelativeTolerance;     // relative tolerance for numerical integrator
  // double fNuclDensMomentCalc_AbsoluteTolerance;     // absolute tolerance for numerical integrator
  // int    fNuclDensMomentCalc_MaxNumOfEvaluations;   // maximum number of integran evaluations in numerical integration

};

}       // genie namespace


#endif // _BERTUZZO_DNu_COH_CROSS_SECTION_H_
