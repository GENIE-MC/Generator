//____________________________________________________________________________
/*!

  \class    genie::BergerSehgalFMCOHPiPXSec2015

  \brief    Computes the triple differential cross section for CC & NC coherent
  pion production according to the \b Berger-Sehgal model.
  v(vbar)A->v(vbar)Api0, vA->l-Api+, vbarA->l+Api-

  Is a concrete implementation of the XSecAlgorithmI interface.

  \ref  PRD 79, 053003 (2009) by Berger and Sehgal  



  \author  G. Perdue, H. Gallagher, D. Cherdack 

  \created 2014 

  \cpright  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE
  */
//____________________________________________________________________________

#ifndef _BERGER_SEHGAL_FM_COHPI_PXSEC_2015_H_
#define _BERGER_SEHGAL_FM_COHPI_PXSEC_2015_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

  class XSecIntegratorI;

  // this class should probably inherit from the infinite mass target version
  class BergerSehgalFMCOHPiPXSec2015 : public XSecAlgorithmI {

    public:
      BergerSehgalFMCOHPiPXSec2015();
      BergerSehgalFMCOHPiPXSec2015(string config);
      virtual ~BergerSehgalFMCOHPiPXSec2015();

      //-- XSecAlgorithmI interface implementation
      double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral        (const Interaction * i) const;
      bool   ValidProcess    (const Interaction * i) const;

      //-- overload the Algorithm::Configure() methods to load private data
      //   members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:
      void LoadConfig(void);

      double ExactKinematicTerm(const Interaction * i) const;
      double PionCOMAbsMomentum(const Interaction * i) const;

      //-- private data members loaded from config Registry or set to defaults
      double fMa;          ///< axial mass
      double fRo;          ///< nuclear size scale parameter
      double fCos8c2;      ///< cos^2(Cabibbo angle)
      bool fRSPionXSec;    ///< Use Rein-Sehgal "style" pion-nucleon xsecs

      const XSecIntegratorI * fXSecIntegrator;
  };

}       // genie namespace

#endif  // _BERGER_SEHGAL_FM_COHPI_PXSEC_2015_H_
