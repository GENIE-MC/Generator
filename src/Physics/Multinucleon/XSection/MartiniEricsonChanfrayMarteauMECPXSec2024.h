//____________________________________________________________________________
/*!

\class    genie::MartiniEricsonChanfrayMarteauMECPXSec2024

\brief    Computes the Martini, Ericson, Chanfray and Marteau MEC model
          differential cross section.
          Uses precomputed hadon tensor tables.

\author   Lavinia Russo <lavinia.russo \at lpnhe.in2p3.fr>
          Laboratoire de physique nucléaire et des hautes énergies (LPNHE) - CNRS - Sorbonne University

\ref      L. Russo, M. Martini, S. Dolan, L. Munteanu, B. Popov, C. Giganti
          Implementation of the Martini-Ericson-Chanfray-Marteau RPA-based
          (anti)neutrino cross-section model in the GENIE neutrino event generator
          arXiv:2508.13939 [hep-ex]

\created  Feb 2024

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _MARTINI_ERICSON_CHANFRAY_MARTEAU_MEC_PXSEC_2016_H_
#define _MARTINI_ERICSON_CHANFRAY_MARTEAU_MEC_PXSEC_2016_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HadronTensors/HadronTensorModelI.h"
#include "Physics/Common/XSecScaleI.h"
#include "Physics/Common/QvalueShifter.h"

namespace genie {

class XSecIntegratorI;

class MartiniEricsonChanfrayMarteauMECPXSec2024 : public XSecAlgorithmI {

public:

  MartiniEricsonChanfrayMarteauMECPXSec2024();
  MartiniEricsonChanfrayMarteauMECPXSec2024(string config);
  virtual ~MartiniEricsonChanfrayMarteauMECPXSec2024();

  // XSecAlgorithmI interface implementation
  double XSec(const Interaction* i, KinePhaseSpace_t k) const;
  double Integral(const Interaction* i) const;
  bool   ValidProcess(const Interaction* i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

  // Method specifically for evaluating np/pp pair probabilities
  double PairRatio(const Interaction* i,
    const std::string& final_state_ratio = "pnFraction") const;

private:

  /// Load algorithm configuration
  void LoadConfig (void);

  // Calculate Qvalue Shift for susa:
  double Qvalue(const Interaction & interaction ) const ;

  /// External scaling factor for this cross section
  double fXSecScale;

  const genie::HadronTensorModelI* fHadronTensorModel;

  // Fermi momentum table used for scaling
  string fKFTable;

  // Binding energies:
  double fEbHe;
  double fEbLi;
  double fEbC;
  double fEbO;
  double fEbMg;
  double fEbAr;
  double fEbCa;
  double fEbFe;
  double fEbNi;
  double fEbSn;
  double fEbAu;
  double fEbPb;

  /// GSL numerical integrator
  const XSecIntegratorI*  fXSecIntegrator;

  const XSecScaleI * fMECScaleAlg ; // Optional algorithm to scale the xsec as a function of W
  const QvalueShifter * fQvalueShifter ; // Optional algorithm to retrieve the qvalue shift for a given target

};

} // genie namespace
#endif // _MARTINI_ERICSON_CHANFRAY_MARTEAU_MEC_PXSEC_2016_H_
