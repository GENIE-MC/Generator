//____________________________________________________________________________
/*!

\class    genie::StrumiaVissaniIBDPXSec

\brief    An implementation of the neutrino - (free) nucleon [inverse beta
          decay] cross section, valid from the threshold energy (1.806MeV)
          up to hundreds of MeV. Currently cut off at 1/2 nucleon mass.
	  Based on the Strumia/Vissani paper Phys.Lett.B564:42-54,2003

\ref      Strumia A., Vissani F., Phys. Lett. B564, pp42-54 (2003)

\author   Corey Reed <cjreed \at nikhef.nl>
          Nikhef

\created  June 22, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SV_QUASIELASTIC_NU_NUCLEON_XSEC_H_
#define _SV_QUASIELASTIC_NU_NUCLEON_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class StrumiaVissaniIBDPXSec : public XSecAlgorithmI {

public:
  StrumiaVissaniIBDPXSec();
  StrumiaVissaniIBDPXSec(string config);
  virtual ~StrumiaVissaniIBDPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
  
  //-- Routines for calculating the scattering amplitute
  double dSigDt(const double sMinusU,
		const double sMinusMnuc,
		const double t) const;
  double MtxElm(const double sMinusU,
		const double t) const;
  static double MAterm(const double t,
		       const double t2,
		       const double f124,
		       const double f22,
		       const double g124,
		       const double g224meM2,
		       const double f1cf2R8,
		       const double g1cg2R16me,
		       const double g1cFsumR);
  static double MBterm(const double t,
		       const double f1cf2,
		       const double g1cg2,
		       const double g1cFsumR,
		       const double f22);
  static double MCterm(const double t,
		       const double f124,
		       const double f22,
		       const double g124);
  double RadiativeCorr(const double Ee) const;
  double FinalStateCorr(const double Ee) const;

  // variables
  double                  fCosCabibbo2;    //  cos^2(cabibbo)
  double                  fg1of0;          //  axial form factor at q2=0
  double                  fMa2;            //  axial mass squared
  double                  fMv2;            //  vector mass squared
  double                  fNucleonMMDiff;  //  nucleon magnetic moment difference
  double                  fEpsilon;        //  small number used to compare floats with 0

  const XSecIntegratorI * fXSecIntegrator; //! the integrator to get total xsec

public:      
  ClassDef(StrumiaVissaniIBDPXSec, 1) // Inverse Beta Decay partial cross section calculation based on a Strumia/Vissani paper
};
}       // genie namespace
#endif  // _SV_QUASIELASTIC_NU_NUCLEON_XSEC_H_
