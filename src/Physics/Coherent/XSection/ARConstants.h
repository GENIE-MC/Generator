//____________________________________________________________________________
/*!

\class    genie::alvarezruso::Constants

\brief    Class containing constants for AlvarezRuso coherent pion production xsec

\ref

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  05/12/2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef AR_CONSTANTS_H
#define AR_CONSTANTS_H

#include "Framework/Registry/Registry.h"
#include "Framework/Conventions/Constants.h"

namespace genie {
namespace alvarezruso {

class ARConstants
{
  public:

    ARConstants();
    ~ARConstants();

    double HBar();
    double Ma_Nucleon();
    double Mv_Nucleon();
    double Ma_Delta();
    double Mv_Delta();
    double GAxial();
    double Rho0();
    double CA4_A();
    double CA5_A();
    double CA4_B();
    double CA5_B();
    double PiDecayConst();
    double DeltaNCoupling();
    double CosCabibboAngle();
    double SinWeinbergAngle();
    double GFermi();
    double ElectronMass();
    double MuonMass();
    double TauMass();
    double ProtonMass();
    double NeutronMass();
    double NucleonMass();
    double NucleonMassSq();
    double DeltaPMass();
    double Delta0Mass();
    double PiPMass();
    double Pi0Mass();
    double cm38Conversion();

    double NCFactor();

  private:
    // unused // const genie::Registry *reg;

    double COHAR_Ma_Nuc      ;
    double COHAR_Mv_Nuc      ;
    double COHAR_Ma_Delta    ;
    double COHAR_Mv_Delta    ;
    double COHAR_GA0         ;
    double COHAR_Rho0        ;
    double COHAR_a4          ;
    double COHAR_a5          ;
    double COHAR_b4          ;
    double COHAR_b5          ;
    double COHAR_fPi_byHbar  ;
    double COHAR_fStar       ;
    double fCosCabibboAngle  ;
    double fSinWeinbergAngle ;

    double massElectron      ;
    double massMuon          ;
    double massTau           ;
    double massProton        ;
    double massNeutron       ;
    double massNucleon       ;
    double massNucleon2      ;
    double massDeltaP        ;
    double massDelta0        ;
    double massPiP           ;
    double massPi0           ;

    double ncFactor;

};

} //namespace alvarezruso

} //namespace genie
#endif
