//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the namespace documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 11, 2008 - CA
   Fixed a bug with the Delta- pdg code. It was incorrectly set to -2214.
   Now set to 1114. The bug affected the final state nubar RES events.
 @ Jun 17, 2009 - CA
   Used resonance codes from PDG/PDGCodes.h
 @ Oct 20, 2009 - CA
   Modified ResonanceCharge() to take into account the probe charge (so as
   to conserve charge in charged lepton scattering)
 @ Jul 23, 2010 - CA
   Moved ResonanceCharge(Interaction) to EVGModules/HadronicSystemGenerator 
   to avoid dependency of BaryonResonance package on the Interaction package.
   Added OrbitalAngularMom(Resonance_t), ResonanceIndex(Resonance_t)
   Width(Resonance_t) and BWNorm(Resonance_t) functions, previously available
   through a BaryonResDataSetI implementation. Simplified BaryonResonance
   package by removing the redundant BaryonResDataPDG, BaryonResDataSetI
   BreitWignerI, BreitWignerRes, BreitWignerLRes and BaryonResParams classes.

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <map> 

#include <TMath.h>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/BWFunc.h"

using namespace genie;

//____________________________________________________________________________
const char * genie::utils::res::AsString(Resonance_t res)
{
    switch(res) {
        case kP33_1232  : return "P33(1232)" ; break;
        case kS11_1535  : return "S11(1535)" ; break;
        case kD13_1520  : return "D13(1520)" ; break;
        case kS11_1650  : return "S11(1650)" ; break;
        case kD13_1700  : return "D13(1700)" ; break;
        case kD15_1675  : return "D15(1675)" ; break;
        case kS31_1620  : return "S31(1620)" ; break;
        case kD33_1700  : return "D33(1700)" ; break;
        case kP11_1440  : return "P11(1440)" ; break;
        case kP33_1600  : return "P33(1600)" ; break;
        case kP13_1720  : return "P13(1720)" ; break;
        case kF15_1680  : return "F15(1680)" ; break;
        case kP31_1910  : return "P31(1910)" ; break;
        case kP33_1920  : return "P33(1920)" ; break;
        case kF35_1905  : return "F35(1905)" ; break;
        case kF37_1950  : return "F37(1950)" ; break;
        case kP11_1710  : return "P11(1710)" ; break;
        case kF17_1970  : return "F17(1970)" ; break;
        default: break;
    }
    return "unknown resonance!";
}
//____________________________________________________________________________
Resonance_t genie::utils::res::FromString(const char * res)
{
    if     ( strcmp( res,"P33(1232)" ) == 0 )  return kP33_1232;
    else if( strcmp( res,"S11(1535)" ) == 0 )  return kS11_1535;
    else if( strcmp( res,"D13(1520)" ) == 0 )  return kD13_1520;
    else if( strcmp( res,"S11(1650)" ) == 0 )  return kS11_1650;
    else if( strcmp( res,"D13(1700)" ) == 0 )  return kD13_1700;
    else if( strcmp( res,"D15(1675)" ) == 0 )  return kD15_1675;
    else if( strcmp( res,"S31(1620)" ) == 0 )  return kS31_1620;
    else if( strcmp( res,"D33(1700)" ) == 0 )  return kD33_1700;
    else if( strcmp( res,"P11(1440)" ) == 0 )  return kP11_1440;
    else if( strcmp( res,"P33(1600)" ) == 0 )  return kP33_1600;
    else if( strcmp( res,"P13(1720)" ) == 0 )  return kP13_1720;
    else if( strcmp( res,"F15(1680)" ) == 0 )  return kF15_1680;
    else if( strcmp( res,"P31(1910)" ) == 0 )  return kP31_1910;
    else if( strcmp( res,"P33(1920)" ) == 0 )  return kP33_1920;
    else if( strcmp( res,"F35(1905)" ) == 0 )  return kF35_1905;
    else if( strcmp( res,"F37(1950)" ) == 0 )  return kF37_1950;
    else if( strcmp( res,"P11(1710)" ) == 0 )  return kP11_1710;
    else if( strcmp( res,"F17(1970)" ) == 0 )  return kF17_1970;
    else return kNoResonance;
}
//____________________________________________________________________________
Resonance_t genie::utils::res::FromPdgCode(int pdgc)
{
    // Delta(1232) code from the standard PDG table.
    // Higher resonance codes come from MINOS extensions to PDG tables.

    switch(pdgc) {

        case  (kPdgP33m1232_DeltaM ) : /* Delta-  */
        case  (kPdgP33m1232_Delta0 ) : /* Delta0  */
        case  (kPdgP33m1232_DeltaP ) : /* Delta+  */
        case  (kPdgP33m1232_DeltaPP) : /* Delta++ */
            return kP33_1232; break;

        case (kPdgS11m1535_N0) :       /* N0  */
        case (kPdgS11m1535_NP) :       /* N+  */
            return kS11_1535; break;

        case (kPdgD13m1520_N0) :       /* N0  */
        case (kPdgD13m1520_NP) :       /* N+  */
            return kD13_1520; break;

        case (kPdgS11m1650_N0) :       /* N0  */
        case (kPdgS11m1650_NP) :       /* N+  */
            return kS11_1650; break;

        case (kPdgD13m1700_N0) :       /* N0  */
        case (kPdgD13m1700_NP) :       /* N+  */
            return kD13_1700; break;

        case (kPdgD15m1675_N0) :       /* N0  */
        case (kPdgD15m1675_NP) :       /* N+  */
            return kD15_1675; break;

        case (kPdgS31m1620_DeltaM ) :  /* Delta-  */
        case (kPdgS31m1620_Delta0 ) :  /* Delta0  */
        case (kPdgS31m1620_DeltaP ) :  /* Delta+  */
        case (kPdgS31m1620_DeltaPP) :  /* Delta++ */
            return kS31_1620; break;

        case (kPdgD33m1700_DeltaM ) :  /* Delta-  */
        case (kPdgD33m1700_Delta0 ) :  /* Delta0  */
        case (kPdgD33m1700_DeltaP ) :  /* Delta+  */
        case (kPdgD33m1700_DeltaPP) :  /* Delta++ */
            return kD33_1700; break;

        case (kPdgP11m1440_N0) :       /* N0  */
        case (kPdgP11m1440_NP) :       /* N+  */
            return kP11_1440; break;

        case (kPdgP13m1720_N0) :       /* N0  */
        case (kPdgP13m1720_NP) :       /* N+  */
            return kP13_1720; break;

        case (kPdgF15m1680_N0) :       /* N0  */
        case (kPdgF15m1680_NP) :       /* N+  */
            return kF15_1680; break;

        case (kPdgP31m1910_DeltaM ) :  /* Delta-  */
        case (kPdgP31m1910_Delta0 ) :  /* Delta0  */
        case (kPdgP31m1910_DeltaP ) :  /* Delta+  */
        case (kPdgP31m1910_DeltaPP) :  /* Delta++ */
            return kP31_1910; break;

        case (kPdgP33m1920_DeltaM ) :  /* Delta-  */
        case (kPdgP33m1920_Delta0 ) :  /* Delta0  */
        case (kPdgP33m1920_DeltaP ) :  /* Delta+  */
        case (kPdgP33m1920_DeltaPP) :  /* Delta++ */
            return kP33_1920; break;

        case (kPdgF35m1905_DeltaM ) :  /* Delta-  */
        case (kPdgF35m1905_Delta0 ) :  /* Delta0  */
        case (kPdgF35m1905_DeltaP ) :  /* Delta+  */
        case (kPdgF35m1905_DeltaPP) :  /* Delta++ */
            return kF35_1905; break;

        case (kPdgF37m1950_DeltaM ) :  /* Delta-  */
        case (kPdgF37m1950_Delta0 ) :  /* Delta0  */
        case (kPdgF37m1950_DeltaP ) :  /* Delta+  */
        case (kPdgF37m1950_DeltaPP) :  /* Delta++ */
            return kF37_1950; break;

        case (kPdgP11m1710_N0) :       /* N0  */
        case (kPdgP11m1710_NP) :       /* N+  */
            return kP11_1710; break;
    }

    return kNoResonance;
}
//____________________________________________________________________________
int genie::utils::res::PdgCode(Resonance_t res, int Q)
{
    // Delta(1232) code from the standard PDG table
    // Higher resonance codes come from MINOS extensions to PDG tables

    switch(res) {

        case kP33_1232:
            if(Q == -1) return  kPdgP33m1232_DeltaM;  /* Delta-  */
            if(Q ==  0) return  kPdgP33m1232_Delta0;  /* Delta0  */
            if(Q ==  1) return  kPdgP33m1232_DeltaP;  /* Delta+  */
            if(Q ==  2) return  kPdgP33m1232_DeltaPP; /* Delta++ */
            break;

        case kS11_1535:
            if(Q ==  0) return  kPdgS11m1535_N0; /* N0  */
            if(Q ==  1) return  kPdgS11m1535_NP; /* N+  */
            break;

        case kD13_1520:
            if(Q ==  0) return  kPdgD13m1520_N0; /* N0  */
            if(Q ==  1) return  kPdgD13m1520_NP; /* N+  */
            break;

        case kS11_1650:
            if(Q ==  0) return  kPdgS11m1650_N0; /* N0  */
            if(Q ==  1) return  kPdgS11m1650_NP; /* N+  */
            break;

        case kD13_1700:
            if(Q ==  0) return  kPdgD13m1700_N0; /* N0  */
            if(Q ==  1) return  kPdgD13m1700_NP; /* N+  */
            break;

        case kD15_1675:
            if(Q ==  0) return  kPdgD15m1675_N0; /* N0  */
            if(Q ==  1) return  kPdgD15m1675_NP; /* N+  */
            break;

        case kS31_1620:
            if(Q == -1) return  kPdgS31m1620_DeltaM;  /* Delta-  */
            if(Q ==  0) return  kPdgS31m1620_Delta0;  /* Delta0  */
            if(Q ==  1) return  kPdgS31m1620_DeltaP;  /* Delta+  */
            if(Q ==  2) return  kPdgS31m1620_DeltaPP; /* Delta++ */
            break;

        case kD33_1700:
            if(Q == -1) return  kPdgD33m1700_DeltaM;  /* Delta-  */
            if(Q ==  0) return  kPdgD33m1700_Delta0;  /* Delta0  */
            if(Q ==  1) return  kPdgD33m1700_DeltaP;  /* Delta+  */
            if(Q ==  2) return  kPdgD33m1700_DeltaPP; /* Delta++ */
            break;

        case kP11_1440:
            if(Q ==  0) return  kPdgP11m1440_N0; /* N0  */
            if(Q ==  1) return  kPdgP11m1440_NP; /* N+  */
            break;

        case kP33_1600:
            return 0;   

        case kP13_1720:
            if(Q ==  0) return  kPdgP13m1720_N0; /* N0  */
            if(Q ==  1) return  kPdgP13m1720_NP; /* N+  */
            break;

        case kF15_1680:
            if(Q ==  0) return  kPdgF15m1680_N0; /* N0  */
            if(Q ==  1) return  kPdgF15m1680_NP; /* N+  */
            break;

        case kP31_1910:
            if(Q == -1) return  kPdgP31m1910_DeltaM;  /* Delta-  */
            if(Q ==  0) return  kPdgP31m1910_Delta0;  /* Delta0  */
            if(Q ==  1) return  kPdgP31m1910_DeltaP;  /* Delta+  */
            if(Q ==  2) return  kPdgP31m1910_DeltaPP; /* Delta++ */
            break;

        case kP33_1920:
            if(Q == -1) return  kPdgP33m1920_DeltaM;  /* Delta-  */
            if(Q ==  0) return  kPdgP33m1920_Delta0;  /* Delta0  */
            if(Q ==  1) return  kPdgP33m1920_DeltaP;  /* Delta+  */
            if(Q ==  2) return  kPdgP33m1920_DeltaPP; /* Delta++ */
            break;

        case kF35_1905:
            if(Q == -1) return  kPdgF35m1905_DeltaM;  /* Delta-  */
            if(Q ==  0) return  kPdgF35m1905_Delta0;  /* Delta0  */
            if(Q ==  1) return  kPdgF35m1905_DeltaP;  /* Delta+  */
            if(Q ==  2) return  kPdgF35m1905_DeltaPP; /* Delta++ */
            break;

        case kF37_1950:
            if(Q == -1) return  kPdgF37m1950_DeltaM;  /* Delta-  */
            if(Q ==  0) return  kPdgF37m1950_Delta0;  /* Delta0  */
            if(Q ==  1) return  kPdgF37m1950_DeltaP;  /* Delta+  */
            if(Q ==  2) return  kPdgF37m1950_DeltaPP; /* Delta++ */
            break;

        case kP11_1710:
            if(Q ==  0) return  kPdgP11m1710_N0; /* N0  */
            if(Q ==  1) return  kPdgP11m1710_NP; /* N+  */
            break;

        case kF17_1970:
            return 0;   
            break;

        default:
            return 0;
    }

    return 0;
}
//____________________________________________________________________________
bool genie::utils::res::IsBaryonResonance(int pdgc)
{
    // Delta(1232) code from the standard PDG table.
    // Higher resonance codes come from MINOS extensions to PDG tables.

    switch(pdgc) {

        /* ------ P33(1232) ------*/
        case  (kPdgP33m1232_DeltaM ) : /* Delta-  */
        case  (kPdgP33m1232_Delta0 ) : /* Delta0  */
        case  (kPdgP33m1232_DeltaP ) : /* Delta+  */
        case  (kPdgP33m1232_DeltaPP) : /* Delta++ */

            /* ------ S11(1535) ------*/
        case (kPdgS11m1535_N0) :       /* N0  */
        case (kPdgS11m1535_NP) :       /* N+  */

            /* ------ D13(1520) ------*/
        case (kPdgD13m1520_N0) :       /* N0  */
        case (kPdgD13m1520_NP) :       /* N+  */

            /* ------ S11(1650) ------*/
        case (kPdgS11m1650_N0) :       /* N0  */
        case (kPdgS11m1650_NP) :       /* N+  */

            /* ------ D13(1700) ------*/
        case (kPdgD13m1700_N0) :       /* N0  */
        case (kPdgD13m1700_NP) :       /* N+  */

            /* ------ D15(1675) ------*/
        case (kPdgD15m1675_N0) :       /* N0  */
        case (kPdgD15m1675_NP) :       /* N+  */

            /* ------ S31(1620) ------*/
        case (kPdgS31m1620_DeltaM ) :  /* Delta-  */
        case (kPdgS31m1620_Delta0 ) :  /* Delta0  */
        case (kPdgS31m1620_DeltaP ) :  /* Delta+  */
        case (kPdgS31m1620_DeltaPP) :  /* Delta++ */

            /* ------ D33(1700) ------*/
        case (kPdgD33m1700_DeltaM ) :  /* Delta-  */
        case (kPdgD33m1700_Delta0 ) :  /* Delta0  */
        case (kPdgD33m1700_DeltaP ) :  /* Delta+  */
        case (kPdgD33m1700_DeltaPP) :  /* Delta++ */

            /* ------ P11(1440) ------*/
        case (kPdgP11m1440_N0) :       /* N0  */
        case (kPdgP11m1440_NP) :       /* N+  */

            /* ------ P33(1600) ------*/
            // are you?

            /* ------ P13(1720) ------*/
        case (kPdgP13m1720_N0) :       /* N0  */
        case (kPdgP13m1720_NP) :       /* N+  */

            /* ------ F15(1680) ------*/
        case (kPdgF15m1680_N0) :       /* N0  */
        case (kPdgF15m1680_NP) :       /* N+  */

            /* ------ P31(1910) ------*/
        case (kPdgP31m1910_DeltaM ) :  /* Delta-  */
        case (kPdgP31m1910_Delta0 ) :  /* Delta0  */
        case (kPdgP31m1910_DeltaP ) :  /* Delta+  */
        case (kPdgP31m1910_DeltaPP) :  /* Delta++ */

            /* ------ P33(1920) ------*/
        case (kPdgP33m1920_DeltaM ) :  /* Delta-  */
        case (kPdgP33m1920_Delta0 ) :  /* Delta0  */
        case (kPdgP33m1920_DeltaP ) :  /* Delta+  */
        case (kPdgP33m1920_DeltaPP) :  /* Delta++ */

            /* ------ F35(1905) ------*/
        case (kPdgF35m1905_DeltaM ) :  /* Delta-  */
        case (kPdgF35m1905_Delta0 ) :  /* Delta0  */
        case (kPdgF35m1905_DeltaP ) :  /* Delta+  */
        case (kPdgF35m1905_DeltaPP) :  /* Delta++ */

            /* ------ F37(1950) ------*/
        case (kPdgF37m1950_DeltaM ) :  /* Delta-  */
        case (kPdgF37m1950_Delta0 ) :  /* Delta0  */
        case (kPdgF37m1950_DeltaP ) :  /* Delta+  */
        case (kPdgF37m1950_DeltaPP) :  /* Delta++ */

            /* ------ P11(1710) ------*/
        case (kPdgP11m1710_N0) :       /* N0  */
        case (kPdgP11m1710_NP) :       /* N+  */

            /* ------ F17(1970) ------*/
            // are you?

            return true;
    }

    return false;
}
//____________________________________________________________________________
bool genie::utils::res::IsDelta(Resonance_t res)
{
    switch(res) {

        case kP33_1232:  return true;  break;
        case kS11_1535:  return false; break;
        case kD13_1520:  return false; break;
        case kS11_1650:  return false; break;
        case kD13_1700:  return false; break;
        case kD15_1675:  return false; break;
        case kS31_1620:  return true;  break;
        case kD33_1700:  return true;  break;
        case kP11_1440:  return false; break;
        case kP33_1600:  return true;  break;
        case kP13_1720:  return false; break;
        case kF15_1680:  return false; break;
        case kP31_1910:  return true;  break;
        case kP33_1920:  return true;  break;
        case kF35_1905:  return true;  break;
        case kF37_1950:  return true;  break;
        case kP11_1710:  return false; break;
        case kF17_1970:  return false; break;
        default:
                         // should not be here - meaningless to return anything
                         assert(false);
    }
    return false;
}
//____________________________________________________________________________
// The values of resonance mass and width is taken from 
// M. Tanabashi et al. (Particle Data Group) Phys. Rev. D 98, 030001
bool genie::utils::res::IsN(Resonance_t res)
{
    return (! utils::res::IsDelta(res) );
}
//____________________________________________________________________________
double genie::utils::res::Mass(Resonance_t res)
{
    switch(res) {
        case kP33_1232  : return 1.232 * units::GeV ; break;
        case kS11_1535  : return 1.535 * units::GeV ; break;
        case kD13_1520  : return 1.515 * units::GeV ; break;
        case kS11_1650  : return 1.655 * units::GeV ; break;
        case kD13_1700  : return 1.700 * units::GeV ; break;
        case kD15_1675  : return 1.675 * units::GeV ; break;
        case kS31_1620  : return 1.630 * units::GeV ; break;
        case kD33_1700  : return 1.700 * units::GeV ; break;
        case kP11_1440  : return 1.430 * units::GeV ; break;
        case kP33_1600  : return 1.600 * units::GeV ; break;
        case kP13_1720  : return 1.720 * units::GeV ; break;
        case kF15_1680  : return 1.685 * units::GeV ; break;
        case kP31_1910  : return 1.890 * units::GeV ; break;
        case kP33_1920  : return 1.920 * units::GeV ; break;
        case kF35_1905  : return 1.880 * units::GeV ; break;
        case kF37_1950  : return 1.930 * units::GeV ; break;
        case kP11_1710  : return 1.710 * units::GeV ; break;
        case kF17_1970  : return 2.190 * units::GeV ; break;
        default: break;
    }
    return -1;
}
//____________________________________________________________________________
double genie::utils::res::Width(Resonance_t res)
{
    switch(res) {
        case kP33_1232  : return 0.117 * units::GeV ; break;
        case kS11_1535  : return 0.150 * units::GeV ; break;
        case kD13_1520  : return 0.115 * units::GeV ; break;
        case kS11_1650  : return 0.140 * units::GeV ; break;
        case kD13_1700  : return 0.150 * units::GeV ; break;
        case kD15_1675  : return 0.150 * units::GeV ; break;
        case kS31_1620  : return 0.140 * units::GeV ; break;
        case kD33_1700  : return 0.300 * units::GeV ; break;
        case kP11_1440  : return 0.350 * units::GeV ; break;
        case kP33_1600  : return 0.320 * units::GeV ; break;
        case kP13_1720  : return 0.250 * units::GeV ; break;
        case kF15_1680  : return 0.130 * units::GeV ; break;
        case kP31_1910  : return 0.280 * units::GeV ; break;
        case kP33_1920  : return 0.260 * units::GeV ; break;
        case kF35_1905  : return 0.330 * units::GeV ; break;
        case kF37_1950  : return 0.285 * units::GeV ; break;
        case kP11_1710  : return 0.100 * units::GeV ; break;
        case kF17_1970  : return 0.500 * units::GeV ; break;
        default: break;
    }
    return -1;
}
//____________________________________________________________________________
double genie::utils::res::BWNorm(Resonance_t res, double N0ResMaxNWidths, double N2ResMaxNWidths, double GnResMaxNWidths)
{
    static genie::utils::res::CacheBWNorm cbwn;
    if (res==kNoResonance)      return -1;
    if (cbwn.cache[res]!=0)   return cbwn.cache[res];


    // Get baryon resonance parameters
    int    IR  = utils::res::ResonanceIndex    (res);
    int    LR  = utils::res::OrbitalAngularMom (res);
    double MR  = utils::res::Mass              (res);
    double WR  = utils::res::Width             (res);

    // imported part of code from src/contrib/misc/bwnorm.C

    double NW = GnResMaxNWidths;
    if(IR==2) NW = N2ResMaxNWidths;
    if(IR==0) NW = N0ResMaxNWidths;

    double Wmin = 0.001;
    double Wmax = MR + NW*WR;
    int N = 1000* TMath::Nint( (Wmax-Wmin)/WR );
    if(N%2==0) N++;

    double dW = (Wmax-Wmin)/(N-1);

    double sum = 0.5 * (genie::utils::bwfunc::BreitWignerL(Wmin,LR,MR,WR,1.0) + genie::utils::bwfunc::BreitWignerL(Wmax,LR,MR,WR,1.0));

    for(int i=1; i<N-1; i++) {
        double W = Wmin + i*dW;
        sum += ( genie::utils::bwfunc::BreitWignerL(W,LR,MR,WR,1.0) * (i%2+1) );
    }
    sum *= (2.*dW/3.);
    cbwn.cache[res]=sum;
    return cbwn.cache[res];
}
//____________________________________________________________________________
int genie::utils::res::OrbitalAngularMom(Resonance_t res)
{
    switch(res) {
        case kP33_1232:  return 1; break;
        case kS11_1535:  return 0; break;
        case kD13_1520:  return 2; break;
        case kS11_1650:  return 0; break;
        case kD13_1700:  return 2; break;
        case kD15_1675:  return 2; break;
        case kS31_1620:  return 0; break;
        case kD33_1700:  return 2; break;
        case kP11_1440:  return 1; break;
        case kP33_1600:  return 1; break;
        case kP13_1720:  return 1; break;
        case kF15_1680:  return 3; break;
        case kP31_1910:  return 1; break;
        case kP33_1920:  return 1; break;
        case kF35_1905:  return 3; break;
        case kF37_1950:  return 3; break;
        case kP11_1710:  return 1; break;
        case kF17_1970:  return 3; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonRes", pFATAL) 
                             << "Unknown resonance " << res;
                         exit(1);
    }
    return 0;
}
//____________________________________________________________________________
int genie::utils::res::ResonanceIndex(Resonance_t res)
{
    switch(res) {
        case kP33_1232:  return 0; break;
        case kS11_1535:  return 1; break;
        case kD13_1520:  return 1; break;
        case kS11_1650:  return 1; break;
        case kD13_1700:  return 1; break;
        case kD15_1675:  return 1; break;
        case kS31_1620:  return 1; break;
        case kD33_1700:  return 1; break;
        case kP11_1440:  return 2; break;
        case kP33_1600:  return 2; break;
        case kP13_1720:  return 2; break;
        case kF15_1680:  return 2; break;
        case kP31_1910:  return 2; break;
        case kP33_1920:  return 2; break;
        case kF35_1905:  return 2; break;
        case kF37_1950:  return 2; break;
        case kP11_1710:  return 2; break;
        case kF17_1970:  return 2; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonRes", pFATAL) 
                             << "Unknown resonance " << res;
                         exit(1);
    }
    return 0;
}
//____________________________________________________________________________
