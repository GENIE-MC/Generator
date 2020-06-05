//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
         Igor Kakorin <kakorin@jinr.ru> (latest updates)
         Joint Institute for Nuclear Research 

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
 @ May 9, 2016 -IK
   BWNorm is calculated on-fly.
 @ Nov 12, 2019 -IK
   Updated resonance masses and widths according to PDG-2018.
   Added previously missing resonances P33(1600) and F17(1970).
   Now mass and widths are taken from PDG table via TDatabasePDG and cached.
   Add extra function needed for MK model.
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
#include "Framework/Conventions/Constants.h"

using namespace genie;
using namespace genie::constants;

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
            
        case (kPdgP33m1600_DeltaM)  :  /* Delta-  */
        case (kPdgP33m1600_Delta0)  :  /* Delta0  */
        case (kPdgP33m1600_DeltaP)  :  /* Delta+  */
        case (kPdgP33m1600_DeltaPP) :  /* Delta++ */
            return kP33_1600; break;

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
       
        case (kPdgF17m1970_N0) :       /* N0  */     
        case (kPdgF17m1970_NP) :       /* N+  */ 
            return kF17_1970; break;
    }

    return kNoResonance;
}
//____________________________________________________________________________
int genie::utils::res::PdgCode(Resonance_t res, int Q)
{
 
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
            if(Q == -1) return  kPdgP33m1600_DeltaM;  /* Delta-  */
            if(Q ==  0) return  kPdgP33m1600_Delta0;  /* Delta0  */
            if(Q ==  1) return  kPdgP33m1600_DeltaP;  /* Delta+  */
            if(Q ==  2) return  kPdgP33m1600_DeltaPP; /* Delta++ */

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
            if(Q ==  0) return  kPdgF17m1970_N0; /* N0  */
            if(Q ==  1) return  kPdgF17m1970_NP; /* N+  */ 
            break;

        default:
            return 0;
    }

    return 0;
}
//____________________________________________________________________________
bool genie::utils::res::IsBaryonResonance(int pdgc)
{

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
        case (kPdgP33m1600_DeltaM ) :  /* Delta-  */
        case (kPdgP33m1600_Delta0 ) :  /* Delta0  */
        case (kPdgP33m1600_DeltaP ) :  /* Delta+  */
        case (kPdgP33m1600_DeltaPP) :  /* Delta++ */

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
        case (kPdgF17m1970_N0) :       /* N0  */
        case (kPdgF17m1970_NP) :       /* N+  */

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
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
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
  // Hardcoded data are removed, now they are taken from PDG table via TDatabasePDG and cached
  static double cm[18];
  if (res == kNoResonance)      
    return -1;
  double mass = cm[res];
  if ( mass != 0)
    return mass;
  
  PDGLibrary * pdglib = PDGLibrary::Instance();
  int pdg = genie::utils::res::PdgCode(res, 0); // the mass doesn't depend on resonance charge
  TParticlePDG * res_pdg = pdglib->Find( pdg );
  if (res_pdg != 0)
  {
    mass = res_pdg->Mass() * units::GeV;
    cm[res] = mass;
    return mass;
  }
  
  // should not be here - meaningless to return anything
  gAbortingInErr = true;
  LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
  exit(1);
}
//____________________________________________________________________________
double genie::utils::res::Width(Resonance_t res)
{
  // Hardcoded data are removed, now they are taken from PDG table via TDatabasePDG and cached
  static double cw[18];
  if (res == kNoResonance)      
    return -1;
  double width = cw[res];
  if ( width != 0)
    return width;

  PDGLibrary * pdglib = PDGLibrary::Instance();
  int pdg = genie::utils::res::PdgCode(res, 0); // the width doesn't depend on resonance charge
  TParticlePDG * res_pdg = pdglib->Find( pdg );
  if (res_pdg != 0)
  {
    width = res_pdg->Width() * units::GeV;
    cw[res] = width;
    return width;
  }
  
  // should not be here - meaningless to return anything
  gAbortingInErr = true;
  LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
  exit(1);

}
//____________________________________________________________________________
double genie::utils::res::BWNorm(Resonance_t res, double N0ResMaxNWidths, double N2ResMaxNWidths, double GnResMaxNWidths)
{
    static double cbwn[18];   // cache resonance BWNorm
    if (res == kNoResonance)      
      return -1;
    
    double norm = cbwn[res];
    if ( norm != 0)   
      return norm;

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

    norm = 0.5 * (genie::utils::bwfunc::BreitWignerL(Wmin,LR,MR,WR,1.0) + genie::utils::bwfunc::BreitWignerL(Wmax,LR,MR,WR,1.0));

    for(int i=1; i<N-1; i++) {
        double W = Wmin + i*dW;
        norm += ( genie::utils::bwfunc::BreitWignerL(W,LR,MR,WR,1.0) * (i%2+1) );
    }
    norm *= (2.*dW/3.);
    cbwn[res]=norm;
    return norm;
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
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
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
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
    return 0;
}
//____________________________________________________________________________
int genie::utils::res::Isospin(Resonance_t res)
{
    switch(res) {
        case kP33_1232:  return 3; break;
        case kS11_1535:  return 1; break;
        case kD13_1520:  return 1; break;
        case kS11_1650:  return 1; break;
        case kD13_1700:  return 1; break;
        case kD15_1675:  return 1; break;
        case kS31_1620:  return 3; break;
        case kD33_1700:  return 3; break;
        case kP11_1440:  return 1; break;
        case kP33_1600:  return 3; break;
        case kP13_1720:  return 1; break;
        case kF15_1680:  return 1; break;
        case kP31_1910:  return 3; break;
        case kP33_1920:  return 3; break;
        case kF35_1905:  return 3; break;
        case kF37_1950:  return 3; break;
        case kP11_1710:  return 1; break;
        case kF17_1970:  return 1; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
    return 0;
}
//____________________________________________________________________________
// The function returns 2*j, j-resonance angular momentum 
int genie::utils::res::AngularMom(Resonance_t res)
{
    switch(res) {
        case kP33_1232:  return 3; break;
        case kS11_1535:  return 1; break;
        case kD13_1520:  return 3; break;
        case kS11_1650:  return 1; break;
        case kD13_1700:  return 3; break;
        case kD15_1675:  return 5; break;
        case kS31_1620:  return 1; break;
        case kD33_1700:  return 3; break;
        case kP11_1440:  return 1; break;
        case kP33_1600:  return 3; break;
        case kP13_1720:  return 3; break;
        case kF15_1680:  return 5; break;
        case kP31_1910:  return 1; break;
        case kP33_1920:  return 3; break;
        case kF35_1905:  return 5; break;
        case kF37_1950:  return 7; break;
        case kP11_1710:  return 1; break;
        case kF17_1970:  return 7; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
    return 0;
}
//____________________________________________________________________________
int genie::utils::res::Cjsgn_plus(Resonance_t res)
// signs of angular momentum Clebsch-Gordon coefficient for RSPP model
{
   
    switch(res) {
        case kP33_1232:  return  1; break;
        case kS11_1535:  return  1; break;
        case kD13_1520:  return -1; break;
        case kS11_1650:  return  1; break;
        case kD13_1700:  return -1; break;
        case kD15_1675:  return  1; break;
        case kS31_1620:  return  1; break;
        case kD33_1700:  return -1; break;
        case kP11_1440:  return -1; break;
        case kP33_1600:  return  1; break;
        case kP13_1720:  return  1; break;
        case kF15_1680:  return -1; break;
        case kP31_1910:  return -1; break;
        case kP33_1920:  return  1; break;
        case kF35_1905:  return -1; break;
        case kF37_1950:  return  1; break;
        case kP11_1710:  return -1; break;
        case kF17_1970:  return  1; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
    return 0;
}
//____________________________________________________________________________
// Rein-Sehgal signs for RSPP model
int genie::utils::res::Dsgn(Resonance_t res)
{
   
    switch(res) {
        case kP33_1232:  return  1; break;
        case kS11_1535:  return  1; break;
        case kD13_1520:  return -1; break;
        case kS11_1650:  return -1; break;
        case kD13_1700:  return -1; break;
        case kD15_1675:  return  1; break;
        case kS31_1620:  return  1; break;
        case kD33_1700:  return  1; break;
        case kP11_1440:  return  1; break;
        case kP33_1600:  return -1; break;
        case kP13_1720:  return  1; break;
        case kF15_1680:  return  1; break;
        case kP31_1910:  return -1; break;
        case kP33_1920:  return  1; break;
        case kF35_1905:  return  1; break;
        case kF37_1950:  return  1; break;
        case kP11_1710:  return -1; break;
        case kF17_1970:  return -1; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
    return 0;
}
/*
//  Not used in the latest version
//____________________________________________________________________________
//phases between resonances and nonresonant helicity amplitudes for MK model
double genie::utils::res::AxialPhase(Resonance_t res)
{
   
    switch(res) {
        case kP33_1232:  return   2.97; break;
        case kS11_1535:  return   0.93; break;
        case kD13_1520:  return   0.00; break;
        case kS11_1650:  return   0.00; break;
        case kD13_1700:  return   0.00; break;
        case kD15_1675:  return   0.00; break;
        case kS31_1620:  return   0.00; break;
        case kD33_1700:  return   0.00; break;
        case kP11_1440:  return   0.64; break;
        case kP33_1600:  return   0.00; break;
        case kP13_1720:  return   0.00; break;
        case kF15_1680:  return   0.00; break;
        case kP31_1910:  return   0.00; break;
        case kP33_1920:  return   0.00; break;
        case kF35_1905:  return   0.00; break;
        case kF37_1950:  return   0.00; break;
        case kP11_1710:  return   0.00; break;
        case kF17_1970:  return   0.00; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
}
//____________________________________________________________________________
//phases between resonances and nonresonant helicity amplitudes for MK model
double genie::utils::res::VectorPhase(Resonance_t res)
{
   
    switch(res) {
        case kP33_1232:  return  -2.860; break;
        case kS11_1535:  return   0.000; break;
        case kD13_1520:  return   2.498; break;
        case kS11_1650:  return   kPi;   break;
        case kD13_1700:  return   0.000; break;
        case kD15_1675:  return   0.000; break;
        case kS31_1620:  return   0.000; break;
        case kD33_1700:  return   0.000; break;
        case kP11_1440:  return   0.000; break;
        case kP33_1600:  return   0.000; break;
        case kP13_1720:  return   0.000; break;
        case kF15_1680:  return   0.000; break;
        case kP31_1910:  return   0.000; break;
        case kP33_1920:  return   0.000; break;
        case kF35_1905:  return   0.000; break;
        case kF37_1950:  return   0.000; break;
        case kP11_1710:  return   0.000; break;
        case kF17_1970:  return   0.000; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
}
//____________________________________________________________________________
//Vector form factors parameter for MK model
double genie::utils::res::CV40(Resonance_t res)
{
   
    switch(res) {
        case kP33_1232:  return   1.12; break;
        case kS11_1535:  return   1.00; break;
        case kD13_1520:  return   1.00; break;
        case kS11_1650:  return   1.00; break;
        case kD13_1700:  return   1.00; break;
        case kD15_1675:  return   1.00; break;
        case kS31_1620:  return   1.00; break;
        case kD33_1700:  return   1.00; break;
        case kP11_1440:  return   1.00; break;
        case kP33_1600:  return   1.00; break;
        case kP13_1720:  return   1.00; break;
        case kF15_1680:  return   1.00; break;
        case kP31_1910:  return   1.00; break;
        case kP33_1920:  return   1.00; break;
        case kF35_1905:  return   1.00; break;
        case kF37_1950:  return   1.00; break;
        case kP11_1710:  return   1.00; break;
        case kF17_1970:  return   1.00; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
}
//____________________________________________________________________________
//Axial form factors parameter for MK model
double genie::utils::res::CA50(Resonance_t res)
{
   
    switch(res) {
        case kP33_1232:  return   1.02; break;
        case kS11_1535:  return   1.53; break;
        case kD13_1520:  return   0.50; break;
        case kS11_1650:  return   1.00; break;
        case kD13_1700:  return   1.00; break;
        case kD15_1675:  return   0.60; break;
        case kS31_1620:  return   1.00; break;
        case kD33_1700:  return   1.00; break;
        case kP11_1440:  return   1.00; break;
        case kP33_1600:  return   1.00; break;
        case kP13_1720:  return   1.10; break;
        case kF15_1680:  return   1.30; break;
        case kP31_1910:  return   0.60; break;
        case kP33_1920:  return   0.60; break;
        case kF35_1905:  return   0.60; break;
        case kF37_1950:  return   1.40; break;
        case kP11_1710:  return   1.20; break;
        case kF17_1970:  return   1.00; break;
        default:
                         // should not be here - meaningless to return anything
                         gAbortingInErr = true;
                         LOG("BaryonResUtils", pFATAL) << "Unknown resonance " << res;
                         exit(1);
    }
}
//____________________________________________________________________________
*/
