//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - November 25, 2004

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
*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Interaction/Interaction.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"

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
double genie::utils::res::Mass(Resonance_t res)
{
  switch(res) {
    case kP33_1232  : return 1.232 ; break;
    case kS11_1535  : return 1.535 ; break;
    case kD13_1520  : return 1.520 ; break;
    case kS11_1650  : return 1.650 ; break;
    case kD13_1700  : return 1.700 ; break;
    case kD15_1675  : return 1.675 ; break;
    case kS31_1620  : return 1.620 ; break;
    case kD33_1700  : return 1.700 ; break;
    case kP11_1440  : return 1.440 ; break;
    case kP33_1600  : return 1.600 ; break;
    case kP13_1720  : return 1.720 ; break;
    case kF15_1680  : return 1.680 ; break;
    case kP31_1910  : return 1.910 ; break;
    case kP33_1920  : return 1.920 ; break;
    case kF35_1905  : return 1.905 ; break;
    case kF37_1950  : return 1.950 ; break;
    case kP11_1710  : return 1.710 ; break;
    case kF17_1970  : return 1.970 ; break;
    default: break;
  }
  return -1;
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
bool genie::utils::res::IsN(Resonance_t res)
{
  return (! utils::res::IsDelta(res) );
}
//____________________________________________________________________________
int genie::utils::res::ResonanceCharge(const Interaction * interaction) 
{
// Figure out what the resonance charge should be to conserve the charge in
// RES interactions 

  const InitialState & init_state = interaction->InitState();

  int prb_pdgc = init_state.ProbePdg();
  int nuc_pdgc   = init_state.Tgt().HitNucPdg();
  int fsl_pdgc   = interaction->FSPrimLeptonPdg();

  int q_prb = int( PDGLibrary::Instance()->Find(prb_pdgc)->Charge() );
  int q_nuc = int( PDGLibrary::Instance()->Find(nuc_pdgc)->Charge() );
  int q_fsl = int( PDGLibrary::Instance()->Find(fsl_pdgc)->Charge() );
  int q_res = (q_prb + q_nuc - q_fsl) /3;

  return q_res;
}
//___________________________________________________________________________

