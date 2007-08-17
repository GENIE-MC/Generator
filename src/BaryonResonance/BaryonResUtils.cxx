//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - November 25, 2004

 For the namespace documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Interaction/Interaction.h"
#include "PDG/PDGLibrary.h"

using namespace genie;

//____________________________________________________________________________
char * genie::utils::res::AsString(Resonance_t res)
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

    case (-2214) : /* Delta-  */
    case  (2114) : /* Delta0  */
    case  (2214) : /* Delta+  */
    case  (2224) : /* Delta++ */
                   return kP33_1232; break;

    case (22112) : /* N0  */
    case (22212) : /* N+  */
                   return kS11_1535; break;

    case  (1214) : /* N0  */
    case  (2124) : /* N+  */
                   return kD13_1520; break;

    case (32112) : /* N0  */
    case (32212) : /* N+  */
                   return kS11_1650; break;

    case (21214) : /* N0  */
    case (22124) : /* N+  */
                   return kD13_1700; break;

    case  (2116) : /* N0  */
    case  (2216) : /* N+  */
                   return kD15_1675; break;

    case  (1112) : /* Delta-  */
    case  (1212) : /* Delta0  */
    case  (2122) : /* Delta+  */
    case  (2222) : /* Delta++ */
                   return kS31_1620; break;

    case (11114) : /* Delta-  */
    case (12114) : /* Delta0  */
    case (12214) : /* Delta+  */
    case (12224) : /* Delta++ */
                   return kD33_1700; break;

    case (12112) : /* N0  */
    case (12212) : /* N+  */
                   return kP11_1440; break;

    case (31214) : /* N0  */
    case (32124) : /* N+  */
                   return kP13_1720; break;

    case (12116) : /* N0  */
    case (12216) : /* N+  */
                   return kF15_1680; break;

    case (21112) : /* Delta-  */
    case (21212) : /* Delta0  */
    case (22122) : /* Delta+  */
    case (22222) : /* Delta++ */
                   return kP31_1910; break;

    case (21114) : /* Delta-  */
    case (22114) : /* Delta0  */
    case (22214) : /* Delta+  */
    case (22224) : /* Delta++ */
                   return kP33_1920; break;

    case  (1116) : /* Delta-  */
    case  (1216) : /* Delta0  */
    case  (2126) : /* Delta+  */
    case  (2226) : /* Delta++ */
                   return kF35_1905; break;

    case  (1118) : /* Delta-  */
    case  (2118) : /* Delta0  */
    case  (2218) : /* Delta+  */
    case  (2228) : /* Delta++ */
                   return kF37_1950; break;

    case (42112) : /* N0  */
    case (42212) : /* N+  */
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
        if(Q == -1) return -2214; /* Delta-  */
        if(Q ==  0) return  2114; /* Delta0  */
        if(Q ==  1) return  2214; /* Delta+  */
        if(Q ==  2) return  2224; /* Delta++ */
        break;

    case kS11_1535:
        if(Q ==  0) return  22112; /* N0  */
        if(Q ==  1) return  22212; /* N+  */
        break;

    case kD13_1520:
        if(Q ==  0) return  1214; /* N0  */
        if(Q ==  1) return  2124; /* N+  */
        break;

    case kS11_1650:
        if(Q ==  0) return  32112; /* N0  */
        if(Q ==  1) return  32212; /* N+  */
        break;

    case kD13_1700:
        if(Q ==  0) return  21214; /* N0  */
        if(Q ==  1) return  22124; /* N+  */
        break;

    case kD15_1675:
        if(Q ==  0) return  2116; /* N0  */
        if(Q ==  1) return  2216; /* N+  */
        break;

    case kS31_1620:
        if(Q == -1) return 1112; /* Delta-  */
        if(Q ==  0) return 1212; /* Delta0  */
        if(Q ==  1) return 2122; /* Delta+  */
        if(Q ==  2) return 2222; /* Delta++ */
        break;

    case kD33_1700:
        if(Q == -1) return 11114; /* Delta-  */
        if(Q ==  0) return 12114; /* Delta0  */
        if(Q ==  1) return 12214; /* Delta+  */
        if(Q ==  2) return 12224; /* Delta++ */
        break;

    case kP11_1440:
        if(Q ==  0) return  12112; /* N0  */
        if(Q ==  1) return  12212; /* N+  */
        break;

    case kP33_1600:
        return 0;   // ???????

    case kP13_1720:
        if(Q ==  0) return  31214; /* N0  */
        if(Q ==  1) return  32124; /* N+  */
        break;

    case kF15_1680:
        if(Q ==  0) return  12116; /* N0  */
        if(Q ==  1) return  12216; /* N+  */
        break;

    case kP31_1910:
        if(Q == -1) return 21112; /* Delta-  */
        if(Q ==  0) return 21212; /* Delta0  */
        if(Q ==  1) return 22122; /* Delta+  */
        if(Q ==  2) return 22222; /* Delta++ */
        break;

    case kP33_1920:
        if(Q == -1) return 21114; /* Delta-  */
        if(Q ==  0) return 22114; /* Delta0  */
        if(Q ==  1) return 22214; /* Delta+  */
        if(Q ==  2) return 22224; /* Delta++ */
        break;

    case kF35_1905:
        if(Q == -1) return 1116; /* Delta-  */
        if(Q ==  0) return 1216; /* Delta0  */
        if(Q ==  1) return 2126; /* Delta+  */
        if(Q ==  2) return 2226; /* Delta++ */
        break;

    case kF37_1950:
        if(Q == -1) return 1118; /* Delta-  */
        if(Q ==  0) return 2118; /* Delta0  */
        if(Q ==  1) return 2218; /* Delta+  */
        if(Q ==  2) return 2228; /* Delta++ */
        break;

    case kP11_1710:
        if(Q ==  0) return  42112; /* N0  */
        if(Q ==  1) return  42212; /* N+  */
        break;

    case kF17_1970:
        return 0;   // ???????
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

    case (-2214) : /* Delta-  */
    case  (2114) : /* Delta0  */
    case  (2214) : /* Delta+  */
    case  (2224) : /* Delta++ */

    /* ------ S11(1535) ------*/

    case (22112) : /* N0  */
    case (22212) : /* N+  */

    /* ------ D13(1520) ------*/

    case  (1214) : /* N0  */
    case  (2124) : /* N+  */

    /* ------ S11(1650) ------*/

    case (32112) : /* N0  */
    case (32212) : /* N+  */

    /* ------ D13(1700) ------*/

    case (21214) : /* N0  */
    case (22124) : /* N+  */

    /* ------ D15(1675) ------*/

    case  (2116) : /* N0  */
    case  (2216) : /* N+  */

    /* ------ S31(1620) ------*/

    case  (1112) : /* Delta-  */
    case  (1212) : /* Delta0  */
    case  (2122) : /* Delta+  */
    case  (2222) : /* Delta++ */

    /* ------ D33(1700) ------*/

    case (11114) : /* Delta-  */
    case (12114) : /* Delta0  */
    case (12214) : /* Delta+  */
    case (12224) : /* Delta++ */

    /* ------ P11(1440) ------*/

    case (12112) : /* N0  */
    case (12212) : /* N+  */

    /* ------ P33(1600) ------*/


    /* ------ P13(1720) ------*/

    case (31214) : /* N0  */
    case (32124) : /* N+  */

    /* ------ F15(1680) ------*/

    case (12116) : /* N0  */
    case (12216) : /* N+  */

    /* ------ P31(1910) ------*/

    case (21112) : /* Delta-  */
    case (21212) : /* Delta0  */
    case (22122) : /* Delta+  */
    case (22222) : /* Delta++ */

    /* ------ P33(1920) ------*/

    case (21114) : /* Delta-  */
    case (22114) : /* Delta0  */
    case (22214) : /* Delta+  */
    case (22224) : /* Delta++ */

    /* ------ F35(1905) ------*/

    case  (1116) : /* Delta-  */
    case  (1216) : /* Delta0  */
    case  (2126) : /* Delta+  */
    case  (2226) : /* Delta++ */

    /* ------ F37(1950) ------*/

    case  (1118) : /* Delta-  */
    case  (2118) : /* Delta0  */
    case  (2218) : /* Delta+  */
    case  (2228) : /* Delta++ */

    /* ------ P11(1710) ------*/

    case (42112) : /* N0  */
    case (42212) : /* N+  */

    /* ------ F17(1970) ------*/


    /* ------ ?         ------*/
    case (-1114) :
    case (1114)  :

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

  int nuc_pdgc = init_state.Tgt().HitNucPdg();
  int fsl_pdgc = interaction->FSPrimLeptonPdg();

  int q_nuc = int( PDGLibrary::Instance()->Find(nuc_pdgc)->Charge() );
  int q_fsl = int( PDGLibrary::Instance()->Find(fsl_pdgc)->Charge() );
  int q_res = (q_nuc - q_fsl) /3;

  return q_res;
}
//___________________________________________________________________________

