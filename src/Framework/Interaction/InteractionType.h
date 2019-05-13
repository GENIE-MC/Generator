//____________________________________________________________________________
/*!

\class    genie::InteractionType

\brief    Enumeration of interaction types: e/m, weak cc, weak nc

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Changes required to implement the GENIE Boosted Dark Matter module
          were installed by Josh Berger (Univ. of Wisconsin)

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/ 
//____________________________________________________________________________

#ifndef _INTERACTION_TYPE_H_
#define _INTERACTION_TYPE_H_

#include <cassert>
#include <string>
#include <cstring>

using std::string;

namespace genie {


typedef enum EInteractionType {

  kIntNull   = 0,
  kIntEM,          //
  kIntWeakCC,      //
  kIntWeakNC,      //
  kIntWeakMix,     // CC + NC + interference
  kIntDarkMatter,  // 
  kIntNDecay,      //
  kIntNOsc         //

} InteractionType_t;


class InteractionType
{
public:

  //__________________________________________________________________________
  static string AsString(InteractionType_t type)
  {
    switch (type) {

      case(kIntEM)         : return "EM";                        break;
      case(kIntWeakCC)     : return "Weak[CC]";                  break;
      case(kIntWeakNC)     : return "Weak[NC]";                  break; 
      case(kIntWeakMix)    : return "Weak[CC+NC+interference]";  break;
      case(kIntDarkMatter) : return "DarkMatter";                break; 
      case(kIntNDecay)     : return "NucleonDecay";              break;
      case(kIntNOsc)       : return "NeutronOsc";                break;
      default :              return "Unknown";                   break;
    }
    return "Unknown";    
  }
  //__________________________________________________________________________
  static InteractionType_t FromString(string type) 
  {
    //-- Make uppercase/lowercase irrelevant

    for(unsigned int i=0; i<type.size(); i++) type[i] = toupper(type[i]);
    
    //-- Figure out the ScatteringType_t from the input string

    const char * t = type.c_str();

    if ( strcmp(t,"EM") == 0 ||
            strcmp(t,"E-M") == 0 ||
                strcmp(t,"E/M") == 0 ||
                   strcmp(t,"ELECTROMAGNETIC") == 0 ||
                      strcmp(t,"ELECTRO-MAGNETIC") == 0 ) return kIntEM;

    else if ( strcmp(t,"WEAK-CC") == 0 ||
                  strcmp(t,"CHARGED-CURRENT") == 0 ||
                      strcmp(t,"CHARGED CURRENT") == 0 ||
                          strcmp(t,"WEAK-CHARGED-CURRENT") == 0 ||
                               strcmp(t,"WEAK CHARGED CURRENT") == 0 ||
                                    strcmp(t,"CC") == 0 ) return kIntWeakCC;

    else if ( strcmp(t,"WEAK-NC") == 0 ||
                 strcmp(t,"NEUTRAL-CURRENT") == 0 ||
                      strcmp(t,"NEUTRAL CURRENT") == 0 ||
                          strcmp(t,"WEAK-NEUTRAL-CURRENT") == 0 ||
                               strcmp(t,"WEAK NEUTRAL CURRENT") == 0 ||
                                     strcmp(t,"NC") == 0 ) return kIntWeakNC;
                                     
    else if ( strcmp(t,"NDECAY") == 0 ) return kIntNDecay;

    else if ( strcmp(t,"NOSC") == 0 ) return kIntNOsc;

    else return kIntNull;
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _INTERACTION_TYPE_H_
