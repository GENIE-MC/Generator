//____________________________________________________________________________
/*!

\class    genie::RegistryItemTypeId

\brief    An enumeration of Registry item types

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  October 20, 2006

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _REGISTRY_ITEM_TYPE_ID_H_
#define _REGISTRY_ITEM_TYPE_ID_H_

#include<string>

using std::string;

namespace genie {

typedef enum ERgType {

  kRgUndefined = 0,
  kRgBool,
  kRgInt,
  kRgDbl,
  kRgStr,
  kRgAlg,
  kRgH1F,
  kRgH2F,
  kRgTree

} RgType_t;

class RgType {

public:
  static string AsString(RgType_t rt)
  {
    switch (rt) {
      case (kRgUndefined) : return "undefined"; break;
      case (kRgBool)      : return "bool";      break;
      case (kRgInt)       : return "int";       break;
      case (kRgDbl)       : return "double";    break;
      case (kRgStr)       : return "string";    break;
      case (kRgAlg)       : return "alg";       break;
      case (kRgH1F)       : return "h1f";       break;
      case (kRgH2F)       : return "h2f";       break;
      case (kRgTree)      : return "tree";      break;
      default             : return "undefined";
    }
  }
};

}        // genie namespace
#endif   // _REGISTRY_ITEM_TYPE_ID_H_
