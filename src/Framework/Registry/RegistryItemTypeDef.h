//____________________________________________________________________________
/*!

\class    genie::RegistryItemTypeDef

\brief    Definition of Registry item types

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  October 20, 2006

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _REGISTRY_ITEM_TYPE_DEF_H_
#define _REGISTRY_ITEM_TYPE_DEF_H_

#include <map>
#include <string>
#include <ostream>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

using std::pair;
using std::string;
using std::ostream;

typedef string       RgKey;
typedef int          RgInt;
typedef bool         RgBool;
typedef double       RgDbl;
typedef string       RgStr;
typedef const char * RgCChAr;
typedef TH1F *       RgH1F;
typedef TH2F *       RgH2F;
typedef TTree *      RgTree;

class RgAlg {
public:
  RgAlg();
  RgAlg(string n, string c);
 ~RgAlg();
  friend ostream & operator << (ostream & stream, const RgAlg & alg);
  RgAlg &          operator =  (const RgAlg & alg);
  string  name;
  string  config;
};

#endif // _REGISTRY_ITEM_TYPE_DEF_H_
