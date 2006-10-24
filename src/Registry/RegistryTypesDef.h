//____________________________________________________________________________
/*!

\class    genie::RegistryItem

\brief    A templated concrete implementation of the RegistryItemI interface.
          Provides an arbitrary basic type (bool, int, double, string) value
          for RegistryI-type containers.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _REGISTRY_TYPES_H_
#define _REGISTRY_TYPES_H_

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

#endif // _REGISTRY_TYPES_H_
