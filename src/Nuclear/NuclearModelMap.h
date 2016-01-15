//____________________________________________________________________________
/*!

\class    genie::NuclearModelMap

\brief    This class is a hook for  nuclear models and allows associating each
          one of them with specific nuclei.
          Is a concrete implementation of the NuclearModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 07, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NUCLEAR_MODEL_MAP_H_
#define _NUCLEAR_MODEL_MAP_H_

#include "Nuclear/NuclearModelI.h"

#include <map>

using std::map;

namespace genie {

class NuclearModelMap : public NuclearModelI {

public:
  NuclearModelMap();
  NuclearModelMap(string config);
  virtual ~NuclearModelMap();

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const;
  double         Prob            (double p, double w, const Target & t) const;
  NuclearModel_t ModelType       (const Target & t) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

private:
  void LoadConfig(void);
  const NuclearModelI * SelectModel(const Target & t) const;

  const NuclearModelI * fDefGlobModel;            ///< default basic model (should work for all nuclei)
  map<int, const NuclearModelI *> fRefinedModels; ///< refinements for specific elements
};

}      // genie namespace
#endif // _NUCLEAR_MODEL_MAP_H_
