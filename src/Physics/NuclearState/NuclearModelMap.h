//____________________________________________________________________________
/*!

\class    genie::NuclearModelMap

\brief    This class is a hook for  nuclear models and allows associating each
          one of them with specific nuclei.
          Is a concrete implementation of the NuclearModelI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 07, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _NUCLEAR_MODEL_MAP_H_
#define _NUCLEAR_MODEL_MAP_H_

#include <map>
#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class NuclearModelMap : public NuclearModelI {

public:
  NuclearModelMap();
  NuclearModelMap(string config);
  virtual ~NuclearModelMap();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- Allow GenerateNucleon to be called with a radius
  virtual bool   GenerateNucleon (const Target & t,
                                  double hitNucleonRadius) const;
  virtual double  Prob           (double p, double w, const Target & t,
                                  double hitNucleonRadius) const;

  //-- implement the NuclearModelI interface
  bool GenerateNucleon (const Target & t) const {
    return GenerateNucleon(t,0.0);
  }
  double Prob (double p, double w, const Target & t) const {
    return Prob(p,w,t,0.0);
  }
  NuclearModel_t ModelType       (const Target & t) const;

  virtual double FermiMomentum( const Target & t, int nucleon_pdg ) const ;
  virtual double LocalFermiMomentum( const Target & t, int nucleon_pdg, double radius ) const ;


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
