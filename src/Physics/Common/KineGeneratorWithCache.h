//____________________________________________________________________________
/*!

\class    genie::KineGeneratorWithCache

\brief    Abstract class. Provides a data caching mechanism for for concrete
          implementations of the EventRecordVisitorI interface, generating
          kinematics and wishing to cache maximum differential xsecs.

          This class provides some common implementation for handling
          (retrieving, creating, searching, adding to) the cache.
          The various super-classes should implement the ComputeMaxXSec(...)
          method for computing the maximum xsec in case it has not already
          being pushed into the cache at a previous iteration.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  December 15, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KINE_GENERATOR_WITH_CACHE_H_
#define _KINE_GENERATOR_WITH_CACHE_H_

#include <string>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Utils/Range1.h"

using std::string;

namespace genie {

class CacheBranchFx;
class XSecAlgorithmI;

class KineGeneratorWithCache : public EventRecordVisitorI {

protected:
  KineGeneratorWithCache();
  KineGeneratorWithCache(string name);
  KineGeneratorWithCache(string name, string config);
  ~KineGeneratorWithCache();

  virtual double ComputeMaxXSec (const Interaction * in) const = 0;
  virtual double MaxXSec        (GHepRecord * evrec) const;
  virtual double FindMaxXSec    (const Interaction * in) const;
  virtual void   CacheMaxXSec   (const Interaction * in, double xsec) const;
  virtual double Energy         (const Interaction * in) const;

  virtual CacheBranchFx * AccessCacheBranch (const Interaction * in) const;

  virtual void AssertXSecLimits (const Interaction * in, double xsec, double xsec_max) const;

  mutable const XSecAlgorithmI * fXSecModel;

  double fSafetyFactor;         ///< maxxsec -> maxxsec * safety_factor
  double fMaxXSecDiffTolerance; ///< max{100*(xsec-maxxsec)/.5*(xsec+maxxsec)} if xsec>maxxsec
  double fEMin;                 ///< min E for which maxxsec is cached - forcing explicit calc.
  bool   fGenerateUniformly;    ///< uniform over allowed phase space + event weight?
};

}      // genie namespace

#endif // _KINE_GENERATOR_WITH_CACHE_H_
