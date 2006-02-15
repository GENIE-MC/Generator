//____________________________________________________________________________
/*!

\class   genie::KineGeneratorWithCache

\brief   Abstract class. Provides a data caching mechanism for for concrete
         implementations of the EventRecordVisitorI interface, generating
         kinematics and wishing to cache maximum differential xsecs.

         This class provides some common implementation for handling
         (retrieving, creating, searching, adding to) the cache.
         The various super-classes should implement the ComputeMaxXSec(...)
         method for computing the maximum xsec in case it has not already
         being pushed into the cache at a previous iteration.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 15, 2004

*/
//____________________________________________________________________________

#ifndef _KINE_GENERATOR_WITH_CACHE_H_
#define _KINE_GENERATOR_WITH_CACHE_H_

#include <string>

#include "Base/XSecAlgorithmI.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "Utils/Range1.h"

class TNtuple;

using std::string;

namespace genie {

class KineGeneratorWithCache : public EventRecordVisitorI {

protected:

  KineGeneratorWithCache();
  KineGeneratorWithCache(string name);
  KineGeneratorWithCache(string name, string config);
  ~KineGeneratorWithCache();

  virtual double ComputeMaxXSec (const Interaction * in) const = 0;

  virtual double   MaxXSec           (GHepRecord * evrec) const;
  virtual double   FindMaxXSec       (const Interaction * in) const;
  virtual void     CacheMaxXSec      (const Interaction * in, double xsec) const;
  virtual double   Energy            (const Interaction * in) const;
  virtual TNtuple* AccessCacheBranch (const Interaction * in) const;
};

}      // genie namespace

#endif // _KINE_GENERATOR_WITH_CACHE_H_
