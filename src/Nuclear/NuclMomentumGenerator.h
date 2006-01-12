//____________________________________________________________________________
/*!

\class    genie::NuclMomentumGenerator

\brief    Describes a nucleon momentum probability distribution (constructed
          from the attached NuclearPDistributionModelI) and can act as a
          nucleon momentum generator.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#ifndef _NUCLEON_MOMENTUM_GENERATOR_H_
#define _NUCLEON_MOMENTUM_GENERATOR_H_

#include <map>

class TVector3;
class TH1D;

using std::map;

namespace genie {

class NuclMomentumModelI;
class Target;

class NuclMomentumGenerator {

public:

  static NuclMomentumGenerator * Instance (void);

  void     UseProbDistribution (const NuclMomentumModelI * m, const Target & t);
  double   Probability         (double p)           const;
  double   Probability         (const TVector3 & p) const;
  double   RandomMomentum      (void) const;
  TVector3 RandomMomentum3     (void) const;

private:

  NuclMomentumGenerator();
  NuclMomentumGenerator(const NuclMomentumGenerator & nmd);
  virtual ~NuclMomentumGenerator();

  string BuildProbDistributionKey(const NuclMomentumModelI * m, const Target & t);

  static NuclMomentumGenerator * fInstance;

  map<string, TH1D*> fProbDistributionMap;
  TH1D*              fCurrProbDistribution;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (NuclMomentumGenerator::fInstance !=0) {
            delete NuclMomentumGenerator::fInstance;
            NuclMomentumGenerator::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}         // genie namespace

#endif    // _NUCLEON_MOMENTUM_GENERATOR_H_

