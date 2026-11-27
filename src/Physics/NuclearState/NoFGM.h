#ifndef _NO_FGM_H_
#define _NO_FGM_H_

#include <map>
#include <TVector3.h>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class NoFGM : public NuclearModelI {

public:
  NoFGM();
  NoFGM(string config);
  virtual ~NoFGM();


  //-- allow methods to be called with a radius
  bool   GenerateNucleon (const Target & t, double) const;
  //-- implement the NuclearModelI interface
  bool GenerateNucleon (const Target & t) const;
  double Prob (double , double , const Target & ) const ;
  double Prob            (double p, double w, const Target & t,  double ) const ;
  NuclearModel_t ModelType       (const Target &) const { return kNucmNoFermiGas;};
  virtual double LocalFermiMomentum( const Target & , int , double  ) const ;
  inline double Momentum (void)  const { return 0; }
  inline void SetMomentum3(const TVector3 & mom) const {fCurrMomentum = mom;}
;
  void Configure (const Registry & /*config*/);
  void Configure (string /*param_set*/) ;

 protected:
  void   LoadConfig (void);


};

}         // genie namespace
#endif    // _LOCAL_FGM_H_

