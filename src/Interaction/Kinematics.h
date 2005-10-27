//____________________________________________________________________________
/*!

\class    genie::Kinematics

\brief    Kinematic variables for an event

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 08, 2004

*/
//____________________________________________________________________________

#ifndef _KINEMATICS_H_
#define _KINEMATICS_H_

#include <map>
#include <iostream>

#include <TObject.h>

#include "Interaction/KineVar.h"

using std::map;
using std::ostream;

namespace genie {

class Kinematics : public TObject {

public:

  Kinematics();
  Kinematics(const Kinematics & kv);
  virtual ~Kinematics();

  double x     (void) const;
  double y     (void) const;
  double Q2    (void) const;
  double q2    (void) const;
  double W     (void) const;
  double lnx   (void) const;
  double lny   (void) const;
  double lnQ2  (void) const;
  double lnW   (void) const;
  double logx  (void) const;
  double logy  (void) const;
  double logQ2 (void) const;
  double logW  (void) const;

  void   Setx  (double x );
  void   Sety  (double y );
  void   SetQ2 (double Q2);
  void   Setq2 (double q2);
  void   SetW  (double W );

  bool   KVSet(KineVar_t kv) const;
  double GetKV(KineVar_t kv) const;
  void   SetKV(KineVar_t kv, double value);

  void Initialize (void);
  void Copy       (const Kinematics & kinematics);
  void Print      (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const Kinematics & kinematics);

private:

  map<KineVar_t, double> fKV;

ClassDef(Kinematics,1)
};

}       // genie namespace

#endif  // _KINEMATICS_H_
