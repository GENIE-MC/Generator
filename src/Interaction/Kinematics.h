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
  ~Kinematics();

  double x       (void) const;
  double y       (void) const;
  double Q2      (void) const;
  double q2      (void) const;
  double W       (void) const;
  double Logx    (void) const;
  double Logy    (void) const;
  double LogQ2   (void) const;
  double LogW    (void) const;
  double Log10x  (void) const;
  double Log10y  (void) const;
  double Log10Q2 (void) const;
  double Log10W  (void) const;

  void   Setx  (double x );
  void   Sety  (double y );
  void   SetQ2 (double Q2);
  void   Setq2 (double q2);
  void   SetW  (double W );

  bool   KVSet(KineVar_t kv) const;
  double GetKV(KineVar_t kv) const;
  void   SetKV(KineVar_t kv, double value);

  //! Copy, reset, compare and print itself
  void Reset    (void);
  void Copy     (const Kinematics & kine);
  void Print    (ostream & stream) const;

  Kinematics &     operator =  (const Kinematics & kine);
  friend ostream & operator << (ostream & stream, const Kinematics & kine);

private:

  //! Private data members
  map<KineVar_t, double> fKV;

ClassDef(Kinematics,1)
};

}       // genie namespace

#endif  // _KINEMATICS_H_
