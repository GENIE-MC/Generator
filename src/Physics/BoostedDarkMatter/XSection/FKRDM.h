//____________________________________________________________________________
/*
Dark Matter model parameters calculated in the same form as
the Feynmann-Kislinger-Ravndall (FKR) baryon excitation model parameters. (header)

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#ifndef _FKRDM_H_
#define _FKRDM_H_

#include <iostream>

using std::ostream;

namespace genie {

class FKRDM;
ostream & operator<< (ostream & stream, const FKRDM & parameters);

class FKRDM {

public:

  friend ostream & operator<< (ostream & stream, const FKRDM & parameters);

  //Unchanged variables
  double Lamda;
  double Tv;
  double Rv;
  double S;
  double Ta;
  double Ra;

  //New variables for DM case that differ from the FKR calculation
  double Bs;
  double Cs;
  double Bz;
  double Cz;

  void Reset (void);
  void Print (ostream & stream) const;

  FKRDM();
  ~FKRDM();
};

}        // genie namespace

#endif   // _FKRDM_H_
