//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmpl

\brief    A class holding the Rein-Seghal's helicity amplitudes.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          RSHelicityAmplModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _RS_HELICITY_AMPL_H_
#define _RS_HELICITY_AMPL_H_

#include <iostream>

#include "Interaction/Interaction.h"
#include "ReinSeghal/FKR.h"
#include "ReinSeghal/RSHelicityAmplModelI.h"

using std::ostream;

namespace genie {

class RSHelicityAmpl {

friend class RSHelicityAmplModelCC;
friend class RSHelicityAmplModelNCp;
friend class RSHelicityAmplModelNCn;
friend class RSHelicityAmplModelEMp;
friend class RSHelicityAmplModelEMn;

public:

  RSHelicityAmpl();
  RSHelicityAmpl(const RSHelicityAmpl & hamp);
  ~RSHelicityAmpl() { }

  double AmpMinus1 (void) const { return fMinus1; }
  double AmpPlus1  (void) const { return fPlus1;  }
  double AmpMinus3 (void) const { return fMinus3; }
  double AmpPlus3  (void) const { return fPlus3;  }
  double Amp0Minus (void) const { return f0Minus; }
  double Amp0Plus  (void) const { return f0Plus;  }

  friend ostream & operator<< (ostream & stream, const RSHelicityAmpl & hamp);

  void Print(ostream & stream) const;

private:

  void   Init(void);

  double fMinus1;
  double fPlus1;
  double fMinus3;
  double fPlus3;
  double f0Minus;
  double f0Plus;
};

}        // genie namespace

#endif   // _RS_HELICITY_AMPL_H_


