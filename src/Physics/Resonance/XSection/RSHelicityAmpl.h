//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmpl

\brief    A class holding the Rein-Sehgal's helicity amplitudes.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          RSHelicityAmplModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RS_HELICITY_AMPL_H_
#define _RS_HELICITY_AMPL_H_

#include <iostream>

#include <TMath.h>

#include "Framework/Interaction/Interaction.h"
#include "Physics/Resonance/XSection/FKR.h"

using std::ostream;

namespace genie {

class RSHelicityAmpl;
ostream & operator<< (ostream & stream, const RSHelicityAmpl & hamp);

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

  //! return helicity amplitude
  double AmpMinus1 (void) const  { return fMinus1; } /* f(-1) */
  double AmpPlus1  (void) const  { return fPlus1;  } /* f(+1) */
  double AmpMinus3 (void) const  { return fMinus3; } /* f(-3) */
  double AmpPlus3  (void) const  { return fPlus3;  } /* f(+3) */
  double Amp0Minus (void) const  { return f0Minus; } /* f(0-) */
  double Amp0Plus  (void) const  { return f0Plus;  } /* f(0+) */

  //! return |helicity amplitude|^2
  double Amp2Minus1 (void) const { return TMath::Power(fMinus1, 2.); } /* |f(-1)|^2 */
  double Amp2Plus1  (void) const { return TMath::Power(fPlus1,  2.); } /* |f(+1)|^2 */
  double Amp2Minus3 (void) const { return TMath::Power(fMinus3, 2.); } /* |f(-3)|^2 */
  double Amp2Plus3  (void) const { return TMath::Power(fPlus3,  2.); } /* |f(+3)|^2 */
  double Amp20Minus (void) const { return TMath::Power(f0Minus, 2.); } /* |f(0-)|^2 */
  double Amp20Plus  (void) const { return TMath::Power(f0Plus,  2.); } /* |f(0+)|^2 */

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


