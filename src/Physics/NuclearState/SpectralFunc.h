//____________________________________________________________________________
/*!

\class    genie::SpectralFunc

\brief    A realistic spectral function - based nuclear model.
          Is a concrete implementation of the NuclearModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 07, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_H_
#define _SPECTRAL_FUNCTION_H_

#include "Physics/NuclearState/NuclearModelI.h"

class TNtupleD;
class TGraph2D;

namespace genie {

class SpectralFunc : public NuclearModelI {

public:
  SpectralFunc();
  SpectralFunc(string config);
  virtual ~SpectralFunc();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const;
  double         Prob            (double p, double w, const Target & t) const;
  NuclearModel_t ModelType       (const Target &) const 
  {
    return kNucmSpectralFunc;
  }

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

private:
  void       LoadConfig             (void);
  TGraph2D * Convert2Graph          (TNtupleD & data) const;
  TGraph2D * SelectSpectralFunction (const Target & target) const; 

  TGraph2D * fSfFe56;   ///< Benhar's Fe56 SF
  TGraph2D * fSfC12;    ///< Benhar's C12 SF
};

}      // genie namespace
#endif // _SPECTRAL_FUNCTION_H_
