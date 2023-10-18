//____________________________________________________________________________
/*!

\class    genie::SpectralFunc

\brief    A realistic spectral function - based nuclear model.
          Is a concrete implementation of the NuclearModelI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  May 07, 2004

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org


*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_H_
#define _SPECTRAL_FUNCTION_H_

#include <map>

#include "Physics/NuclearState/NuclearModelI.h"

class TGraph2D;
class TH2D;
class TNtupleD;

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

  TH2D* SelectSpectralFunction (const Target& target) const;

 protected:
  void       LoadConfig             (void);
  TGraph2D*  Convert2Graph          (TNtupleD& data) const;
  TH2D* LoadSFDataFile( const std::string& full_file_name, int targetN) const;

  /// The path to the folder containing the spectral function data files
  std::string fDataPath;
  std::string fDataPathProton;
  std::string fDataPathNeutron;

  /// Map storing cached spectral functions. Keys are pair 
  /// <isospin 1/-1,nuclear PDG codes>,
  /// values are 2D histograms representing the probability distribution
  mutable std::map<std::pair<int,int>, TH2D*> fSpectralFunctionMap;

  /// The number of nucleon momentum bins to use when making spectral function
  /// histograms
  int fNumXBins;

  /// The number of removal energy bins to use when making spectral function
  /// histograms
  int fNumYBins;
};

}      // genie namespace
#endif // _SPECTRAL_FUNCTION_H_
