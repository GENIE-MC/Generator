//____________________________________________________________________________
/*!

\class    genie::flux::GPowerLawFlux

\brief    A simple GENIE flux driver for neutrinos following a power law
          spectrum. Can handle a mix of neutrinos with their corresponding
          weight.

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC

\created  May 02, 2023

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _G_POWERLAW_FLUX_H_
#define _G_POWERLAW_FLUX_H_

#include <string>
#include <map>

#include <TLorentzVector.h>

#include "Framework/EventGen/GFluxI.h"

using std::string;
using std::map;

namespace genie {
namespace flux  {

class GPowerLawFlux: public GFluxI {

public :
  GPowerLawFlux();
  GPowerLawFlux(double alpha, double emin, double emax, int pdg);
  GPowerLawFlux(double alpha, double emin, double emax, const map<int,double> & numap /* pdg -> weight*/);
 ~GPowerLawFlux();

  // methods implementing the GENIE GFluxI interface
  const PDGCodeList &    FluxParticles (void) { return *fPdgCList; }
  double                 MaxEnergy     (void) { return  fMaxEv;    }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return  fgPdgC;    }
  double                 Weight        (void) { return  1.0;       }
  const TLorentzVector & Momentum      (void) { return  fgP4;      }
  const TLorentzVector & Position      (void) { return  fgX4;      }
  bool                   End           (void) { return  false;     }
  long int               Index         (void) { return -1;         }
  void                   Clear            (Option_t * opt);
  void                   GenerateWeighted (bool gen_weighted);

  // special setters for this class
  void                   SetDirectionCos (double dx, double dy, double dz);
  void                   SetRayOrigin    (double x,  double y,  double z);
  // setters consistent w/ GCylindTH1Flux naming
  void                   SetNuDirection  (const TVector3 & direction);
  void                   SetBeamSpot     (const TVector3 & spot);

  // allow re-initialization, and/or initialization after default ctor
  void   Initialize (double alpha, double emin, double emax, int pdg);
  void   Initialize (double alpha, double emin, double emax, const map<int,double> & numap);

private:

  // private methods
  void   CleanUp    (void);

  // private data members
  double           fSpectralIndex; ///< spectral index (E^{-alpha})
  double           fMinEv;         ///< minimum energy
  double           fMaxEv;         ///< maximum energy
  PDGCodeList *    fPdgCList;      ///< list of neutrino pdg-codes
  int              fgPdgC;         ///< running generated nu pdg-code
  TLorentzVector   fgP4;           ///< running generated nu 4-momentum
  TLorentzVector   fgX4;           ///< running generated nu 4-position
  map<int, double> fProb;          ///< cumulative probability of neutrino types
  double           fProbMax;
};

} // flux namespace
} // genie namespace

#endif // _G_POWERLAW_FLUX_H_
