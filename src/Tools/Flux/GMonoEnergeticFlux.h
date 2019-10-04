//____________________________________________________________________________
/*!

\class    genie::flux::GMonoEnergeticFlux

\brief    A simple GENIE flux driver for monoenergetic neutrinos along the
          z direction. Can handle a mix of neutrinos with their corresponding
          weight. This trivial flux driver was added in so that single energy 
          neutrinos can be easily used with the event generation driver that 
          can handle a target mix or detailed geometries.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Feb 08, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_MONOENERGETIC_FLUX_H_
#define _G_MONOENERGETIC_FLUX_H_

#include <string>
#include <map>

#include <TLorentzVector.h>

#include "Framework/EventGen/GFluxI.h"

using std::string;
using std::map;

namespace genie {
namespace flux  {

class GMonoEnergeticFlux: public GFluxI {

public :
  GMonoEnergeticFlux();
  GMonoEnergeticFlux(double Ev, int pdg);
  GMonoEnergeticFlux(double Ev, const map<int,double> & numap /* pdg -> weight*/);
 ~GMonoEnergeticFlux();

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
  void   Initialize (double Ev, int pdg);
  void   Initialize (double Ev, const map<int,double> & numap);

private:

  // private methods
  void   CleanUp    (void);

  // private data members
  double           fMaxEv;       ///< maximum energy
  PDGCodeList *    fPdgCList;    ///< list of neutrino pdg-codes
  int              fgPdgC;       ///< running generated nu pdg-code
  TLorentzVector   fgP4;         ///< running generated nu 4-momentum
  TLorentzVector   fgX4;         ///< running generated nu 4-position
  map<int, double> fProb;
  double           fProbMax;
};

} // flux namespace
} // genie namespace

#endif // _G_MONOENERGETIC_FLUX_H_
