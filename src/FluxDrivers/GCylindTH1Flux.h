//____________________________________________________________________________
/*!

\class   genie::flux::GCylindTH1Flux

\brief   A simple GENIE flux driver. Generates a 'cylindrical' neutrino beam
         along the input direction, with the input transverse radius and
         centered at the input beam spot position.
         The energies are generated from the input energy spectrum (TH1D)
         Multiple neutrino species can be generated (you will need to supply
         an energy spectrum for each).

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 4, 2005

*/
//____________________________________________________________________________

#ifndef _G_TH1_CYLICDRICAL_FLUX_H_
#define _G_TH1_CYLICDRICAL_FLUX_H_

#include <vector>

#include <TLorentzVector.h>

#include "EVGDrivers/GFluxI.h"

class TH1D;
class TVector3;

using std::vector;

namespace genie {
namespace flux  {

class GCylindTH1Flux: public GFluxI {

public :

  GCylindTH1Flux();
  ~GCylindTH1Flux();

  //-- methods specific to this flux object
  void SetNuDirection      (const TVector3 & direction);
  void SetBeamSpot         (const TVector3 & spot);
  void SetTransverseRadius (double Rt);
  void AddEnergySpectrum   (int nu_pdgc, TH1D * spectrum);

  //-- methods implementing the GENIE GFluxI interface
  const PDGCodeList &    FluxParticles (void);
  double                 MaxEnergy     (void);
  bool                   GenerateNext  (void);
  int                    PdgCode       (void);
  const TLorentzVector & Momentum      (void);
  const TLorentzVector & Position      (void);

private:

  //-- private methods
  void   Initialize        (void);
  void   CleanUp           (void);
  void   ResetSelection    (void);
  void   AddAllFluxes      (void);
  int    SelectNeutrino    (double Ev);

  //-- private data members
  double         fMaxEv;       ///< maximum energy
  PDGCodeList *  fPdgCList;    ///< list of neutrino pdg-codes
  int            fgPdgC;       ///< running generated nu pdg-code
  TLorentzVector fgP4;         ///< running generated nu 4-momentum
  TLorentzVector fgX4;         ///< running generated nu 4-position
  vector<TH1D *> fSpectrum;    ///< flux = f(Ev), 1/neutrino species
  TH1D *         fTotSpectrum; ///< combined flux = f(Ev)
  TVector3 *     fDirVec;      ///< neutrino direction
  TVector3 *     fBeamSpot;    ///< beam spot position
  double         fRt;          ///< transverse size of neutrino beam
};

} // flux namespace
} // genie namespace

#endif // _G_TH1_CYLICDRICAL_FLUX_H_
