//____________________________________________________________________________
/*!

\class    genie::flux::GAtmoFlux

\brief    A base class for the concrete FLUKA and Bartol atmospheric neutrino 
          flux drivers.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  January 26, 2008

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GATMO_FLUX_H_
#define _GATMO_FLUX_H_

#include <string>

#include <TLorentzVector.h>

#include "EVGDrivers/GFluxI.h"

class TH2D;

using std::string;

namespace genie {
namespace flux  {

// Number of neutrino species considered by concrete flux drivers
const unsigned int kNNu = 4;

class GAtmoFlux: public GFluxI {

public :
  virtual ~GAtmoFlux();

  // methods implementing the GENIE GFluxI interface
  virtual const PDGCodeList &    FluxParticles (void) { return *fPdgCList; }
  virtual double                 MaxEnergy     (void) { return  fMaxEv;    }
  virtual bool                   GenerateNext  (void);
  virtual int                    PdgCode       (void) { return fgPdgC;     }
  virtual double                 Weight        (void) { return fWeight;    }
  virtual const TLorentzVector & Momentum      (void) { return fgP4;       }
  virtual const TLorentzVector & Position      (void) { return fgX4;       }

  // methods specific to the atmospheric flux drivers
  void SetRadii           (double Rlongitudinal, double Rtransverse);
  void SetCosThetaBins    (unsigned int nbins, const double * bins);
  void SetEnergyBins      (unsigned int nbins, const double * bins);
  void SetNuMuFluxFile    (string filename) { fFluxFile[0] = filename; }
  void SetNuMuBarFluxFile (string filename) { fFluxFile[1] = filename; }
  void SetNuEFluxFile     (string filename) { fFluxFile[2] = filename; }
  void SetNuEBarFluxFile  (string filename) { fFluxFile[3] = filename; }
  bool LoadFluxData       (void);

protected:

  //
  GAtmoFlux();

  // protected methods
  void    Initialize        (void);
  void    CleanUp           (void);
  void    ResetSelection    (void);
  TH2D *  CreateFluxHisto2D (string name, string title);
  void    ZeroFluxHisto2D   (TH2D * h2);
  void    AddAllFluxes      (void);
  int     SelectNeutrino    (double Ev, double costheta);

  //
  virtual bool FillFluxHisto2D   (TH2D * h2, string filename) = 0;

  // protected data members
  double         fMaxEv;          ///< (declared) maximum energy
  PDGCodeList *  fPdgCList;       ///< (declared) list of neutrino pdg-codes
  int            fgPdgC;          ///< (current) generated nu pdg-code
  TLorentzVector fgP4;            ///< (current) generated nu 4-momentum
  TLorentzVector fgX4;            ///< (current) generated nu 4-position
  string         fFluxFile[kNNu]; ///< (config) input flux files
  double         fRl;             ///< (config) flux neutrino generation surface: longitudinal radius
  double         fRt;             ///< (config) flux neutrino generation surface: transverse radius
  unsigned int   fkNCosBins;      ///< (config) number of cos(theta) bins in input flux data files
  unsigned int   fkNEvBins;       ///< (config) number of energy bins in input flux data files
  const double * fkCosBins;       ///< (config) cos(theta) bins in input flux data files
  const double * fkEvBins;        ///< (config) energy bins in input flux data files
  TH2D *         fFlux2D[kNNu];   ///< flux = f(Ev,cos8), 1/neutrino species
  TH2D *         fFluxSum2D;      ///< combined flux = f(Ev,cos8)
  int            fNSkipped;       ///< number of skipped fluxes
  double         fWeight;         ///< integral of fFluxSum2D used for flux normalization
};

} // flux namespace
} // genie namespace

#endif // _GATMO_FLUX_H_

