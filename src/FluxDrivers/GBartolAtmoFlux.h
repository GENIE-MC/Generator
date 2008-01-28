//____________________________________________________________________________
/*!

\class   genie::flux::GBartolAtmoFlux

\brief   A flux driver for the Bartol Atmospheric Neutrino Flux

\ref     published shortly as
	 G. Barr, T.K. Gaisser, P. Lipari, S. Robbins and T. Stanev

         To be able to use this flux driver you will need to download the
         flux data from:
	 http://www-pnp.physics.ox.ac.uk/~barr/fluxfiles/0408i/index.html

\author  Christopher Backhouse <c.backhouse1@physics.ox.ac.uk>
         Oxford University

\created January 26, 2008

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GBARTOL_ATMO_FLUX_I_H_
#define _GBARTOL_ATMO_FLUX_I_H_

#include <string>

#include <TLorentzVector.h>

#include "EVGDrivers/GFluxI.h"

// Get hold of kNNu - kind of a hack
#include "FluxDrivers/GFlukaAtmo3DFlux.h"

class TH2D;

using std::string;

namespace genie {
namespace flux  {

// Number of cos8 and energy bins in flux simulation
const unsigned int kNGBrtCos = 20;
const unsigned int kNGBrtEv  = 30;

// Number of cos8 and energies for which the neutrino flux is given
const double kGBrtCos[kNGBrtCos] = {
  -.95, -.85, -.75, -.65, -.55, -.45, -.35, -.25, -.15, -.05,
  +.05, +.15, +.25, +.35, +.45, +.55, +.65, +.75, +.85, +.95
};

const double kGBrtEv[kNGBrtEv] = {
      11.220,   14.125,   17.783,   22.387,   28.184,   35.481,   44.668,
      56.234,   70.795,   89.125,  112.202,  141.254,  177.828,  223.872,
     281.838,  354.813,  446.684,  562.341,  707.946,  891.251, 1122.018,
    1412.538, 1778.279, 2238.721, 2818.383, 3548.134, 4466.836, 5623.413,
    7079.458, 8912.509
};

class GBartolAtmoFlux: public GFluxI {

public :
  GBartolAtmoFlux();
 ~GBartolAtmoFlux();

  //-- methods specific to this flux object
  bool LoadFluxData       (void);
  void SetRadii           (double Rlongitudinal, double Rtransverse);
  void SetNuMuFluxFile    (string filename) { fFluxFile[0] = filename; }
  void SetNuMuBarFluxFile (string filename) { fFluxFile[1] = filename; }
  void SetNuEFluxFile     (string filename) { fFluxFile[2] = filename; }
  void SetNuEBarFluxFile  (string filename) { fFluxFile[3] = filename; }

  //-- methods implementing the GENIE GFluxI interface
  const PDGCodeList &    FluxParticles (void) { return *fPdgCList; }
  double                 MaxEnergy     (void) { return  fMaxEv;    }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return fgPdgC;     }
  double                 Weight        (void) { return 1.0;        }
  const TLorentzVector & Momentum      (void) { return fgP4;       }
  const TLorentzVector & Position      (void) { return fgX4;       }

private:

  //-- private methods
  void   Initialize        (void);
  void   CleanUp           (void);
  void   ResetSelection    (void);
  TH2D * CreateFluxHisto2D (string name, string title);
  bool   FillFluxHisto2D   (TH2D * h2, string filename);
  void   ZeroFluxHisto2D   (TH2D * h2);
  void   AddAllFluxes      (void);
  int    SelectNeutrino    (double Ev, double costheta);

  //-- private data members
  double         fMaxEv;          ///< maximum energy
  PDGCodeList *  fPdgCList;       ///< list of neutrino pdg-codes
  int            fgPdgC;          ///< running generated nu pdg-code
  TLorentzVector fgP4;            ///< running generated nu 4-momentum
  TLorentzVector fgX4;            ///< running generated nu 4-position
  TH2D *         fFlux2D[kNNu];   ///< flux = f(Ev,cos8), 1/neutrino species
  TH2D *         fFluxSum2D;      ///< combined flux = f(Ev,cos8)
  string         fFluxFile[kNNu]; ///< flux file
  int            fNSkipped;       ///< number of skipped fluxes
  double         fRl;             ///< longitudinal radius
  double         fRt;             ///< transverse radius
};


} // flux namespace
} // genie namespace

#endif // _GBARTOL_ATMO_FLUX_I_H_

