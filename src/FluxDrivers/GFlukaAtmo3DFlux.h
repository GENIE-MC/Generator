//____________________________________________________________________________
/*!

\class   genie::GFlukaAtmo3DFlux

\brief   A flux driver for the FLUKA 3-D Atmospheric Neutrino Flux

\ref     Astrop.Phys.19 (2003) p.269; hep-ph/0207035; hep-ph/9907408
         Alfredo.Ferrari     <Alfredo.Ferrari@cern.ch>
         Paola.Sala          <Paola.Sala@cern.ch>
         Giuseppe Battistoni <Giuseppe.Battistoni@mi.infn.it>
         Teresa Montaruli    <Teresa.Montaruli@ba.infn.it>

         To be able to use this flux driver you will need to download the
         flux data from:  http://lxmi.mi.infn.it/~battist/neutrino.html

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 3, 2005 [during the most boring MINOS shift ever!]

*/
//____________________________________________________________________________

#ifndef _GFLUKA_ATMO_3D_FLUX_I_H_
#define _GFLUKA_ATMO_3D_FLUX_I_H_

#include <string>
#include <TLorentzVector.h>
#include "GMCJob/GFluxI.h"

using std::string;

class TH2D;

namespace genie {

class PDGCodeList;

namespace flux  {

// Number of cos8 and energy bins in flux simulation
const unsigned int kNGFlk3DCos = 40;
const unsigned int kNGFlk3DEv  = 61;

// Number of neutrino species
const unsigned int kNNu = 4;

// Number of cos8 and energis for which the neutrino flux is given
const double kGFlk3DCos[kNGFlk3DCos] = {
  -0.975, -0.925, -0.875, -0.825, -0.775, -0.725, -0.675, -0.625,
  -0.575, -0.525, -0.475, -0.425, -0.375, -0.325, -0.275, -0.225,
  -0.175, -0.125, -0.075, -0.025,  0.025,  0.075,  0.125,  0.175,
   0.225,  0.275,  0.325,  0.375,  0.425,  0.475,  0.525,  0.575,
   0.625,  0.675,  0.725,  0.775,  0.825,  0.875,  0.925,  0.975
};
const double kGFlk3DEv[kNGFlk3DEv]   = {
    0.106,  0.119,  0.133,  0.150,   0.168,  0.188,  0.211,
    0.237,  0.266,  0.299,  0.335,   0.376,  0.422,  0.473,
    0.531,  0.596,  0.668,  0.750,   0.841,  0.944,  1.059,
    1.189,  1.334,  1.496,  1.679,   1.884,  2.113,  2.371,
    2.661,  2.985,  3.350,  3.758,   4.217,  4.732,  5.309,
    5.957,  6.683,  7.499,  8.414,   9.441, 10.593, 11.885,
   13.335, 14.962, 16.788, 18.837,  21.135, 23.714, 26.607,
   29.854, 33.497, 37.584, 42.170,  47.315, 53.088, 59.566,
   66.834, 74.989, 84.140, 94.406, 105.925
};

class GFlukaAtmo3DFlux: public GFluxI {

public :

  GFlukaAtmo3DFlux();
  ~GFlukaAtmo3DFlux();

  //-- methods specific to this flux object
  bool LoadFluxData       (void);
  void SetRadii           (double Rlongitudinal, double Rtransverse);
  void SetNuMuFluxFile    (string filename) { fFluxFile[0] = filename; }
  void SetNuMuBarFluxFile (string filename) { fFluxFile[1] = filename; }
  void SetNuEFluxFile     (string filename) { fFluxFile[2] = filename; }
  void SetNuEBarFluxFile  (string filename) { fFluxFile[3] = filename; }

  //-- methods implementing the GENIE GFluxI interface
  const PDGCodeList &    FluxParticles (void) = 0;
  double                 MaxEnergy     (void) = 0;
  bool                   GenerateNext  (void) = 0;
  int                    PdgCode       (void) = 0;
  const TLorentzVector & Momentum      (void) = 0;
  const TLorentzVector & Position      (void) = 0;

private:

  //-- private methods for initializing / loading flux
  void Initialize        (void);
  bool FillFluxHisto2D   (TH2D * h2, string filename);
  void CreateFluxHisto2D (TH2D * h2, string name, string title);
  void ZeroFluxHisto2D   (TH2D * h2);
  void AddAllFluxes      (void);
  int  SelectNeutrino    (double Ev, double costheta);

  //-- private data members
  double         fMaxEv;          ///< maximum energy
  PDGCodeList *  fPdgCList;       ///< list of neutrino pdg-codes
  int            fgPdgC;          ///< running generated nu pdg-code
  TLorentzVector fgP4;            ///< running generated nu 4-momentum
  TLorentzVector fgX4;            ///< running generated nu 4-position
  TH2D *         fFlux2D[kNNu];   ///< flux = f(Ev,cos8)
  TH2D *         fFluxSum2D;      ///< combined flux = f(Ev,cos8)
  string         fFluxFile[kNNu]; ///< flux file
  int            fNSkipped;       ///< number of skipped fluxes
  double         fRl;             ///< longitudinal radius
  double         fRt;             ///< transverse radius
};

} // flux namespace
} // genie namespace

#endif // _GFLUKA_ATMO_3D_FLUX_I_H_
