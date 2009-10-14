//____________________________________________________________________________
/*!

\class   genie::flux::GFlukaAtmo3DFlux

\brief   A flux driver for the FLUKA 3-D Atmospheric Neutrino Flux

\ref     Astrop.Phys.19 (2003) p.269; hep-ph/0207035; hep-ph/9907408
         Alfredo.Ferrari     <Alfredo.Ferrari@cern.ch>
         Paola.Sala          <Paola.Sala@cern.ch>
         Giuseppe Battistoni <Giuseppe.Battistoni@mi.infn.it>
         Teresa Montaruli    <Teresa.Montaruli@ba.infn.it>

         To be able to use this flux driver you will need to download the
         flux data from:  http://lxmi.mi.infn.it/~battist/neutrino.html

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created July 3, 2005 [during the most boring MINOS shift ever!]

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GFLUKA_ATMO_3D_FLUX_I_H_
#define _GFLUKA_ATMO_3D_FLUX_I_H_

#include "FluxDrivers/GAtmoFlux.h"

namespace genie {
namespace flux  {

// Number of cos8 and energy bins in flux simulation
const unsigned int kNGFlk3DCos = 40;
const unsigned int kNGFlk3DEv  = 61;

const double kGFlk3DCos[kNGFlk3DCos+1] = {
  -1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65,
  -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25,
  -0.20, -0.15, -0.10, -0.05,  0.00,  0.05,  0.10,  0.15,
   0.20,  0.25,  0.30,  0.35,  0.40,  0.45,  0.50,  0.55,
   0.60,  0.65,  0.70,  0.75,  0.80,  0.85,  0.90,  0.95, 1.00
};

const double kGFlk3DEv[kNGFlk3DEv+1]   = {
     0.100,   0.112,  0.126,  0.141,  0.158,  0.178,  0.200,  0.224,  0.251,  0.282,
     0.316,   0.355,  0.398,  0.447,  0.501,  0.562,  0.631,  0.708,  0.794,  0.891,
     1.000,   1.120,  1.260,  1.410,  1.580,  1.780,  2.000,  2.240,  2.510,  2.820,
     3.160,   3.550,  3.980,  4.470,  5.010,  5.620,  6.310,  7.080,  7.940,  8.910,
    10.000,  11.200, 12.600, 14.100, 15.800, 17.800, 20.000, 22.400, 25.100, 28.200,
    31.600,  35.500, 39.800, 44.700, 50.100, 56.200, 63.100, 70.800, 79.400, 89.100,
   100.000, 112.000
};

class GFlukaAtmo3DFlux: public GAtmoFlux {

public :
  GFlukaAtmo3DFlux();
 ~GFlukaAtmo3DFlux();

  // Override one of the default interface methods to take into account
  // ad-hoc factor of 100 required for bringing agreement with Bartol
  // flux simulation - to be investigated.
  double Weight (void) { return 100*fWeight; }

  // Most implementation is derived from the base GAtmoFlux
  // The concrete driver is only required to implement a function for
  // loading the input data files
private:
  bool FillFluxHisto2D(TH2D * h2, string filename);
};

} // flux namespace
} // genie namespace

#endif // _GFLUKA_ATMO_3D_FLUX_I_H_
