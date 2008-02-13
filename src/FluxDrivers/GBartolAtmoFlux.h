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

#include "FluxDrivers/GAtmoFlux.h"

namespace genie {
namespace flux  {

// Number of cos8 and energy bins in flux simulation
const unsigned int kNGBrtCos = 20;
const unsigned int kNGBrtEv  = 30;

const double kGBrtCos[kNGBrtCos+1] = {
  -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1,
   0, +.1, +.2, +.3, +.4, +.5, +.6, +.7, +.8, +.9, +1
};

const double kGBrtEv[kNGBrtEv+1] = {
  10.00, 12.59, 15.85, 19.95, 25.12, 31.62, 39.81, 50.12, 63.10, 79.43,
  100.0, 125.9, 158.5, 199.5, 251.2, 316.2, 398.1, 501.2, 631.0, 794.3,
  1000., 1259., 1585., 1995., 2512., 3162., 3981., 5012., 6310., 7943.,
  10000
};

class GBartolAtmoFlux: public GAtmoFlux {

public :
  GBartolAtmoFlux();
 ~GBartolAtmoFlux();

  // Most implementation is derived from the base GAtmoFlux
  // The concrete driver is only required to implement a function for
  // loading the input data files
private:
  bool FillFluxHisto2D(TH2D * h2, string filename);
};

} // flux namespace
} // genie namespace

#endif // _GBARTOL_ATMO_FLUX_I_H_

