//____________________________________________________________________________
/*!

\class   genie::flux::GBGLRSAtmoFlux

\brief   A flux driver for the Bartol Atmospheric Neutrino Flux

\ref     G. Barr, T.K. Gaisser, P. Lipari, S. Robbins and T. Stanev,
         astro-ph/0403630

         To be able to use this flux driver you will need to download the
         flux data from:  http://www-pnp.physics.ox.ac.uk/~barr/fluxfiles/

         Please note that this class expects to read flux files formatted as 
         described in the above BGLRS flux page.
         Each file contains 5 columns:
         - neutrino energy (GeV) at bin centre
         - neutrino cos(zenith angle) at bin centre
         - neutrino flux dN/dlnE (#neutrinos /m^2 /sec /sr)
         - MC statistical error on the flux (not used here)
         - Number of unweighted events entering in the bin (not used here)
         The flux is given in 20 bins of cos(zenith angle) from -1.0 to 1.0 
         (bin width = 0.1) and 30 equally log-spaced energy bins (10 bins per 
         decade), with Emin = 10.00 GeV.

         Note that in the BGLRS input files the flux is defined as dN/dlnE, 
         while in the FLUKA files the flux is defined as dN/dE.
         We compensate for logarithmic units (dlnE = dE/E) as we read-in the
         BGLRS files.

\author  Christopher Backhouse <c.backhouse1@physics.ox.ac.uk>
         Oxford University

\created January 26, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_BGLRS_ATMO_FLUX_H_
#define _GENIE_BGLRS_ATMO_FLUX_H_

#include "Tools/Flux/GAtmoFlux.h"

namespace genie {
namespace flux  {

// Number of cos(zenith) and energy bins in flux simulation
const unsigned int kBGLRS3DNumCosThetaBins           = 20;
const double       kBGLRS3DCosThetaMin               = -1.0;
const double       kBGLRS3DCosThetaMax               =  1.0;
const unsigned int kBGLRS3DNumLogEvBinsLow           = 40;
const unsigned int kBGLRS3DNumLogEvBinsPerDecadeLow  = 20;
const unsigned int kBGLRS3DNumLogEvBinsHigh          = 30;
const unsigned int kBGLRS3DNumLogEvBinsPerDecadeHigh = 10;
const double       kBGLRS3DEvMin                     = 0.1; // GeV

class GBGLRSAtmoFlux: public GAtmoFlux {

public :
  GBGLRSAtmoFlux();
 ~GBGLRSAtmoFlux();

  //     
  // Most implementation is derived from the base GAtmoFlux
  // The concrete driver is only required to implement a function for
  // loading the input data files
  //

private:

  void SetBinSizes   (void);
  bool FillFluxHisto (int nu_pdg, string filename);
};

} // flux namespace
} // genie namespace

#endif // _GBARTOL_ATMO_FLUX_I_H_

