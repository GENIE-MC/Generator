//____________________________________________________________________________
/*!

\class   genie::flux::GFLUKAAtmoFlux

\brief   A flux driver for the FLUKA 3-D Atmospheric Neutrino Flux

\ref     Astrop.Phys.19 (2003) p.269; hep-ph/0207035; hep-ph/9907408
         Alfredo.Ferrari     <Alfredo.Ferrari@cern.ch>
         Paola.Sala          <Paola.Sala@cern.ch>
         Giuseppe Battistoni <Giuseppe.Battistoni@mi.infn.it>
         Teresa Montaruli    <Teresa.Montaruli@ba.infn.it>

         To be able to use this flux driver you will need to download the
         flux data from:  http://pcbat1.mi.infn.it/~battist/neutrino.html

         Please note that this class expects to read flux files formatted as 
         described in the above FLUKA flux page.
         Each file contains 3 columns:
	 - neutrino energy (GeV) at bin centre
	 - neutrino cos(zenith angle) at bin centre
         - neutrino flux (#neutrinos /GeV /m^2 /sec /sr)
         The flux is given in 40 bins of cos(zenith angle) from -1.0 to 1.0 
         (bin width = 0.05) and 61 equally log-spaced energy bins (20 bins per 
         decade), with Emin = 0.100 GeV.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created July 3, 2005 

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_FLUKA_ATMO_FLUX_H_
#define _GENIE_FLUKA_ATMO_FLUX_H_

#include "Tools/Flux/GAtmoFlux.h"

namespace genie {
namespace flux  {

// Number of cos(zenith) and energy bins in flux simulation
const unsigned int kGFlk3DNumCosThetaBins       = 40;
const double       kGFlk3DCosThetaMin           = -1.0;
const double       kGFlk3DCosThetaMax           =  1.0;
const unsigned int kGFlk3DNumLogEvBins          = 61;
const unsigned int kGFlk3DNumLogEvBinsPerDecade = 20;
const double       kGFlk3DEvMin                 = 0.100; // GeV

class GFLUKAAtmoFlux: public GAtmoFlux {

public :
  GFLUKAAtmoFlux();
 ~GFLUKAAtmoFlux();

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

#endif // _GFLUKA_ATMO_3D_FLUX_I_H_
