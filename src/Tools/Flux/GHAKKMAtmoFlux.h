//____________________________________________________________________________
/*!

\class   genie::flux::GHAKKMAtmoFlux

\brief   A driver for the HAKKM 3-D atmospheric neutrino flux (commonly known 
         as the `Honda flux') 

\ref     M. Honda, M. Sajjad Athar, T. Kajita, K. Kasahara, and S. Midorikawa
         Phys. Rev. D 92 (2015) 023004 

         The flux files necessary for running this flux driver can be obtained
         from:​http://www.icrr.u-tokyo.ac.jp/~mhonda/

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created July 9, 2015

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         for the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_HAKKM_ATMO_FLUX_H_
#define _GENIE_HAKKM_ATMO_FLUX_H_

#include "Tools/Flux/GAtmoFlux.h"

namespace genie {
namespace flux  {

// Number of cos(zenith), azimuthat, and log(energy) bins in flux simulation

const unsigned int kGHnd3DNumCosThetaBins       =  20;
const double       kGHnd3DCosThetaMin           =  -1.0;
const double       kGHnd3DCosThetaMax           =   1.0;
const unsigned int kGHnd3DNumPhiBins       	=  12;
const double       kGHnd3DPhiMin           	=   0.0;
const double       kGHnd3DPhiMax           	= 360.0;
const unsigned int kGHnd3DNumLogEvBins          = 101;	
const unsigned int kGHnd3DNumLogEvBinsPerDecade =  20;
const double       kGHnd3DEvMin                 =   0.1; // GeV

class GHAKKMAtmoFlux: public GAtmoFlux {

public :
  GHAKKMAtmoFlux();
 ~GHAKKMAtmoFlux();

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

#endif // _GENIE_ATMNC_ATMO_3D_FLUX_I_H_
