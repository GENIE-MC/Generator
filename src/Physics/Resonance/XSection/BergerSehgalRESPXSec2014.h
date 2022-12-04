//____________________________________________________________________________
/*!

\class    genie::BergerSehgalRESPXSec2014

\brief    Computes the double differential cross section for resonance 
          electro- or neutrino-production according to the Berger Sehgal model.

          The computed cross section is the d^2 xsec/ dQ^2 dW \n
          where \n
            \li \c Q^2 : momentum transfer ^ 2
            \li \c W   : invariant mass of the final state hadronic system

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      Berger, Sehgal Phys. Rev. D76, 113004 (2007)

          Modifications within original format of 
	  D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

          Modifications based on a MiniBooNE tune courtesy of J. Nowak
          and S. Dytman.

\author   Steve Dytman
          University of Pittsburgh

          Jarek Nowak 
          University of Lancaster

	  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Sep 15, 2015

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BERGER_SEHGAL_RES_PXSEC_2014_H_
#define _BERGER_SEHGAL_RES_PXSEC_2014_H_

#include "ReinSehgal/BSKLNBaseRESPXSec2014.h"

namespace genie {

 class BergerSehgalRESPXSec2014: public BSKLNBaseRESPXSec2014
 {
   public:
     BergerSehgalRESPXSec2014();
     BergerSehgalRESPXSec2014(string config);
     virtual ~BergerSehgalRESPXSec2014();
 };

}       // genie namespace

#endif  // _BERGER_SEHGAL_RES_PXSEC_2014_H_
