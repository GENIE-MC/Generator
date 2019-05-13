//____________________________________________________________________________
/*!

  \namespace genie::utils::res

  \brief     Baryon Resonance utilities.

  \author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
  University of Liverpool & STFC Rutherford Appleton Lab

  \created   November 25, 2004

  \cpright   Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE
  */
//____________________________________________________________________________

#ifndef _BARYON_RESONANCE_UTILS_H_
#define _BARYON_RESONANCE_UTILS_H_

#include <string>
#include <map> 

#include "Framework/ParticleData/BaryonResonance.h"

using std::string;

namespace genie {

    namespace utils {
        namespace res {

            const char* AsString   (Resonance_t res);  ///< resonance id -> string
            Resonance_t FromString (const char * res); ///< string -> resonance id

            int         PdgCode     (Resonance_t res, int Q); ///< (resonance id, charge) -> PDG code
            Resonance_t FromPdgCode (int pdgc);               ///< PDG code -> resonance id

            bool        IsBaryonResonance (int pdgc);         		 ///< is input a baryon resonance?
            bool        IsDelta           (Resonance_t res);  		 ///< is it a Delta resonance?
            bool        IsN               (Resonance_t res);  		 ///< is it an N resonance?
            double      Mass              (Resonance_t res); 			 ///< resonance mass (GeV)
            double      Width             (Resonance_t res); 			 ///< resonance width (GeV)
            double      BWNorm            (Resonance_t res, 
                    double N0ResMaxNWidths=6, 
                    double N2ResMaxNWidths=2, 
                    double GnResMaxNWidths=4);  ///< breit-wigner normalization factor
            int         OrbitalAngularMom (Resonance_t res);  		///< orbital angular momentum
            int         ResonanceIndex    (Resonance_t res);  		///< resonance idx, quark model / SU(6)
            struct      CacheBWNorm		{
                CacheBWNorm()
                {
                    cache[kP33_1232] = 0.0;
                    cache[kS11_1535] = 0.0;
                    cache[kD13_1520] = 0.0;
                    cache[kS11_1650] = 0.0;
                    cache[kD13_1700] = 0.0;
                    cache[kD15_1675] = 0.0;
                    cache[kS31_1620] = 0.0;
                    cache[kD33_1700] = 0.0;
                    cache[kP11_1440] = 0.0;
                    cache[kP33_1600] = 0.0;
                    cache[kP13_1720] = 0.0;
                    cache[kF15_1680] = 0.0;
                    cache[kP31_1910] = 0.0;
                    cache[kP33_1920] = 0.0;
                    cache[kF35_1905] = 0.0;
                    cache[kF37_1950] = 0.0;
                    cache[kP11_1710] = 0.0;
                    cache[kF17_1970] = 0.0;
                }
                std::map <genie::EResonance, double> cache;
            };															///< cached breit-wigner normalization factor

        }        // res   namespace
    }        // utils namespace
}        // genie namespace

#endif   // _BARYON_RESONANCE_UTILS_H_
