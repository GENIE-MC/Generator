//____________________________________________________________________________
/*!

\namespace genie::utils::res

\brief     Baryon Resonance utilities.


  \authors    Costas Andreopoulos <c.andreopoulos \at cern.ch>
              University of Liverpool \n
              Updates were made by 
              Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research 
  
  \created   November 25, 2004
  
  \update   November 12, 2019
            Added extra functions for MK model. \n
            Updated resonance masses and widths according to PDG-2018. \n
            Added previously missing resonances P33(1600) and F17(1970). \n
            Now mass and widths are taken from PDG table via TDatabasePDG and cached. \n


\cpright   Copyright (c) 2003-2025, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
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
					   double      N0ResMaxNWidths=6, 
					   double      N2ResMaxNWidths=2, 
					   double      GnResMaxNWidths=4);  ///< breit-wigner normalization factor
            int         OrbitalAngularMom (Resonance_t res);  		///< orbital angular momentum
            int         ResonanceIndex    (Resonance_t res);  		///< resonance idx, quark model / SU(6)
            int         Isospin           (Resonance_t res);
            int         AngularMom        (Resonance_t res);
            int         Cjsgn_plus        (Resonance_t res);
            int         Dsgn              (Resonance_t res);
/*
            // Not used in the latest version
            double      VectorPhase       (Resonance_t res);
            double      AxialPhase        (Resonance_t res);
            double      CV40              (Resonance_t res);
            double      CA50              (Resonance_t res);
*/ 

        }        // res   namespace
    }        // utils namespace
}        // genie namespace

#endif   // _BARYON_RESONANCE_UTILS_H_
