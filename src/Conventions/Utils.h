/*!___________________________________________________________________________

\class    Utils

\package  GENIE

\brief    Etc. utils. -- remove -- remove -- remove

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Aug 20, 2004
 
____________________________________________________________________________*/

#ifndef _UTILS_H_
#define _UTILS_H_

#include <TParticlePDG.h>

#include "Conventions/Constants.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Messenger/Messenger.h"

using namespace genie::constants;

namespace genie {

class Utils
{
public:

  //__________________________________________________________________________
  static double QuarkCharge(int pdgc)
  {
    if( pdg::IsQuark(pdgc) ) {

      return PDGLibrary::Instance()->Find(pdgc)->Charge();

    }  else {
        LOG("Utils", pERROR)
               << "Input PDG code " << pdgc << " does not belong to a quark";
    }
    return 0;
  }
  //__________________________________________________________________________
  static double CkmElement(int init_pdgc, int fin_pdgc)
  {
     if( pdg::IsUQuark(init_pdgc) ) {

         if( pdg::IsDQuark(fin_pdgc) ) return kVud;
         if( pdg::IsSQuark(fin_pdgc) ) return kVus;

     } else

     if( pdg::IsDQuark(init_pdgc) ) {

         if( pdg::IsUQuark(fin_pdgc) ) return kVud;
         if( pdg::IsCQuark(fin_pdgc) ) return kVcd;

     } else

     if( pdg::IsSQuark(init_pdgc) ) {

         if( pdg::IsUQuark(fin_pdgc) ) return kVus;
         if( pdg::IsCQuark(fin_pdgc) ) return kVcs;

     } else {

       LOG("Utils", pERROR)
                << "CKM element for " << init_pdgc 
                                     << " --> " << fin_pdgc << " was not set";
    }
     return 0;
  }
  //__________________________________________________________________________
};

}        // genie namespace
#endif   // _UTILS_H_
