//____________________________________________________________________________
/*!

\class    genie::rew::GSystUncertainty

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Sep 1, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_SYST_UNCERTAINTY_H_
#define _G_SYST_UNCERTAINTY_H_

#include "ReWeight/GSyst.h"

namespace genie {
namespace rew   {

class GSystUncertainty {

public:  
  static GSystUncertainty * Instance (void);

  double OneSigmaErr(GSyst_t syst) const;

private:
  GSystUncertainty();
  GSystUncertainty(const GSystUncertainty & err);
 ~GSystUncertainty();
  
  static GSystUncertainty * fInstance;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (GSystUncertainty::fInstance !=0) {
            delete GSystUncertainty::fInstance;
            GSystUncertainty::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

} // rew   namespace
} // genie namespace

#endif 

