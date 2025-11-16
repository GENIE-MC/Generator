//____________________________________________________________________________
/*!

\class    genie::masterclass::FastSimCherenkov

\brief    Fast simulation of the response of a Cherenkov detector

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Sep 22, 2010

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org           
*/
//____________________________________________________________________________

#ifndef _CHERENKOV_MC_H_
#define _CHERENKOV_MC_H_

#include <TRootEmbeddedCanvas.h>

namespace genie {

 class EventRecord;

 namespace masterclass {

    class FastSimCherenkov {
    public:
      FastSimCherenkov();
     ~FastSimCherenkov();
      void SetEmbeddedCanvas (TRootEmbeddedCanvas * ec);
      void Draw              (EventRecord * event);
    private:
      TRootEmbeddedCanvas * fEmbeddedCanvas;
    };

 }  // masterclass namespace
}  // genie namespace

#endif  // _CHERENKOV_MC_H_
