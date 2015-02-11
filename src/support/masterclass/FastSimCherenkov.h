//____________________________________________________________________________
/*!

\class    genie::masterclass::FastSimCherenkov

\brief    Fast simulation of the response of a Cherenkov detector

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Sep 22, 2010

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

