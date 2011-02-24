//____________________________________________________________________________
/*!

\class    genie::gview::fastsim::FastSimCherenkov

\brief    Fast simulation of the response of a Cherenkov detector

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep 22, 2010

*/
//____________________________________________________________________________

#ifndef _CHERENKOV_MC_H_
#define _CHERENKOV_MC_H_

class TRootEmbeddedCanvas;

namespace genie {

 class EventRecord;

 namespace gview {
  namespace fastsim {

    class FastSimCherenkov {
    public:
      FastSimCherenkov();
     ~FastSimCherenkov();
      void SetEmbeddedCanvas (TRootEmbeddedCanvas * ec);
      void Draw              (EventRecord * event);
    private:
      TRootEmbeddedCanvas * fEmbeddedCanvas;
    };

  }  // fastsim namespace
 }  // gview namespace
}  // genie namespace

#endif  // _CHERENKOV_MC_H_	

