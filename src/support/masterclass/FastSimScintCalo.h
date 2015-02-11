//____________________________________________________________________________
/*!

\class    genie::masterclass::FastSimScintCalo

\brief    Fast simulation of the response of a scintillator calorimeter.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Sep 22, 2010

*/
//____________________________________________________________________________

#ifndef _SCINT_CALO_SIM_H_
#define _SCINT_CALO_SIM_H_

#include <TRootEmbeddedCanvas.h>

namespace genie {

 class EventRecord;

 namespace masterclass {

    class FastSimScintCalo {
    public:
      FastSimScintCalo();
     ~FastSimScintCalo();
      void SetEmbeddedCanvas (TRootEmbeddedCanvas * ec);
  //  void SetNumOfScintLayers()
  //  void SetNumOfScintStripsPerLayer()
      void Draw              (EventRecord * event);
    private:
      TRootEmbeddedCanvas * fEmbeddedCanvas;
    };

 }  // masterclass namespace
}  // genie namespace

#endif  // _SCINT_CALO_SIM_H_

