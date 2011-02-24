//____________________________________________________________________________
/*!

\class    genie::gview::fastsim::FastSimScintCalo

\brief    Fast simulation of the response of a scintillator calorimeter.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep 22, 2010

*/
//____________________________________________________________________________

#ifndef _SCINT_CALO_SIM_H_
#define _SCINT_CALO_SIM_H_

class TRootEmbeddedCanvas;

namespace genie {

 class EventRecord;

 namespace gview {
  namespace fastsim {

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

  }  // fastsim namespace
 }  // gview namespace
}  // genie namespace

#endif  // _SCINT_CALO_SIM_H_

