//____________________________________________________________________________
/*!

\class    genie::masterclass::FastSimScintCalo

\brief    Fast simulation of the response of a scintillator calorimeter.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Sep 22, 2010

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
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
