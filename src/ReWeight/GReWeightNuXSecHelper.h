//____________________________________________________________________________
/*!

\class   genie::rew::GReWeightNuXSecHelper

\brief   Helper class for cross section model reweighting

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NEUTRINO_CROSS_SECTION_HELPER_H_
#define _G_REWEIGHT_NEUTRINO_CROSS_SECTION_HELPER_H_

#include <map>

#include "Conventions/KinePhaseSpace.h"
#include "Interaction/ScatteringType.h"
#include "EVGDrivers/GEVGPool.h"

using std::map;

namespace genie {

class EventRecord;

namespace rew {

class GReWeightNuXSecHelper {

public :
  GReWeightNuXSecHelper();
 ~GReWeightNuXSecHelper();

  void   HandleInitState  (const InitialState & init_state);            
  void   DiffCrossSecType (ScatteringType_t sct, KinePhaseSpace_t kps); 
  double NewWeight        (const EventRecord & event, bool shape_only = false);

private:

   void Initialize();

   GEVGPool                                fGPool;             ///<
   map<ScatteringType_t, KinePhaseSpace_t> fCrossSecModelPhSp; ///<
};

}      // rew   namespace
}      // genie namespace

#endif // _G_REWEIGHT_NEUTRINO_CROSS_SECTION_H_
