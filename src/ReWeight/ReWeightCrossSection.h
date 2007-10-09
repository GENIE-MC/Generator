//____________________________________________________________________________
/*!

\class   genie::ReWeightCrossSection

\brief   A cross section model reweighting engine.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#ifndef _REWEIGHT_CROSS_SECTION_H_
#define _REWEIGHT_CROSS_SECTION_H_

#include <map>

#include "Conventions/KinePhaseSpace.h"
#include "Interaction/ScatteringType.h"
#include "EVGCore/InteractionList.h"
#include "EVGDrivers/GEVGPool.h"

using std::map;

namespace genie {

class EventRecord;

class ReWeightCrossSection {

public :
  ReWeightCrossSection();
 ~ReWeightCrossSection();

  // various configuration options for the cross section 
  // re-weighting engine

  void HandleInitState  (const InitialState & init_state);
  void DiffCrossSecType (ScatteringType_t sct, KinePhaseSpace_t kps);
  void DontReweight     (const Interaction & interaction);

  // the actual re-weighting code

  double NewWeight (const EventRecord & event);

private:

   void Initialize();

   GEVGPool                                fGPool;             ///<
   map<ScatteringType_t, KinePhaseSpace_t> fCrossSecModelPhSp; ///<
   InteractionList                         fNoRewProc;         ///<
};

}      // genie namespace

#endif // _REWEIGHT_CROSS_SECTION_H_
