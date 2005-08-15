//____________________________________________________________________________
/*!

\class    genie::GHepOrder

\brief    A very simple utility class to make more transparent the usage of
          the standard GHEP record ordering conventions

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 16, 2005

*/
//____________________________________________________________________________

#ifndef _GHEP_ORDER_H_
#define _GHEP_ORDER_H_

#include "Interaction/Interaction.h"

namespace genie {

class GHepOrder {

public:

  //_________________________________________________________________________
  static int ProbePosition(void)
  {
    // the probe is *always* at slot 0
    return 0;
  }
  //_________________________________________________________________________
  static int TargetNucleusPosition(const Interaction * const interaction)
  {
    // the target nucleus (if the interaction happens at a nuclear target)
    // will always be at 1

    const InitialState & init_state = interaction->GetInitialState();
    bool tgt_is_nucleus = init_state.GetTarget().IsNucleus();

    int target_nucleus_pos = (tgt_is_nucleus) ? 1 : -1 /*no nucl. tgt*/;

    return target_nucleus_pos;
  }
  //_________________________________________________________________________
  static int StruckNucleonPosition(const Interaction * const interaction)
  {
    // the struck nucleon (if any) is either at slot 2 (for nuclear targets)
    // or slot 1 (for free nucleon targets)

    const ProcessInfo & proc_info = interaction->GetProcessInfo();
    if(proc_info.IsCoherent()) return -1; // no struck nucleon

    const InitialState & init_state = interaction->GetInitialState();
    bool tgt_is_nucleus = init_state.GetTarget().IsNucleus();

    int struck_nucleon_pos = (tgt_is_nucleus) ? 2 : 1;

    return struck_nucleon_pos;
  }
  //_________________________________________________________________________
};

}         // genie namespace
#endif    // _GHEP_ORDER_H_

