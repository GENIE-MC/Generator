//____________________________________________________________________________
/*!

\class   MerenyiAnalysisD2

\brief   Merenyi D2 analysis class

\author  Pauli Kehayias (Tufts Univ)

\created Jan 13, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MERENYI_ANALYSIS_D2_H_
#define _MERENYI_ANALYSIS_D2_H_

#include "validation/Merenyi/MerenyiAnalysis.h"

namespace genie       {
namespace vld_merenyi {

class MerenyiAnalysisD2 : public MerenyiAnalysis {
public:
  MerenyiAnalysisD2() {}
 ~MerenyiAnalysisD2() {}

private:
  void classifyHadrons(MerenyiNuEvent& myEvent);
};

} // vld_merenyi
} // genie

#endif


