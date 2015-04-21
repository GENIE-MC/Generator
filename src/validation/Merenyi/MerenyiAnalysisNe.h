//____________________________________________________________________________
/*!

\class   MerenyiAnalysisNe

\brief   Merenyi Neon analysis class

\author  Pauli Kehayias (Tufts Univ)

\created Jan 13, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MERENYI_ANALYSIS_NEON_H_
#define _MERENYI_ANALYSIS_NEON_H_

#include "validation/Merenyi/MerenyiAnalysis.h"

namespace genie       {
namespace vld_merenyi {

class MerenyiAnalysisNe : public MerenyiAnalysis {

public:
  MerenyiAnalysisNe() { }
 ~MerenyiAnalysisNe() { }

private:
  void classifyHadrons (MerenyiNuEvent& myEvent);
  void checkCandP      (MerenyiNuEvent& myEvent);
};


} // vld_merenyi
} // genie

#endif


