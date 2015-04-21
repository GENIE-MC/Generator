//____________________________________________________________________________
/*!

\namespace genie::rew::margin

\brief     Utilities for nuisance parameter marginalization

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Lab

\created   Oct 20, 2009

\cpright   Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RW_MRG_H_
#define _RW_MRG_H_

#include "Conventions/GBuild.h"
#include "EVGCore/EventRecord.h"
#include "ReWeight/GSyst.h"
#include "ReWeight/GReWeight.h"

namespace genie  {
namespace rew    {
namespace margin {

  double MarginalizeFates (
      double fixed_dial, genie::rew::GSyst_t fixed_syst, 
      const genie::EventRecord * ev, int n=200);

}  // margin namespace
}  // rew    namespace
}  // genie  namespace

#endif // _RW_MRG_H_
