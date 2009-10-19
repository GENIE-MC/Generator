
//____________________________________________________________________________
/*!

\namespace genie::rew::margin

\brief     

\author    Jim Dobson <J.Dobson07 \at imperial.ac.uk>
           Imperial College London

           Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           STFC, Rutherford Appleton Laboratory

\created   Sep 09, 2009

\cpright   Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RW_MRG_H_
#define _RW_MRG_H_

//#include <Math/IFunction.h>

#include "Conventions/GBuild.h"
#include "EVGCore/EventRecord.h"
#include "ReWeight/GSyst.h"
#include "ReWeight/GReWeight.h"

namespace genie  {
namespace rew    {
namespace margin {


  double MarginalizeFates (
     double fixed_dial, genie::rew::GSyst_t fixed_syst, const genie::EventRecord * ev, int n=200);

/*           
  //
  // utility class for integration (nuisance param marginalization) using GSL
  //
  class INukeFateIntgrd: public ROOT::Math::IBaseFunctionMultiDim 
  {
   public:
    INukeFateIntgrd(
      double fixed_dial, genie::rew::GSyst_t fixed_syst, const genie::EventRecord * ev);
   ~INukeFateIntgrd();

    unsigned int                        NDim   (void)               const;
    double                              DoEval (const double * xin) const;
    ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

   private:
    double                     fFixedTwKDial;
    genie::rew::GSyst_t        fFixedFateSyst;
    const genie::EventRecord * fEvent;
    genie::rew::GReWeight *    fRW;
  };
*/

}  // margin namespace
}  // rew    namespace
}  // genie  namespace

#endif // _RW_MRG_H_
