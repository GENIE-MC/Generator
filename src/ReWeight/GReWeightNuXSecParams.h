//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecParams

\brief    

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep 1, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NEUTRINO_CROSS_SECTION_PARAMS_H_
#define _G_REWEIGHT_NEUTRINO_CROSS_SECTION_PARAMS_H_

#include <map>

#include "ReWeight/GSyst.h"

using std::map;

namespace genie {
namespace rew   {

 class GReWeightNuXSecParams
 {
 public:
   GReWeightNuXSecParams();
  ~GReWeightNuXSecParams();

  double DefValue   (GSyst_t) const; ///< default value of a physics parameter
  double CurValue   (GSyst_t) const; ///< current value of a physics parameter
  double CurTwkDial (GSyst_t) const; ///< current phys param tweaking dial value
  bool   IsIncluded (GSyst_t) const; ///< is included?
  bool   IsTweaked  (GSyst_t) const; ///< is included & tweaked to non-default value?
  bool   IsTweaked  (void)    const; ///< is any parameter tweaked?

  void Reset          (GSyst_t);     ///< set cur=def, twk_dial=0, tweaked=false
  void SetDefValue    (GSyst_t syst, double value);
  void SetCurValue    (GSyst_t syst, double value);
  void SetCurTwkDial  (GSyst_t syst, double value);
  void SetTweakedFlag (GSyst_t syst, bool   value);
      
 private:

   void Init(void);

   map<GSyst_t, double> fDefParams;
   map<GSyst_t, double> fCurParams;
   map<GSyst_t, double> fCurTwkDial;
   map<GSyst_t, bool>   fIsTweaked;
 };

} // rew   namespace
} // genie namespace

#endif

