//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecParams

\brief    

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep 10, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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

  double DefValue       (GSyst_t s) const;       ///< default value of a physics parameter
  double CurValue       (GSyst_t s) const;       ///< current value of a physics parameter
  double CurTwkDial     (GSyst_t s) const;       ///< current phys param tweaking dial value
  bool   IsIncluded     (GSyst_t s) const;       ///< is included?
  bool   IsTweaked      (GSyst_t s) const;       ///< is included & tweaked to non-default value?
  bool   IsTweaked      (void)      const;       ///< is any parameter tweaked?
  double ChisqPenalty   (void)      const;       ///< chi^2_{penalty} for the current parameter shifts
  void   Reconfigure    (void)      const;       ///< propagate physics parameter changes to GENIE
  void   LoadDefaults   (void);                  ///< load default params from GENIE
  void   Reset          (GSyst_t s);             ///< set cur=def, twk_dial=0, tweaked=false
  void   Reset          (void);                  ///< reset all cross section params
  void   SetCurTwkDial  (GSyst_t s, double val); ///<
      
 private:

  void   SetDefValue    (GSyst_t s, double val); ///<
  void   SetCurValue    (GSyst_t s, double val); ///<
  void   SetTweakedFlag (GSyst_t s, bool   val); ///<
  void   Init           (void);                  ///<

  map<GSyst_t, double> fDefParams;  ///< default values for all physics params supported by ReWeight
  map<GSyst_t, double> fCurParams;  ///< current values for all physics params included by the user
  map<GSyst_t, double> fCurTwkDial; ///< corresponding tweaking dial, current = default * (1 + dial * fractional_err)
  map<GSyst_t, bool>   fIsTweaked;  ///<
 };

} // rew   namespace
} // genie namespace

#endif

