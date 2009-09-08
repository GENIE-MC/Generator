//____________________________________________________________________________
/*!

\class    genie::rew::GSystType_t

\brief    An enumeration of systematic parameter types

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_SYSTEMATIC_TYPE_H_
#define _G_SYSTEMATIC_TYPE_H_

namespace genie {
namespace rew   {

typedef enum GSystType {

  kSystType_Null = 0, 
  kSystType_NuXSec,      ///< neutrino cross section systematics
  kSystType_INuke,       ///< intranuclear rescattering systematics

} GSystType_t;

} // rew   namespace
} // genie namespace

#endif 

