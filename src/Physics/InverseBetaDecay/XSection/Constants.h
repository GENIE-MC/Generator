//____________________________________________________________________________
/*!

\namespace genie::VLEConstants

\brief     Constants used in VLE package

\author    Corey Reed <cjreed \at nikhef.nl>
           Nikhef

\created   June 22, 2009

\cpright   Copyright (c) 2003-2019, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _VLECONSTANTS_H_
#define _VLECONSTANTS_H_

#include <TMath.h>

#include "Framework/Conventions/Constants.h"

namespace genie {
namespace constants {

// masses
static const double k4NucMass2      = 4.000 * kNucleonMass2;
static const double k4EleMass2      = 4.000 * kElectronMass2;
static const double kNucMassDiff    = kNeutronMass - kProtonMass;
static const double kNucMassDiff2   = kNucMassDiff * kNucMassDiff;

} // namespace constants
} // namespace genie

#endif // _VLECONSTANTS_H_
