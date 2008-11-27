//____________________________________________________________________________
/*!

\namespace  genie::utils::xgen

\brief      Utilities for cross-generator comparisons

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            STFC, Rutherford Appleton Laboratory

\created    Nov 30, 2008

\cpright    Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XGEN_UTILS_H_
#define _XGEN_UTILS_H_

namespace genie {

class GHepRecord;

namespace utils {
namespace xgen  {

  int NeutReactionCode   (const GHepRecord * evrec); 
  int NuanceReactionCode (const GHepRecord * evrec); 

} // xgen  namespace
} // utils namespace
} // genie namespace

#endif // _XGEN_UTILS_H_
