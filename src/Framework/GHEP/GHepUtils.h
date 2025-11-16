//____________________________________________________________________________
/*!

\namespace  genie::utils::ghep

\brief      GHEP event record utilities

\author     Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

\created    Nov 30, 2008

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _GHEP_UTILS_H_
#define _GHEP_UTILS_H_

namespace genie {

class GHepRecord;

namespace utils {
namespace ghep  {

  int NeutReactionCode   (const GHepRecord * evrec);
  int NuanceReactionCode (const GHepRecord * evrec);

} // ghep  namespace
} // utils namespace
} // genie namespace

#endif // _GHEP_UTILS_H_
