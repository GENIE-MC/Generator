//____________________________________________________________________________
/*!

\namespace  genie::utils::ghep

\brief      GHEP event record utilities

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    Nov 30, 2008

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
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
