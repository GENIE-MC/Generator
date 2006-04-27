//____________________________________________________________________________
/*!

\namespace  genie::utils::fragmrec

\brief      Simple utilities for the Fragmentation Event Record.

            The Fragmentation event record is a TClonesArray of TMCParticles -
            equivalent to PYTHIA's PYJETS.

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    November 26, 2004

\cpright    Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
            All rights reserved.
            For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _FRAGM_REC_UTILS_H_
#define _FRAGM_REC_UTILS_H_

#include <TClonesArray.h>

namespace genie {
namespace utils {

namespace fragmrec
{
  int NParticles(int pdg_code, const TClonesArray * const particle_list);
  int NParticles(int pdg_code, int status, const TClonesArray * const particle_list);
  int NPositives(const TClonesArray * const particle_list);
  int NNegatives(const TClonesArray * const particle_list);

  void Print(const TClonesArray * const part_list);

} // fragmrec namespace
} // utils    namespace
} // genie    namespace

#endif // _FRAGM_REC_UTILS_H_
