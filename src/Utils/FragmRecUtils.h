//____________________________________________________________________________
/*!

\namespace  genie::fragm_rec_utils

\brief      Simple utilities for the Fragmentation Event Record.

            The Fragmentation event record is a TClonesArray of TMCParticles -
            equivalent to PYTHIA's PYJETS.

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    November 26, 2004

*/
//____________________________________________________________________________

#ifndef _FRAGM_REC_UTILS_H_
#define _FRAGM_REC_UTILS_H_

#include <TClonesArray.h>

namespace genie {

namespace fragm_rec_utils
{
  int NParticles(int pdg_code, const TClonesArray * const particle_list);
  int NParticles(int pdg_code, int status, const TClonesArray * const particle_list);
  int NPositives(const TClonesArray * const particle_list);
  int NNegatives(const TClonesArray * const particle_list);

  void Print(const TClonesArray * const part_list);

}      // fragm_rec_utils namespace

}      // genie namespace

#endif // _FRAGM_REC_UTILS_H_
