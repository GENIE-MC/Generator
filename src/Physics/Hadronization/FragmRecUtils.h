//____________________________________________________________________________
/*!

\namespace  genie::utils::fragmrec

\brief      Simple utilities for the Fragmentation Event Record.

            The Fragmentation event record is a TClonesArray of TMCParticles -
            equivalent to PYTHIA's PYJETS.

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    November 26, 2004

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
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
