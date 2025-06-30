//____________________________________________________________________________
/*!

\namespace  genie::utils

\brief      Common functions used for handling generation of the primary
            lepton, regardless of whether the relevant class inherits from
            PrimaryLeptonGenerator or not.

\author     Steven Gardiner <gardiner \at fnal.gov>
            Fermi National Accelerator Laboratory

\created    May 01, 2020

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PRIMARY_LEPTON_UTILS_H
#define _PRIMARY_LEPTON_UTILS_H

namespace genie {

class GHepRecord;

namespace utils {

  void SetPrimaryLeptonPolarization( GHepRecord* ev );

} // utils   namespace
} // genie   namespace

#endif // _PRIMARY_LEPTON_UTILS_H
