//____________________________________________________________________________
/*!

\class    genie::SuSAv2MECHadronTensorModel

\brief    Creates hadron tensor objects for calculations of MEC
          cross sections using the SuSAv2 approach

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  April 26, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SUSAV2_MEC_HADRON_TENSOR_MODEL_H_
#define _SUSAV2_MEC_HADRON_TENSOR_MODEL_H_

// GENIE includes
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"

namespace genie {

class SuSAv2MECHadronTensorModel : public TabulatedHadronTensorModelI {

public:

  SuSAv2MECHadronTensorModel();
  SuSAv2MECHadronTensorModel(std::string config);

  virtual ~SuSAv2MECHadronTensorModel();

protected:

  // Implementation of TabulatedHadronTensorModelI interface
  virtual HadronTensorI* ParseTensorFile( const std::string& full_file_name ) const;

};

} // namespace genie

#endif // _SUSAV2_MEC_HADRON_TENSOR_MODEL_H_
