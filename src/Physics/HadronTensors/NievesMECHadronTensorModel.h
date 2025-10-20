//____________________________________________________________________________
/*!

\class    genie::NievesMECHadronTensorModel

\brief    Creates hadron tensor objects for calculations of MEC
          cross sections using the Valencia model

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  April 26, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NIEVES_MEC_HADRON_TENSOR_MODEL_H_
#define _NIEVES_MEC_HADRON_TENSOR_MODEL_H_

// GENIE includes
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"

namespace genie {

class NievesMECHadronTensorModel : public TabulatedHadronTensorModelI {

public:

  NievesMECHadronTensorModel();
  NievesMECHadronTensorModel(std::string config);

  virtual ~NievesMECHadronTensorModel();

protected:

  // Implementation of TabulatedHadronTensorModelI interface
  virtual HadronTensorI* ParseTensorFile( const std::string& full_file_name ) const;

};

} // namespace genie

#endif // _NIEVES_MEC_HADRON_TENSOR_MODEL_H_
