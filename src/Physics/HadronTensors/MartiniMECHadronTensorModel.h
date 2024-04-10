//____________________________________________________________________________
/*!

\class    genie::MartiniMECHadronTensorModel

\brief    Creates hadron tensor objects for calculations of MEC
          cross sections using the Martini approach

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  April 26, 2019

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MARTINI_MEC_HADRON_TENSOR_MODEL_H_
#define _MARTINI_MEC_HADRON_TENSOR_MODEL_H_

// GENIE includes
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"

namespace genie {

class MartiniMECHadronTensorModel : public TabulatedHadronTensorModelI {

public:

  MartiniMECHadronTensorModel();
  MartiniMECHadronTensorModel(std::string config);

  virtual ~MartiniMECHadronTensorModel();

protected:

  // Implementation of TabulatedHadronTensorModelI interface
  virtual HadronTensorI* ParseTensorFile( const std::string& full_file_name ) const;

};

} // namespace genie

#endif // _Martini_MEC_HADRON_TENSOR_MODEL_H_
