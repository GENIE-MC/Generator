//____________________________________________________________________________
/*!

\class    genie::MartiniQELHadronTensorModel

\brief    Creates hadron tensor objects for calculations of quasielastic
          cross sections using the SuSAv2 approach

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  April 26, 2019

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MARTINI_QEL_HADRON_TENSOR_MODEL_H_
#define _MARTINI_QEL_HADRON_TENSOR_MODEL_H_

// GENIE includes
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"

namespace genie {

class MartiniQELHadronTensorModel : public TabulatedHadronTensorModelI {

public:

  MartiniQELHadronTensorModel();
  MartiniQELHadronTensorModel(std::string config);

  virtual ~MartiniQELHadronTensorModel();

protected:

  // Implementation of TabulatedHadronTensorModelI interface
  virtual HadronTensorI* ParseTensorFile( const std::string& full_file_name ) const;

};

} // namespace genie

#endif // _MARTINI_QEL_HADRON_TENSOR_MODEL_H_
