//____________________________________________________________________________
/*!

\class    genie::HadronTensorModelI

\brief    Creates hadron tensor objects for use in cross section calculations

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  April 26, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HADRON_TENSOR_MODEL_H_
#define _HADRON_TENSOR_MODEL_H_

// standard library includes
#include <map>
#include <string>
#include <vector>

// GENIE includes
#include "Framework/Algorithm/Algorithm.h"
#include "Physics/HadronTensors/HadronTensorI.h"

namespace genie {

class HadronTensorModelI : public Algorithm {

public:
  virtual ~HadronTensorModelI();

  /// Retrieves a pointer to a hadron tensor object appropriate for this model
  /// \param[in] tensor_pdg The PDG code for the nuclide described by the tensor
  /// \param[in] type The desired kind of hadron tensor
  /// \returns A pointer to the requested hadron tensor, or NULL if a match
  /// could not be found/created
  virtual const HadronTensorI* GetTensor(int tensor_pdg, HadronTensorType_t type) const = 0;

protected:
  HadronTensorModelI();
  HadronTensorModelI(std::string name);
  HadronTensorModelI(std::string name, std::string config);
};

} // namespace genie

#endif // _HADRON_TENSOR_MODEL_H_
