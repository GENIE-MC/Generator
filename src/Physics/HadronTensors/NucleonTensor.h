//____________________________________________________________________________
/*

\class    genie::NucleonTensor

\brief    Pure abstract base class. Defines the NucleonTensor interface
          to be implemented by any algorithmic class computing the hadronic 
	        response tensor

\author   Noah Steinberg <nsteinbe \at fnal.gov>
\author   Liang Liu <liangliu \at fnal.gov>

\created  Oct 20, 2023

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _NUCLEON_TENSOR_H_
#define _NUCLEON_TENSOR_H_

#include "Physics/HadronTensors/Rank2LorentzTensor.h"

namespace genie {

class NucleonTensor : public Rank2LorentzTensor {

public:
  virtual ~NucleonTensor() = default;

  //! Compute individual elements of tensor
  virtual std::complex<double> operator()(genie::TensorIndex_t mu,
	 genie::TensorIndex_t nu) const override;  

  //! Initialize the tensor object, pure virtual 
  virtual void initialize_tensor(std::complex<double> (&hadron_tensor)[4][4]) const = 0;

private:  
  mutable bool fCreated_Tensor = false;
  mutable std::complex<double> fhadron_tensor[4][4];
};

}         // genie namespace
#endif    // _NUCLEON_TENSOR_H_
