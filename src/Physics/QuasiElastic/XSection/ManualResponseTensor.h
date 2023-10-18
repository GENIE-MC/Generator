//____________________________________________________________________________
/*!
\class    genie::ManualResponseTensor
\brief    Constructs full Nuclear Response Tensor W^{\mu\nu} from an input 4x4 array
\author   Noah Steinberg <nsteinbe \at fnal.gov>
          Fermi National Accelerator Laboratory
\created  Nov 23, 2021
\cpright  Copyright (c) 2003-2021, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef MANUAL_RESPONSE_TENSOR_H
#define MANUAL_RESPONSE_TENSOR_H

#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InteractionType.h"
#include "Physics/QuasiElastic/XSection/Rank2LorentzTensorI.h"

namespace genie {

class ManualResponseTensor : public Rank2LorentzTensorI {

public:
	std::complex<double> RespTensor[4][4];
	/*https://stackoverflow.com/questions/13054243/set-one-array-equal-to-another-without-a-loop/13054295*/
	ManualResponseTensor(const std::complex<double> in_RespArray[4][4]); 

	inline virtual ~ManualResponseTensor() {}

	 /// Retrieves a tensor element corresponding to the given indices
  	virtual std::complex<double> operator()(TensorIndex_t mu, TensorIndex_t nu)
  	const /*override*/;

}; // class ManualResponseTensor
 
} // genie namespace
#endif
