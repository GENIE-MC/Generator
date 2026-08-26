//____________________________________________________________________________
/*
\class    genie::NucleonTensor

\brief    Pure abstract base class. Defines the NucleonTensor interface
          to be implemented by any algorithmic class computing the hadronic 
	        response tensor


 Noah Steinberg <nsteinbe \at fnal.gov>
 Liang Liu <liangliu \at fnal.gov>
 Fermi National Accelerator Laboratory 

 Copyright (c) 2003-2026, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include "Physics/HadronTensors/NucleonTensor.h"

//____________________________________________________________________________
std::complex<double> genie::NucleonTensor::operator()(genie::TensorIndex_t mu,
  genie::TensorIndex_t nu) const 
{

  // Check to see if the hadron tensor has already been filled
  if (fCreated_Tensor == false) {
    
    // Create tensor if it hasn't been
    // filled before
    initialize_tensor(fhadron_tensor);
    fCreated_Tensor = true;
  }

  // Return tensor elements
  std::complex<double> result = fhadron_tensor[mu][nu];
  
  return result;

}

