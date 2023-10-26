//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Noah Steinberg <nsteinbe \at fnal.gov>
 Fermi National Accelerator Laboratory 
*/
//____________________________________________________________________________

#include "Physics/Fortran/XSection/HadronTensorInterfaceI.h"

//____________________________________________________________________________
std::complex<double> genie::HadronTensorInterfaceI::operator()(genie::TensorIndex_t mu,
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

