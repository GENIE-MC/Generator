//____________________________________________________________________________
/*

   \class    genie::IASingleNucleonTensor

   \brief    Concrete implementation of HadronTensorI interface

   Impulse approximation single nucleon response tensor

   \created  Oct 20, 2023
   \author   Noah Steinberg <nsteinbe \at fnal.gov>
   Fermi National Accelerator Laboratory 

   Native C++ implementation, replacing the previous Fortran wrapper.

   \updated May 02, 2026
   \author Liang Liu <liangliu \at fnal.gov>
   Fermi National Accelerator Laboratory

   \cpright  Copyright (c) 2003-2026, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _IA_SINGLE_NUCLEON_TENSOR_H_
#define _IA_SINGLE_NUCLEON_TENSOR_H_


// GENIE includes
#include "Physics/HadronTensors/NucleonTensor.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/HadronTensors/Rank2LorentzTensor.h"

namespace genie {

  class IASingleNucleonTensor : public NucleonTensor {

    public:

      // Constructor takes in hadron information and form factors
      // as well as name of model for interface
      IASingleNucleonTensor(double w, double mNi, 
          TLorentzVector p4Ni, TLorentzVector p4Nf, 
          const genie::QELFormFactors& FormFactors, 
          std::string model);

      // Overridden initialization function
      virtual void initialize_tensor(std::complex<double> (&hadron_tensor)[4][4]) const override;

      ~IASingleNucleonTensor() {}

    protected:
      TLorentzVector fp4Ni;
      TLorentzVector fp4Nf;
      TLorentzVector fq4t;
      double fmNi;
      double fw;
      const QELFormFactors fFormFactors;
      std::string fModel;
  };

}
#endif	
