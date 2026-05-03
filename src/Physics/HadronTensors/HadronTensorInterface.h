//____________________________________________________________________________
/*

   \class    genie::HadronTensorInterface

   \brief    Concrete implementation of HadronTensorInterfaceI interface
   that calls fortran subroutines to compute response tensor/

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

#ifndef _HADRON_TENSOR_FORT_INTERFACE_H
#define _HADRON_TENSOR_FORT_INTERFACE_H


// GENIE includes
#include "Physics/HadronTensors/HadronTensorInterfaceI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/HadronTensors/Rank2LorentzTensorI.h"

namespace genie {

  class HadronTensorInterface : public HadronTensorInterfaceI {

    public:

      // Constructor takes in hadron information and form factors
      // as well as name of model for interface
      HadronTensorInterface(double w, double mNi, TLorentzVector p4Ni, TLorentzVector p4Nf, const genie::QELFormFactors& FormFactors, std::string model);

      // Overridden initialization function
      virtual void initialize_tensor(std::complex<double> (&hadron_tensor)[4][4]) const override;

      ~HadronTensorInterface() {}

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
