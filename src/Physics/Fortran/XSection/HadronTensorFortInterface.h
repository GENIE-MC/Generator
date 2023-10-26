//____________________________________________________________________________
/*

\class    genie::HadronTensorFortInterface

\brief    Concrete implementation of HadronTensorInterfaceI interface
          that calls fortran subroutines to compute response tensor/

\author   Noah Steinberg <nsteinbe \at fnal.gov>
	  Fermi National Accelerator Laboratory 

\created  Oct 20, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
 	   For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _HADRON_TENSOR_FORT_INTERFACE_H
#define _HADRON_TENSOR_FORT_INTERFACE_H


// GENIE includes
#include "Physics/Fortran/XSection/HadronTensorInterfaceI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/Rank2LorentzTensorI.h"

namespace genie {

class HadronTensorFortInterface : public HadronTensorInterfaceI {

public:
	
        // Constructor takes in hadron information and form factors
        // as well as name of model for interface
	HadronTensorFortInterface(double w, double mNi, TLorentzVector p4Ni, TLorentzVector p4Nf, const genie::QELFormFactors& FormFactors, std::string model);

	// Overridden initialization function
	virtual void initialize_tensor(std::complex<double> (&hadron_tensor)[4][4]) const override;

	~HadronTensorFortInterface() {}

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
