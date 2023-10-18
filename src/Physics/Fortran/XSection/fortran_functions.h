//____________________________________________________________________________
//
/*
  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Noah Steinberg <nsteinbe \at fnal.gov>
  Fermi National Acclerator Laboratory


*/
//____________________________________________________________________________ 

// For each additional Fortran file added to compute response tensor
// Must declare the corresponding function which returns the tensor
// as 'extern "C"' to expose fortran subroutine to c++
extern "C" {

  void compute_hadron_tensor_SF(double *mn,double *w, double *wt,double *pNi_x,double *pNi_y,
        double *pNi_z, double *q_x, double *q_y, double *q_z, double *f1v, double *f2v,
        double *ffa, double *ffp, std::complex<double> HadronTensor[4][4]);
}

extern "C" { 

  void compute_hadron_tensor_SF_CC(double *mn,double *w, double *wt,double *pNi_x,double *pNi_y,
        double *pNi_z, double *q_x, double *q_y, double *q_z, double *f1v, double *f2v,
        double *ffa, double *ffp, std::complex<double> HadronTensor[4][4]);
}

// This function calls the above fortran subroutines
// based on an identifying string
// string is filled in UnifiedQELPXsec config file
void Get_hadrontensor_from_fortran(std::string TensorModel, double xmn, double w, double wt, double pNi_x, double pNi_y, double pNi_z, double q_x, double q_y, double q_z, double f1v, double f2v, double ffa, double ffp, std::complex<double> HadronTensor[4][4]) {

        if(TensorModel == "compute_hadron_tensor_SF") {
                compute_hadron_tensor_SF(&xmn, &w, &wt, &pNi_x, &pNi_y, &pNi_z, &q_x, &q_y, &q_z, &f1v, &f2v, &ffa, &ffp, HadronTensor);
        }

        else if(TensorModel == "compute_hadron_tensor_SF_CC") {
                compute_hadron_tensor_SF_CC(&xmn, &w, &wt, &pNi_x, &pNi_y, &pNi_z, &q_x, &q_y, &q_z, &f1v, &f2v, &ffa, &ffp, HadronTensor);
        }

        else std::cout << "Unrecognized Fortran Tensor Model in genie::UnifiedQELPXSec::XSec()" << "\n";
}

	
