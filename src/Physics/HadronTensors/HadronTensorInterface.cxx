//____________________________________________________________________________
/*
   Copyright (c) 2003-2023, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org

   Modify the Fortran wrapper (by Noah Steinberg) to C++ native implemetation of one body current
   Native C++ implementation, replacing the previous Fortran wrapper.

   By Liang Liu and ChatGPT

   Noah Steinberg <nsteinbe \at fnal.gov>
   Fermi National Accelerator Laboratory 


   Liang Liu <liangliu \at fnal.gov>
   Fermi National Accelerator Laboratory 

*/    
//____________________________________________________________________________
#include "Physics/HadronTensors/HadronTensorInterface.h"
#include "Physics/HadronTensors/onebody_currents_sf.h"

genie::HadronTensorInterface::HadronTensorInterface(double w, double mNi, TLorentzVector p4Ni, TLorentzVector p4Nf, const genie::QELFormFactors& FormFactors, std::string model) : fFormFactors(FormFactors) 
{

  // Get Kinematics fo interaction
  fmNi = mNi;
  fw = w;

  // On shell energy of struck particle
  double E_NiOnShell = std::sqrt(p4Ni.P()*p4Ni.P() + mNi*mNi);
  //double epsilon_B = E_NiOnShell - p4Ni.E();

  // On shell four vector
  TLorentzVector p4NiOnShell = TLorentzVector(p4Ni.Vect(), E_NiOnShell);
  fp4Ni = p4NiOnShell;

  // Final state particle
  fp4Nf = p4Nf;

  // 4-momentum transfer
  fq4t = p4Nf - p4NiOnShell; 

  // Set hadron response tensor model
  fModel = model;
}

void genie::HadronTensorInterface::initialize_tensor(std::complex<double> (&hadron_tensor)[4][4]) const
{

  // Gather all kinematic and form factor info
  // TO DO: TLorentzVector and FormFactor objects rvalues can't
  // be passed to f90 as they don't have an lvalue
  double mNi, w;
  double pNix, pNiy, pNiz;
  double qtx, qty, qtz, wt;
  double f1v, xif2v, fa, fp;

  w = fw; mNi = fmNi;
  pNix = fp4Ni.X(); pNiy = fp4Ni.Y(); pNiz = fp4Ni.Z();
  qtx = fq4t.X(); qty = fq4t.Y(); qtz = fq4t.Z(); wt = fq4t.E();
  f1v = fFormFactors.F1V(); xif2v = fFormFactors.xiF2V(); fa = fFormFactors.FA(); fp = fFormFactors.Fp();

  // Call appropriate fortran subroutine based on name of model provided
  if (fModel == "Noemi-hadron-tensor") {
    genie::onebody_currents_sf::compute_hadron_tensor_eigen(mNi, w, wt, pNix, pNiy, pNiz, qtx, qty,
        qtz, f1v, xif2v, fa, fp, hadron_tensor);
  }
  else if (fModel == "Noemi-hadron-tensor-cc") {  // cc stands for current conservation ??
    genie::onebody_currents_sf::compute_hadron_tensor_eigen_cc(mNi, w, wt, pNix, pNiy, pNiz, qtx, qty,
        qtz, f1v, xif2v, fa, fp, hadron_tensor);
  }
  else {	/* What do we want to say if there isn't a model? */
    LOG("UnifiedQELPXSec", pERROR) 
      << "Unknown hadron tensor model: " << fModel
      << ". The options are Noemi-hadron-tensor or Noemi-hadron-tensor-cc";
    exit(1);
  }

}

