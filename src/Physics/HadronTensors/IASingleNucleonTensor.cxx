//____________________________________________________________________________
/*

   Impulse approximation single nucleon response tensor

   Modify the Fortran wrapper (by Noah Steinberg) to C++ native implemetation of one body current
   Native C++ implementation, replacing the previous Fortran wrapper.

   By Liang Liu and ChatGPT and Claude

   Noah Steinberg <nsteinbe \at fnal.gov>
   Fermi National Accelerator Laboratory 


   Liang Liu <liangliu \at fnal.gov>
   Fermi National Accelerator Laboratory 
   Copyright (c) 2003-2026, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org

*/    
//____________________________________________________________________________
#include "Physics/HadronTensors/IASingleNucleonTensor.h"
#include "Physics/HadronTensors/onebody_currents_sf.h"

genie::IASingleNucleonTensor::IASingleNucleonTensor(double w, double mNi, 
    TLorentzVector p4Ni, TLorentzVector p4Nf, 
    const genie::QELFormFactors& FormFactors, 
    std::string model) : fFormFactors(FormFactors) 
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

void genie::IASingleNucleonTensor::initialize_tensor(std::complex<double> (&hadron_tensor)[4][4]) const
{

  // Gather all kinematic and form factor info
  // TO DO: TLorentzVector and FormFactor objects rvalues can't
  // be passed to f90 as they don't have an lvalue

  double mNi  = fmNi;
  double w    = fw;
  double pNix = fp4Ni.X();
  double pNiy = fp4Ni.Y();
  double pNiz = fp4Ni.Z();
  double qtx  = fq4t.X();
  double qty  = fq4t.Y();
  double qtz  = fq4t.Z();
  double wt   = fq4t.E();
  double f1v  = fFormFactors.F1V();
  double xif2v = fFormFactors.xiF2V();
  double fa   = fFormFactors.FA();
  double fp   = fFormFactors.Fp();

  // *Call appropriate fortran subroutine based on name of model provided*
  // translate fortron code to native c++ code
  if (fModel == "Noemi-hadron-tensor") {
    genie::onebody_currents_sf::ComputeNucleonTensor(mNi, w, wt, pNix, pNiy, pNiz, qtx, qty,
        qtz, f1v, xif2v, fa, fp, hadron_tensor);
  }
  else if (fModel == "Noemi-hadron-tensor-cc") {  // cc stands for current conservation ??
    genie::onebody_currents_sf::ComputeNucleonTensorCC(mNi, w, wt, pNix, pNiy, pNiz, qtx, qty,
        qtz, f1v, xif2v, fa, fp, hadron_tensor);
  }
  else {	/* What do we want to say if there isn't a model? */
    LOG("HadronTensors", pFATAL) 
      << "Unknown hadron tensor model: " << fModel
      << ". The options are Noemi-hadron-tensor or Noemi-hadron-tensor-cc";
    exit(1);
  }

}

