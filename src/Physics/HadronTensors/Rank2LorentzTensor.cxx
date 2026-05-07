//____________________________________________________________________________
/*!

\class    genie::Rank2LorentzTensor

\brief    Abstract interface for an object that computes the elements of a rank
          2 Lorentz tensor \f$A^{\mu\nu}\f$.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 17, 2019

\updated  May 06, 2026

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#include "Physics/HadronTensors/Rank2LorentzTensor.h"
#include "Framework/Messenger/Messenger.h"

std::complex<double> genie::Rank2LorentzTensor::Contract(
  const genie::Rank2LorentzTensor& other) const
{
  std::complex<double> sum(0., 0.);
  for (int m = 0; m < 4; ++m) {
    for (int n = 0; n < 4; ++n) {
      genie::TensorIndex_t mu = static_cast<genie::TensorIndex_t>( m );
      genie::TensorIndex_t nu = static_cast<genie::TensorIndex_t>( n );
      sum += this->operator()(mu, nu) * other(mu, nu) * this->Metric(mu, mu)
        * this->Metric(nu, nu);
    }
  }
  return sum;
}

double genie::Rank2LorentzTensor::LeviCivitaProduct(genie::TensorIndex_t mu,
  genie::TensorIndex_t nu, const TLorentzVector& p1, const TLorentzVector& p2)
  const
{
  if (mu == nu) return 0.;

  // This should really be a constexpr std::array once GENIE upgrades to C++11
  // TODO: is there a cleaner way to do this?
  static const genie::TensorIndex_t indices[]
    = { kTIdx_SpaceX, kTIdx_SpaceY, kTIdx_SpaceZ, kTIdx_Time };

  bool set_alpha = false;
  bool set_beta = false;
  genie::TensorIndex_t alpha = genie::kTIdx_Time, beta = genie::kTIdx_Time; // initialize
  for (int k = 0; k < 4; ++k) {
    genie::TensorIndex_t index = indices[k];
    if (mu != index && nu != index) {
      if ( !set_alpha ) {
        alpha = index;
        set_alpha = true;
      }
      else {
        beta = index;
        set_beta = true;
      }
    }
  }

  if ( !(set_alpha && set_beta) ) {
    LOG("HadronTensors", pFATAL) << "alpha and beta not set";
    exit(1);
  }

  // See http://www.cpluscode.org/p/physics-c.html
  // TODO: check that this actually works and that epsilon(0,1,2,3) has the
  // right sign!
  int i = mu + 1;
  int j = nu + 1;
  int k = alpha + 1;
  int l = beta + 1;

  int levi_civita_sign = (i - j) * (i - k) * (i - l) * (j - k) * (j - l)
    * (k - l) / 12;

  double result = levi_civita_sign * (p1[alpha]*p2[beta] - p1[beta]*p2[alpha]);

  return result;
}


// Slight convention change with the fortran wrapper SF model, need different
// convention for the levi civita tensor
double genie::Rank2LorentzTensor::LeviCivitaProductSF(genie::TensorIndex_t mu,
  genie::TensorIndex_t nu, const TLorentzVector& p1, const TLorentzVector& p2)
  const
{
  if (mu == nu) return 0.;

    // Need to use contravariant 4 vectors for this piece of the leptonic tensor
    // TODO: Test
    TLorentzVector p1_con(p1);
    TLorentzVector p2_con(p2);
    
    p1_con.SetVect(-1.*p1.Vect());
    p2_con.SetVect(-1.*p2.Vect());
  
  

  static const genie::TensorIndex_t indices[]
    = { kTIdx_SpaceX, kTIdx_SpaceY, kTIdx_SpaceZ, kTIdx_Time };

  bool set_alpha = false;
  bool set_beta = false;
  genie::TensorIndex_t alpha = genie::kTIdx_Time, beta = genie::kTIdx_Time;  // initialize
  for (int k = 0; k < 4; ++k) {
    genie::TensorIndex_t index = indices[k];
    if (mu != index && nu != index) {
      if ( !set_alpha ) {
        alpha = index;
        set_alpha = true;
      }
      else {
        beta = index;
        set_beta = true;
      }
    }
  }

  if ( !(set_alpha && set_beta) ) {
    LOG("HadronTensors", pFATAL) << "alpha and beta not set";
    exit(1);
  }
  
  int levi_civita_sign = (mu - nu) * (mu - alpha) * (mu - beta) * (nu - alpha) * (nu - beta) * (alpha - beta) / 12;
  double result = levi_civita_sign * (p1_con[alpha]*p2_con[beta] - p1_con[beta]*p2_con[alpha]);

  return result;

}
