//____________________________________________________________________________
/*!

\namespace  genie::utils::math

\brief      Simple mathematical utilities not found in ROOT's TMath

\author     Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

\created    May 06, 2004

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org            
*/
//____________________________________________________________________________

#ifndef _MATH_UTILS_H_
#define _MATH_UTILS_H_

#include <vector>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TLorentzVector.h>
#include "Framework/Utils/Range1.h"
#include "cmath"

using std::vector;

namespace genie {
namespace utils {

namespace math
{

  // This class has been created to perform several operations with long 
  // doubles. It is needed in HEDIS because the kinematics of the outgoing
  // particles can be so large that the on-shell feature is not fulfilled 
  // many times due to the precission of double. 
  class LongLorentzVector {

    public :
      LongLorentzVector(double px, double py, double pz, double e) {
        fPx = (long double) px;
        fPy = (long double) py;
        fPz = (long double) pz;
        fE  = (long double) e;
      }
      LongLorentzVector(const TLorentzVector & p4) { 
        fPx = (long double) p4.Px();  
        fPy = (long double) p4.Py();  
        fPz = (long double) p4.Pz();  
        fE  = (long double) p4.E();  
      }
     ~LongLorentzVector() {}

      long double Px (void) { return fPx; }
      long double Py (void) { return fPy; }
      long double Pz (void) { return fPz; }
      long double E  (void) { return fE;  }
      long double P  (void) { return sqrtl(fPx*fPx+fPy*fPy+fPz*fPz);     }
      long double M  (void) { return sqrtl(fE*fE-fPx*fPx-fPy*fPy-fPz*fPz); }
      long double M2 (void) { return fE*fE-fPx*fPx-fPy*fPy-fPz*fPz; }
      long double Dx (void) { return fPx/sqrtl(fPx*fPx+fPy*fPy+fPz*fPz); }
      long double Dy (void) { return fPy/sqrtl(fPx*fPx+fPy*fPy+fPz*fPz); }
      long double Dz (void) { return fPz/sqrtl(fPx*fPx+fPy*fPy+fPz*fPz); }

      void Rotate    (LongLorentzVector axis) {
        long double up = axis.Dx()*axis.Dx() + axis.Dy()*axis.Dy();
        if (up) {
          up = sqrtl(up);
          long double pxaux = fPx,  pyaux = fPy,  pzaux = fPz;
          fPx = (axis.Dx()*axis.Dz()*pxaux - axis.Dy()*pyaux + axis.Dx()*up*pzaux)/up;
          fPy = (axis.Dy()*axis.Dz()*pxaux + axis.Dx()*pyaux + axis.Dy()*up*pzaux)/up;
          fPz = (axis.Dz()*axis.Dz()*pxaux -           pxaux + axis.Dz()*up*pzaux)/up;
        } 
        else if (axis.Dz() < 0.) { // phi=0  teta=pi
          fPx = -fPx; 
          fPz = -fPz; 
        }
      }    

    void BoostZ    (long double bz) {
      long double b2 = bz*bz;
      long double gamma = 1.0 / sqrtl(1.0 - b2);
      long double bp = bz*fPz;
      long double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;
      fPz = fPz + gamma2*bp*bz + gamma*bz*fE;
      fE  = gamma*(fE + bp);    
    }    

    void BoostY    (long double by) {
      long double b2 = by*by;
      long double gamma = 1.0 / sqrtl(1.0 - b2);
      long double bp = by*fPy;
      long double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;
      fPy = fPy + gamma2*bp*by + gamma*by*fE;
      fE  = gamma*(fE + bp);    
    }
    
    private :

      long double fPx;
      long double fPy;
      long double fPz;
      long double fE;
  };

  // Cholesky decomposition. Returns lower triangular matrix.
  TMatrixD CholeskyDecomposition (const TMatrixD& cov);
  // Generates a vector of correlated parameters.
  TVectorD CholeskyGenerateCorrelatedParams (const TMatrixD& Lch, TVectorD& mean);
  TVectorD CholeskyGenerateCorrelatedParams (const TMatrixD& Lch, TVectorD& mean, TVectorD& g_uncorrelated);
  // Generates a vector of correlated parameter variations.
  TVectorD CholeskyGenerateCorrelatedParamVariations (const TMatrixD& Lch);
  TVectorD CholeskyCalculateCorrelatedParamVariations(const TMatrixD& Lch, TVectorD& g_uncorrelated);

  double KahanSummation (double x[], unsigned int n);
  double KahanSummation (const vector<double> & x);

  bool   AreEqual       (double x1, double x2);
  bool   AreEqual       (float  x1, float  x2);

  bool   IsWithinLimits (double x, Range1D_t range);
  bool   IsWithinLimits (float  x, Range1F_t range);
  bool   IsWithinLimits (int    i, Range1I_t range);

  double NonNegative    (double x);
  double NonNegative    (float  x);

} // math  namespace
} // utils namespace
} // genie namespace

#endif // _MATH_UTILS_H_
