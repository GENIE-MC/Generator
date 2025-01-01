//____________________________________________________________________________
/*!

\namespace  genie::utils

\brief      Common functions used for handling generation of the primary
            lepton, regardless of whether the relevant class inherits from
            PrimaryLeptonGenerator or not.

\author     Steven Gardiner <gardiner \at fnal.gov>
            Fermi National Accelerator Laboratory

\created    May 01, 2020

\cpright    Copyright (c) 2003-2023, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PRIMARY_LEPTON_UTILS_H
#define _PRIMARY_LEPTON_UTILS_H

#include <vector>
#include <complex>
#include <stdexcept>

#include "TVector3.h"

namespace genie {

class GHepRecord;

namespace utils {
  class HermitianMatrix;

  void SetPrimaryLeptonPolarization( GHepRecord* ev );
  
  void CalculatePolarizationVectorWithNuclearTensor(
                        TVector3 & polarization,
                        const TLorentzVector & neutrinoMom,
                        const TLorentzVector & leptonMom, 
                        bool isLeftPolarized,
                        const HermitianMatrix & NTensor);
  
  void CalculatePolarizationVectorWithStructureFunctions(
                        TVector3 & polarization,
                        const TLorentzVector & neutrinoMom,
                        const TLorentzVector & leptonMom, 
                        const TLorentzVector & inNucleonMom,
                        const TLorentzVector & q4,
                        bool isLeftPolarized,
                        double M,
                        double W1,
                        double W2,
                        double W3,
                        double W4,
                        double W5,
                        double W6);
                        
  void CalculatePolarizationVectorInTargetRestFrame(
                      TVector3 & polarization,
                      const TLorentzVector & neutrinoMomTRF,
                      const TLorentzVector & leptonMomTRF, 
                      bool isLeftPolarized,
                      double M,
                      double W1,
                      double W2,
                      double W3,
                      double W4,
                      double W5,
                      double W6);
                        

  inline int g(int a, int b) ///< metric g^{ab}=g_{ab}=diag(1,-1,-1,1)
  {
      return (a==b)*(2*(a==0) - 1);
  }
  inline int e(int a, int b, int c, int d) ///< Levi-Chevita symbol, where e_{0123}=+1
  {
      return (b - a)*(c - a)*(d - a)*(c - b)*(d - b)*(d - c)/12;
  }

  
  class HermitianMatrix 
  {
    public:
        // Constructor
        explicit HermitianMatrix(size_t size) : n(size), data(size * (size + 1) / 2) {}
    
        // Set an element
        void set(size_t i, size_t j, const std::complex<double>& value) 
        {
            if (i >= n || j >= n) 
            {
                throw std::out_of_range("Indices out of range.");
            }
            if (i == j) 
            {
                data[index(i, j)] = std::real(value);
            }
            data[index(i, j)] = (i < j) ? value : std::conj(value);
        }
    
        // Get an element
        std::complex<double> get(size_t i, size_t j) const 
        {
            if (i >= n || j >= n) 
            {
                throw std::out_of_range("Indices out of range.");
            }
            return (i <= j) ? data[index(i, j)] : std::conj(data[index(j, i)]);
        }
    
        // Operator () for accessing elements
        std::complex<double> operator()(size_t i, size_t j) const 
        {
            return get(i, j);
        }
        
    private:
        size_t n; // Matrix size
        std::vector<std::complex<double>> data; // Vector for storing elements of the upper triangular part
    
        // Convert (i, j) indices to a single index in the 1D array
        size_t index(size_t i, size_t j) const 
        {
            if (i > j) std::swap(i, j); // Use only the upper triangular part
            return i * n - i * (i + 1) / 2 + j;
        }
  };
  
  

} // utils   namespace
} // genie   namespace

#endif // _PRIMARY_LEPTON_UTILS_H
