//____________________________________________________________________________
/*!

\class    genie::MKSPPPXSec

\brief    Class calculate differental cross-section d^4(sig)/d(Q2)d(W)d(cost)d(phi),
          for specific W, Q2, neutrino energy(in lab frame) & pion angles in the Adler frame, where \n
          Q2          : Sqaured 4-momentum transfer, Q2 = -k*k         \n                        
          W           : Invariant mass                                 \n                      
          cost        : Cosine of pion polar angle in Adler frame      \n                      
          phi         : Pion azimuthal angle in Adler frame            \n

for the following channels:      

          1       nu + p -> l      + p + pi+
          2       nu + n -> l      + p + pi0
          3       nu + n -> l      + n + pi+
          4   antinu + n -> l+     + n + pi-
          5   antinu + p -> l+     + n + pi0
          6   antinu + p -> l+     + p + pi-
          7       nu + p -> nu     + p + pi0
          8       nu + p -> nu     + n + pi+
          9       nu + n -> nu     + n + pi0
          10      nu + n -> nu     + p + pi-
          11  antinu + p -> antinu + p + pi0
          12  antinu + p -> antinu + n + pi+
          13  antinu + n -> antinu + n + pi0
          14  antinu + n -> antinu + p + pi-
                                                                                           
                                                          
\ref      
          1.  Kabirnezhad M., Phys.Rev.D 97 (2018) 013002 (Erratim: arXiv:1711.02403)
          2.  Kabirnezhad M., Ph. D. Thesis (https://www.ncbj.gov.pl/sites/default/files/m.kabirnezhad_thesis_0.pdf , 
                                        https://www.ncbj.gov.pl/sites/default/files/m.kabirnezhad-thesis_final_0.pdf) 
          3.  Rein D., Sehgal L., Ann. of Phys. 133 (1981) 79-153
          4.  Rein D., Z.Phys.C 35 (1987) 43-64
          5.  Adler S.L., Ann. Phys. 50 (1968) 189
          6.  Graczyk K., Sobczyk J., Phys.Rev.D 77 (2008) 053001 [Erratum: ibid.D 79 (2009) 079903]
          7.  Jacob M., Wick G.C., Ann. of Phys. 7 (1959) 404-428
          8.  Hernandez E., Nieves J., Valverde M., Phys.Rev.D 76 (2007) 033005
          9.  Adler S.L., Nussinov S., Paschos E.A., Phys. Rev. D 9 (1974) 2125-2143 [Erratum: ibid D 10 (1974) 1669]
          10. Paschos E.A., Yu J.Y., Sakuda M., Phys. Rev. D 69 (2004) 014013 [arXiv: hep-ph/0308130] 
          11. Yu J.Y., "Neutrino interactions and  nuclear  effects in oscillation experiments and the 
             nonperturbative dispersive  sector in strong (quasi-)abelian  fields", Ph. D.Thesis, Dortmund U., Dortmund, 2002 (unpublished)
          12. Kakorin I., Kuzmin K., Naumov V. "Report on implementation of the MK-model for resonance single-pion production into GENIE"
                                               (https://genie-docdb.pp.rl.ac.uk/cgi-bin/private/ShowDocument?docid=181, 
                                                http://theor.jinr.ru/NeutrinoOscillations/Papers/Report_MK_implementation_GENIE.pdf)              

\authors  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>,  Joint Institute for Nuclear Research \n
          Vadim Naumov <vnaumov@theor.jinr.ru>,  Joint Institute for Nuclear Research \n
          adapted from code provided by \n
          Minoo Kabirnezhad <minoo.kabirnezhad@physics.ox.ac.uk>
          University of Oxford, Department of Physics \n
          based on code of \n
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Nov 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MK_SPP_PXSEC_2014_H_
#define _MK_SPP_PXSEC_2014_H_

#include <vector>
#include <complex>
#include <functional>
#include <algorithm>
#include <type_traits>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResList.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"
#include "Physics/Resonance/XSection/FKR.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"

namespace genie {

  
  class XSecIntegratorI;

  class MKSPPPXSec: public XSecAlgorithmI {

    
    public:
      MKSPPPXSec();
      MKSPPPXSec(string config);
      virtual ~MKSPPPXSec();

      // implement the XSecAlgorithmI interface 
      double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral     (const Interaction * i) const;
      bool   ValidProcess (const Interaction * i) const;
      bool   ValidKinematics(const Interaction * interaction) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:
      
      // Helicities \~{F}^{\lambda_k}_{\lambda_2 \lambda_1} or \~{G}^{\lambda_k}_{\lambda_2 \lambda_1} 
      // Definition are given in eq. 16 from ref. 1
      // The structure:
      // Number of resonance
      // F(Vector) or G(Axial)
      // \lambda_k (boson polarization eL(M), eR(P), e-(OM), e+(OP))
      // \lambda_2 (final nucleon polarization -1/2(M), +1/2(P)
      // \lambda_1 (initial nucleon polarization -1/2(M), +1/2(P)
      enum Current { VECTOR, AXIAL };
      enum BosonPolarization   { LEFT, RIGHT, MINUS0, PLUS0 };
      enum NucleonPolarization { MINUS, PLUS };
      
      //#ifndef DOXYGEN_SHOULD_SKIP_THIS
      template < typename C, C beginVal, C endVal>
      class Iterator {
        typedef typename std::underlying_type<C>::type val_t;
        int val;
      public:
        Iterator(const C & f) : val(static_cast<val_t>(f)) {}
        Iterator() : val(static_cast<val_t>(beginVal)) {}
        Iterator operator++() {
          ++val;
          return *this;
        }
        C operator*() { return static_cast<C>(val); }
        Iterator begin() { return *this; } //default ctor is good
        Iterator end() {
            static const Iterator endIter=++Iterator(endVal); // cache it
            return endIter;
        }
        bool operator!=(const Iterator& i) { return val != i.val; }
      };
      
      using  CurrentIterator             = Iterator<Current,             Current::VECTOR,            Current::AXIAL>            ;
      using  BosonPolarizationIterator   = Iterator<BosonPolarization,   BosonPolarization::LEFT,    BosonPolarization::PLUS0>   ;
      using  NucleonPolarizationIterator = Iterator<NucleonPolarization, NucleonPolarization::MINUS, NucleonPolarization::PLUS> ;
       
      template<typename T> 
      struct is_complex : std::false_type {};

      template<typename T> 
      struct is_complex<std::complex<T>> : std::true_type {};

      template<bool C, typename T = void>
      using enable_if_t = typename std::enable_if<C, T>::type;
      
      template <typename T>
      class HelicityBkgAmp {
              
        public:

          HelicityBkgAmp() : array(32) {}
          HelicityBkgAmp(const HelicityBkgAmp& ha)
          {
            array = ha.array;
          }
          ~HelicityBkgAmp(){}
          T& operator() (Current hatype, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1)
          {
            int indx = 2*(2*(4*hatype+lambda_k)+lambda_2)+lambda_1;
            return array[indx];
          }
          template<typename S = T, enable_if_t<is_complex<S>{}>* = nullptr>
          auto Re(Current hatype, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> typename S::value_type
          {
            return this->operator()(hatype, lambda_k, lambda_2, lambda_1).real();
          }
          template<typename S = T, enable_if_t<is_complex<S>{}>* = nullptr>
          auto Im(Current hatype, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> typename S::value_type
          {
            return this->operator()(hatype, lambda_k, lambda_2, lambda_1).imag();
          }
          template<typename S = T, enable_if_t<!is_complex<S>{}>* = nullptr>
          auto Re(Current hatype, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> S
          {
            return this->operator()(hatype, lambda_k, lambda_2, lambda_1);
          }
          template<typename S = T, enable_if_t<!is_complex<S>{}>* = nullptr>
          auto Im(Current hatype, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> S
          {
            return 0;
          }
          HelicityBkgAmp& operator*= (double factor)
          {
            std::transform(array.begin(), array.end(), array.begin(), std::bind(std::multiplies<T>(), std::placeholders::_1, factor));
            return *this;
          }
          
          HelicityBkgAmp& operator/= (double factor)
          {
            std::transform(array.begin(), array.end(), array.begin(), std::bind(std::multiplies<T>(), std::placeholders::_1, 1./factor));
            return *this;
          }
    
          HelicityBkgAmp& operator+= (const HelicityBkgAmp& ha)
          {
            std::transform(ha.array.begin(), ha.array.end(), array.begin(), array.begin(), std::bind(std::plus<T>(), std::placeholders::_1, std::placeholders::_2));
            return *this;
          }
          
          HelicityBkgAmp& operator-= (const HelicityBkgAmp& ha)
          {
            std::transform(ha.array.begin(), ha.array.end(), array.begin(), array.begin(), std::bind(std::minus<T>(), std::placeholders::_2, std::placeholders::_1));
            return *this;
          }
    
          HelicityBkgAmp& operator= (const HelicityBkgAmp& ha)
          {
            if (this != &ha)
              array = ha.array;
            return *this;
          }
          
          friend HelicityBkgAmp operator+(HelicityBkgAmp lhs, const HelicityBkgAmp& rhs)
          {
             lhs += rhs;
             return lhs;
          }
          
          friend HelicityBkgAmp operator-(HelicityBkgAmp lhs, const HelicityBkgAmp& rhs)
          {
             lhs -= rhs;
             return lhs;
          }
          
          friend HelicityBkgAmp operator*(HelicityBkgAmp ha, double factor)
          {
             ha *= factor;
             return ha;
          }
          
          friend HelicityBkgAmp operator*(double factor, HelicityBkgAmp ha)
          {
             ha *= factor;
             return ha;
          }
          
          friend HelicityBkgAmp operator/(HelicityBkgAmp ha, double factor)
          {
             ha /= factor;
             return ha;
          }
          
          friend HelicityBkgAmp operator/(double factor, HelicityBkgAmp ha)
          {
             ha /= factor;
             return ha;
          }
          
        private:
          std::vector<T> array; //2*4*2*2; 
      
      };
      
      template <typename T>
      class HelicityAmpVminusARes {
      
        public:

          HelicityAmpVminusARes() : array(1)
          {
            array.reserve(288);
          }
          ~HelicityAmpVminusARes(){}
          T& operator() (Resonance_t res, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1)
          {
            if (res == kNoResonance)
            {
              // meaningless to return anything
              gAbortingInErr = true;
              LOG("MKSPPPXSec", pFATAL) << "Unknown resonance " << res;
              exit(1);
            }
            int indx = 2*(2*(4*res+lambda_k)+lambda_2)+lambda_1;
            return array[indx];
          }
          
          template<typename S = T, enable_if_t<is_complex<S>{}>* = nullptr>
          auto Re(Resonance_t res, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> typename S::value_type
          {
            return this->operator()(res, lambda_k, lambda_2, lambda_1).real();
          }
          template<typename S = T, enable_if_t<is_complex<S>{}>* = nullptr>
          auto Im(Resonance_t res, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> typename S::value_type
          {
            return this->operator()(res, lambda_k, lambda_2, lambda_1).imag();
          }
          template<typename S = T, enable_if_t<!is_complex<S>{}>* = nullptr>
          auto Re(Resonance_t res, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> S
          {
            return this->operator()(res, lambda_k, lambda_2, lambda_1);
          }
          template<typename S = T, enable_if_t<!is_complex<S>{}>* = nullptr>
          auto Im(Resonance_t res, BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> S
          {
            return 0;
          }
          
        private:
          std::vector<T> array; //nres*4*2*2, nres=18
      
      };
      
      template <typename T>
      class SumHelicityAmpVminusARes {
      
        public:

          SumHelicityAmpVminusARes() : array(16){}
          ~SumHelicityAmpVminusARes(){}
          T& operator() (BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1)
          {
            int indx = 2*(2*lambda_k+lambda_2)+lambda_1;
            return array[indx];
          }
          
          template<typename S = T, enable_if_t<is_complex<S>{}>* = nullptr>
          auto Re(BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> typename S::value_type
          {
            return this->operator()(lambda_k, lambda_2, lambda_1).real();
          }
          template<typename S = T, enable_if_t<is_complex<S>{}>* = nullptr>
          auto Im(BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> typename S::value_type
          {
            return this->operator()(lambda_k, lambda_2, lambda_1).imag();
          }
          template<typename S = T, enable_if_t<!is_complex<S>{}>* = nullptr>
          auto Re(BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> S
          {
            return this->operator()(lambda_k, lambda_2, lambda_1);
          }
          template<typename S = T, enable_if_t<!is_complex<S>{}>* = nullptr>
          auto Im(BosonPolarization lambda_k, NucleonPolarization lambda_2, NucleonPolarization lambda_1) -> S
          {
            return 0;
          }
          
        private:
          std::vector<T> array; //4*2*2
      
      };
        
      //#endif
      int Lambda (NucleonPolarization l) const;
      int Lambda (BosonPolarization l) const;
      int PhaseFactor(BosonPolarization lk, NucleonPolarization l1, NucleonPolarization l2) const;
      
      void LoadConfig (void);
      
/*
      // Not used n the latest version 
      mutable FKR_MK fFKR;
      
      const MKHelicityAmplModelI * fHAmplModelCC;
      const MKHelicityAmplModelI * fHAmplModelNCp;
      const MKHelicityAmplModelI * fHAmplModelNCn;
*/      
      
      mutable FKR fFKR;
      const RSHelicityAmplModelI * fHAmplModelCC;
      const RSHelicityAmplModelI * fHAmplModelNCp;
      const RSHelicityAmplModelI * fHAmplModelNCn;


      // configuration data
      double   fFermiConstant ;
      double   fCA50;              ///< FKR parameter Zeta
      double   fOmega;             ///< FKR parameter Omega
      double   fMa2;               ///< (axial mass)^2
      double   fMv2;               ///< (vector mass)^2
      double   fSin2Wein;          ///< sin^2(Weingberg angle)
      double   fVud;               ///< |Vud| (magnitude ud-element of CKM-matrix)
      double   fXSecScaleCC;       ///< external CC xsec scaling factor
      double   fXSecScaleNC;       ///< external NC xsec scaling factor
      const ELFormFactorsModelI  * fElFFModel;          ///< Elastic form facors model for background contribution
      const QELFormFactorsModelI * fFormFactorsModel;   ///< Quasielastic form facors model for background contribution
      const QELFormFactorsModelI * fEMFormFactorsModel; ///< Electromagnetic form factors model for background contribution
      
      string fKFTable;             ///< table of Fermi momentum (kF) constants for various nuclei
      bool fUseRFGParametrization; ///< use parametrization for fermi momentum insted of table?
      bool fUsePauliBlocking;      ///< account for Pauli blocking?

      mutable QELFormFactors  fFormFactors;      ///<  Quasielastic form facors for background contribution
      mutable QELFormFactors  fEMFormFactors;    ///<  Electromagnetic form facors for background contribution
      double  f_pi;                              ///<  Constant for pion-nucleon interaction
      double  FA0;                               ///<  Axial coupling (value of axial form factor at Q2=0)
      double  Frho0;                             ///<  Value of form fator F_rho at t=0 
      /// Parameters for vector virtual form factor
      /// for background contribution, which equal to:  
      /// 1,                                              W<VWmin   
      /// V3*W^3+V2*W^2+V1*W+V0                     VWmin<W<VWmax
      /// 0                                               W>VWmax
      double fBkgVWmin;
      double fBkgVWmax; 
      double fBkgV3;  
      double fBkgV2;   
      double fBkgV1;   
      double fBkgV0;   
      /*
      // Not used n the latest version 
      /// Parameters for axial-vector virtual form factor
      /// for background contribution, which equal to:  
      /// 1,                                              W<VWmin   
      /// A3*W^3+A2*W^2+A1*W+A0                     VWmin<W<VWmax
      /// 0                                               W>VWmax
      double fBkgAWmin;
      double fBkgAWmax; 
      double fBkgA3;  
      double fBkgA2;  
      double fBkgA1;  
      double fBkgA0; 
      */
      double fRho770Mass;                        ///< Mass of rho(770) meson 
      
      const XSecIntegratorI * fXSecIntegrator;
      
      BaryonResList  fResList;
                 
  };
  
  
}       // genie namespace

#endif  // _MK_SPP_PXSEC_2014_H_
