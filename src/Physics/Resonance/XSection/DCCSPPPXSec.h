//____________________________________________________________________________
/*!

\class    genie::DCCSPPPXSec

\brief    
Class calculate differental cross-sections  
\f[
\frac{d^4\sigma}{dQ^2dWd\cos\theta_\pi d\phi_\pi}
\f]
or
\f[
\frac{d^3\sigma}{dQ^2dWd\cos\theta_\pi}
\f]
for specific neutrino energy (in lab frame), where:


Variable             | Description
---------------------|-----------------------------------------------------
\f$W\f$              | Invariant mass
\f$Q^2\f$            | Sqaured 4-momentum transfer
\f$\cos\theta_\pi\f$ | Cosine of pion polar angle in \f$\pi\f$N rest frame
\f$\phi_\pi\f$       | Pion azimuthal angle in \f$\pi\f$N rest frame
for the following channels:
-#  \f$\ell + p \to \ell + p + \pi^0\f$
-#  \f$\ell + p \to \ell + n + \pi^+\f$
-#  \f$\ell + n \to \ell + n + \pi^0\f$
-#  \f$\ell + n \to \ell + p + \pi^-\f$
                                                                                      
\ref      
          1. T. Sato , T.-S. H. Lee, Phys. Rev. C 54, 2660(1996)
          2. A. Matsuyama , T.-S. H. Lee, T. Sato, Phys. Rept. 439, 193(2007)
          3. https://www.phy.anl.gov/theory/research/anl-osaka-pwa/ (and other references in it)
          4. https://www.phy.anl.gov/theory/research/anl-osaka-pwa/crst-eepi.pdf
    

\authors  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n

\created  Sep 06, 2022

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DCC_SPP_PXSEC_H_
#define _DCC_SPP_PXSEC_H_

#include <vector>
#include <complex>
#include <map>
#include <memory>
#include <string>

#include "Framework/Interaction/SppChannel.h"
#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

  
  class XSecIntegratorI;

  class DCCSPPPXSec: public XSecAlgorithmI {

    
    public:
      DCCSPPPXSec();
      DCCSPPPXSec(std::string config);
      virtual ~DCCSPPPXSec();

      // implement the XSecAlgorithmI interface 
      double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral     (const Interaction * i) const;
      bool   ValidProcess (const Interaction * i) const;
      bool   ValidKinematics(const Interaction * interaction) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(std::string config);

    private:
            
      using VAmpl = std::vector < std::complex<double> >;
      using CPtrDT  = const std::vector<std::vector<double> > *;
      using UPtrDT = std::unique_ptr < std::vector<std::vector<double> > >;
      
       
      struct TablePos
      {
        unsigned int lo_row_W; 
        unsigned int hi_row_W; 
        unsigned int lo_row_Q2; 
        unsigned int hi_row_Q2; 
        double lo_W;
        double hi_W; 
        double lo_Q2; 
        double hi_Q2;
      };
      
      void LoadConfig (void);
      VAmpl Amplitudes(double W, double Q2, unsigned int L, SppChannel_t spp_chn) const;
      CPtrDT GetDataTable(SppChannel_t spp_chn) const;
      CPtrDT BuildDataTable(SppChannel_t spp_chn) const;
      UPtrDT ParseDataTableFile( std::string full_file_name ) const;
      void GetTablePos(double W, double Q2, unsigned int L, TablePos & tabpos) const;
      std::string FindDataTableFile(const std::string &basename, bool &ok) const;
      std::string GetDataTableFileBasename(SppChannel_t spp_chn) const;
      double dPdx (int L, double x) const;
      double d2Pdx2 (int L, double x) const;
            
      
      // configuration data
      /// Cache of tables with amplitudes that have been fully loaded into memory
      ///
      /// Keys are SPP channels IDs, values are pointers to table with amplitudes
      mutable std::map<SppChannel_t, UPtrDT> fTables;
      /// If true, logging messages will be issued when a requested hadron tensor
      /// file cannot be found
      bool fWarnIfMissing;
      /// Paths to check when searching for hadron tensor data files
      std::vector<std::string> fDataPaths;
      /// Table of Fermi momentum (kF) constants for various nuclei
      std::string fKFTable;
      /// Use parametrization for fermi momentum insted of table?
      bool fUseRFGParametrization;
      /// Account for Pauli blocking?
      bool fUsePauliBlocking;
      /// External EM xsec scaling factor
      double fXSecScaleEM;
      const XSecIntegratorI * fXSecIntegrator;
                 
  };
  
  
}       // genie namespace

#endif  // _DCC_SPP_PXSEC_H_
