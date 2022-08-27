//____________________________________________________________________________
/*!

\class    genie::DCCSPPPXSec

\brief    Class calculate differental cross-section d^4(sig)/d(Q2)d(W)d(cost)d(phi),
          for specific W, Q2, neutrino energy(in lab frame) & pion angles in the Adler frame, where \n
          Q2          : Sqaured 4-momentum transfer, Q2 = -k*k         \n                        
          W           : Invariant mass                                 \n                      
          cost        : Cosine of pion polar angle in N\pi rest frame  \n                      
          phi         : Pion azimuthal angle in N\pi rest frame        \n

for the following channels:      

          1       l + p -> l + p + pi0
          2       l + p -> l + n + pi+
          3       l + n -> l + n + pi0
          4       l + n -> l + p + pi-
                                                                                           
                                                          
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
//#include "Framework/EventGen/XSecAlgorithmI.h"
//#include "Framework/ParticleData/BaryonResonance.h"

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
      
      void LoadConfig (void);
      VAmpl Amplitudes(double W, double Q2, unsigned int L, SppChannel_t spp_chn) const;
      CPtrDT GetDataTable(SppChannel_t spp_chn) const;
      CPtrDT BuildDataTable(SppChannel_t spp_chn) const;
      UPtrDT ParseDataTableFile( std::string full_file_name ) const;
      void GetTablePos(double W, double Q2, unsigned int L, TablePos & tabpos) const;
      std::string FindDataTableFile(const std::string &basename, bool &ok) const;
      std::string GetDataTableFileBasename(SppChannel_t spp_chn) const;
      double dPdx (unsigned int L, double x) const;
      double d2Pdx2 (unsigned int L, double x) const;
      
      // configuration data
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

      const XSecIntegratorI * fXSecIntegrator;
                 
  };
  
  
}       // genie namespace

#endif  // _DCC_SPP_PXSEC_H_
