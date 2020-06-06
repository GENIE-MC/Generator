//____________________________________________________________________________
/*!

\class    genie::RSPPEventGenerator

\brief    Generates resonance single pion production event for the following channels:      

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
          
Is a concrete implementation of the EventRecordVisitorI interface.

\authors  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>,  Joint Institute for Nuclear Research \n
          Vadim Naumov <vnaumov@theor.jinr.ru>,  Joint Institute for Nuclear Research \n

\created  May 9, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RSPP_EVENT_GENERATOR_H_
#define _RSPP_EVENT_GENERATOR_H_

#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>

#include "Framework/Utils/Range1.h"
#include "Physics/Common/KineGeneratorWithCache.h"


namespace genie {

class RSPPEventGenerator : public EventRecordVisitorI {

public :
  RSPPEventGenerator();
  RSPPEventGenerator(string config);
 ~RSPPEventGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig      (void);
  int GetRecoilNucleonPdgCode(Interaction * interaction) const;
  int GetFinalPionPdgCode(Interaction * interaction) const;
  string CreateKey(int num_grid, Interaction * interaction) const;
  void OpenFile (void);
  
  mutable TFoam * fFoam = 0;        ///< MC adaptive simulator
  mutable TFile * fFoam_file;       ///< File to store MC Foam grid
  mutable const XSecAlgorithmI * fXSecModel;
  bool   fGenerateUniformly;    ///< Uniform over allowed phase space + event weight?
  double fNextGriddE;           ///< Interval between reference energies of neutrino (in hit nucleon rest frame) at which the foam grid is built
  int    fFOAMNCells;           ///< No of allocated number of cells  
  int    fFOAMNSampl;           ///< No. of MC events in the cell MC exploration        
  int    fFOAMNBin;             ///< No. of bins in edge-histogram in cell exploration        
  int    fFOAMOptRej;           ///< OptRej = 0, weighted; OptRej=1, wt=1 MC events        
  int    fFOAMOptDrive;         ///< Maximum weight reduction, =1 for variance reduction         
  int    fFOAMEvPerBin;         ///< Maximum number of the effective wt=1 events/bin, EvPerBin=0 deactivates this option
  int    fFOAMChat;             ///< =0,1,2 is the `‘chat level’' in the standard output      
  double fFOAMMaxWtRej;         ///< Maximum weight used to get w=1 MC events       
  
  
  
  
};

class XSecAlgorithmI;
class Interaction;

namespace utils {
namespace gsl   {



//.....................................................................................
//
// genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E
// A 4-D cross section function: d4XSecMK_dWQ2CosThetaPhi_E = f(W, Q2, CosTheta, Phi)|(fixed E)
//
class d4XSecMK_dWQ2CosThetaPhi_E: public ROOT::Math::IBaseFunctionMultiDim, public TFoamIntegrand
{
public:
  d4XSecMK_dWQ2CosThetaPhi_E(const XSecAlgorithmI * m, const Interaction * i);
 ~d4XSecMK_dWQ2CosThetaPhi_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
  
  // TFoamIntegrand interface
  double Density (int ndim, double * xin);

private:
  const XSecAlgorithmI * fModel;
  Interaction *    fInteraction;
  Range1D_t Wl;
  bool isZero;
  KPhaseSpace * kps;
};


} // gsl   namespace
} // utils namespace

}      // genie namespace
#endif // _RSPP_EVENT_GENERATOR_H_
