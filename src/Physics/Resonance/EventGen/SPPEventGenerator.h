//____________________________________________________________________________
/*!

\class    genie::SPPEventGenerator

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
          Vadim Naumov <vnaumov@theor.jinr.ru>,  Joint Institute for Nuclear Research \n

\created  May 9, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SPP_EVENT_GENERATOR_H_
#define _SPP_EVENT_GENERATOR_H_

#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Framework/Utils/Range1.h"
#include "Physics/Common/KineGeneratorWithCache.h"


namespace genie {

class SPPEventGenerator : public KineGeneratorWithCache {

public :
  SPPEventGenerator();
  SPPEventGenerator(string config);
 ~SPPEventGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig      (void);
  double ComputeMaxXSec  (const Interaction * interaction) const;
  int GetRecoilNucleonPdgCode(Interaction * interaction) const;
  int GetFinalPionPdgCode(Interaction * interaction) const;
  double fWcut;
  struct Vertex
  {
    Vertex () : Vertex (0., 0., 0., 0.)
    {};
    Vertex (double px1, double px2, double px3, double px4) : x1(px1), x2(px2), x3(px3), x4(px4)
    {};
    ~Vertex(){};
    double x1, x2, x3, x4;
    void Print (std::ostream& os)
    {
      os << "(" << x1 << "," << x2 << "," << x3 << "," << x4 << ")";
    };
    bool operator == (const Vertex &v) const
    {
       double epsilon = 1e-5;
       return (TMath::Abs(this->x1 - v.x1) < epsilon || TMath::Abs(this->x2 - v.x2) < epsilon || TMath::Abs(this->x3 - v.x3) < epsilon || TMath::Abs(this->x4 - v.x4) < epsilon);
    };
   
  };
  
  struct Cell
  {
     Cell(){};
     ~Cell(){};
     Vertex Vertex1;
     Vertex Vertex2;
     void Print (std::ostream& os)
     {
        os << std::endl;
        os << "vertex1 = ";
        Vertex1.Print(os);
        os << std::endl;
        os << "vertex2 = ";
        Vertex2.Print(os);
        os << std::endl;
     };
  };

  int fMaxDepth;  ///< Maximum depth of dividing parent cell
  
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
class d4XSecMK_dWQ2CosThetaPhi_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d4XSecMK_dWQ2CosThetaPhi_E(const XSecAlgorithmI * m, const Interaction * i, double wcut);
 ~d4XSecMK_dWQ2CosThetaPhi_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  Interaction *    fInteraction;
  Range1D_t Wl;
  double fWcut;
  bool isZero;
  KPhaseSpace * kps;
};


} // gsl   namespace
} // utils namespace

}      // genie namespace
#endif // _SPP_EVENT_GENERATOR_H_
