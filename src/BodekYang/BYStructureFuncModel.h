//____________________________________________________________________________
/*!

\class    genie::BYStructureFuncModel

\brief    Abstract class.
          Implements part of the DISFormFactorsModelI interface and
          provides some common implementation for concrete Bodek-Yang
          form factor algorithms.

\ref      U.K.Yang and A.Bodek,
          Modeling Deep Inelastic Cross Sections in the Few GeV Region,
          NuINT-01 Proceedings

          U.K.Yang and A.Bodek,
          Parton distributions, d/u and higher twist effect at high x,
          PRL 82, 2467 (1999), hep-ph/9809480,
          PRL 84, 5456 (2000), hep-ph/9912543

          U.K.Yang and A.Bodek,
          Studies of Hugher Twist and Higher Order Effects in NLO and
          NNLO QCD analysis of lepton-nucleon scattering data on F2 and R,
          Eur.Phys.J C13 245, 2000, hep-ex/9908058

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 28, 2004

*/
//____________________________________________________________________________

#ifndef _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_
#define _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_

#include "PartonModel/DISStructureFuncModel.h"
#include "Interaction/Interaction.h"
#include "PDF/PDFModelI.h"

namespace genie {

class BYStructureFuncModel : public DISStructureFuncModel {

public:

  virtual ~BYStructureFuncModel();

  // Overload Algorithm::Configure() to read the config. registry
  // at the algorithm initialization and set private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

  // Overload DISStructureFuncModel::CalculatePDFs() to get the private
  // PDF data member with an attached with a BYPDFModel algorithm
  // computed at x = xw (BY scaling var) and apply the K sea/valence
  // scaling factors.
  bool CalculatePDFs(const Interaction * interaction) const;

protected:

  // protected constructors - abstract class
  BYStructureFuncModel();
  BYStructureFuncModel(const char * param_set);

  // Bodek-Yang model parameters A,B,Csea,Cv1,Cv2
  void   Init(void);
  void   ReadBYParams(void);
  double fA;
  double fB;
  double fCs;
  double fCv1;
  double fCv2;
};

}         // genie namespace

#endif    // _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_
