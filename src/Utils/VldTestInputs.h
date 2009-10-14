//__________________________________________________________________________
/*!

\class    VldTestInputs

\brief    Hold typical inputs for validation programs

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Oct 12, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//__________________________________________________________________________

#ifndef _VLD_TEST_INPUTS_H_
#define _VLD_TEST_INPUTS_H_

#include <TChain.h>
#include <TTree.h>
#include <TFile.h>

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace genie {

class VldTestInputs
{
public:
  VldTestInputs(const int nmaxmodels=10);
 ~VldTestInputs(void);

  int     NModels  (void)       const;
  string  ModelTag (int imodel) const;
  TFile*  XSecFile (int imodel) const;
  TChain* EvtChain (int imodel) const;

  bool LoadFromFile(string xmlfile);

private:

  void Init(const int nmaxmodels);


  int             fNModels;
  vector<string>  *fModelTag;
  vector<TFile*>  *fXSecFile;
  vector<TChain*> *fEvtChain;
};

}      // genie namespace

#endif // _VLD_TEST_INPUTS_H_

