//__________________________________________________________________________
/*!

\class    VldTestInputs

\brief    Hold typical inputs for validation programs

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Oct 12, 2009

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//__________________________________________________________________________

#ifndef _VLD_TEST_INPUTS_H_
#define _VLD_TEST_INPUTS_H_

#include <iostream>
#include <string>
#include <vector>

#include <TChain.h>
#include <TTree.h>
#include <TFile.h>

using std::ostream;
using std::string;
using std::vector;

namespace genie {
namespace utils {
namespace vld   {

class VldTestInputs
{
public:
  VldTestInputs(bool chain=true, const int nmaxmodels=10);
 ~VldTestInputs(void);

  int              NModels       (void)             const;
  string           ModelTag      (int imodel)       const;
  TFile *          XSecFile      (int imodel)       const;
  string           XSecFileName  (int imodel)       const;
  TChain *         EvtChain      (int imodel)       const;
  vector<string> & EvtFileNames  (int imodel)       const;
  void             Print         (ostream & stream) const;
  bool             LoadFromFile  (string xmlfile);

  friend ostream & operator << (ostream & stream, const VldTestInputs & inp);  

private:

  void Init    (const int nmaxmodels);
  void CleanUp (void);

  bool                      fDoChain;
  int                       fNModels;
  vector<string>          * fModelTag;
  vector<TFile*>          * fXSecFile;
  vector<string>          * fXSecFileName;
  vector<TChain*>         * fEvtChain;
  vector<vector<string> > * fEvtFileNames;
};

}      // vld   namespace
}      // utils namespace
}      // genie namespace

#endif // _VLD_TEST_INPUTS_H_

