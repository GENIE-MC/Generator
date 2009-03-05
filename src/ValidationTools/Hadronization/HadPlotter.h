//____________________________________________________________________________
/*!

\class   HadPlots

\brief   Class to make data MC plots

\author  Tingjun Yang (Stanford Univ)

\created Feb 28, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HADPLOTTER_H_
#define _HADPLOTTER_H_

#include "ValidationTools/Hadronization/HadPlots.h"

#include <vector>
#include <string>

using namespace std;

namespace genie {
namespace vld   {

class HadPlotter {
      
public:
  
  HadPlotter() {};

  virtual ~HadPlotter() {}

  void AddPlots(HadPlots hp);
  
  void ShowPlots();

private:

  vector<HadPlots> hadPlots;


};
}//vld
}//genie

#endif
