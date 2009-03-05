//____________________________________________________________________________
/*!

\class   HadPlotter

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

#include <vector>
#include <string>

#include "ValidationTools/Hadronization/HadPlots.h"

class TH2D;
class TGraph;
class TCanvas;
class TLegend;

using namespace std;

namespace genie {
namespace vld_hadronization {

class HadPlotter {
      
public:
  HadPlotter(string data_file_directory = "");
  virtual ~HadPlotter();

  void AddPlots  (HadPlots hp);
  void ShowPlots (void);

private:

  TGraphErrors* MakeGraph (string file);
  TH2D*         DrawFrame (TCanvas *c1, int dir, 
                           double fr_xmin, double fr_xmax, double fr_ymin, double fr_ymax, 
                           string frtit, string frxtit, string frytit, bool logx, bool logy);
  void          DrawGraph (TGraph* gr, int mstyle, int mcol, double msize=0.8, string opt = "def");
  void          SetLegend (TLegend *leg);

  string fDataDir; ///< top level dir for data files
  vector<HadPlots> hadPlots;
};

} // vld_hadronization
} // genie

#endif
