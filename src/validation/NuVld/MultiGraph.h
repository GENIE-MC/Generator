//_____________________________________________________________________________
/*!

\class    genie::nuvld::MultiGraph

\brief    

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  January 2004
_______________________________________________________________________________*/

#ifndef _MULTI_GRAPH_H_
#define _MULTI_GRAPH_H_

#include <string>
#include <vector>

#include <TGraphAsymmErrors.h>
#include <TLegend.h>

using std::string;
using std::vector;

namespace genie {
namespace nuvld {

class MultiGraph 
{
public:

  MultiGraph();
  ~MultiGraph();

  void AddGraph(string legend_entry, TGraphAsymmErrors * graph);

  unsigned int        NGraphs         (void)                const;
  TGraphAsymmErrors * GetGraph        (unsigned int igraph) const;
  string              GetLegendEntry  (unsigned int igraph) const;

  TLegend * GetLegend  (const char * option) const;
  void      FillLegend (const char * option, TLegend * legend) const;
                                  
private:

  void  FormatGraph(unsigned int igraph);

  vector<TGraphAsymmErrors *> _mgraph;
  vector<string>              _legend_entries;
};

} // nuvld namespace
} // genie namespace

#endif
