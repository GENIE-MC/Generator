//_____________________________________________________________________________
/*!

\class    genie::nuvld::GraphUtils

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _GRAPH_UTILS_H_
#define _GRAPH_UTILS_H_

#include <TGraph.h>
#include <TMath.h>

namespace genie {
namespace nuvld {

class GraphUtils {

public:

  //___________________________________________________________________________
  static void range(TGraph * graph, 
                    double & xmin, double & xmax, double & ymin, double & ymax)
  {
    const int N = graph->GetN();

    xmin = graph->GetX() [TMath::LocMin(N, graph->GetX())];
    xmax = graph->GetX() [TMath::LocMax(N, graph->GetX())];
    ymin = graph->GetY() [TMath::LocMin(N, graph->GetY())];
    ymax = graph->GetY() [TMath::LocMax(N, graph->GetY())];
  }
  //___________________________________________________________________________
};

} // nuvld namespace
} // genie namespace

#endif


