//____________________________________________________________________________
/*!

\class   genie::geometry::PointGeomAnalyzer

\brief   The PointGeomAnalyzer class is the simplest implementation of the
         GeomAnalyserI interface and defines a simple 'point-like' geometry.

         Use this geometry analyzer to generate events when you do not want
         to use a detailed GEANT/ROOT geometry description but you only need
         to generate events for a 'single' nuclear target while you still want
         to use the GENIE MC job driver 'loaded' with a GENIE flux driver.
         The geometry can also support a mix of targets, each with its 
         corresponding weight.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created July 14, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _POINT_GEOMETRY_ANALYZER_H_
#define _POINT_GEOMETRY_ANALYZER_H_

#include <map>

#include "Framework/EventGen/GeomAnalyzerI.h"

using std::map;

namespace genie    {
namespace geometry {

class PointGeomAnalyzer : public GeomAnalyzerI {

public :
  PointGeomAnalyzer(int tgtpdgc);
  PointGeomAnalyzer(unsigned int n, const int tgt_pdg[], const double weight[]);
  PointGeomAnalyzer(const map<int,double> & tgtmap /* pdg -> weight*/);
 ~PointGeomAnalyzer();

  // implement the GeomAnalyzerI interface

  const PDGCodeList &    ListOfTargetNuclei    (void);
  const PathLengthList & ComputeMaxPathLengths (void);

  const PathLengthList &
           ComputePathLengths
             (const TLorentzVector & x, const TLorentzVector & p);
  const TVector3 &
           GenerateVertex
             (const TLorentzVector & x, const TLorentzVector & p, int tgtpdg);
private:

  void Initialize (const map<int,double> & tgtmap);
  void CleanUp    (void);

  TVector3 *       fCurrVertex;          ///< current generated vertex
  PathLengthList * fCurrPathLengthList;  ///< current list of path-lengths
  PDGCodeList *    fCurrPDGCodeList;     ///< current list of target nuclei
};

}      // geometry namespace
}      // genie    namespace

#endif // _POINT_GEOMETRY_ANALYZER_H_
