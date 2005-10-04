//____________________________________________________________________________
/*!

\class   genie::PointGeomAnalyzer

\brief   The PointGeomAnalyzer class is the simplest implementation of the
         GeomAnalyserI interface and defines a simple 'point-like' geometry.

         Use this geometry analyzer to generate events when you do not want
         to use a detailed GEANT/ROOT geometry description but you only need
         to generate events for a 'single' nuclear target while you still want
         to use the GENIE MC job driver 'loaded' with a GENIE flux driver.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 14, 2005

*/
//____________________________________________________________________________

#ifndef _POINT_GEOMETRY_ANALYZER_H_
#define _POINT_GEOMETRY_ANALYZER_H_

#include "EVGDrivers/GeomAnalyzerI.h"

namespace genie {

class PointGeomAnalyzer : public GeomAnalyzerI {

public :

  PointGeomAnalyzer(int tgtpdgc);
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

  TVector3 *       fCurrVertex;          ///< current generated vertex
  PathLengthList * fCurrPathLengthList;  ///< current list of path-lengths
  PDGCodeList *    fCurrPDGCodeList;     ///< current list of target nuclei
};

}      // genie namespace

#endif // _POINT_GEOMETRY_ANALYZER_H_
