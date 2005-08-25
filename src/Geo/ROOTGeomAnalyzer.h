//____________________________________________________________________________
/*!

\class   genie::ROOTGeomAnalyzer

\brief   A ROOT/GEANT Geometry Analyzer

\author  Anselmo Meregaglia <anselmo.meregaglia@cern.ch>
         ETH Zurich

\created May 24, 2005

*/
//____________________________________________________________________________

#ifndef _ROOT_GEOMETRY_ANALYZER_H_
#define _ROOT_GEOMETRY_ANALYZER_H_

#include <string>

#include <TGeoManager.h>

#include "Geo/GeomAnalyzerI.h"
#include "PDG/PDGUtils.h"

class TGeoVolume;

using std::string;

namespace genie {

class ROOTGeomAnalyzer : public GeomAnalyzerI {

public :

  ROOTGeomAnalyzer(string filename);
 ~ROOTGeomAnalyzer();
 double SetVtxMaterial(int pdgc);
 double ComputeMaxPathLength(double* XYZ,double* direction,int pdgc);

  // implement the GeomAnalyzerI interface

  const PDGCodeList & ListOfTargetNuclei (void);

  const PathLengthList &
           ComputePathLengths
             (const TLorentzVector & x, const TLorentzVector & p); 

  const TVector3 &
           GenerateVertex
             (const TLorentzVector & x, const TLorentzVector & p, int tgtpdg);

   void test(void);

private:

  void Initialize              (string filename);
  void BuildListOfTargetNuclei (void);
 

  int              fMaterial;            ///< [input] selected material for vertex
  
  TGeoManager *    fGeometry;            ///< [input] detector geometry
  TVector3 *       fCurrVertex;          ///< current generated vertex
  PathLengthList * fCurrPathLengthList;  ///< current list of path-lengths
  PDGCodeList *    fCurrPDGCodeList;     ///< current list of target nuclei
};

}      // genie namespace

#endif // _ROOT_GEOMETRY_ANALYZER_H_
