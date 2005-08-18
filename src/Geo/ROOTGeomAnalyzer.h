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


class TGeoVolume;

#include "GeomAnalyzerI.h"
#include <TGeoManager.h>

namespace genie {

class ROOTGeomAnalyzer : public GeomAnalyzerI {

public :

  ROOTGeomAnalyzer();
 ~ROOTGeomAnalyzer();
 void Load(char* filename);
 int SetVtxMaterial(char* material);

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

  void Initialize              (void);
  void BuildListOfTargetNuclei (void);
  TGeoVolume* GetWorldVolume(void);

  char*            fMaterial;            ///< [input] selected material for vertex
  
  TGeoManager *    fGeometry;            ///< [input] detector geometry
  TVector3 *       fCurrVertex;          ///< current generated vertex
  PathLengthList * fCurrPathLengthList;  ///< current list of path-lengths
  PDGCodeList *    fCurrPDGCodeList;     ///< current list of target nuclei
};

}      // genie namespace

#endif // _ROOT_GEOMETRY_ANALYZER_H_
