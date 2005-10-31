//____________________________________________________________________________
/*!

\class   genie::geometry::ROOTGeomAnalyzer

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

#include "EVGDrivers/GeomAnalyzerI.h"
#include "PDG/PDGUtils.h"

class TGeoVolume;
class TGeoMaterial;
class TGeoElement;

using std::string;

namespace genie    {
namespace geometry {

class ROOTGeomAnalyzer : public GeomAnalyzerI {

public :

  ROOTGeomAnalyzer(string geometry_filename);
  ROOTGeomAnalyzer(TGeoManager * gm);
 ~ROOTGeomAnalyzer();

  // set or enquire for analyzer configuration options

  void SetScannerNPoints    (int    np) { fNPoints    = np; }
  void SetScannerNRays      (int    nr) { fNRays      = nr; }
  void SetWeightWithDensity (bool   wt) { fDensWeight = wt; }
  void SetUnits             (double lu);
  void SetTopVolName        (string nm);

  int    ScannerNPoints    (void) const { return fNPoints;       }
  int    ScannerNRays      (void) const { return fNRays;         }
  bool   WeightWithDensity (void) const { return fDensWeight;    }
  double Units             (void) const { return fScale;         }
  string TopVolName        (void) const { return fTopVolumeName; }

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

  void   Initialize              (void);
  void   Load                    (string geometry_filename);
  void   Load                    (TGeoManager * gm);
  void   CleanUp                 (void);
  void   BuildListOfTargetNuclei (void);
  int    GetTargetPdgCode        (const TGeoMaterial * const m) const;
  int    GetTargetPdgCode        (const TGeoElement  * const e) const;
  void   ScalePathLengths        (PathLengthList & pl);
  double ComputeMaxPathLengthPDG (double* XYZ, double* direction, int pdgc);

  int              fMaterial;               ///< input selected material for vertex generation
  TGeoManager *    fGeometry;               ///< input detector geometry
  string           fTopVolumeName;          ///< input top volume name (if other that TGeoManager::GetTopVolume()]
  int              fNPoints;                ///< max path length scanner: points/surface [def:200]
  int              fNRays;                  ///< max path length scanner: rays/point [def:200]
  bool             fDensWeight;             ///< if true pathlengths are weighted with density [def:true]
  double           fScale;                  ///< conversion factor: input geometry units -> meters
  TVector3 *       fCurrVertex;             ///< current generated vertex
  PathLengthList * fCurrPathLengthList;     ///< current list of path-lengths
  PathLengthList * fCurrMaxPathLengthList;  ///< current list of max path-lengths
  PDGCodeList *    fCurrPDGCodeList;        ///< current list of target nuclei
  TGeoVolume *     fTopVolume;              ///< top volume 
};

}      // geometry namespace
}      // genie    namespace

#endif // _ROOT_GEOMETRY_ANALYZER_H_
