//____________________________________________________________________________
/*!

\class    genie::geometry::ROOTGeomAnalyzer

\brief    A ROOT/GEANT4 geometry driver

\author   Anselmo Meregaglia <anselmo.meregaglia@cern.ch>, ETH Zurich
          Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, STFC, Rutherford Lab

\created  May 24, 2005

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
class TGeoMixture;
class TGeoElement;
class TVector3;

using std::string;

namespace genie    {

class GFluxI;

namespace geometry {

class ROOTGeomAnalyzer : public GeomAnalyzerI {

public :
  ROOTGeomAnalyzer(string geometry_filename);
  ROOTGeomAnalyzer(TGeoManager * gm);
 ~ROOTGeomAnalyzer();

  // set or enquire for analyzer configuration options

  void SetScannerNPoints    (int    np) { fNPoints    = np; } /* box  scanner */
  void SetScannerNRays      (int    nr) { fNRays      = nr; } /* box  scanner */
  void SetScannerNParticles (int    np) { fNParticles = np; } /* flux scanner */
  void SetScannerFlux       (GFluxI* f) { fFlux       = f;  } /* flux scanner */
  void SetWeightWithDensity (bool   wt) { fDensWeight = wt; }
  void SetMixtureWeightsSum (double sum);
  void SetLengthUnits       (double lu);
  void SetDensityUnits      (double du);
  void SetMaxPlSafetyFactor (double sf);
  void SetTopVolName        (string nm);

  int     ScannerNPoints    (void) const { return fNPoints;           }
  int     ScannerNRays      (void) const { return fNRays;             }
  int     ScannerNParticles (void) const { return fNParticles;        }
  bool    WeightWithDensity (void) const { return fDensWeight;        }
  double  LengthUnits       (void) const { return fLengthScale;       }
  double  DensityUnits      (void) const { return fDensityScale;      }
  double  MixtureWeightsSum (void) const { return fMixtWghtSum;       }
  double  MaxPlSafetyFactor (void) const { return fMaxPlSafetyFactor; }
  string  TopVolName        (void) const { return fTopVolumeName;     }
  TGeoManager * GetGeometry (void) const { return fGeometry;          }

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
  void   MaxPathLengthsBoxMethod (void);
  void   MaxPathLengthsFluxMethod(void);
  int    GetTargetPdgCode        (const TGeoMaterial * const m) const;
  int    GetTargetPdgCode        (const TGeoElement  * const e) const;
  void   ScalePathLengths        (PathLengthList & pl);
  double ComputePathLengthPDG    (const TVector3 & r, const TVector3 & udir, int pdgc);
  double GetWeight               (TGeoMaterial * mat, int pdgc);
  double GetWeight               (TGeoMixture * mixt, int pdgc);
  double GetWeight               (TGeoMixture * mixt, int ielement, int pdgc);
  bool   WillNeverEnter          (double step);
  double StepToNextBoundary      (void);
  double Step                    (void);
  double StepUntilEntering       (void);

  int              fMaterial;              ///< input selected material for vertex generation
  TGeoManager *    fGeometry;              ///< input detector geometry
  string           fTopVolumeName;         ///< input top vol [other than TGeoManager::GetTopVolume()]
  int              fNPoints;               ///< max path length scanner (box method): points/surface [def:200]
  int              fNRays;                 ///< max path length scanner (box method): rays/point [def:200]
  int              fNParticles;            ///< max path length scanner (flux method): particles in [def:10000]
  GFluxI *         fFlux;                  ///< a flux objects that can be used to scan the max path lengths
  bool             fDensWeight;            ///< if true pathlengths are weighted with density [def:true]
  double           fLengthScale;           ///< conversion factor: input geometry length units -> meters
  double           fDensityScale;          ///< conversion factor: input geometry density units -> kgr/meters^3
  double           fMaxPlSafetyFactor;     ///< factor that can multiply the computed max path lengths
  double           fMixtWghtSum;           ///< norm of relative weights (<0 if explicit summing required)
  TVector3 *       fCurrVertex;            ///< current generated vertex
  PathLengthList * fCurrPathLengthList;    ///< current list of path-lengths
  PathLengthList * fCurrMaxPathLengthList; ///< current list of max path-lengths
  PDGCodeList *    fCurrPDGCodeList;       ///< current list of target nuclei
  TGeoVolume *     fTopVolume;             ///< top volume
};

}      // geometry namespace
}      // genie    namespace

#endif // _ROOT_GEOMETRY_ANALYZER_H_
