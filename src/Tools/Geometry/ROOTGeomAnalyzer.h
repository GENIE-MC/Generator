//____________________________________________________________________________
/*!

\class    genie::geometry::ROOTGeomAnalyzer

\brief    A ROOT/GEANT4 geometry driver

\author   Anselmo Meregaglia <anselmo.meregaglia \at cern.ch>, ETH Zurich
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>, STFC, Rutherford Lab
          Robert Hatcher <rhatcher \at fnal.gov>, Fermilab

\created  May 24, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ROOT_GEOMETRY_ANALYZER_H_
#define _ROOT_GEOMETRY_ANALYZER_H_

#include <string>
#include <algorithm>

#include <TGeoManager.h>
#include <TVector3.h>

#include "Framework/EventGen/GeomAnalyzerI.h"
#include "Framework/ParticleData/PDGUtils.h"

class TGeoVolume;
class TGeoMaterial;
class TGeoMixture;
class TGeoElement;
class TGeoHMatrix;

using std::string;

namespace genie    {

class GFluxI;

namespace geometry {

class PathSegmentList;
class GeomVolSelectorI;

class ROOTGeomAnalyzer : public GeomAnalyzerI {

public :
  ROOTGeomAnalyzer(string geometry_filename);
  ROOTGeomAnalyzer(TGeoManager * gm);
  ROOTGeomAnalyzer() : GeomAnalyzerI() { ; } // used ONLY for derived class overloading
 ~ROOTGeomAnalyzer();

  /// implement the GeomAnalyzerI interface

  virtual const  PDGCodeList &    ListOfTargetNuclei    (void);
  virtual const  PathLengthList & ComputeMaxPathLengths (void);

  virtual const  PathLengthList & ComputePathLengths(const TLorentzVector & x, 
                                                     const TLorentzVector & p);
  virtual const  TVector3 &       GenerateVertex(const TLorentzVector & x, 
                                                 const TLorentzVector & p, int tgtpdg);

  /// set geometry driver's configuration options

  virtual void SetScannerNPoints    (int    np) { fNPoints    = np; } /* box  scanner */
  virtual void SetScannerNRays      (int    nr) { fNRays      = nr; } /* box  scanner */
  virtual void SetScannerNParticles (int    np) { fNParticles = np; } /* flux scanner */
  virtual void SetScannerFlux       (GFluxI* f) { fFlux       = f;  } /* flux scanner */
  virtual void SetWeightWithDensity (bool   wt) { fDensWeight = wt; }
  virtual void SetMixtureWeightsSum (double sum);
  virtual void SetLengthUnits       (double lu);
  virtual void SetDensityUnits      (double du);
  virtual void SetMaxPlSafetyFactor (double sf);
  virtual void SetTopVolName        (string nm);
  virtual void SetKeepSegPath       (bool keep) { fKeepSegPath = keep; }
  virtual void SetDebugFlags        (int  flgs) { fDebugFlags  = flgs; }

  /// retrieve geometry driver's configuration options

  virtual int           ScannerNPoints    (void) const { return fNPoints;           }
  virtual int           ScannerNRays      (void) const { return fNRays;             }
  virtual int           ScannerNParticles (void) const { return fNParticles;        }
  virtual bool          WeightWithDensity (void) const { return fDensWeight;        }
  virtual double        LengthUnits       (void) const { return fLengthScale;       }
  virtual double        DensityUnits      (void) const { return fDensityScale;      }
  virtual double        MixtureWeightsSum (void) const { return fMixtWghtSum;       }
  virtual double        MaxPlSafetyFactor (void) const { return fMaxPlSafetyFactor; }
  virtual string        TopVolName        (void) const { return fTopVolumeName;     }
  virtual TGeoManager * GetGeometry       (void) const { return fGeometry;          }
  virtual bool          GetKeepSegPath    (void) const { return fKeepSegPath;       }
  virtual const PathLengthList& GetMaxPathLengths(void) const { return *fCurrMaxPathLengthList; } // call only after ComputeMaxPathLengths() has been called 

  /// access to geometry coordinate/unit transforms for validation/test purposes

  virtual void   Local2SI      (PathLengthList & pl) const;
  virtual void   Local2SI      (TVector3 & v) const;
  virtual void   SI2Local      (TVector3 & v) const;
  virtual void   Master2Top    (TVector3 & v) const;
  virtual void   Master2TopDir (TVector3 & v) const;
  virtual void   Top2Master    (TVector3 & v) const;
  virtual void   Top2MasterDir (TVector3 & v) const;

  /// configure processing to perform path segment trimming

  virtual GeomVolSelectorI* AdoptGeomVolSelector (GeomVolSelectorI* selector) /// take ownership, return old
  { std::swap(selector,fGeomVolSelector); return selector; }


protected:

  virtual void   Initialize              (void);
  virtual void   CleanUp                 (void);
  virtual void   Load                    (string geometry_filename);
  virtual void   Load                    (TGeoManager * gm);
  virtual void   BuildListOfTargetNuclei (void);

  virtual int    GetTargetPdgCode        (const TGeoMaterial * const m) const;
  virtual int    GetTargetPdgCode        (const TGeoMixture * const m, int ielement) const;
  virtual double GetWeight               (const TGeoMaterial * mat, int pdgc);
  virtual double GetWeight               (const TGeoMixture * mixt, int pdgc);
  virtual double GetWeight               (const TGeoMixture * mixt, int ielement, int pdgc);

  virtual void   MaxPathLengthsFluxMethod(void);
  virtual void   MaxPathLengthsBoxMethod (void);
  virtual bool   GenBoxRay               (int indx, TLorentzVector& x4, TLorentzVector& p4);

  virtual double ComputePathLengthPDG    (const TVector3 & r, const TVector3 & udir, int pdgc);
  virtual void   SwimOnce                (const TVector3 & r, const TVector3 & udir);

  virtual bool   FindMaterialInCurrentVol(int pdgc);
  virtual bool   WillNeverEnter          (double step);
  virtual double StepToNextBoundary      (void);
  virtual double Step                    (void);
  virtual double StepUntilEntering       (void);



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
  TGeoHMatrix *    fMasterToTop;           ///< matrix connecting master coordinates to top volume coordinates
  bool             fMasterToTopIsIdentity; ///< is fMasterToTop matrix the identity matrix?

  bool             fKeepSegPath;           ///< need to fill path segment "path"
  PathSegmentList* fCurrPathSegmentList;   ///< current list of path-segments
  GeomVolSelectorI* fGeomVolSelector;      ///< optional path seg trimmer (owned)

  // used by GenBoxRay to retain history between calls
  TVector3         fGenBoxRayPos;
  TVector3         fGenBoxRayDir;
  int              fiface, fipoint, firay;
  bool             fnewpnt;
  double           fdx, fdy, fdz, fox, foy, foz;  ///< top vol size/origin (top vol units)
  
  // test purposes
  double           fmxddist, fmxdstep;   ///< max errors in pathsegmentlist
  int              fDebugFlags;

};

}      // geometry namespace
}      // genie    namespace

#endif // _ROOT_GEOMETRY_ANALYZER_H_
