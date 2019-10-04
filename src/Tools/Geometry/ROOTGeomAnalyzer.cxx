//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Anselmo Meregaglia <anselmo.meregaglia \at cern.ch>, ETH Zurich
         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>, STFC - Rutherford Lab
         Robert Hatcher <rhatcher \at fnal.gov>, Fermilab

         May 24, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 28, 2008 - CA
   Slight code restructuring to make it easier keeping track of conversions
   between the SI and the current geometry system of units.
   Fixed a bug in density unit conversions reported by Elaine Schulte.
 @ Mar 28, 2008 - Pawel Guzowksi, Jim Dobson
   Fixed bug preventing the code crossing more than a fixed number of volume
   boundaries. The problem was not seen before as the code was tested with
   relatively simpler geometries (small number of volumes).
   Fixed problem with SetTopVolume() not actually setting it, so that the 
   option to specify sub-volumes of the master volume wasn't working properly.
 @ Nov 12, 2008 - CA, Jim Dobson
   Cleaned-up geometry navigation code. Improve the vertex placement method
   GenerateVertex() to avoid geo volume overshots.
 @ May 30, 2008 - CA, Jim Dobson
   Fixed path-lengths units problem in Local2SI(). The problem went undetected
   for long as the relevant quantities were used in ratios & the extra 
   multiplicative factor was cancelled out. It was spotted when we tried to
   work out T2K event sample normalization in terms of POTs.
 @ June 4, 2008 - CA
   Modified code forming pdg codes in case of mixtures: Getting A from
   TGeoMixture::GetAmixt()[ielement] rather than TGeoElement::GetA().
   Updated code enables GENIE to see all isotopes in the fixed nd280 geometry.
 @ August 29, 2008 - Pawel Guzowksi
   Fixed a long-standing limitation: Now the top volume can be any geometry
   volume and not just the master geometry volume. ROOT is placing the origin
   of the coordinatee system at the centre of the top volume. On the other
   hand GENIE assumed the the origin of the coordinate system is at the centre 
   of the master geometry volume. So, if the top volume is to be set to anything
   other than the master geometry volume, an appropriate transformation is 
   required for the input flux neutrino positions/directions and the generated 
   neutrino interaction vertices.
   Added Master2Top(v), Master2TopDir(v) and Top2Master(v) methods and the
   fMasterToTop TGeoHMatrix to store the transformation matrix between the 
   master and top volume coordinates.
 @ June 4, 2009 - RWH
   Refactorized code so swimming the same start point and direction only
   ever happens once for all the target PDG codes; do so generates a 
   PathSegmentList from whence the individual PDG lengths can be derived.  
   This also allows for future culling or trimming of path segments by an 
   externally supplied function as well as being more efficient.
 @ July 17, 2009 - RWH
   StepUntilEntering() wasn't accumulating total distance stepped if the 
   StepToNextBoudary() didn't leave the position IsEntering().  
   Also, small tweaks to LOG messages.
 @ July 27, 2009 - RWH
   The method StepToNextBoundary() doesn't actually move the current point 
   but just finds the boundary so the returned "step" size shouldn't be part 
   of the sum in StepUntilEntering().
 @ August 3, 2009 - RWH
   Handle case where initial step is from outside the geometry into the top
   volume.  It appears that the volume name is given as the top volume 
   (perhaps incorrectly), but the material (correctly) is null. Such steps 
   don't contribute to any material path length but document the total path.          
   Provision for keeping track of maximum ray distance and step size errors.
 @ August 28, 2009 - RWH
   Adjust to PathSegmentList move to genie::geometry namespace.  
   Add AdoptGeomVolSelector() function with takes ownership of GeomVolSelectorI*.  
   This selector can generate an trimmed PathSegmentList* which can be swapped 
   in for the original to limit what material is considered.
 @ February 4, 2010 - RWH
   Allow forcing fetch of geometry hierarchy path (or'ed with desire of the 
   GeomVolSelector, if any) but to get it by default as it is costly.  
   Fill the PathSegment medium and material upon entry along with volume name, 
   not at exit.
 @ February 4, 2011 - JD
   Change the way the topvolume is matched in SetTopVolName(string name). 
   Previously used TString::Contains("vol2match") which did not require the string
   length to be the same and sometime lead to degeneracies and selection of 
   incorrect top volume. Bug and fix were found by Kevin Connolly.   

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <set>

#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TList.h>
#include <TSystem.h>
#include <TMath.h>
#include <TPolyMarker3D.h>
#include <TGeoBBox.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Tools/Geometry/PathSegmentList.h"
#include "Framework/EventGen/PathLengthList.h"
#include "Framework/EventGen/GFluxI.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"
#include "Tools/Geometry/GeomVolSelectorI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

using namespace genie;
using namespace genie::geometry;
using namespace genie::controls;

//#define RWH_DEBUG
//#define RWH_DEBUG_2
//#define RWH_COUNTVOLS

#ifdef RWH_COUNTVOLS
// keep some statistics about how many volumes traversed for each box face
long int mxsegments = 0; //rwh
long int curface = 0; //rwh
long int nswims[6]  = { 0, 0, 0, 0, 0, 0}; //rwh
long int nnever[6]  = { 0, 0, 0, 0, 0, 0}; //rwh
double   dnvols[6]  = { 0, 0, 0, 0, 0, 0}; //rwh
double   dnvols2[6] = { 0, 0, 0, 0, 0, 0}; //rwh
bool  accum_vol_stat = false;
#endif

//___________________________________________________________________________
ROOTGeomAnalyzer::ROOTGeomAnalyzer(string geometry_filename)
  : GeomAnalyzerI()
{
///
/// Constructor from a geometry file
///
  LOG("GROOTGeom", pDEBUG)
    << "ROOTGeomAnalyzer ctor \"" << geometry_filename << "\"";
  this->Initialize();
  this->Load(geometry_filename);
}

//___________________________________________________________________________
ROOTGeomAnalyzer::ROOTGeomAnalyzer(TGeoManager * gm)
  : GeomAnalyzerI()
{
///
/// Constructor from a TGeoManager
///
  LOG("GROOTGeom", pDEBUG)
    << "ROOTGeomAnalyzer ctor passed TGeoManager*";
  this->Initialize();
  this->Load(gm);
}

//___________________________________________________________________________
ROOTGeomAnalyzer::~ROOTGeomAnalyzer()
{
  this->CleanUp();

  if ( fmxddist > 0 || fmxdstep > 0 )
    LOG("GROOTGeom",pNOTICE)
      << "ROOTGeomAnalyzer " 
      << " mxddist " << fmxddist
      << " mxdstep " << fmxdstep; 
}

//===========================================================================
// Geometry driver interface implementation:

//___________________________________________________________________________
const PDGCodeList & ROOTGeomAnalyzer::ListOfTargetNuclei(void)
{
  return *fCurrPDGCodeList;
}

//___________________________________________________________________________
const PathLengthList & ROOTGeomAnalyzer::ComputeMaxPathLengths(void)
{
/// Computes the maximum path lengths for all materials in the input 
/// geometry. The computed path lengths are in SI units (kgr/m^2, if 
/// density weighting is enabled)

  LOG("GROOTGeom", pNOTICE)
     << "Computing the maximum path lengths for all materials";

  if (!fGeometry) {
      LOG("GROOTGeom", pFATAL) << "No ROOT geometry is loaded!!";
      exit(1);
  }

  //-- initialize max path lengths
  fCurrMaxPathLengthList->SetAllToZero();

  //-- select maximum path length calculation method
  if ( fFlux ) {
    this->MaxPathLengthsFluxMethod();
    // clear any accumulated exposure accounted generated 
    // while exploring the geometry
    fFlux->Clear("CycleHistory");
  } else {
    this->MaxPathLengthsBoxMethod();
  }

  return *fCurrMaxPathLengthList;
}

//___________________________________________________________________________
const PathLengthList & ROOTGeomAnalyzer::ComputePathLengths(
                          const TLorentzVector & x, const TLorentzVector & p)
{
/// Computes the path-length within each detector material for a 
/// neutrino starting from point x (master coord) and travelling along 
/// the direction of p (master coord).
/// The computed path lengths are in SI units (kgr/m^2, if density 
/// weighting is enabled)

  //LOG("GROOTGeom", pDEBUG)
  //     << "Computing path-lengths for the input neutrino";

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
       << "\nInput nu: 4p (GeV) = " << utils::print::P4AsShortString(&p)
       << ", 4x (m,s) = " << utils::print::X4AsString(&x);
#endif

  // if trimming configure with neutrino ray's info
  if ( fGeomVolSelector ) {
    fGeomVolSelector->SetCurrentRay(x,p);
    fGeomVolSelector->SetSI2Local(1/this->LengthUnits());
  }

  TVector3 udir = p.Vect().Unit(); // unit vector along direction
  TVector3 pos = x.Vect();         // initial position
  this->SI2Local(pos);             // SI -> curr geom units
  
  if (!fMasterToTopIsIdentity) {
    this->Master2Top(pos);         // transform position (master -> top)
    this->Master2TopDir(udir);     // transform direction (master -> top)
  }

  // reset current list of path-lengths
  fCurrPathLengthList->SetAllToZero();

  //loop over materials & compute the path-length
  vector<int>::iterator itr;
  for (itr=fCurrPDGCodeList->begin();itr!=fCurrPDGCodeList->end();itr++) {

    int pdgc = *itr;

    Double_t pl = this->ComputePathLengthPDG(pos,udir,pdgc);
    fCurrPathLengthList->AddPathLength(pdgc,pl);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pINFO)
      <<"Calculated path length for material: " << pdgc << " = " << pl;
#endif

  } // loop over materials

  this->Local2SI(*fCurrPathLengthList); // curr geom units -> SI

  return *fCurrPathLengthList;
}

//___________________________________________________________________________
const TVector3 & ROOTGeomAnalyzer::GenerateVertex(
              const TLorentzVector & x, const TLorentzVector & p, int tgtpdg)
{
/// Generates a random vertex, within the detector material with the input
/// PDG code, for a neutrino starting from point x (master coord) and 
/// travelling along the direction of p (master coord).

  LOG("GROOTGeom", pNOTICE)
       << "Generating vtx in material: " << tgtpdg
       << " along the input neutrino direction";

  int nretry = 0;
  retry:  // goto label in case of abject failure
  nretry++;

  // reset current interaction vertex
  fCurrVertex->SetXYZ(0.,0.,0.);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
       << "\nInput nu: 4p (GeV) = " << utils::print::P4AsShortString(&p)
       << ", 4x (m,s) = " << utils::print::X4AsString(&x);
#endif

  if (!fGeometry) {
      LOG("GROOTGeom", pFATAL) << "No ROOT geometry is loaded!!";
      exit(1);
  }

  // calculate the max path length for the selected material starting from
  // x and looking along the direction of p
  TVector3 udir = p.Vect().Unit();
  TVector3 pos = x.Vect();
  this->SI2Local(pos);           // SI -> curr geom units

  if (!fMasterToTopIsIdentity) {
    this->Master2Top(pos);       // transform position (master -> top)
    this->Master2TopDir(udir);   // transform direction (master -> top)
  }

  double maxwgt_dist = this->ComputePathLengthPDG(pos,udir,tgtpdg);
  if ( maxwgt_dist <= 0 ) {
    LOG("GROOTGeom", pERROR)
     << "The current trajectory does not cross the selected material!!";
    return *fCurrVertex;
  }

  // generate random number between 0 and max_dist
  RandomGen * rnd = RandomGen::Instance();
  double genwgt_dist(maxwgt_dist * rnd->RndGeom().Rndm());

  LOG("GROOTGeom", pINFO)
    << "Swim mass: Top Vol dir = " << utils::print::P3AsString(&udir)
    << ", pos = " << utils::print::Vec3AsString(&pos);
  LOG("GROOTGeom", pINFO)
     << "Max {L x Density x Weight} given (init,dir) = " << maxwgt_dist;
  LOG("GROOTGeom", pINFO)
       << "Generated 'distance' in selected material = " << genwgt_dist;
#ifdef RWH_DEBUG
  if ( ( fDebugFlags & 0x01 ) ) {
    fCurrPathSegmentList->SetDoCrossCheck(true);       //RWH
    LOG("GROOTGeom", pINFO) << *fCurrPathSegmentList;  //RWH
    double mxddist = 0, mxdstep = 0;
    fCurrPathSegmentList->CrossCheck(mxddist,mxdstep);
    fmxddist = TMath::Max(fmxddist,mxddist);
    fmxdstep = TMath::Max(fmxdstep,mxdstep);
  }
#endif

  // compute the pdg weight for each material just once, then use a stl map 
  PathSegmentList::MaterialMap_t wgtmap;
  PathSegmentList::MaterialMapCItr_t mitr     = 
    fCurrPathSegmentList->GetMatStepSumMap().begin();
  PathSegmentList::MaterialMapCItr_t mitr_end = 
    fCurrPathSegmentList->GetMatStepSumMap().end();
  // loop over map to get tgt weight for each material (once)
  // steps outside the geometry may have no assigned material
  for ( ; mitr != mitr_end; ++mitr ) {
    const TGeoMaterial* mat = mitr->first;
    double wgt = ( mat ) ? this->GetWeight(mat,tgtpdg) : 0;
    wgtmap[mat] = wgt;
#ifdef RWH_DEBUG
    if ( ( fDebugFlags & 0x02 ) ) {
      LOG("GROOTGeom", pINFO)
        << " wgtmap[" << mat->GetName() << "] pdg " << tgtpdg << " wgt " << Form("%.6f",wgt);
    }
#endif
  }

  // walk down the path to pick the vertex
  const genie::geometry::PathSegmentList::PathSegmentV_t& segments = 
    fCurrPathSegmentList->GetPathSegmentV();
  genie::geometry::PathSegmentList::PathSegVCItr_t sitr;
  double walked = 0;
  for ( sitr = segments.begin(); sitr != segments.end(); ++sitr) {
    const genie::geometry::PathSegment& seg = *sitr;
    const TGeoMaterial* mat = seg.fMaterial;
    double trimmed_step = seg.GetSummedStepRange();
    double wgtstep = trimmed_step * wgtmap[mat];
    double beyond = walked + wgtstep;
#ifdef RWH_DEBUG
    if ( ( fDebugFlags & 0x04 ) ) {
      LOG("GROOTGeom", pINFO)
        << " beyond " << beyond << " genwgt_dist " << genwgt_dist
        << " trimmed_step " << trimmed_step << " wgtstep " << wgtstep;
    }
#endif
    if ( beyond > genwgt_dist ) {
      // the end of this segment is beyond our generation point
#ifdef RWH_DEBUG
      if ( ( fDebugFlags & 0x08 ) ) {
        LOG("GROOTGeom", pINFO)
          << "Choose vertex pos walked=" << walked 
          << " beyond=" << beyond 
          << " wgtstep " << wgtstep
          << " ( " << trimmed_step << "*" << wgtmap[mat] << ")"
          << " look for " << genwgt_dist
          << " in " << seg.fVolume->GetName() << " "
          << mat->GetName();
      }
#endif
      // choose a vertex in this segment (possibly multiple steps)
      double frac = ( genwgt_dist - walked ) / wgtstep;
      if ( frac > 1.0 ) {
        LOG("GROOTGeom", pWARN)
          << "Hey, frac = " << frac << " ( > 1.0 ) "
          << genwgt_dist << " " << walked << " " << wgtstep;
      }
      pos = seg.GetPosition(frac);
      fGeometry -> SetCurrentPoint (pos[0],pos[1],pos[2]);
      fGeometry -> FindNode();
      LOG("GROOTGeom", pINFO)
        << "Choose vertex position in " << seg.fVolume->GetName() << " "
         << utils::print::Vec3AsString(&pos);
      break;
    }
    walked = beyond;
  }

  LOG("GROOTGeom", pNOTICE)
     << "The vertex was placed in volume: " 
     << fGeometry->GetCurrentVolume()->GetName()
     << ", path: " << fGeometry->GetPath();

  // warn for any volume overshoots
  bool ok = this->FindMaterialInCurrentVol(tgtpdg);
  if (!ok) {
    LOG("GROOTGeom", pWARN)
       << "Geometry volume was probably overshot";
    LOG("GROOTGeom", pWARN)
       << "No material with code = " << tgtpdg << " could be found at genwgt_dist="
       << genwgt_dist << " (maxwgt_dist=" << maxwgt_dist << ")";
    if ( nretry < 10 ) {
      LOG("GROOTGeom", pWARN)
        << "retry placing vertex";
      goto retry;  // yeah, I know! MyBad.
    }
  }

  if (!fMasterToTopIsIdentity) {
     this->Top2Master(pos); // transform position (top -> master)
  }

  this->Local2SI(pos);   // curr geom units -> SI

  fCurrVertex->SetXYZ(pos[0],pos[1],pos[2]);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG) 
      << "Vtx (m) = " << utils::print::Vec3AsString(&pos);
#endif

  return *fCurrVertex;
}

//===========================================================================
// Driver configuration methods:

//___________________________________________________________________________
void ROOTGeomAnalyzer::SetLengthUnits(double u)
{
/// Use the units of the input geometry, 
///    e.g. SetLengthUnits(genie::units::centimeter)
/// GENIE uses the physical system of units (hbar=c=1) almost throughtout 
/// so everything is expressed in GeV but when analyzing detector geometries 
/// use meters. Setting input geometry units will allow the code to compute 
/// the conversion factor.
/// As input, use one of the constants in $GENIE/src/Conventions/Units.h

  fLengthScale = u/units::meter;
  LOG("GROOTGeom", pNOTICE)
     << "Geometry length units scale factor (geom units -> m): " 
     << fLengthScale;
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::SetDensityUnits(double u)
{
/// Like SetLengthUnits, but for density (default units = kgr/m3)

  fDensityScale = u / (units::kilogram / units::meter3);
  LOG("GROOTGeom", pNOTICE)
    << "Geometry density units scale factor (geom units -> kgr/m3): " 
    << fDensityScale;
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::SetMaxPlSafetyFactor(double sf)
{
/// Set a factor that can multiply the computed max path lengths.
/// The maximum path lengths are computed by performing an MC scanning of 
/// the input geometry. If you configure the scanner with a low number of 
/// points or rays you might understimate the path lengths, so you might 
/// want to 'inflate' them a little bit using this method.
/// Do not set this number too high, because the max interaction probability
/// will be grossly overestimated and you would need lots of attempts before
/// getting a flux neutrino to interact...

  fMaxPlSafetyFactor = sf;

  LOG("GROOTGeom", pNOTICE)
    << "Max path length safety factor: " << fMaxPlSafetyFactor;
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::SetMixtureWeightsSum(double sum)
{
/// Set it to x, if the relative weight proportions of elements in a mixture
/// add up to x (eg x=1, 100, etc). Set it to a negative value to explicitly 
/// compute the correct weight normalization.

  fMixtWghtSum = sum;
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::SetTopVolName(string name)
{
/// Set the name of the top volume.
/// This driver would ask the TGeoManager::GetTopVolume() for the top volume.
/// Use this method for changing this if for example you want to set a smaller
/// volume as the top one so as to generate events only in a specific part of
/// your detector.

  if (name.size() == 0) return;

  fTopVolumeName = name;
  LOG("GROOTGeom",pNOTICE) << "Geometry Top Volume name: " << fTopVolumeName;

  TGeoVolume * gvol = fGeometry->GetVolume(fTopVolumeName.c_str());
  if (!gvol) {
     LOG("GROOTGeom",pWARN) << "Could not find volume: " << name.c_str();
     LOG("GROOTGeom",pWARN) << "Will not change the current top volume";
     fTopVolumeName = "";
     return;
  }

  // Get a matrix connecting coordinates of master and top volumes.
  // The matrix will be used for transforming the coordinates of incoming 
  // flux neutrinos & generated interaction vertices.
  // This is needed (in case that the input top volume != master volume) 
  // because ROOT always sets the coordinate system origin at the centre of 
  // the specified top volume (whereas GENIE assumes that the global reference 
  // frame is that of the master volume)

  TGeoIterator next(fGeometry->GetMasterVolume());
  TGeoNode *node;
  TString nodeName, volNameStr;
  const char* volName = fTopVolumeName.c_str();
  while ((node = next())) {
    nodeName = node->GetVolume()->GetName();
    if (nodeName == volName) {
      if (fMasterToTop) delete fMasterToTop;
      fMasterToTop = new TGeoHMatrix(*next.GetCurrentMatrix());
      fMasterToTopIsIdentity = fMasterToTop->IsIdentity();
      break;
    }
  }

  // set volume name
  fTopVolume = gvol;
  fGeometry->SetTopVolume(fTopVolume);
}

//===========================================================================
// Geometry/Unit transforms:

//___________________________________________________________________________
void ROOTGeomAnalyzer::Local2SI(PathLengthList & pl) const
{
/// convert path lengths from current geometry units to SI units
///

  double scaling_factor = this->LengthUnits();
  if (this->WeightWithDensity()) { scaling_factor *= this->DensityUnits(); }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
    << "Scaling path-lengths from local units -> meters "
    << ((this->WeightWithDensity()) ? "* kgr/m^3" : "")
    << " - scale = " << scaling_factor;
#endif

  PathLengthList::iterator pliter;
  for(pliter = pl.begin(); pliter != pl.end(); ++pliter)
  {
    int pdgc = pliter->first;
    pl.ScalePathLength(pdgc, scaling_factor);
  }
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::Local2SI(TVector3 & vec) const
{
/// convert position vector from current geometry units to SI units
///

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Position (loc): " << utils::print::Vec3AsString(&vec);
#endif

    vec *= (this->LengthUnits()); 

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Position (SI): " << utils::print::Vec3AsString(&vec);
#endif
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::SI2Local(TVector3 & vec) const
{
/// convert position vector from SI units to current geometry units
///
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Position (SI): " << utils::print::Vec3AsString(&vec);
#endif

    vec *= (1./this->LengthUnits()); 

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Position (loc): " << utils::print::Vec3AsString(&vec);
#endif
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::Master2Top(TVector3 & vec) const
{
/// transform the input position vector from the master volume coordinate
/// system to the specified top volume coordinate system, but not 
/// change the units.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Position (coord:master): " << utils::print::Vec3AsString(&vec);
#endif

    Double_t mast[3], top[3];
    vec.GetXYZ(mast);
    fMasterToTop->MasterToLocal(mast, top);
    vec.SetXYZ(top[0], top[1], top[2]);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Position (coord:top): " << utils::print::Vec3AsString(&vec);
#endif
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::Master2TopDir(TVector3 & vec) const
{
/// transform the input direction vector from the master volume coordinate
/// system to the specified top volume coordinate system, but not
/// change the units.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Direction (coord:master): " << utils::print::Vec3AsString(&vec);
#endif

    Double_t mast[3], top[3];
    vec.GetXYZ(mast);
    fMasterToTop->MasterToLocalVect(mast, top);
    vec.SetXYZ(top[0], top[1], top[2]);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Direction (coord:top): " << utils::print::Vec3AsString(&vec);
#endif
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::Top2Master(TVector3 & vec) const
{
/// transform the input position vector from the specified top volume 
/// coordinate system to the master volume coordinate system, but not
/// change the units.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Position (coord:top): " << utils::print::Vec3AsString(&vec);
#endif

    Double_t mast[3], top[3];
    vec.GetXYZ(top);
    fMasterToTop->LocalToMaster(top, mast);
    vec.SetXYZ(mast[0], mast[1], mast[2]);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Position (coord:master): " << utils::print::Vec3AsString(&vec);
#endif
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::Top2MasterDir(TVector3 & vec) const
{
/// transform the input direction vector from the specified top volume 
/// coordinate system to the master volume coordinate system.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Direction (coord:top): " << utils::print::Vec3AsString(&vec);
#endif

    Double_t mast[3], top[3];
    vec.GetXYZ(top);
    fMasterToTop->LocalToMasterVect(top, mast);
    vec.SetXYZ(mast[0], mast[1], mast[2]);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Direction (coord:master): " << utils::print::Vec3AsString(&vec);
#endif
}

//===========================================================================
// Private methods:

//___________________________________________________________________________
void ROOTGeomAnalyzer::Initialize(void)
{
  LOG("GROOTGeom", pNOTICE)
                << "Initializing ROOT geometry driver & setting defaults";

  fCurrMaxPathLengthList = 0;
  fCurrPathLengthList    = 0;
  fCurrPathSegmentList   = 0;
  fGeomVolSelector       = 0;
  fCurrPDGCodeList       = 0;
  fTopVolume             = 0;
  fTopVolumeName         = "";
  fKeepSegPath           = false;

  // some defaults:
  this -> SetScannerNPoints    (200);
  this -> SetScannerNRays      (200);
  this -> SetScannerNParticles (10000);
  this -> SetScannerFlux       (0);
  this -> SetMaxPlSafetyFactor (1.1);
  this -> SetLengthUnits       (genie::units::meter);
  this -> SetDensityUnits      (genie::units::kilogram/genie::units::meter3);
  this -> SetWeightWithDensity (true);
  this -> SetMixtureWeightsSum (-1.);

  fMasterToTopIsIdentity = true;

  fmxddist = 0;
  fmxdstep = 0;
  fDebugFlags = 0;
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::CleanUp(void)
{
  LOG("GROOTGeom", pNOTICE) << "Cleaning up...";

  if ( fCurrPathSegmentList   ) delete fCurrPathSegmentList;
  if ( fCurrPathLengthList    ) delete fCurrPathLengthList;
  if ( fCurrMaxPathLengthList ) delete fCurrMaxPathLengthList;
  if ( fCurrPDGCodeList       ) delete fCurrPDGCodeList;
  if ( fMasterToTop           ) delete fMasterToTop;
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::Load(string filename)
{
/// Load the detector geometry from the input ROOT file
///
  LOG("GROOTGeom", pNOTICE) << "Loading geometry from: " << filename;

  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
  if (!is_accessible) {
     LOG("GROOTGeom", pFATAL)
       << "The ROOT geometry doesn't exist! Initialization failed!";
     exit(1);
  }
  TGeoManager * gm = TGeoManager::Import(filename.c_str());

  this->Load(gm);
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::Load(TGeoManager * gm)
{
/// Load the detector geometry from the input TGeoManager

  LOG("GROOTGeom", pNOTICE)
         << "A TGeoManager is being loaded to the geometry driver";
  fGeometry = gm;

  if (!fGeometry) {
    LOG("GROOTGeom", pFATAL) << "Null TGeoManager! Aborting";
  }
  assert(fGeometry);

  this->BuildListOfTargetNuclei();

  const PDGCodeList & pdglist = this->ListOfTargetNuclei();

  fTopVolume             = 0;
  fCurrPathSegmentList   = new PathSegmentList();
  fCurrPathLengthList    = new PathLengthList(pdglist);
  fCurrMaxPathLengthList = new PathLengthList(pdglist);
  fCurrVertex            = new TVector3(0.,0.,0.);

  // ask geometry manager for its top volume
  fTopVolume = fGeometry->GetTopVolume();
  if (!fTopVolume) {
      LOG("GROOTGeom", pFATAL) << "Could not get top volume!!!";
  }
  assert(fTopVolume);

  // load matrix (identity) of top volume
  fMasterToTop = new TGeoHMatrix(*fGeometry->GetCurrentMatrix());
  fMasterToTopIsIdentity = true;

//#define PRINT_MATERIALS
#ifdef PRINT_MATERIALS
  fGeometry->GetListOfMaterials()->Print();
  fGeometry->GetListOfMedia()->Print();
#endif

}

//___________________________________________________________________________
void ROOTGeomAnalyzer::BuildListOfTargetNuclei(void)
{
/// Determine possible target PDG codes.
/// Note:  If one is using a top volume other than the master level
/// then the final list might include PDG codes that will never
/// be seen during swimming through the volumes if those code are only
/// found in materials outside the top volume.

  fCurrPDGCodeList = new PDGCodeList;

  if (!fGeometry) {
    LOG("GROOTGeom", pFATAL) << "No ROOT geometry is loaded!!";
    exit(1);
  }

  TObjArray * volume_list = fGeometry->GetListOfVolumes();
  if (!volume_list) {
     LOG("GROOTGeom", pERROR)
        << "Null list of geometry volumes. Can not find build target list!";
     return;
  }

  std::set<Int_t> seen_mat; // list of materials we've already handled
  std::vector<TGeoVolume*> volvec; //RWH

  int numVol = volume_list->GetEntries();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG) << "Number of volumes found: " << numVol;
#endif

  for (int ivol = 0; ivol < numVol; ivol++) {
      TGeoVolume * volume = dynamic_cast <TGeoVolume *>(volume_list->At(ivol));
      if (!volume) {
         LOG("GROOTGeom", pWARN)
           << "Got a null geometry volume!! Skiping current list element";
         continue;
      }
      TGeoMaterial * mat = volume->GetMedium()->GetMaterial();

      // shortcut if we've already seen this material
      Int_t mat_indx = mat->GetIndex();
      if ( seen_mat.find(mat_indx) != seen_mat.end() ) continue;
      seen_mat.insert(mat_indx);
      volvec.push_back(volume); //RWH

      if (mat->IsMixture()) {
         TGeoMixture * mixt = dynamic_cast <TGeoMixture*> (mat);
         int Nelements = mixt->GetNelements();
         for (int i=0; i<Nelements; i++) {
            int ion_pdgc = this->GetTargetPdgCode(mixt,i);
            fCurrPDGCodeList->push_back(ion_pdgc);
         }
      } else {
          int ion_pdgc = this->GetTargetPdgCode(mat);
          fCurrPDGCodeList->push_back(ion_pdgc);
      }
  }
  // sort the list
  // we don't calculate this list but once per geometry and a sorted
  // list is easier to read so this doesn't cost much
  std::sort(fCurrPDGCodeList->begin(),fCurrPDGCodeList->end());

}

//___________________________________________________________________________
int ROOTGeomAnalyzer::GetTargetPdgCode(const TGeoMaterial * const m) const
{
  int A = TMath::Nint(m->GetA());
  int Z = TMath::Nint(m->GetZ());

  int pdgc = pdg::IonPdgCode(A,Z);

  return pdgc;
}

//___________________________________________________________________________
int ROOTGeomAnalyzer::GetTargetPdgCode(
                     const TGeoMixture * const m, int ielement) const
{
  int A = TMath::Nint(m->GetAmixt()[ielement]);
  int Z = TMath::Nint(m->GetZmixt()[ielement]);

  int pdgc = pdg::IonPdgCode(A,Z);

  return pdgc;
}

//___________________________________________________________________________
double ROOTGeomAnalyzer::GetWeight(const TGeoMaterial * mat, int pdgc)
{
/// Get the weight of the input material.
/// Return the weight only if the material's pdg code matches the input code.
/// If the material is found to be a mixture, call the corresponding method
/// for mixtures.
/// Weight is in the curr geom density units.

  if (!mat) {
    LOG("GROOTGeom", pERROR) << "Null input material. Return weight = 0.";
    return 0;
  }

  bool exists = fCurrPDGCodeList->ExistsInPDGCodeList(pdgc);
  if (!exists) {
    LOG("GROOTGeom", pERROR) << "Target doesn't exist. Return weight = 0.";
    return 0;
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom",pDEBUG)
       << "Curr. material: A/Z = " << mat->GetA() << " / " << mat->GetZ();
#endif

  // if the input material is a mixture, get a the sum of weights for
  // all matching elements
  double weight = 0.;
  if (mat->IsMixture()) {
    const TGeoMixture * mixt = dynamic_cast <const TGeoMixture*> (mat);
    if (!mixt) {
     LOG("GROOTGeom", pERROR) << "Null input mixture. Return weight = 0.";
     return 0;
    }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Material : " << mat->GetName()
          << " is a mixture with " << mixt->GetNelements() << " elements";
#endif
    // loop over elements & sum weights of matching elements
    weight = this->GetWeight(mixt,pdgc);
    return weight;
  } // is mixture?

  // pure material
  int ion_pdgc = this->GetTargetPdgCode(mat);
  if (ion_pdgc != pdgc) return 0.;

  if (this->WeightWithDensity()) 
    weight = mat->GetDensity(); // material density (curr geom units)
  else                           
    weight = 1.0;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
       << "Weight[mat:" << mat->GetName() << "] = " << weight;
#endif

  return weight;
}

//___________________________________________________________________________
double ROOTGeomAnalyzer::GetWeight(const TGeoMixture * mixt, int pdgc)
{
/// Loop over the mixture elements, find the one matching the input pdgc
/// and  return its weight.
/// Weight is in the curr geom density units.

  double weight = 0;

  int nm = 0;
  for (int i = 0; i < mixt->GetNelements(); i++) {
     double dw = (this->GetWeight(mixt,i,pdgc));
     if (dw>0) nm++;
     weight += dw;
  }

  if (nm>1) {
     for (int j = 0; j < mixt->GetNelements(); j++) {
           LOG("GROOTGeom", pWARN)
              << "[" << j << "] Z = " << mixt->GetZmixt()[j] 
              << ", A = " << mixt->GetAmixt()[j]
              << " (pdgc = " << this->GetTargetPdgCode(mixt,j)
              << "), w = " << mixt->GetWmixt()[j];
     }
     LOG("GROOTGeom", pERROR)
        << "Material pdgc = " << pdgc << " appears " << nm
        << " times (>1) in mixture = " << mixt->GetName();
     LOG("GROOTGeom", pFATAL)
        << "Your geometry must be incorrect - Aborting";
     exit(1);
  }

  // if we are not weighting with the density then the weight=1 if the pdg
  // code was matched for any element of this mixture
  if ( !this->WeightWithDensity() && weight>0. ) weight=1.0;

  return weight;
}

//___________________________________________________________________________
double ROOTGeomAnalyzer::GetWeight(const TGeoMixture* mixt, int ielement, int pdgc)
{
/// Get the weight of the input ith element of the input material.
/// Return the weight only if the element's pdg code matches the input code.
/// Weight is in the curr geom density units.

//  int ion_pdgc = this->GetTargetPdgCode(mixt->GetElement(ielement));
  int ion_pdgc = this->GetTargetPdgCode(mixt, ielement);
  if (ion_pdgc != pdgc) return 0.;

  double d = mixt->GetDensity();         // mixture density (curr geom units)
  double w = mixt->GetWmixt()[ielement]; // relative proportion by mass

  double wtot = this->MixtureWeightsSum();

  // <0 forces explicit calculation of relative proportion normalization
  if (wtot < 0) {
    wtot = 0;
    for (int i = 0; i < mixt->GetNelements(); i++) {
      wtot += (mixt->GetWmixt()[i]);
    }
  }
  assert(wtot>0);

  w /= wtot;
  double weight = d*w;

  return weight;
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::MaxPathLengthsFluxMethod(void)
{
/// Use the input flux driver to generate "rays", and then follow them through
/// the detector and figure out the maximum path length for each material

  LOG("GROOTGeom", pNOTICE)
               << "Computing the maximum path lengths using the FLUX method";

  int iparticle = 0;
  PathLengthList::const_iterator pl_iter;

  const int nparticles = abs(this->ScannerNParticles());

  // if # scanner particles is negative, this signals that the user
  // desires to force rays to have the maximum energy (useful if the
  // volume considered changes size with neutrino energy)
  bool rescale_e = (this->ScannerNParticles() < 0 );
  double emax = fFlux->MaxEnergy();
  if ( rescale_e ) {
    LOG("GROOTGeom", pNOTICE)
      << "max path lengths with FLUX method forcing Enu=" << emax;
  }

  while (iparticle < nparticles ) {

    bool ok = fFlux->GenerateNext();
    if (!ok) {
       LOG("GROOTGeom", pWARN) << "Couldn't generate a flux neutrino";
       continue;
    }

    TLorentzVector   nup4  = fFlux->Momentum();
    if ( rescale_e ) {
      double ecurr = nup4.E();
      if ( ecurr > 0 ) nup4 *= (emax/ecurr);
    }
    const TLorentzVector & nux4  = fFlux->Position();

    //LOG("GMCJDriver", pNOTICE)
    //   << "\n [-] Generated flux neutrino: "
    //   << "\n  |----o 4-momentum : " << utils::print::P4AsString(&nup4)
    //   << "\n  |----o 4-position : " << utils::print::X4AsString(&nux4);

    const PathLengthList & pl = this->ComputePathLengths(nux4, nup4);

    bool enters = false;

    for (pl_iter = pl.begin(); pl_iter != pl.end(); ++pl_iter) {
       int    pdgc       = pl_iter->first;
       double pathlength = pl_iter->second;

       if ( pathlength > 0 ) {
          pathlength *= (this->MaxPlSafetyFactor());

          pathlength = TMath::Max(pathlength, fCurrMaxPathLengthList->PathLength(pdgc));
          fCurrMaxPathLengthList->SetPathLength(pdgc,pathlength);
          enters = true;
       }
    }
    if (enters) iparticle++;  
  }
}

//___________________________________________________________________________
void ROOTGeomAnalyzer::MaxPathLengthsBoxMethod(void)
{
/// Generate points in the geometry's bounding box and for each point 
/// generate random rays, follow them through the detector and figure out 
/// the maximum path length for each material

  LOG("GROOTGeom", pNOTICE)
    << "Computing the maximum path lengths using the BOX method";
#ifdef RWH_COUNTVOLS
  accum_vol_stat = true;
#endif

  int  iparticle = 0;
  bool ok = true;
  TLorentzVector nux4;
  TLorentzVector nup4;

  PathLengthList::const_iterator pl_iter;

  while ( (ok = this->GenBoxRay(iparticle++,nux4,nup4)) ) {

    //LOG("GMCJDriver", pNOTICE)
    //  << "\n [-] Generated flux neutrino: "
    //  << "\n  |----o 4-momentum : " << utils::print::P4AsString(&nup4)
    //  << "\n  |----o 4-position : " << utils::print::X4AsString(&nux4);

    const PathLengthList & pllst = this->ComputePathLengths(nux4, nup4);

    for (pl_iter = pllst.begin(); pl_iter != pllst.end(); ++pl_iter) {
       int    pdgc = pl_iter->first;
       double pl   = pl_iter->second;

       if (pl>0) {
          pl *= (this->MaxPlSafetyFactor());

          pl = TMath::Max(pl, fCurrMaxPathLengthList->PathLength(pdgc));
          fCurrMaxPathLengthList->SetPathLength(pdgc,pl);
       }
    }
  }

  // print out the results
  LOG("GROOTGeom", pDEBUG)
    << "DensWeight \"" << (fDensWeight?"true":"false") 
    << "\" MixtWghtSum " << fMixtWghtSum;
  LOG("GROOTGeom", pDEBUG) << "CurrMaxPathLengthList: "
    << *fCurrMaxPathLengthList;

#ifdef RWH_COUNTVOLS
  // rwh
  // print statistics for average,rms of number of volumes seen for
  // various rays for each face
  for (int j = 0; j < 6; ++j ) {
    long int ns = nswims[j];
    double   x  = dnvols[j];
    double   x2 = dnvols2[j];
    if ( ns == 0 ) ns = 1;
    double avg = x / (double)ns;
    double rms = TMath::Sqrt((x2/(double)ns) - avg*avg);
    LOG("GROOTGeom", pNOTICE)
      << "RWH: nswim after BOX face " << j << " is " << ns
      << " avg " << avg << " rms " << rms 
      << " never " << nnever[j];
  }
  LOG("GROOTGeom", pNOTICE)
    << "RWH: Max PathSegmentList size " << mxsegments;
  accum_vol_stat = false;
#endif

}

//___________________________________________________________________________
bool ROOTGeomAnalyzer::GenBoxRay(int indx, TLorentzVector& x4, TLorentzVector& p4)
{
/// Generate points in the geometry's bounding box and for each point generate
/// random rays -- a pseudo-flux -- in master coordinates and SI units

  firay++;
  fnewpnt = false;

  // first time through ... special case
  if ( indx == 0 )           { fiface = 0; fipoint = 0; firay = 0; fnewpnt = true; }

  if ( firay   >= fNRays   ) {             fipoint++;   firay = 0; fnewpnt = true; }
  if ( fipoint >= fNPoints ) { fiface++;   fipoint = 0; firay = 0; fnewpnt = true; }
  if ( fiface  >= 3 ) {
    LOG("GROOTGeom",pINFO) << "Box surface scanned: " << indx << " points/rays";
    return false; // signal end
  }

  if ( indx == 0 ) {
    // get geometry's bounding box
    //LOG("GROOTGeom", pNOTICE) << "Getting a TGeoBBox enclosing the detector";
    TGeoShape * TS  = fTopVolume->GetShape();
    TGeoBBox *  box = (TGeoBBox *)TS;
    //get box origin and dimensions (in the same units as the geometry)
    fdx = box->GetDX(); // half-length
    fdy = box->GetDY(); // half-length
    fdz = box->GetDZ(); // half-length
    fox = (box->GetOrigin())[0];
    foy = (box->GetOrigin())[1];
    foz = (box->GetOrigin())[2];
    
    LOG("GROOTGeom",pNOTICE) 
      << "Box size (GU)   :"
      << " x = " << 2*fdx << ", y = " << 2*fdy << ", z = " << 2*fdz;
    LOG("GROOTGeom",pNOTICE) 
      << "Box origin (GU) :"
      << " x = " << fox   << ", y = " << foy   << ", z = " <<   foz;
    LOG("GROOTGeom",pNOTICE)
      << "Will generate [" << fNPoints << "] random points / box surface";
    LOG("GROOTGeom",pNOTICE)
      << "Will generate [" << fNRays   << "] rays / point";

#ifdef VALIDATE_CORNERS
    // RWH: test that we know the coordinate transforms for the box corners
    for (int sz = -1; sz <= +1; ++sz) {
      for (int sy = -1; sy <= +1; ++sy) {
        for (int sx = -1; sx <= +1; ++sx) {
          if (sx == 0 || sy == 0 || sz == 0 ) continue;
          TVector3& pos = fGenBoxRayPos;
          pos.SetXYZ(fox+(sx*fdx),foy+(sy*fdy),foz+(sz*fdz));
          TVector3 master(fGenBoxRayPos);
          this->Top2Master(master);   // transform position (top -> master)
          this->Local2SI(master);
          TVector3 pos2(master);
          this->Master2Top(pos2);
          this->SI2Local(pos2);
          LOG("GROOTGeom", pNOTICE)  
            << "TopVol corner "
            << " [" << pos[0] << "," << pos[1] << "," << pos[2] << "] "
            << "Master corner "
            << " [" << master[0] << "," << master[1] << "," << master[2] << "] "
            << " top again" 
            << " [" << pos2[0] << "," << pos2[1] << "," << pos2[2] << "] ";
        }
      }
    }
#endif

  }

  RandomGen * rnd = RandomGen::Instance();

  double eps = 0.0; //1.0e-12; // tweak to make sure we're inside the box

  switch ( fiface ) {
  case 0: 
    {
      //top:
      if ( firay == fNRays-1 && fipoint == fNPoints-1 )
        LOG("GROOTGeom",pINFO) << "Box surface scanned: [TOP]";
      fGenBoxRayDir.SetXYZ(-0.5+rnd->RndGeom().Rndm(), 
                           -rnd->RndGeom().Rndm(), 
                           -0.5+rnd->RndGeom().Rndm());
      if ( fnewpnt )
        fGenBoxRayPos.SetXYZ(fox-fdx+eps+2*(fdx-eps)*rnd->RndGeom().Rndm(), 
                             foy+fdy-eps, 
                             foz-fdz+eps+2*(fdz-eps)*rnd->RndGeom().Rndm());
      break;
    } 
  case 1: 
    {
      //left:
      if ( firay == fNRays-1 && fipoint == fNPoints-1 )
        LOG("GROOTGeom",pINFO) << "Box surface scanned: [LEFT]";
      fGenBoxRayDir.SetXYZ(rnd->RndGeom().Rndm(), 
                           -0.5+rnd->RndGeom().Rndm(), 
                           -0.5+rnd->RndGeom().Rndm());
      if ( fnewpnt )
        fGenBoxRayPos.SetXYZ(fox-fdx+eps, 
                             foy-fdy+eps+2*(fdy-eps)*rnd->RndGeom().Rndm(),
                             foz-fdz+eps+2*(fdz-eps)*rnd->RndGeom().Rndm());
      break;
    } 
  case 2: 
    {
      //front: (really what I, RWH, would call the back)
      if ( firay == fNRays-1 && fipoint == fNPoints-1 )
        LOG("GROOTGeom",pINFO) << "Box surface scanned: [FRONT]";
      fGenBoxRayDir.SetXYZ(-0.5+rnd->RndGeom().Rndm(), 
                           -0.5+rnd->RndGeom().Rndm(), 
                           -rnd->RndGeom().Rndm());
      if ( fnewpnt )
        fGenBoxRayPos.SetXYZ(fox-fdx+eps+2*(fdx-eps)*rnd->RndGeom().Rndm(), 
                             foy-fdy+eps+2*(fdy-eps)*rnd->RndGeom().Rndm(), 
                             foz+fdz+eps);
      break;
    }
  }
/*
    //bottom:
      pos.SetXYZ(ox-dx+2*dx*rnd->RndGeom().Rndm(), 
                 oy-dy, 
                 oz-dz+2*dz*rnd->RndGeom().Rndm());
      dir.SetXYZ(-0.5+rnd->RndGeom().Rndm(), 
                 rnd->RndGeom().Rndm(), 
                 -0.5+rnd->RndGeom().Rndm());
    //right:
      pos.SetXYZ(ox+dx, 
                 oy-dy+2*dy*rnd->RndGeom().Rndm(), 
                 oz-dz+2*dz*rnd->RndGeom().Rndm());
      dir.SetXYZ(-rnd->RndGeom().Rndm(), 
                 -0.5+rnd->RndGeom().Rndm(),
                 -0.5+rnd->RndGeom().Rndm());
    //back:
      pos.SetXYZ(ox-dx+2*dx*rnd->RndGeom().Rndm(), 
                 oy-dy+2*dy*rnd->RndGeom().Rndm(), 
                 oz-dz);
      dir.SetXYZ(-0.5+rnd->RndGeom().Rndm(), 
                 -0.5+rnd->RndGeom().Rndm(), 
                  rnd->RndGeom().Rndm());
*/

#ifdef RWH_COUNTVOLS
  curface = fiface;
#endif

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  if ( fnewpnt )
    LOG("GROOTGeom", pNOTICE)  
      << "GenBoxRay(topvol) "
      << " iface " << fiface << " ipoint " << fipoint << " iray " << firay
      << " newpnt " << (fnewpnt?"true":"false")
      << " x " << utils::print::Vec3AsString(&fGenBoxRayPos)
      << " p " << utils::print::P3AsString(&fGenBoxRayDir);
#endif

  if ( fnewpnt ) {
    if ( ! fMasterToTopIsIdentity) {
      this->Top2Master(fGenBoxRayPos); // transform position (top -> master)
    }
    this->Local2SI(fGenBoxRayPos);
  }
  this->Top2MasterDir(fGenBoxRayDir);  // transform direction (top -> master)

  x4.SetVect(fGenBoxRayPos);
  p4.SetVect(fGenBoxRayDir.Unit());

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pNOTICE)  
    << "GenBoxRay(master) "
    << " iface " << fiface << " ipoint " << fipoint << " iray " << firay
    << " newpnt " << (fnewpnt?"true":"false")
    << " x " << utils::print::Vec3AsString(&fGenBoxRayPos)
    << " p " << utils::print::P3AsString(&fGenBoxRayDir);
#endif

  return true;
}

//________________________________________________________________________
double ROOTGeomAnalyzer::ComputePathLengthPDG(
                  const TVector3 & r0, const TVector3 & udir, int pdgc)
{
/// Compute the path length for the material with pdg-code = pdc, staring 
/// from the input position r (top vol coord & units) and moving along the 
/// direction of the unit vector udir (top vol coord).

  double pl = 0; // path-length (x density, if density-weighting is ON)

  this->SwimOnce(r0,udir);

  double step   = 0;
  double weight = 0;

  //  const TGeoVolume   * vol = 0;
  //  const TGeoMedium   * med = 0;
  const TGeoMaterial * mat = 0;

  // loop over independent materials, which is shorter or equal to # of volumes
  PathSegmentList::MaterialMapCItr_t itr     = 
    fCurrPathSegmentList->GetMatStepSumMap().begin();
  PathSegmentList::MaterialMapCItr_t itr_end = 
    fCurrPathSegmentList->GetMatStepSumMap().end();
  for ( ; itr != itr_end; ++itr ) {
    mat  = itr->first;
    if ( ! mat ) continue;  // segment outside geometry has no material
    step = itr->second;
    weight = this->GetWeight(mat,pdgc);
    pl += (step*weight);
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG) 
       << "PathLength[" << pdgc << "] = " << pl << " in curr geom units";
#endif

  return pl;
}

//________________________________________________________________________
void ROOTGeomAnalyzer::SwimOnce(const TVector3 & r0, const TVector3 & udir)
{
/// Swim through the geometry from the from the input position 
/// r0 (top vol coord & units) and moving along the direction of the 
/// unit vector udir (topvol coord) to create a filled PathSegmentList

  int nvolswim = 0; //rwh

  if ( ! fCurrPathSegmentList ) fCurrPathSegmentList = new PathSegmentList();

  // don't swim if the current PathSegmentList is up-to-date
  if ( fCurrPathSegmentList->IsSameStart(r0,udir) ) return;

  // start fresh
  fCurrPathSegmentList->SetAllToZero();

  // set start info so next time we don't swim for the same ray 
  fCurrPathSegmentList->SetStartInfo(r0,udir);
 
  PathSegment ps_curr;

  bool found_vol (false);
  bool keep_on   (true);

  double step    = 0;
  double raydist = 0;

  const TGeoVolume *   vol = 0;
  const TGeoMedium *   med = 0;
  const TGeoMaterial * mat = 0;

  // determining the geometry path is expensive, do it only if necessary
  bool selneedspath = ( fGeomVolSelector && fGeomVolSelector->GetNeedPath() );
  const bool fill_path = fKeepSegPath || selneedspath;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pNOTICE)  
    << "SwimOnce x [" << r0[0] << "," << r0[1] << "," << r0[2]
    << "] udir [" << udir[0] << "," << udir[1] << "," << udir[2];
#endif

  fGeometry -> SetCurrentDirection (udir[0],udir[1],udir[2]);
  fGeometry -> SetCurrentPoint     (r0[0],  r0[1],  r0[2]  );

  while (!found_vol || keep_on) {
     keep_on = true;

     fGeometry->FindNode();

     ps_curr.SetEnter( fGeometry->GetCurrentPoint() , raydist );
     vol = fGeometry->GetCurrentVolume();
     med = vol->GetMedium();
     mat = med->GetMaterial();
     ps_curr.SetGeo(vol,med,mat);
#ifdef PATHSEG_KEEP_PATH
     if (fill_path) ps_curr.SetPath(fGeometry->GetPath());
#endif

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
#ifdef DUMP_SWIM
       LOG("GROOTGeom", pDEBUG) << "Current volume: " << vol->GetName()
                             << " pos " << fGeometry->GetCurrentPoint()[0]
                             << " "     << fGeometry->GetCurrentPoint()[1]
                             << " "     << fGeometry->GetCurrentPoint()[2]
                             << " dir " << fGeometry->GetCurrentDirection()[0]
                             << " "     << fGeometry->GetCurrentDirection()[1]
                             << " "     << fGeometry->GetCurrentDirection()[2]
                             << "[path: " << fGeometry->GetPath() << "]";
#endif
#endif

     // find the start of top
     if (fGeometry->IsOutside() || !vol) {
        keep_on = false;
        if (found_vol) break;
        step = 0;
          this->StepToNextBoundary();
        //rwh//raydist += step;  // STNB doesn't actually "step"

#ifdef RWH_DEBUG
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("GROOTGeom", pDEBUG) << "Outside ToNextBoundary step: " << step 
                                 << " raydist: " << raydist;
#endif
#endif

        while (!fGeometry->IsEntering()) {
          step = this->Step();
          raydist += step;
#ifdef RWH_DEBUG
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("GROOTGeom", pDEBUG) 
              << "Stepping... [step size = " << step << "]";
          LOG("GROOTGeom", pDEBUG) << "Outside step: " << step 
                                   << " raydist: " << raydist;
#endif
#endif
          if (this->WillNeverEnter(step)) {
#ifdef RWH_COUNTVOLS
            if ( accum_vol_stat ) {
              // this really shouldn't happen for the box exploration... 
              // if coord transforms done right
              // could happen for neutrinos on a flux window
              nnever[curface]++;  //rwh
              if ( nnever[curface]%21 == 0 )
                LOG("GROOTGeom", pNOTICE)  
                  << "curface " << curface << " " << nswims[curface]
                  << " never " << nnever[curface]
                  << " x [" << r0[0] << "," << r0[1] << "," << r0[2] << "] "
                  << " p [" << udir[0] << "," << udir[1] << "," << udir[2] << "]";
            }
#endif
            fCurrPathSegmentList->SetAllToZero();            
            return;
          }
        } // finished while

        ps_curr.SetExit(fGeometry->GetCurrentPoint());
        ps_curr.SetStep(step);
        if ( ( fDebugFlags & 0x10 ) ) {
          // In general don't add the path segments from the start point to
          // the top volume (here for debug purposes)
          // Clear out the step range even if we keep it
          ps_curr.fStepRangeSet.clear();
          LOG("GROOTGeom", pNOTICE)
            << "debug: step towards top volume: " << ps_curr;
          fCurrPathSegmentList->AddSegment(ps_curr);
        }

     }  // outside or !vol

     if (keep_on) {
       if (!found_vol) found_vol = true;

       step   = this->StepUntilEntering();
       raydist += step;

       ps_curr.SetExit(fGeometry->GetCurrentPoint());
       ps_curr.SetStep(step);
       fCurrPathSegmentList->AddSegment(ps_curr);

       nvolswim++; //rwh

#ifdef DUMP_SWIM
       LOG("GROOTGeom", pDEBUG) << "Current volume: " << vol->GetName() 
                                << " step " << step << " in " << mat->GetName();
#endif
         
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("GROOTGeom", pDEBUG)
         << "Cur med.: " << med->GetName() << ", mat.: " << mat->GetName();
       LOG("GROOTGeom", pDEBUG)
         << "Step = " << step; // << ", weight = " << weight;
#endif
     } //keep on
  }

#ifdef RWH_COUNTVOLS
  if ( accum_vol_stat ) {
    nswims[curface]++;   //rwh
    dnvols[curface]  += (double)nvolswim;
    dnvols2[curface] += (double)nvolswim * (double)nvolswim;
    long int ns = fCurrPathSegmentList->size();
    if ( ns > mxsegments ) mxsegments = ns;
  }
#endif

//rwh:debug
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
    << "PathSegmentList size " << fCurrPathSegmentList->size();
#endif

#ifdef RWH_DEBUG_2
  if ( ( fDebugFlags & 0x20 ) ) {
    fCurrPathSegmentList->SetDoCrossCheck(true);       //RWH
    LOG("GROOTGeom", pNOTICE) << "Before trimming" << *fCurrPathSegmentList;
    double mxddist = 0, mxdstep = 0;
    fCurrPathSegmentList->CrossCheck(mxddist,mxdstep);
    fmxddist = TMath::Max(fmxddist,mxddist);
    fmxdstep = TMath::Max(fmxdstep,mxdstep);
  }
#endif

  // PathSegmentList trimming occurs here!
  if ( fGeomVolSelector ) {
    PathSegmentList* altlist = 
      fGeomVolSelector->GenerateTrimmedList(fCurrPathSegmentList);
    std::swap(altlist,fCurrPathSegmentList);
    delete altlist;  // after swap delete original
  }

  fCurrPathSegmentList->FillMatStepSum();

#ifdef RWH_DEBUG_2
  if ( fGeomVolSelector) { 
    // after FillMatStepSum() so one can see the summed mass
    if ( ( fDebugFlags & 0x40 ) ) {
      fCurrPathSegmentList->SetPrintVerbose(true);
      LOG("GROOTGeom", pNOTICE) << "After  trimming" << *fCurrPathSegmentList;
      fCurrPathSegmentList->SetPrintVerbose(false);
    }
  }
#endif


  return;
}

//___________________________________________________________________________
bool ROOTGeomAnalyzer::FindMaterialInCurrentVol(int tgtpdg)
{
  TGeoVolume * vol = fGeometry -> GetCurrentVolume();
  if(vol) {
    TGeoMaterial * mat = vol->GetMedium()->GetMaterial();
    if(mat->IsMixture()) {
      TGeoMixture * mixt = dynamic_cast <TGeoMixture*> (mat);
      for(int i = 0; i < mixt->GetNelements(); i++) {
         int pdg = this->GetTargetPdgCode(mixt, i);
         if(tgtpdg == pdg) return true;
      }
    } else {
       int pdg = this->GetTargetPdgCode(mat);
       if(tgtpdg == pdg) return true;
    }
  } else {
     LOG("GROOTGeom", pWARN) << "Current volume is null!";
     return false;
  }
  return false;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::StepToNextBoundary(void)
{
  fGeometry->FindNextBoundary();
  double step=fGeometry->GetStep();
  return step;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::Step(void)
{
  fGeometry->Step();
  double step=fGeometry->GetStep();
  return step;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::StepUntilEntering(void)
{
  this->StepToNextBoundary();  // doesn't actually step, so don't include in sum
  double step = 0; // 

  while(!fGeometry->IsEntering()) {
    step += this->Step();
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__

  bool isen = fGeometry->IsEntering();
  bool isob = fGeometry->IsOnBoundary();

  LOG("GROOTGeom",pDEBUG)
      << "IsEntering = "     << utils::print::BoolAsYNString(isen)
      << ", IsOnBoundary = " << utils::print::BoolAsYNString(isob)
      << ", Step = " << step;
#endif

  return step;
}
//___________________________________________________________________________
bool ROOTGeomAnalyzer::WillNeverEnter(double step)
{
// If the neutrino trajectory would never enter the detector, then the
// TGeoManager::GetStep returns the maximum step (1E30).
// Compare surrent step with max step and figure out whether the particle
// would never enter the detector

  if (step > 9.99E29) {

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GROOTGeom", pINFO) << "Wow! Current step is dr = " << step;
     LOG("GROOTGeom", pINFO) << "This trajectory isn't entering the detector";
#endif
     return true;

  } else 
    return false;
}

//___________________________________________________________________________
