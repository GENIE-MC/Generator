//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include "Physics/NeutralHeavyLepton/NHLDecayVolume.h"

using namespace genie;
using namespace genie::NHL;
using namespace genie::units;

double NHLDecayVolume::lunits = genie::units::mm;
double NHLDecayVolume::aunits = genie::units::rad;
double NHLDecayVolume::tunits = genie::units::ns;

double NHLDecayVolume::fSx = 0.0, NHLDecayVolume::fSy = 0.0, NHLDecayVolume::fSz = 0.0; // start point
double NHLDecayVolume::fPx = 0.0, NHLDecayVolume::fPy = 0.0, NHLDecayVolume::fPz = 0.0; // momentum direction
double NHLDecayVolume::fEx = 0.0, NHLDecayVolume::fEy = 0.0, NHLDecayVolume::fEz = 0.0; // entry point
double NHLDecayVolume::fXx = 0.0, NHLDecayVolume::fXy = 0.0, NHLDecayVolume::fXz = 0.0; // exit  point

double NHLDecayVolume::fDx = 0.0, NHLDecayVolume::fDy = 0.0, NHLDecayVolume::fDz = 0.0; // decay point
double NHLDecayVolume::fOx = 0.0, NHLDecayVolume::fOy = 0.0, NHLDecayVolume::fOz = 0.0; // origin of DV
double NHLDecayVolume::fLx = 0.0, NHLDecayVolume::fLy = 0.0, NHLDecayVolume::fLz = 0.0; // DV bounding box lengths

double NHLDecayVolume::fAx = 0.0, NHLDecayVolume::fAy = 0.0, NHLDecayVolume::fAz = 0.0; // axis of rotation
double NHLDecayVolume::fAlpha = 0.0; // rotation angle;
double NHLDecayVolume::ft = 0.0; // elapsed time

double NHLDecayVolume::kNewSpeedOfLight = genie::units::kSpeedOfLight 
  * (genie::units::m / genie::units::mm)
  / (genie::units::s / genie::units::ns);

//____________________________________________________________________________
void NHLDecayVolume::EnforceUnits( std::string length_units, std::string angle_units, std::string time_units ){
  
  LOG( "NHL", pDEBUG )
    << "Switching units to " << length_units.c_str() << " , " << angle_units.c_str() << " , " << time_units.c_str();

  double old_lunits = lunits;
  double old_aunits = aunits;
  double old_tunits = tunits;

  lunits = utils::units::UnitFromString( length_units );
  aunits = utils::units::UnitFromString( angle_units );
  tunits = utils::units::UnitFromString( time_units );

  // convert to new units
  fSx /= lunits/old_lunits; fSy /= lunits/old_lunits; fSz /= lunits/old_lunits;
  fPx /= lunits/old_lunits; fPy /= lunits/old_lunits; fPz /= lunits/old_lunits;
  fEx /= lunits/old_lunits; fEy /= lunits/old_lunits; fEz /= lunits/old_lunits;
  fXx /= lunits/old_lunits; fXy /= lunits/old_lunits; fXz /= lunits/old_lunits;

  fDx /= lunits/old_lunits; fDy /= lunits/old_lunits; fDz /= lunits/old_lunits;
  fOx /= lunits/old_lunits; fOy /= lunits/old_lunits; fOz /= lunits/old_lunits;

  fAx /= lunits/old_lunits; fAy /= lunits/old_lunits; fAz /= lunits/old_lunits;
  fAlpha /= aunits/old_aunits;
  ft /= tunits/old_tunits;

  kNewSpeedOfLight /= (lunits / old_lunits) / (tunits / old_tunits);

  LOG( "NHL", pDEBUG )
    << "kNewSpeedOfLight = " << kNewSpeedOfLight << " [mm/ns]";
}
//____________________________________________________________________________
double NHLDecayVolume::CalcTravelLength( double betaMag, double CoMLifetime, double maxLength )
{
  // decay probability P0(t) = 1 - exp( -t/tau ) where:
  // t   = time-of-flight (in rest frame)
  // tau = CoMLifetime

  double maxLabTime = maxLength / kNewSpeedOfLight;
  assert( betaMag > 0.0 && betaMag < 1.0 ); // massive moving particle
  double gamma = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
  double maxRestTime = maxLabTime / ( betaMag * gamma ); // this is how "wide" the detector looks like

  // if P(DL=0) = 1, P(DL = LMax) = exp( - LMax / c * 1/( beta * gamma ) * 1 / CoMLifetime )
  double PExit = std::exp( - maxRestTime / CoMLifetime );

  // from [0,1] we'd reroll anything in [0, PExit] and keep (PExit, 1]. That's expensive.
  // Instead, let 1 ==> maxRestTime, 0 ==> 0, exponential decay
  
  RandomGen * rnd = RandomGen::Instance();
  double ranthrow = rnd->RndGen().Uniform();

  double S0 = (1.0 - PExit) * ranthrow + PExit; 
  double elapsed_time = CoMLifetime * std::log( 1.0 / S0 );
  double elapsed_length = elapsed_time * kNewSpeedOfLight;

  LOG( "NHL", pDEBUG )
    << "betaMag, maxLength, CoMLifetime = " << betaMag << ", " << maxLength << ", " << CoMLifetime
    << "\nbetaMag = " << betaMag << " ==> gamma = " << gamma
    << "\n==> maxLength = " << maxRestTime << " (rest frame) = " << maxLabTime << " (lab frame)"
    << "\nranthrow = " << ranthrow << ", PExit = " << PExit
    << "\n==> S0 = " << S0 << " ==> elapsed_time = " << elapsed_time << " ==> elapsed_length = " << elapsed_length;

  return elapsed_length;
}
//____________________________________________________________________________
TVector3 NHLDecayVolume::GetDecayPoint( double travelLength, TVector3 & entryPoint, TVector3 & momentum )
{
  double ex = entryPoint.X(); double ey = entryPoint.Y(); double ez = entryPoint.Z();
  double px = momentum.X(); double py = momentum.Y(); double pz = momentum.Z();
  double p2 = px*px + py*py + pz*pz; double p = std::sqrt(p2);
  px *= 1./p; py *= 1./p; pz *= 1./p;

  double dx = ex + travelLength * px; fDx = dx;
  double dy = ey + travelLength * py; fDy = dy;
  double dz = ez + travelLength * pz; fDz = dz;

  LOG( "NHL", pDEBUG )
    << "decayPoint = (" << dx << ", " << dy << ", " << dz << ")";

  TVector3 decayPoint( dx, dy, dz );
  return decayPoint;
}
//____________________________________________________________________________
double NHLDecayVolume::GetMaxLength( TVector3 & entryPoint, TVector3 & exitPoint )
{
  double ex = entryPoint.X(); double ey = entryPoint.Y(); double ez = entryPoint.Z();
  double xx = exitPoint.X(); double xy = exitPoint.Y(); double xz = exitPoint.Z();

  return std::sqrt( (ex-xx)*(ex-xx) + (ey-xy)*(ey-xy) + (ez-xz)*(ez-xz) );
}
//____________________________________________________________________________
void NHLDecayVolume::MakeSDV()
{
  fOx = 0.0; fOy = 0.0; fOz = 0.0;
  fLx = 1.0; fLy = 1.0; fLz = 1.0;

  LOG("NHL", pDEBUG)
    << "Setting simple decay volume with unit-m side."
    << "\nSetting units to \"m\", \"rad\", \"ns\"";

  EnforceUnits("m","rad","ns");
}
//____________________________________________________________________________
// if entry and exit points, populate TVector3's with their coords. If not, return false
bool NHLDecayVolume::SDVEntryAndExitPoints( TVector3 & startPoint, TVector3 & momentum,
					    TVector3 & entryPoint, TVector3 & exitPoint )
{
  assert( fOx == 0.0 && fOy == 0.0 && fOz == 0.0 && fLx == 1.0 && fLy == 0.0 && fLz == 0.0 ); // SDV
  fSx = startPoint.X(); fSy = startPoint.Y(); fSz = startPoint.Z();
  fPx = momentum.X(); fPy = momentum.Y(); fPz = momentum.Z();
  double fP2 = fPx*fPx + fPy*fPy + fPz*fPz; double fP = std::sqrt(fP2);
  fPx *= 1.0/fP; fPy *= 1.0/fP; fPz *= 1.0/fP;

  // calc parameter for line at each face
  double txP = (  fLx - fSx ) / fPx;
  double txM = ( -fLx - fSx ) / fPx;
  double tyP = (  fLy - fSy ) / fPy;
  double tyM = ( -fLy - fSy ) / fPy;
  double tzP = (  fLz - fSz ) / fPz;
  double tzM = ( -fLz - fSz ) / fPz;

  // do we have an entry or exit anywhere?
  // entry from face Q = const <==> {  pr(momentum, Q) points to origin && within bounding square }
  // exit  from face Q = const <==> { -pr(momentum, Q) points to origin && within bounding square }
  double q1t = 0.0, q2t = 0.0;
  bool pointsOnX = false, pointsOnY = false, pointsOnZ = false;

  // case x = +fLx
  q1t = fSy + txP * fPy; q2t = fSz + txP * fPz;
  if( std::abs( q1t ) <= fLy && std::abs( q2t ) <= fLz ){ // within bounding square
    pointsOnX = true;
    if( fSx * fPx < 0 ){ // pointing towards origin
      fEx = fLx; fEy = q1t; fEz = q2t;
    } else if( fSx * fPx > 0 ){ // pointing away from origin
      fXx = fLx; fXy = q1t; fXz = q2t;
    } else return false; // treat tangent as no entry
  }
  // case x = -fLx
  q1t = fSy + txM * fPy; q2t = fSz + txM * fPz;
  if( std::abs( q1t ) <= fLy && std::abs( q2t ) <= fLz ){ // within bounding square
    pointsOnX = true;
    if( fSx * fPx < 0 ){ // pointing towards origin
      fEx = -fLx; fEy = q1t; fEz = q2t;
    } else if( fSx * fPx > 0 ){ // pointing away from origin
      fXx = -fLx; fXy = q1t; fXz = q2t;
    } else return false; // treat tangent as no entry
  }

  // case y = +fLy
  q1t = fSz + tyP * fPz; q2t = fSx + tyP * fPx;
  if( std::abs( q1t ) <= fLz && std::abs( q2t ) <= fLx ){ // within bounding square
    pointsOnY = true;
    if( fSy * fPy < 0 ){ // pointing towards origin
      fEx = q2t; fEy = fLy; fEz = q1t;
    } else if( fSy * fPy > 0 ){ // pointing away from origin
      fXx = q2t; fXy = fLy; fXz = q1t;
    } else return false; // treat tangent as no entry
  }
  // case y = -fLy
  q1t = fSz + tyM * fPz; q2t = fSx + tyM * fPx;
  if( std::abs( q1t ) <= fLz && std::abs( q2t ) <= fLx ){ // within bounding square
    pointsOnY = true;
    if( fSy * fPy < 0 ){ // pointing towards origin
      fEx = q2t; fEy = -fLy; fEz = q1t;
    } else if( fSy * fPy > 0 ){ // pointing away from origin
      fXx = q2t; fXy = -fLy; fXz = q1t;
    } else return false; // treat tangent as no entry
  }

  // case z = +fLz
  q1t = fSx + tzP * fPx; q2t = fSy + tzP * fPy;
  if( std::abs( q1t ) <= fLx && std::abs( q2t ) <= fLy ){ // within bounding square
    pointsOnZ = true;
    if( fSz * fPz < 0 ){ // pointing towards origin
      fEx = q1t; fEy = q2t; fEz = fLz;
    } else if( fSz * fPz > 0 ){ // pointing away from origin
      fXx = q1t; fXy = q2t; fXz = fLz;
    } else return false; // treat tangent as no entry
  }
  // case z = -fLz
  q1t = fSx + tzM * fPx; q2t = fSy + tzM * fPy;
  if( std::abs( q1t ) <= fLx && std::abs( q2t ) <= fLy ){ // within bounding square
    pointsOnZ = true;
    if( fSz * fPz < 0 ){ // pointing towards origin
      fEx = q1t; fEy = q2t; fEz = -fLz;
    } else if( fSz * fPz > 0 ){ // pointing away from origin
      fXx = q1t; fXy = q2t; fXz = -fLz;
    } else return false; // treat tangent as no entry
  }

  bool finalPoints = ( pointsOnX || pointsOnY || pointsOnZ );
  if( finalPoints ){
    entryPoint.SetXYZ( fEx, fEy, fEz );
    exitPoint.SetXYZ( fXx, fXy, fXz );
    return true;
  }
  
  // missed detector
  return false;
}
//____________________________________________________________________________
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
void NHLDecayVolume::ImportBoundingBox( TGeoBBox * box )
{
  fLx = 2.0 * box->GetDX();
  fLy = 2.0 * box->GetDY();
  fLz = 2.0 * box->GetDZ();
  fOx = (box->GetOrigin())[0];
  fOy = (box->GetOrigin())[1];
  fOz = (box->GetOrigin())[2];

  LOG( "NHL", pDEBUG )
    << "Imported bounding box with origin at ( " << fOx << ", " << fOy << ", " << fOz << " ) and sides " << fLx << " x " << fLy << " x " << fLz;
}
//____________________________________________________________________________
bool NHLDecayVolume::VolumeEntryAndExitPoints( TVector3 & startPoint, TVector3 & momentum,
					       TVector3 & entryPoint, TVector3 & exitPoint,
					       TGeoManager * gm, TGeoVolume * vol )
{
  double sx = startPoint.X(); double sy = startPoint.Y(); double sz = startPoint.Z();
  double px = momentum.X(); double py = momentum.Y(); double pz = momentum.Z();
  double p2 = px*px + py*py + pz*pz; double p = std::sqrt(p2);
  px *= 1./p; py *= 1./p; pz *= 1./p;

  fSx = sx; fSy = sy; fSz = sz;
  fPx = px; fPy = py; fPz = pz;

  // RETHERE I am hacking to set MINERvA ID tracker as top volume.
  /*
  TGeoVolume * tracker = gm->FindVolumeFast("DetectorlvTracker");
  assert(tracker);
  gm->SetTopVolume(tracker);
  */

  // put first point at z = const 0.5m behind the bounding box
  // RETHERE, this is a placeholder
  double firstZOffset = 0.1; // m
  firstZOffset *= genie::units::m / lunits;

  double firstZ = fOx - fLz/2.0 - firstZOffset;

  LOG( "NHL", pDEBUG )
    << "firstZ = " << firstZ;

  // now find which point the line would hit this z at
  double dz = firstZ - sz;
  double tz = dz / pz;
  double dx = tz * px;
  double dy = tz * py;
  double firstX = sx + dx;
  double firstY = sy + dy;

  gm->SetCurrentPoint( firstX, firstY, firstZ );
  gm->SetCurrentDirection( px, py, pz );

  LOG( "NHL", pDEBUG )
    << "\nCurrent point     is: ( " << firstX << ", " << firstY << ", " << firstZ << " )"
    << "\nCurrent direction is: ( " << px << ", " << py << ", " << pz << " )";

  assert( gm->FindNode() == NULL || gm->FindNode() == gm->GetTopNode() ); // need to be outside volume!

  // RETHERE - the root file units are cm! Gotta implement or use -L option
  double stepmax = 1.0e+6; // cm 
  stepmax *= genie::units::cm / lunits;
  double stepreg = 1.0; // mm
  stepreg *= genie::units::mm / lunits;
  double stepsmall = 1.0e-2; // cm
  stepsmall *= genie::units::cm / lunits;

  int ibound = 0;
  const int imax = 10;

  LOG( "NHL", pDEBUG )
    << "Starting to search for intersections...";
  
  // enter the volume.
  TGeoNode * nextNode = gm->FindNextBoundaryAndStep( stepmax ); 

  if( nextNode == NULL ) return false;

  // entered the detector, let's save this point
  fEx = ( gm->GetCurrentPoint() )[0] * genie::units::cm / genie::units::mm; // RETHERE fix this conversion!
  fEy = ( gm->GetCurrentPoint() )[1] * genie::units::cm / genie::units::mm;
  fEz = ( gm->GetCurrentPoint() )[2] * genie::units::cm / genie::units::mm;
  entryPoint.SetXYZ( fEx, fEy, fEz );

  LOG( "NHL", pDEBUG )
    << "Entry point found at ( " << fEx << ", " << fEy << ", " << fEz << " ) [mm]"; 

  // now propagate until we exit again
  
  int bdIdx = 0;
  const int bdIdxMax = 1e+4;

  double sfx = 0.0, sfy = 0.0, sfz = 0.0; // coords of the "safe" points

  // do one big step first, half of largest BBox dimension
  // then if not outside yet, step by ever smaller steps until some threshold
  Double_t sNext = std::max( fLx, std::max( fLy, fLz ) ) / 2.0;
  gm->SetStep( sNext );
  LOG( "NHL", pINFO )
    << "fLx, fLy, fLz = " << fLx << ", " << fLy << ", " << fLz << " ==> sNext = " << sNext;
  gm->Step();
  
  // FindNextBoundaryAndStep() sets step size to distance to next boundary and executes that step
  // so one "step" here is actually one big step + one small step
  while( gm->FindNextBoundaryAndStep() && bdIdx < bdIdxMax ){
    const Double_t * currPoint = gm->GetCurrentPoint();
    LOG( "NHL", pINFO )
      << "Step " << bdIdx << " : ( " << currPoint[0] << ", " << currPoint[1] << ", " << currPoint[2] << " )";
    sNext *= 0.5;
    gm->SetStep( sNext );
    gm->Step();
    bdIdx++;
  }
  if( bdIdx == bdIdxMax ){
    LOG( "NHL", pWARN )
      << "Failed to exit this volume. Dropping this trajectory.";
    return false;
  }

  // check to see if we're outside the bounding box.
  // If yes, we've overstepped, return to detector!
  const Double_t * ffPoint = gm->GetCurrentPoint();
  if( std::abs(ffPoint[0] - fOx) > fLx/2.0 || std::abs(ffPoint[1] - fOy) > fLy/2.0 || std::abs(ffPoint[2] - fOz) > fLz/2.0 ){
    LOG( "NHL", pDEBUG )
      << "Overstepped bounding box: we're at ( " << ffPoint[0] << ", " << ffPoint[1] << ", " << ffPoint[2] << " )";
    const Double_t * sfDir = gm->GetCurrentDirection();
    gm->SetCurrentDirection( -sfDir[0], -sfDir[1], -sfDir[2] );
    TGeoNode * tmpNode = gm->FindNextBoundaryAndStep();
    LOG( "NHL", pDEBUG )
      << "We turned back with new step = " << gm->GetStep();
    // and set direction back to normal
    gm->SetCurrentDirection( sfDir[0], sfDir[1], sfDir[2] );
  }
  const Double_t * sfPoint = gm->GetCurrentPoint();
  sfx = sfPoint[0]; sfy = sfPoint[1]; sfz = sfPoint[2];

  // exited the detector, let's save this point
  fXx = sfx * genie::units::cm / genie::units::mm; // RETHERE fix this conversion!
  fXy = sfy * genie::units::cm / genie::units::mm;
  fXz = sfz * genie::units::cm / genie::units::mm;
  exitPoint.SetXYZ( fXx, fXy, fXz );

  LOG( "NHL", pDEBUG )
    << "Exit point found at ( " << fXx << ", " << fXy << ", " << fXz << " ) [mm]"; 

  return true;
  
}
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__
//____________________________________________________________________________
