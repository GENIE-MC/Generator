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

double NHLDecayVolume::lunits = genie::units::mm; string NHLDecayVolume::lunitString = "mm";
double NHLDecayVolume::aunits = genie::units::rad;
double NHLDecayVolume::tunits = genie::units::ns; string NHLDecayVolume::tunitString = "ns";

double NHLDecayVolume::fSx = 0.0, NHLDecayVolume::fSy = 0.0, NHLDecayVolume::fSz = 0.0; // start point
double NHLDecayVolume::fPx = 0.0, NHLDecayVolume::fPy = 0.0, NHLDecayVolume::fPz = 0.0; // momentum direction
double NHLDecayVolume::fEx = 0.0, NHLDecayVolume::fEy = 0.0, NHLDecayVolume::fEz = 0.0; // entry point
double NHLDecayVolume::fXx = 0.0, NHLDecayVolume::fXy = 0.0, NHLDecayVolume::fXz = 0.0; // exit  point

double NHLDecayVolume::fSxROOT = 0.0, NHLDecayVolume::fSyROOT = 0.0, NHLDecayVolume::fSzROOT = 0.0;
double NHLDecayVolume::fExROOT = 0.0, NHLDecayVolume::fEyROOT = 0.0, NHLDecayVolume::fEzROOT = 0.0;
double NHLDecayVolume::fXxROOT = 0.0, NHLDecayVolume::fXyROOT = 0.0, NHLDecayVolume::fXzROOT = 0.0;

double NHLDecayVolume::fDx = 0.0, NHLDecayVolume::fDy = 0.0, NHLDecayVolume::fDz = 0.0; // decay point
double NHLDecayVolume::fOx = 0.0, NHLDecayVolume::fOy = 0.0, NHLDecayVolume::fOz = 0.0; // origin of DV
double NHLDecayVolume::fLx = 0.0, NHLDecayVolume::fLy = 0.0, NHLDecayVolume::fLz = 0.0; // DV bounding box lengths

double NHLDecayVolume::fDxROOT = 0.0, NHLDecayVolume::fDyROOT = 0.0, NHLDecayVolume::fDzROOT = 0.0;
double NHLDecayVolume::fOxROOT = 0.0, NHLDecayVolume::fOyROOT = 0.0, NHLDecayVolume::fOzROOT = 0.0;
double NHLDecayVolume::fLxROOT = 0.0, NHLDecayVolume::fLyROOT = 0.0, NHLDecayVolume::fLzROOT = 0.0;

double NHLDecayVolume::fAx = 0.0, NHLDecayVolume::fAy = 0.0, NHLDecayVolume::fAz = 0.0; // axis of rotation
double NHLDecayVolume::fAlpha = 0.0; // rotation angle;
double NHLDecayVolume::ft = 0.0; // elapsed time

double NHLDecayVolume::kNewSpeedOfLight = genie::units::kSpeedOfLight 
  * (genie::units::m / lunits)
  / (genie::units::s / tunits);

//____________________________________________________________________________
void NHLDecayVolume::EnforceUnits( std::string length_units, std::string angle_units, std::string time_units ){
  
  LOG( "NHL", pDEBUG )
    << "Switching units to " << length_units.c_str() << " , " << angle_units.c_str() << " , " << time_units.c_str();

  double old_lunits = lunits;
  double old_aunits = aunits;
  double old_tunits = tunits;

  lunits = utils::units::UnitFromString( length_units ); lunitString = length_units;
  aunits = utils::units::UnitFromString( angle_units );
  tunits = utils::units::UnitFromString( time_units ); tunitString = time_units;

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
    << "kNewSpeedOfLight = " << kNewSpeedOfLight << " [" << lunitString.c_str() << "/"
    << tunitString.c_str() << "]";
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
  double maxRestTime = maxLabTime / gamma ; // this is how "wide" the detector looks like

  // if P(DL=0) = 1, P(DL = LMax) = exp( - LMax / c * 1/( beta * gamma ) * 1 / CoMLifetime )
  double PExit = std::exp( - maxRestTime / CoMLifetime );

  // from [0,1] we'd reroll anything in [0, PExit] and keep (PExit, 1]. That's expensive.
  // Instead, let 1 ==> maxRestTime, 0 ==> 0, exponential decay
  
  RandomGen * rnd = RandomGen::Instance();
  double ranthrow = rnd->RndGen().Uniform();

  double S0 = (1.0 - PExit) * ranthrow + PExit; 
  double rest_time = CoMLifetime * std::log( 1.0 / S0 );
  double elapsed_time = rest_time * gamma;
  double elapsed_length = elapsed_time * betaMag * kNewSpeedOfLight;

  LOG( "NHL", pDEBUG )
    << "betaMag, maxLength, CoMLifetime = " << betaMag << ", " << maxLength << ", " << CoMLifetime
    << "\nbetaMag = " << betaMag << " ==> gamma = " << gamma
    << "\n==> maxLength [" << tunitString.c_str()
    << "] = " << maxRestTime << " (rest frame) = " << maxLabTime << " (lab frame)"
    << "\nranthrow = " << ranthrow << ", PExit = " << PExit
    << "\n==> S0 = " << S0 << " ==> rest_time [" << lunitString.c_str() << "] = " << rest_time
    << " ==> elapsed_time [" << tunitString.c_str()
    << "] = " << elapsed_time << " ==> elapsed_length [" << lunitString.c_str()
    << "] = " << elapsed_length;

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

  fDxROOT = fDx * lunits / units::cm;
  fDyROOT = fDy * lunits / units::cm;
  fDzROOT = fDz * lunits / units::cm;

  LOG( "NHL", pDEBUG )
    << "\ndecayPoint = ( " << dx << ", " << dy << ", " << dz << " ) ["
    << lunitString.c_str() << "]"
    << "\ndecayPoint(ROOT) = ( " << fDxROOT << ", " << fDyROOT << ", " << fDzROOT << " ) [cm]";

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
  fLx = 1.0 * units::m; fLy = 1.0 * units::m; fLz = 1.0 * units::m;

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
  assert( fOx == 0.0 && fOy == 0.0 && fOz == 0.0 && 
	  fLx == 1.0 && fLy == 1.0 && fLz == 1.0 ); // SDV
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
  fLx = 2.0 * box->GetDX() * units::cm / lunits;
  fLy = 2.0 * box->GetDY() * units::cm / lunits;
  fLz = 2.0 * box->GetDZ() * units::cm / lunits;
  fOx = (box->GetOrigin())[0] * units::cm / lunits;
  fOy = (box->GetOrigin())[1] * units::cm / lunits;
  fOz = (box->GetOrigin())[2] * units::cm / lunits;

  fLxROOT = 2.0 * box->GetDX();
  fLyROOT = 2.0 * box->GetDY();
  fLzROOT = 2.0 * box->GetDZ();

  fOxROOT = (box->GetOrigin())[0];
  fOyROOT = (box->GetOrigin())[1];
  fOzROOT = (box->GetOrigin())[2];

  LOG( "NHL", pDEBUG )
    << "\nImported bounding box with origin at ( " << fOx << ", " << fOy << ", " << fOz << " ) and sides " << fLx << " x " << fLy << " x " << fLz << " [units: " << lunitString.c_str() << "]"
    << "\nIn ROOT units this is origin at ( " << fOxROOT << ", " << fOyROOT << ", " << fOzROOT << " ) and sides " << fLxROOT << " x " << fLyROOT << " x " << fLzROOT << " [cm]";
}
//____________________________________________________________________________
bool NHLDecayVolume::VolumeEntryAndExitPoints( TVector3 & startPoint, TVector3 & momentum,
					       TVector3 & entryPoint, TVector3 & exitPoint,
					       TGeoManager * gm, TGeoVolume * /* vol */ )
{
  const double mmtolunits = units::mm / lunits;

  double sx = startPoint.X(); double sy = startPoint.Y(); double sz = startPoint.Z();
  sx *= mmtolunits; sy *= mmtolunits; sz *= mmtolunits;
  double px = momentum.X(); double py = momentum.Y(); double pz = momentum.Z();
  double p2 = px*px + py*py + pz*pz; double p = std::sqrt(p2);
  px *= 1./p; py *= 1./p; pz *= 1./p;

  fSx = sx; fSy = sy; fSz = sz;
  fSxROOT = fSx * lunits / units::cm; fSyROOT = fSy * lunits / units::cm; fSzROOT = fSz * lunits / units::cm;
  fPx = px; fPy = py; fPz = pz;

  // put first point slightly inside the bounding box
  double firstZOffset = -0.1; // m
  firstZOffset *= units::m / lunits;

  double firstZ = fOx - fLz/2.0 - firstZOffset;

  LOG( "NHL", pDEBUG )
    << "\nfirstZ = " << firstZ << " [" << lunitString.c_str() << "]";

  // now find which point the line would hit this z at
  double dz = firstZ - sz;
  double tz = dz / pz;
  double dx = tz * px;
  double dy = tz * py;
  double firstX = sx + dx;
  double firstY = sy + dy;

  // now we gotta return everything to cm for ROOT to work its magic.
  double firstXROOT = firstX * lunits / units::cm,
    firstYROOT = firstY * lunits / units::cm, firstZROOT = firstZ * lunits / units::cm;

  gm->SetCurrentPoint( firstXROOT, firstYROOT, firstZROOT );
  gm->SetCurrentDirection( px, py, pz );

  LOG( "NHL", pINFO )
    << "\nCurrent point     is: ( " << firstX << ", " << firstY << ", " << firstZ << " ) [" << lunitString.c_str() << "]"
    << "\nFrom start point    : ( " << sx << ", " << sy << ", " << sz << " ) [" << lunitString.c_str() << "]"
    << "\nIn ROOT, current is : ( " << firstXROOT << ", " << firstYROOT << ", " << firstZROOT << " ) [cm]"
    << "\nIn ROOT, start is   : ( " << fSxROOT << ", " << fSxROOT << ", " << fSzROOT << " ) [cm]"
    << "\nCurrent direction is: ( " << px << ", " << py << ", " << pz << " ) [GeV/GeV]";

  assert( gm->FindNode() == NULL || gm->FindNode() == gm->GetTopNode() ); // need to be outside volume!

  double stepmax = 1.0e+6; // cm 
  stepmax *= genie::units::cm / lunits;

  // int ibound = 0;
  // const int imax = 10;

  LOG( "NHL", pDEBUG )
    << "Starting to search for intersections...";
  
  // enter the volume.
  /*
  // Do 10 cm steps. Once we hit the detector, reverse that step. Then we can use ROOT's geom stuff.
  double stepOffset = 1.0; // cm
  stepOffset *= units::cm / lunits; gm->SetStep( stepOffset );
  int iOffset = 0;
  TGeoVolume * gvol = gm->GetCurrentNode()->GetVolume();
  while( gvol->IsTopVolume() && iOffset < controls::kRjMaxIterations ){
    gm->Step( stepOffset );
    gvol = gm->GetCurrentNode()->GetVolume();
    iOffset++;
  }
  LOG( "NHL", pDEBUG )
    << "\nReached detector element at ( " << (gm->GetCurrentPoint())[0] << ", "
    << (gm->GetCurrentPoint())[1] << ", " << (gm->GetCurrentPoint())[2] << " ) after "
    << iOffset << " steps. Going back one...";

  gm->SetCurrentDirection( -px, -py, -pz ); gm->Step(); gm->SetCurrentDirection( px, py, pz );
  LOG( "NHL", pDEBUG )
    << "\nCurrent point     is: ( " << (gm->GetCurrentPoint())[0] << ", "
    << (gm->GetCurrentPoint())[1] << ", " << (gm->GetCurrentPoint())[2] 
    << " ) [" << lunitString.c_str() << "]"
    << "\nFrom start point    : ( " << sx << ", " << sy << ", " << sz << " ) [" << lunitString.c_str() << "]"
    << "\nCurrent direction is: ( " << px << ", " << py << ", " << pz << " ) [GeV/GeV]";
  */

  TGeoNode * nextNode = gm->FindNextBoundaryAndStep( stepmax );

  if( nextNode == NULL ) return false;

  // entered the detector, let's save this point
  fEx = ( gm->GetCurrentPoint() )[0] * genie::units::cm / lunits;
  fEy = ( gm->GetCurrentPoint() )[1] * genie::units::cm / lunits;
  fEz = ( gm->GetCurrentPoint() )[2] * genie::units::cm / lunits;
  entryPoint.SetXYZ( fEx, fEy, fEz );

  fExROOT = ( gm->GetCurrentPoint() )[0];
  fEyROOT = ( gm->GetCurrentPoint() )[1];
  fEzROOT = ( gm->GetCurrentPoint() )[2];

  LOG( "NHL", pDEBUG )
    << "\nEntry point found at ( " << fEx << ", " << fEy << ", " << fEz << " ) [" << lunitString.c_str() << "]"
    << "\nIn ROOT, entry at    ( " << fExROOT << ", " << fEyROOT << ", " << fEzROOT << " ) [cm]"; 

  // now propagate until we exit again
  
  int bdIdx = 0;
  const int bdIdxMax = 1e+4;

  double sfx = 0.0, sfy = 0.0, sfz = 0.0; // coords of the "safe" points in user units
  double sfxROOT = 0.0, sfyROOT = 0.0, sfzROOT = 0.0; // same, in cm

  // do one big step first
  // then if not outside yet, step by ever smaller steps until some threshold
  //Double_t sNext = std::max( fLx, std::max( fLy, fLz ) ) / 2.0;
  Double_t sNext = std::min( std::max( fLx, std::max( fLy, fLz ) ), 100.0 * lunits / units::cm ) / 2.0;
  Double_t sNextROOT = sNext * lunits / units::cm;
  gm->SetStep( sNextROOT );
  LOG( "NHL", pINFO )
    << "fLx, fLy, fLz = " << fLx << ", " << fLy << ", " << fLz << " ==> sNextROOT = " << sNextROOT;
  gm->Step();
  
  // FindNextBoundaryAndStep() sets step size to distance to next boundary and executes that step
  // so one "step" here is actually one big step + one small step
  while( gm->FindNextBoundaryAndStep() && bdIdx < bdIdxMax ){
    const Double_t * currPoint = gm->GetCurrentPoint();
    if( bdIdx % 100 == 0 ){
      LOG( "NHL", pDEBUG )
	<< "Step " << bdIdx << " : ( " << currPoint[0] << ", " << currPoint[1] << ", " << currPoint[2] << " ) [cm]";
    }
    sfxROOT = currPoint[0]; sfyROOT = currPoint[1]; sfzROOT = currPoint[2];
    sNextROOT *= 0.5;
    gm->SetStep( sNextROOT );
    gm->Step();
    bdIdx++;
  }
  if( bdIdx == bdIdxMax ){
    LOG( "NHL", pWARN )
      << "Failed to exit this volume. Dropping this trajectory.";
    return false;
  }

  // Always step back one step
  /*
  const Double_t * ffPoint = gm->GetCurrentPoint();
  if( std::abs(ffPoint[0] - fOxROOT) > fLxROOT/2.0 || 
      std::abs(ffPoint[1] - fOyROOT) > fLyROOT/2.0 || 
      std::abs(ffPoint[2] - fOzROOT) > fLzROOT/2.0 ){
    LOG( "NHL", pDEBUG )
      << "Overstepped bounding box: we're at ( " << ffPoint[0] << ", " << ffPoint[1] << ", " << ffPoint[2] << " ) [cm]";
    const Double_t * sfDir = gm->GetCurrentDirection();
    gm->SetCurrentDirection( -sfDir[0], -sfDir[1], -sfDir[2] );
    __attribute__((unused)) TGeoNode * tmpNode = gm->FindNextBoundaryAndStep();
    LOG( "NHL", pDEBUG )
      << "We turned back with new step = " << gm->GetStep();
    // and set direction back to normal
    gm->SetCurrentDirection( sfDir[0], sfDir[1], sfDir[2] );
  }
  const Double_t * sfPoint = gm->GetCurrentPoint();
  sfxROOT = sfPoint[0]; sfyROOT = sfPoint[1]; sfzROOT = sfPoint[2];
  */
  sfx = sfxROOT * units::cm / lunits; sfy = sfyROOT * units::cm / lunits; sfz = sfzROOT * units::cm / lunits;

  // exited the detector, let's save this point
  fXx = sfx; fXxROOT = sfxROOT;
  fXy = sfy; fXyROOT = sfyROOT;
  fXz = sfz; fXzROOT = sfzROOT;
  exitPoint.SetXYZ( fXx, fXy, fXz );

  LOG( "NHL", pINFO )
    << "\nExit point found at ( " << fXx << ", " << fXy << ", " << fXz << " ) ["
    << lunitString.c_str() << "]"
    << "\nIn ROOT, exit at    ( " << fXxROOT << ", " << fXyROOT << ", " << fXzROOT << " ) [cm]"; 

  return true;
  
}
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__
//____________________________________________________________________________
