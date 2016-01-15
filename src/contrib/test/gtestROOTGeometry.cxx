//____________________________________________________________________________
/*!

\program gtestROOTGeometry

\brief   Tests the GENIE ROOT geometry driver by generating random rays,
         following them through the detector, computing density-weighted
         path-lengths for each material & generating vertices

\syntax  gtestROOTGeometry [-f geom] [-n nvtx] [-d dx,dy,dz] [-s x,y,z] 
                           [-r size] [-p pdg] 
  
         Options:

          -f  A ROOT file containing a ROOT/GEANT geometry description
              [default: $GENIE/data/geo/samples/BoxWithLArPbLayers.root]
          -n  Number of vertices to generate
          -d  Direction of generated rays (dirx,diry,dirz) 
              [default: 1,0,0 - along x]
          -s  Specifies the ray generation surface.
              That surface is perpendicular to the beam direction and
              contains the point specified here.
              Use the same units as your input geometry.
              [default: 0,0,0]
          -r  Specifies the radius of the area over which rays will be
              generated. That event generation area is the area contained 
              within the circle centered at the point specified with the
              -s option and a radius specified here.
              Use the same units as your input geometry.
              [default: 100 in your geometry unit]
          -p  Can be used to specify a target pdg code. If that option is set
              then the vertex is generated only at that material. If not set
              then vertices will be generated in all detector materials (by
              weight).
          -v  Can specify specific volumes of the input geometry. If not set
              will use the master volume.

          Examples:
            
          gtestROOTGeometry -n 20000 -d 1,0,0 -s -10,0,0 -r 10 -p 1000180400

          will take a test geom ($GENIE/data/geo/samples/BoxWithLArPbLayers.root)
          and generate 20k vertices in Ar40 (pdg=1000180400) only, using rays 
          (flux) generated uniformly (in area) within a circle of radius = 10
          (in units of that input geometry, m in this case), centered at 
          (-10,0,0) (see prev comment on units) and perpendicular to the ray 
          direction (1,0,0).

\Author  Anselmo Meregaglia <anselmo.meregaglia@cern.ch>
         ETH Zurich

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Lab

\created August 11, 2005

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <vector>
#include <cassert>

#include <TFile.h>
#include <TNtupleD.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TApplication.h>
#include <TPolyMarker3D.h>

#include "Conventions/Constants.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "EVGDrivers/PathLengthList.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLnArgParser.h"

using std::string;
using std::vector;

using namespace genie;
using namespace genie::constants;

void GetCommandLineArgs (int argc, char ** argv);
void GetRandomRay       (TLorentzVector & x, TLorentzVector & p);
int  GetTargetMaterial  (const PathLengthList & pl);

string   gOptGeomFile;       // input ROOT geom file
string   gOptRootGeomTopVol; // top volume / can be used to override the master volume
TVector3 gOptRayDirection;   // ray direction
TVector3 gOptRaySurf;        // ray generation surface
double   gOptRayR;           // ray generation area radius
int      gOptNVtx;           // number of vertices to generate
int      gOptTgtPdg;         // 

double   kDefOptRayR         = 100;
TVector3 kDefOptRayDirection (1,0,0);
TVector3 kDefOptRaySurf      (0,0,0);

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);
 
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__

  TApplication theApp("App", &argc, argv);
 
  // Create & configure the geometry driver
  //
  LOG("test", pINFO) 
     << "Creating a geometry driver for ROOT geometry at: " 
     << gOptGeomFile;
  geometry::ROOTGeomAnalyzer * geom_driver = 
               new geometry::ROOTGeomAnalyzer(gOptGeomFile);
  geom_driver ->SetTopVolName(gOptRootGeomTopVol);

  // Draw the geometry
  // & define TPolyMarker3D for drawing vertices later on
  //
  LOG("test", pINFO) 
      << "Drawing the ROOT geometry";
  geom_driver->GetGeometry()->GetTopVolume()->Draw();

  TPolyMarker3D * vtxp = new TPolyMarker3D();
  vtxp->SetMarkerColor(kRed);
  vtxp->SetMarkerStyle(8);
  vtxp->SetMarkerSize(0.4);

  // Compute & Printout the (density weighted) max path lengths
  //
  LOG("test", pINFO) 
       << "Computing max {density-weighted path lengths}";
  const PathLengthList & maxpl = geom_driver->ComputeMaxPathLengths();

  LOG("test", pINFO) << "Maximum math lengths: " << maxpl;

  TFile f("geomtest.root","recreate");
  TNtupleD vtxnt("vtxnt","","x:y:z:A:Z");

  TLorentzVector x(0,0,0,0);
  TLorentzVector p(0,0,0,0);

  int n = 0;
  while (n < gOptNVtx) {

    // generate a random ray 
    GetRandomRay(x,p);

    // compute density-weighted path lengths for each geometry
    // material for the current ray
    const PathLengthList & pl = geom_driver->ComputePathLengths(x,p);
    LOG("test",pINFO)        
       << "Current path lengths: " << pl;

    // select detector material (amongst all materials defined in the 
    // detector geometry -- do so based on density-weighted path lengths) 
    // or force it to the user-selected material
    int tpdg = GetTargetMaterial(pl);
    if (tpdg == -1) continue;
    LOG("test",pINFO) << "Selected target material: " << tpdg;

    // generate an 'interaction vertex' in the selected material
    const TVector3 & vtx = geom_driver->GenerateVertex(x,p,tpdg);
    LOG("test",pINFO) 
      << "Generated vtx: (x = " << vtx.X() 
      << ", y = " << vtx.Y() << ", z = " <<vtx.Z() << ")";

    // add it at the ntuple & at the vtx marker
    vtxnt.Fill(vtx.X(),vtx.Y(),vtx.Z(), 
               pdg::IonPdgCodeToA(tpdg), pdg::IonPdgCodeToZ(tpdg));
    vtxp->SetNextPoint(vtx.X(),vtx.Y(),vtx.Z()); 

    n++;
    LOG("test", pNOTICE) 
      << " *** Vertices generated so far: " << n;
  }
 
  // draw vertices
  vtxp->Draw("same");
  
  vtxnt.Write();
  f.Close();
  
  theApp.Run(kTRUE);

#else
    LOG("test", pERROR) 
       << "*** You should have enabled the geometry drivers first!";
#endif

  return 0;
}
//____________________________________________________________________________
void GetRandomRay(TLorentzVector & x, TLorentzVector & p)
{
// generate a random ray (~flux neutrino)
//
  RandomGen * rnd = RandomGen::Instance();

  TVector3 vec0(gOptRayDirection);
  TVector3 vec = vec0.Orthogonal();

  double phi = 2.*kPi * rnd->RndFlux().Rndm();
  double Rt  = -1;
  bool accept = false;
  while(!accept) {
    double r = gOptRayR * rnd->RndFlux().Rndm();
    double y = gOptRayR * rnd->RndFlux().Rndm();
    if(y<r) {
      accept = true;
      Rt = r;
    }
  }

  vec.Rotate(phi,vec0);
  vec.SetMag(Rt);

  vec = vec + gOptRaySurf;

  TLorentzVector xx(vec, 0.);
  TLorentzVector pp(gOptRayDirection, gOptRayDirection.Mag());

  x = xx;
  p = pp;

  LOG("test", pNOTICE) 
   << "** Curr ray:";
  LOG("test", pNOTICE) 
   << "    x = " << x.X() << ",  y = " << x.Y() << ",  z = " << x.Z();
  LOG("test", pNOTICE) 
   << "    px = " << p.X() << ", py = " << p.Y() << ", pz = " << p.Z();

}
//____________________________________________________________________________
int GetTargetMaterial(const PathLengthList & pl)
{
  if(pl.AreAllZero()) return -1;

  if(gOptTgtPdg > 0) {
    if(pl.PathLength(gOptTgtPdg) > 0) return gOptTgtPdg;
  } 
  else {
    RandomGen * rnd = RandomGen::Instance();

    PathLengthList::const_iterator pliter;
    double sum = 0;
    for(pliter = pl.begin(); pliter != pl.end(); ++pliter) {
        sum += pliter->second;
    }
    double cpl = sum * rnd->RndFlux().Rndm();
    sum = 0;
    for(pliter = pl.begin(); pliter != pl.end(); ++pliter) {
       sum += pliter->second;
       if(cpl < sum) {
          return pliter->first;
       }
    }
  }
  return -1;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  CmdLnArgParser parser(argc,argv);

  // geometry file
  if( parser.OptionExists('f') ) {
    LOG("test", pINFO) << "Getting ROOT geometry filename";
    gOptGeomFile = parser.ArgAsString('f');
  } else {
    string base_dir = string( gSystem->Getenv("GENIE") );
    string filename = base_dir + 
          string("/data/geo/samples/BoxWithLArPbLayers.root");
    gOptGeomFile = filename;
  }

  // check whether an event generation volume name has been 
  // specified -- default is the 'top volume'
  if( parser.OptionExists('v') ) {
    LOG("test", pDEBUG) << "Checking for input volume name";
    gOptRootGeomTopVol = parser.ArgAsString('v');
  } else {
    LOG("test", pDEBUG) << "Using the <master volume>";
  } 

  // direction
  if( parser.OptionExists('d') ) {
    LOG("test", pINFO) << "Reading ray direction";
    // split the comma separated list
    vector<double> dirv = parser.ArgAsDoubleTokens('d', ",");
    assert(dirv.size() == 3);
    gOptRayDirection.SetXYZ(dirv[0],dirv[1],dirv[2]);
  } else {
    LOG("test", pINFO) << "No input ray direction - Using default";
    gOptRayDirection = kDefOptRayDirection;
  }

  // ray surface:
  if( parser.OptionExists('s') ) {
    LOG("test", pINFO) << "Reading ray generation surface";
    // split the comma separated list
    vector<double> rsv = parser.ArgAsDoubleTokens('s', ",");
    assert(rsv.size() == 3);
    gOptRaySurf.SetXYZ(rsv[0],rsv[1],rsv[2]);
  } else {
    LOG("test", pINFO) << "No input ray generation surface - Using default";
    gOptRaySurf = kDefOptRaySurf;
  }

  // ray generation area radius:
  if( parser.OptionExists('r') ) {
    LOG("test", pINFO) << "Reading radius of ray generation area";
    gOptRayR = parser.ArgAsDouble('r');
  } else {
    LOG("test", pINFO) << "No input radius of ray generation area - Using default";
    gOptRayR = kDefOptRayR;
  }
  gOptRayR = TMath::Abs(gOptRayR); // must be positive

  // number of vertices to generate:
  if( parser.OptionExists('n') ) {
    LOG("test", pINFO) << "Getting number of vertices to generate";
    gOptNVtx = parser.ArgAsInt('n');
  } else {
    gOptNVtx = 0;
  }

  // 'forced' target pdg:
  if( parser.OptionExists('p') ) {
    LOG("test", pINFO) << "Getting 'forced' target pdg";
    gOptTgtPdg = parser.ArgAsInt('p');
  } else {
    gOptTgtPdg = -1;
  }

  LOG("test", pNOTICE) 
    << "\n Options: "
    << "\n ROOT geometry file: " << gOptGeomFile
    << "\n ROOT geometry top volume: " << gOptRootGeomTopVol
    << "\n Ray direction: (" 
           << gOptRayDirection.X() << ", "
           << gOptRayDirection.Y() << ", "
           << gOptRayDirection.Z() << ") "
    << "\n Ray generation surface : (" 
           << gOptRaySurf.X() << ", "
           << gOptRaySurf.Y() << ", "
           << gOptRaySurf.Z() << ") "
    << "\n Ray generation area radius : " << gOptRayR 
    << "\n Number of vertices : " << gOptNVtx
    << "\n Forced targer PDG : "  << gOptTgtPdg;

}
//____________________________________________________________________________
