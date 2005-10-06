//____________________________________________________________________________
/*!

\class   genie::geometry::ROOTGeomAnalyzer

\brief   A ROOT/GEANT Geometry Analyzer

\author  Anselmo Meregaglia <anselmo.meregaglia@cern.ch>
         ETH Zurich

\created May 24, 2005

*/
//____________________________________________________________________________

#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TSystem.h>

#include "Conventions/Units.h"
#include "EVGDrivers/PathLengthList.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

#include <TPolyMarker3D.h>
#include <TGeoBBox.h>
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::geometry;

//___________________________________________________________________________
ROOTGeomAnalyzer::ROOTGeomAnalyzer(string filename) :
GeomAnalyzerI()
{
  this->Initialize(filename);
}
//___________________________________________________________________________
ROOTGeomAnalyzer::~ROOTGeomAnalyzer()
{
  if( fCurrPathLengthList    ) delete fCurrPathLengthList;
  if( fCurrMaxPathLengthList ) delete fCurrMaxPathLengthList;
  if( fCurrPDGCodeList       ) delete fCurrPDGCodeList;
  if( fGeometry              ) delete fGeometry;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::SetUnits(double u)
{
// Use the units of your input geometry, eg
//               geom.SetUnits(genie::units::centimeter)
// GENIE uses the physical system of units (hbar=c=1) almost throughtout so
// everything is expressed in GeV but when analyzing detector geometries we 
// use meters. Setting your input geometry units will allow us to figure the
// conversion factor. 
// As input, use one of the constants in $GENIE/src/Conventions/Units.h

  fScale = u/units::meter;
}
//___________________________________________________________________________
const PathLengthList & ROOTGeomAnalyzer::ComputeMaxPathLengths(void)
{
  LOG("GROOTGeom", pINFO)
                  << "Computing the maximum path lengths for all materials";
  if(!fGeometry)
    {
      LOG("GROOTGeom",pERROR) << "No ROOT geometry is loaded!";
      return *fCurrMaxPathLengthList;
    }

  fCurrMaxPathLengthList->SetAllToZero();

  //select World volume
  TObjArray *LV =0;
  LV=fGeometry->GetListOfVolumes();

  int numVol;
  numVol=(LV->GetEntries());

  TGeoVolume *TV =0;
  TGeoVolume *TVWorld =0;
  TGeoShape *TS=0;;

  char *name;
  char *str;
  str="World";
  int FlagFound(0);

  for(Int_t i=0;i<numVol;i++)
    {
      TV= dynamic_cast <TGeoVolume *> (LV->At(i));
      name=const_cast<char*>(TV->GetName());
      if(!strcmp(str,name))
        {
          FlagFound=1;
          TVWorld=TV;
          break;
        }
    }

  if(!FlagFound)
    {
      LOG("GROOTGeom",pERROR) << "The World Volume does not exist in your geometry!!! ";
      return *fCurrMaxPathLengthList;
    }

  LOG("GROOTGeom", pINFO) << "Getting a TGeoBBox enclosing the detector";
  TS=TVWorld->GetShape();
  TGeoBBox *box=(TGeoBBox *)TS;

  //get box origin and dimensions
  double dx = box->GetDX(); // half-length
  double dy = box->GetDY(); // half-length
  double dz = box->GetDZ(); // half-length
  double ox = (box->GetOrigin())[0];
  double oy = (box->GetOrigin())[1];
  double oz = (box->GetOrigin())[2];
  LOG("GROOTGeom",pINFO)
     << "Box dimensions : x = "<< 2*dx << ", y = "<< 2*dy << ", z = "<< 2*dz;
  LOG("GROOTGeom",pINFO)
     << "Box origin     : x = "<< ox   << ", y = "<< oy   << ", z = "<<   oz;

  //generate 200 random points on each surface, use 200 rays to
  //calculate maximum path for each material

  RandomGen* rand=RandomGen::Instance();
  TRandom & r3=rand->Random3();

  LOG("GROOTGeom",pINFO)
        << "Will generate [" << fNPoints << "] random points on each box surface";
  LOG("GROOTGeom",pINFO)
        << "Will generate [" << fNRays   << "] rays for each point";

  //loop on materials
  int pdgc(0);
  vector<int>::iterator itrPDG;
  for(itrPDG=fCurrPDGCodeList->begin();itrPDG!=fCurrPDGCodeList->end();itrPDG++)
    {
      pdgc=*itrPDG;
      LOG("GROOTGeom", pINFO)
             <<" Calculating max path length for material: " << pdgc;

      int igen(0);
      int Rgen(0);
      int maxPoints(fNPoints);
      int maxRays(fNRays);
      double xyz[3];
      double direction[3];
      double MaxPath(0);
      double Length(0);
      double dirTot(0);

      //top
      LOG("GROOTGeom",pINFO) << "Box surface scanned: [TOP]";
      igen=0;
      while(igen<maxPoints)
        {
          igen++;
          xyz[0] = ox-dx+2*dx*r3.Rndm();
          xyz[1] = oy+dy;
          xyz[2] = oz-dz+2*dz*r3.Rndm();

          Rgen=0;
          while(Rgen<maxRays)
            {
              Rgen++;
              direction[0]=-0.5+r3.Rndm();
              direction[1]=-r3.Rndm();
              direction[2]=-0.5+r3.Rndm();

              dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

              direction[0]/=dirTot;
              direction[1]/=dirTot;
              direction[2]/=dirTot;

              Length=ComputeMaxPathLengthPDG(xyz,direction,pdgc);
              if(Length>MaxPath)
                MaxPath=Length;

            }
        }


      //bottom
      LOG("GROOTGeom",pINFO) << "Box surface scanned: [BOTTOM]";
      igen=0;
      while(igen<maxPoints)
        {
          igen++;
          xyz[0] = ox-dx+2*dx*r3.Rndm();
          xyz[1] = oy-dy;
          xyz[2] = oz-dz+2*dz*r3.Rndm();

          Rgen=0;
          while(Rgen<maxRays)
            {
              Rgen++;
              direction[0]=-0.5+r3.Rndm();
              direction[1]=r3.Rndm();
              direction[2]=-0.5+r3.Rndm();

              dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

              direction[0]/=dirTot;
              direction[1]/=dirTot;
              direction[2]/=dirTot;

              Length=ComputeMaxPathLengthPDG(xyz,direction,pdgc);
              if(Length>MaxPath)
                MaxPath=Length;
            }
        }

      //left
      LOG("GROOTGeom",pINFO) << "Box surface scanned: [LEFT]";
      igen=0;
      while(igen<maxPoints)
        {
          igen++;
          xyz[0] = ox-dx;
          xyz[1] = oy-dy+2*dy*r3.Rndm();
          xyz[2] = oz-dz+2*dz*r3.Rndm();

          Rgen=0;
          while(Rgen<maxRays)
            {
              Rgen++;
              direction[0]=r3.Rndm();
              direction[1]=-0.5+r3.Rndm();
              direction[2]=-0.5+r3.Rndm();

              dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

              direction[0]/=dirTot;
              direction[1]/=dirTot;
              direction[2]/=dirTot;

              Length=ComputeMaxPathLengthPDG(xyz,direction,pdgc);
              if(Length>MaxPath)
                MaxPath=Length;
            }
        }

      //right
      LOG("GROOTGeom",pINFO) << "Box surface scanned: [RIGHT]";
      igen=0;
      while(igen<maxPoints)
        {
          igen++;
          xyz[0] = ox+dx;
          xyz[1] = oy-dy+2*dy*r3.Rndm();
          xyz[2] = oz-dz+2*dz*r3.Rndm();

          Rgen=0;
          while(Rgen<maxRays)
            {
              Rgen++;
              direction[0]=-r3.Rndm();
              direction[1]=-0.5+r3.Rndm();
              direction[2]=-0.5+r3.Rndm();

              dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

              direction[0]/=dirTot;
              direction[1]/=dirTot;
              direction[2]/=dirTot;

              Length=ComputeMaxPathLengthPDG(xyz,direction,pdgc);
              if(Length>MaxPath)
                MaxPath=Length;
            }
        }

      //back
      LOG("GROOTGeom",pINFO) << "Box surface scanned: [BACK]";
      igen=0;
      while(igen<maxPoints)
        {
          igen++;
          xyz[0] = ox-dx+2*dx*r3.Rndm();
          xyz[1] = oy-dy+2*dy*r3.Rndm();
          xyz[2] = oz-dz;

          Rgen=0;
          while(Rgen<maxRays)
            {
              Rgen++;
              direction[0]=-0.5+r3.Rndm();
              direction[1]=-0.5+r3.Rndm();
              direction[2]=r3.Rndm();

              dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

              direction[0]/=dirTot;
              direction[1]/=dirTot;
              direction[2]/=dirTot;

              Length=ComputeMaxPathLengthPDG(xyz,direction,pdgc);
              if(Length>MaxPath)
                MaxPath=Length;
            }
        }

      //front
      LOG("GROOTGeom",pINFO) << "Box surface scanned: [FRONT]";
      igen=0;
      while(igen<maxPoints)
        {
          igen++;
          xyz[0] = ox-dx+2*dx*r3.Rndm();
          xyz[1] = oy-dy+2*dy*r3.Rndm();
          xyz[2] = oz+dz;

          Rgen=0;
          while(Rgen<maxRays)
            {
              Rgen++;
              direction[0]=-0.5+r3.Rndm();
              direction[1]=-0.5+r3.Rndm();
              direction[2]=-r3.Rndm();

              dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

              direction[0]/=dirTot;
              direction[1]/=dirTot;
              direction[2]/=dirTot;

              Length=ComputeMaxPathLengthPDG(xyz,direction,pdgc);
              if(Length>MaxPath)
                MaxPath=Length;
            }
        }

      LOG("GROOTGeom", pINFO) << "Max path length found = " << MaxPath;
      fCurrMaxPathLengthList->AddPathLength(pdgc,MaxPath);
    }

  this->ScalePathLengths(*fCurrMaxPathLengthList);

  return *fCurrMaxPathLengthList;
}
//________________________________________________________________________
void ROOTGeomAnalyzer::Initialize(string filename)
{
  LOG("GROOTGeom", pINFO)
               << "Initializing with input geometry from: " << filename;

  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
  if (!is_accessible) {
     LOG("GROOTGeom", pERROR)
       << "The ROOT geometry doesn't exist! Initialization failed!";
     return;
  }
  fGeometry = TGeoManager::Import(filename.c_str());

  fCurrMaxPathLengthList = 0;
  fCurrPathLengthList    = 0;
  fCurrPDGCodeList       = 0;

  this->BuildListOfTargetNuclei();

  const PDGCodeList & pdglist = this->ListOfTargetNuclei();

  fCurrPathLengthList    = new PathLengthList(pdglist);
  fCurrMaxPathLengthList = new PathLengthList(pdglist);
  fCurrVertex            = new TVector3(0.,0.,0.);

  // some defaults:
  this -> SetScannerNPoints (200);
  this -> SetScannerNRays   (200);
  this -> SetUnits          (genie::units::meter);
}
//___________________________________________________________________________
const PDGCodeList & ROOTGeomAnalyzer::ListOfTargetNuclei(void)
{
  return *fCurrPDGCodeList;
}
//___________________________________________________________________________
const PathLengthList & ROOTGeomAnalyzer::ComputePathLengths(
                          const TLorentzVector & x, const TLorentzVector & p)
{
  // Computes the path-length within each detector material for a neutrino
  // starting from point x and travelling along the direction of p
  // reset current list of path-lengths

  LOG("GROOTGeom", pINFO) << "Computing path-lengths for the input neutrino";

  fCurrPathLengthList->SetAllToZero();

  LOG("GROOTGeom", pDEBUG)
       << "\nInput nu: 4p = " << utils::print::P4AsShortString(&p)
       << ", 4x = " << utils::print::X4AsString(&x);

  TGeoVolume *current =0;

  int FlagNotInYet(0);
  bool condition(kTRUE);

  float xx,yy,zz;
  double xyz[3];
  double step(0);
  xx=x.X();
  yy=x.Y();
  zz=x.Z();

  fGeometry->SetCurrentDirection(p.Px()/p.P(),p.Py()/p.P(),p.Pz()/p.P());

  while((!FlagNotInYet) || condition)
    {
      condition=kTRUE;

      LOG("GROOTGeom",pDEBUG)<<" x "<<xx<<" y "<<yy<<" z "<<zz<<" flag not yet in "<<FlagNotInYet;
      xyz[0]=xx;
      xyz[1]=yy;
      xyz[2]=zz;

      fGeometry->SetCurrentPoint(xyz);
      fGeometry->FindNode(xyz[0],xyz[1],xyz[2]);
      current =  fGeometry->GetCurrentVolume();
      LOG("GROOTGeom",pDEBUG)<<" current volume "<<current->GetName();
      TGeoMedium *med;
      TGeoMaterial *mat;

      if (fGeometry->IsOutside() || !current)
        {
          condition=kFALSE;

          if(FlagNotInYet)
            break;

          fGeometry->FindNextBoundary();
          step=fGeometry->GetStep();
          while(!fGeometry->IsEntering())
            {
              fGeometry->Step();
              step=fGeometry->GetStep();
            }

          xx+=step * p.Px()/p.P();
          yy+=step * p.Py()/p.P();
          zz+=step * p.Pz()/p.P();
        }

      if(condition)
        {
          if(!FlagNotInYet)
            FlagNotInYet=1;

          med = current->GetMedium();
          LOG("GROOTGeom",pDEBUG)<<" current medium "<<med->GetName();
          if (!med)
            condition=kFALSE;
        }

      if(condition)
        {
          mat = med->GetMaterial();
          LOG("GROOTGeom",pDEBUG)<<" current material "<<mat->GetName();
          LOG("GROOTGeom",pDEBUG)<<" current material A"<<mat->GetA();
          LOG("GROOTGeom",pDEBUG)<<" current material Z"<<mat->GetZ();
          LOG("GROOTGeom",pDEBUG)<<" current material is mix "<<mat->IsMixture();
          if (!mat)
            condition=kFALSE;
        }

      if(condition)
        {
          if(mat->IsMixture())
            {
              int Nelements((dynamic_cast <TGeoMixture*> (mat))->GetNelements());
              TGeoElement *ele;
              for(int i=0;i<Nelements;i++)
                {
                  LOG("GROOTGeom",pDEBUG)<<" number of elements "<<Nelements;
                  ele=(dynamic_cast <TGeoMixture*> (mat))->GetElement(i);
                  int ion_pdgc = this->GetTargetPdgCode(ele);
                  fGeometry->FindNextBoundary();
                  step=fGeometry->GetStep();
                  while(!fGeometry->IsEntering())
                    {
                      fGeometry->Step();
                      step=fGeometry->GetStep();
                    }
                  LOG("GROOTGeom",pDEBUG) <<" IsEntering   = "
                        << utils::print::BoolAsYNString(fGeometry->IsEntering());
                  LOG("GROOTGeom",pDEBUG) <<" IsOnBoundary = "
                        << utils::print::BoolAsYNString(fGeometry->IsOnBoundary());
                  LOG("GROOTGeom",pDEBUG)
                      <<" PDG-Code = " << ion_pdgc << ", Step = "<<step;

                  fCurrPathLengthList->AddPathLength(ion_pdgc,step);
                  xx+=step * p.Px()/p.P();
                  yy+=step * p.Py()/p.P();
                  zz+=step * p.Pz()/p.P();
                }
            }
          else
            {
              int ion_pdgc = this->GetTargetPdgCode(mat);
              fGeometry->FindNextBoundary();
              step=fGeometry->GetStep();
              while(!fGeometry->IsEntering())
                {
                  fGeometry->Step();
                  step=fGeometry->GetStep();
                }
              LOG("GROOTGeom",pDEBUG) <<" IsEntering   = "
                   << utils::print::BoolAsYNString(fGeometry->IsEntering());
              LOG("GROOTGeom",pDEBUG) <<" IsOnBoundary = "
                   << utils::print::BoolAsYNString(fGeometry->IsOnBoundary());
              LOG("GROOTGeom",pDEBUG)
                   <<" PDG-Code = " << ion_pdgc << ", Step = "<<step;

              fCurrPathLengthList->AddPathLength(ion_pdgc,step);
              xx+=step * p.Px()/p.P();
              yy+=step * p.Py()/p.P();
              zz+=step * p.Pz()/p.P();
            }

        }
    }

  this->ScalePathLengths(*fCurrPathLengthList);

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const TVector3 & ROOTGeomAnalyzer::GenerateVertex(
               const TLorentzVector & x, const TLorentzVector & p, int tgtpdg)
{
// Generates a random vertex, within the detector material with the input
// PDG code, for a neutrino starting from point x and travelling along the
// direction of p

  LOG("GROOTGeom", pINFO)
           << "Generating vtx in material: " << tgtpdg
                                     << " along the input neutrino direction";
  // reset current interaction vertex
  fCurrVertex->SetXYZ(0.,0.,0.);

  LOG("GROOTGeom", pDEBUG)
       << "\nInput nu: 4p = " << utils::print::P4AsShortString(&p)
       << ", 4x = " << utils::print::X4AsString(&x);

  if(!fGeometry)
    {
      LOG("GROOTGeom",pERROR) << "No ROOT geometry is loaded!";
      return *fCurrVertex;
    }

  //calculate length weighted with density
  double dist(0);

  TGeoVolume *current =0;

  int FlagNotInYet(0);
  bool condition(kTRUE);

  float xx,yy,zz;
  double xyz[3];
  double step(0);
  xx=x.X();
  yy=x.Y();
  zz=x.Z();

  fGeometry->SetCurrentDirection(p.Px()/p.P(),p.Py()/p.P(),p.Pz()/p.P());

  while((!FlagNotInYet) || condition)
    {
      condition=kTRUE;

      LOG("GROOTGeom",pDEBUG)<<" x "<<xx<<" y "<<yy<<" z "<<zz<<" flag not yet in "<<FlagNotInYet;
      xyz[0]=xx;
      xyz[1]=yy;
      xyz[2]=zz;

      fGeometry->SetCurrentPoint(xyz);
      fGeometry->FindNode(xyz[0],xyz[1],xyz[2]);
      current =  fGeometry->GetCurrentVolume();
      LOG("GROOTGeom",pDEBUG)<<" current volume "<<current->GetName();
      TGeoMedium *med;
      TGeoMaterial *mat;

      if (fGeometry->IsOutside() || !current)
        {
          condition=kFALSE;

          if(FlagNotInYet)
            break;

          fGeometry->FindNextBoundary();
          step=fGeometry->GetStep();
          while(!fGeometry->IsEntering())
            {
              fGeometry->Step();
              step=fGeometry->GetStep();
            }

          xx+=step * p.Px()/p.P();
          yy+=step * p.Py()/p.P();
          zz+=step * p.Pz()/p.P();
        }

      if(condition)
        {
          if(!FlagNotInYet)
            FlagNotInYet=1;

          med = current->GetMedium();
          LOG("GROOTGeom",pDEBUG)<<" current medium "<<med->GetName();
          if (!med)
            condition=kFALSE;
        }

      if(condition)
        {
          mat = med->GetMaterial();
          LOG("GROOTGeom",pDEBUG)<<" current material "<<mat->GetName();
          LOG("GROOTGeom",pDEBUG)<<" current material A"<<mat->GetA();
          LOG("GROOTGeom",pDEBUG)<<" current material Z"<<mat->GetZ();
          LOG("GROOTGeom",pDEBUG)<<" current material is mix "<<mat->IsMixture();
          if (!mat)
            condition=kFALSE;
        }

      if(condition)
        {
          int ion_pdgc = this->GetTargetPdgCode(mat);
          double density(mat->GetDensity());
          fGeometry->FindNextBoundary();
          step=fGeometry->GetStep();
          while(!fGeometry->IsEntering())
            {
              fGeometry->Step();
              step=fGeometry->GetStep();
            }

          if(ion_pdgc == tgtpdg)
            dist+=(step*density);

          xx+=step * p.Px()/p.P();
          yy+=step * p.Py()/p.P();
          zz+=step * p.Pz()/p.P();

        }
    }

  if(dist==0)
    {
      LOG("GROOTGeom",pERROR)
           <<"No material selected along this direction from set point!!! ";
      return *fCurrVertex;
    }

  LOG("GROOTGeom",pDEBUG)<<" dist times density "<<dist;
  //generate random number between 0 and dist
  RandomGen* rand=RandomGen::Instance();
  TRandom & r3=rand->Random3();
  double distVertex(r3.Rndm()*dist);
  LOG("GROOTGeom",pDEBUG)<<" Random distance in selected material "<<distVertex;

  //-- generate the vertex

  xx=x.X();
  yy=x.Y();
  zz=x.Z();
  double StepIncrease(0.001);
  double distToVtx(0);
  FlagNotInYet=0;
  condition=kTRUE;

  while(((!FlagNotInYet) || condition) && distToVtx<distVertex)
    {
      condition=kTRUE;

      LOG("GROOTGeom",pDEBUG)<<" x "<<xx<<" y "<<yy<<" z "<<zz<<" flag not yet in "<<FlagNotInYet;
      xx+=StepIncrease * p.Px()/p.P();
      yy+=StepIncrease * p.Py()/p.P();
      zz+=StepIncrease * p.Pz()/p.P();
      xyz[0]=xx;
      xyz[1]=yy;
      xyz[2]=zz;

      fGeometry->SetCurrentPoint(xyz);
      fGeometry->FindNode(xyz[0],xyz[1],xyz[2]);
      current =  fGeometry->GetCurrentVolume();
      LOG("GROOTGeom",pDEBUG)<<" current volume "<<current->GetName();
      TGeoMedium *med;
      TGeoMaterial *mat;

      if (fGeometry->IsOutside() || !current)
        {
          condition=kFALSE;

          if(FlagNotInYet)
            break;
        }

      if(condition)
        {
          if(!FlagNotInYet)
            FlagNotInYet=1;

          med = current->GetMedium();
          LOG("GROOTGeom",pDEBUG)<<" current medium "<<med->GetName();
          if (!med)
            condition=kFALSE;
        }

      if(condition)
        {
          mat = med->GetMaterial();
          LOG("GROOTGeom",pDEBUG)<<" current material "<<mat->GetName();
          LOG("GROOTGeom",pDEBUG)<<" current material A"<<mat->GetA();
          LOG("GROOTGeom",pDEBUG)<<" current material Z"<<mat->GetZ();
          LOG("GROOTGeom",pDEBUG)<<" current material is mix "<<mat->IsMixture();
          if (!mat)
            condition=kFALSE;
        }

      if(condition)
        {
          int ion_pdgc = this->GetTargetPdgCode(mat);
          double density(mat->GetDensity());
          if(ion_pdgc == tgtpdg)
            distToVtx+=(StepIncrease*density);
        }
    }

  xx-=StepIncrease * p.Px()/p.P();
  yy-=StepIncrease * p.Py()/p.P();
  zz-=StepIncrease * p.Pz()/p.P();

  fCurrVertex->SetXYZ(xx,yy,zz);

  LOG("GROOTGeom",pDEBUG)<<" Vtx : x "<<xx<<" y "<<yy<<" z "<<zz;
  return *fCurrVertex;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::BuildListOfTargetNuclei(void)
{
  fCurrPDGCodeList = new PDGCodeList;

  if(!fGeometry)
    {
      LOG("GROOTGeom", pERROR) << "No ROOT geometry is loaded!";
      return;
    }

  TObjArray *LV = new TObjArray();

  LV=fGeometry->GetListOfVolumes();
  int numVol;

  numVol=(LV->GetEntries());
  TGeoVolume *TV =0;

  for(Int_t i=0;i<numVol;i++)
    {
      TV= dynamic_cast <TGeoVolume *>(LV->At(i));
      TGeoMedium *med;
      TGeoMaterial *mat;
      med = TV->GetMedium();
      mat = med->GetMaterial();
      if(mat->IsMixture())
        {
          int Nelements((dynamic_cast <TGeoMixture*> (mat))->GetNelements());
          TGeoElement *ele;
          for(int i=0;i<Nelements;i++)
            {
              ele=(dynamic_cast <TGeoMixture*> (mat))->GetElement(i);
              int ion_pdgc = this->GetTargetPdgCode(ele);
              fCurrPDGCodeList->push_back(ion_pdgc);
            }
        }
      else
        {
          int ion_pdgc = this->GetTargetPdgCode(mat);
          fCurrPDGCodeList->push_back(ion_pdgc);
        }
    }
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::ComputeMaxPathLengthPDG(double* XYZ,double* direction,int pdgc)
{

  TGeoVolume *current =0;
  int counterloop(0);
  double Length(0);

  int FlagNotInYet(0);
  bool condition(kTRUE);

  float xx,yy,zz;
  double xyz[3];
  double step(0);
  xx=XYZ[0];
  yy=XYZ[1];
  zz=XYZ[2];

  fGeometry->SetCurrentDirection(direction);

  while(((!FlagNotInYet) || condition) && counterloop <100)
    {
      counterloop++;
      condition=kTRUE;

      xyz[0]=xx;
      xyz[1]=yy;
      xyz[2]=zz;

      fGeometry->SetCurrentPoint(xyz);
      fGeometry->FindNode(xyz[0],xyz[1],xyz[2]);
      current =  fGeometry->GetCurrentVolume();
      TGeoMedium *med;
      TGeoMaterial *mat;

      if (fGeometry->IsOutside() || !current)
        {
          condition=kFALSE;

          if(FlagNotInYet)
            break;

          fGeometry->FindNextBoundary();
          step=fGeometry->GetStep();
          while(!fGeometry->IsEntering())
            {
              fGeometry->Step();
              step=fGeometry->GetStep();
            }

          xx+=step * direction[0];
          yy+=step * direction[1];
          zz+=step * direction[2];
        }

      if(condition)
        {
          if(!FlagNotInYet)
            FlagNotInYet=1;

          med = current->GetMedium();
          if (!med)
            condition=kFALSE;
        }

      if(condition)
        {
          mat = med->GetMaterial();
          if (!mat)
            condition=kFALSE;
        }

      if(condition)
        {
          int ion_pdgc = this->GetTargetPdgCode(mat);
          fGeometry->FindNextBoundary();
          step=fGeometry->GetStep();
          while(!fGeometry->IsEntering())
            {
              fGeometry->Step();
              step=fGeometry->GetStep();
            }

          if(ion_pdgc == pdgc)
            Length+=step;

          xx+=step * direction[0];
          yy+=step * direction[1];
          zz+=step * direction[2];
        }
    }
  return Length;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::ScalePathLengths(PathLengthList & pl)
{
// convert path lenghts to default GENIE length scale
//
  LOG("GROOTGeom", pDEBUG) 
              << "Scaling path-lengths -> meters (scale = " << fScale << ")";

  PathLengthList::iterator pliter;
  for(pliter = pl.begin(); pliter != pl.end(); ++pliter) 
  {
    int pdgc = pliter->first;
    pl.ScalePathLength(pdgc,fScale);
  }
}
//___________________________________________________________________________
int ROOTGeomAnalyzer::GetTargetPdgCode(const TGeoMaterial * const m) const
{
  int A = int(m->GetA());
  int Z = int(m->GetZ());

  int pdgc = pdg::IonPdgCode(A,Z);

  return pdgc;
}
//___________________________________________________________________________
int ROOTGeomAnalyzer::GetTargetPdgCode(const TGeoElement * const e) const
{
  int A = int(e->A());
  int Z = int(e->Z());

  int pdgc = pdg::IonPdgCode(A,Z);

  return pdgc;
}
//___________________________________________________________________________

