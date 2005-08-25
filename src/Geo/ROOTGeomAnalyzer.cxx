//____________________________________________________________________________
/*!

\class   genie::ROOTGeomAnalyzer

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

#include "Geo/ROOTGeomAnalyzer.h"
#include "Geo/PathLengthList.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

#include <TPolyMarker3D.h>
#include <TGeoBBox.h>
#include <TRandom3.h>

using namespace genie;

//___________________________________________________________________________
ROOTGeomAnalyzer::ROOTGeomAnalyzer(string filename) :
GeomAnalyzerI()
{
  this->Initialize(filename);
}
//___________________________________________________________________________
ROOTGeomAnalyzer::~ROOTGeomAnalyzer()
{
  if( fCurrPathLengthList ) delete fCurrPathLengthList;
  if( fCurrPDGCodeList    ) delete fCurrPDGCodeList;
  if( fGeometry) delete fGeometry;

}
//___________________________________________________________________________
double ROOTGeomAnalyzer::SetVtxMaterial(int pdgc)
{

  if(!fGeometry)
    {
      std::cout<<" ERROR!!! Load geometry before setting the material!!! "<<std::endl;
      fMaterial=-1;
      return 0;
    }

  fMaterial=pdgc;

  if(!fCurrPDGCodeList->CheckMaterial(fMaterial))
    {
      std::cout<<" ERROR!!! The selected material does not exist!!! "<<std::endl;
      fMaterial=-1;
      return 0;
    }
  else
    {
      std::cout<<" Material selected : "<<fMaterial<<std::endl;
    }

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
      std::cout<<" ERROR!!! The World Volume does not exist in your geometry!!! "<<std::endl;
      fMaterial=-1;
      return 0;
    }

  //generate 200 random points on each surface, use 200 rays to calculate maximum path for selected material
   
  TS=TVWorld->GetShape();
  TGeoBBox *box=(TGeoBBox *)TS;

  double dx = box->GetDX();
  double dy = box->GetDY();
  double dz = box->GetDZ();
  double ox = (box->GetOrigin())[0];
  double oy = (box->GetOrigin())[1];
  double oz = (box->GetOrigin())[2];

  std::cout<<" max dimensions : x = "<<dx<<" ; y = "<<dy<<" ; z = "<<dz<<std::endl;
  std::cout<<" origin : x = "<<ox<<" ; y = "<<oy<<" ; z = "<<oz<<std::endl;

  gRandom = new TRandom3();
  int igen(0);
  int Rgen(0);
  int maxPoints(200);
  int maxRays(200);
  double xyz[3];
  double direction[3];
  double MaxPath(0);
  double Length(0);
  double dirTot(0);

  //top
  igen=0;
  while(igen<maxPoints)
    {
      igen++;
      xyz[0] = ox-dx+2*dx*gRandom->Rndm();
      xyz[1] = oy+dy;
      xyz[2] = oz-dz+2*dz*gRandom->Rndm();
      
      Rgen=0;
      while(Rgen<maxRays)
	{
	  Rgen++;
	  direction[0]=-0.5+gRandom->Rndm();
	  direction[1]=-gRandom->Rndm();
	  direction[2]=-0.5+gRandom->Rndm();
	  
	  dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
	  
	  direction[0]/=dirTot;
	  direction[1]/=dirTot;
	  direction[2]/=dirTot;
	   
	  Length=ComputeMaxPathLength(xyz,direction,pdgc); 
	  if(Length>MaxPath)
	    MaxPath=Length;
	  
	}
    }
  
 
  //bottom
  igen=0;
  while(igen<maxPoints)
    { 
      igen++;
      xyz[0] = ox-dx+2*dx*gRandom->Rndm();
      xyz[1] = oy-dy;
      xyz[2] = oz-dz+2*dz*gRandom->Rndm();

      Rgen=0;
      while(Rgen<maxRays)
	{ 
	  Rgen++;
	  direction[0]=-0.5+gRandom->Rndm();
	  direction[1]=gRandom->Rndm();
	  direction[2]=-0.5+gRandom->Rndm();
	  
	  dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
	  
	  direction[0]/=dirTot;
	  direction[1]/=dirTot;
	  direction[2]/=dirTot;
	  
	  Length=ComputeMaxPathLength(xyz,direction,pdgc); 
	  if(Length>MaxPath)
	    MaxPath=Length;
	}
    }

  //left  
  igen=0;
  while(igen<maxPoints)
    {
      igen++;
      xyz[0] = ox-dx;
      xyz[1] = oy-dy+2*dy*gRandom->Rndm();
      xyz[2] = oz-dz+2*dz*gRandom->Rndm(); 
      
      Rgen=0;
      while(Rgen<maxRays)
	{
	  Rgen++;
	  direction[0]=gRandom->Rndm();
	  direction[1]=-0.5+gRandom->Rndm();
	  direction[2]=-0.5+gRandom->Rndm();
	  
	  dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
	  
	  direction[0]/=dirTot;
	  direction[1]/=dirTot;
	  direction[2]/=dirTot;
	  
	  Length=ComputeMaxPathLength(xyz,direction,pdgc); 
	  if(Length>MaxPath)
	    MaxPath=Length;
	}
    }
  
  //right  
  igen=0;
  while(igen<maxPoints)
    {
      igen++;
      xyz[0] = ox+dx;
      xyz[1] = oy-dy+2*dy*gRandom->Rndm();
      xyz[2] = oz-dz+2*dz*gRandom->Rndm();
      
      Rgen=0;
      while(Rgen<maxRays)
	{
	  Rgen++;
	  direction[0]=-gRandom->Rndm();
	  direction[1]=-0.5+gRandom->Rndm();
	  direction[2]=-0.5+gRandom->Rndm();
	  
	  dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
	  
	  direction[0]/=dirTot;
	  direction[1]/=dirTot;
	  direction[2]/=dirTot;
	  
	  Length=ComputeMaxPathLength(xyz,direction,pdgc); 
	  if(Length>MaxPath)
	    MaxPath=Length;
	}
    }
  
  //back
  igen=0;
  while(igen<maxPoints)
    {
      igen++;
      xyz[0] = ox-dx+2*dx*gRandom->Rndm();
      xyz[1] = oy-dy+2*dy*gRandom->Rndm();
      xyz[2] = oz-dz; 

      Rgen=0;
      while(Rgen<maxRays)
	{
	  Rgen++;
	  direction[0]=-0.5+gRandom->Rndm();
	  direction[1]=-0.5+gRandom->Rndm();
	  direction[2]=gRandom->Rndm();
	  
	  dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
	  
	  direction[0]/=dirTot;
	  direction[1]/=dirTot;
	  direction[2]/=dirTot;
	  
	  Length=ComputeMaxPathLength(xyz,direction,pdgc); 
	  if(Length>MaxPath)
	    MaxPath=Length;
	}
    }  

  //front
  igen=0;
  while(igen<maxPoints)
    {
      igen++;
      xyz[0] = ox-dx+2*dx*gRandom->Rndm();
      xyz[1] = oy-dy+2*dy*gRandom->Rndm();
      xyz[2] = oz+dz; 

      Rgen=0;
      while(Rgen<maxRays)
	{
	  Rgen++;
	  direction[0]=-0.5+gRandom->Rndm();
	  direction[1]=-0.5+gRandom->Rndm();
	  direction[2]=-gRandom->Rndm();
	  
	  dirTot=sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
	  
	  direction[0]/=dirTot;
	  direction[1]/=dirTot;
	  direction[2]/=dirTot;
	  
	  Length=ComputeMaxPathLength(xyz,direction,pdgc); 
	  if(Length>MaxPath)
	    MaxPath=Length;
	}
    }
 

  return MaxPath;

}
//________________________________________________________________________
void ROOTGeomAnalyzer::Initialize(string filename)
{
  fCurrPathLengthList = 0;
  fCurrPDGCodeList    = 0;

  fGeometry=TGeoManager::Import(filename.c_str());

  this->BuildListOfTargetNuclei();

  const PDGCodeList & pdglist = this->ListOfTargetNuclei();

  fCurrPathLengthList = new PathLengthList(pdglist);
  fCurrVertex = new TVector3(0.,0.,0.);
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
  
  fCurrPathLengthList->SetAllToZero();
  
  LOG("GROOTGeom", pINFO)
       << "\nComputing path-lengths for neutrino: "
       << "\n  with 4-momentum : " << print_utils::P4AsString(&p)
       << "\n  starting from   : " << print_utils::X4AsString(&x);

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

      std::cout<<" x "<<xx<<" y "<<yy<<" z "<<zz<<" flag not yet in "<<FlagNotInYet<<std::endl;
      xyz[0]=xx;
      xyz[1]=yy;
      xyz[2]=zz;
      
      fGeometry->SetCurrentPoint(xyz);
      fGeometry->FindNode(xyz[0],xyz[1],xyz[2]);
      current =  fGeometry->GetCurrentVolume();
      std::cout<<" current volume "<<current->GetName()<<std::endl;
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
          std::cout<<" current medium "<<med->GetName()<<std::endl;
          if (!med)
            condition=kFALSE;
        }

      if(condition)
        {
          mat = med->GetMaterial();
          std::cout<<" current material "<<mat->GetName()<<std::endl;
	  std::cout<<" current material A"<<mat->GetA()<<std::endl;
	  std::cout<<" current material Z"<<mat->GetZ()<<std::endl;
	  std::cout<<" current material is mix "<<mat->IsMixture()<<std::endl;
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
		  std::cout<<" number of elements "<<Nelements<<std::endl;
		  ele=(dynamic_cast <TGeoMixture*> (mat))->GetElement(i);
		  std::cout<<" test "<<ele->A()<<std::endl;
		  int A=(int)(ele->A());
		  int Z(ele->Z());
		  int ion_pdgc = pdg::IonPdgCode(A,Z);
		  fGeometry->FindNextBoundary();
		  step=fGeometry->GetStep();
		  while(!fGeometry->IsEntering())
		    {
		      fGeometry->Step(); 
		      step=fGeometry->GetStep();
		    }
		  std::cout<<" isentering? "<< fGeometry->IsEntering()<<std::endl;
		  std::cout<<" isonboundary? "<< fGeometry->IsOnBoundary()<<std::endl;
		  std::cout<<" A "<<A<<" Z "<<Z<<" code "<<ion_pdgc<<" step "<<step<<std::endl;
		  
		  fCurrPathLengthList->AddPathLength(ion_pdgc,step);
		  xx+=step * p.Px()/p.P();
		  yy+=step * p.Py()/p.P();
		  zz+=step * p.Pz()/p.P();
		}
	    }
	  else
	    {
	      int A=(int)(mat->GetA());
	      int Z(mat->GetZ());
	      int ion_pdgc = pdg::IonPdgCode(A,Z);
	      fGeometry->FindNextBoundary();
	      step=fGeometry->GetStep(); 
	      while(!fGeometry->IsEntering())
		{
		  fGeometry->Step();
		  step=fGeometry->GetStep();
		}
	      std::cout<<" isentering? "<< fGeometry->IsEntering()<<std::endl;
	      std::cout<<" isonboundary? "<< fGeometry->IsOnBoundary()<<std::endl;
	      std::cout<<" A "<<A<<" Z "<<Z<<" code "<<ion_pdgc<<" step "<<step<<std::endl;
	      
	      fCurrPathLengthList->AddPathLength(ion_pdgc,step);
	      xx+=step * p.Px()/p.P();
	      yy+=step * p.Py()/p.P();
	      zz+=step * p.Pz()/p.P();
	    }
	  
	}
    }
  return *fCurrPathLengthList;
}

//___________________________________________________________________________
const TVector3 & ROOTGeomAnalyzer::GenerateVertex(
              const TLorentzVector & x, const TLorentzVector & p, int tgtpdg)
{
// Generates a random vertex, within the detector material with the input
// PDG code, for a neutrino starting from point x and travelling along the
// direction of p

  // reset current interaction vertex
  fCurrVertex->SetXYZ(0.,0.,0.);

  LOG("GROOTGeom", pINFO)
       << "\nGenerating an vertex at material with PDG code = " << tgtpdg
       << "\nfor a neutrino: "
       << "\n  with 4-momentum : " << print_utils::P4AsString(&p)
       << "\n  starting from   : " << print_utils::X4AsString(&x);

  //-- generate the vertex
  //
  // ... ... ...

  return *fCurrVertex;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::BuildListOfTargetNuclei(void)
{
  fCurrPDGCodeList = new PDGCodeList;

  if(!fGeometry)
    {
      std::cout<<" ERROR!!! Load geometry Fisrt!!! "<<std::endl;
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
	      int A=(int)(ele->A());
	      int Z(ele->Z());
	      int ion_pdgc = pdg::IonPdgCode(A,Z);
	      fCurrPDGCodeList->push_back(ion_pdgc);
	      
	    }
	}
      else
	{
	  int A=(int)(mat->GetA());
	  int Z(mat->GetZ());
	  int ion_pdgc = pdg::IonPdgCode(A,Z);
	  fCurrPDGCodeList->push_back(ion_pdgc);
	}
    }
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::ComputeMaxPathLength(double* XYZ,double* direction,int pdgc)
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

      //std::cout<<" x "<<xx<<" y "<<yy<<" z "<<zz<<" flag not yet in "<<FlagNotInYet<<std::endl;
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
          //std::cout<<" current medium "<<med->GetName()<<std::endl;
          if (!med)
            condition=kFALSE;
        }

      if(condition)
        {
          mat = med->GetMaterial();
          //std::cout<<" current material "<<mat->GetName()<<std::endl;
	  //std::cout<<" current material A"<<mat->GetA()<<std::endl;
	  //std::cout<<" current material Z"<<mat->GetZ()<<std::endl;
	  //std::cout<<" current material is mix "<<mat->IsMixture()<<std::endl;
          if (!mat)
            condition=kFALSE;
        }

      if(condition)
        {
	  int A=(int)(mat->GetA());
	  int Z(mat->GetZ());
	  int ion_pdgc = pdg::IonPdgCode(A,Z);
	  fGeometry->FindNextBoundary();
	  step=fGeometry->GetStep();  
	  while(!fGeometry->IsEntering())
	    {
	      fGeometry->Step();
	      step=fGeometry->GetStep();
	    }
	  //std::cout<<" isentering? "<< fGeometry->IsEntering()<<std::endl;
	  //std::cout<<" isonboundary? "<< fGeometry->IsOnBoundary()<<std::endl;
	  //std::cout<<" A "<<A<<" Z "<<Z<<" code "<<ion_pdgc<<" step "<<step<<std::endl;
	  
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
void ROOTGeomAnalyzer::test(void)
{
  TGeoManager *TGM = new TGeoManager("TGM","test");
  //TGM->Import("/eth/store_6/home_6/amerega/TEST/vgm.2.03/examples/E01/N03/bin/Linux-g++/Geometry.root");
  TGM->Import("$GENIE/src/test/TestGeometry.root");

  TObjArray *LV = new TObjArray();

  //TObjArray *LV=TGM->GetListOfVolumes();
  LV=gGeoManager->GetListOfVolumes();

  int numVol;
  numVol=(LV->GetEntries());

  std::cout<<numVol<<" volumes found  "<<std::endl;

  TGeoVolume *TV = new TGeoVolume();
  TGeoVolume *TVWorld = new TGeoVolume();

  TGeoShape *TS;

  char *name;
  char *str;
  str="World";

  Double_t  point [3];
  point[0]=0;
  point[1]=0;
  point[2]=0;

  Double_t  dir [3];
  dir[0]=1;
  dir[1]=0;
  dir[2]=0;

  for(Int_t i=0;i<numVol;i++)
    {
      TV= dynamic_cast <TGeoVolume *> (LV->At(i));
      std::cout<<i<<"  "<<TV->GetName()<<" made of "<<TV->GetMaterial()->GetName()<<std::endl;
      name=const_cast<char*>(TV->GetName());
      TS=TV->GetShape();
      TS->ComputeBBox();
      std::cout<<i<< " test "<<TV->GetNdaughters()<<" contains "<<TV->Contains(point)<<" distance out "<<TS->DistFromOutside(point,dir)<<std::endl;
      if(!strcmp(str,name))
        TVWorld=TV;
      std::cout<<TVWorld->GetName()<<" FOUND "<<std::endl;
    }

  //help Andrei

  TGM->SetVisOption(0);
  TVWorld->Draw();

  TPolyMarker3D *marker = new TPolyMarker3D();
  marker->SetMarkerColor(kRed);
  marker->SetMarkerStyle(8);
  marker->SetMarkerSize(0.5);

  const TGeoShape *TS1=TVWorld->GetShape();
  TGeoBBox *box=(TGeoBBox *)TS1;

  Double_t dx = box->GetDX();
  Double_t dy = box->GetDY();
  Double_t dz = box->GetDZ();
  Double_t ox = (box->GetOrigin())[0];
  Double_t oy = (box->GetOrigin())[1];
  Double_t oz = (box->GetOrigin())[2];

  std::cout<<" max dimensions : x = "<<dx<<" ; y = "<<dy<<" ; z = "<<dz<<std::endl;
  std::cout<<" origin : x = "<<ox<<" ; y = "<<oy<<" ; z = "<<oz<<std::endl;

  gRandom = new TRandom3();
  Int_t igen(0);
  Double_t xyz[3];

  Int_t found(0);
  //while(igen<1)
 char *nameMat;
 char *strmat;

  while(found==0 && igen<100)
    {
      // xyz[0] = ox-dx+2*dx*gRandom->Rndm();
      //xyz[1] = oy-dy+2*dy*gRandom->Rndm();
      //xyz[2] = oz-dz+2*dz*gRandom->Rndm();

      xyz[0] = 0.5;
      xyz[1] = 0;
      xyz[2] = 0;

      std::cout<<" random generated point: x = "<<xyz[0]<<" ; y = "<<xyz[1]<<" ; z = "<<xyz[2]<<std::endl;

      //TGM->SetCurrentPoint(xyz);
      gGeoManager->SetCurrentPoint(xyz);
      igen++;
      //TGM->FindNode(xyz[0],xyz[1],xyz[2]);
      gGeoManager->FindNode(xyz[0],xyz[1],xyz[2]);

      Bool_t condition;
      condition=kTRUE;

      TGeoVolume *current = gGeoManager->GetCurrentVolume();
      std::cout<<" current volume "<<current->GetName()<<std::endl;
      if (TGM->IsOutside() || !current)
        condition=kFALSE;

      TGeoMedium *med;
      TGeoMaterial *mat;

      if(condition)
        {
          med = current->GetMedium();
          std::cout<<" current medium "<<med->GetName()<<std::endl;
          if (!med)
            condition=kFALSE;
        }
      if(condition)
        {
          mat = med->GetMaterial();
          std::cout<<" current material "<<mat->GetName()<<std::endl;
          if (!mat)
            condition=kFALSE;
        }

      if(condition)
        {

          //strmat="Lead";
          //strmat="liquidArgon";
          strmat="Galactic";
          nameMat=const_cast<char*>(mat->GetName());
          if(!strcmp(strmat,nameMat))
            found=1;
        }

      marker->SetNextPoint(xyz[0],xyz[1],xyz[2]);

    }

  if(found==1)
    std::cout<<" found point : x = "<<xyz[0]<<" ; y = "<<xyz[1]<<" ; z = "<<xyz[2]<<" ; in material : "<<strmat<<std::endl;
  else
    std::cout<<" point not found!!!!"<<std::endl;
  marker->Draw("same");
}
//___________________________________________________________________________
