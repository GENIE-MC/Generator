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
ROOTGeomAnalyzer::ROOTGeomAnalyzer(void) :
GeomAnalyzerI()
{
  this->Initialize();
}
//___________________________________________________________________________
ROOTGeomAnalyzer::~ROOTGeomAnalyzer()
{
  if( fCurrPathLengthList ) delete fCurrPathLengthList;
  if( fCurrPDGCodeList    ) delete fCurrPDGCodeList;
  if( fGeometry) delete fGeometry;

}
//___________________________________________________________________________
void ROOTGeomAnalyzer::Load(char* filename)
{
  fGeometry=TGeoManager::Import(filename);
}
//___________________________________________________________________________
int ROOTGeomAnalyzer::SetVtxMaterial(char* material)
{
  if(!fGeometry)
    {
      std::cout<<" WARNING!!! Load geometry before setting the material!!! "<<std::endl;
      return 0;
    }

  fMaterial=material;

  TObjArray *LV = new TObjArray();

  LV=fGeometry->GetListOfVolumes();
  int numVol;

  numVol=(LV->GetEntries());
  TGeoVolume *TV = new TGeoVolume();
  int MaterialFound(0);
  char* MaterialName;

  for(int i=0;i<numVol;i++)
    {
      TV= dynamic_cast <TGeoVolume *>(LV->At(i));

      // std::cout<<i<<"  "<<TV->GetName()<<" made of "<<TV->GetMaterial()->GetName()<<std::endl;

      MaterialName=const_cast<char*>(TV->GetMaterial()->GetName());

      if(!strcmp(material,MaterialName))
        {
          MaterialFound=1;
          //std::cout<<material<<" FOUND "<<std::endl;
        }
    }

  if(!MaterialFound)
    {
      std::cout<<" WARNING!!! No volume made of the selected material "<<std::endl;
      delete LV;
      delete TV;
      return 0;
    }

  delete LV;
  delete TV;
  return 1;
}
//________________________________________________________________________
TGeoVolume* ROOTGeomAnalyzer::GetWorldVolume(void)
{
  if(!fGeometry)
    {
      std::cout<<" WARNING!!! Load geometry before looking for the Wrld!!! "<<std::endl;
      return 0;
    }
  TObjArray *LV = new TObjArray();

  LV=fGeometry->GetListOfVolumes();
  int numVol;

  numVol=(LV->GetEntries());
  TGeoVolume *TV = new TGeoVolume();
  TGeoVolume *TVWorld = new TGeoVolume();

  char *name;
  char *str;
  str="World";

  for(Int_t i=0;i<numVol;i++)
    {
      TV= dynamic_cast <TGeoVolume *>(LV->At(i));
      name=const_cast<char*>(TV->GetName());
      if(!strcmp(str,name))
        TVWorld=TV;
    }

  return TVWorld;
}
//________________________________________________________________________
void ROOTGeomAnalyzer::Initialize(void)
{
  fCurrPathLengthList = 0;
  fCurrPDGCodeList    = 0;

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
  float step(0.1);

  xx=x.X();
  yy=x.Y();
  zz=x.Z();

  while((!FlagNotInYet) || condition)
    {
      xx+=step * p.Px()/p.P();
      yy+=step * p.Py()/p.P();
      zz+=step * p.Pz()/p.P();

      std::cout<<" x "<<xx<<" y "<<yy<<" z "<<zz<<std::endl;
      xyz[0]=xx;
      xyz[1]=yy;
      xyz[2]=zz;
      fGeometry->SetCurrentPoint(xyz);
      fGeometry->FindNode(xyz[0],xyz[1],xyz[2]);
      //fGeometry->FindNode();
      //current =  dynamic_cast <TGeoVolume *>(gGeoManager->GetCurrentVolume());
      
      //gGeoManager->SetCurrentPoint(xyz);
      //gGeoManager->FindNode(xyz[0],xyz[1],xyz[2]);
      //current =  gGeoManager->GetCurrentVolume();
      current =  fGeometry->GetCurrentVolume();
      std::cout<<" current volume "<<current->GetName()<<std::endl;
      TGeoMedium *med;
      TGeoMaterial *mat;

      if (fGeometry->IsOutside() || !current)
        condition=kFALSE;

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
		  ele=mat->GetElement(i);
		  int A(ele->A());
		  int Z(ele->Z());
		  char strA[3];
		  char strZ[3];
		  char strMat[11];
		  if(A>99)
		    sprintf( strA, "%d", A );
		  else if(A>9)
		    sprintf( strA, "0%d", A );
		  else
		    sprintf( strA, "00%d", A );
		  
		  if(Z>99)
		    sprintf( strZ, "%d", Z );
		  else if(Z>9)
		    sprintf( strZ, "0%d", Z );
		  else
		    sprintf( strZ, "00%d", Z );

		  sprintf( strMat, "1%s%s0000", strA,strZ);
		  std::cout<<" a "<<strA<<" z "<<strZ<<" mat "<<strMat<<std::endl;
		  
		  fCurrPathLengthList->AddPathLength(atoi(strMat),step);
		
		}
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

  //-- fill in list using the input TGeoVolume
  //
  //... ...
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::test(void)
{
  TGeoManager *TGM = new TGeoManager("TGM","test");
  //TGM->Import("/eth/store_6/home_6/amerega/TEST/vgm.2.03/examples/E01/N03/bin/Linux-g++/Geometry.root");
  TGM->Import("$GENIE/src/test/TestGeometry.root");

  TObjArray *LV = new TObjArray();
  TObjArray *LN = new TObjArray();

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
