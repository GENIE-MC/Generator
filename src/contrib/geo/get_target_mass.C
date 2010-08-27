//****************************************************************************
//*
//* get_target_mass()
//*
//* Reads-in a ROOT geometry and prints-out its mass composition.
//*
//* Input arguments:
//*    - geometry_file
//*    - topvolname    (default: "")
//*    - length_unit   (default: units::cm)
//*    - density_unit  (default: units::g_cm3)
//*    - checkoverlaps (default: kFALSE)
//*
//*
//* T. Ferber (University Hamburg)
//* Luruper Chaussee 149
//* 20259 Hamburg
//* torben.ferber@desy.de
//*
//****************************************************************************

#include <iostream>
#include <string>
#include <iostream>
#include <iomanip>

#include "../../Conventions/Units.h"

//#define _debug_

using namespace std;
using namespace genie;

Int_t set_top_volume(string name);
void  get_mass      (double length_unit, double density_unit);

//____________________________________________________________________________
void get_target_mass ( 
  string geometry_file, string topvolname = "", 
  double length_unit=units::cm, double density_unit=units::g_cm3, 
  Bool_t checkoverlaps = kFALSE)
{
   // import geometry
   TGeoManager *tgeo = new TGeoManager("TGeo","TGeo");
   tgeo->Import( geometry_file.c_str() );

   // set top volume 
   if( ! set_top_volume(topvolname) ) return;

   // draw
   gGeoManager->GetTopVolume()->Draw();

   // method fails if there are volume overlaps, check it if the user wants to... (takes some time)
   if(checkoverlaps) {
      Double_t overlap_acc = 0.05;
      gGeoManager->CheckOverlaps( overlap_acc );
   }

   get_mass(length_unit, density_unit);
}
//____________________________________________________________________________
Int_t set_top_volume(string topvolname)
{	
   // no user input, set to overal top volume 
   if( topvolname=="" ) {
      topvolname = gGeoManager->GetTopVolume()->GetName();
   }
   TGeoVolume * topvol = gGeoManager->GetVolume( topvolname.c_str() );
   if (!topvol) {
     cout << "top volume does not exist" << endl;
     return 0;
   }
   gGeoManager->SetTopVolume(topvol);
   return 1;
}
//____________________________________________________________________________
void get_mass(Double_t length_unit, Double_t density_unit)
{
   // calc unit conversion factors
   Double_t density_unit_to_SI = density_unit / units::kg_m3;
   Double_t length_unit_to_SI  = length_unit  / units::m;
   Double_t volume_unit_to_SI  = TMath::Power(length_unit_to_SI, 3.);
#ifdef _debug_
   cout << "Input density unit --> kg/m^3 : x" << density_unit_to_SI << endl;
   cout << "Input length  unit --> m      : x" << length_unit_to_SI  << endl;
#endif

   // get materials in geometry
   TList *matlist = gGeoManager->GetListOfMaterials();
   if (!matlist ) { 
     cout << "Null list of materials!" << endl; 
     return; 
   } else {
#ifdef _debug_
     matlist->Print();
#endif
   }
   int max_idx = 0; // number of mixtures in geometry
   Int_t nmat = matlist->GetEntries();
   for( Int_t imat = 0; imat < nmat; imat++ )
   {
      Int_t idx = gGeoManager->GetMaterial(imat)->GetIndex();
      max_idx = TMath::Max(max_idx, idx);
   }

#ifdef _debug_
   cout << "max_idx = " << max_idx << endl;
   cout << "nmat    = " << nmat    << endl;
#endif

   TGeoVolume * topvol = gGeoManager->GetTopVolume(); //get top volume
   if (!topvol) {
     cout << "volume does not exist" << endl;
     return;
   }
   TGeoIterator NodeIter(topvol);
   TGeoNode *node=NodeIter();
   NodeIter.SetType(0); // include  all daughters

   Double_t * volume = new Double_t[max_idx+1];
   Double_t * mass   = new Double_t[max_idx+1];

   volume[ topvol->GetMaterial()->GetIndex() ] = topvol->Capacity(); //iterator does not include topvolume

   Int_t first = 1;
   while ( node=NodeIter() ) 
   {
      if( first )
      {
         NodeIter.Reset();
         first = 0;
      }	

      Int_t momidx = node->GetMotherVolume()->GetMaterial()->GetIndex() ;
      Int_t idx    = node->GetVolume()      ->GetMaterial()->GetIndex() ;

      Double_t node_vol = node->GetVolume()->Capacity() * volume_unit_to_SI;

      volume[ momidx ] -= node_vol; //substract subvolume from mother
      volume[ idx    ] += node_vol;
   }

   //
   // print out volume/mass for each `material'
   //

   for( Int_t i=0; i<gGeoManager->GetListOfMaterials()->GetEntries(); i++ )
   {
      Int_t    idx     = gGeoManager->GetMaterial(i)->GetIndex();
      Double_t density = gGeoManager->GetMaterial(i)->GetDensity() * density_unit_to_SI;

      mass[idx] = density * volume[idx];

      if( volume[idx] > 1e-20) {
        cout << setw(6) << i << "  " 
             << setw(15) << gGeoManager->GetMaterial(i)->GetName() << ": "
             << setprecision(6) 
             << setw(20) << volume[idx] << " m^3 " 
             << setw(20) << mass[idx] << " kg" 
             <<  endl;
      }
   }

   delete [] volume;
   delete [] mass; 

   //
   // print out mass contribution for each nuclear target
   //

   // ...
   // ...
   // ...
   // ...

}
//____________________________________________________________________________

