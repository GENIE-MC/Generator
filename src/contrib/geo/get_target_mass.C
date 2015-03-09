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
//* @ Aug 30, 2010 - TF
//* added list of isotope masses, added index check, added arrary initialization, corrected iterator initialization
//*
//****************************************************************************


#include <iostream>
#include <string>
#include <iostream>
#include <iomanip>

#include "Conventions/Units.h"

#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

// #define _debug_

using namespace std;
using namespace genie;

Int_t set_top_volume(string name);
void  get_mass      (double length_unit, double density_unit);

string gFileName;

//____________________________________________________________________________
void get_target_mass ( 
  string geometry_file, string topvolname = "", 
  double length_unit=units::cm, double density_unit=units::g_cm3, 
  Bool_t checkoverlaps = kFALSE)
{
  PDGLibrary* pdglib = PDGLibrary::Instance(); // get message out of the way

   gFileName = geometry_file;
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
   //tables of Z and A
   const Int_t lcin_Z = 150;
   const Int_t lcin_A = 300;

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

   //check if material index is unique
   Int_t * checkindex = new Int_t[max_idx+1];
   for( Int_t i = 0; i<max_idx+1; i++ ) checkindex[i] = 0;
   for( Int_t imat = 0; imat < nmat; imat++ )
   {
      if( !checkindex[imat] ) checkindex[imat] = 1;
      else 
      {
         cout << "material index is not unique" << endl;
        return;
      }
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
   TGeoNode *node;
   NodeIter.SetType(0); // include  all daughters

   Double_t * volume = new Double_t[max_idx+1];
   Double_t * mass   = new Double_t[max_idx+1];

   for( Int_t i = 0; i<max_idx+1; i++ ){ volume[i]=0.; mass[i]=0.; } // IMPORTANT! force empty arrays, allows repated calls without ending ROOT session

   volume[ topvol->GetMaterial()->GetIndex() ] = topvol->Capacity() * volume_unit_to_SI; //iterator does not include topvolume  

   while ( (node=NodeIter()) )
   {
      Int_t momidx = node->GetMotherVolume()->GetMaterial()->GetIndex() ;
      Int_t idx    = node->GetVolume()      ->GetMaterial()->GetIndex() ;

      Double_t node_vol = node->GetVolume()->Capacity() * volume_unit_to_SI;

      volume[ momidx ] -= node_vol; //substract subvolume from mother
      volume[ idx    ] += node_vol;
   }


   Double_t larr_MassIsotopes[lcin_Z][lcin_A] = {0.}; //[Z][A], no map in pure ROOT
   Double_t larr_VolumeIsotopes[lcin_Z][lcin_A] = {0.}; //[Z][A], no map in pure ROOT

   for( Int_t i=0; i<gGeoManager->GetListOfMaterials()->GetEntries(); i++ )
   {
      TGeoMaterial *lgeo_Mat = gGeoManager->GetMaterial(i);
      Int_t    idx     = gGeoManager->GetMaterial(i)->GetIndex();

      if( lgeo_Mat->IsMixture() )
      {
         TGeoMixture * lgeo_Mix = dynamic_cast <TGeoMixture*> ( lgeo_Mat );
         Int_t lint_Nelements = lgeo_Mix->GetNelements();

         for ( Int_t j=0; j<lint_Nelements; j++) 
         {
            Int_t lint_Z = TMath::Nint( (Double_t) lgeo_Mix->GetZmixt()[j] );
            Int_t lint_A = TMath::Nint( (Double_t) lgeo_Mix->GetAmixt()[j] );
            Double_t ldou_Fraction = lgeo_Mix->GetWmixt()[j];
            Double_t ldou_Density = lgeo_Mix->GetDensity() * density_unit_to_SI;

            larr_MassIsotopes[ lint_Z ][ lint_A ] += volume[idx] * ldou_Fraction * ldou_Density;
            larr_VolumeIsotopes[ lint_Z ][ lint_A ] += volume[idx] * ldou_Fraction;
         }
      }
   }

   //
   // print out volume/mass for each `material'
   //

   Double_t ldou_MinimumVolume = 1e-20;

   cout << endl
        << " Geometry: \"" <<  gFileName << "\"" << endl
        << " TopVolume: \"" << topvol->GetName() << "\"" 
        << endl;

   cout <<endl << "materials:" << endl;
   cout << setw(5) << "index"
        << setw(15) << "name"
        << setprecision(6) 
        << setw(14) << "volume (m^3)"
        << setw(14) << "mass (kg)"
        << setw(14) << "mass (%)"
        <<  endl;

   double total_mass_materials = 0;
   for( Int_t i=0; i<gGeoManager->GetListOfMaterials()->GetEntries(); i++ )
   {
     Int_t    idx     = gGeoManager->GetMaterial(i)->GetIndex();
     Double_t density = gGeoManager->GetMaterial(i)->GetDensity() * density_unit_to_SI;
     Double_t mass_material = density * volume[idx];
     if ( volume[idx] > ldou_MinimumVolume ) {
       total_mass_materials += mass_material;
     }
   }


   for( Int_t i=0; i<gGeoManager->GetListOfMaterials()->GetEntries(); i++ )
   {
      Int_t    idx     = gGeoManager->GetMaterial(i)->GetIndex();
      Double_t density = gGeoManager->GetMaterial(i)->GetDensity() * density_unit_to_SI;

      mass[idx] = density * volume[idx];

      if( volume[idx] > ldou_MinimumVolume ) {
        cout << setw(5) << i 
             << setw(15) << gGeoManager->GetMaterial(i)->GetName() 
             << setprecision(6) 
             << setw(14) << volume[idx] 
             << setw(14) << mass[idx] 
             << setw(14) << mass[idx]*100./total_mass_materials
             <<  endl;
      }
   }


   //
   // print out mass contribution for each nuclear target
   //
   PDGLibrary* pdglib = PDGLibrary::Instance();

   cout <<endl << "isotopes:" << endl;
   cout << setw(4) << "Z" 
        << setw(4) << "A"
        << setw(14) << "PDG isotope"
        << setw(6) << "      "
        << setprecision(6)
        << setw(14) << "volume (m^3)"
        << setw(14) << "mass (kg)"
        << setw(14) << "mass (%)"
        <<  endl;

   double total_mass_isotopes = 0;
   for( Int_t i=0; i<lcin_Z; i++ ) {
     for( Int_t j=0; j<lcin_A; j++ ) {
       if( larr_VolumeIsotopes[ i ][ j ] > ldou_MinimumVolume ) {
         total_mass_isotopes += larr_MassIsotopes[ i ][ j ];
       }
     }
   }

   for( Int_t i=0; i<lcin_Z; i++ )
   {
      for( Int_t j=0; j<lcin_A; j++ )
      {
         if( larr_VolumeIsotopes[ i ][ j ] > ldou_MinimumVolume ) {
           int pdgcode = 1000000000 + i*10000 + j*10;
              cout << setw(4) << i
             << setw(4)<< j
             << setw(14) << pdgcode
             << setw(6) << pdglib->Find(pdgcode)->GetName()
             << setprecision(6) 
             << setw(14) << larr_VolumeIsotopes[ i ][ j ]
             << setw(14) << larr_MassIsotopes[ i ][ j ] 
             << setw(14) << larr_MassIsotopes[ i ][ j ]*100.0/total_mass_isotopes
             <<  endl;
         }
         else if ( larr_VolumeIsotopes[ i ][ j ] < -ldou_MinimumVolume ) {
            cout << "negative volume, check geometry " << larr_VolumeIsotopes[ i ][ j ] << endl;
         }
      }
   }

   cout << endl << " mass totals: " << total_mass_materials << " " << total_mass_isotopes 
        << endl << endl;

   delete [] volume;
   delete [] mass;

}
//____________________________________________________________________________

