//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 06, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 26, 2009 - CA
   That file was added in 2.5.1 - Added RecursiveExhaust() method contributed
   by Jacek Holeczek.

*/
//____________________________________________________________________________

#include <TGeoVolume.h>
#include <TGeoNode.h>

#include "Tools/Geometry/GeoUtils.h"

//___________________________________________________________________________
void genie::utils::geometry::RecursiveExhaust(
            TGeoVolume *topvol, string volnames, bool exhaust)
{
// contributed by Jacek Holeczek
//
//
    if (!topvol) return; // just a precaution
            
    // check if the "current volume" is in the list of volumes
    if (topvol->GetName()) // non-null pointer to volume name?
    {
        const char *name = topvol->GetName();
        size_t length = strlen(name);
        if (length) // non-empty volume name?
        {
            size_t ind = 0;
            while ( (ind = volnames.find_first_of("+-", ind)) != std::string::npos )
            {
                ind += 1;
                if (ind == volnames.length()) break; // just a precaution
                if ( (!(volnames.compare(ind, length, name))) &&
                     ( ((ind + length) == volnames.length()) ||
                       (volnames[(ind + length)] == '+') ||
                       (volnames[(ind + length)] == '-') ||
                       (volnames[(ind + length)] == ' ') ) )
                {
                    // a "match" is found ... now check what to do with it
                    if ( volnames[(ind - 1)] == '+')
                        exhaust = false;
                    else if ( volnames[(ind - 1)] == '-')
                        exhaust = true;
                    break; // we are done
                }
            }
        }
    }

#if defined(DEBUG_RECURSIVE_EXHAUST)
    std::cout << topvol->GetName()
              << " <" << topvol->GetMedium()->GetName() << ">"
              << " : " << exhaust << " :";
#endif /* defined(DEBUG_RECURSIVE_EXHAUST) */
   
    // "exhaust" the "current top volume" if requested 
    if (exhaust)
    {
        static TGeoMaterial *matVacuum = ((TGeoMaterial *)0);
        static TGeoMedium *Vacuum = ((TGeoMedium *)0);
   
        if (!Vacuum)
        {  
#if defined(DEBUG_RECURSIVE_EXHAUST)
            std::cout << " Creating the Vaccum material and medium :";
#endif /* defined(DEBUG_RECURSIVE_EXHAUST) */
            // Actually ... one should check if the "Vacuum" TGeoMaterial and
            // TGeoMedium are already defined in the geometry and, if found,
            // re-use them ... but I was too lazy to implement it here, sorry.
            if (!matVacuum) matVacuum = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);
            if (matVacuum) Vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
        }
        
        if (Vacuum) topvol->SetMedium(Vacuum); // "exhaust" volume
    }
         
#if defined(DEBUG_RECURSIVE_EXHAUST)
    std::cout << " <" << topvol->GetMedium()->GetName() << ">"
              << std::endl;                                 
#endif /* defined(DEBUG_RECURSIVE_EXHAUST) */

    // proceed with all daughters of the "current volume"
    Int_t nd = topvol->GetNdaughters();
    for (Int_t i = 0; i < nd; i++)
    {
       // non-null pointer to node?
       if (topvol->GetNode(i)) {
            genie::utils::geometry::RecursiveExhaust(
              topvol->GetNode(i)->GetVolume(), volnames, exhaust);
       }
    }
}
//___________________________________________________________________________

