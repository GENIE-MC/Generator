//____________________________________________________________________________
/*!

\program gtestNaturalIsotopes

\brief   A simple test program to test the natural isotopes table

\author  Jim Dobson <j.dobson07@imperial.ac.uk>
         Imperial College London

\created May 30, 2008

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "Utils/NaturalIsotopes.h"

using namespace genie;

int main()
{
  LOG("Main", pNOTICE) 
     << "Testing the NaturalIsotopes utility";  

  NaturalIsotopes * ni = NaturalIsotopes::Instance();

  for(int Z = 1; Z < 104; Z++){
    int nel = ni->NElements(Z);
    LOG("Main", pNOTICE) << "** Z = " << Z;
    LOG("Main", pNOTICE) << "   Number of elements = " << nel;
                       
    for(int n = 0; n < nel; n++){
      const NaturalIsotopeElementData * el = ni->ElementData(Z,n);   
      LOG("Main", pNOTICE) 
        << "   -- Element: " << n 
        << ", PdgCode = " << el->PdgCode() 
        << ", Abundance = " << el->Abundance();
    } // n
  } // Z  
   
  return 0;
}


