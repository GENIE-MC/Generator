//____________________________________________________________________________
/*!

\program gtestNaturalIsotopes

\brief   A simple test program to test the natural isotopes table

\author  Jim Dobson <j.dobson07@imperial.ac.uk>
         Imperial College London

\created May 30, 2008

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <iomanip>
#include <string>

#include "Messenger/Messenger.h"
#include "Utils/NaturalIsotopes.h"
#include "PDG/PDGLibrary.h"
#include "TParticlePDG.h"
#include "TMath.h"

using namespace genie;
using namespace std;

int getA(int pdg);
int setA(int pdg, int a);
std::string PDGcheck(int pdg);

int main(int argc, char** argv)
{
  bool docheckpdg = false;
  bool calcavg    = false;

  for (int iarg=1; iarg < argc; ++iarg) {
    // argv[0] is program name, skip it
    //  std::cout << "[" << iarg << "] = \"" 
    //     << argv[iarg] << "\"" << std::endl;
    std::string argvstr = std::string(argv[iarg]);
    if      ( argvstr == "--check-pdg" ) docheckpdg = true;
    else if ( argvstr == "--avg-a"     ) calcavg = true;
    else {
      cout << "Usage: " << argv[0] << " [--check-pdg] [--avg-a]" << endl;
      exit(1);
    } 
  }

  LOG("test", pNOTICE) 
     << "Testing the NaturalIsotopes utility";  

  /* PDGLibrary * pdglib  =  */ PDGLibrary::Instance();  // get msg out
  NaturalIsotopes * ni = NaturalIsotopes::Instance();


  if ( docheckpdg ) {
    LOG("test", pNOTICE) 
      << "Check on PDG in PDGLibrary: [PDG,PDG-proton,PDG-neutron]";
  }

  for(int Z = 1; Z < 104; Z++){
    int nel = ni->NElements(Z);
    LOG("test", pNOTICE) << "** Z = " << Z
                         << "   Number of elements = " << nel;
    
    int pdg;
    double abund, avg_a = 0;
    std::string extra;
    for(int n = 0; n < nel; n++){
      const NaturalIsotopeElementData * el = ni->ElementData(Z,n);   
      pdg = el->PdgCode();
      abund = el->Abundance();
      avg_a += double(getA(pdg)) * abund;
      extra = ( docheckpdg ) ? PDGcheck(pdg) : "";
      LOG("test", pNOTICE) 
        << "   -- Element: " << n 
        << ", PdgCode = " << pdg
        << extra
        << ", Abundance = " << abund;
    } // n
    if ( calcavg && ( nel > 1 ) ) {
      avg_a /= 100.; // abundances were in %
      pdg = setA(pdg,TMath::Nint(avg_a));
      extra = ( docheckpdg ) ? PDGcheck(pdg) : "";
      LOG("test", pNOTICE)
        << "                 "
        << " PdgCode = " << pdg
        << extra
        << ", average <A> " << avg_a ;
    }
  } // Z  
   
  return 0;
}

int getA(int pdg)
{ return (int(pdg/10)%1000); }

int setA(int pdg, int a)
{
  return int(pdg/10000)*10000 + 10*a;
}

std::string PDGcheck(int pdg)
{
  const int minus_p = 10010;  // Z-1, A-1
  const int minus_n =    10;  // A-1

  PDGLibrary * pdglib  = PDGLibrary::Instance();
  int pdg_minus_p = pdg - minus_p;
  int pdg_minus_n = pdg - minus_n;
  const TParticlePDG * part         = pdglib->Find(pdg);
  const TParticlePDG * part_minus_p = pdglib->Find(pdg_minus_p);
  const TParticlePDG * part_minus_n = pdglib->Find(pdg_minus_n);
  string codeck = " [";
  codeck += ( (part) ? "ok" : "NO" );
  if ( pdg != 1000010010 ) {
  codeck += ",";
    codeck += ( (part_minus_p) ? "ok" : "NO" );
    codeck += ",";
    codeck += ( (part_minus_n) ? "ok" : "NO" );
    codeck += "]";
  } else { // special case for H1
    codeck += ",--,--]";
  }

  return codeck;
}
