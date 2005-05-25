//____________________________________________________________________________
/*!

\class   genie::PDGCodeList

\brief   A list of PDG codes

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 13, 2005

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"

using std::endl;
using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const PDGCodeList & list)
 {
   list.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
PDGCodeList::PDGCodeList() : list<int>()
{

}
//___________________________________________________________________________
PDGCodeList::PDGCodeList(size_type n) : list<int>(n)
{

}
//___________________________________________________________________________
PDGCodeList::~PDGCodeList()
{

}
//___________________________________________________________________________
void PDGCodeList::push_front(int pdg_code)
{
  bool exists = this->ExistsInPDGLibrary(pdg_code);

  if(!exists) {
    LOG("PDG", pERROR)
                   << "\n*** Can't push_front non-existent particle in list";
  } else {
    list<int>::push_front(pdg_code);
  }  
}
//___________________________________________________________________________
void PDGCodeList::push_back(int pdg_code)
{
  bool exists = this->ExistsInPDGLibrary(pdg_code);

  if(!exists) {
    LOG("PDG", pERROR)
                   << "\n*** Can't push_back non-existent particle in list";
  } else {
    list<int>::push_back(pdg_code);
  }
}
//___________________________________________________________________________
bool PDGCodeList::ExistsInPDGLibrary(int pdg_code)
{
  PDGLibrary * pdglib = PDGLibrary::Instance();

  TParticlePDG * particle = pdglib->Find(pdg_code);

  if(!particle) return false;

  return true;
}
//___________________________________________________________________________
void PDGCodeList::Print(ostream & stream) const
{
  stream << "\n[-] PDG Code List" << endl;
  
  PDGCodeList::const_iterator iter;

  for(iter = this->begin(); iter != this->end(); ++iter) {

    int pdg_code = *iter;
    stream << " |-----o code: " << pdg_code << endl;
  }
}
//___________________________________________________________________________
