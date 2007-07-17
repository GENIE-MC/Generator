//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <algorithm>
#include <iomanip>
#include <string>

#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"

using std::setw;
using std::setfill;
using std::endl;
using std::string;

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
PDGCodeList::PDGCodeList(bool allowdup) : 
vector<int>()
{
  fAllowDuplicateEntries = allowdup;
}
//___________________________________________________________________________
PDGCodeList::PDGCodeList(size_type n, bool allowdup) : 
vector<int>(n)
{
  fAllowDuplicateEntries = allowdup;
}
//___________________________________________________________________________
PDGCodeList::PDGCodeList(const PDGCodeList & list) :
vector<int>()
{
  this->Copy(list);
}
//___________________________________________________________________________
PDGCodeList::~PDGCodeList()
{

}
//___________________________________________________________________________
void PDGCodeList::push_back(int pdg_code)
{
  if(this->CheckPDGCode(pdg_code)) vector<int>::push_back(pdg_code);
}
//___________________________________________________________________________
void PDGCodeList::insert(iterator pos, size_type n, const int& pdg_code)
{
  if(this->CheckPDGCode(pdg_code)) {
     if(n>1) n = 1;
     vector<int>::insert(pos,n,pdg_code);
  }
}
//___________________________________________________________________________
bool PDGCodeList::CheckPDGCode(int pdg_code)
{
// check whether the PDG code can be inserted

  bool exists = this->ExistsInPDGLibrary(pdg_code);
  if(!exists) {
    LOG("PDG", pERROR)
           << "Can't add non-existent particle [pdgc = " << pdg_code << "]";
    return false;
  }

  if(!fAllowDuplicateEntries) {
    bool added = this->ExistsInPDGCodeList(pdg_code);
    if(added) {
      LOG("PDG", pDEBUG)
                << "Particle [pdgc = " << pdg_code << "] was already added";
      return false;
    }
  }
  return true;
}
//___________________________________________________________________________
bool PDGCodeList::ExistsInPDGLibrary(int pdg_code)
{
// check whether the PDG code is a valid one (exists in PDGLibrary)

  PDGLibrary * pdglib = PDGLibrary::Instance();
  TParticlePDG * particle = pdglib->Find(pdg_code);
  if(!particle) return false;
  return true;
}
//___________________________________________________________________________
bool PDGCodeList::ExistsInPDGCodeList(int pdg_code)
{
// check whether the PDG code already exists in the list

  int n = count(this->begin(), this->end(), pdg_code);
  if(n!=0) return true;
  return false;
}
//___________________________________________________________________________
void PDGCodeList::Print(ostream & stream) const
{
  stream << "\n[-]" << endl;

  PDGLibrary * pdglib = PDGLibrary::Instance();

  PDGCodeList::const_iterator iter;
  size_t nc = this->size();

  for(iter = this->begin(); iter != this->end(); ++iter) {
    int pdg_code = *iter;
    TParticlePDG * p = pdglib->Find(pdg_code);

    if(!p) {
      stream << " |---o ** ERR: no particle with PDG code: " << pdg_code;
    } else {
      string name = p->GetName();
      stream << " |---o code: " << pdg_code
                       << " [" << setfill(' ') << setw(5) << name << "]";
    }
    if( (--nc) > 0) stream << endl;
  }
}
//___________________________________________________________________________
void PDGCodeList::Copy(const PDGCodeList & list)
{
  this->clear();

  PDGCodeList::const_iterator iter;
  for(iter = list.begin(); iter != list.end(); ++iter) {
    int code = *iter;
    this->push_back(code);
  }

  fAllowDuplicateEntries = list.fAllowDuplicateEntries;
}
//___________________________________________________________________________
PDGCodeList & PDGCodeList::operator = (const PDGCodeList & list)
{
  this->Copy(list);
  return (*this);
}
//___________________________________________________________________________

