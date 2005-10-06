//____________________________________________________________________________
/*!

\class   genie::PathLengthList

\brief   Object to be filled with the neutrino path-length, for all detector
         geometry materials, when starting from a position x and travelling
         along the direction of the neutrino 4-momentum.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 24, 2005

*/
//____________________________________________________________________________

#include <string>
#include <iomanip>

#include <TLorentzVector.h>

#include "EVGDrivers/PathLengthList.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

using std::string;
using std::setw;
using std::setfill;
using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const PathLengthList & list)
 {
   list.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
PathLengthList::PathLengthList(const PDGCodeList & pdg_list) :
map<int, double>()
{
  PDGCodeList::const_iterator pdg_iter;

  for(pdg_iter = pdg_list.begin(); pdg_iter != pdg_list.end(); ++pdg_iter) {

     int pdgc = *pdg_iter;
     this->insert( map<int, double>::value_type(pdgc, 0.) );
  }
}
//___________________________________________________________________________
PathLengthList::~PathLengthList()
{

}
//___________________________________________________________________________
void PathLengthList::AddPathLength(int pdgc, double pl)
{
// Adds pl to the total path length for material with code = pdgc

  if (this->count(pdgc) == 1) { (*this)[pdgc] += pl; } 
  else {
     LOG("PathL", pWARN)
         << "No material with PDG code = " << pdgc << " in path length list";
  }
}
//___________________________________________________________________________
void PathLengthList::SetPathLength(int pdgc, double pl)
{
// Sets the total path length for material with code = pdgc to be pl

  if (this->count(pdgc) == 1) { (*this)[pdgc] = pl; }
  else {
     LOG("PathL", pWARN)
         << "No material with PDG code = " << pdgc << " in path length list";
  }
}
//___________________________________________________________________________
void PathLengthList::ScalePathLength(int pdgc, double scale)
{
// Scales pl for material with code = pdgc with the input scale factor

  if (this->count(pdgc) == 1) {
     double pl = (*this)[pdgc];
     pl *= scale;
     (*this)[pdgc] = pl;
  } else {
     LOG("PathL", pWARN)
         << "No material with PDG code = " << pdgc << " in path length list";
  }
}
//___________________________________________________________________________
double PathLengthList::PathLength(int pdgc) const
{
// Gets the total path length for material with code = pdgc to be pl

  if ( this->count(pdgc) == 1 ) {
     map<int, double>::const_iterator pl_iter = this->find(pdgc);
     return pl_iter->second;
  } else {
     LOG("PathL", pWARN)
         << "No material with PDG code = " << pdgc << " in path length list";
  }
  return 0;
}
//___________________________________________________________________________
void PathLengthList::SetAllToZero(void)
{
  PathLengthList::const_iterator pl_iter;

  for(pl_iter = this->begin(); pl_iter != this->end(); ++pl_iter) {
    int pdgc = pl_iter->first;
    (*this)[pdgc] = 0.;
  }
}
//___________________________________________________________________________
void PathLengthList::Print(ostream & stream) const
{
  stream << "\n[-]" << endl;

  PDGLibrary * pdglib = PDGLibrary::Instance();

  PathLengthList::const_iterator pl_iter;
  size_t nc = this->size();

  for(pl_iter = this->begin(); pl_iter != this->end(); ++pl_iter) {

    int    pdgc = pl_iter->first;
    double pl   = pl_iter->second; // path length

    TParticlePDG * p = pdglib->Find(pdgc);

    if(!p) {
      stream << " |---o ** ERR: no particle with PDG code: " << pdgc;
    } else {
      string name = p->GetName();
      stream << " |---o code: " << pdgc << " [" << setfill(' ')
         << setw(5) << name << "] " << "-----> path-length = " << pl;
    }
    if( (--nc) > 0) stream << endl;
  }
}
//___________________________________________________________________________

