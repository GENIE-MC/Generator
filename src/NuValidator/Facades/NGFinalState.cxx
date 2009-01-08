//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGFinalState

\brief    NeuGEN's Final State information

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#include "Facades/NGFinalState.h"

using std::endl;
using namespace genie::nuvld::facades;

ClassImp(NGFinalState)

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
  namespace facades {

     ostream & operator << (ostream & stream, const NGFinalState & final)
     {
       final.Print(stream);
       return stream;
     }
  }
 }
}       
//____________________________________________________________________________
NGFinalState::NGFinalState() :
_name("default")
{
  _proton  = 0;
  _neutron = 0;
  _piplus  = 0;
  _piminus = 0;
  _pizero  = 0;
}
//____________________________________________________________________________
NGFinalState::NGFinalState(const char * name) :
_name(name)
{
  _proton  = 0;
  _neutron = 0;
  _piplus  = 0;
  _piminus = 0;
  _pizero  = 0;
}
//____________________________________________________________________________
NGFinalState::NGFinalState(
               int proton, int neutron, int piplus, int piminus, int pizero) :
_proton  (proton  ),
_neutron (neutron ),
_piplus  (piplus  ),
_piminus (piminus ),
_pizero  (pizero  )
{
  _name = "default";
}
//____________________________________________________________________________
NGFinalState::~NGFinalState()
{

}
//____________________________________________________________________________
void NGFinalState::SetFinalState(
                 int proton, int neutron, int piplus, int piminus, int pizero)
{
  _proton  = proton;
  _neutron = neutron;
  _piplus  = piplus;
  _piminus = piminus;
  _pizero  = pizero;
}
//____________________________________________________________________________
void NGFinalState::Print(ostream & stream) const
{
  stream << endl;
  stream << "Protons:..........." << _proton  << endl;
  stream << "Neutrons:.........." << _neutron << endl;
  stream << "Piplus:............" << _piplus  << endl;
  stream << "Piminus:..........." << _piminus << endl;
  stream << "Pizero:............" << _pizero  << endl;
}
//____________________________________________________________________________


