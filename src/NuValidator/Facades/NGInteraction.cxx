//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGInteraction

\brief    NeuGEN's Interaction

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#include <iostream>

#include "Facades/NGInteraction.h"

using std::endl;
using namespace genie::nuvld::facades;

ClassImp(NGInteraction)

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
  namespace facades {

     ostream & operator << (ostream & stream, const NGInteraction & in)
     {
       in.Print(stream);
       return stream;
     }
  }
 }
}        
//____________________________________________________________________________
NGInteraction::NGInteraction() :
_name("default")
{

}
//____________________________________________________________________________
NGInteraction::NGInteraction(const char * name) :
_name(name)
{

}
//____________________________________________________________________________
NGInteraction::NGInteraction(const NGInteraction * inter) 
{
  _name       = inter->_name;
  _flavor     = inter->_flavor;
  _nucleus    = inter->_nucleus;
  _ccnc       = inter->_ccnc;
  _init_state = inter->_init_state;
}
//____________________________________________________________________________
NGInteraction::NGInteraction(
                  NGFlavor_t f, NGNucleus_t n, NGCcNc_t c, NGInitState_t in) :
_flavor(f),
_nucleus(n),
_ccnc(c),
_init_state(in)
{

}
//____________________________________________________________________________
NGInteraction::~NGInteraction()
{

}
//____________________________________________________________________________
int NGInteraction::GetProcess(void) const
{
  int icode = 0;
  int pa    = 0;
  int im    = _init_state;
  int it    = _flavor;
  int cn    = _ccnc;

  makestate_(&pa, &pa, &pa, &pa, &pa, &it, &im, &cn, &icode);

  return icode;
}
//____________________________________________________________________________
int NGInteraction::GetProcess(NGFinalState * final) const
{
  int proton =  final->GetProton();
  int neutron = final->GetNeutron();
  int piplus =  final->GetPiplus();
  int piminus = final->GetPiminus();
  int pizero =  final->GetPizero();
  
  int icode = 0;
  int im    = _init_state;
  int it    = _flavor;
  int cn    = _ccnc;

  makestate_(&piplus, &piminus, &pizero, &proton, &neutron,
             &it, &im, &cn, &icode);

  return icode;
}
//____________________________________________________________________________
void NGInteraction::Print(ostream & stream) const
{
  stream << "Interaction:..." << _name                              << endl;
  stream << "Flavor:........" << NGFlavor::AsString(_flavor)        << endl;
  stream << "Nucleus:......." << NGNucleus::AsString(_nucleus)      << endl;
  stream << "CCNC:.........." << NGCcNc::AsString(_ccnc)            << endl;
  stream << "InitState:....." << NGInitState::AsString(_init_state) << endl;  
}
//____________________________________________________________________________

