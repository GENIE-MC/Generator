//_____________________________________________________________________________
/*!

\class    genie::nuvld::BrowserSingleton

\brief

\author  Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created January 20, 2004
*/
//_____________________________________________________________________________

#include "NuVldGUI/BrowserSingleton.h"

using namespace genie::nuvld;

ClassImp(BrowserSingleton)

//_____________________________________________________________________________
BrowserSingleton * BrowserSingleton::_self = 0;
//_____________________________________________________________________________
BrowserSingleton * BrowserSingleton::Instance()
{
  if(_self == 0) _self = new BrowserSingleton;

  return _self;
}
//_____________________________________________________________________________
BrowserSingleton::BrowserSingleton()
{
  _self = 0;

  _e_canvas  = 0;
  _text_edit = 0;
}
//_____________________________________________________________________________
BrowserSingleton::~BrowserSingleton()
{

}
//_____________________________________________________________________________
