//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 20, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/GuiBrowserSingleton.h"

using namespace genie::nuvld;

ClassImp(GuiBrowserSingleton)

//_____________________________________________________________________________
GuiBrowserSingleton * GuiBrowserSingleton::_self = 0;
//_____________________________________________________________________________
GuiBrowserSingleton * GuiBrowserSingleton::Instance()
{
  if(_self == 0) _self = new GuiBrowserSingleton;

  return _self;
}
//_____________________________________________________________________________
GuiBrowserSingleton::GuiBrowserSingleton()
{
  _self = 0;

  _e_canvas  = 0;
  _text_edit = 0;
}
//_____________________________________________________________________________
GuiBrowserSingleton::~GuiBrowserSingleton()
{

}
//_____________________________________________________________________________
