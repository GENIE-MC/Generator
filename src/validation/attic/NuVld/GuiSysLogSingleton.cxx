//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/GuiSysLogSingleton.h"

using namespace genie::nuvld;

ClassImp(GuiSysLogSingleton)

//_____________________________________________________________________________
GuiSysLogSingleton * GuiSysLogSingleton::fSelf = 0;
//_____________________________________________________________________________
GuiSysLogSingleton * GuiSysLogSingleton::Instance()
{
  if(fSelf == 0) 
     fSelf = new GuiSysLogSingleton;

  return fSelf;     
}
//_____________________________________________________________________________
GuiSysLogSingleton::GuiSysLogSingleton()
{
  fSelf = 0;

  fStatusBar   = 0;
  fProgressBar = 0;
  fLog         = 0;  
}
//_____________________________________________________________________________
GuiSysLogSingleton::~GuiSysLogSingleton()
{

}
//_____________________________________________________________________________

