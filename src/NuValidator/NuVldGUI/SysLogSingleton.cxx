//_____________________________________________________________________________
/*!

\class    genie::nuvld::SysLogSingleton

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "NuVldGUI/SysLogSingleton.h"

using namespace genie::nuvld;

ClassImp(SysLogSingleton)

//_____________________________________________________________________________
SysLogSingleton * SysLogSingleton::fSelf = 0;
//_____________________________________________________________________________
SysLogSingleton * SysLogSingleton::Instance()
{
  if(fSelf == 0) 
     fSelf = new SysLogSingleton;

  return fSelf;     
}
//_____________________________________________________________________________
SysLogSingleton::SysLogSingleton()
{
  fSelf = 0;

  fStatusBar   = 0;
  fProgressBar = 0;
  fLog         = 0;  
}
//_____________________________________________________________________________
SysLogSingleton::~SysLogSingleton()
{

}
//_____________________________________________________________________________

