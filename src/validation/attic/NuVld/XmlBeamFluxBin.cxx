//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 01, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/XmlBeamFluxBin.h"

namespace genie {
namespace nuvld {
  
//______________________________________________________________________________________
ostream & operator << (ostream & stream, const XmlBeamFluxBin & bin)
{
  stream << "[Emin = " << bin._Emin << ", Emax = " << bin._Emax << "] " 
         << " - <E> = " << bin._E << " ---> flux = " << bin._Flux 
         << " (+" << bin._Flux_perr << ", -" << bin._Flux_nerr << ")" << endl;
	 
  return stream;
}
//______________________________________________________________________________________
XmlBeamFluxBin::XmlBeamFluxBin() :
_E(""),
_Emin(""),
_Emax(""),
_Flux(""),
_Flux_perr(""),
_Flux_nerr("")
{

}
//______________________________________________________________________________________
XmlBeamFluxBin::XmlBeamFluxBin(
                 string E, string Emin, string Emax, string F, string dFp, string dFn) :
_E(E),
_Emin(Emin),
_Emax(Emax),
_Flux(F),
_Flux_perr(dFp),
_Flux_nerr(dFn)
{

}
//______________________________________________________________________________________
XmlBeamFluxBin::XmlBeamFluxBin(const XmlBeamFluxBin & /*bin*/)
{

}
//______________________________________________________________________________________

} // nuvld namespace
} // genie namespace
