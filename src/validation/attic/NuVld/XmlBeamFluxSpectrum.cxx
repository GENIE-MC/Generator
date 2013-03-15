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

#include "ValidationTools/NuVld/XmlBeamFluxSpectrum.h"

namespace genie {
namespace nuvld {
  
//______________________________________________________________________________________
ostream & operator << (ostream & stream, const XmlBeamFluxSpectrum & spectrum)
{
  stream << endl;
  stream << "----------------- Printing beam flux spectrum -------------------" << endl;
  
  stream << " --> flux units   = " << spectrum._flux_units << endl;
  stream << " --> energy units = " << spectrum._E_units    << endl;
  stream << " --> energy frame = " << spectrum._E_frame    << endl;

  vector<XmlBeamFluxBin *>::iterator bin_iter;
  for(bin_iter = spectrum._spectrum->begin(); 
                         bin_iter != spectrum._spectrum->end(); ++bin_iter)
                                                                 stream << *(*bin_iter);
  return stream;
}
//______________________________________________________________________________________
XmlBeamFluxSpectrum::XmlBeamFluxSpectrum()
{
  _flux_units = "cm-2 sec-1 sr-1 GeV-1";
  _E_units    = "GeV";
  _E_frame    = "lab";
  
  _spectrum = new vector<XmlBeamFluxBin *>;
}
//______________________________________________________________________________________
XmlBeamFluxSpectrum::XmlBeamFluxSpectrum(const XmlBeamFluxSpectrum & /*spectrum*/)
{

}
//______________________________________________________________________________________
void XmlBeamFluxSpectrum::Add(XmlBeamFluxBin * bin)
{
  _spectrum->push_back(bin);

  sort( _spectrum->begin(), _spectrum->end(), less_energy() );
}
//______________________________________________________________________________________

} // nuvld namespace
} // genie namespace
