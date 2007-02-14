//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - November 20, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <algorithm>

#include "EVGCore/GVldContext.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"

using std::count;

using namespace genie;

//___________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const GVldContext & vldc)
 {
   vldc.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
GVldContext::GVldContext()
{
  this->Init();
}
//___________________________________________________________________________
GVldContext::~GVldContext()
{

}
//___________________________________________________________________________
void GVldContext::Decode(string encoded_vld_context)
{
//example:
// ENERGY:0-100;

  vector<string> fields = utils::str::Split(encoded_vld_context, ";");

  vector<string>::const_iterator field_iter;

  for(field_iter = fields.begin(); field_iter != fields.end(); ++field_iter){

     SLOG("VldContext", pINFO) << " ************ " << *field_iter;

     vector<string> curr_field = utils::str::Split(*field_iter, ":");

     assert(curr_field.size() == 2);

     string name   = curr_field[0];
     string values = curr_field[1];

     //-- Make lowercase/uppercase irrelevant
     name = utils::str::ToUpper(name);

     //-- send the string to an appropriate decoder
     if (name.find("ENERGY") != string::npos) DecodeENERGY (values);
     else {
       SLOG("VldContext", pWARN)
                      << "**** Unknown field named: "
                                        << name << " in encoded vld-context";
     }
  }
}
//___________________________________________________________________________
void GVldContext::DecodeENERGY(string encoded_energy)
{
  SLOG("VldContext", pDEBUG) << "Decoding ENERGY: " << encoded_energy;

  vector<string> energy = utils::str::Split(encoded_energy, "-");

  assert (energy.size() == 2);

  fEmin = atof( energy[0].c_str() );
  fEmax = atof( energy[1].c_str() );
}
//___________________________________________________________________________
void GVldContext::Init(void)
{
  fEmin    = -1.0;
  fEmax    = -1.0;
}
//___________________________________________________________________________
void GVldContext::Print(ostream & stream) const
{
  stream << "Energy range:..." << "[" << fEmin << ", " << fEmax << "]";
  stream << "\n";
}
//___________________________________________________________________________
