//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 20, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <algorithm>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/GVldContext.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"

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
//Example:
// energy:0-100;

  string vldc = utils::str::ToUpper(encoded_vld_context);

  // set defauts for missing entries
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->CommonList("Param", "Validation");
  
  if(vldc.find("ENERGY") == string::npos) {
    fEmin = gc->GetDouble("GVLD-Emin");
    fEmax = gc->GetDouble("GVLD-Emax");
  }

  LOG("VldContext", pDEBUG) << "Validity context: " << vldc;

  vector<string> fields = utils::str::Split(vldc, ";");
  if(fields.size()==0) return;

  vector<string>::const_iterator field_iter;

  for(field_iter = fields.begin(); field_iter != fields.end(); ++field_iter){

     string curr_field = *field_iter;
     SLOG("VldContext", pINFO) << " ************ " << curr_field;
     if(curr_field.size()==0) continue;

     vector<string> curr_fieldv = utils::str::Split(curr_field, ":");
     assert(curr_fieldv.size() == 2);

     string name   = curr_fieldv[0];
     string values = curr_fieldv[1];

     //-- send the string to an appropriate decoder
     if (name.find("ENERGY") != string::npos) DecodeENERGY (values);
     else {
       SLOG("VldContext", pWARN)
            << "**** Unknown field named: " << name << " in vld context";
     }
  }
}
//___________________________________________________________________________
void GVldContext::DecodeENERGY(string encoded_energy)
{
  SLOG("VldContext", pDEBUG) << "Decoding energy range: " << encoded_energy;

  vector<string> energy = utils::str::Split(encoded_energy, "-");
  assert (energy.size() == 2);
  fEmin = atof( energy[0].c_str() );
  fEmax = atof( energy[1].c_str() );
}
//___________________________________________________________________________
void GVldContext::Init(void)
{
  fEmin = -1.0;
  fEmax = -1.0;
}
//___________________________________________________________________________
void GVldContext::Print(ostream & stream) const
{
  stream << "Energy range:..." << "[" << fEmin << ", " << fEmax << "]";
  stream << "\n";
}
//___________________________________________________________________________
