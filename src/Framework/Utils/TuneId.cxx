//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda <Marco.Roda \at liverpool.ac.uk>
         University of Liverpool

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

//#include <sstream>

#include "Framework/Utils/TuneId.h"

//using std::ostringstream;
using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const TuneId & id)
  {
    id.Print(stream);
    return stream;
  }
  //..........................................................................
  bool operator == (const TuneId & id1, const TuneId & id2)
  {
    return id1.Compare(id2);
  }
  //..........................................................................
  bool operator != (const TuneId & id1, const TuneId & id2)
  {
    return !id1.Compare(id2);
  }
}
//____________________________________________________________________________
TuneId::TuneId() :
TObject()
{
  this->Init();
}
//____________________________________________________________________________
TuneId::TuneId(string id_str) :
TObject()
{
  this->Decode(id_str);
}
//____________________________________________________________________________
TuneId::TuneId(const TuneId & id) :
TObject()
{
  this->Copy(id);
}
//____________________________________________________________________________
TuneId::~TuneId()
{

}
//____________________________________________________________________________
void TuneId::Decode(string id_str)
{
  this->fName = id_str;
/*
  this->fPrefix          = 
  this->fYear            = 
  this->fMajorModelId    = 
  this->fMinorModelId    = 
  this->fTunedParamSetId = 
  this->fFitDataSetId    = 
*/
}
//____________________________________________________________________________
void TuneId::Copy(const TuneId & id)
{
  this->fName            = id.Name();
  this->fPrefix          = id.Prefix();
  this->fYear            = id.Year();
  this->fMajorModelId    = id.MajorModelId();
  this->fMinorModelId    = id.MinorModelId();
  this->fTunedParamSetId = id.TunedParamSetId();
  this->fFitDataSetId    = id.FitDataSetId();
}
//____________________________________________________________________________
bool TuneId::Compare(const TuneId & id) const
{
  return (this->Name() == id.Name());
}
//____________________________________________________________________________
void TuneId::Print(ostream & stream) const
{
  stream << " GENIE tune: " << endl;
  stream << " - Prefix ............... : " << this->Prefix()          << endl; 
  stream << " - Year ................. : " << this->Year()            << endl; 
  stream << " - Major model ID ....... : " << this->MajorModelId()    << endl; 
  stream << " - Minor model ID ....... : " << this->MinorModelId()    << endl; 
  stream << " - Tuned param set ID ... : " << this->TunedParamSetId() << endl;
  stream << " - Fit dataset ID ....... : " << this->FitDataSetId()    << endl;
}
//____________________________________________________________________________
void TuneId::Init(void)
{
  this->fName            = "";
  this->fPrefix          = "";
  this->fYear            = "";
  this->fMajorModelId    = "";
  this->fMinorModelId    = "";
  this->fTunedParamSetId = "";
  this->fFitDataSetId    = "";
}
//____________________________________________________________________________
