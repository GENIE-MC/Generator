//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

\author   Aaron Meyer <asmeyer2012 \at uchicago.edu>
          based off AxialFormFactor by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <string>

#include "Physics/QuasiElastic/XSection/AxialFormFactor.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/MathUtils.h"

using std::endl;
using std::string;

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const AxialFormFactor & ff)
  {
     ff.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
AxialFormFactor::AxialFormFactor()
{
  this->Reset();
}
//____________________________________________________________________________
AxialFormFactor::AxialFormFactor(const AxialFormFactor & ff)
{
  this->Copy(ff);
}
//____________________________________________________________________________
void AxialFormFactor::SetModel(const AxialFormFactorModelI * model)
{
  this->Reset();
  this->fModel = model;
}
//____________________________________________________________________________
void AxialFormFactor::Calculate(const Interaction * interaction)
{
  if(!this->fModel)
  {
    LOG("AxialFormFactor", pERROR)
                   << "No AxialFormFactorModelI algorithm was defined!";
    this->Reset("D");
  }
  else {
    this->fFA = this->fModel->FA(interaction);
  }
}
//____________________________________________________________________________
void AxialFormFactor::Reset(Option_t * opt)
{
// Reset the AxialFormFactor object (data & attached model). If the input
// option = D it resets the data only and not the attached model.

  this->fFA = 0.;

  string option(opt);
  if(option.find("D") == string::npos) {this->fModel = 0;}
}
//____________________________________________________________________________
void AxialFormFactor::Copy(const AxialFormFactor & ff)
{
  this->fModel = ff.fModel;
  this->fFA    = ff.fFA;
}
//____________________________________________________________________________
bool AxialFormFactor::Compare(const AxialFormFactor & ff) const
{
  return math::AreEqual(this->fFA, ff.fFA);
}
//____________________________________________________________________________
void AxialFormFactor::Print(ostream & stream) const
{
  stream<< endl;
  stream<< "(FA  = " << this->fFA << ") "<< endl;
}
//____________________________________________________________________________
bool AxialFormFactor::operator == (const AxialFormFactor & ff) const
{
  return this->Compare(ff);
}
//___________________________________________________________________________
AxialFormFactor & AxialFormFactor::operator = (const AxialFormFactor & ff)
{
  this->Copy(ff);
  return (*this);
}
//___________________________________________________________________________

