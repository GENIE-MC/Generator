//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   Moved into the ElFF package from its previous location               

*/
//____________________________________________________________________________

#include <string>

#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
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
  ostream & operator << (ostream & stream, const ELFormFactors & ff)
  {
     ff.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
ELFormFactors::ELFormFactors()
{
  this->Reset();
}
//____________________________________________________________________________
ELFormFactors::ELFormFactors(const ELFormFactors & ff)
{
  this->Copy(ff);
}
//____________________________________________________________________________
void ELFormFactors::SetModel(const ELFormFactorsModelI * model)
{
  this->Reset();
  this->fModel = model;
}
//____________________________________________________________________________
void ELFormFactors::Calculate(const Interaction * interaction)
{
  if(!this->fModel)
  {
    LOG("ELFormFactors", pERROR)
                   << "No ELFormFactorModelI algorithm was defined!";
    this->Reset("D");
  }
  else {
    this->fGep = this->fModel->Gep(interaction);
    this->fGmp = this->fModel->Gmp(interaction);
    this->fGen = this->fModel->Gen(interaction);
    this->fGmn = this->fModel->Gmn(interaction);
  }
}
//____________________________________________________________________________
void ELFormFactors::Reset(Option_t * opt)
{
// Reset the ELFormFactors object (data & attached model). If the input
// option = D it resets the data only and not the attached model.

  this->fGep = 0.;
  this->fGmp = 0.;
  this->fGen = 0.;
  this->fGmn = 0.;

  string option(opt);
  if(option.find("D") == string::npos) {this->fModel = 0;}
}
//____________________________________________________________________________
void ELFormFactors::Copy(const ELFormFactors & ff)
{
  this->fModel = ff.fModel;
  this->fGep   = ff.fGep;
  this->fGmp   = ff.fGmp;
  this->fGen   = ff.fGen;
  this->fGmn   = ff.fGmn;
}
//____________________________________________________________________________
bool ELFormFactors::Compare(const ELFormFactors & ff) const
{
  bool equal =
          math::AreEqual(this->fGep, ff.fGep) &&
          math::AreEqual(this->fGmp, ff.fGmp) &&
          math::AreEqual(this->fGen, ff.fGen) &&
          math::AreEqual(this->fGmn, ff.fGmn);
  return equal;
}
//____________________________________________________________________________
void ELFormFactors::Print(ostream & stream) const
{
  stream<< endl;
  stream<< "(Gep = " << this->fGep << ", Gmp = " << this->fGmp << ")" << endl;
  stream<< "(Gen = " << this->fGen << ", Gmn = " << this->fGmn << ")" << endl;
}
//____________________________________________________________________________
bool ELFormFactors::operator == (const ELFormFactors & ff) const
{
  return this->Compare(ff);
}
//___________________________________________________________________________
ELFormFactors & ELFormFactors::operator = (const ELFormFactors & ff)
{
  this->Copy(ff);
  return (*this);
}
//___________________________________________________________________________

