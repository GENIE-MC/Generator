//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <string>

#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/MathUtils.h"

using namespace genie;
using namespace genie::utils;

using std::endl;
using std::string;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const QELFormFactors & ff)
  {
     ff.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
QELFormFactors::QELFormFactors()
{
  this->Reset();
}
//____________________________________________________________________________
QELFormFactors::QELFormFactors(const QELFormFactors & form_factors)
{
  this->Copy(form_factors);
}
//____________________________________________________________________________
void QELFormFactors::SetModel(const QELFormFactorsModelI * model)
{
  this->Reset();
  this->fModel = model;
}
//____________________________________________________________________________
void QELFormFactors::Calculate(const Interaction * interaction)
{
  if(!this->fModel) {
    LOG("QELFF",pERROR)
             << "No QELFormFactorsModelI attached. Can not calculate FF's";
    this->Reset("D");
    return;
  }

  this -> fF1V   = fModel -> F1V   (interaction);
  this -> fxiF2V = fModel -> xiF2V (interaction);
  this -> fFA    = fModel -> FA    (interaction);
  this -> fFp    = fModel -> Fp    (interaction);
}
//____________________________________________________________________________
void QELFormFactors::Reset(Option_t * opt)
{
// Reset the QELFormFactors object (data & attached model). If the input
// option = D it resets the data only and not the attached model.

  this->fF1V   = 0;
  this->fxiF2V = 0;
  this->fFA    = 0;
  this->fFp    = 0;

  string option(opt);
  if(option.find("D") == string::npos) {this->fModel = 0;}
}
//____________________________________________________________________________
void QELFormFactors::Copy(const QELFormFactors & ff)
{
  this->fModel = ff.fModel;

  this->fF1V   = ff.fF1V;
  this->fxiF2V = ff.fxiF2V;
  this->fFA    = ff.fFA;
  this->fFp    = ff.fFp;
}
//____________________________________________________________________________
bool QELFormFactors::Compare(const QELFormFactors & ff) const
{
  bool equal =
          math::AreEqual(this->fF1V,   ff.fF1V)   &&
          math::AreEqual(this->fxiF2V, ff.fxiF2V) &&
          math::AreEqual(this->fFA,    ff.fFA)    &&
          math::AreEqual(this->fFp,    ff.fFp);
  return equal;
}
//____________________________________________________________________________
void QELFormFactors::Print(ostream & stream) const
{
  stream << endl;
  stream << "F1V    = " << this->fF1V   << endl;
  stream << "xi*F2V = " << this->fxiF2V << endl;
  stream << "FA     = " << this->fFA    << endl;
  stream << "Fp     = " << this->fFp    << endl;
}
//____________________________________________________________________________
bool QELFormFactors::operator == (const QELFormFactors & ff) const
{
  return this->Compare(ff);
}
//___________________________________________________________________________
QELFormFactors & QELFormFactors::operator = (const QELFormFactors & ff)
{
  this->Copy(ff);
  return (*this);
}
//___________________________________________________________________________


