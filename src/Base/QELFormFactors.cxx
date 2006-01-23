//____________________________________________________________________________
/*!

\class    genie::QELFormFactors

\brief    A class holding Quasi Elastic (QEL) Form Factors.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          QELFormFactorsModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 20, 2004

*/
//____________________________________________________________________________

#include <string>

#include "Base/QELFormFactors.h"
#include "Messenger/Messenger.h"
#include "Utils/MathUtils.h"

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


