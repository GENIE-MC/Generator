//____________________________________________________________________________
/*!

\class    genie::ELFormFactors

\brief    A class holding the Elastic Form Factors Ge,Gm

          This class is using the \b Strategy Pattern. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 20, 2004

*/
//____________________________________________________________________________

#include <string>

#include "Base/ELFormFactors.h"
#include "Messenger/Messenger.h"
#include "Utils/MathUtils.h"

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
void ELFormFactors::Calculate(double q2)
{
  if(!this->fModel)
  {
    LOG("ELFormFactors", pERROR)
                   << "No ELFormFactorModelI algorithm was defined!";
    this->Reset("D");
  }
  else {
    this->fGep = this->fModel->Gep(q2);
    this->fGmp = this->fModel->Gmp(q2);
    this->fGen = this->fModel->Gen(q2);
    this->fGmn = this->fModel->Gmn(q2);
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

