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

#include "Base/ELFormFactors.h"
#include "Messenger/Messenger.h"

using std::endl;
using namespace genie;

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
  this->fModel = 0;
  this->InitFormFactors();
}
//____________________________________________________________________________
ELFormFactors::ELFormFactors(const ELFormFactors & ff)
{
  this->fModel = ff.fModel;
  this->fGep   = ff.fGep;
  this->fGmp   = ff.fGmp;
  this->fGen   = ff.fGen;
  this->fGmn   = ff.fGmn;
}
//____________________________________________________________________________
void ELFormFactors::SetModel(const ELFormFactorsModelI * model)
{
  this->fModel = model;
  this->InitFormFactors();
}
//____________________________________________________________________________
void ELFormFactors::Calculate(double q2)
{
  if(!this->fModel)
  {
    LOG("ELFormFactors", pERROR)
                   << "No ELFormFactorModelI algorithm was defined!";
    this->InitFormFactors();
  }
  else {
    this->fGep = this->fModel->Gep(q2);
    this->fGmp = this->fModel->Gmp(q2);
    this->fGen = this->fModel->Gen(q2);
    this->fGmn = this->fModel->Gmn(q2);
  }
}
//____________________________________________________________________________
void ELFormFactors::InitFormFactors()
{
  this->fGep = 0.;
  this->fGmp = 0.;
  this->fGen = 0.;
  this->fGmn = 0.;
}
//____________________________________________________________________________
void ELFormFactors::Print(ostream & stream) const
{
  stream<< endl;
  stream<< "(Gep = " << this->fGep << ", Gmp = " << this->fGmp << ")" << endl;
  stream<< "(Gen = " << this->fGen << ", Gmn = " << this->fGmn << ")" << endl;
}
//____________________________________________________________________________

