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
  this->fGe    = ff.fGe;
  this->fGm    = ff.fGm;
}
//____________________________________________________________________________
void ELFormFactors::SetModel(const ELFormFactorsModelI * model)
{
  this->fModel = model;
  this->InitFormFactors();
}
//____________________________________________________________________________
void ELFormFactors::Calculate(const Interaction * interaction)
{
  if(!this->fModel)
  {
    LOG("ELFormFactors", pERROR)
                   << "No ELFormFactorModelI algorithm was defined!";
    this->InitFormFactors();
  }
  else {
    this->fGe = this->fModel->Ge(interaction);
    this->fGm = this->fModel->Gm(interaction);
  }
}
//____________________________________________________________________________
void ELFormFactors::InitFormFactors()
{
  this->fGe = 0.;
  this->fGm = 0.;
}
//____________________________________________________________________________
void ELFormFactors::Print(ostream & stream) const
{
  stream << "(Ge = " << this->fGe << ", Gm = " << this->fGm << ")" << endl;
}
//____________________________________________________________________________

