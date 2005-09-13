//____________________________________________________________________________
/*!

\class    genie::DISStructureFunc

\brief    A class holding Deep Inelastic Scattering (DIS) Form Factors
          (invariant structure funstions)

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          DISStructureFuncModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#include "Base/DISStructureFunc.h"
#include "Messenger/Messenger.h"

using namespace genie;

using std::endl;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const DISStructureFunc & ff)
  {
     ff.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
DISStructureFunc::DISStructureFunc()
{
  this->fModel = 0;
  this->InitFormFactors();
}
//____________________________________________________________________________
DISStructureFunc::DISStructureFunc(const DISStructureFunc & form_factors)
{
  this->fModel = form_factors.fModel;
  this->fF1    = form_factors.fF1;
  this->fF2    = form_factors.fF2;
  this->fF3    = form_factors.fF3;
  this->fF4    = form_factors.fF4;
  this->fF5    = form_factors.fF5;
  this->fF6    = form_factors.fF6;
}
//____________________________________________________________________________
void DISStructureFunc::SetModel(const DISStructureFuncModelI * model)
{
  this->InitFormFactors();

  this->fModel = model;
}
//____________________________________________________________________________
void DISStructureFunc::Calculate(const Interaction * interaction)
{
  if(!this->fModel) {
    LOG("DISSF",pERROR)
             << "No DISStructureFuncModelI attached. Can not calculate SF's";
    return;
  }

  fModel->Calculate(interaction);

  this->fF1 = fModel->F1();
  this->fF2 = fModel->F2();
  this->fF3 = fModel->F3();
  this->fF4 = fModel->F4();
  this->fF5 = fModel->F5();
  this->fF6 = fModel->F6();
}
//____________________________________________________________________________
void DISStructureFunc::InitFormFactors(void)
{
  this->fF1 = 0.0;
  this->fF2 = 0.0;
  this->fF3 = 0.0;
  this->fF4 = 0.0;
  this->fF5 = 0.0;
  this->fF6 = 0.0;
}
//____________________________________________________________________________
void DISStructureFunc::Print(ostream & stream) const
{
  stream << "F1  = " << this->fF1 << endl;
  stream << "F2  = " << this->fF2 << endl;
  stream << "F3  = " << this->fF3 << endl;
  stream << "F4  = " << this->fF4 << endl;
  stream << "F5  = " << this->fF5 << endl;
  stream << "F6  = " << this->fF6 << endl;
}
//____________________________________________________________________________

