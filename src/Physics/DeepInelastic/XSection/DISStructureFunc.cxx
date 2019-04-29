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

#include "Physics/DeepInelastic/XSection/DISStructureFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/MathUtils.h"

using namespace genie;
using namespace genie::utils;

using std::endl;
using std::string;

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
  this->Reset();
}
//____________________________________________________________________________
DISStructureFunc::DISStructureFunc(const DISStructureFunc & sf)
{
  this->Copy(sf);
}
//____________________________________________________________________________
void DISStructureFunc::SetModel(const DISStructureFuncModelI * model)
{
  this->Reset();
  this->fModel = model;
}
//____________________________________________________________________________
void DISStructureFunc::Calculate(const Interaction * interaction)
{
  if(!this->fModel) {
    LOG("DISSF",pERROR)
             << "No DISStructureFuncModelI attached. Can not calculate SF's";
    this->Reset("D");
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
void DISStructureFunc::Reset(Option_t * opt)
{
// Reset the DISStructureFunc object (data & attached model). If the input
// option = D it resets the data only and not the attached model.

  this->fF1 = 0.0;
  this->fF2 = 0.0;
  this->fF3 = 0.0;
  this->fF4 = 0.0;
  this->fF5 = 0.0;
  this->fF6 = 0.0;

  string option(opt);
  if(option.find("D") == string::npos) {this->fModel = 0;}
}
//____________________________________________________________________________
void DISStructureFunc::Copy(const DISStructureFunc & sf)
{
  this->fF1 = sf.fF1;
  this->fF2 = sf.fF2;
  this->fF3 = sf.fF3;
  this->fF4 = sf.fF4;
  this->fF5 = sf.fF5;
  this->fF6 = sf.fF6;

  this->fModel = sf.fModel;
}
//____________________________________________________________________________
bool DISStructureFunc::Compare(const DISStructureFunc & sf) const
{
  bool equal =
          math::AreEqual(this->fF1, sf.fF1) &&
          math::AreEqual(this->fF2, sf.fF2) &&
          math::AreEqual(this->fF3, sf.fF3) &&
          math::AreEqual(this->fF4, sf.fF4) &&
          math::AreEqual(this->fF5, sf.fF5) &&
          math::AreEqual(this->fF6, sf.fF6);
  return equal;
}
//____________________________________________________________________________
void DISStructureFunc::Print(ostream & stream) const
{
  stream << "(F1-F6) = (" 
         << this->fF1 << ", " << this->fF2 << ", "
         << this->fF3 << ", " << this->fF4 << ", "
         << this->fF5 << ", " << this->fF6 << ")" << endl;
/*
  stream << "F1  = " << this->fF1 << endl;
  stream << "F2  = " << this->fF2 << endl;
  stream << "F3  = " << this->fF3 << endl;
  stream << "F4  = " << this->fF4 << endl;
  stream << "F5  = " << this->fF5 << endl;
  stream << "F6  = " << this->fF6 << endl;
*/
}
//____________________________________________________________________________
bool DISStructureFunc::operator == (const DISStructureFunc & sf) const
{
  return this->Compare(sf);
}
//___________________________________________________________________________
DISStructureFunc & DISStructureFunc::operator = (const DISStructureFunc & sf)
{
  this->Copy(sf);
  return (*this);
}
//___________________________________________________________________________

