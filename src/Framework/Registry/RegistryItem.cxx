//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 06, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Oct 11, 2007 - CA
 Added 'bool IsLocal()' & 'SetLocal(bool)' implementations now required from
 the RegistryItemI interface. Indicate 'local' or 'global' status in Print()

*/
//____________________________________________________________________________

#include <string>
#include <iostream>

#include "Framework/Registry/RegistryItem.h"

using std::string;
using std::endl;

namespace genie {

//____________________________________________________________________________
template ostream & operator << 
               (ostream & stream, const RegistryItem<RgBool> & r);
template ostream & operator << 
               (ostream & stream, const RegistryItem<RgInt>  & r);
template ostream & operator << 
               (ostream & stream, const RegistryItem<RgDbl>  & r);
template ostream & operator << 
               (ostream & stream, const RegistryItem<RgStr>  & r);
template ostream & operator << 
               (ostream & stream, const RegistryItem<RgAlg>  & r);
template ostream & operator << 
               (ostream & stream, const RegistryItem<RgH1F>  & r);
template ostream & operator << 
               (ostream & stream, const RegistryItem<RgH2F>  & r);
template ostream & operator << 
               (ostream & stream, const RegistryItem<RgTree> & r);
//____________________________________________________________________________
template <typename T>
   ostream & operator << (ostream & stream, const RegistryItem<T> & rec)
{
   rec.Print(stream);
   return stream;
}
//____________________________________________________________________________
template<typename T> 
             RegistryItem<T>::RegistryItem(T item, bool locked, bool local)
{
  fItem     = item;
  fIsLocked = locked;
  fIsLocal  = local;
}
//____________________________________________________________________________
template<typename T> RegistryItem<T>::RegistryItem(const RegistryItem * ri)
{
  fItem     = ri->fItem;
  fIsLocked = ri->fIsLocked;
  fIsLocal  = ri->fIsLocal;
}
//____________________________________________________________________________
template<typename T> RegistryItem<T>::~RegistryItem()
{

}
//____________________________________________________________________________
template<> RegistryItem<RgH1F>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
template<> RegistryItem<RgH2F>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
template<> RegistryItem<RgTree>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
template<typename T> RegistryItemI * RegistryItem<T>::Clone(void) const
{
  RegistryItemI * item = new RegistryItem<T>(fItem,fIsLocked);
  return item;
}
//____________________________________________________________________________
template<> RgType_t RegistryItem<RgBool>::TypeInfo(void) const 
{ 
  return kRgBool; 
}
//____________________________________________________________________________
template<> RgType_t RegistryItem<RgInt>::TypeInfo (void) const 
{ 
  return kRgInt;  
}
//____________________________________________________________________________
template<> RgType_t RegistryItem<RgDbl>::TypeInfo (void) const 
{ 
  return kRgDbl;  
}
//____________________________________________________________________________
template<> RgType_t RegistryItem<RgStr>::TypeInfo (void) const 
{ 
  return kRgStr;  
}
//____________________________________________________________________________
template<> RgType_t RegistryItem<RgAlg>::TypeInfo (void) const 
{ 
  return kRgAlg;  
}
//____________________________________________________________________________
template<> RgType_t RegistryItem<RgH1F>::TypeInfo (void) const 
{ 
  return kRgH1F;  
}
//____________________________________________________________________________
template<> RgType_t RegistryItem<RgH2F>::TypeInfo (void) const 
{ 
  return kRgH2F;  
}
//____________________________________________________________________________
template<> RgType_t RegistryItem<RgTree>::TypeInfo(void) const 
{ 
  return kRgTree; 
}
//____________________________________________________________________________
template<typename T> void RegistryItem<T>::Print(ostream & stream) const
{
  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") 
         << " "
         << ((fIsLocal)  ? "[l]" : "[g]") 
         << " : " 
         << (fItem);
}
//____________________________________________________________________________
template<> void RegistryItem<RgAlg>::Print(ostream & stream) const
{
  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") 
         << " "
         << ((fIsLocal)  ? "[l]" : "[g]") 
         << " : " 
         << (fItem);
}
//____________________________________________________________________________
template<> void RegistryItem<RgH1F>::Print(ostream & stream) const
{
  TH1F * histo = dynamic_cast<TH1F *>(fItem);
  if(!histo) stream << "*** NULL RgH1F ***";

  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") 
         << " "
         << ((fIsLocal)  ? "[l]" : "[g]") 
         << " : " 
         << (histo ? histo->GetName() : "NULL");
}
//____________________________________________________________________________
template<> void RegistryItem<RgH2F>::Print(ostream & stream) const
{
  TH2F * histo = dynamic_cast<TH2F *>(fItem);
  if(!histo) stream << "*** NULL RgH2F ***";

  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") 
         << " "
         << ((fIsLocal)  ? "[l]" : "[g]") 
         << " : " 
         << (histo ? histo->GetName() : "NULL");
}
//____________________________________________________________________________
template<> void RegistryItem<RgTree>::Print(ostream & stream) const
{
  TTree * tree = dynamic_cast<TTree *>(fItem);
  if(!tree) stream << "*** NULL RgTree ***";

  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") 
         << " "
         << ((fIsLocal)  ? "[l]" : "[g]") 
         << " : " 
         << (tree ? tree->GetName() : "NULL");
}
//____________________________________________________________________________

// Declare all template specializations.
// They need to be added at the end of the implementation file as described at
// http://www.comeaucomputing.com/techtalk/templates/#whylinkerror
// The need for them arises from the export C++ keyword not being properly 
// implemented at most C++ compilers. See http://gcc.gnu.org/bugs.html#known
// 
template class RegistryItem<RgBool>;
template class RegistryItem<RgInt>;
template class RegistryItem<RgDbl>;
template class RegistryItem<RgStr>;
template class RegistryItem<RgAlg>;
template class RegistryItem<RgH1F>;
template class RegistryItem<RgH2F>;
template class RegistryItem<RgTree>;

//____________________________________________________________________________

} // genie namespace
