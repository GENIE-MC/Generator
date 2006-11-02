//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 06, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <string>
#include <iostream>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include "Registry/RegistryItem.h"

using std::string;
using std::endl;

using namespace genie;

//____________________________________________________________________________

template class RegistryItem<RgBool>;
template class RegistryItem<RgInt>;
template class RegistryItem<RgDbl>;
template class RegistryItem<RgStr>;
template class RegistryItem<RgAlg>;
template class RegistryItem<RgH1F>;
template class RegistryItem<RgH2F>;
template class RegistryItem<RgTree>;

namespace genie {
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
}
//____________________________________________________________________________
namespace genie {
 template <class T>
    ostream & operator << (ostream & stream, const RegistryItem<T> & rec)
    {
      rec.Print(stream);
      return stream;
    }
}
//____________________________________________________________________________
template<class T> RegistryItem<T>::RegistryItem(T item, bool locked)
{
  fItem     = item;
  fIsLocked = locked;
}
//____________________________________________________________________________
template<class T> RegistryItem<T>::RegistryItem(const RegistryItem * ri)
{
  fItem     = ri->fItem;
  fIsLocked = ri->fIsLocked;
}
//____________________________________________________________________________
template<class T> RegistryItem<T>::~RegistryItem()
{

}
//____________________________________________________________________________
RegistryItem<RgH1F>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
RegistryItem<RgH2F>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
RegistryItem<RgTree>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
template<class T> RegistryItemI * RegistryItem<T>::Clone(void) const
{
  RegistryItemI * item = new RegistryItem<T>(fItem,fIsLocked);
  return item;
}
//____________________________________________________________________________
RgType_t RegistryItem<RgBool>::TypeInfo(void) const { return kRgBool; }
RgType_t RegistryItem<RgInt>::TypeInfo (void) const { return kRgInt;  }
RgType_t RegistryItem<RgDbl>::TypeInfo (void) const { return kRgDbl;  }
RgType_t RegistryItem<RgStr>::TypeInfo (void) const { return kRgStr;  }
RgType_t RegistryItem<RgAlg>::TypeInfo (void) const { return kRgAlg;  }
RgType_t RegistryItem<RgH1F>::TypeInfo (void) const { return kRgH1F;  }
RgType_t RegistryItem<RgH2F>::TypeInfo (void) const { return kRgH2F;  }
RgType_t RegistryItem<RgTree>::TypeInfo(void) const { return kRgTree; }
//____________________________________________________________________________
template<class T> void RegistryItem<T>::Print(ostream & stream) const
{
  if(fIsLocked) stream << "[  locked] : " << fItem;
  else          stream << "[unlocked] : " << fItem;
}
//____________________________________________________________________________
void RegistryItem<RgAlg>::Print(ostream & stream) const
{
  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") << " : " << (fItem);
}
//____________________________________________________________________________
void RegistryItem<RgH1F>::Print(ostream & stream) const
{
  TH1F * histo = dynamic_cast<TH1F *>(fItem);
  if(!histo) stream << "*** NULL RgH1F ***";

  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") 
         << " : " << fItem->GetName();
}
//____________________________________________________________________________
void RegistryItem<RgH2F>::Print(ostream & stream) const
{
  TH2F * histo = dynamic_cast<TH2F *>(fItem);
  if(!histo) stream << "*** NULL RgH2F ***";

  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") 
          << " : " << fItem->GetName();
}
//____________________________________________________________________________
void RegistryItem<RgTree>::Print(ostream & stream) const
{
  TTree * tree = dynamic_cast<TTree *>(fItem);
  if(!tree) stream << "*** NULL RgTree ***";

  stream << ((fIsLocked) ? "[  locked]" : "[unlocked]") 
         << " : " << fItem->GetName();
}
//____________________________________________________________________________
