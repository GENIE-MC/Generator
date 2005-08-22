//____________________________________________________________________________
/*!

\class    genie::RegistryItem

\brief    A templated concrete implementation of the RegistryItemI interface.
          Provides an arbitrary basic type (bool, int, double, string) value
          for RegistryI-type containers.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

*/
//____________________________________________________________________________

#include <string>
#include <iostream>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include "RegistryItem.h"

using std::string;
using std::endl;

using namespace genie;

//____________________________________________________________________________

template class RegistryItem<bool>;
template class RegistryItem<int>;
template class RegistryItem<double>;
template class RegistryItem<string>;
template class RegistryItem<TH1F*>;
template class RegistryItem<TH2F*>;
template class RegistryItem<TTree*>;

namespace genie {
 template
    ostream & operator << (ostream & stream, const RegistryItem<bool> &   rb);
 template
    ostream & operator << (ostream & stream, const RegistryItem<int> &    ri);
 template
    ostream & operator << (ostream & stream, const RegistryItem<double> & rd);
 template
    ostream & operator << (ostream & stream, const RegistryItem<string> & rs);
 template
    ostream & operator << (ostream & stream, const RegistryItem<TH1F*> &  rh);
 template
    ostream & operator << (ostream & stream, const RegistryItem<TH2F*> &  rh);
 template
    ostream & operator << (ostream & stream, const RegistryItem<TTree*> & rt);
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
RegistryItem<TH1F*>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
RegistryItem<TH2F*>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
RegistryItem<TTree*>::~RegistryItem()
{
  if (fItem) delete fItem;
}
//____________________________________________________________________________
template<class T> void RegistryItem<T>::Print(ostream & stream) const
{
  if(fIsLocked) stream << "[  locked] : " << fItem << endl;
  else          stream << "[unlocked] : " << fItem << endl;
}
//____________________________________________________________________________
void RegistryItem<TH1F*>::Print(ostream & stream) const
{
  if(fIsLocked) stream << "[  locked] : TH1F = "
                            << ( (fItem) ? fItem->GetName() : "NULL") << endl;
  else          stream << "[unlocked] : TH1F = "
                            << ( (fItem) ? fItem->GetName() : "NULL") << endl;
}
//____________________________________________________________________________
void RegistryItem<TH2F*>::Print(ostream & stream) const
{
  if(fIsLocked) stream << "[  locked] : TH2F = "
                            << ( (fItem) ? fItem->GetName() : "NULL") << endl;
  else          stream << "[unlocked] : TH2F = "
                            << ( (fItem) ? fItem->GetName() : "NULL") << endl;
}
//____________________________________________________________________________
void RegistryItem<TTree*>::Print(ostream & stream) const
{
  if(fIsLocked) stream << "[  locked] : TTree = "
                            << ( (fItem) ? fItem->GetName() : "NULL") << endl;
  else          stream << "[unlocked] : TTree = "
                            << ( (fItem) ? fItem->GetName() : "NULL") << endl;
}
//____________________________________________________________________________
