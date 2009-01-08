//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGCardPairList

\brief    A list of NGCardPair objects

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#include <iostream>
#include <TMath.h>

#include "Facades/NGCardPairList.h"

using std::endl;
using std::cout;

using namespace genie::nuvld::facades;

ClassImp(NGCardPairList)

//____________________________________________________________________________
namespace genie   {
 namespace nuvld   {
  namespace facades {
    
    ostream & operator << (ostream & stream, const NGCardPairList & pairs)
    {
      map<string, NGCardPair *>::const_iterator pair_iter;

      stream << "Contained Card Pairs: " << endl;

      for(pair_iter = pairs._card_pairs_map.begin();
                      pair_iter != pairs._card_pairs_map.end(); ++pair_iter) {

         stream << "NGCardPair: " << pair_iter->first << ": " << endl;
         stream << *(pair_iter->second);
      }
      return stream;
    }
  }  
 }  
}
//______________________________________________________________________________
NGCardPairList::NGCardPairList()
{

}
//______________________________________________________________________________
NGCardPairList::NGCardPairList(const NGCardPairList * pairs)
{
  map<string, NGCardPair *>::const_iterator pair_iter;

  for(pair_iter = pairs->_card_pairs_map.begin();
                       pair_iter != pairs->_card_pairs_map.end(); ++pair_iter) {

     string        name  =  pair_iter->first;
     NGCardPair *  cards =  pair_iter->second;

     this->AddCardPair(name, cards);
  }
}
//______________________________________________________________________________
NGCardPairList::~NGCardPairList()
{

}
//______________________________________________________________________________
void NGCardPairList::AddCardPair(string name, NGCardPair * pair)
{
  NGCardPair * cloned_pair = new NGCardPair(pair);

  _card_pairs_map.insert(
                      map<string, NGCardPair *>::value_type(name, cloned_pair) );
}
//______________________________________________________________________________
void NGCardPairList::Merge(const NGCardPairList * pairs)
{
  map<string, NGCardPair *>::const_iterator pair_iter;

  for(pair_iter = pairs->_card_pairs_map.begin();
                       pair_iter != pairs->_card_pairs_map.end(); ++pair_iter) {

     string       name  =  pair_iter->first;
     NGCardPair * cards =  pair_iter->second;

     this->AddCardPair(name, cards);
  }
}
//______________________________________________________________________________
void NGCardPairList::Erase(string name)
{
   _card_pairs_map.erase(name);
}
//______________________________________________________________________________
const vector<string> * NGCardPairList::GetListOfNames(void) const
{
   map<string, NGCardPair *>::const_iterator pair_iter;

   vector<string> * names = new vector<string>;

   for(pair_iter = _card_pairs_map.begin();
                    pair_iter != _card_pairs_map.end(); ++pair_iter)
                                           names->push_back( pair_iter->first );
   return names;
}
//______________________________________________________________________________
NGCardPair * NGCardPairList::GetCardPair(string name) const
{
  if( _card_pairs_map.count(name) == 1 ) {

   map<string, NGCardPair *>::const_iterator pair_iter = _card_pairs_map.find(name);
   return pair_iter->second;

  } else return 0;
}
//______________________________________________________________________________
unsigned int NGCardPairList::NPairs(void) const
{
  return _card_pairs_map.size();
}
//______________________________________________________________________________
bool NGCardPairList::Exists(string name) const
{
  return (_card_pairs_map.count(name) == 1);
}
//______________________________________________________________________________


