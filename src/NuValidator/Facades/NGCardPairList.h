//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGCardPairList

\brief    A list of NGCardPair objects

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _CARD_PAIR_LIST_H_
#define _CARD_PAIR_LIST_H_

#include <map>
#include <vector>
#include <string>

#include <TObject.h>

#include "Facades/NGCardPair.h"

using std::map;
using std::vector;
using std::string;

namespace genie   {
namespace nuvld   {
namespace facades {

class NGCardPairList : public TObject
{
public:

  friend ostream & operator << (ostream & stream, const NGCardPairList & pairs);

  NGCardPairList();
  NGCardPairList(const NGCardPairList * pairs);
  ~NGCardPairList();

  void                   AddCardPair    (string name, NGCardPair * pair);
  void                   Merge          (const NGCardPairList * pairs);
  void                   Erase          (string name);
  const vector<string> * GetListOfNames (void)        const;
  NGCardPair *           GetCardPair    (string name) const;
  unsigned int           NPairs         (void)        const;
  bool                   Exists         (string name) const;

private:

  map<string, NGCardPair *> _card_pairs_map;

ClassDef(NGCardPairList, 1)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif
