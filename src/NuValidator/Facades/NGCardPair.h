//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGCardPair

\brief    A NeuGenInputs, NeuGenConfig pair

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _CARD_PAIR_H_
#define _CARD_PAIR_H_

#include <iostream>

#include <TObject.h>

#include "Facades/NeuGenInputs.h"
#include "Facades/NeuGenConfig.h"

using std::ostream;

namespace genie   {
namespace nuvld   {
namespace facades {

class NGCardPair : public TObject
{
public:

  NGCardPair();
  NGCardPair(const NGCardPair * pair);
  ~NGCardPair();

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NGCardPair & pair);

  NeuGenInputs * GetInputs (void) const  { return _inputs; }
  NeuGenConfig * GetConfig (void) const  { return _config; }

  void SetInputs (const NeuGenInputs * inputs);
  void SetConfig (const NeuGenConfig * config);

private:

  void Init(void);
  
  NeuGenConfig * _config;
  NeuGenInputs * _inputs;

  ClassDef(NGCardPair, 1)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif // _CARD_PAIR_H_
