//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenCards

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "NeuGenCards.h"

using namespace genie::nuvld;

ClassImp(NeuGenCards)

//_____________________________________________________________________________
NeuGenCards * NeuGenCards::_self = 0;
//_____________________________________________________________________________
NeuGenCards * NeuGenCards::Instance(void)
{
  if(_self == 0) _self = new NeuGenCards;

  return _self;
}
//_____________________________________________________________________________
NeuGenCards::NeuGenCards() :
TObject()
{
  _self = 0;

  _config = new NeuGenConfig;
  _inputs = new NeuGenInputs;
}
//_____________________________________________________________________________
NeuGenCards::NeuGenCards(const NeuGenCards & /*cards*/) :
TObject()
{

}
//_____________________________________________________________________________
NeuGenCards::~NeuGenCards()
{
  _self = 0;

  delete _config;
  delete _inputs;
}
//_____________________________________________________________________________
void NeuGenCards::SetInputs(const NeuGenInputs * inputs)
{
  if(_inputs) delete _inputs;

  _inputs = new NeuGenInputs(inputs);  
}
//_____________________________________________________________________________
void NeuGenCards::SetConfig(const NeuGenConfig * config)
{
  if(_config) delete _config;

  _config = new NeuGenConfig(config);  
}
//_____________________________________________________________________________
void NeuGenCards::SetDefaultConfig(void)
{
  if(_config) delete _config;

  _config = new NeuGenConfig;  
}
//_____________________________________________________________________________

