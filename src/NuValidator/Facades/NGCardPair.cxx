//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGCardPair

\brief    A NeuGenInputs, NeuGenConfig pair

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004          
*/
//_____________________________________________________________________________

#include "Facades/NGCardPair.h"

using std::endl;
using namespace genie::nuvld::facades;

ClassImp(NGCardPair)

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
  namespace facades {
    
    ostream & operator << (ostream & stream, const NGCardPair & pair)
    {
      pair.Print(stream);
      return stream;
    }
  }
 }
}    
//____________________________________________________________________________
NGCardPair::NGCardPair()
{
  this->Init();
}
//____________________________________________________________________________
NGCardPair::NGCardPair(const NGCardPair * pair)
{
  _inputs = pair->_inputs;
  _config = pair->_config;
}
//____________________________________________________________________________
NGCardPair::~NGCardPair()
{
  if(_inputs) delete _inputs;
  if(_config) delete _config;
}
//____________________________________________________________________________
void NGCardPair::SetInputs(const NeuGenInputs * inputs)
{
  if(_inputs) delete _inputs;

  _inputs = new NeuGenInputs(inputs);
}
//____________________________________________________________________________
void NGCardPair::SetConfig(const NeuGenConfig * config)
{
  if(_config) delete _config;

  _config = new NeuGenConfig(config);
}
//____________________________________________________________________________
void NGCardPair::Init(void) 
{
  _inputs = 0;
  _config = 0;
}
//____________________________________________________________________________
void NGCardPair::Print(ostream & stream) const
{
  stream << *_inputs << endl;
  stream << *_config << endl;
}
//____________________________________________________________________________

