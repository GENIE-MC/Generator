//____________________________________________________________________________
/*!

\class    genie::Algorithm

\brief    Algorithm abstract base class. 

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 02, 2004
 
*/
//____________________________________________________________________________

#include "Algorithm/Algorithm.h"
#include "Config/ConfigPool.h"
#include "Messenger/Messenger.h"

using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const Algorithm & alg)
  {
     // print algorithm name & parameter-set
     stream << alg.fName << "/" << alg.fParamSet << endl;

     // print algorithm configuration
     stream << *(alg.fConfig);

     return stream;
  }
}
//____________________________________________________________________________
Algorithm::Algorithm()
{
  fConfigIsOwned = false; 
  fConfig        = 0;
} 
//____________________________________________________________________________
Algorithm::Algorithm(const char * param_set) :
fParamSet(param_set)
{ 
  fConfigIsOwned = false; 
  fConfig        = 0;
}
//____________________________________________________________________________
Algorithm::~Algorithm()
{
  //if(fConfigIsOwned && fConfig) delete fConfig; -- check
} 
//____________________________________________________________________________
void Algorithm::Configure(const Registry & config)
{
// Configure the Algorithm using the input configuration Registry

  // check if its already owns a configuration Registry & delete it

  if(fConfigIsOwned && fConfig) {

    LOG("Algorithm", pINFO) << "Deleting old owned configuration";
        
    //delete fConfig; -- check
  }

  // create a new configuration Registry & raise the "owned" flag

  fConfig        = new Registry(config);
  fConfigIsOwned = true; 

  LOG("Algorithm", pINFO) << *fConfig;
}
//____________________________________________________________________________
void Algorithm::Configure(string param_set)
{
// Configure the Algorithm looking up at the ConfigPool singleton for a
// configuration Registry corresponding to the input named parameter set.

  fParamSet = param_set;

  FindConfig();
}
//____________________________________________________________________________
void Algorithm::FindConfig(void)
{
// Finds its configration Registry from the ConfigPool and gets a pointer to
// it. If the Registry comes from the ConfigPool then the Algorithm does not
// own its configuration (the ConfigPool singleton keeps the ownership and the 
// responsibility to -eventually- delete all the Registries it instantiates
// by parsing the XML config files).

  ConfigPool * pool = ConfigPool::Instance();

  Registry * config = pool->FindRegistry( this );

  if(!config)

        // notify & keep whatever config Registry was used before.

        LOG("Algorithm", pWARN)         
               << "No Configuration available for " 
                      << this->Name() << "/" << this->ParamSet() 
                                                    << " at the ConfigPool";
  else {
        // check if its already owns a configuration Registry & delete it;

        if(fConfigIsOwned && fConfig) delete fConfig;

        // point to the configuration Registry retrieved from the ConfigPool
        // and raise the "not-owned" flag.

        fConfig        = config;
        fConfigIsOwned = false;

        LOG("Algorithm", pDEBUG) << ENDL << *fConfig;
  }
}
//____________________________________________________________________________
