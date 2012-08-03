////////////////////////////////////////////////////////////////////////
/// \file  GFlavorMixerFactory.cxx
/// \brief factory for generating GENIE GFlavorMixerI class objects
///
/// \version 
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \update  2012-08-03 initial version
////////////////////////////////////////////////////////////////////////

#include "FluxDrivers/GFlavorMixerFactory.h"

#include "Messenger/Messenger.h"

namespace genie {
namespace flux {

// Define static variable which holds the one-and-only instance
GFlavorMixerFactory* GFlavorMixerFactory::fgTheInstance;

GFlavorMixerFactory::GFlavorMixerFactory() 
{
  fgTheInstance = this;   // record created self in static pointer
}

GFlavorMixerFactory::~GFlavorMixerFactory()
{
  fgTheInstance = 0;
}

GFlavorMixerFactory& GFlavorMixerFactory::Instance()
{
  // Cleaner dtor calls GFlavorMixerFactory dtor at job end
  static Cleaner cleaner;

  if ( ! fgTheInstance ) {
    // need to create one
    cleaner.UseMe();   // dummy call to quiet compiler warnings
    fgTheInstance = new GFlavorMixerFactory();
  }
  
  return *fgTheInstance;
}

GFlavorMixerI* 
GFlavorMixerFactory::GetFlavorMixer(const std::string& name)
{
  GFlavorMixerI* p = 0;
  
  // we don't want map creating an entry if it doesn't exist
  // so use map::find() not map::operator[]
  std::map<std::string, GFlavorMixerICtorFuncPtr_t>::iterator itr
    = fFunctionMap.find(name);
  if ( fFunctionMap.end() != itr ) { 
    // found an appropriate entry in the list
    GFlavorMixerICtorFuncPtr_t foo = itr->second;  // this is the function
    p = (*foo)();  // use function to create the GFlavorMixerI
  }
  if ( ! p ) {
    LOG("Flux",pWARN) << "### GFlavorMixerFactory WARNING: "
                      << "GFlavorMixerI " << name << " is not known";
  }
  return p;
}
  
bool GFlavorMixerFactory::IsKnownFlavorMixer(const std::string& name)
{
  //  check if we know the name
  bool res = false;
  std::map<std::string, GFlavorMixerICtorFuncPtr_t>::iterator itr
    = fFunctionMap.find(name);
  if ( fFunctionMap.end() != itr ) res = true;
  return res;
}

const std::vector<std::string>& 
GFlavorMixerFactory::AvailableFlavorMixers() const
{
  // list of names might be out of date due to new registrations
  // rescan the std::map on each call (which won't be frequent)
  listnames.clear();

  // scan map for registered names
  std::map<std::string, GFlavorMixerICtorFuncPtr_t>::const_iterator itr;
  for ( itr = fFunctionMap.begin(); itr != fFunctionMap.end(); ++itr )
    listnames.push_back(itr->first);

  return listnames;
}

bool GFlavorMixerFactory::RegisterCreator(std::string name, 
                                          GFlavorMixerICtorFuncPtr_t foo,
                                          bool* boolptr)
{
  // record new functions for creating processes
  fFunctionMap[name] = foo;
  fBoolPtrMap[name]  = boolptr;
  return true;
}

} // namespace flux
} // namespace genie
