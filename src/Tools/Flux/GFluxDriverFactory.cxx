////////////////////////////////////////////////////////////////////////
/// \file  GFluxDriverFactory.cxx
/// \brief factory for generating GENIE GFluxI class objects
///
/// \version 
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \update  2012-08-03 initial version
////////////////////////////////////////////////////////////////////////

#include "FluxDrivers/GFluxDriverFactory.h"

#include "Messenger/Messenger.h"
#include <iostream>
#include <iomanip>

namespace genie {
namespace flux {

// Define static variable which holds the one-and-only instance
GFluxDriverFactory* GFluxDriverFactory::fgTheInstance;

GFluxDriverFactory::GFluxDriverFactory() 
{
  fgTheInstance = this;   // record created self in static pointer
}

GFluxDriverFactory::~GFluxDriverFactory()
{
  fgTheInstance = 0;
}

GFluxDriverFactory& GFluxDriverFactory::Instance()
{
  // Cleaner dtor calls GFluxDriverFactory dtor at job end
  static Cleaner cleaner;

  if ( ! fgTheInstance ) {
    // need to create one
    cleaner.UseMe();   // dummy call to quiet compiler warnings
    fgTheInstance = new GFluxDriverFactory();
  }
  
  return *fgTheInstance;
}

GFluxI* 
GFluxDriverFactory::GetFluxDriver(const std::string& name)
{
  GFluxI* p = 0;
  
  // we don't want map creating an entry if it doesn't exist
  // so use map::find() not map::operator[]
  std::map<std::string, GFluxICtorFuncPtr_t>::iterator itr
    = fFunctionMap.find(name);
  if ( fFunctionMap.end() != itr ) { 
    // found an appropriate entry in the list
    GFluxICtorFuncPtr_t foo = itr->second;  // this is the function
    p = (*foo)();  // use function to create the GFluxI
  }
  if ( ! p ) {
    LOG("Flux",pWARN) << "### GFluxDriverFactory WARNING: "
                      << "GFluxI " << name << " is not known";
  }
  return p;
}
  
bool GFluxDriverFactory::IsKnownFluxDriver(const std::string& name)
{
  //  check if we know the name
  bool res = false;
  std::map<std::string, GFluxICtorFuncPtr_t>::iterator itr
    = fFunctionMap.find(name);
  if ( fFunctionMap.end() != itr ) res = true;
  return res;
}

const std::vector<std::string>& 
GFluxDriverFactory::AvailableFluxDrivers() const
{
  // list of names might be out of date due to new registrations
  // rescan the std::map on each call (which won't be frequent)
  listnames.clear();

  // scan map for registered names
  std::map<std::string, GFluxICtorFuncPtr_t>::const_iterator itr;
  for ( itr = fFunctionMap.begin(); itr != fFunctionMap.end(); ++itr )
    listnames.push_back(itr->first);

  return listnames;
}

bool GFluxDriverFactory::RegisterCreator(std::string name, 
                                          GFluxICtorFuncPtr_t foo,
                                          bool* boolptr)
{
  // record new functions for creating processes
  fFunctionMap[name] = foo;
  fBoolPtrMap[name]  = boolptr;
  return true;
}

void GFluxDriverFactory::PrintConfig() const
{
  const std::vector<std::string>& avail = AvailableFluxDrivers();
  size_t n = avail.size();
  std::cout << "GFluxDriverFactory has the following drivers registered:"
            << std::endl;
  for (size_t i=0; i<n; ++i) {
    std::cout << "  [" << std::setw(3) << i << "]  " 
              << avail[i] << std::endl;
  }
}

} // namespace flux
} // namespace genie
