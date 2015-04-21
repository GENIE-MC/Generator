////////////////////////////////////////////////////////////////////////
/// \file  GFluxDriverFactory.h
/// \brief A class for generating concrete GFluxI derived classes
///        based on the factory pattern.  This code supplies a CPP
///        macro which allows the classes to self-register and thus
///        no modification of this class is needed in order to expand
///        the list of classes it knows about.
///
///        Implemented as a singleton holding a map between names and
///        pointers-to-functions (that call a class default constructor).
///        The functions pointers must return GFluxI*.
///
/// \version 
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
////////////////////////////////////////////////////////////////////////
#ifndef GENIE_FLUX_GFLUXDRIVERFACTORY_H
#define GENIE_FLUX_GFLUXDRIVERFACTORY_H

#include <string>
#include <vector>
#include <map>

#include "EVGDrivers/GFluxI.h"

namespace genie {
namespace flux {

// while most drivers are defined genie::flux::MySpecificDriver
// the base interface is only genie::GFluxI

// define a type for the pointer to a function that returns a 
//    genie::GFluxI* 
// i.e. calls the (typically default) ctor for the class.
typedef genie::GFluxI* (*GFluxICtorFuncPtr_t)();

class GFluxDriverFactory
{
public:
  static GFluxDriverFactory& Instance();
  // no public ctor for singleton, all user access is through Instance()

  genie::GFluxI* GetFluxDriver(const std::string&);
  // instantiate a GFluxI driver by name

  bool IsKnownFluxDriver(const std::string&);
  // check if the name is in the list of names

  const std::vector<std::string>& AvailableFluxDrivers() const;
  // return a list of available names

  bool RegisterCreator(std::string name, 
                       GFluxICtorFuncPtr_t ctorptr, bool* ptr);
  // register a new GFluxI type by passing pointer to creator function

  void PrintConfig() const;

private:
  static GFluxDriverFactory* fgTheInstance;
  // the one and only instance

  std::map<std::string, GFluxICtorFuncPtr_t> fFunctionMap;
  // mapping between known class names and a registered ctor function

  std::map<std::string, bool*> fBoolPtrMap;

  mutable std::vector<std::string> listnames;
  // copy of list of names, used solely due to AvailableFluxDrivers() 
  // method returning a const reference rather than a vector object.
  // mutable because AvailableFluxDrivers() is const, but list might 
  // need recreation if new entries have been registered.

private:
  GFluxDriverFactory();
  // private ctor, users access class via Instance()

  virtual ~GFluxDriverFactory();

  GFluxDriverFactory(const GFluxDriverFactory&);
  // method private and not implement, declared to prevent copying

  void operator=(const GFluxDriverFactory&);
  // method private and not implement, declared to prevent assignment

  // sub-class Cleaner struct is used to clean up singleton at the end of job.
  struct Cleaner {
     void UseMe() { }  // Dummy method to quiet compiler
    ~Cleaner() {
       if (GFluxDriverFactory::fgTheInstance != 0) {
         delete GFluxDriverFactory::fgTheInstance;
         GFluxDriverFactory::fgTheInstance = 0;
  } } };
  friend struct Cleaner; 

};

} // namespace flux
} // namespcae genie

// Define macro to create a function to call the class' ctor
// and then registers this function with the factory instance for later use
// Users should have in their  MyFluxClass.cc two lines that look like:
//     #include "GFluxDriverFactory.h"
//     FLUXDRIVERREG(MyFluxClass)  // no semicolon
// where "MyFluxClass" is the name of the class (assuming no special namespace)
// If the class is defined in a namespace (or two) use:
//     #include "GFluxDriverFactory.h"
//     FLUXDRIVERREG3(myspace,myAltFlux,myspace::myAltFlux) // no semicolon
//     FLUXDRIVERREG4(genie,flux,YAFlux,genie::flux::YAFlux) // no semicolon
// and either can then be retrieved from the factory using:
//     genie::flux::GFluxDriverFactory& factory =
//         genie::flux::GFluxDriverFactory::Instance();
//     genie::GFluxI* p = 0;
//     p = factory.GetFluxDriver("MyFluxClass");
//     p = factory.GetFluxDriver("myspace::myAltFlux");
//     p = factory.GetFluxDriver("genie::flux::YAFlux");
//
// The expanded code looks like:
//   genie::GFluxI* MyFluxClass_ctor_function () { return new MyFluxClass; }
//   static bool MyFluxClass_creator_registered = 
//     GFluxDriverFactory::Instance().RegisterCreator("MyFluxClass",
//                                               & MyFluxClass_ctor_function );
//   namespace myspace {
//     genie::GFluxI* myAltAltFlux_ctor_function () { return new myspace::myAltAltFlux; }
//     static bool myAltFlux_creator_registered = 
//       GFluxDriverFactory::Instance().RegisterCreator("myspace::myAltAltFlux",
//                                                 & myspace::myAltAltFlux_ctor_function ); }

#define FLUXDRIVERREG( _name ) \
  genie::GFluxI* _name ## _ctor_function () { return new _name; } \
  static bool _name ## _creator_registered =                            \
    genie::flux::GFluxDriverFactory::Instance().RegisterCreator(# _name, \
                                        & _name ## _ctor_function,        \
                                        & _name ## _creator_registered ); 

#define FLUXDRIVERREG3( _ns, _name, _fqname ) \
namespace _ns { \
  genie::GFluxI* _name ## _ctor_function () { return new _fqname; }   \
  static bool _name ## _creator_registered =                                \
    genie::flux::GFluxDriverFactory::Instance().RegisterCreator(# _fqname, \
                                        & _fqname ## _ctor_function,          \
                                        & _fqname ## _creator_registered );}

#define FLUXDRIVERREG4( _nsa, _nsb, _name, _fqname )  \
namespace _nsa { \
 namespace _nsb { \
  genie::GFluxI* _name ## _ctor_function () { return new _fqname; }   \
  static bool _name ## _creator_registered =                                \
    genie::flux::GFluxDriverFactory::Instance().RegisterCreator(# _fqname, \
                                        & _fqname ## _ctor_function,          \
                                        & _fqname ## _creator_registered );}}
#endif
