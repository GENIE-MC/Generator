////////////////////////////////////////////////////////////////////////
/// \file  GFlavorMixerFactory.h
/// \brief A class for generating concrete GFlavorMixerI derived classes
///        based on the factory pattern.  This code supplies a CPP
///        macro which allows the classes to self-register and thus
///        no modification of this class is needed in order to expand
///        the list of classes it knows about.
///
///        Implemented as a singleton holding a map between names and
///        pointers-to-functions (that call a class default constructor).
///        The functions pointers must return GFlavorMixerI*.
///
/// \version 
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
////////////////////////////////////////////////////////////////////////
#ifndef GENIE_FLUX_GFLAVORMIXERFACTORY_H
#define GENIE_FLUX_GFLAVORMIXERFACTORY_H

#include <string>
#include <vector>
#include <map>

#include "FluxDrivers/GFlavorMixerI.h"

namespace genie {
namespace flux {

// define a type for the pointer to a function that returns a 
//    genie::flux::GFlavorMixerI* 
// i.e. calls the (typically default) ctor for the class.
typedef genie::flux::GFlavorMixerI* (*GFlavorMixerICtorFuncPtr_t)();

class GFlavorMixerFactory
{
public:
  static GFlavorMixerFactory& Instance();
  // no public ctor for singleton, all user access is through Instance()

  genie::flux::GFlavorMixerI* GetFlavorMixer(const std::string&);
  // instantiate a PhysProc by name

  bool IsKnownFlavorMixer(const std::string&);
  // check if the name is in the list of names

  const std::vector<std::string>& AvailableFlavorMixers() const;
  // return a list of available names

  bool RegisterCreator(std::string name, 
                       GFlavorMixerICtorFuncPtr_t ctorptr, bool* ptr);
  // register a new GFlavorMixerI type by passing pointer to creator function

private:
  static GFlavorMixerFactory* fgTheInstance;
  // the one and only instance

  std::map<std::string, GFlavorMixerICtorFuncPtr_t> fFunctionMap;
  // mapping between known class names and a registered ctor function

  std::map<std::string, bool*> fBoolPtrMap;

  mutable std::vector<std::string> listnames;
  // copy of list of names, used solely due to AvailableFlavorMixers() 
  // method returning a const reference rather than a vector object.
  // mutable because AvailableFlavorMixers() is const, but list might 
  // need recreation if new entries have been registered.

private:
  GFlavorMixerFactory();
  // private ctor, users access class via Instance()

  virtual ~GFlavorMixerFactory();

  GFlavorMixerFactory(const GFlavorMixerFactory&);
  // method private and not implement, declared to prevent copying

  void operator=(const GFlavorMixerFactory&);
  // method private and not implement, declared to prevent assignment

  // sub-class Cleaner struct is used to clean up singleton at the end of job.
  struct Cleaner {
     void UseMe() { }  // Dummy method to quiet compiler
    ~Cleaner() {
       if (GFlavorMixerFactory::fgTheInstance != 0) {
         delete GFlavorMixerFactory::fgTheInstance;
         GFlavorMixerFactory::fgTheInstance = 0;
  } } };
  friend struct Cleaner; 

};

} // namespace flux
} // namespcae genie

// Define macro to create a function to call the class' ctor
// and then registers this function with the factory instance for later use
// Users should have in their  myPhyList.cc two lines that look like:
//     #include "GFlavorMixerFactory.h"
//     FLAVORMIXREG(MyMixerClass)  // no semicolon
// where "MyMixerClass" is the name of the class (assuming no special namespace)
// If the class is defined in a namespace (or two) use:
//     #include "GFlavorMixerFactory.h"
//     FLAVORMIXREG3(myspace,myAltMixer,myspace::myAltMixer) // no semicolon
//     FLAVORMIXREG4(genie,flux,YAMixer,genie::flux::YAMixer) // no semicolon
// and either can then be retrieved from the factory using:
//     GFlavorMixerFactory& factory =
//         GFlavorMixerFactory::Instance();
//     genie::flux::GFlavorMixerI* p = 0;
//     p = factory.GetFlavorMixer("MyMixerClass");
//     p = factory.GetFlavorMixer("myspace::myAltMixer");
//     p = factory.GetFlavorMixer("genie::flux::YAMixer");
//
// The expanded code looks like:
//   genie::flux::GFlavorMixerI* MyMixerClass_ctor_function () { return new MyMixerClass; }
//   static bool MyMixerClass_creator_registered = 
//     GFlavorMixerFactory::Instance().RegisterCreator("MyMixerClass",
//                                               & MyMixerClass_ctor_function );
//   namespace myspace {
//     genie::flux::GFlavorMixerI* myAltAltMixer_ctor_function () { return new myspace::myAltAltMixer; }
//     static bool myAltMixer_creator_registered = 
//       GFlavorMixerFactory::Instance().RegisterCreator("myspace::myAltAltMixer",
//                                                 & myspace::myAltAltMixer_ctor_function ); }

#define FLAVORMIXREG( _name ) \
  genie::flux::GFlavorMixerI* _name ## _ctor_function () { return new _name; } \
  static bool _name ## _creator_registered =                            \
    genie::flux::GFlavorMixerFactory::Instance().RegisterCreator(# _name, \
                                        & _name ## _ctor_function,        \
                                        & _name ## _creator_registered ); 

#define FLAVORMIXREG3( _ns, _name, _fqname ) \
namespace _ns { \
  genie::flux::GFlavorMixerI* _name ## _ctor_function () { return new _fqname; }   \
  static bool _name ## _creator_registered =                                \
    genie::flux::GFlavorMixerFactory::Instance().RegisterCreator(# _fqname, \
                                        & _fqname ## _ctor_function,          \
                                        & _fqname ## _creator_registered );}

#define FLAVORMIXREG4( _nsa, _nsb, _name, _fqname )  \
namespace _nsa { \
 namespace _nsb { \
  genie::flux::GFlavorMixerI* _name ## _ctor_function () { return new _fqname; }   \
  static bool _name ## _creator_registered =                                \
    genie::flux::GFlavorMixerFactory::Instance().RegisterCreator(# _fqname, \
                                        & _fqname ## _ctor_function,          \
                                        & _fqname ## _creator_registered );}}
#endif
