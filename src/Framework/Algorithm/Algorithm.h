//____________________________________________________________________________
/*!

\class    genie::Algorithm

\brief    Algorithm abstract base class.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Marco Roda <mroda \at liverpool.ac.uk>
          University of Liverpool

\created  May 02, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ALGORITHM_H_
#define _ALGORITHM_H_

#include <string>
#include <sstream>
#include <iostream>
#include <cassert>
#include <map>

#include "Framework/Algorithm/AlgStatus.h"
#include "Framework/Algorithm/AlgCmp.h"
#include "Framework/Algorithm/AlgId.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Framework/Messenger/Messenger.h"

using std::string;
using std::ostream;
using std::map;

namespace genie {

class Algorithm;

ostream & operator << (ostream & stream, const Algorithm & alg);

typedef map <string, Algorithm *>                 AlgMap;
typedef map <string, Algorithm *>::iterator       AlgMapIter;
typedef map <string, Algorithm *>::const_iterator AlgMapConstIter;
typedef pair<string, Algorithm *>                 AlgMapPair;

class Algorithm {

public:
  virtual ~Algorithm();

  //! Configure the algorithm with an external registry
  //!   The registry is merged with the top level registry if it is owned,
  //!   Otherwise a copy of it is added with the highest priority
  virtual void Configure (const Registry & config);

  //! Configure the algorithm from the AlgoConfigPool
  //!   based on param_set string given in input
  //! An algorithm contains a vector of registries coming from different 
  //!   xml configuration files, which are loaded according a very precise prioriy
  //! This methods will load a number registries in order of priority:
  //!   1) "Tunable" parameter set from CommonParametes. This is loaded with the
  //!      highest prioriry and it is designed to be used for tuning procedure
  //!      Usage not expected from the user.
  //!   2) For every string defined in "CommonParame" the corresponding parameter set will be loaded 
  //!      from CommonParameter.xml
  //!   3) parameter set specified by the config string and defined in the xml file of the algorithm
  //!   4) if config is not "Default" also the Default parameter set from the same xml file will be loaded
  //!      Effectively this avoids the repetion of a parameter when it is not changed in the requested configuration
  virtual void Configure (string config);

  //! Lookup configuration from the config pool
  //!   Similar logic from void Configure(string)
  virtual void FindConfig (void);

  //! Get configuration registry
  //!  Evaluate the summary of the configuration and returns it
  //!  The summary of a configuration is a merge of all the registries 
  //!  known to the algorithm (see Configure() methods) but every parameter is appearing only
  //!  once and in case of repetitions, only the parameter from the registry with the highest prioriry 
  //!  is considered.
  virtual const Registry & GetConfig(void) const ;

  //! Returns the pointer of the summary registry, see previous method
  //!  Gives access to the summary so it could be changed. 
  //!  The usage of this method is deprecated as it is mantained only for back compatibility. 
  //! If you need to add or chage a parter (or more), use the AddTopRegistry() instead
  Registry * GetOwnedConfig(void);

  //! Get algorithm ID
  virtual const AlgId & Id(void) const { return fID; }

  //! Get algorithm status
  virtual AlgStatus_t GetStatus(void) const { return fStatus; }

  //! Allow reconfigration after initializaton?
  //! Algorithms may opt-out, if reconfiguration is not necessary,
  //! to improve event reweighting speed.
  virtual bool AllowReconfig(void) const { return fAllowReconfig; }

  //! Compare with input algorithm
  virtual AlgCmp_t Compare(const Algorithm * alg) const;

  //! Set algorithm ID
  virtual void SetId(const AlgId & id);
  virtual void SetId(string name,  string config);

  //! Access the sub-algorithm pointed to by the input key, either from the
  //! local pool or from AlgFactory's pool
  const Algorithm * SubAlg(const RgKey & registry_key) const;

  //! Clone the configuration registry looked up from the configuration pool
  //! and take its ownership
  void AdoptConfig (void);

  //! Take ownership of the algorithms subtructure (sub-algorithms,...)
  //! by copying them from the AlgFactory pool to the local pool
  //! Also bring all the configuration variables to the top level config Registry.
  //! This can be used to group together a series of algorithms & their
  //! configurations and extract (a clone of) this group from the shared
  //! pools. Having a series of algorithms/configurations behaving as a
  //! monolithic block, with a single point of configuration (the top level)
  //! is to be used when bits & pieces of GENIE are used in isolation for
  //! data fitting or reweighting
  void AdoptSubstructure (void);

  //! Print algorithm info
  virtual void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const Algorithm & alg);

protected:
  Algorithm();
  Algorithm(string name);
  Algorithm(string name, string config);

  void Initialize         (void);
  void DeleteConfig       (void);
  void DeleteSubstructure (void);

  //! Split an incoming configuration Registry into a block valid for this algorithm
  //! Ownership of the returned registry belongs to the algo
  Registry * ExtractLocalConfig( const Registry & in ) const ;
  //! Split an incoming configuration Registry into a block valid for the sub-algo identified by alg_key
  Registry * ExtractLowerConfig( const Registry & in, const string & alg_key ) const ;

  bool         fAllowReconfig; ///<
  //  bool         fOwnsConfig;    ///< true if it owns its config. registry
  bool         fOwnsSubstruc;  ///< true if it owns its substructure (sub-algs,...)
  AlgId        fID;            ///< algorithm name and configuration set


  /// ideally these members should go private
  /// Registry will be access only through the GetParam method
  vector<Registry*>  fConfVect ;   ///< configurations registries from various sources
                                   ///<  the order of the vector is the precedence in case of repeated parameters
                                   ///<  position 0 -> Highest precedence
  vector<bool>       fOwnerships ; ///< ownership for every registry in fConfVect

  AlgStatus_t  fStatus;        ///< algorithm execution status
  AlgMap *     fOwnedSubAlgMp; ///< local pool for owned sub-algs (taken out of the factory pool)

  //! Ideal access to a parameter value from the vector of registries
  //! Returns true if the value is found and the parameters is set
  template<class T>
    bool GetParam( const RgKey & name, T & p, bool is_top_call = true ) const ;

  //! Ideal access to a parameter value from the vector of registries,
  //! With default value. Returns true if the value is set from the
  //! registries, false if the value is the default
  template<class T>
     bool GetParamDef( const RgKey & name, T & p, const T & def ) const ;

  //! Handle to load vectors of parameters
  //! It looks for different registry item with name comm_name0, comm_name1, etc...
  template<class T>
    bool GetParamVect( const std::string & comm_name, std::vector<T> & v,
    		           unsigned int max, bool is_top_call = true ) const ;

  int   AddTopRegistry( Registry * rp, bool owns = true );  ///< add registry with top priority, also update ownership
  int   AddLowRegistry( Registry * rp, bool owns = true );  ///< add registry with lowest priority, also update ownership
  int   MergeTopRegistry( const Registry & r ) ;            ///< Merge with top level registry if first reg of the vector is owned
                                                            ///< Otherwise an owned copy is added as a top registry
  int   AddTopRegisties( const vector<Registry*> & rs, bool owns = false ) ; ///< Add registries with top priority, also udated Ownerships  

private:

  Registry *   fConfig;        ///< Summary configuration derived from fConvVect, not necessarily allocated

};

}       // genie namespace

#ifndef __CINT__  // don't even try for ROOT 5
#include "Framework/Algorithm/Algorithm.icc"
#endif

#endif  // _ALGORITHM_H_
