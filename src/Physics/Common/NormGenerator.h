//____________________________________________________________________________
/*!

\class    genie::NormGenerator

\brief    Normalization channel. Its main property is a constant cross section 
          per nucleon over the whole energy range. For nu/charged probes this produces 
          NC/EM events with the probe & target "echoed" back as final state particles.

\ref      [1] GENIE docdb 297


\author   Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research

\created  May 16, 2022

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org      
    
*/
//____________________________________________________________________________

#ifndef _NORM_GENERATOR_H_
#define _NORM_GENERATOR_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/EventGen/EventRecordVisitorI.h"


namespace genie {

  class NormGenerator: public EventRecordVisitorI {

  public :
    NormGenerator();
    NormGenerator(string config);
    ~NormGenerator();

    // implement the EventRecordVisitorI interface
    void ProcessEventRecord(GHepRecord * event_rec) const;

    // overload the Algorithm::Configure() methods to load private data
    // members from configuration options
    void Configure(const Registry & config);
    void Configure(string config);

    // methods to load sub-algorithms and config data from the Registry
    void LoadConfig (void);



  private:
    mutable const XSecAlgorithmI * fXSecModel;
  };

}      // genie namespace
#endif // _NORM_GENERATOR_H_
