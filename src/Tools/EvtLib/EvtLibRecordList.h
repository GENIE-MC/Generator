////////////////////////////////////////////////////////////////////////
// \author Chris Backhouse -- c.backhouse@ucl.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef _EVTLIBRECORDLIST_H
#define _EVTLIBRECORDLIST_H

#include <vector>
#include <string>

class TFile;
class TTree;

namespace genie{
namespace evtlib{

  //---------------------------------------------------------------------------
  struct EvtLibParticle
  {
    int pdg;
    float E, px, py, pz;
  };

  //---------------------------------------------------------------------------
  struct EvtLibRecord
  {
    EvtLibRecord();
    EvtLibRecord(float _E, int _prod_id,
                 const std::vector<EvtLibParticle>& _ps);

    /// Order by energy - this allows OnDemandRecordList to work efficiently
    bool operator<(const EvtLibRecord& rhs) const;

    float E;
    int prod_id;
    std::vector<EvtLibParticle> parts;
  };

  //---------------------------------------------------------------------------
  class IEvtLibRecordList
  {
  public:
    virtual ~IEvtLibRecordList(){}
    virtual const EvtLibRecord* GetRecord(float E) const = 0;
  };

  /// Maximum number of particles supported in a single library event record
  const int kEvtLibMaxParts = 1024;

  //---------------------------------------------------------------------------
  /// Helper for \ref SimpleRecordList and \ref OnDemandRecordList
  class RecordLoader
  {
  public:
    RecordLoader(TTree* tree);
    ~RecordLoader();

    long NRecords() const;
    EvtLibRecord GetRecord(int i) const;
  protected:
    TTree* fTree;

    float Enu;
    int prod_id;
    int nparts;
    int pdgs[kEvtLibMaxParts];
    float Es[kEvtLibMaxParts];
    float px[kEvtLibMaxParts], py[kEvtLibMaxParts], pz[kEvtLibMaxParts];
  };

  //---------------------------------------------------------------------------
  class SimpleRecordList: public IEvtLibRecordList
  {
  public:
    SimpleRecordList(TTree* tree, const std::string& prettyName);
    virtual ~SimpleRecordList(){}

    const EvtLibRecord* GetRecord(float E) const override;
  protected:
    std::vector<EvtLibRecord> fRecs;
  };

  //---------------------------------------------------------------------------
  class OnDemandRecordList: public IEvtLibRecordList
  {
  public:
    OnDemandRecordList(TTree* tree, const std::string& prettyName);
    virtual ~OnDemandRecordList(){delete fLoader;}

    const EvtLibRecord* GetRecord(float E) const override;
  protected:
    void LoadIndex() const;

    TTree* fTree;
    std::string fPrettyName;
    mutable RecordLoader* fLoader;

    mutable std::vector<std::pair<float, int>> fEnergies;

    mutable EvtLibRecord fRecord;
  };

}} // namespaces

#endif
