////////////////////////////////////////////////////////////////////////
// \author Christopher Backhouse -- c.backhouse@ucl.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef _RECORDLIST_H
#define _RECORDLIST_H

#include <vector>

class TFile;
class TTree;

namespace genie{
namespace vmc{
  //---------------------------------------------------------------------------
  struct Particle
  {
    int pdg;
    float E, px, py, pz;
  };

  //---------------------------------------------------------------------------
  struct Record
  {
    Record();
    Record(float _E, /*float _w,*/ int _prod_id, const std::vector<Particle>& _ps);

    /// Order by energy
    bool operator<(const Record& rhs) const;

    float E;
    //    float weight;
    int prod_id;
    std::vector<Particle> parts;
  };

  //---------------------------------------------------------------------------
  class IRecordList
  {
  public:
    virtual ~IRecordList(){}
    virtual const Record* GetRecord(float E) const = 0;
  };

  //---------------------------------------------------------------------------
  /// Helper for \ref SimpleRecordList and \ref OnDemandRecordList
  class RecordLoader
  {
  public:
    RecordLoader(const char* fname);
    ~RecordLoader();

    long NRecords() const;
    Record GetRecord(int i) const;
  protected:
    TFile* fFile;
    TTree* fTree;

    float Enu;//, weight;
    int prod_id;
    int nparts;
    int pdgs[1024];
    float Es[1024], px[1024], py[1024], pz[1024];
  };

  //---------------------------------------------------------------------------
  class SimpleRecordList: public IRecordList
  {
  public:
    SimpleRecordList(const char* fname);
    virtual ~SimpleRecordList(){}

    const Record* GetRecord(float E) const override;
  protected:
    std::vector<Record> fRecs;
  };

  //---------------------------------------------------------------------------
  class OnDemandRecordList: public IRecordList
  {
  public:
    OnDemandRecordList(const char* fname);
    virtual ~OnDemandRecordList(){}

    const Record* GetRecord(float E) const override;
  protected:
    RecordLoader fLoader;

    std::vector<std::pair<float, int>> fEnergies;

    mutable Record fRecord;
  };

}} // namespaces

#endif
