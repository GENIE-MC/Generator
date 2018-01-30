#ifndef GRWIORECORD_H_
#define GRWIORECORD_H_

#include <vector>

#include <TObject.h>

class TRootIOCtor;

namespace genie {
namespace rew   {

class GReWeightInfo : public TObject {

   public:
      GReWeightInfo() : TObject(), fTweak(0.), fWeight(1.) {}
      GReWeightInfo( double twk, double wt ) : TObject() { fTweak=twk; fWeight=wt; }
      GReWeightInfo( TRootIOCtor* ) : TObject(), fTweak(0.), fWeight(1.) {}
      GReWeightInfo( const GReWeightInfo& r ) : TObject() { fTweak=r.fTweak; fWeight=r.fWeight; }
      ~GReWeightInfo() {}
      
      double GetTweak()  const { return fTweak; }
      double GetWeight() const { return fWeight; }

   private:
      
      double fTweak;
      double fWeight;
      
// NOTE (JY): the 2nd parameter does have a "sense" of versioning, but it's more than that;
//            in particular, 2 is the "minimal level" for reasonable I/O (e.g.0 would mean 
//            "never store it");
//            one would increase it when chaning/adding to the data member list, etc.
//
ClassDef(GReWeightInfo,2)

};

} // end namespace rew
} // end namespace genie

namespace genie
{
namespace rew
{

class GReWeightIORecord : public TObject {

   public:
      using TObject::Copy; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
      
      GReWeightIORecord();
      GReWeightIORecord( const GReWeightIORecord& );
      GReWeightIORecord( TRootIOCtor* );
      
      ~GReWeightIORecord() {}
      
      void Reset();
      void Copy( const GReWeightIORecord& );
             
      void SetOriginalEvtNumber( const int ievt ) { fOrigEvtNum = ievt; return; }
      void Insert( const double, const double );
      
      int    GetOriginalEvtNumber()   const { return fOrigEvtNum; }
      int    GetNumOfRWResults()      const { return fRWResults.size(); }
      double GetTweak( const int i )  const { return ( ( fRWResults.size() > 0 && i >= 0 && i < (int)fRWResults.size() ) ? fRWResults[i].GetTweak() : 0. ); }
      double GetWeight( const int i ) const { return ( ( fRWResults.size() > 0 && i >= 0 && i < (int)fRWResults.size() ) ? fRWResults[i].GetWeight() : 1. ); }

   private:
                  
      int                 fOrigEvtNum;
      std::vector<GReWeightInfo> fRWResults;
      
ClassDef(GReWeightIORecord,2)

};

} // end namespace rew
} // end namespace genie

#endif // RWOBJECT_H_
