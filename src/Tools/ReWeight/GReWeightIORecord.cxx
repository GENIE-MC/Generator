
#include <TRootIOCtor.h>

#include "Tools/ReWeight/GReWeightIORecord.h"

#include <cassert>

using namespace genie;
using namespace genie::rew;

ClassImp(GReWeightIORecord)

GReWeightIORecord::GReWeightIORecord()
   : TObject(),
     fOrigEvtNum(-1)
{

   fRWResults.clear();

}

GReWeightIORecord::GReWeightIORecord( const GReWeightIORecord& rwrec )
   : TObject()
{

   this->Copy(rwrec);

}

GReWeightIORecord::GReWeightIORecord( TRootIOCtor* )
   : TObject(),
     fOrigEvtNum(-1)
{

   fRWResults.clear();

}

void GReWeightIORecord::Reset()
{

   fOrigEvtNum=-1;
   fRWResults.clear();

   return;

}


void GReWeightIORecord::Copy( const GReWeightIORecord& rwrec )
{

   this->Reset();
   fOrigEvtNum    = rwrec.fOrigEvtNum; 
   fRWResults     = rwrec.fRWResults;
   
   return;

}

void GReWeightIORecord::Insert( const double twk, const double wt )
{
   
   fRWResults.push_back( GReWeightInfo(twk,wt) );

   return;

}
