
#include <TRootIOCtor.h>

#include "Tools/ReWeight/GReWeightIOBranchDesc.h"

#include <cassert>

using namespace genie;
using namespace genie::rew;

ClassImp(GReWeightIOBranchDesc)

GReWeightIOBranchDesc::GReWeightIOBranchDesc( const std::string& name, 
                                              const double mean, 
					      const double sigpls, const double sigmin )
   : TObject()
{

   SetParameter( name, mean, sigpls, sigmin );

}

GReWeightIOBranchDesc::GReWeightIOBranchDesc( const GReWeightIOBranchDesc& brdesc )
   : TObject()
{

   fParameterName  = brdesc.fParameterName;
   fParameterMean  = brdesc.fParameterMean;
   fParameterSigmaPlus = brdesc.fParameterSigmaPlus;
   fParameterSigmaMinus = brdesc.fParameterSigmaMinus;

}

GReWeightIOBranchDesc::GReWeightIOBranchDesc( TRootIOCtor* )
   : TObject(), fParameterName(), fParameterMean(0.), fParameterSigmaPlus(0.), fParameterSigmaMinus(0.)
{
}

void GReWeightIOBranchDesc::SetParameter( const std::string& name, 
                                          const double mean, 
					  const double sigpls, const double sigmin )
{

   fParameterName = name;
   fParameterMean = mean;
   fParameterSigmaPlus = sigpls;
   fParameterSigmaMinus = sigmin;

   return;

}
