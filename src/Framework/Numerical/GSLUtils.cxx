//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 03, 2008 - CA
   Was first added in GENIE v2.5.1.

*/
//____________________________________________________________________________

#include <RVersion.h>

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,18,0)
#define _OLD_GSL_INTEGRATION_ENUM_TYPES_
#endif

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Numerical/GSLUtils.h"

//____________________________________________________________________________
ROOT::Math::IntegrationOneDim::Type 
     genie::utils::gsl::Integration1DimTypeFromString (string type)
{
// Returns the appropriate IntegrationOneDim type based on the input string

  string t = genie::utils::str::ToLower(type);


#ifdef _OLD_GSL_INTEGRATION_ENUM_TYPES_

  if      (t=="adaptive")          return ROOT::Math::IntegrationOneDim::ADAPTIVE;
  else if (t=="adaptive_singular") return ROOT::Math::IntegrationOneDim::ADAPTIVESINGULAR;
  else if (t=="non_adaptive")      return ROOT::Math::IntegrationOneDim::NONADAPTIVE;

  LOG("GSL", pWARN) 
       << "Unknown 1-dim GSL integration type = " << type 
       << ". Setting it to default [adaptive].";

  return ROOT::Math::IntegrationOneDim::ADAPTIVE;

#else

  if      (t=="gauss")             return ROOT::Math::IntegrationOneDim::kGAUSS;
  else if (t=="adaptive")          return ROOT::Math::IntegrationOneDim::kADAPTIVE;
  else if (t=="adaptive_singular") return ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR;
  else if (t=="non_adaptive")      return ROOT::Math::IntegrationOneDim::kNONADAPTIVE;

  LOG("GSL", pWARN) 
       << "Unknown 1-dim GSL integration type = " << type 
       << ". Setting it to default [adaptive].";

  return ROOT::Math::IntegrationOneDim::kADAPTIVE;

#endif
}
//____________________________________________________________________________
ROOT::Math::IntegrationMultiDim::Type 
     genie::utils::gsl::IntegrationNDimTypeFromString (string type)
{
// Returns the appropriate IntegrationMultiDim type based on the input string

  string t = genie::utils::str::ToLower(type);

#ifdef _OLD_GSL_INTEGRATION_ENUM_TYPES_

  if      (t=="adaptive") return ROOT::Math::IntegrationMultiDim::ADAPTIVE;
  else if (t=="plain")    return ROOT::Math::IntegrationMultiDim::PLAIN;
  else if (t=="vegas")    return ROOT::Math::IntegrationMultiDim::VEGAS;
  else if (t=="miser")    return ROOT::Math::IntegrationMultiDim::MISER;

  LOG("GSL", pWARN) 
       << "Unknown N-dim GSL integration type = " << type 
       << ". Setting it to default [adaptive].";

  return ROOT::Math::IntegrationMultiDim::ADAPTIVE;

#else

  if      (t=="adaptive") return ROOT::Math::IntegrationMultiDim::kADAPTIVE;
  else if (t=="plain")    return ROOT::Math::IntegrationMultiDim::kPLAIN;
  else if (t=="vegas")    return ROOT::Math::IntegrationMultiDim::kVEGAS;
  else if (t=="miser")    return ROOT::Math::IntegrationMultiDim::kMISER;

  LOG("GSL", pWARN) 
       << "Unknown N-dim GSL integration type = " << type 
       << ". Setting it to default [adaptive].";

  return ROOT::Math::IntegrationMultiDim::kADAPTIVE;

#endif

}
//____________________________________________________________________________

