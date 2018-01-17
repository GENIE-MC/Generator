////////////////////////////////////////////////////////////////////////
/// \file  GFlavorMap.cxx
/// \brief GENIE interface for flavor modification
///
/// \version $Id: GFlavorMap.cxx,v 1.1.1.1 2010/12/22 16:18:52 p-nusoftart Exp $
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \update  2010-10-31 initial version
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>

#include "FluxDrivers/GFlavorMap.h"
#include "FluxDrivers/GFlavorMixerFactory.h"
// self register with the factory
FLAVORMIXREG4(genie,flux,GFlavorMap,genie::flux::GFlavorMap)

#include "Messenger/Messenger.h"
#define  LOG_BEGIN(a,b)   LOG(a,b)
#define  LOG_END ""

// GENIE includes
#include "Utils/StringUtils.h"

namespace genie {
namespace flux {
//____________________________________________________________________________
GFlavorMap::GFlavorMap() 
{
  // Initialize with identity matrix
  size_t jin, jout;
  for (jin=0; jin<7; ++jin) {
    for (jout=0; jout<7; ++jout ) {
      fProb[jin][jout] = ( (jin==jout) ? 1. : 0. );
    }
  }
}

GFlavorMap::~GFlavorMap() { ; }

//____________________________________________________________________________
void GFlavorMap::Config(std::string configIn)
{
  std::string config = genie::utils::str::TrimSpaces(configIn);
  LOG_BEGIN("FluxBlender", pINFO) 
    << "GFlavorMap::Config \"" << config << "\"" << LOG_END;
  
  if ( config.find("swap") == 0 || 
       config.find("map")  == 0 || 
       config.find("genie::flux::GFlavorMap") == 0 ) {
    ParseMapString(config);
  } else if ( config.find("fixedfrac") == 0 ) {
    ParseFixedfracString(config);
  } else {
    LOG_BEGIN("FluxBlender", pWARN) 
      << "GFlavorMap::Config don't know how to parse \"" 
      << config << "\"" << LOG_END;
    LOG_BEGIN("FluxBlender", pWARN) 
      << " ... will attempt \"map\" strategy" << LOG_END;
    
  }

}

//____________________________________________________________________________
void GFlavorMap::ParseMapString(std::string config)
{
  LOG_BEGIN("FluxBlender", pINFO) 
    << "GFlavorMap::ParseMapString \"" << config << "\"" << LOG_END;
  vector<string> tokens = genie::utils::str::Split(config," ");
  for (unsigned int jtok = 0; jtok < tokens.size(); ++jtok ) {
    string tok1 = tokens[jtok];
    if ( tok1 == "" ) continue;
    if ( tok1 == "swap" || tok1 == "map" ) continue;
    if ( tok1 == "genie::flux::GFlavorMap" ) continue;
    // should have the form  <int>:<int>
    vector<string> pair = genie::utils::str::Split(tok1,":");
    if ( pair.size() != 2 ) {
      LOG_BEGIN("FluxBlender", pWARN) 
        << "could not parse " << tok1 << " split size=" << pair.size() 
        << LOG_END;
      continue;
    }
    int pdg_in   = strtol(pair[0].c_str(),NULL,0);
    int indx_in  = PDG2Indx(pdg_in);
    int pdg_out  = strtol(pair[1].c_str(),NULL,0);
    int indx_out = PDG2Indx(pdg_out);
    for (int jout = 0; jout < 7; ++jout ) {
      fProb[indx_in][jout] = ( ( jout == indx_out ) ? 1 : 0 );
    }
    
  }
}

//____________________________________________________________________________
void GFlavorMap::ParseFixedfracString(std::string config)
{
  LOG_BEGIN("FluxBlender", pINFO) 
    << "GFlavorMap::ParseFixedFracString \"" << config << "\"" << LOG_END;
  vector<string> tokens = genie::utils::str::Split(config,"{}");
  for (unsigned int jtok = 0; jtok< tokens.size(); ++jtok ) {
    string tok1 = genie::utils::str::TrimSpaces(tokens[jtok]);
    if ( tok1 == "" ) continue;
    if ( tok1 == "fixedfrac" ) continue;
    // should have the form pdg:f0,f12,f14,f16,f-12,f-14,f-16
    vector<string> pair = genie::utils::str::Split(tok1,":");
    if ( pair.size() != 2 ) {
      LOG_BEGIN("FluxBlender", pWARN) 
        << "could not parse \"" << tok1 << "\" split size=" << pair.size() 
        << LOG_END;
      continue;
    }
    int pdg_in   = strtol(pair[0].c_str(),NULL,0);
    int indx_in  = PDG2Indx(pdg_in);
    vector<string> fracs = genie::utils::str::Split(pair[1],",");
    if ( fracs.size() != 7 ) {
      LOG_BEGIN("FluxBlender", pWARN) 
        << "could not parse frac list \"" << pair[1] << "\" split size=" << fracs.size() 
        << LOG_END;
      continue;
    }
    // set each value
    double psum = 0;
    for (int indx_out = 0; indx_out < 7; ++indx_out ) {
      double p = strtod(fracs[indx_out].c_str(),NULL);
      psum += p;
      fProb[indx_in][indx_out] = p;
    }
    if ( psum > 0 ) {
      // normalize to 1.0
      for (int indx_out = 0; indx_out < 7; ++indx_out )
        fProb[indx_in][indx_out] /= psum;
    }
  }
}

//____________________________________________________________________________
double GFlavorMap::Probability(int pdg_initial, int pdg_final, 
                                double /* energy */ , double /* dist */ )
{
  double prob = fProb[PDG2Indx(pdg_initial)][PDG2Indx(pdg_final)];
  if ( false ) {
    LOG_BEGIN("FluxBlender", pINFO) 
      << "Probability " << pdg_initial << "=>" << pdg_final
      << " = " << prob << LOG_END;
  }
  return prob;
}

//____________________________________________________________________________
void GFlavorMap::PrintConfig(bool /* verbose */)
{
  size_t jin, jout;
  LOG_BEGIN("FluxBlender", pINFO) 
    << "GFlavorMap::PrintConfig():" << LOG_END;

  //            1234567890[xxx]:
  std::cout << "       in      \\ out  ";
  for (jout=0; jout<7; ++jout ) 
    std::cout << "  " << std::setw(3) << Indx2PDG(jout) << "    ";
  std::cout << std::endl;

  std::cout << "----------------+";
  for (jout=0; jout<7; ++jout ) std::cout << "----------";
  std::cout << std::endl;

  for (jin=0; jin<7; ++jin) {
    std::cout << std::setw(10) << IndxName(jin) 
              << "[" << std::setw(3) << Indx2PDG(jin) << "] |  ";
    for (jout=0; jout<7; ++jout ) 
      std::cout <<  std::setw(8) << fProb[jin][jout] << " ";
    std::cout << std::endl;
  }
  
}
//____________________________________________________________________________
const char* GFlavorMap::IndxName(int indx)
{
  static const char* name[] = { "sterile", 
                                "nu_e", "nu_mu", "nu_tau",
                                "nu_e_bar", "nu_mu_bar", "nu_tau_bar" };
  return name[indx];

}

//____________________________________________________________________________
} // namespace flux
} // namespace genie
