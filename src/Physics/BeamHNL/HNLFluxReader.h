//----------------------------------------------------------------------------
/*!

  This reads in flux histograms and constructs SimpleHNL objects
  with appropriate energy + momentum + vertex position + pol

  The couplings to SM are not relevant at this stage
  However, labelling the flux with its couplings should be so we don't lose track

  Solution: Read in couplings directly from the *file name*. 

\namespace  genie::HNL::HNLFluxReader

\brief      Reads in user inout flux from histograms

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    January 5th, 2022

\cpright    ??? - TBD

 */
//----------------------------------------------------------------------------
// TODO: Calculate polarisation vector magnitude in 3-body HNL production
//----------------------------------------------------------------------------

#ifndef _HNL_FLUXREADER_H_
#define _HNL_FLUXREADER_H_

// -- ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"

// -- GENIE includes

#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
//#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/BeamHNL/HNLEnums.h"
#include "Physics/BeamHNL/HNLKinUtils.h"
#include "Physics/BeamHNL/SimpleHNL.h"

namespace genie{
namespace HNL{

  class SimpleHNL;

  namespace HNLFluxReader{

    extern std::string fPath; // flux file to poke

    /// perform selections from input files
    int selectMass( const double mN );
    void selectFile( std::string fin, const double mN ); // find but don't open the file

    /// get the histogram and energy from it
    TH1F * getFluxHist1F( std::string fin, std::string hName, int HType );
    // and overloaded method for the newer interface
    TH1F * getFluxHist1F( std::string fin, int masspoint, bool isParticle );

    TH3D * getFluxHist3D( std::string fin, std::string dirName, std::string hName );

    std::vector< double > * generateVtx3X( TH3D * prodVtxHist );

  } // namespace HNLFluxReader

} // namespace HNL
} // namespace genie

#endif // #ifndef _HNL_FLUXREADER_H_
