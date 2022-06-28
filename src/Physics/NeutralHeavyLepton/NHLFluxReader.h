//----------------------------------------------------------------------------
/*!

  This reads in flux histograms and constructs SimpleNHL objects
  with appropriate energy + momentum + vertex position + pol

  The couplings to SM are not relevant at this stage
  However, labelling the flux with its couplings should be so we don't lose track

  Solution: Read in couplings directly from the *file name*. 

\namespace  genie::NHL::NHLFluxReader

\brief      Reads in user inout flux from histograms

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    January 5th, 2022

\cpright    ??? - TBD

 */
//----------------------------------------------------------------------------
// TODO: Calculate polarisation vector magnitude in 3-body NHL production
//----------------------------------------------------------------------------

#ifndef _NHL_FLUXREADER_H_
#define _NHL_FLUXREADER_H_

// -- ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"

// -- GENIE includes

#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
//#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/NeutralHeavyLepton/NHLEnums.h"
#include "Physics/NeutralHeavyLepton/NHLKinUtils.h"
#include "Physics/NeutralHeavyLepton/SimpleNHL.h"

namespace genie{
namespace NHL{

  class SimpleNHL;

  namespace NHLFluxReader{

    extern std::string fPath; // flux file to poke
    extern double fmN; // NHL mass

    // RETHERE: will the following line behave?
    //static const double kKaonMass = ( PDGLibrary::Instance() )->Find( genie::kPdgKP )->Mass();
    static const double kKaonMass = 0.4936 * genie::units::GeV; // to match pdg table
    
    /// initialise bare NHL parameters
    inline void bareInit( ){ fmN = 0.; }
    
    inline void setMass( const double mN ){ fmN = mN; }

    /// perform selections from input files
    int selectMass( const double mN );
    void selectFile( std::string fin, const double mN ); // find but don't open the file

    /// get the histogram and energy from it
    TH1F * getFluxHist1F( std::string fin, std::string hName, int HType );
    // and overloaded method for the newer interface
    TH1F * getFluxHist1F( std::string fin, int masspoint, bool isParticle );

    TH3D * getFluxHist3D( std::string fin, std::string dirName, std::string hName );

    std::vector< double > * generateVtx3X( TH3D * prodVtxHist );
    
    /// interface to SimpleNHL constructor
    genie::NHL::SimpleNHL generateNHL( const int PDG, const int parPDG, const double mN,
				       const double Ue42, const double Um42, const double Ut42 );

  } // namespace NHLFluxReader

} // namespace NHL
} // namespace genie

#endif // #ifndef _NHL_FLUXREADER_H_
