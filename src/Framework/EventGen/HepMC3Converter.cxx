#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_HEPMC3_INTERFACE_ENABLED__
//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// ROOT includes
#include "TCollection.h"

// GENIE includes
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GVersion.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/HepMC3Converter.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Framework/Utils/RunOpt.h"

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/WriterAscii.h"

// GHepParticles return units of GeV/c for p.  the V_i are all in fermis
// and are relative to the center of the struck nucleus.
// add the vertex X/Y/Z to the V_i for status codes 0 and 1

// Definitions unique to this source file
namespace {
  constexpr int DUMMY_PARTICLE_INDEX = -1;

  // Implemented version of the NuHepMC standard
  // (https://github.com/NuHepMC/Spec)
  constexpr int NUHEPMC_MAJOR_VERSION = 0;
  constexpr int NUHEPMC_MINOR_VERSION = 1;
  constexpr int NUHEPMC_PATCH_VERSION = 0;

  // E.C.1
  // Mapping of GENIE EventGeneratorList names to NuHepMC3 process ID codes.
  // Note that these are distinct from GENIE's own process labeling scheme.
  //
  // NOTE: NuHepMC 0.1.0 does not designate ID conventions for EM processes.
  // They are assigned negative codes here for the time being.
  const std::map< std::string, int > NUHEPMC3_PROC_MAP = {
    { "AM-NUGAMMA", 751 }, // "kScAMNuGamma, kIntWeakNC"
    { "COH-CC-PION", 100 },
    { "COH-NC-PION", 150 },
    { "DFR-CC", 700 },
    { "DFR-NC", 750 },
    { "DIS-CC", 600 },
    { "DIS-CC-CHARM", 601 },
    { "DIS-CC-SINGLEK", 602 },
    { "DIS-EM", -600 },
    { "DIS-NC", 650 },
    { "IMD", 701 },
    { "IMD-ANH", 702 },
    { "MEC-CC", 300 },
    { "MEC-EM", -300 },
    { "MEC-NC", 350 },
    { "NNBarOsc", 900 },
    { "NUE-EL", 703 },
    { "NucleonDecay", 901 },
    { "QEL-CC", 200 },
    { "QEL-CC-CHARM", 201 },
    { "QEL-CC-LAMBDA", 202 },
    { "QEL-EM", -200 },
    { "QEL-NC", 250 },
    { "RES-CC", 400 },
    { "RES-EM", -400 },
    { "RES-NC", 450 },
  };

  // P.R.1
  // Mapping of GENIE particle status codes to NuHepMC-compliant codes. The
  // correspondence isn't purely one-to-one.
  struct PartStatusInfo {
    PartStatusInfo() {}
    PartStatusInfo( int code, const std::string& name,
      const std::string& desc ) : code_( code ), name_( name ), desc_( desc ) {}

    int code_;
    std::string name_;
    std::string desc_;
  };

  const std::map< genie::GHepStatus_t, PartStatusInfo >
    NUHEPMC_PARTICLE_STATUS_MAP
  {
    { genie::EGHepStatus::kIStUndefined,
      { 0, "Not defined", "Not meaningful" } },

    // NuHepMC divides the initial state into "incoming beam" and "target"
    // particles
    //{ genie::EGHepStatus::kIStInitialState,
    //  { 4, "Incoming beam particle", "Not meaningful" } },

    { genie::EGHepStatus::kIStStableFinalState,
      { 1, "Final state", "Undecayed physical particle" } },
    { genie::EGHepStatus::kIStIntermediateState,
      { 21, "Intermediate state", "Temporary particle for internal use" } },
    { genie::EGHepStatus::kIStDecayedState,
      { 2, "Decayed state", "Decayed physical particle" } },
    { genie::EGHepStatus::kIStCorrelatedNucleon,
      { 22, "Correlated nucleon", "Spectator nucleon in a correlated pair" } },
    { genie::EGHepStatus::kIStNucleonTarget,
      { 23, "Target nucleon", "Struck nucleon in the initial state" } },
    { genie::EGHepStatus::kIStDISPreFragmHadronicState,
      { 24, "Prefragmentation", "Temporary prefragmentation hadronic state"
        " for deep inelastic scattering" } },
    { genie::EGHepStatus::kIStPreDecayResonantState,
      { 25, "Pre-decay resonance", "Temporary hadronic resonance"
        " before decay" } },
    { genie::EGHepStatus::kIStHadronInTheNucleus,
      { 26, "Hadron in the nucleus",
        "Input particle for intranuclear cascade" } },
    { genie::EGHepStatus::kIStFinalStateNuclearRemnant,
      { 27, "Hadronic blob", "Pseudoparticle representing the"
        " final-state remnant nucleus" } },
    { genie::EGHepStatus::kIStNucleonClusterTarget,
      { 28, "Nucleon cluster", "Temporary multi-nucleon system for"
        " internal use" } },
  };

}

//____________________________________________________________________________
genie::HepMC3Converter::HepMC3Converter()
{
  fWriter = std::make_shared< HepMC3::WriterAscii >( "example.hepmc3" );
}
//____________________________________________________________________________
std::shared_ptr< HepMC3::GenEvent > genie::HepMC3Converter::ConvertToHepMC3(
  const genie::EventRecord& gevrec )
{
  // TODO: check GENIE's unit conventions
  // E.R.3
  auto evt = std::make_shared< HepMC3::GenEvent >( HepMC3::Units::GEV,
    HepMC3::Units::CM );

  // Create the primary vertex
  // E.R.5
  // TODO: add vertex positions
  auto prim_vtx = std::make_shared< HepMC3::GenVertex >();
  prim_vtx->set_status( 1 );

  evt->add_vertex( prim_vtx );

  // Add the particles from the GENIE event
  TIter g_part_iter( &gevrec );
  genie::GHepParticle* g_part = nullptr;

  int g_part_idx = 0;
  while( (g_part = dynamic_cast< genie::GHepParticle* >(g_part_iter.Next())) ) {

    HepMC3::FourVector mom4( g_part->Px(), g_part->Py(), g_part->Pz(),
      g_part->Energy() );

    int hepmc3_status = this->GetHepMC3ParticleStatus( g_part, gevrec );

    auto part = std::make_shared< HepMC3::GenParticle >( mom4, g_part->Pdg(),
      hepmc3_status );

    // Primary particles have no mother
    int first_mommy = g_part->FirstMother();
    if ( first_mommy == DUMMY_PARTICLE_INDEX ) {
      // Initial-state primary particles go into the primary vertex
      if ( g_part->Status() == genie::EGHepStatus::kIStInitialState ) {
        prim_vtx->add_particle_in( part );
      }
      // In unusual cases where they are not labeled as initial-state,
      // have them come out of the primary vertex
      else {
        // TODO: revisit this
        prim_vtx->add_particle_out( part );
      }
    }
    // Other particles come out of a vertex
    else {
      const auto& part_vec = evt->particles();
      auto mother_part = part_vec.at( first_mommy );
      auto production_vtx = mother_part->end_vertex();

      production_vtx->add_particle_out( part );
      if ( production_vtx != prim_vtx ) {
        // TODO: set vertex status based on the current particle status
        production_vtx->set_status( 3 );
      }
    }

    // All particles except those without daughters have an end vertex.
    int first_daughter = g_part->FirstDaughter();
    if ( first_daughter != DUMMY_PARTICLE_INDEX
      // Don't create a new end vertex for the primary particles since the
      // primary vertex already exists
      && first_mommy != DUMMY_PARTICLE_INDEX )
    {
      // TODO: add vertex positions

      // Handle cases with multiple mothers. If the first daughter of the
      // current GENIE particle has a first mother which appears before the
      // current GENIE particle in the event record, then retrieve the
      // end vertex which was already made. Otherwise, make a new end vertex.
      genie::GHepParticle* fd_part = gevrec.Particle( first_daughter );
      int fd_fm = fd_part->FirstMother();

      auto end_vtx = std::make_shared< HepMC3::GenVertex >();
      if ( fd_fm < g_part_idx ) {
        end_vtx = evt->particles().at( fd_fm )->end_vertex();
      }
      else {
        evt->add_vertex( end_vtx );
      }

      end_vtx->add_particle_in( part );
    }

    // Add the completed HepMC3::GenParticle object to the event
    evt->add_particle( part );
    // TODO: add particle attributes

    // Move to the next particle in the GENIE event record
    ++g_part_idx;
  }

  // Create a HepMC3::GenRunInfo object if needed, then associate it with the
  // current event
  if ( !fRunInfo ) this->PrepareRunInfo();
  evt->set_run_info( fRunInfo );

  // Return the completed HepMC3::GenEvent object
  return evt;
}
//--------------------------------------------------------------------------
int genie::HepMC3Converter::GetHepMC3ParticleStatus(
  const genie::GHepParticle* /*gpart*/, const genie::EventRecord& /*gevrec*/ ) const
{
  // TODO: write this
  return 4;
}
//--------------------------------------------------------------------------
void genie::HepMC3Converter::WriteEvent( const HepMC3::GenEvent& evt ) const
{
  fWriter->write_event( evt );
}
//--------------------------------------------------------------------------
void genie::HepMC3Converter::PrepareRunInfo()
{
  // G.R.1
  fRunInfo = std::make_shared<HepMC3::GenRunInfo>();

  // G.R.2
  fRunInfo->add_attribute( "NuHepMC.Version.Major",
    std::make_shared<HepMC3::IntAttribute>(NUHEPMC_MAJOR_VERSION) );

  fRunInfo->add_attribute( "NuHepMC.Version.Minor",
    std::make_shared<HepMC3::IntAttribute>(NUHEPMC_MINOR_VERSION) );

  fRunInfo->add_attribute( "NuHepMC.Version.Patch",
    std::make_shared<HepMC3::IntAttribute>(NUHEPMC_PATCH_VERSION) );

  // G.R.3
  fRunInfo->tools().emplace_back(
    HepMC3::GenRunInfo::ToolInfo{ "GENIE", __GENIE_RELEASE__,
      __GENIE_GIT_REVISION__ }
  );

  // Store the CMC identifier ("tune") as a GENIE-specific attribute of the run
  // info
  genie::RunOpt* ro = genie::RunOpt::Instance();
  // NOTE: This assumes that genie::RunOpt::BuildTune() has been called
  // previously
  std::string tune_name = ro->Tune()->Name();
  fRunInfo->add_attribute( "GENIE.XSecTune",
    std::make_shared<HepMC3::StringAttribute>(tune_name) );

  // G.R.4
  std::vector< int > proc_IDs;

  // Build the list of process IDs from the global list of xsec models
  genie::AlgConfigPool* acp = genie::AlgConfigPool::Instance();
  genie::Registry* gpl = acp->GlobalParameterList();
  const auto& gpl_item_map = gpl->GetItemMap();
  const std::string prefix( "XSecModel@genie::EventGenerator/" );
  for ( const auto& pair : gpl_item_map ) {
    // Skip keys that aren't cross-section model settings
    const std::string& key = pair.first;
    if ( key.find(prefix) == std::string::npos ) continue;

    // Get the interaction mode string (the part of the key after the prefix)
    std::string mode = key;
    mode.erase( 0, prefix.size() );

    // Look up the NuHepMC3 process ID code for the current kind of interaction
    int id_code = NUHEPMC3_PROC_MAP.at( mode );
    std::string id_code_str = std::to_string( id_code );

    // Add it to the vector of active process ID codes
    proc_IDs.push_back( id_code );

    // Save the name of the process as a string attribute
    fRunInfo->add_attribute( "NuHepMC.ProcessInfo[" + id_code_str + "].Name",
      std::make_shared< HepMC3::StringAttribute >(mode) );

    // Get the registry algorithm ID for the cross-section model used by GENIE
    RgAlg alg = gpl->GetAlg( key );

    // Look up its XML configuration settings
    genie::Registry* xsec_config = acp->FindRegistry( alg.name, alg.config );

    // G.C.5
    // Check the number of citations listed for the given algorithm
    //int num_citations = xsec_config->GetIntDef( "Citation-Count", 0 );

    // Get the DOI for the first citation for this algorithm, if it exists
    std::string doi = xsec_config->GetStringDef( "Citation-0-DOI", "" );

    // Try using the default configuration of the algorithm if a DOI wasn't
    // found and the algorithm has a non-default configuration
    if ( doi.empty() && alg.config != "Default" ) {
      genie::Registry* xsec_config_def = acp->FindRegistry( alg.name,
        "Default" );
      doi = xsec_config_def->GetStringDef( "Citation-0-DOI", "" );
    }

    std::string cite;
    if ( !doi.empty() ) {
      cite = " -- http://doi.org/" + doi;
    }

    std::string desc =  alg.name + '/' + alg.config + cite;

    // Save the description of the process as a string attribute
    fRunInfo->add_attribute( "NuHepMC.ProcessInfo[" + id_code_str
      + "].Description", std::make_shared< HepMC3::StringAttribute >(desc) );
  }

  fRunInfo->add_attribute( "NuHepMC.ProcessIDs",
    std::make_shared< HepMC3::VectorIntAttribute >(proc_IDs) );

  // G.R.5
  std::vector< int > vertex_IDs;

  const std::map< int, std::pair<std::string, std::string> > vtx_status_map = {
    { 1, { "Primary", "The primary vertex or hard scatter" } },
    { 3, { "Secondary", "Secondary vertex" } },
  };

  for ( const auto& vtx_pair : vtx_status_map ) {
    int vtx_code = vtx_pair.first;
    vertex_IDs.push_back( vtx_code );

    std::string vtx_prefix = "NuHepMC.VertexStatusInfo["
      + std::to_string( vtx_code );
    fRunInfo->add_attribute( vtx_prefix + "].Name",
      std::make_shared< HepMC3::StringAttribute >(vtx_pair.second.first) );
    fRunInfo->add_attribute( vtx_prefix + "].Description",
      std::make_shared< HepMC3::StringAttribute >(vtx_pair.second.second) );
  }

  fRunInfo->add_attribute( "NuHepMC.VertexStatusIDs",
    std::make_shared< HepMC3::VectorIntAttribute >(vertex_IDs) );

  // G.R.6
  std::vector< int > particle_statuses = { 4, 11 };

  fRunInfo->add_attribute( "NuHepMC.ParticleStatusInfo[4].Name",
    std::make_shared< HepMC3::StringAttribute >("Beam") );

  fRunInfo->add_attribute( "NuHepMC.ParticleStatusInfo[4].Description",
    std::make_shared< HepMC3::StringAttribute >("Incoming beam particle") );

  fRunInfo->add_attribute( "NuHepMC.ParticleStatusInfo[11].Name",
    std::make_shared< HepMC3::StringAttribute >("Target") );

  fRunInfo->add_attribute( "NuHepMC.ParticleStatusInfo[11].Description",
    std::make_shared< HepMC3::StringAttribute >("Target particle") );

  for ( const auto& pstatus_pair : NUHEPMC_PARTICLE_STATUS_MAP ) {
    const auto& ps = pstatus_pair.second;
    std::string ps_code_str = std::to_string( ps.code_ );

    particle_statuses.push_back( ps.code_ );

    fRunInfo->add_attribute( "NuHepMC.ParticleStatusInfo[" + ps_code_str
      + "].Name", std::make_shared< HepMC3::StringAttribute >(ps.name_) );

    fRunInfo->add_attribute( "NuHepMC.ParticleStatusInfo[" + ps_code_str
      + "].Description",
      std::make_shared< HepMC3::StringAttribute >(ps.desc_) );
  }

  fRunInfo->add_attribute( "NuHepMC.ParticleStatusIDs",
    std::make_shared< HepMC3::VectorIntAttribute >(particle_statuses) );

  // G.R.7
  fRunInfo->set_weight_names( { "CV" } );

  // G.C.1
  fRunInfo->add_attribute( "NuHepMC.Conventions",
      std::make_shared<HepMC3::VectorStringAttribute>(
          std::vector<std::string>{"G.C.1", "G.C.2", "G.C.4", "G.C.5", "E.C.1",
                                   "E.C.2", "E.C.3", "E.C.5", "E.C.6"}));

  //// G.C.2
  //fRunInfo->add_attribute("NuHepMC.Exposure.NEvents",
  //                        std::make_shared<HepMC3::IntAttribute>(3));

  //// G.C.4
  //fRunInfo->add_attribute(
  //    "NuHepMC.FluxAveragedTotalCrossSection",
  //    std::make_shared<HepMC3::DoubleAttribute>(1.234E-38 * cm2_to_pb));
}
//--------------------------------------------------------------------------

#endif //__GENIE_HEPMC3_INTERFACE_ENABLED__
