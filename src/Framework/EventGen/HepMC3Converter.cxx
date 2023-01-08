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

// Standard library includes
#include <algorithm>
#include <sstream>

// ROOT includes
#include "TCollection.h"

// GENIE includes
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GVersion.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/HepMC3Converter.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Framework/Utils/RunOpt.h"

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

// Definitions unique to this source file
namespace {
  constexpr int DEFAULT_RESCATTER_CODE = -1;
  constexpr int DEFAULT_DECAY_MODE = -1;
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
      { 23, "Intermediate state", "Temporary particle for internal use" } },
    { genie::EGHepStatus::kIStDecayedState,
      { 2, "Decayed state", "Decayed physical particle" } },
    { genie::EGHepStatus::kIStCorrelatedNucleon,
      { 22, "Correlated nucleon", "Spectator nucleon in a correlated pair" } },
    { genie::EGHepStatus::kIStNucleonTarget,
      { 21, "Target nucleon", "Struck nucleon in the initial state" } },
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
    { genie::EGHepStatus::kIStFormZone,
      { 29, "Formation zone", "Hadron before formation zone free step" } },
  };

  // Convert the contents of a TBits object into a string that can be stored in
  // a HepMC3::StringAttribute
  std::string tbits_to_string( const TBits& bits ) {
    std::stringstream temp_ss;
    temp_ss << bits;
    return temp_ss.str();
  }

  // Set the contents of a TBits object from a string. This is used to retrieve
  // a HepMC3::StringAttribute originally converted from a TBits object via
  // tbits_to_string()
  void set_tbits_from_string( const std::string& str, TBits& bits ) {
    std::stringstream temp_ss( str );
    size_t num_bits = str.size();
    bool bit_val;
    for ( size_t b = 0u; b < num_bits; ++b ) {
      temp_ss >> bit_val;
      bits.SetBitNumber( b, bit_val );
    }
  }

  // Converts a TLorentzVector to a form suitable for storage as a HepMC3
  // attribute
  std::shared_ptr< HepMC3::VectorDoubleAttribute > four_vector_to_attribute(
    const TLorentzVector& vec4 )
  {
    std::vector< double > temp_vec = { vec4.X(), vec4.Y(), vec4.Z(), vec4.T() };
    return std::make_shared< HepMC3::VectorDoubleAttribute >( temp_vec );
  }

  // Retrieves a TLorentzVector stored in a HepMC3::VectorDoubleAttribute
  TLorentzVector attribute_to_four_vector(
    const HepMC3::VectorDoubleAttribute& attr )
  {
    const std::vector< double >& vec = attr.value();
    TLorentzVector result( vec.at(0), vec.at(1), vec.at(2), vec.at(3) );
    return result;
  }

}
//____________________________________________________________________________
genie::HepMC3Converter::HepMC3Converter()
{

}
//____________________________________________________________________________
std::shared_ptr< HepMC3::GenEvent > genie::HepMC3Converter::ConvertToHepMC3(
  const genie::EventRecord& gevrec )
{
  // TODO: check GENIE's unit conventions
  // E.R.3
  auto evt = std::make_shared< HepMC3::GenEvent >( HepMC3::Units::GEV,
    HepMC3::Units::CM );

  // Set the overall event 4-position using the GHepRecord vertex
  const TLorentzVector* g_vtx = gevrec.Vertex();
  HepMC3::FourVector evt_pos4( g_vtx->X(), g_vtx->Y(), g_vtx->Z(), g_vtx->T() );
  evt->shift_position_to( evt_pos4 );

  // Create the primary vertex
  // E.R.5
  auto prim_vtx = std::make_shared< HepMC3::GenVertex >();
  prim_vtx->set_status( 1 );

  evt->add_vertex( prim_vtx );

  // Indices for especially important particles in the GENIE event (these are
  // equal to DUMMY_PARTICLE_INDEX in cases where the relevant particle doesn't
  // exist)
  int probe_idx = gevrec.ProbePosition();

  // DUMMY_PARTICLE_INDEX for a free nucleon target
  int tgt_idx = gevrec.TargetNucleusPosition();

  // DUMMY_PARTICLE_INDEX for COH, nu-e, etc.
  int hit_nucleon_idx = gevrec.HitNucleonPosition();

  // Set the 4-position of the primary vertex using a reference particle from
  // the GENIE event. Prefer taking coordinates from the probe, hit nucleon, and
  // target nucleus (in that order). Complain if none of those are available.
  genie::GHepParticle* ref_part = nullptr;
  if ( probe_idx != DUMMY_PARTICLE_INDEX ) {
    ref_part = gevrec.Probe();
  }
  else if ( hit_nucleon_idx != DUMMY_PARTICLE_INDEX ) {
    ref_part = gevrec.HitNucleon();
  }
  else if ( tgt_idx != DUMMY_PARTICLE_INDEX ) {
    ref_part = gevrec.TargetNucleus();
  }
  else {
    LOG( "HepMC3Converter", pFATAL ) << "Could not set a valid 4-position"
      << " for the primary vertex!";
    std::exit( 1 );
  }

  HepMC3::FourVector prim_vtx_pos4( ref_part->Vx(), ref_part->Vy(),
    ref_part->Vz(), ref_part->Vt() );

  prim_vtx->set_position( prim_vtx_pos4 );

  // In events with both a complex nuclear target and a defined initial-state
  // struck nucleon, define a "nuclear vertex" separating the nucleon target
  // from the spectator remnant nucleus
  std::shared_ptr< HepMC3::GenVertex > nuclear_vtx;
  if ( tgt_idx != DUMMY_PARTICLE_INDEX &&
    hit_nucleon_idx != DUMMY_PARTICLE_INDEX )
  {
    nuclear_vtx = std::make_shared< HepMC3::GenVertex >();
    nuclear_vtx->set_status( 2 );
    evt->add_vertex( nuclear_vtx );

    // Get the 4-position of the nuclear vertex from the target nucleus
    // GHepParticle
    genie::GHepParticle* tgt = gevrec.TargetNucleus();
    HepMC3::FourVector nuclear_vtx_pos4( tgt->Vx(), tgt->Vy(),
      tgt->Vz(), tgt->Vt() );
    nuclear_vtx->set_position( nuclear_vtx_pos4 );
  }

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
      // Initial-state primary particles go into either the primary vertex
      // or, in the case of a complex nuclear target with a defined
      // initial-state struck nucleon, the nuclear vertex
      if ( g_part->Status() == genie::EGHepStatus::kIStInitialState ) {
        if ( nuclear_vtx && g_part_idx == tgt_idx ) {
          nuclear_vtx->add_particle_in( part );
        }
        else {
          prim_vtx->add_particle_in( part );
        }
      }
      // In unusual cases where the primary particle is not labeled as
      // initial-state, have it come out of the nuclear vertex (if it
      // exists) or the primary vertex (if it does not)
      // NOTE: This takes care of the "binding energy" pseudoparticles
      // made by GENIE for some interaction modes
      // TODO: revisit this
      else {
        if ( nuclear_vtx ) {
          nuclear_vtx->add_particle_out( part );
        }
        else {
          prim_vtx->add_particle_out( part );
        }
      }
    }
    // Other particles come out of a vertex
    else {
      const auto& part_vec = evt->particles();
      auto mother_part = part_vec.at( first_mommy );
      auto production_vtx = mother_part->end_vertex();

      production_vtx->add_particle_out( part );

      // Set the vertex status if it hasn't already been defined
      // TODO: set vertex status based on the current particle status
      if ( production_vtx->status() != 0 ) {
        production_vtx->set_status( 3 );
      }

      // In the case of an initial-state struck nucleon, send it into the
      // primary vertex
      if ( g_part_idx == hit_nucleon_idx ) {
        prim_vtx->add_particle_in( part );
      }
    }

    // All particles except those without daughters have an end vertex.
    int first_daughter = g_part->FirstDaughter();
    if ( first_daughter != DUMMY_PARTICLE_INDEX
      // Don't create a new end vertex for the primary particles since the
      // primary vertex already exists
      && first_mommy != DUMMY_PARTICLE_INDEX )
    {
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

        // This is a new vertex, so get its 4-position from the first daughter
        // of the current particle
        HepMC3::FourVector end_vtx_pos4( fd_part->Vx(), fd_part->Vy(),
          fd_part->Vz(), fd_part->Vt() );
        end_vtx->set_position( end_vtx_pos4 );
      }

      end_vtx->add_particle_in( part );
    }

    // Add the completed HepMC3::GenParticle object to the event
    evt->add_particle( part );

    // Set the rescatter code and "is bound" attributes if they differ from
    // their default values
    int resc = g_part->RescatterCode();
    if ( resc != DEFAULT_RESCATTER_CODE ) {
      part->add_attribute( "GENIE.RescatterCode",
        std::make_shared< HepMC3::IntAttribute >(resc) );
    }

    // bool <--> int conversions are implicit
    bool bound = g_part->IsBound();
    if ( bound ) {
      part->add_attribute( "GENIE.IsBound",
        std::make_shared< HepMC3::IntAttribute >(bound) );
    }

    // If a polarization is defined for the particle, then set extra attributes
    // to store it
    if ( g_part->PolzIsSet() ) {
      double polz_polar_ang = g_part->PolzPolarAngle();
      part->add_attribute( "GENIE.PolzPolarAngle",
        std::make_shared< HepMC3::DoubleAttribute >(polz_polar_ang) );

      double polz_azim_ang = g_part->PolzAzimuthAngle();
      part->add_attribute( "GENIE.PolzAzimuthAngle",
        std::make_shared< HepMC3::DoubleAttribute >(polz_azim_ang) );
    }

    // If a removal energy is defined for the particle, then set an extra
    // attribute to store it
    double E_rem = g_part->RemovalEnergy();
    if ( E_rem != 0. ) {
      part->add_attribute( "GENIE.RemovalEnergy",
        std::make_shared< HepMC3::DoubleAttribute >(E_rem) );
    }

    // Move to the next particle in the GENIE event record
    ++g_part_idx;
  }

  // Set the overall event weight
  double wgt = gevrec.Weight();
  evt->weights().push_back( wgt );

  // Set string attributes to store the boolean event flags and mask
  std::string flags = tbits_to_string( *gevrec.EventFlags() );
  std::string mask = tbits_to_string( *gevrec.EventMask() );

  evt->add_attribute( "GENIE.EventFlags",
    std::make_shared< HepMC3::StringAttribute >(flags) );
  evt->add_attribute( "GENIE.EventMask",
    std::make_shared< HepMC3::StringAttribute >(mask) );

  // Set attributes for other GHepRecord data members
  double prob = gevrec.Probability();
  double xsec = gevrec.XSec();
  double diff_xsec = gevrec.DiffXSec();
  int phase_space = static_cast< int >( gevrec.DiffXSecVars() );

  evt->add_attribute( "GENIE.Probability",
    std::make_shared< HepMC3::DoubleAttribute >(prob) );
  evt->add_attribute( "GENIE.XSec",
    std::make_shared< HepMC3::DoubleAttribute >(xsec) );
  evt->add_attribute( "GENIE.DiffXSec",
    std::make_shared< HepMC3::DoubleAttribute >(diff_xsec) );
  evt->add_attribute( "GENIE.DiffXSecVars",
    std::make_shared< HepMC3::IntAttribute >(phase_space) );

  // Store the interaction summary in the HepMC3 event
  this->StoreInteraction( *gevrec.Summary(), *evt );

  // Create a HepMC3::GenRunInfo object if needed, then associate it with the
  // current event
  if ( !fRunInfo ) this->PrepareRunInfo();
  evt->set_run_info( fRunInfo );

  // Return the completed HepMC3::GenEvent object
  return evt;
}
//____________________________________________________________________________
int genie::HepMC3Converter::GetHepMC3ParticleStatus(
  const genie::GHepParticle* gpart, const genie::EventRecord& gevrec ) const
{
  // The initial state status is split by the NuHepMC standard into "beam"
  // (probe) and "target" particles. Decide which to use here if needed.
  genie::GHepStatus_t status = gpart->Status();
  if ( status == genie::EGHepStatus::kIStInitialState ) {
    genie::GHepParticle* probe = gevrec.Probe();
    if ( gpart == probe ) return 4; // NuHepMC beam particle
    else return 11; // NuHepMC target particle
  }

  // Otherwise, there is a one-to-one mapping of GENIE codes to NuHepMC
  // codes. Look up the conversion here.
  auto end = NUHEPMC_PARTICLE_STATUS_MAP.end();
  auto iter = NUHEPMC_PARTICLE_STATUS_MAP.find( status );

  // If the lookup was unsuccessful, then complain
  if ( iter == end ) {
    LOG( "HepMC3Converter", pFATAL ) << "Could not convert GENIE particle"
      " status code!";
    std::exit( 1 );
  }

  return iter->second.code_;
}
//____________________________________________________________________________
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
//____________________________________________________________________________
std::shared_ptr< genie::EventRecord > genie::HepMC3Converter::RetrieveGHEP(
  const HepMC3::GenEvent& evt )
{
  auto gevrec = std::make_shared< genie::EventRecord >();

  // Retrieve and store the overall event weight
  double wgt = evt.weight();
  gevrec->SetWeight( wgt );

  const auto& part_vec = evt.particles();
  for ( const auto& part : part_vec ) {
    int pdg = part->pid();
    genie::GHepStatus_t status = this->GetGHepParticleStatus( part->status() );
    const HepMC3::FourVector& p4 = part->momentum();

    // The particle 4-position will be retrieved from the production or end
    // vertex below and stored in this 4-vector
    HepMC3::FourVector x4;

    int mommy1 = DUMMY_PARTICLE_INDEX;
    int mommy2 = DUMMY_PARTICLE_INDEX;
    const auto& prod_vtx = part->production_vertex();
    if ( prod_vtx && prod_vtx->id() != 0 ) {
      x4 = prod_vtx->position();
      const auto& mommy_vec = prod_vtx->particles_in();
      size_t mommy_count = mommy_vec.size();
      if ( mommy_count > 0u ) mommy1 = mommy_vec.front()->id() - 1;
      if ( mommy_count > 1u ) mommy2 = mommy_vec.back()->id() - 1;

      // Nuclear binding energy pseudoparticles are recorded in the GENIE event
      // record as if they were primary (motherless). Ignore the vertex
      // relationships in this special case.
      if ( pdg == genie::kPdgBindino ) {
        mommy1 = DUMMY_PARTICLE_INDEX;
        mommy2 = DUMMY_PARTICLE_INDEX;
      }
    }

    int dau1 = DUMMY_PARTICLE_INDEX;
    int dau2 = DUMMY_PARTICLE_INDEX;
    const auto& end_vtx = part->end_vertex();
    if ( end_vtx ) {
      // If we don't have a production vertex, then use the end vertex as a
      // backup for assigning a 4-position to the current GHepParticle. Also
      // use the end vertex (which will be the primary vertex) to assign a
      // position to the initial hit nucleon when scattering on a complex
      // target.
      if ( !prod_vtx || prod_vtx->id() == 0 || status == kIStNucleonTarget ) {
        x4 = end_vtx->position();
      }
      const auto& dau_vec = end_vtx->particles_out();
      size_t dau_count = dau_vec.size();
      for ( size_t d = 0u; d < dau_count; ++d ) {
        const auto& daughter = dau_vec.at( d );
        // Skip "bindinos" since GENIE treats them as motherless
        if ( daughter->pid() != genie::kPdgBindino ) {
          if ( dau1 == DUMMY_PARTICLE_INDEX) {
            dau1 = daughter->id() - 1;
            dau2 = dau1;
          }
          else dau2 = daughter->id() - 1;
        }
      }
    }

    gevrec->AddParticle( pdg, status, mommy1, mommy2, dau1, dau2, p4.px(),
      p4.py(), p4.pz(), p4.e(), x4.x(), x4.y(), x4.z(), x4.t() );

    // Get a pointer to the newly-created GHepParticle object so that we can
    // set additional data members using HepMC3 attributes if they are present
    genie::GHepParticle* g_part = gevrec->Particle( gevrec->GetEntries() - 1 );

    auto resc_ptr = part->attribute< HepMC3::IntAttribute >(
      "GENIE.RescatterCode" );
    if ( resc_ptr ) g_part->SetRescatterCode( resc_ptr->value() );

    // bool <--> int conversions are implicit
    auto bound_ptr = part->attribute< HepMC3::IntAttribute >(
      "GENIE.IsBound" );
    if ( bound_ptr ) g_part->SetBound( bound_ptr->value() );

    // TODO: add error handling for when only one of these is set
    auto polz_pol_ptr = part->attribute< HepMC3::DoubleAttribute >(
      "GENIE.PolzPolarAngle" );
    auto polz_azm_ptr = part->attribute< HepMC3::DoubleAttribute >(
      "GENIE.PolzAzimuthAngle" );

    if ( polz_pol_ptr && polz_azm_ptr ) {
      g_part->SetPolarization( polz_pol_ptr->value(), polz_azm_ptr->value() );
    }

    auto E_rem_ptr = part->attribute< HepMC3::DoubleAttribute >(
      "GENIE.RemovalEnergy" );
    if ( E_rem_ptr ) g_part->SetRemovalEnergy( E_rem_ptr->value() );
  }

  // Retrieve event attributes used to store extra data members
  auto flags_ptr = evt.attribute< HepMC3::StringAttribute >(
    "GENIE.EventFlags" );
  if ( flags_ptr ) set_tbits_from_string( flags_ptr->value(),
    *gevrec->EventFlags() );

  auto mask_ptr = evt.attribute< HepMC3::StringAttribute >(
    "GENIE.EventMask" );
  if ( mask_ptr ) set_tbits_from_string( mask_ptr->value(),
    *gevrec->EventMask() );

  auto prob_ptr = evt.attribute< HepMC3::DoubleAttribute >(
    "GENIE.Probability" );
  if ( prob_ptr ) gevrec->SetProbability( prob_ptr->value() );

  auto xsec_ptr = evt.attribute< HepMC3::DoubleAttribute >(
    "GENIE.XSec" );
  if ( xsec_ptr ) gevrec->SetXSec( xsec_ptr->value() );

  // TODO: add error handling if only one of these is set
  auto diff_xsec_ptr = evt.attribute< HepMC3::DoubleAttribute >(
    "GENIE.DiffXSec" );
  auto phase_space_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.DiffXSecVars" );

  if ( diff_xsec_ptr && phase_space_ptr ) {
    auto ps = static_cast< genie::KinePhaseSpace_t >(
      phase_space_ptr->value() );
    gevrec->SetDiffXSec( diff_xsec_ptr->value(), ps );
  }

  genie::Interaction* itr = this->RetrieveInteraction( evt );
  gevrec->AttachSummary( itr );

  // For reactions that do not simultaneously contain a complex target nucleus
  // and an initial struck nucleon, prune the daughters of the primary particles
  // to match the GHepRecord conventions
  genie::GHepParticle* tgt = gevrec->TargetNucleus();
  genie::GHepParticle* hit_nuc = gevrec->HitNucleon();
  if ( !tgt || !hit_nuc ) {
    genie::GHepParticle* probe = gevrec->Probe();
    genie::GHepParticle* fs_lep = gevrec->FinalStatePrimaryLepton();
    if ( probe && fs_lep ) {
      int pr_pos = gevrec->ProbePosition();
      fs_lep->SetFirstMother( pr_pos );
      fs_lep->SetLastMother( DUMMY_PARTICLE_INDEX );

      int fsl_pos = gevrec->FinalStatePrimaryLeptonPosition();
      int pr_dau1 = probe->FirstDaughter();
      int pr_dau2 = probe->LastDaughter();
      if ( pr_dau2 != DUMMY_PARTICLE_INDEX ) {
        for ( int d = pr_dau1; d <= pr_dau2; ++d ) {
          genie::GHepParticle* dau = gevrec->Particle( d );
          if ( d == fsl_pos ) continue;
          dau->SetFirstMother( dau->LastMother() );
          dau->SetLastMother( DUMMY_PARTICLE_INDEX );
        }
      }

      probe->SetFirstDaughter( fsl_pos );
      probe->SetLastDaughter( fsl_pos );

      genie::GHepParticle* struck = nullptr;
      genie::GHepParticle* hit_e = gevrec->HitElectron();

      if ( hit_nuc ) struck = hit_nuc;
      else if ( tgt ) struck = tgt;
      else if ( hit_e ) struck = hit_e;

      if ( struck ) {
        int struck_dau1 = struck->FirstDaughter();
        if ( struck_dau1 == fsl_pos ) {
          struck->SetFirstDaughter( struck_dau1 + 1 );
        }

        if ( struck == tgt ) {
          tgt->SetPosition( TLorentzVector() );

          // For COH
          int tgt_dau1 = tgt->FirstDaughter();
          if ( tgt_dau1 != DUMMY_PARTICLE_INDEX ) {
            genie::GHepParticle* tgt_dau = gevrec->Particle( tgt_dau1 );
            if ( tgt_dau ) {
              int tgt_dau_pdg = tgt_dau->Pdg();
              if ( genie::pdg::IsIon(tgt_dau_pdg) ) {
                tgt_dau->SetPosition( TLorentzVector() );
              }
            }
          }
        }
      }
    }

  }

  // To match the GHepRecord convention, ignore mothers beyond the first for
  // strings from Pythia
  TIter str_iter( gevrec.get() );
  genie::GHepParticle* sp = nullptr;
  while ( (sp = dynamic_cast< genie::GHepParticle* >(str_iter.Next())) ) {
    if ( sp->Pdg() == genie::kPdgString ) {
      sp->SetLastMother( DUMMY_PARTICLE_INDEX );
    }
  }

  return gevrec;
}
//
genie::GHepStatus_t genie::HepMC3Converter::GetGHepParticleStatus(
  int nuhepmc_status ) const
{
  // Both the NuHepMC "beam" and "target" particle status codes correspond to
  // the initial state status used by GENIE
  if ( nuhepmc_status == 4 || nuhepmc_status == 11 ) {
    return genie::EGHepStatus::kIStInitialState;
  }

  // Otherwise, there is a one-to-one mapping of NuHepMC codes to GENIE codes.
  // Look up the conversion here.
  auto end = NUHEPMC_PARTICLE_STATUS_MAP.end();
  auto iter = std::find_if( NUHEPMC_PARTICLE_STATUS_MAP.begin(), end,
    [nuhepmc_status](const auto& pair)
    { return pair.second.code_ == nuhepmc_status; }
  );

  // If the lookup was unsuccessful, then complain
  if ( iter == end ) {
    LOG( "HepMC3Converter", pFATAL ) << "Could not convert NuHepMC particle"
      " status code!";
    std::exit( 1 );
  }

  return iter->first;
}
//____________________________________________________________________________
void genie::HepMC3Converter::StoreInteraction( const genie::Interaction& inter,
  HepMC3::GenEvent& evt )
{
  const genie::InitialState& istate = inter.InitState();
  const genie::Target& tgt = istate.Tgt();

  evt.add_attribute( "GENIE.Interaction.ProbePDG",
    std::make_shared< HepMC3::IntAttribute >(istate.ProbePdg()) );

  TLorentzVector* temp_probe_p4 = istate.GetProbeP4( genie::kRfLab );
  evt.add_attribute( "GENIE.Interaction.ProbeP4",
    four_vector_to_attribute(*temp_probe_p4) );
  delete temp_probe_p4;

  TLorentzVector* temp_tgt_p4 = istate.GetTgtP4( genie::kRfLab );
  evt.add_attribute( "GENIE.Interaction.TargetP4",
    four_vector_to_attribute(*temp_tgt_p4) );
  delete temp_tgt_p4;

  evt.add_attribute( "GENIE.Interaction.TargetPDG",
    std::make_shared< HepMC3::IntAttribute >(tgt.Pdg()) );
  if ( tgt.HitNucIsSet() ) {
    evt.add_attribute( "GENIE.Interaction.HitNucleonPDG",
      std::make_shared< HepMC3::IntAttribute >(tgt.HitNucPdg()) );
    evt.add_attribute( "GENIE.Interaction.HitNucleonP4",
      four_vector_to_attribute(tgt.HitNucP4()) );
    evt.add_attribute( "GENIE.Interaction.HitNucleonRadius",
      std::make_shared< HepMC3::DoubleAttribute >(tgt.HitNucPosition()) );
  }
  if ( tgt.HitQrkIsSet() ) {
    evt.add_attribute( "GENIE.Interaction.HitQuarkPDG",
      std::make_shared< HepMC3::IntAttribute >(tgt.HitQrkPdg()) );
    evt.add_attribute( "GENIE.Interaction.HitSeaQuark",
      std::make_shared< HepMC3::IntAttribute >(tgt.HitSeaQrk()) );
  }

  const genie::ProcessInfo& proc_info = inter.ProcInfo();
  int s_type = static_cast< int >( proc_info.ScatteringTypeId() );
  int i_type = static_cast< int >( proc_info.InteractionTypeId() );

  evt.add_attribute( "GENIE.Interaction.ScatteringType",
    std::make_shared< HepMC3::IntAttribute >(s_type) );
  evt.add_attribute( "GENIE.Interaction.InteractionType",
    std::make_shared< HepMC3::IntAttribute >(i_type) );

  const genie::Kinematics& kine = inter.Kine();
  evt.add_attribute( "GENIE.Interaction.FSLeptonP4",
    four_vector_to_attribute(kine.FSLeptonP4()) );
  evt.add_attribute( "GENIE.Interaction.HadSystP4",
    four_vector_to_attribute(kine.HadSystP4()) );

  // Convert the map of kinematic variable values into two vectors of the same
  // size
  std::vector< int > kvar_labels;
  std::vector< double > kvar_values;
  for ( const auto kvar_pair : kine.GetMap() ) {
    int label = static_cast< int >( kvar_pair.first );
    kvar_labels.push_back( label );
    kvar_values.push_back( kvar_pair.second );
  }

  // Only bother to store the map contents if there is at least one kinematic
  // variable set
  if ( !kvar_labels.empty() ) {
    evt.add_attribute( "GENIE.Interaction.KineVarLabels",
      std::make_shared< HepMC3::VectorIntAttribute >(kvar_labels) );
    evt.add_attribute( "GENIE.Interaction.KineVarValues",
      std::make_shared< HepMC3::VectorDoubleAttribute >(kvar_values) );
  }

  // Store data members of the exclusive tag if they differ from their default
  // values
  const genie::XclsTag& xt = inter.ExclTag();
  if ( xt.IsStrangeEvent() ) {
    evt.add_attribute( "GENIE.Interaction.StrangeHadronPDG",
      std::make_shared< HepMC3::IntAttribute >(xt.StrangeHadronPdg()) );
  }
  if ( xt.IsCharmEvent() ) {
    evt.add_attribute( "GENIE.Interaction.CharmHadronPDG",
      std::make_shared< HepMC3::IntAttribute >(xt.CharmHadronPdg()) );
  }
  if ( xt.IsFinalLeptonEvent() ) {
    evt.add_attribute( "GENIE.Interaction.FinalLeptonPDG",
      std::make_shared< HepMC3::IntAttribute >(xt.FinalLeptonPdg()) );
  }
  if ( xt.IsFinalQuarkEvent() ) {
    evt.add_attribute( "GENIE.Interaction.FinalQuarkPDG",
      std::make_shared< HepMC3::IntAttribute >(xt.FinalQuarkPdg()) );
  }
  if ( xt.KnownResonance() ) {
    int res_id = static_cast< int >( xt.Resonance() );
    evt.add_attribute( "GENIE.Interaction.Resonance",
      std::make_shared< HepMC3::IntAttribute >(res_id) );
  }

  int decay_mode = xt.DecayMode();
  if ( decay_mode != DEFAULT_DECAY_MODE ) {
    evt.add_attribute( "GENIE.Interaction.DecayMode",
      std::make_shared< HepMC3::IntAttribute >(decay_mode) );
  }

  int num_protons = xt.NProtons();
  if ( num_protons != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NProtons",
      std::make_shared< HepMC3::IntAttribute >(num_protons) );
  }

  int num_neutrons = xt.NNeutrons();
  if ( num_neutrons != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NNeutrons",
      std::make_shared< HepMC3::IntAttribute >(num_neutrons) );
  }

  int num_pi0 = xt.NPi0();
  if ( num_pi0 != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NPi0",
      std::make_shared< HepMC3::IntAttribute >(num_pi0) );
  }

  int num_pips = xt.NPiPlus();
  if ( num_pips != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NPiPlus",
      std::make_shared< HepMC3::IntAttribute >(num_pips) );
  }

  int num_pims = xt.NPiMinus();
  if ( num_pims != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NPiMinus",
      std::make_shared< HepMC3::IntAttribute >(num_pims) );
  }

  int num_sgam = xt.NSingleGammas();
  if ( num_sgam != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NSingleGammas",
      std::make_shared< HepMC3::IntAttribute >(num_sgam) );
  }

  int num_rho0 = xt.NRho0();
  if ( num_rho0 != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NRho0",
      std::make_shared< HepMC3::IntAttribute >(num_rho0) );
  }

  int num_rhop = xt.NRhoPlus();
  if ( num_rhop != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NRhoPlus",
      std::make_shared< HepMC3::IntAttribute >(num_rhop) );
  }

  int num_rhom = xt.NRhoMinus();
  if ( num_rhom != 0 ) {
    evt.add_attribute( "GENIE.Interaction.NRhoMinus",
      std::make_shared< HepMC3::IntAttribute >(num_rhom) );
  }

}
//____________________________________________________________________________
genie::Interaction* genie::HepMC3Converter::RetrieveInteraction(
  const HepMC3::GenEvent& evt )
{
  // Start by loading the initial-state information
  int probe_pdg = 0;
  int tgt_pdg = 0;

  auto tgt_pdg_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.TargetPDG" );
  if ( tgt_pdg_ptr ) tgt_pdg = tgt_pdg_ptr->value();

  auto probe_pdg_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.ProbePDG" );
  if ( probe_pdg_ptr ) probe_pdg = probe_pdg_ptr->value();

  genie::InitialState istate( tgt_pdg, probe_pdg );
  genie::Target& tgt = *istate.TgtPtr();

  auto probe_p4_ptr = evt.attribute< HepMC3::VectorDoubleAttribute >(
   "GENIE.Interaction.ProbeP4" );
  if ( probe_p4_ptr ) {
    TLorentzVector temp_p4 = attribute_to_four_vector( *probe_p4_ptr );
    istate.SetProbeP4( temp_p4 );
  }

  auto tgt_p4_ptr = evt.attribute< HepMC3::VectorDoubleAttribute >(
    "GENIE.Interaction.TargetP4" );
  if ( tgt_p4_ptr ) {
    TLorentzVector temp_p4 = attribute_to_four_vector( *tgt_p4_ptr );
    istate.SetTgtP4( temp_p4 );
  }

  auto hit_nuc_pdg_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.HitNucleonPDG" );

  if ( hit_nuc_pdg_ptr ) {

    tgt.SetHitNucPdg( hit_nuc_pdg_ptr->value() );

    auto hit_nuc_p4_ptr = evt.attribute< HepMC3::VectorDoubleAttribute >(
     "GENIE.Interaction.HitNucleonP4" );
    if ( hit_nuc_p4_ptr ) {
      TLorentzVector temp_p4 = attribute_to_four_vector( *hit_nuc_p4_ptr );
      tgt.SetHitNucP4( temp_p4 );
    }

    auto hit_nuc_radius_ptr = evt.attribute< HepMC3::DoubleAttribute >(
      "GENIE.Interaction.HitNucleonRadius" );
    if ( hit_nuc_radius_ptr ) {
      tgt.SetHitNucPosition( hit_nuc_radius_ptr->value() );
    }

  }

  auto hit_q_pdg_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.HitQuarkPDG" );
  if ( hit_q_pdg_ptr ) {
    tgt.SetHitQrkPdg( hit_q_pdg_ptr->value() );

    auto sea_q_ptr = evt.attribute< HepMC3::IntAttribute >(
      "GENIE.Interaction.HitSeaQuark" );
    if ( sea_q_ptr ) tgt.SetHitSeaQrk( sea_q_ptr->value() );
  }

  // TODO: add error handling for when only one of these is set
  auto s_type_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.ScatteringType" );
  auto i_type_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.InteractionType" );

  genie::ProcessInfo proc_info;
  if ( s_type_ptr && i_type_ptr ) {
    auto s_type = static_cast< genie::ScatteringType_t >(
      s_type_ptr->value() );
    auto i_type = static_cast< genie::InteractionType_t >(
      i_type_ptr->value() );

    proc_info.Set( s_type, i_type );
  }

  auto inter = new genie::Interaction( istate, proc_info );

  genie::Kinematics& kine = *inter->KinePtr();

  auto fsl_p4_ptr = evt.attribute< HepMC3::VectorDoubleAttribute >(
   "GENIE.Interaction.FSLeptonP4" );
  if ( fsl_p4_ptr ) {
    TLorentzVector temp_p4 = attribute_to_four_vector( *fsl_p4_ptr );
    kine.SetFSLeptonP4( temp_p4 );
  }

  auto hs_p4_ptr = evt.attribute< HepMC3::VectorDoubleAttribute >(
   "GENIE.Interaction.HadSystP4" );
  if ( hs_p4_ptr ) {
    TLorentzVector temp_p4 = attribute_to_four_vector( *hs_p4_ptr );
    kine.SetHadSystP4( temp_p4 );
  }

  // TODO: add error handling for when only one of these is set
  auto kv_label_ptr = evt.attribute< HepMC3::VectorIntAttribute >(
    "GENIE.Interaction.KineVarLabels" );
  auto kv_value_ptr = evt.attribute< HepMC3::VectorDoubleAttribute >(
    "GENIE.Interaction.KineVarValues" );
  if ( kv_label_ptr && kv_value_ptr ) {
    const auto& label_vec = kv_label_ptr->value();
    const auto& value_vec = kv_value_ptr->value();

    size_t num_KVs = label_vec.size();
    assert( num_KVs == value_vec.size() );
    for ( size_t kv = 0u; kv < num_KVs; ++kv ) {
      auto label = static_cast< genie::KineVar_t >( label_vec.at(kv) );
      kine.SetKV( label, value_vec.at(kv) );
    }
  }

  genie::XclsTag& xt = *inter->ExclTagPtr();

  auto strange_pdg_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.StrangeHadronPDG" );
  if ( strange_pdg_ptr ) {
    xt.SetStrange( strange_pdg_ptr->value() );
  }

  auto charm_pdg_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.CharmHadronPDG" );
  if ( charm_pdg_ptr ) {
    xt.SetCharm( charm_pdg_ptr->value() );
  }

  auto fslep_pdg_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.FinalLeptonPDG" );
  if ( fslep_pdg_ptr ) {
    xt.SetFinalLepton( fslep_pdg_ptr->value() );
  }

  auto fsq_pdg_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.FinalQuarkPDG" );
  if ( fsq_pdg_ptr ) {
    xt.SetFinalQuark( fsq_pdg_ptr->value() );
  }

  auto res_id_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.Resonance" );
  if ( res_id_ptr ) {
    auto res_id = static_cast< genie::Resonance_t >( res_id_ptr->value() );
    xt.SetResonance( res_id );
  }

  auto dm_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.DecayMode" );
  if ( dm_ptr ) {
    xt.SetDecayMode( dm_ptr->value() );
  }

  auto num_p_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NProtons" );
  if ( num_p_ptr ) {
    xt.SetNProtons( num_p_ptr->value() );
  }

  auto num_n_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NNeutrons" );
  if ( num_n_ptr ) {
    xt.SetNNeutrons( num_n_ptr->value() );
  }

  int num_pi0 = 0;
  int num_pip = 0;
  int num_pim = 0;

  auto num_pi0_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NPi0" );
  if ( num_pi0_ptr ) num_pi0 = num_pi0_ptr->value();

  auto num_pim_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NPiMinus" );
  if ( num_pim_ptr ) num_pim = num_pim_ptr->value();

  auto num_pip_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NPiPlus" );
  if ( num_pip_ptr ) num_pip = num_pip_ptr->value();

  xt.SetNPions( num_pip, num_pi0, num_pim );

  auto num_sgam_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NSingleGammas" );
  if ( num_sgam_ptr ) {
    xt.SetNSingleGammas( num_sgam_ptr->value() );
  }

  int num_rho0 = 0;
  int num_rhop = 0;
  int num_rhom = 0;

  auto num_rho0_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NRho0" );
  if ( num_rho0_ptr ) num_rho0 = num_rho0_ptr->value();

  auto num_rhop_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NRhoPlus" );
  if ( num_rhop_ptr ) num_rhop = num_rhop_ptr->value();

  auto num_rhom_ptr = evt.attribute< HepMC3::IntAttribute >(
    "GENIE.Interaction.NRhoMinus" );
  if ( num_rhom_ptr ) num_rhom = num_rhom_ptr->value();

  xt.SetNRhos( num_rhop, num_rho0, num_rhom );

  return inter;
}
#endif //__GENIE_HEPMC3_INTERFACE_ENABLED__
