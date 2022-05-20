//----------------------------------------------------------------------------
/*!

  Implementation of NHLFluxCreator

 */
//----------------------------------------------------------------------------

#include "Physics/NeutralHeavyLepton/NHLFluxCreator.h"

using namespace genie;
using namespace genie::NHL;

// extern definitions
std::string NHLFluxCreator::fCurrPath;

std::map< NHLProd_t, double > NHLFluxCreator::dynamicScores;
std::map< NHLProd_t, double > NHLFluxCreator::dynamicScores_pion;
std::map< NHLProd_t, double > NHLFluxCreator::dynamicScores_kaon;
std::map< NHLProd_t, double > NHLFluxCreator::dynamicScores_muon;
std::map< NHLProd_t, double > NHLFluxCreator::dynamicScores_neuk;
double NHLFluxCreator::BR_pi2mu, NHLFluxCreator::BR_pi2e;
double NHLFluxCreator::BR_K2mu,  NHLFluxCreator::BR_K2e;
double NHLFluxCreator::BR_K3mu,  NHLFluxCreator::BR_K3e;
double NHLFluxCreator::BR_K03mu, NHLFluxCreator::BR_K03e;

bool NHLFluxCreator::doProduceHists = false;
bool NHLFluxCreator::isParentOnAxis = true;

double NHLFluxCreator::fLx, NHLFluxCreator::fLy, NHLFluxCreator::fLz;
double NHLFluxCreator::fCx, NHLFluxCreator::fCy, NHLFluxCreator::fCz;
double NHLFluxCreator::fAx1, NHLFluxCreator::fAz, NHLFluxCreator::fAx2;
double NHLFluxCreator::fDx, NHLFluxCreator::fDy, NHLFluxCreator::fDz;

double NHLFluxCreator::parentMass, NHLFluxCreator::parentMomentum, NHLFluxCreator::parentEnergy;

TGenPhaseSpace NHLFluxCreator::fPhaseSpaceGenerator;

TF1 * NHLFluxCreator::fNHL, * NHLFluxCreator::fSMv;

double NHLFluxCreator::potnum;
int    NHLFluxCreator::decay_ptype;
double NHLFluxCreator::decay_vx, NHLFluxCreator::decay_vy, NHLFluxCreator::decay_vz;
double NHLFluxCreator::decay_pdpx, NHLFluxCreator::decay_pdpy, NHLFluxCreator::decay_pdpz;
double NHLFluxCreator::decay_necm;
double NHLFluxCreator::decay_nimpwt;

int    NHLFluxCreator::job;
double NHLFluxCreator::pots;

TFile * NHLFluxCreator::fin = 0;
TTree * NHLFluxCreator::tree = 0, * NHLFluxCreator::meta = 0;
bool NHLFluxCreator::isTreeInit = false, NHLFluxCreator::isMetaInit = false, NHLFluxCreator::isBoxInit = false;

//----------------------------------------------------------------------------
void NHLFluxCreator::MakeTupleFluxEntry( int iEntry, flux::GNuMIFluxPassThroughInfo * gnmf, std::string finpath, int run )
{
  // This method creates 1 NHL from the flux info and saves the information
  // Essentially, it replaces a SMv with an NHL

  // Open flux input and initialise trees
  if( !isTreeInit || !isMetaInit || !isBoxInit ){
    OpenFluxInput( finpath );
    InitialiseTree();
    InitialiseMeta();
    MakeBBox();
  }

  // All these in m
  // Beam (0,0,0) is user (0.2486, 60.35,   -1022.74) // = -(fCx, fCy, fCz)
  // Beam (0,1,0) is user (0.2486, 61.3483, -1022.68)
  // Beam (0,0,1) is user (0.2486, 60.2917, -1021.74)

  TVector3 fCvec_beam( fCx, fCy, fCz );
  TVector3 fCvec = ApplyUserRotation( fCvec_beam );

  tree->GetEntry(iEntry);
    
  // turn cm to m and make origin wrt detector 
  fDx = decay_vx * units::cm / units::m;
  fDy = decay_vy * units::cm / units::m;
  fDz = decay_vz * units::cm / units::m;
  TVector3 fDvec( fDx, fDy, fDz );
  TVector3 fDvec_beam = ApplyUserRotation( fDvec, true );

  TVector3 detO_beam( fCvec_beam.X() - fDvec_beam.X(),
		      fCvec_beam.Y() - fDvec_beam.Y(),
		      fCvec_beam.Z() - fDvec_beam.Z() ); // separation in beam coords
  TVector3 detO( fCvec.X() - fDvec.X(),
		 fCvec.Y() - fDvec.Y(),
		 fCvec.Z() - fDvec.Z() ); // separation in rotated coords
  TVector3 detO_user( -detO.X(), -detO.Y(), -detO.Z() );

  LOG( "NHL", pDEBUG )
    << "\n\n\t***** In BEAM coords: *****"
    << "\nCentre     = " << utils::print::Vec3AsString( &fCvec_beam )
    << "\nDecay      = " << utils::print::Vec3AsString( &fDvec_beam )
    << "\nSeparation = " << utils::print::Vec3AsString( &detO_beam )
    << "\n\t***** In ROTATED coords: *****"
    << "\nCentre     = " << utils::print::Vec3AsString( &fCvec )
    << "\nDecay      = " << utils::print::Vec3AsString( &fDvec )
    << "\nSeparation = " << utils::print::Vec3AsString( &detO )
    << "\n\t***** In USER coords: *****"
    << "\nSeparation = " << utils::print::Vec3AsString( &detO_user )
    << "\n\n";
  
  double acc_saa = CalculateDetectorAcceptanceSAA( detO_user );
  //double acc_drc = CalculateDetectorAcceptanceDRC( detO, fLx, fLy, fLz );
  
  // set parent mass
  switch( std::abs( decay_ptype ) ){
  case kPdgPiP: case kPdgKP: case kPdgMuon: case kPdgK0L:
    parentMass = PDGLibrary::Instance()->Find(decay_ptype)->Mass(); break;
  default:
    LOG( "NHL", pERROR ) << "Parent with PDG code " << decay_ptype << " not handled!"
			 << "\n\tProceeding, but results are possibly unphysical.";
    parentMass = PDGLibrary::Instance()->Find(decay_ptype)->Mass(); break;
  }
  parentMomentum = std::sqrt( decay_pdpx*decay_pdpx + decay_pdpy*decay_pdpy + decay_pdpz*decay_pdpz );
  parentEnergy = std::sqrt( parentMass*parentMass + parentMomentum*parentMomentum );

  isParentOnAxis = utils::nhl::GetCfgBool( "NHL", "FluxCalc", "IsParentOnAxis" );
  // CAUTION: p4par is in USER coords
  TLorentzVector p4par = ( isParentOnAxis ) ? 
    TLorentzVector( parentMomentum * (detO.Unit()).X(), 
		    parentMomentum * (detO.Unit()).Y(),
		    parentMomentum * (detO.Unit()).Z(),
		    parentEnergy ) :
    TLorentzVector( decay_pdpx, decay_pdpy, decay_pdpz, parentEnergy );

  if( !isParentOnAxis ){
    // rotate p4par to user coordinates
    TVector3 tmpv3 = ApplyUserRotation( p4par.Vect() );
    p4par.SetPxPyPzE( tmpv3.Px(), tmpv3.Py(), tmpv3.Pz(), p4par.E() );
  }


  TVector3 boost_beta = p4par.BoostVector();

  // explicitly check if there are any allowed decays for this parent
  bool canGoForward = true;
  switch( std::abs( decay_ptype ) ){
  case kPdgPiP:
    canGoForward = 
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdPion2Muon ) ||
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdPion2Electron ); break;
  case kPdgKP:
    canGoForward =
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdKaon2Muon ) || 
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdKaon2Electron ) ||
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdKaon3Muon ) ||
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdKaon3Electron ); break;
  case kPdgMuon:
    canGoForward =
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdMuon3Nue ); break; // SM nus are massless so doesn't matter
  case kPdgK0L:
    canGoForward =
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdNeuk3Muon ) ||
      utils::nhl::IsProdKinematicallyAllowed( kNHLProdNeuk3Electron ); break;
  }

  if( !canGoForward ){ FillNonsense( iEntry, gnmf, run ); return; }
  // now calculate which decay channel produces the NHL.
  dynamicScores = GetProductionProbs( decay_ptype );
  assert( dynamicScores.size() > 0 );
  
  RandomGen * rnd = RandomGen::Instance();
  double score = rnd->RndGen().Uniform( 0.0, 1.0 );
  NHLProd_t prodChan;
  // compare with cumulative prob. If < 1st in map, pick 1st chan. If >= 1st and < (1st+2nd), pick 2nd, etc
  
  unsigned int imap = 0; double s1 = 0.0;
  std::map< NHLProd_t, double >::iterator pdit = dynamicScores.begin();
  while( score >= s1 && pdit != dynamicScores.end() ){
    s1 += (*pdit).second;
    if( parentMass > 0.495 ){
      LOG( "NHL", pDEBUG )
	<< "(*pdit).first = " << utils::nhl::ProdAsString( (*pdit).first )
	<< " : (*pdit).second = " << (*pdit).second;
    }
    if( score >= s1 ){
      imap++; pdit++;
    }
  }
  assert( imap < dynamicScores.size() ); // should have decayed to *some* NHL
  prodChan = (*pdit).first;

  LOG( "NHL", pDEBUG )
    << "Selected channel: " << utils::nhl::ProdAsString( prodChan );
  
  // decay channel specified, now time to make kinematics
  TLorentzVector p4NHL_rest = NHLEnergy( prodChan, p4par ); // this is a random direction rest-frame NHL. We don't care about where it's pointing

  /* 
   * it is NOT sufficient to boost this into lab frame! 
   * Only a small portion of the CM decays can possibly reach the detector, 
   * imposing a constraint on the allowed directions of p4NHL_rest. 
   * You will miscalculate the NHL energy if you just Boost here. 
   */
  // explicitly calculate the boost correction to lab-frame energy
  // in a dk2nu-like fashion. See bsim::CalcEnuWgt()
  double betaMag = boost_beta.Mag();
  double gamma   = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
  double betaNHL = p4NHL_rest.P() / p4NHL_rest.E();
  double boost_correction = 0.0;
  double costh_pardet = 0.0;
  if( parentMomentum > 0.0 ){
    costh_pardet = ( decay_pdpx * detO.X() +
		     decay_pdpy * detO.Y() +
		     decay_pdpz * detO.Z() ) / ( parentMomentum * detO.Mag() );
    if( costh_pardet < -1.0 ) costh_pardet = -1.0;
    if( costh_pardet > 1.0 ) costh_pardet = 1.0;
    // assume boost is on z' direction where z' = parent momentum direction, subbing betaMag ==> betaMag * costh_pardet
    boost_correction = 1.0 / ( gamma * ( 1.0 - betaMag * betaNHL * costh_pardet ) );
  }

  LOG( "NHL", pFATAL )
    << "\n *** Event " << iEntry << " *** :"
    << "\nbetaMag = " << betaMag
    << "\ngamma   = " << gamma
    << "\nbetaNHL = " << betaNHL
    << "\ncosth_pardet     = " << costh_pardet
    << "\nboost_correction = " << boost_correction
    << "\nChannel = " << utils::nhl::ProdAsString( prodChan )
    << "\nCM ENHL = " << p4NHL_rest.E()
    << "\nENHL = " << p4NHL_rest.E() * boost_correction
    << "\nacceptance = " << acc_saa * boost_correction * boost_correction;

  assert( boost_correction > 0.0 );

  // so now we have the random decay. Direction = parent direction, energy = what we calculated
  double ENHL = p4NHL_rest.E() * boost_correction;
  TVector3 pdu = ( p4par.Vect() ).Unit();
  TLorentzVector p4NHL_rand( ENHL * pdu.X(), ENHL * pdu.Y(), ENHL * pdu.Z(), ENHL );

  // find random point in BBox and force momentum to point to that point
  // first, separation in beam frame
  TVector3 fRVec_beam = PointToRandomPointInBBox( detO_beam );
  // rotate it and get unit
  TVector3 fRVec_unit = (ApplyUserRotation( fRVec_beam )).Unit();
  // force NHL to point along this direction
  TLorentzVector p4NHL( p4NHL_rand.P() * fRVec_unit.X(),
			p4NHL_rand.P() * fRVec_unit.Y(),
			p4NHL_rand.P() * fRVec_unit.Z(),
			p4NHL_rand.E() );

  TVector3 pNHL_beam = ApplyUserRotation( p4NHL.Vect(), true );
  TLorentzVector p4NHL_beam( pNHL_beam.X(), pNHL_beam.Y(), pNHL_beam.Z(), p4NHL.E() );
  
  LOG( "NHL", pDEBUG )
    << "\nRandom:  " << utils::print::P4AsString( &p4NHL_rand )
    << "\nPointed: " << utils::print::P4AsString( &p4NHL )
    << "\nRest:    " << utils::print::P4AsString( &p4NHL_rest );
  
  // calculate acceptance correction
  // first, get minimum and maximum deviation from parent momentum to hit detector in degrees
  // RETHERE generalise condition in case momentum hits detector

  double zm = ( isParentOnAxis ) ? 0.0 : GetAngDeviation( p4par, detO, false );
  double zp = GetAngDeviation( p4par, detO, true );
  
  // now get the actual acceptance correction
  LOG( "NHL", pDEBUG )
    << "Entering CalculateAcceptanceCorrection with:"
    << "\np4par      = " << utils::print::P4AsString(&p4par)
    << "\np4NHL_rest = " << utils::print::P4AsString(&p4NHL_rest)
    << "\ndecay_necm = " << decay_necm
    << "\nzm = " << zm
    << "\nzp = " << zp;
  double accCorr = CalculateAcceptanceCorrection( p4par, p4NHL_rest, decay_necm, zm, zp );
  
  // also have to factor in boost correction itself... that's same as energy boost correction squared
  // which means a true acceptance of...
  double acceptance = acc_saa * boost_correction * boost_correction * accCorr;

  LOG( "NHL", pDEBUG )
    << "\nacc_saa    = " << acc_saa
    << "\nboost_correction^2 = " << boost_correction * boost_correction
    << "\naccCorr    = " << accCorr
    << "\nacceptance = " << acceptance;

  // RETHERE add delay! 
  // write 4-position all this happens at
  TLorentzVector x4NHL_beam( decay_vx, decay_vy, decay_vz, 0.0 ); // in cm
  TLorentzVector x4NHL( -detO.X(), -detO.Y(), -detO.Z(), 0.0 ); // in m
  TLorentzVector x4NHL_cm( units::m / units::cm * ( -detO.X() ),
			   units::m / units::cm * ( -detO.Y() ),
			   units::m / units::cm * ( -detO.Z() ), 0.0 ); // in cm

  LOG( "NHL", pDEBUG )
    << "\nx4NHL_beam    = " << utils::print::X4AsString( &x4NHL_beam )
    << "\nx4NHL_user    = " << utils::print::X4AsString( &x4NHL )
    << "\nx4NHL_user_cm = " << utils::print::X4AsString( &x4NHL_cm );

  // fill all the GNuMIFlux stuff
  // comments as seeon on https://www.hep.utexas.edu/~zarko/wwwgnumi/v19/v19/output_gnumi.html

  gnmf->pcodes = 1;                          ///< converted to PDG
  gnmf->units = 0;                           ///< cm
  
  int typeMod = ( decay_ptype > 0 ) ? 1 : -1;
  gnmf->fgPdgC = typeMod * kPdgNHL;          ///< PDG code

  gnmf->fgXYWgt = acceptance;                ///< geometrical * collimation correction

  gnmf->fgP4 = p4NHL_beam;                   ///< generated 4-momentum, beam coord [GeV]
  gnmf->fgX4 = x4NHL_beam;                   ///< generated 4-position, beam coord [cm]
  gnmf->fgP4User = p4NHL;                    ///< generated 4-momentum, user coord [GeV]
  gnmf->fgX4User = x4NHL_cm;                  ///< generated 4-position, user coord [cm]

  gnmf->run      = run;                      ///< Run number
  gnmf->evtno    = iEntry;                   ///< Event number (proton on target) 
                                                 // RETHERE which is it?
  gnmf->ndxdz    = p4NHL.Px() / p4NHL.Pz();  ///< Neutrino direction slope for a random decay
  gnmf->ndydz    = p4NHL.Py() / p4NHL.Pz();  ///< See above
  gnmf->npz      = p4NHL.Pz();               ///< Neutrino momentum [GeV] along z direction (beam axis)
  gnmf->nenergy  = p4NHL.E();                ///< Neutrino energy [GeV] for a random decay
  gnmf->ndxdznea = -9999.9;                  ///< Neutrino direction slope for a decay forced to ND
  gnmf->ndydznea = -9999.9;                  ///< See above
  gnmf->nenergyn = -9999.9;                  ///< Neutrino energy for decay forced to ND
  gnmf->nwtnear  = accCorr;                  ///< weight for decay forced to ND / now acceptance correction
  gnmf->ndxdzfar = -9999.9;                  ///< Same as ND but FD
  gnmf->ndydzfar = -9999.9;                  ///< See above
  gnmf->nenergyf = -9999.9;                  ///< See above
  gnmf->nwtfar   = -9999.9;                  ///< See above
  gnmf->norig    = -9999.9;                  ///< Obsolete...
  
  int iNdecay = -1, iNtype = -1;
  switch( prodChan ){
  case kNHLProdNeuk3Electron: iNdecay = 1; iNtype = 53; break;
  case kNHLProdNeuk3Muon: iNdecay = 3; iNtype = 56; break;
  case kNHLProdKaon2Electron: iNdecay = 5; iNtype = 53; break;
  case kNHLProdKaon2Muon: iNdecay = 7; iNtype = 56; break;
  case kNHLProdKaon3Electron: iNdecay = 9; iNtype = 53; break;
  case kNHLProdKaon3Muon: iNdecay = 11; iNtype = 56; break;
  case kNHLProdMuon3Nue: iNdecay = 13; iNtype = 53; break;
  case kNHLProdMuon3Numu: iNdecay = 15; iNtype = 53; break;
  case kNHLProdMuon3Nutau: iNdecay = 17; iNtype = 53; break;
  case kNHLProdPion2Electron: iNdecay = 19; iNtype = 53; break;
  case kNHLProdPion2Muon: iNdecay = 21; iNtype = 56; break;
  default: iNdecay = -2; iNtype = -2; break;
  }
  if( decay_ptype < 0 ){ iNdecay++; iNtype--; }
  if( iNdecay > 0 ) iNdecay += 30; // horrendous hack to fit all the 22 decay channels...
  if( iNtype > 0 ) iNtype += 10; // and a hack for "NHL type" (referring to co-produced leptons)

  gnmf->ndecay = iNdecay;                    ///< Decay mode that produced neutrino
  gnmf->ntype  = iNtype;                     ///< Neutrino "flavour" (i.e. of co-produced lepton)

  gnmf->vx = decay_vx;                       ///< X position of hadron/muon decay
  gnmf->vy = decay_vy;                       ///< Y position of hadron/muon decay
  gnmf->vz = decay_vz;                       ///< Z position of hadron/muon decay

  gnmf->pdpx = decay_pdpx;                   ///< Parent X momentum at decay point
  gnmf->pdpy = decay_pdpy;                   ///< Parent Y momentum at decay point
  gnmf->pdpz = decay_pdpz;                   ///< Parent Z momentum at decay point

  gnmf->ppdxdz = -9999.9;                    ///< Parent dxdz direction at production
  gnmf->ppdydz = -9999.9;                    ///< Parent dydz direction at production
  gnmf->pppz = -9999.9;                      ///< Parent energy at production

  gnmf->ppmedium = -9999;                    ///< Tracking medium number where parent was produced
  gnmf->ptype = decay_ptype;                 ///< Parent GEANT code particle ID converted to PDG

  gnmf->ppvx = -9999.9;                      ///< Parent production vertex X (cm)
  gnmf->ppvy = -9999.9;                      ///< Parent production vertex Y (cm)
  gnmf->ppvz = -9999.9;                      ///< Parent production vertex Z (cm)

  gnmf->necm = p4NHL_rest.E();               ///< Neutrino energy in COM frame
  gnmf->nimpwt = decay_nimpwt;               ///< Weight of neutrino parent

  gnmf->xpoint = p4par.Px();                 ///< Used here to store parent px in user coords
  gnmf->ypoint = p4par.Py();                 ///< Used here to store parent py in user coords
  gnmf->zpoint = p4par.Pz();                 ///< Used here to store parent pz in user coords

  gnmf->tvx = -9999.9;                       ///< X exit point of parent particle at the target
  gnmf->tvy = -9999.9;                       ///< Y exit point of parent particle at the target
  gnmf->tvz = -9999.9;                       ///< Z exit point of parent particle at the target

  gnmf->tpx = -9999.9;                       ///< Parent momentum exiting the target (X)
  gnmf->tpy = -9999.9;                       ///< Parent momentum exiting the target (Y)
  gnmf->tpz = -9999.9;                       ///< Parent momentum exiting the target (Z)

  gnmf->tptype = -9999;                      ///< Parent particle ID exiting the target conv to PDG
  gnmf->tgen = -9999;                        ///< Parent generation in cascade

  gnmf->tgptype = -9999;                     ///< Type of particle that created a particle...
  
  gnmf->tgppx = -9999.9;                     ///< Momentum of particle that created particle at IP
  gnmf->tgppy = -9999.9;                     ///< Momentum of particle that created particle at IP
  gnmf->tgppz = -9999.9;                     ///< Momentum of particle that created particle at IP

  gnmf->tprivx = -9999.9;                    ///< Primary particle interaction vertex
  gnmf->tprivy = -9999.9;                    ///< Primary particle interaction vertex
  gnmf->tprivz = -9999.9;                    ///< Primary particle interaction vertex

  gnmf->beamx = -9999.9;                     ///< Primary proton origin
  gnmf->beamy = -9999.9;                     ///< Primary proton origin
  gnmf->beamz = -9999.9;                     ///< Primary proton origin

  gnmf->beampx = -9999.9;                    ///< Primary proton momentum
  gnmf->beampy = -9999.9;                    ///< Primary proton momentum
  gnmf->beampz = -9999.9;                    ///< Primary proton momentum

#ifndef SKIP_MINERVA_MODS
  gnmf->ntrajectory = -9;
  gnmf->overflow = false;

  for( unsigned int i = 0; i < 10; i++ ){
    gnmf->pdgcode[i] = -9;
    gnmf->trackId[i] = -9;
    gnmf->parentId[i] = -9;
    
    gnmf->startx[i] = -9999.9;
    gnmf->starty[i] = -9999.9;
    gnmf->startz[i] = -9999.9;
    gnmf->startpx[i] = -9999.9;
    gnmf->startpy[i] = -9999.9;
    gnmf->startpz[i] = -9999.9;
    gnmf->stopx[i] = -9999.9;
    gnmf->stopy[i] = -9999.9;
    gnmf->stopz[i] = -9999.9;
    gnmf->pprodpx[i] = -9999.9;
    gnmf->pprodpy[i] = -9999.9;
    gnmf->pprodpz[i] = -9999.9;

    gnmf->proc[i] = -9;
    gnmf->ivol[i] = -9;
    gnmf->fvol[i] = -9;
  }
#endif

  LOG( "NHL", pDEBUG )
    << "Finished MakeTupleFluxEntry(). Returning to App";
  
}
//----------------------------------------------------------------------------
void NHLFluxCreator::FillNonsense( int iEntry, flux::GNuMIFluxPassThroughInfo * gnmf, int run )
{
  gnmf->pcodes = 1;                          ///< converted to PDG
  gnmf->units = 0;                           ///< cm
  
  gnmf->fgPdgC = -9999;                      ///< PDG code

  gnmf->fgXYWgt = -9999.9;                   ///< geometrical * collimation correction

  TLorentzVector dv( -9999.9, -9999.9, -9999.9, -9999.9 );
  gnmf->fgP4 = dv;                           ///< generated 4-momentum, beam coord
  gnmf->fgX4 = dv;                           ///< generated 4-position, beam coord
  gnmf->fgP4User = dv;                       ///< generated 4-momentum, user coord
  gnmf->fgX4User = dv;                       ///< generated 4-position, user coord

  gnmf->run      = run;                      ///< Run number
  gnmf->evtno    = iEntry;                   ///< Event number (proton on target) 
                                                 // RETHERE which is it?
  gnmf->ndxdz    = -9999.9;                  ///< Neutrino direction slope for a random decay
  gnmf->ndydz    = -9999.9;                  ///< See above
  gnmf->npz      = -9999.9;                  ///< Neutrino momentum [GeV] along z direction (beam axis)
  gnmf->nenergy  = -9999.9;                  ///< Neutrino energy [GeV] for a random decay
  gnmf->ndxdznea = -9999.9;                  ///< Neutrino direction slope for a decay forced to ND
  gnmf->ndydznea = -9999.9;                  ///< See above
  gnmf->nenergyn = -9999.9;                  ///< Neutrino energy for decay forced to ND
  gnmf->nwtnear  = -9999.9;                  ///< weight for decay forced to ND
  gnmf->ndxdzfar = -9999.9;                  ///< Same as ND but FD
  gnmf->ndydzfar = -9999.9;                  ///< See above
  gnmf->nenergyf = -9999.9;                  ///< See above
  gnmf->nwtfar   = -9999.9;                  ///< See above
  gnmf->norig    = -9999;                    ///< Obsolete...
  
  int iNdecay = -1, iNtype = -1;
  gnmf->ndecay = iNdecay;                    ///< Decay mode that produced neutrino
  gnmf->ntype  = iNtype;                     ///< Neutrino "flavour" (i.e. of co-produced lepton)

  gnmf->vx = decay_vx;                       ///< X position of hadron/muon decay
  gnmf->vy = decay_vy;                       ///< Y position of hadron/muon decay
  gnmf->vz = decay_vz;                       ///< Z position of hadron/muon decay

  gnmf->pdpx = decay_pdpx;                   ///< Parent X momentum at decay point
  gnmf->pdpy = decay_pdpy;                   ///< Parent Y momentum at decay point
  gnmf->pdpz = decay_pdpz;                   ///< Parent Z momentum at decay point

  gnmf->ppdxdz = -9999.9;                    ///< Parent dxdz direction at production
  gnmf->ppdydz = -9999.9;                    ///< Parent dydz direction at production
  gnmf->pppz = -9999.9;                      ///< Parent energy at production

  gnmf->ppmedium = -9999;                    ///< Tracking medium number where parent was produced
  gnmf->ptype = decay_ptype;                 ///< Parent GEANT code particle ID converted to PDG

  gnmf->ppvx = -9999.9;                      ///< Parent production vertex X (cm)
  gnmf->ppvy = -9999.9;                      ///< Parent production vertex Y (cm)
  gnmf->ppvz = -9999.9;                      ///< Parent production vertex Z (cm)
  
  gnmf->necm = -9999.9;                      ///< Neutrino energy in COM frame
  gnmf->nimpwt = decay_nimpwt;               ///< Weight of neutrino parent

  gnmf->xpoint = -9999.9;                    ///< Debugging hook (unused)
  gnmf->ypoint = -9999.9;                    ///< Debugging hook (unused)
  gnmf->zpoint = -9999.9;                    ///< Debugging hook (unused)

  gnmf->tvx = -9999.9;                       ///< X exit point of parent particle at the target
  gnmf->tvy = -9999.9;                       ///< Y exit point of parent particle at the target
  gnmf->tvz = -9999.9;                       ///< Z exit point of parent particle at the target

  gnmf->tpx = -9999.9;                       ///< Parent momentum exiting the target (X)
  gnmf->tpy = -9999.9;                       ///< Parent momentum exiting the target (Y)
  gnmf->tpz = -9999.9;                       ///< Parent momentum exiting the target (Z)

  gnmf->tptype = -9999;                      ///< Parent particle ID exiting the target conv to PDG
  gnmf->tgen = -9999;                        ///< Parent generation in cascade

  gnmf->tgptype = -9999;                     ///< Type of particle that created a particle...
  
  gnmf->tgppx = -9999.9;                     ///< Momentum of particle that created particle at IP
  gnmf->tgppy = -9999.9;                     ///< Momentum of particle that created particle at IP
  gnmf->tgppz = -9999.9;                     ///< Momentum of particle that created particle at IP

  gnmf->tprivx = -9999.9;                    ///< Primary particle interaction vertex
  gnmf->tprivy = -9999.9;                    ///< Primary particle interaction vertex
  gnmf->tprivz = -9999.9;                    ///< Primary particle interaction vertex

  gnmf->beamx = -9999.9;                     ///< Primary proton origin
  gnmf->beamy = -9999.9;                     ///< Primary proton origin
  gnmf->beamz = -9999.9;                     ///< Primary proton origin

  gnmf->beampx = -9999.9;                    ///< Primary proton momentum
  gnmf->beampy = -9999.9;                    ///< Primary proton momentum
  gnmf->beampz = -9999.9;                    ///< Primary proton momentum

#ifndef SKIP_MINERVA_MODS
  gnmf->ntrajectory = -9;
  gnmf->overflow = false;

  for( unsigned int i = 0; i < 10; i++ ){
    gnmf->pdgcode[i] = -9;
    gnmf->trackId[i] = -9;
    gnmf->parentId[i] = -9;
    
    gnmf->startx[i] = -9999.9;
    gnmf->starty[i] = -9999.9;
    gnmf->startz[i] = -9999.9;
    gnmf->startpx[i] = -9999.9;
    gnmf->startpy[i] = -9999.9;
    gnmf->startpz[i] = -9999.9;
    gnmf->stopx[i] = -9999.9;
    gnmf->stopy[i] = -9999.9;
    gnmf->stopz[i] = -9999.9;
    gnmf->pprodpx[i] = -9999.9;
    gnmf->pprodpy[i] = -9999.9;
    gnmf->pprodpz[i] = -9999.9;

    gnmf->proc[i] = -9;
    gnmf->ivol[i] = -9;
    gnmf->fvol[i] = -9;
  }
#endif
}
//----------------------------------------------------------------------------
void NHLFluxCreator::OpenFluxInput( std::string finpath )
{
  if( std::strcmp( finpath.c_str(), fCurrPath.c_str() ) == 0 ) return;

  fCurrPath = finpath;

  LOG( "NHL", pDEBUG )
    << "Getting flux input from finpath = " << finpath.c_str();

  fin = TFile::Open( finpath.c_str() );
  assert( fin );

  // show some stuff from dkmeta and dk2nu as proof
  tree = dynamic_cast<TTree *>( fin->Get( "dkRootTree" ) );
  meta = dynamic_cast<TTree *>( fin->Get( "dkRootMeta" ) );

  if( !tree ){ LOG( "NHL", pFATAL ) << "Could not open flux tree!"; }
  if( !meta ){ LOG( "NHL", pFATAL ) << "Could not open meta tree!"; }
  assert( tree && meta );

  const int nFilesInMeta = meta->GetEntries();
  int nEntries = tree->GetEntries();

  LOG( "NHL", pDEBUG )
    << "There were " << nFilesInMeta << " files in meta with " << nEntries << " total nus";
}
//----------------------------------------------------------------------------
void NHLFluxCreator::InitialiseTree()
{
  if( isTreeInit ) return;

  potnum = 0.0;
  decay_ptype = 0;
  decay_vx = 0.0; decay_vy = 0.0; decay_vz = 0.0;
  decay_pdpx = 0.0; decay_pdpy = 0.0; decay_pdpz = 0.0;
  decay_nimpwt = 0.0;
  
  tree->SetBranchAddress( "potnum",       &potnum       );
  tree->SetBranchAddress( "decay_ptype",  &decay_ptype  );
  tree->SetBranchAddress( "decay_vx",     &decay_vx     );
  tree->SetBranchAddress( "decay_vy",     &decay_vy     );
  tree->SetBranchAddress( "decay_vz",     &decay_vz     );
  tree->SetBranchAddress( "decay_pdpx",   &decay_pdpx   );
  tree->SetBranchAddress( "decay_pdpy",   &decay_pdpy   );
  tree->SetBranchAddress( "decay_pdpz",   &decay_pdpz   );
  tree->SetBranchAddress( "decay_necm",   &decay_necm   );
  tree->SetBranchAddress( "decay_nimpwt", &decay_nimpwt );
}
//----------------------------------------------------------------------------
void NHLFluxCreator::InitialiseMeta()
{
  if( isMetaInit ) return;
  
  job = 0;
  pots = 0.0;

  meta->SetBranchAddress( "job",  &job  );
  meta->SetBranchAddress( "pots", &pots );
}
//----------------------------------------------------------------------------
void NHLFluxCreator::ReadBRs()
{
  TParticlePDG * pionParticle = PDGLibrary::Instance()->Find( kPdgPiP );
  TParticlePDG * kaonParticle = PDGLibrary::Instance()->Find( kPdgKP  );
  TParticlePDG * neukParticle = PDGLibrary::Instance()->Find( kPdgK0L );

  TObjArray * pionDecayList = pionParticle->DecayList();
  TObjArray * kaonDecayList = kaonParticle->DecayList();
  TObjArray * neukDecayList = neukParticle->DecayList();

  TDecayChannel * pion2muChannel = ( TDecayChannel * ) pionDecayList->At(0);
  TDecayChannel * pion2elChannel = ( TDecayChannel * ) pionDecayList->At(1);

  TDecayChannel * kaon2muChannel = ( TDecayChannel * ) kaonDecayList->At(0);
  //TDecayChannel * kaon2elChannel = 0; // tiny BR, not in genie_pdg_table.txt
  TDecayChannel * kaon3muChannel = ( TDecayChannel * ) kaonDecayList->At(5);
  TDecayChannel * kaon3elChannel = ( TDecayChannel * ) kaonDecayList->At(4);

  TDecayChannel * neuk3muChannel = ( TDecayChannel * ) neukDecayList->At(4);
  TDecayChannel * neuk3elChannel = ( TDecayChannel * ) neukDecayList->At(2);

  BR_pi2mu = pion2muChannel->BranchingRatio();
  BR_pi2e  = pion2elChannel->BranchingRatio();
  
  BR_K2mu  = kaon2muChannel->BranchingRatio();
  BR_K2e   = 1.6e-5; // RETHERE - add to pdg table? From PDG 2021
  BR_K3mu  = kaon3muChannel->BranchingRatio();
  BR_K3e   = kaon3elChannel->BranchingRatio();

  BR_K03mu = 2.0 * neuk3muChannel->BranchingRatio(); // one from K0L-->mu+ and one from -->mu-
  BR_K03e  = 2.0 * neuk3elChannel->BranchingRatio();
}
//----------------------------------------------------------------------------
TVector3 NHLFluxCreator::GetBoostBetaVec( TLorentzVector parp4 )
{
  double px = parp4.Px(), py = parp4.Py(), pz = parp4.Pz(), pE = parp4.E();
  double p = parp4.P();
  double beta_mag = p/pE;
  TVector3 bvec( beta_mag * px/p, beta_mag * py/p, beta_mag * pz/p );
  return bvec;
}
//----------------------------------------------------------------------------
std::map< NHLProd_t, double > NHLFluxCreator::GetProductionProbs( int parPDG )
{
  // check if we've calculated scores before
  switch( std::abs( parPDG ) ){
  case kPdgPiP : if( dynamicScores_pion.size() > 0 ) return dynamicScores_pion; break;
  case kPdgKP  : if( dynamicScores_kaon.size() > 0 ) return dynamicScores_kaon; break;
  case kPdgMuon: if( dynamicScores_muon.size() > 0 ) return dynamicScores_muon; break;
  case kPdgK0L : if( dynamicScores_neuk.size() > 0 ) return dynamicScores_neuk; break;
  default: LOG( "NHL", pWARN ) << "Unknown parent. Proceeding, but results may be unphysical"; break;
  }

  std::map< NHLProd_t, double > dynScores;

  // first get branching ratios to SM
  ReadBRs();
  // then get NHL parameter space
  double M    = utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Mass" );
  double Ue42 = utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Ue42" );
  double Um42 = utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Um42" );
  double Ut42 = utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Ut42" );

  // now get parent mass
  //double mP = PDGLibrary::Instance()->Find( std::abs( parPDG ) )->Mass();
  
  // first get pure kinematic part of the BRs
  double KScale[4] = { -1.0, -1.0, -1.0, -1.0 }, mixScale[4] = { -1.0, -1.0, -1.0, -1.0 };
  double totalMix = 0.0;
  switch( std::abs( parPDG ) ){
  case genie::kPdgMuon:
    KScale[0] = NHLSelector::KScale_Global( kNHLProdMuon3Numu, M );
    KScale[1] = NHLSelector::KScale_Global( kNHLProdMuon3Nue, M ); // same, convenience for later
    KScale[2] = NHLSelector::KScale_Global( kNHLProdMuon3Nutau, M ); // same, convenience for later
    mixScale[0] = 1.0 * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = 1.0 * Ue42 * KScale[1]; totalMix += mixScale[1];
    mixScale[2] = 1.0 * Ut42 * KScale[2]; totalMix += mixScale[2];

    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdMuon3Numu,  mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdMuon3Nue,   mixScale[1] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdMuon3Nutau, mixScale[2] / totalMix } ) );

    dynamicScores_muon = dynScores;
    break;
  case genie::kPdgKP:
    KScale[0] = NHLSelector::KScale_Global( kNHLProdKaon2Muon, M );
    KScale[1] = NHLSelector::KScale_Global( kNHLProdKaon2Electron, M );
    KScale[2] = NHLSelector::KScale_Global( kNHLProdKaon3Muon, M );
    KScale[3] = NHLSelector::KScale_Global( kNHLProdKaon3Electron, M );
    mixScale[0] = BR_K2mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_K2e  * Ue42 * KScale[1]; totalMix += mixScale[1];
    mixScale[2] = BR_K3mu * Um42 * KScale[2]; totalMix += mixScale[2];
    mixScale[3] = BR_K3e  * Ue42 * KScale[3]; totalMix += mixScale[3];

    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdKaon2Muon,     mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdKaon2Electron, mixScale[1] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdKaon3Muon,     mixScale[2] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdKaon3Electron, mixScale[3] / totalMix } ) );

    dynamicScores_kaon = dynScores;
    break;
  case genie::kPdgPiP:

    KScale[0] = NHLSelector::KScale_Global( kNHLProdPion2Muon, M );
    KScale[1] = NHLSelector::KScale_Global( kNHLProdPion2Electron, M );
    mixScale[0] = BR_pi2mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_pi2e  * Ue42 * KScale[1]; totalMix += mixScale[1];

    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdPion2Muon,     mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdPion2Electron, mixScale[1] / totalMix } ) );

    dynamicScores_pion = dynScores;
    break;
  case genie::kPdgK0L:

    KScale[0] = NHLSelector::KScale_Global( kNHLProdNeuk3Muon, M );
    KScale[1] = NHLSelector::KScale_Global( kNHLProdNeuk3Electron, M );
    mixScale[0] = BR_K03mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_K03e  * Ue42 * KScale[1]; totalMix += mixScale[1];

    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdNeuk3Muon,     mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdNeuk3Electron, mixScale[1] / totalMix } ) );

    dynamicScores_neuk = dynScores;
    break;
  default:
    LOG( "NHL", pERROR )
      << "Unknown parent particle. Cannot make scales, exiting."; exit(1);
  }

  LOG( "NHL", pDEBUG )
    << "Score map now has " << dynScores.size() << " elements. Returning.";
  return dynScores;

}
//----------------------------------------------------------------------------
TLorentzVector NHLFluxCreator::NHLEnergy( NHLProd_t nhldm, TLorentzVector p4par )
{
  // first boost to parent rest frame
  TLorentzVector p4par_rest = p4par;
  TVector3 boost_beta = p4par.BoostVector();
  p4par_rest.Boost( -boost_beta );

  LOG( "NHL", pDEBUG )
    << "Attempting to decay rest-system p4 = " << utils::print::P4AsString(&p4par_rest)
    << " as " << utils::nhl::ProdAsString( nhldm );
  
  // get PDGCodeList and truncate 1st member
  PDGCodeList fullList  = utils::nhl::ProductionProductList( nhldm );
  bool        allow_duplicate = true;
  PDGCodeList decayList( allow_duplicate );
  double * mass = new double[decayList.size()];
  double   sum  = 0.0;

  for( std::vector<int>::const_iterator pdg_iter = fullList.begin(); pdg_iter != fullList.end(); ++pdg_iter )
    {
      if( pdg_iter != fullList.begin() ){
	int pdgc = *pdg_iter;
	decayList.push_back( pdgc );
      }
    }

  int iv = 0;
  for( std::vector<int>::const_iterator pdg_iter = decayList.begin(); pdg_iter != decayList.end(); ++pdg_iter )
    {
      int pdgc = *pdg_iter;
      double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
      mass[iv++] = m; sum += m;
    }
  
  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay( p4par_rest, decayList.size(), mass );
  if(!permitted) {
    LOG("NHL", pERROR)
      << " *** Phase space decay is not permitted \n"
      << " Total particle mass = " << sum << "\n"
      << " Decaying system p4 = " << utils::print::P4AsString(&p4par_rest);
    // clean-up
    delete [] mass;
    // throw exception
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Decay not permitted kinematically");
    exception.SwitchOnFastForward();
    throw exception;
  }

  // Get the maximum weight
  //double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int idec=0; idec<200; idec++) {
     double w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  wmax *= 2;

  LOG("NHL", pNOTICE)
     << "Max phase space gen. weight @ current NHL system: " << wmax;

  // Generate an unweighted decay
  RandomGen * rnd = RandomGen::Instance();
  
  bool accept_decay=false;
  unsigned int itry=0;
  while(!accept_decay)
    {
      itry++;
      
      if(itry > controls::kMaxUnweightDecayIterations) {
	// report, clean-up and return
	LOG("NHL", pWARN)
	  << "Couldn't generate an unweighted phase space decay after "
	  << itry << " attempts";
	// clean up
	delete [] mass;
	// throw exception
	genie::exceptions::EVGThreadException exception;
	exception.SetReason("Couldn't select decay after N attempts");
	exception.SwitchOnFastForward();
	throw exception;
      }
      double w  = fPhaseSpaceGenerator.Generate();
      if(w > wmax) {
        LOG("NHL", pWARN)
	  << "Decay weight = " << w << " > max decay weight = " << wmax;
      }
      double gw = wmax * rnd->RndHadro().Rndm();
      accept_decay = (gw<=w);
      
      LOG("NHL", pINFO)
        << "Decay weight = " << w << " / R = " << gw
        << " - accepted: " << accept_decay;
      
    } //!accept_decay
  
  // [DON'T] Insert final state products into a TClonesArray of GHepParticle's
  // Grab 0th entry energy and return that
  int idp = 0; TLorentzVector p4NHL, p4NHL_rest;
  for(std::vector<int>::const_iterator pdg_iter = decayList.begin(); pdg_iter != decayList.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);
     if( std::abs( pdgc ) == kPdgNHL ) p4NHL = *p4fin;
     idp++;
  }
  
  return p4NHL; // rest frame momentum!
}
//----------------------------------------------------------------------------
TVector3 NHLFluxCreator::PointToRandomPointInBBox( TVector3 detO_beam )
{
  RandomGen * rnd = RandomGen::Instance();
  double ox = detO_beam.X(), oy = detO_beam.Y(), oz = detO_beam.Z();
  double rx = (rnd->RndGen()).Uniform( ox - fLx/2.0, ox + fLx/2.0 ), 
         ry = (rnd->RndGen()).Uniform( oy - fLy/2.0, oy + fLy/2.0 ),
         rz = (rnd->RndGen()).Uniform( oz - fLz/2.0, oz + fLz/2.0 );
  TVector3 vec( rx, ry, rz );
  LOG( "NHL", pDEBUG )
    << "Pointing to this point in BBox (beam coords): " << utils::print::Vec3AsString( &vec );
  return vec;
}
//----------------------------------------------------------------------------
double NHLFluxCreator::GetAngDeviation( TLorentzVector p4par, TVector3 detO, bool seekingMax )
{
  TVector3 ppar = p4par.Vect(); assert( ppar.Mag() > 0.0 );
  TVector3 pparUnit = ppar.Unit();
  // px . (x-x0) + py . (y-y0) + pz . (z-z0) = 0
  // ==> px . x + py . y + pz . z = ppar . detO
  double dterm = pparUnit.X() * detO.X() + pparUnit.Y() * detO.Y() + pparUnit.Z() * detO.Z();
  // trajectory parallel to ppar and passes through (0,0,0)
  // ==> v(t) = ( px, py, pz ) t
  // ==> (px^2 + py^2 + pz^2) t = ppar . detO
  // ==> t = ppar . detO / ppar.Mag2()
  double t = dterm / 1.0;
  double x_incp = t * pparUnit.X(), y_incp = t * pparUnit.Y(), z_incp = t * pparUnit.Z();

  // sweep over plane perp to ppar, passing through centre, and calc intersection point
  // special case: parent is perfectly on axis so hits detector centre
  TVector3 IPdev( detO.X() - x_incp, detO.Y() - y_incp, detO.Z() - z_incp );
  bool parentHitsCentre = ( IPdev.Mag() < controls::kASmallNum );

  // RETHERE: Approximating detector dimension as BBox z dimension (NHL running on almost-z)
  double detRadius = fLz / 2.0;

  if( parentHitsCentre ){
    // calculate angles for four points and return largest (smallest) of them

    // randomly select a phi to go with, make 2 perpendicular vectors from it
    double phi = ( RandomGen::Instance()->RndGen() ).Uniform( 0., 2.0 * constants::kPi );
    TVector3 r1VecPrim( -pparUnit.Y(), pparUnit.X(), 0.0 ); // perp to ppar == on plane
    TVector3 r1Vec = r1VecPrim.Unit();
    r1Vec.Rotate( phi, pparUnit );
    TVector3 r2Vec( r1Vec ); r2Vec.Rotate( 0.5 * constants::kPi, pparUnit );

    double rprod = r1Vec.X() * r2Vec.X() + r1Vec.Y() * r2Vec.Y() + r1Vec.Z() * r2Vec.Z();

    assert( std::abs( rprod ) < controls::kASmallNum );

    // four IP with det. All have distance detRadius from centre.
    TVector3 p1( detO.X() + detRadius*r1Vec.X(),
		 detO.Y() + detRadius*r1Vec.Y(),
		 detO.Z() + detRadius*r1Vec.Z() );
    TVector3 p2( detO.X() - detRadius*r1Vec.X(),
		 detO.Y() - detRadius*r1Vec.Y(),
		 detO.Z() - detRadius*r1Vec.Z() );
    TVector3 p3( detO.X() + detRadius*r2Vec.X(),
		 detO.Y() + detRadius*r2Vec.Y(),
		 detO.Z() + detRadius*r2Vec.Z() );
    TVector3 p4( detO.X() - detRadius*r2Vec.X(),
		 detO.Y() - detRadius*r2Vec.Y(),
		 detO.Z() - detRadius*r2Vec.Z() );

    // Return largest(smallest) angle using inner product magic
    double thLarge = -999.9; double thSmall = 999.9;
    double th1 = TMath::ACos( ( p1.X()*pparUnit.X() + p1.Y()*pparUnit.Y() + p1.Z()*pparUnit.Z() ) / p1.Mag() ); if( th1 > thLarge ){ thLarge = th1; } else if( th1 < thSmall ){ thSmall = th1; }
    double th2 = TMath::ACos( ( p2.X()*pparUnit.X() + p2.Y()*pparUnit.Y() + p2.Z()*pparUnit.Z() ) / p2.Mag() ); if( th2 > thLarge ){ thLarge = th2; } else if( th2 < thSmall ){ thSmall = th2; }
    double th3 = TMath::ACos( ( p3.X()*pparUnit.X() + p3.Y()*pparUnit.Y() + p3.Z()*pparUnit.Z() ) / p3.Mag() ); if( th3 > thLarge ){ thLarge = th3; } else if( th3 < thSmall ){ thSmall = th3; }
    double th4 = TMath::ACos( ( p4.X()*pparUnit.X() + p4.Y()*pparUnit.Y() + p4.Z()*pparUnit.Z() ) / p4.Mag() ); if( th4 > thLarge ){ thLarge = th4; } else if( th4 < thSmall ){ thSmall = th4; }

    return ( seekingMax ) ? thLarge * 180.0 / constants::kPi : thSmall * 180.0 / constants::kPi;
  } else {
    // find direction from IP to det centre.
    TVector3 rVec = IPdev.Unit(); double dist = IPdev.Mag();
    // two IP with det. Closer(farther) has distance detRadius -(+) d( IP, detO )
    double dh = detRadius + dist, dl = detRadius - dist;
    // get those vectors and do inner product magic
    TVector3 ph( x_incp + dh * rVec.X(), y_incp + dh * rVec.Y(), z_incp + dh * rVec.Z() );
    TVector3 pl( x_incp - dl * rVec.X(), y_incp - dl * rVec.Y(), z_incp - dl * rVec.Z() );
    double thh = TMath::ACos( ( ph.X()*pparUnit.X() + ph.Y()*pparUnit.Y() + ph.Z()*pparUnit.Z() ) / ph.Mag() );
    double thl = TMath::ACos( ( pl.X()*pparUnit.X() + pl.Y()*pparUnit.Y() + pl.Z()*pparUnit.Z() ) / pl.Mag() );

    return ( seekingMax ) ? thh * 180.0 / constants::kPi : thl * 180.0 / constants::kPi;
  }

  LOG( "NHL", pERROR )
    << "Could not calculate the angle range for detector at ( " << detO.X() << ", " 
    << detO.Y() << ", " << detO.Z() << " ) [m] from NHL production point with parent momentum = ( "
    << ppar.X() << ", " << ppar.Y() << ", " << ppar.Z() << " ) [GeV]. Returning zero.";
  return 0.0;
}
//----------------------------------------------------------------------------
double NHLFluxCreator::CalculateAcceptanceCorrection( TLorentzVector p4par, TLorentzVector p4NHL,
						      double SMECM, double zm, double zp )
{
  /*
   * This method calculates NHL acceptance by taking into account the collimation effect
   * NHL are massive so Lorentz boost from parent CM ==> lab is more effective
   * This means that, given a desired range of lab-frame emission angles, there are
   * more rest-frame emission angles that map into this range. 
   * Find the measure of the rest-frame that maps onto the allowed lab-frame angles
   * and return the ratio over the relevant measure for a SM neutrino
   */
  
  assert( zm >= 0.0 && zp >= zm );
  if( zp == zm ) return 0.0;

  double M = p4NHL.M();
  if( M == 0.0 ) return 1.0;

  fNHL = ( TF1* ) gROOT->GetListOfFunctions()->FindObject( "fNHL" );
  if( !fNHL ){ fNHL = new TF1( "fNHL", labangle, 0.0, 180.0, 6 ); }
  fNHL->SetParameter( 0, p4par.E()  );
  fNHL->SetParameter( 1, p4par.Px() );
  fNHL->SetParameter( 2, p4par.Py() );
  fNHL->SetParameter( 3, p4par.Pz() );
  fNHL->SetParameter( 4, p4NHL.P()  );
  fNHL->SetParameter( 5, p4NHL.E()  );

  double ymax = fNHL->GetMaximum(), xmax = fNHL->GetMaximumX();
  double range1 = 0.0;

  if( fNHL->GetMinimum() == fNHL->GetMaximum() ) return 1.0; // bail on constant function

  if( zm < fNHL->GetMinimum() ){ // really good collimation, ignore checks on zm
    double z0 = fNHL->GetMinimum();
    if( ymax > zp && xmax < 180.0 ){ // >=2 pre-images, add them together
      int nPreim = 0;

      // RETHERE: Make this more sophisticated! Assumes 2 preimages, 1 before and 1 after max
      double xl1 = fNHL->GetX( z0, 0.0, xmax );
      double xh1 = fNHL->GetX( zp, 0.0, xmax );
      double xl2 = fNHL->GetX( z0, xmax, 180.0 );
      double xh2 = fNHL->GetX( zp, xmax, 180.0 );

      range1 += std::abs( xl1 - xh1 ) + std::abs( xh2 - xl2 ); nPreim = 2;
    } else if( ymax > zp && xmax == 180.0 ){ // 1 pre-image, SMv-like case
      double xl = fNHL->GetX( z0 ), xh = fNHL->GetX( zp );
      range1 = std::abs( xh - xl );
    } else if( ymax <= zp ){ // 1 pre-image but all emissions reach detector
      range1 = 180.0;
    }
  } else { // not so good collimation, enforce checks on zm
    if( ymax <= zm ){ // 0 pre-images
      return 0.0;
    } else if( ymax > zp && xmax < 180.0 ){ // >=2 pre-images, add them together
      int nPreim = 0;

      // RETHERE: Make this more sophisticated! Assumes 2 preimages, 1 before and 1 after max
      double xl1 = fNHL->GetX( zm, 0.0, xmax );
      double xh1 = fNHL->GetX( zp, 0.0, xmax );
      double xl2 = fNHL->GetX( zm, xmax, 180.0 );
      double xh2 = fNHL->GetX( zp, xmax, 180.0 );

      range1 += std::abs( xl1 - xh1 ) + std::abs( xh2 - xl2 ); nPreim = 2;
    } else if( ymax > zp && xmax == 180.0 ){ // 1 pre-image, SMv-like case
      double xl = fNHL->GetX( zm ), xh = fNHL->GetX( zp );
      range1 = std::abs( xh - xl );
    } else if( zm < ymax && ymax <= zp ){ // 1 pre-image
      double xl = fNHL->GetX( zm, 0., xmax ), xh = fNHL->GetX( zm, xmax, 180.0 );
      range1 = std::abs( xh - xl );
    }
  }

  fSMv = ( TF1* ) gROOT->GetListOfFunctions()->FindObject( "fSMv" );
  if( !fSMv ){
    fSMv = new TF1( "fSMv", labangle, 0.0, 180.0, 6 );
  }
  fSMv->SetParameter( 0, p4par.E()  );
  fSMv->SetParameter( 1, p4par.Px() );
  fSMv->SetParameter( 2, p4par.Py() );
  fSMv->SetParameter( 3, p4par.Pz() );
  fSMv->SetParameter( 4, SMECM      );
  fSMv->SetParameter( 5, SMECM      );
  double range2 = -1.0;

  if( fSMv->GetMaximum() == fSMv->GetMinimum() ) return 1.0; // bail

  // SMv deviates more from parent than NHL due to masslessness. This means a larger minimum of labangle
  // Sometimes the angle is so small, that this calculation fails as there is not SMv preimage to compare
  // with. Default to estimating dx/dz * dy/dz ratio in that case.

  if( fSMv->GetMinimum() < zp ){
    if( fSMv->GetMinimum() < zm ){
      range2 = std::abs( fSMv->GetX( zp ) - fSMv->GetX( zm ) );
    } else { // due to monotonicity all of [0.0, fSMv->GetX(zp)] is good
      range2 = fSMv->GetX( zp );
      if( range2 == 0 ) return 1.0; // sometimes this happens, gotta bail!
    }
  } else { // can't decide based on SMv analytically.

    TVector3 bv = p4par.BoostVector();
    
    TLorentzVector vcx( p4NHL.P(), 0.0, 0.0, p4NHL.E() ), 
      vcy( 0.0, p4NHL.P(), 0.0, p4NHL.E() ), 
      vcz( 0.0, 0.0, p4NHL.P(), p4NHL.E() );
    vcx.Boost( bv ); vcy.Boost( bv ); vcz.Boost( bv );
    
    TLorentzVector vsx( SMECM, 0.0, 0.0, SMECM ),
      vsy( 0.0, SMECM, 0.0, SMECM ),
      vsz( 0.0, 0.0, SMECM, SMECM );
    vsx.Boost( bv ); vsy.Boost( bv ); vsz.Boost( bv );
    
    double xpart = std::abs( ( vcx.X() / vcz.Z() ) / ( vsx.X() / vsz.Z() ) );
    double ypart = std::abs( ( vcy.Y() / vcz.Z() ) / ( vsy.Y() / vsz.Z() ) );
    return 1.0 / ( xpart * ypart );
  }

  assert( range2 > 0.0 );
  return range1 / range2;

}
//----------------------------------------------------------------------------
double NHLFluxCreator::labangle( double * x, double * par )
{
  double xrad = x[0] * TMath::DegToRad();
  double Ehad = par[0], pxhad = par[1], pyhad = par[2], pzhad = par[3];
  double pnhl = par[4], Enhl = par[5];

  TLorentzVector p4had( pxhad, pyhad, pzhad, Ehad );
  TVector3 boost_vec = p4had.BoostVector(); // beta of parent in lab frame

  // assume phi invariance so create NHL rest-frame momentum along y'z' plane
  TLorentzVector pncm( 0.0, pnhl * TMath::Sin( xrad ), pnhl * TMath::Cos( xrad ), Enhl );

  // boost into lab frame
  pncm.Boost( boost_vec );
  
  // return lab frame theta wrt parent momentum in deg
  // double theta = TMath::ACos( pncm.Pz() / pncm.P() ) * 180.0 / constants::kPi;
  double num = pxhad * pncm.X() + pyhad * pncm.Y() + pzhad * pncm.Z();
  double den = p4had.P() * pncm.P();
  double theta = TMath::ACos( num / den ) * 180.0 / constants::kPi;
  return theta;
}
//----------------------------------------------------------------------------
void NHLFluxCreator::MakeBBox()
{
  if( isBoxInit ) return;

  LOG( "NHL", pWARN )
    << "WARNING: This is a dummy (==unit-side) bounding box centred at config-given point";

  // read config
  fCx = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "DetCentreXInBeam" );
  fCy = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "DetCentreYInBeam" );
  fCz = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "DetCentreZInBeam" );

  fAx1 = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "EulerExtrinsicX1" );
  fAz  = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "EulerExtrinsicZ"  );
  fAx2 = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "EulerExtrinsicX2" );

  fLx = 1.0; fLy = 1.0; fLz = 1.0;

  isBoxInit = true;
}
//----------------------------------------------------------------------------
TVector3 NHLFluxCreator::ApplyUserRotation( TVector3 vec, bool doBackwards )
{
  double vx = vec.X(), vy = vec.Y(), vz = vec.Z();

  double Ax2 = ( doBackwards ) ? -fAx2 : fAx2;
  double Az  = ( doBackwards ) ? -fAz  : fAz;
  double Ax1 = ( doBackwards ) ? -fAx1 : fAx1;

  // Ax2 first
  double x = vx, y = vy, z = vz;
  vy = y * std::cos( Ax2 ) - z * std::sin( Ax2 );
  vz = y * std::sin( Ax2 ) + z * std::cos( Ax2 );
  y = vy; z = vz;
  // then Az
  vx = x * std::cos( Az )  - y * std::sin( Az );
  vy = x * std::sin( Az )  + y * std::cos( Az );
  x = vx; y = vy;
  // Ax1 last
  vy = y * std::cos( Ax1 ) - z * std::sin( Ax1 );
  vz = y * std::sin( Ax1 ) + z * std::cos( Ax1 );

  TVector3 nvec( vx, vy, vz );
  return nvec;
}
//----------------------------------------------------------------------------
double NHLFluxCreator::CalculateDetectorAcceptanceSAA( TVector3 detO )
{
  // sang is solid-angle / 4pi
  double rad = std::sqrt( detO.X() * detO.X() + detO.Y() * detO.Y() + detO.Z() * detO.Z() );
  double sang = 1.0 - TMath::Cos( TMath::ATan( kRDET / rad ) ); sang *= 0.5;
  return sang;
}
//----------------------------------------------------------------------------
double NHLFluxCreator::CalculateDetectorAcceptanceDRC( TVector3 detO, double Lx, double Ly, double Lz )
{
  /* This method does two steps to calculate the acceptance
   * 1 : Partition the BBox into a number of smaller boxes.
   *     Project the centre of each of these boxes onto the unit sphere as P = (\theta, \phi)
   * 2 : Partition the unit sphere into parallels and meridians
   *     For each area that contains at least one P_i, calculate the volume element
   *     A = \int_{\theta_min}^{\theta_{max}}d\theta \int_{\phi_{min}}^{\phi_{max}} d\phi sin\theta
   *     = \Delta\phi * ( cos(\theta_{min}) - cos(\theta_{max}) ) 
   *     and add that to the solid angle (max = 4pi)
   * // RETHERE STEP 3 SHOULD BE ITS OWN METHOD (need to know flux area-norm)
   *    See CalculateAreaNormalisation()
   * 3 : Project the entire BBox onto the plane tangent to the unit sphere at the point
   *     O' := P( O ) where O the centre of the BBox.
   *     Calculate its area iteratively, and scale up by r(O).
   *     This gives an estimate for the effective area of the BBox.
   */
  LOG( "NHL", pDEBUG )
    << "Starting the calculation with DRC using O = ( "
    << detO.X() << ", " << detO.Y() << ", " << detO.Z() << " ) [m] and lengths = ( " 
    << Lx << ", " << Ly << ", " << Lz << " ) [m]. . .";

  // the lower face
  double xA = detO.X() - Lx/2.0, yA = detO.Y() - Ly/2.0, zA = detO.Z() - Lz/2.0;
  double xB = detO.X() + Lx/2.0, yB = detO.Y() - Ly/2.0, zB = detO.Z() - Lz/2.0;
  double xC = detO.X() + Lx/2.0, yC = detO.Y() + Ly/2.0, zC = detO.Z() - Lz/2.0;
  double xD = detO.X() - Lx/2.0, yD = detO.Y() + Ly/2.0, zD = detO.Z() - Lz/2.0;
  // the upper face is symmetric
  double xE = detO.X() - Lx/2.0, yE = detO.Y() - Ly/2.0, zE = detO.Z() + Lz/2.0;
  double xF = detO.X() + Lx/2.0, yF = detO.Y() - Ly/2.0, zF = detO.Z() + Lz/2.0;
  double xG = detO.X() + Lx/2.0, yG = detO.Y() + Ly/2.0, zG = detO.Z() + Lz/2.0;
  double xH = detO.X() - Lx/2.0, yH = detO.Y() + Ly/2.0, zH = detO.Z() + Lz/2.0;

  LOG( "NHL", pDEBUG )
    << "\nEdge points for this BBox are:"
    << "\nA = ( " << xA << ", " << yA << ", " << zA << " )"
    << "\nB = ( " << xB << ", " << yB << ", " << zB << " )"
    << "\nC = ( " << xC << ", " << yC << ", " << zC << " )"
    << "\nD = ( " << xD << ", " << yD << ", " << zD << " )"
    << "\nE = ( " << xE << ", " << yE << ", " << zE << " )"
    << "\nF = ( " << xF << ", " << yF << ", " << zF << " )"
    << "\nG = ( " << xG << ", " << yG << ", " << zG << " )"
    << "\nH = ( " << xH << ", " << yH << ", " << zH << " )";

  // hard-code this partition
  // double nx = 10.0, ny = 10.0, nz = 10.0;
  double nx = utils::nhl::GetCfgDouble( "NHL", "FluxCalc", "Raycast-PartitionX" );
  double ny = utils::nhl::GetCfgDouble( "NHL", "FluxCalc", "Raycast-PartitionY" );
  double nz = utils::nhl::GetCfgDouble( "NHL", "FluxCalc", "Raycast-PartitionZ" );
  int nix = nx, niy = ny, niz = nz; int nia = nix * niy * niz;

  LOG( "NHL", pDEBUG )
    << "Step 1: Partition the BBox into " << nix << " x " << niy << " x " << niz << " sub-boxes.";
  double xs[nix], ys[niy], zs[niz];
  double centres[nia][3];
  double thetas[nia], phis[nia]; // avoid vectors for speed
  std::array< double, 2 > points[nia]; // 2D array with sort capability. Note it's in degrees!
  for( int ix = 0; ix < nix; ix++ ){ xs[ix] = xA + (ix+0.5) * Lx / nx; }
  for( int iy = 0; iy < niy; iy++ ){ ys[iy] = yA + (iy+0.5) * Ly / ny; }
  for( int iz = 0; iz < niz; iz++ ){ zs[iz] = zA + (iz+0.5) * Lz / nz; }
  LOG( "NHL", pDEBUG )
    << "The small-(x,y,z) centre is at ( " << xs[0] << ", " << ys[0] << ", " << zs[0] << " )";

  for( int ix = 0; ix < nix; ix++ ){
    for( int iy = 0; iy < niy; iy++ ){
      for( int iz = 0; iz < niz; iz++ ){
	int ia = ix*niy*niz + iy*niz +iz;
	centres[ia][0] = xs[ix];
	centres[ia][1] = ys[iy];
	centres[ia][2] = zs[iz];
      }
    }
  }

  double minTheta = 9999.9, maxTheta = -9999.9;
  double minPhi = 9999.0, maxPhi = -9999.9;

  for( int ia = 0; ia < nia; ia++ ){ points[ia] = { 0.0, 0.0 }; }

  for( int ia = 0; ia < nia; ia++ ){
    double rad = TMath::Sqrt(
			     TMath::Power( centres[ia][0], 2.0 ) +
			     TMath::Power( centres[ia][1], 2.0 ) +
			     TMath::Power( centres[ia][2], 2.0 ) );
    thetas[ia] = TMath::ACos( centres[ia][2] / rad );
    phis[ia] = TMath::ACos( centres[ia][0] / ( rad * TMath::Sin( thetas[ia] ) ) );
    points[ia] = { thetas[ia] * 180.0 / constants::kPi, phis[ia] * 180.0 / constants::kPi };
    
    if( thetas[ia] > maxTheta ) maxTheta = thetas[ia];
    if( thetas[ia] < minTheta ) minTheta = thetas[ia];
    if( phis[ia] > maxPhi ) maxPhi = phis[ia];
    if( phis[ia] < minPhi ) minPhi = phis[ia];
  }

  // sort by theta
  double tmpmin = points[0][0], tmpmax = points[0][0]; int tmpi = -1;
  for( int ia = 1; ia < nia; ia++ ){
    //std::ostringstream osts, psts;
    tmpi = -1; // to avoid nonsense
    //LOG( "NHL", pDEBUG ) << "tmpmin, tmpmax = " << tmpmin << ", " << tmpmax;
    if( points[ia][0] >= tmpmax ){ // if greater than max, then must be end point! Don't reorganise
      tmpmax = points[ia][0]; 
      //LOG( "NHL", pDEBUG ) << "tmpmax is now " << tmpmax << " and tmpi = -1";
    } else if( points[ia][0] < tmpmin ){ // if smaller than min, then must be first point! Set tmpi to 0
      tmpmin = points[ia][0]; tmpi = 0; 
      //LOG( "NHL", pDEBUG ) << "tmpmin is now " << tmpmin << " and tmpi = 0";
    } else { // find where the push must go and push
      tmpi = ia;
      while( points[ia][0] < points[tmpi-1][0] && tmpi > 0 ){ 
	//psts << "\n@tmpi = " << tmpi << ": " << points[tmpi-1][0]; 
	tmpi--; 
      }
      //psts << "\nFINAL: tmpi = " << tmpi;
      //LOG( "NHL", pDEBUG ) << (psts.str()).c_str();
    }
    //LOG( "NHL", pDEBUG ) << "ia = " << ia << ", tmpi = " << tmpi << ", theta = " << points[ia][0];
    // now push every entry in [tmpi, ia-1] to +1 and insert ia in tmpi
    double tmptheta = points[ia][0]; double tmpphi = points[ia][1];
    for( int ib = ia; ib > tmpi; ib-- ){
      points[ib] = points[ib-1];
    }
    points[tmpi] = { tmptheta, tmpphi };

    /*
    for( int is = 0; is <= ia; is++ ){
      osts << "\n[ " << is << " ] " << points[is][0]; 
    }
    LOG( "NHL", pDEBUG ) << (osts.str()).c_str();
    */
  }

  /*
  std::ostringstream asts;
  // check sorting
  for( int ia = 0; ia < nia; ia++ ){
    asts << "\n[" << ia << "]: " << points[ia][0];
  }
  LOG( "NHL", pDEBUG ) << (asts.str()).c_str();
  */

  // sort these arrays by theta
  // See https://stackoverflow.com/questions/20931669/sort-a-2d-array-in-c-using-built-in-functionsor-any-other-method
  //std::sort( std::begin( points ), std::end( points ) );

  // now the spherical partition
  //double ntheta = 18000.0, nphi = 36000.0;
  double ntheta = utils::nhl::GetCfgDouble( "NHL", "FluxCalc", "Raycast-NParallels" );
  double nphi   = utils::nhl::GetCfgDouble( "NHL", "FluxCalc", "Raycast-NMeridians" );
  int nitheta = ntheta, niphi = nphi;

  LOG( "NHL", pDEBUG )
    << "Step 2: Partition the unit sphere into " << nitheta << " parallels and " << niphi << " meridians";
  LOG( "NHL", pDEBUG )
    << "\nThe minimum theta is " << minTheta * 180.0 / constants::kPi << " deg,"
    << " and the maximum theta is " << maxTheta * 180.0 / constants::kPi << " deg"
    << "\nThe minimum phi is " << minPhi * 180.0 / constants::kPi << " deg,"
    << " and the maximum phi is " << maxPhi * 180.0 / constants::kPi << " deg";

  double areaElem = 0.0;
  for( int ith = 0; ith < nitheta; ith++ ){
    double thetaSmall = ith / ntheta * 180.0;
    double thetaLarge = (ith+1) / ntheta * 180.0;
    if( minTheta * 180.0 / constants::kPi > thetaLarge ) continue;
    if( maxTheta * 180.0 / constants::kPi < thetaSmall ) break; // done!
    for( int iph = 0; iph < niphi; iph++ ){
      double phiSmall = iph / nphi * 360.0;
      double phiLarge = (iph+1) / nphi * 360.0;
      if( minPhi * 180.0 / constants::kPi > phiLarge ) continue;
      if( maxPhi * 180.0 / constants::kPi < phiSmall ) break; // done, on to next theta
      LOG( "NHL", pDEBUG ) << "theta in [ " << thetaSmall << ", " << thetaLarge
			   << " ]; phi in [ " << phiSmall << ", " << phiLarge << " ]";
      // tag this volume-element only if we haven't seen this region before to avoid double-counting
      bool foundPoint = false;
      int npoints = 0;
      while( !foundPoint && npoints <= nia ){
	foundPoint = ( points[npoints][0] <= thetaLarge &&
		       points[npoints][0] >= thetaSmall &&
		       points[npoints][1] <= phiLarge   &&
		       points[npoints][1] >= phiSmall );
	npoints++;
      }
      if( foundPoint ){ // calculate area element in sterad and add it
	double dPhi = ( phiLarge - phiSmall ) * constants::kPi / 180.0 ;
	double dTheta = TMath::Cos( thetaSmall * constants::kPi / 180.0 ) - TMath::Cos( thetaLarge * constants::kPi / 180.0 );
	areaElem += dPhi * dTheta;
	LOG( "NHL", pDEBUG )
	  << "Tagging element: area = " << dPhi * dTheta;
	continue; // exit out of this region and onto the next
      }
      if( foundPoint ){ // you should never be seeing this, double counting WILL happen here.
	LOG( "NHL", pERROR )
	  << "Double-counting is occurring here. Crashing now...";
	assert( false );
      }
    }
  }

  LOG( "NHL", pDEBUG )
    << "Solid angle = " << areaElem << " sterad";

  // return this / 4pi == acceptance
  return areaElem / (4.0 * constants::kPi);
}
//----------------------------------------------------------------------------
double NHLFluxCreator::CalculateAreaNormalisation()
{
  // for now this is just a square of length kRDET
  // returns 1 / area
  return 1.0 / ( kRDET * kRDET );
}
