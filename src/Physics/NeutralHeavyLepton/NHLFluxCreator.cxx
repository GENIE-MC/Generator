//----------------------------------------------------------------------------
/*!

  Implementation of NHLFluxCreator

 */
//----------------------------------------------------------------------------

#include "Physics/NeutralHeavyLepton/NHLFluxCreator.h"

using namespace genie;
using namespace genie::NHL;
using namespace genie::NHL::NHLenums;

//----------------------------------------------------------------------------
NHLFluxCreator::NHLFluxCreator() :
  EventRecordVisitorI("genie::NHL::NHLFluxCreator")
{

}
//----------------------------------------------------------------------------
NHLFluxCreator::NHLFluxCreator(string config) :
  EventRecordVisitorI("genie::NHL::NHLFluxCreator", config)
{

}
//----------------------------------------------------------------------------
NHLFluxCreator::~NHLFluxCreator()
{

}
//----------------------------------------------------------------------------
void NHLFluxCreator::ProcessEventRecord(GHepRecord * evrec) const
{
  // Adds the inital state NHL at the event record.
  // Also assigns the production vertex to evrec (this will be overwritten by subsequent modules)
  // Also adds (acceptance*nimpwt)^(-1) component of weight

  if( iCurrEntry < fFirstEntry ) iCurrEntry = fFirstEntry;

  LOG( "NHL", pDEBUG ) << "Flux record processing now, current entry = " << iCurrEntry;

  if( iCurrEntry >= fFirstEntry ) {

    //flux::GNuMIFluxPassThroughInfo * gnmf = new flux::GNuMIFluxPassThroughInfo();
    if( !fGnmf ) fGnmf = new flux::GNuMIFluxPassThroughInfo();
    this->MakeTupleFluxEntry( iCurrEntry, fGnmf, fCurrPath );
  
    if( std::abs(fGnmf->fgPdgC) == genie::kPdgNHL ){ // only add particle if parent is valid
      
      double invAccWeight = fGnmf->nimpwt * fGnmf->fgXYWgt;
      evrec->SetWeight( evrec->Weight() / invAccWeight );
      
      // scale by how many POT it takes to make the appropriate parent
      evrec->SetWeight( evrec->Weight() * POTScaleWeight );
      
      evrec->SetVertex( fGnmf->fgX4User ); // NHL production vertex. NOT where NHL decays to visible FS.
      // construct Particle(0). Don't worry about daughter links at this stage.
      TLorentzVector probeP4 = ( fUseBeamMomentum ) ? fGnmf->fgP4 : fGnmf->fgP4User;
      if( fUseBeamMomentum ){
	TVector3 probeP3 = probeP4.Vect(); TVector3 dumor( 0.0, 0.0, 0.0 );
	probeP3 = this->ApplyUserRotation( probeP3, dumor, fDetRotation, false ); // active transformation
	probeP4.SetPxPyPzE( probeP3.X(), probeP3.Y(), probeP3.Z(), probeP4.E() );
      }
      GHepParticle ptNHL( fGnmf->fgPdgC, kIStInitialState, -1, -1, -1, -1, probeP4, fGnmf->fgX4User );
      evrec->AddParticle( ptNHL );
    }
    
    // clean up
    //delete gnmf;

  }
  
  // update iCurrEntry. This should be the last thing to happen before returning.
  // iCurrEntry++;
}
//----------------------------------------------------------------------------
void NHLFluxCreator::SetInputPath(std::string finpath) const
{
  LOG( "NHL", pDEBUG ) << "Setting input path to " << finpath;
  LOG( "NHL", pDEBUG ) << "Before setting, fCurrPath = " << fCurrPath;
  fCurrPath = finpath;
}
//----------------------------------------------------------------------------
int NHLFluxCreator::GetNEntries() const
{
  LOG( "NHL", pDEBUG ) << "fNEntries = " << fNEntries;
  if( fNEntries <= 0 ){
    this->OpenFluxInput( fCurrPath );
  }
  return fNEntries;
}
//----------------------------------------------------------------------------
void NHLFluxCreator::SetFirstEntry( int iFirst ) const
{
  fFirstEntry = iFirst;
}
//----------------------------------------------------------------------------
void NHLFluxCreator::SetCurrentEntry( int iCurr ) const
{
  iCurrEntry = iCurr;
}
//----------------------------------------------------------------------------
flux::GNuMIFluxPassThroughInfo * NHLFluxCreator::RetrieveGNuMIFluxPassThroughInfo() const
{
  return fGnmf;
}
//----------------------------------------------------------------------------
void NHLFluxCreator::MakeTupleFluxEntry( int iEntry, flux::GNuMIFluxPassThroughInfo * gnmf, std::string finpath ) const
{
  // This method creates 1 NHL from the flux info and saves the information
  // Essentially, it replaces a SMv with an NHL

  // Open flux input and initialise trees
  if( iEntry == fFirstEntry ){
    this->OpenFluxInput( finpath );
    this->InitialiseTree();
    this->InitialiseMeta();
    this->MakeBBox();
  } else if( iEntry < fFirstEntry ){
    return;
  }

  // All these in m
  TVector3 fCvec_beam( fCx, fCy, fCz );
  TVector3 fCvec = this->ApplyUserRotation( fCvec_beam );
  fLepPdg = 0;
  
  LOG( "NHL", pDEBUG ) << "Getting entry " << iEntry;
  ctree->GetEntry(iEntry);

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

  LOG( "NHL", pDEBUG )
    << "\niEntry       = " << iEntry
    << "\nnimpwt       = " << decay_nimpwt
    //<< "\npotnum       = " << potnum
    //<< "\ndeltaPotnum  = " << deltaPotnum
    << "\ncanGoForward ? " << canGoForward;
    //<< "\nmultiplicity = " << 1.0 / multiplicity;

  if( !canGoForward ){ // return, but update the NPOT it took to make this rejected parent
    this->FillNonsense( iEntry, gnmf ); return; 
  }
    
  // turn cm to m and make origin wrt detector 
  fDx = decay_vx * units::cm / units::m;
  fDy = decay_vy * units::cm / units::m;
  fDz = decay_vz * units::cm / units::m;
  TVector3 fDvec( fDx, fDy, fDz );
  TVector3 fDvec_beam = this->ApplyUserRotation( fDvec, true );

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
  
  double acc_saa = this->CalculateDetectorAcceptanceSAA( detO_user );
  
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

  // CAUTION: p4par is in USER coords
  TLorentzVector p4par = ( isParentOnAxis ) ? 
    TLorentzVector( parentMomentum * (detO.Unit()).X(), 
		    parentMomentum * (detO.Unit()).Y(),
		    parentMomentum * (detO.Unit()).Z(),
		    parentEnergy ) :
    TLorentzVector( decay_pdpx, decay_pdpy, decay_pdpz, parentEnergy );

  TLorentzVector p4par_beam( decay_pdpx, decay_pdpy, decay_pdpz, parentEnergy );
  if( !isParentOnAxis ){
    // rotate p4par to user coordinates
    TVector3 tmpv3 = ApplyUserRotation( p4par.Vect() );
    p4par.SetPxPyPzE( tmpv3.Px(), tmpv3.Py(), tmpv3.Pz(), p4par.E() );
  }


  TVector3 boost_beta = p4par.BoostVector();

  // now calculate which decay channel produces the NHL.
  dynamicScores = this->GetProductionProbs( decay_ptype );
  assert( dynamicScores.size() > 0 );
  
  if( dynamicScores.find( kNHLProdNull ) != dynamicScores.end() ){ // exists kin allowed channel but 0 coupling
    this->FillNonsense( iEntry, gnmf ); return;
  }
  
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

  // we will now boost detO into rest frame, force rest to point to the new direction, boost the result, and compare the boost corrections
  double boost_correction_two = 0.0;
  
  // 17-Jun-22: Notice the time component needs to be nonzero to get this to work!

  // first guess: betaNHL ~= 1 . Do the Lorentz boosts knocking betaNHL downwards until we hit det centre
  double betaStar  = p4NHL_rest.P() / p4NHL_rest.E();
  double betaMag = boost_beta.Mag();
  double gamma   = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
  double betaLab = 1.0; // first guess

  // now make a TLorentzVector in lab frame to boost back to rest. 
  double timeBit = detO.Mag() / units::kSpeedOfLight ; // s
  TLorentzVector detO_4v( detO.X(), detO.Y(), detO.Z(), timeBit ); detO_4v.Boost( -boost_beta );
  TVector3 detO_rest_unit = (detO_4v.Vect()).Unit();
  TLorentzVector p4NHL_rest_good( p4NHL_rest.P() * detO_rest_unit.X(),
				  p4NHL_rest.P() * detO_rest_unit.Y(),
				  p4NHL_rest.P() * detO_rest_unit.Z(),
				  p4NHL_rest.E() );

  // boost that into lab frame!
  TLorentzVector p4NHL_good = p4NHL_rest_good;
  p4NHL_good.Boost( boost_beta );
  boost_correction_two = p4NHL_good.E() / p4NHL_rest.E();

  TVector3 detO_unit = detO.Unit();

  TVector3 p4NHL_good_vect = p4NHL_good.Vect();
  TVector3 p4NHL_good_unit = p4NHL_good_vect.Unit();

  // now calculate how far away from det centre we are.
  // dist = || ( rC - rD ) x p4NHL || / || p4NHL_good || where x == cross product
  TVector3 distNum = detO.Cross( p4NHL_good_unit );
  double dist = distNum.Mag(); // m

  while( dist >= 1.0e-3 && betaLab > 0.0 ){ // 1mm tolerance
    // that didn't work. Knock betaLab down a little bit and try again.
    betaLab -= 1.0e-3;
    timeBit = detO.Mag() / ( betaLab * units::kSpeedOfLight );
    detO_4v.SetXYZT( detO.X(), detO.Y(), detO.Z(), timeBit );
    detO_4v.Boost( -boost_beta );
    detO_rest_unit = (detO_4v.Vect()).Unit();
    p4NHL_rest_good.SetPxPyPzE( p4NHL_rest.P() * detO_rest_unit.X(),
				p4NHL_rest.P() * detO_rest_unit.Y(),
				p4NHL_rest.P() * detO_rest_unit.Z(),
				p4NHL_rest.E() );

    // boost into lab frame
    p4NHL_good = p4NHL_rest_good;
    p4NHL_good.Boost( boost_beta );
    
    detO_unit = detO.Unit();
    p4NHL_good_vect = p4NHL_good.Vect();
    p4NHL_good_unit = p4NHL_good_vect.Unit();

    distNum = detO.Cross( p4NHL_good_unit );
    dist = distNum.Mag(); // m
  }

  LOG( "NHL", pDEBUG )
    << "\nTimelike bit = " << timeBit
    << "\ndetO = " << utils::print::Vec3AsString( &detO ) << " [m]"
    << "\np4NHL_good = " << utils::print::P4AsString( &p4NHL_good ) << " GeV"
    << "\ndistNum = " << utils::print::Vec3AsString( &distNum ) << " [m GeV]"
    << "\ndist = " << dist << " [m]"
    << "\ndetO_rest_unit = " << utils::print::Vec3AsString( &detO_rest_unit )
    << "\ndetO_unit = " << utils::print::Vec3AsString( &detO_unit )
    << "\np4NHL_good_unit = " << utils::print::Vec3AsString( &p4NHL_good_unit );

  // but we don't care about that. We just want to obtain a proxy for betaNHL in lab frame.
  // Then we can use the dk2nu-style formula modified for betaNHL!

  /* 
   * it is NOT sufficient to boost this into lab frame! 
   * Only a small portion of the CM decays can possibly reach the detector, 
   * imposing a constraint on the allowed directions of p4NHL_rest. 
   * You will miscalculate the NHL energy if you just Boost here. 
   */
  // explicitly calculate the boost correction to lab-frame energy
  // in a dk2nu-like fashion. See bsim::CalcEnuWgt()
  //double betaNHL = p4NHL_rest.P() / p4NHL_rest.E();
  double betaNHL = p4NHL_good.P() / p4NHL_good.E();
  double boost_correction = 0.0;
  double costh_pardet = 0.0;
  if( parentMomentum > 0.0 ){
    costh_pardet = ( decay_pdpx * detO.X() +
		     decay_pdpy * detO.Y() +
		     decay_pdpz * detO.Z() ) / ( parentMomentum * detO.Mag() );
    if( costh_pardet < -1.0 ) costh_pardet = -1.0;
    if( costh_pardet > 1.0 ) costh_pardet = 1.0;
    // assume boost is on z' direction where z' = parent momentum direction, subbing betaMag ==> betaMag * costh_pardet
    //boost_correction = gamma * ( 1.0 + betaNHL * betaMag * costh_pardet );
    if( std::abs( costh_pardet ) >= 0.9 && boost_correction * p4NHL_rest.E() > p4NHL_rest.M() ){
      boost_correction = 1.0 / ( gamma * ( 1.0 - betaMag * betaNHL * costh_pardet ) );
    } else {
      boost_correction = p4NHL_good.E() / p4NHL_rest_good.E();
    }
  }

  LOG( "NHL", pDEBUG )
    << "\ndetO            = " << utils::print::Vec3AsString( &detO )
    << "\ndetO_rest_unit  = " << utils::print::Vec3AsString( &detO_rest_unit )
    << "\np4NHL_rest_good = " << utils::print::P4AsString( &p4NHL_rest_good )
    << "\np4NHL_good      = " << utils::print::P4AsString( &p4NHL_good )
    << "\nbetaNHL         = " << betaNHL
    << "\nboost_corr_two  = " << boost_correction_two
    << "\nboost_corr_one  = " << boost_correction
    << "\ncostheta = " << costh_pardet << ", betaMag = " << betaMag
    << "\nratio           = " << boost_correction_two / boost_correction;

  assert( boost_correction > 0.0 && boost_correction_two > 0.0 );

  // so now we have the random decay. Direction = parent direction, energy = what we calculated
  double ENHL = p4NHL_rest.E() * boost_correction;
  double MNHL = p4NHL_rest.M();
  double PNHL = std::sqrt( ENHL * ENHL - MNHL * MNHL );
  TVector3 pdu = ( p4par.Vect() ).Unit();
  TLorentzVector p4NHL_rand( PNHL * pdu.X(), PNHL * pdu.Y(), PNHL * pdu.Z(), ENHL );

  // find random point in BBox and force momentum to point to that point
  // first, separation in beam frame
  TVector3 fRVec_beam = this->PointToRandomPointInBBox( detO_beam );
  // rotate it and get unit
  TVector3 fRVec_unit = (this->ApplyUserRotation( fRVec_beam )).Unit();
  // force NHL to point along this direction
  TLorentzVector p4NHL( p4NHL_rand.P() * fRVec_unit.X(),
			p4NHL_rand.P() * fRVec_unit.Y(),
			p4NHL_rand.P() * fRVec_unit.Z(),
			p4NHL_rand.E() );

  TVector3 pNHL_beam = this->ApplyUserRotation( p4NHL.Vect(), true );
  TLorentzVector p4NHL_beam( pNHL_beam.X(), pNHL_beam.Y(), pNHL_beam.Z(), p4NHL.E() );
  
  LOG( "NHL", pDEBUG )
    << "\nRandom:  " << utils::print::P4AsString( &p4NHL_rand )
    << "\nPointed: " << utils::print::P4AsString( &p4NHL )
    << "\nRest:    " << utils::print::P4AsString( &p4NHL_rest );
  
  // calculate acceptance correction
  // first, get minimum and maximum deviation from parent momentum to hit detector in degrees

  /*
  double zm = ( isParentOnAxis ) ? 0.0 : this->GetAngDeviation( p4par, detO, false );
  double zp = this->GetAngDeviation( p4par, detO, true );
  */

  if( !fGeoManager )
    fGeoManager = TGeoManager::Import(fGeomFile.c_str());
  double zm = 0.0, zp = 0.0;
  if( fUseBeamMomentum )
    this->GetAngDeviation( p4par_beam, detO_beam, fGeoManager, zm, zp );
  else
    this->GetAngDeviation( p4par, detO, fGeoManager, zm, zp );

  if( zm == -999.9 && zp == 999.9 ){
    this->FillNonsense( iEntry, gnmf ); return;
  }

  if( isParentOnAxis ){ 
    double tzm = zm, tzp = zp;
    zm = 0.0;
    zp = (tzp - tzm)/2.0; // 1/2 * angular opening
  }

  double accCorr = this->CalculateAcceptanceCorrection( p4par, p4NHL_rest, decay_necm, zm, zp );
  
  // also have to factor in boost correction itself... that's same as energy boost correction squared
  // which means a true acceptance of...
  double acceptance = acc_saa * boost_correction * boost_correction * accCorr;

  // finally, a delay calculation
  // if SMv arrives at t=0, then NHL arrives at t = c * ( 1 - beta_NHL ) / L

  double detDist = std::sqrt( detO.X() * detO.X() +
			      detO.Y() * detO.Y() +
			      detO.Z() * detO.Z() ); // m
  const double kSpeedOfLightNs = units::kSpeedOfLight * units::ns / units::s; // m / ns
  double delay = detDist / kSpeedOfLightNs * ( 1.0 / betaNHL - 1.0 );
  delay *= units::ns / units::s;

  LOG( "NHL", pDEBUG )
    << "\ndetDist = " << detDist << " [m]"
    << "\nbetaNHL = " << betaNHL
    << "\ndelay = " << delay << " [ns]";
  
  // write 4-position all this happens at
  TLorentzVector x4NHL_beam( decay_vx, decay_vy, decay_vz, delay ); // in cm, ns
  TLorentzVector x4NHL( -detO.X(), -detO.Y(), -detO.Z(), delay ); // in m, ns
  TLorentzVector x4NHL_cm( units::m / units::cm * ( -detO.X() ),
			   units::m / units::cm * ( -detO.Y() ),
			   units::m / units::cm * ( -detO.Z() ), delay ); // in cm, ns

  LOG( "NHL", pDEBUG )
    << "\nnimpwt                 = " << decay_nimpwt
    << "\nacc_saa                = " << acc_saa
    << "\nboost_correction       = " << boost_correction
    << "\naccCorr                = " << accCorr;

  // fill all the GNuMIFlux stuff
  // comments as seeon on https://www.hep.utexas.edu/~zarko/wwwgnumi/v19/v19/output_gnumi.html

  gnmf->pcodes = 1;                          ///< converted to PDG
  gnmf->units = 0;                           ///< cm
  
  int typeMod = 1;
  if( !fIsMajorana ) typeMod = ( decay_ptype > 0 ) ? 1 : -1;
  gnmf->fgPdgC = typeMod * kPdgNHL;          ///< PDG code

  gnmf->fgXYWgt = acceptance;                ///< geometrical * collimation correction

  gnmf->fgP4 = p4NHL_beam;                   ///< generated 4-momentum, beam coord [GeV]
  gnmf->fgX4 = x4NHL_beam;                   ///< generated 4-position, beam coord [cm]
  gnmf->fgP4User = p4NHL;                    ///< generated 4-momentum, user coord [GeV]
  gnmf->fgX4User = x4NHL_cm;                 ///< generated 4-position, user coord [cm]

  gnmf->evtno    = iEntry;                   ///< Event number (proton on target) 
                                                 // RETHERE which is it?
  gnmf->ndxdz    = p4NHL.Px() / p4NHL.Pz();  ///< Neutrino direction slope for a random decay
  gnmf->ndydz    = p4NHL.Py() / p4NHL.Pz();  ///< See above
  gnmf->npz      = p4NHL.Pz();               ///< Neutrino momentum [GeV] along z direction (beam axis)
  gnmf->nenergy  = p4NHL.E();                ///< Neutrino energy [GeV] for a random decay
  gnmf->ndxdznea = -9999.9;                  ///< Neutrino direction slope for a decay forced to ND
  gnmf->ndydznea = -9999.9;                  ///< See above
  gnmf->nenergyn = boost_correction;         ///< Neutrino energy for decay forced to ND // now houses ratio of boost calcs
  gnmf->nwtnear  = accCorr;                  ///< weight for decay forced to ND / now acceptance correction
  gnmf->ndxdzfar = -9999.9;                  ///< Same as ND but FD
  gnmf->ndydzfar = -9999.9;                  ///< See above
  gnmf->nenergyf = -9999.9;                  ///< See above
  gnmf->nwtfar   = -9999.9;                  ///< See above
  gnmf->norig    = potnum;                  ///< Obsolete...
  
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

  gnmf->ppmedium = fLepPdg;                  ///< Co-produced lepton PDG code
  gnmf->ptype = decay_ptype;                 ///< Parent GEANT code particle ID converted to PDG

  gnmf->ppvx = fLPx;                         ///< Used here to store NHL polarisation x in NHL rest
  gnmf->ppvy = fLPy;                         ///< Used here to store NHL polarisation y in NHL rest
  gnmf->ppvz = fLPz;                         ///< Used here to store NHL polarisation z in NHL rest

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
    << "Finished MakeTupleFluxEntry()";
  
}
//----------------------------------------------------------------------------
void NHLFluxCreator::FillNonsense( int iEntry, flux::GNuMIFluxPassThroughInfo * gnmf ) const
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
void NHLFluxCreator::OpenFluxInput( std::string finpath ) const
{
  //if( std::strcmp( finpath.c_str(), fCurrPath.c_str() ) == 0 ) return;

  iCurrEntry = fFirstEntry;
  fCurrPath = finpath;
  finpath.append("/");

  LOG( "NHL", pDEBUG )
    << "Getting flux input from finpath = " << finpath.c_str();

  // recurse over files in this directory and add to chain
  ctree = new TChain( "dkTree" ); // "dkRootTree"
  cmeta = new TChain( "dkMeta" ); // "dkRootMeta"

  TSystemDirectory dir( finpath.c_str(), finpath.c_str() );
  TList * files = dir.GetListOfFiles(); int nFiles = 0;
  assert( files );
  files->Sort(); // for sanity!

  TSystemFile * file;
  TString fname;
  TIter next(files);
  
  while( (file=( TSystemFile * ) next()) ){
    fname = file->GetName();
    if( !file->IsDirectory() ){
      TString fullpath = TString( finpath.c_str() ) + fname;
      nFiles++;
      ctree->Add( fullpath );
      cmeta->Add( fullpath );
    }
  }

  if( !ctree ){ LOG( "NHL", pFATAL ) << "Could not open flux tree!"; }
  if( !cmeta ){ LOG( "NHL", pFATAL ) << "Could not open meta tree!"; }
  assert( ctree && cmeta );

  const int nEntriesInMeta = cmeta->GetEntries();
  int nEntries = ctree->GetEntries();

  fNEntries = nEntries;

  LOG( "NHL", pDEBUG )
    << "\nThere were " << nEntriesInMeta << " entries in meta with " << nEntries << " total nus"
    << "\n got from " << nFiles << " files";

  delete file;
  delete files;
}
//----------------------------------------------------------------------------
void NHLFluxCreator::InitialiseTree() const
{
  LOG("NHL", pDEBUG) << "Tree initialised";
  potnum = 0.0;
  decay_ptype = 0;
  decay_vx = 0.0; decay_vy = 0.0; decay_vz = 0.0;
  decay_pdpx = 0.0; decay_pdpy = 0.0; decay_pdpz = 0.0;
  decay_nimpwt = 0.0;
  
  ctree->SetBranchAddress( "potnum",       &potnum       );
  ctree->SetBranchAddress( "decay_ptype",  &decay_ptype  );
  ctree->SetBranchAddress( "decay_vx",     &decay_vx     );
  ctree->SetBranchAddress( "decay_vy",     &decay_vy     );
  ctree->SetBranchAddress( "decay_vz",     &decay_vz     );
  ctree->SetBranchAddress( "decay_pdpx",   &decay_pdpx   );
  ctree->SetBranchAddress( "decay_pdpy",   &decay_pdpy   );
  ctree->SetBranchAddress( "decay_pdpz",   &decay_pdpz   );
  ctree->SetBranchAddress( "decay_necm",   &decay_necm   );
  ctree->SetBranchAddress( "decay_nimpwt", &decay_nimpwt );
}
//----------------------------------------------------------------------------
void NHLFluxCreator::InitialiseMeta() const
{ 
  job = 0;
  pots = 0.0;

  cmeta->SetBranchAddress( "job",  &job  );
  cmeta->SetBranchAddress( "pots", &pots );
}
//----------------------------------------------------------------------------
void NHLFluxCreator::ReadBRs() const
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
std::map< NHLProd_t, double > NHLFluxCreator::GetProductionProbs( int parPDG ) const
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
  
  double M    = fMass;
  double Ue42 = fU4l2s.at(0);
  double Um42 = fU4l2s.at(1);
  double Ut42 = fU4l2s.at(2);

  // now get parent mass
  //double mP = PDGLibrary::Instance()->Find( std::abs( parPDG ) )->Mass();

  // also, construct an NHLBRFunctions * object to handle the scalings.
  const Algorithm * algBRFunc = AlgFactory::Instance()->GetAlgorithm("genie::NHL::NHLBRFunctions", "Default");
  const NHLBRFunctions * BRFunc = dynamic_cast< const NHLBRFunctions * >( algBRFunc );
  
  // first get pure kinematic part of the BRs
  double KScale[4] = { -1.0, -1.0, -1.0, -1.0 }, mixScale[4] = { -1.0, -1.0, -1.0, -1.0 };
  double totalMix = 0.0;
  switch( std::abs( parPDG ) ){
  case genie::kPdgMuon:
    KScale[0] = BRFunc->KScale_Global( kNHLProdMuon3Numu, M );
    KScale[1] = BRFunc->KScale_Global( kNHLProdMuon3Nue, M ); // same, convenience for later
    KScale[2] = BRFunc->KScale_Global( kNHLProdMuon3Nutau, M ); // same, convenience for later
    mixScale[0] = 1.0 * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = 1.0 * Ue42 * KScale[1]; totalMix += mixScale[1];
    mixScale[2] = 1.0 * Ut42 * KScale[2]; totalMix += mixScale[2];

    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdMuon3Numu,  mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdMuon3Nue,   mixScale[1] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdMuon3Nutau, mixScale[2] / totalMix } ) );

    // it can happen that NHL is not coupled to the only kinematically available channel.
    // Return bogus map if that happens
    if( totalMix <= 0.0 ){
      dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdNull, -999.9 } ) );
      dynamicScores_muon = dynScores;
      return dynScores;
    }

    dynamicScores_muon = dynScores;
    break;
  case genie::kPdgKP:
    KScale[0] = BRFunc->KScale_Global( kNHLProdKaon2Muon, M );
    KScale[1] = BRFunc->KScale_Global( kNHLProdKaon2Electron, M );
    KScale[2] = BRFunc->KScale_Global( kNHLProdKaon3Muon, M );
    KScale[3] = BRFunc->KScale_Global( kNHLProdKaon3Electron, M );
    mixScale[0] = BR_K2mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_K2e  * Ue42 * KScale[1]; totalMix += mixScale[1];
    mixScale[2] = BR_K3mu * Um42 * KScale[2]; totalMix += mixScale[2];
    mixScale[3] = BR_K3e  * Ue42 * KScale[3]; totalMix += mixScale[3];

    // it can happen that NHL is not coupled to the only kinematically available channel.
    // Return bogus map if that happens
    if( totalMix <= 0.0 ){
      dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdNull, -999.9 } ) );
      dynamicScores_pion = dynScores;
      return dynScores;
    }

    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdKaon2Muon,     mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdKaon2Electron, mixScale[1] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdKaon3Muon,     mixScale[2] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdKaon3Electron, mixScale[3] / totalMix } ) );

    dynamicScores_kaon = dynScores;
    break;
  case genie::kPdgPiP:

    KScale[0] = BRFunc->KScale_Global( kNHLProdPion2Muon, M );
    KScale[1] = BRFunc->KScale_Global( kNHLProdPion2Electron, M );
    mixScale[0] = BR_pi2mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_pi2e  * Ue42 * KScale[1]; totalMix += mixScale[1];

    // it can happen that NHL is not coupled to the only kinematically available channel.
    // Return bogus map if that happens
    if( totalMix <= 0.0 ){
      dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdNull, -999.9 } ) );
      dynamicScores_pion = dynScores;
      return dynScores;
    }

    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdPion2Muon,     mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdPion2Electron, mixScale[1] / totalMix } ) );

    dynamicScores_pion = dynScores;
    break;
  case genie::kPdgK0L:

    KScale[0] = BRFunc->KScale_Global( kNHLProdNeuk3Muon, M );
    KScale[1] = BRFunc->KScale_Global( kNHLProdNeuk3Electron, M );
    mixScale[0] = BR_K03mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_K03e  * Ue42 * KScale[1]; totalMix += mixScale[1];
    
    // it can happen that NHL is not coupled to the only kinematically available channel.
    // Return bogus map if that happens
    if( totalMix <= 0.0 ){
      dynScores.insert( std::pair< NHLProd_t, double >( { kNHLProdNull, -999.9 } ) );
      dynamicScores_neuk = dynScores;
      return dynScores;
    }

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
TLorentzVector NHLFluxCreator::NHLEnergy( NHLProd_t nhldm, TLorentzVector p4par ) const
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
  TGenPhaseSpace fPhaseSpaceGenerator;
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
  
  // Grab 0th entry energy and return that
  int idp = 0; TLorentzVector p4NHL, p4NHL_rest;
  // search for the charged lepton in this stack, get the 4-vector in parent rest frame
  TLorentzVector p4Lep, p4Lep_NHLRest;

  for(std::vector<int>::const_iterator pdg_iter = decayList.begin(); pdg_iter != decayList.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);

     if( std::abs( pdgc ) == kPdgNHL ) p4NHL = *p4fin;
     if( std::abs( pdgc ) == kPdgElectron || 
	 std::abs( pdgc ) == kPdgMuon ||
	 std::abs( pdgc ) == kPdgTau ){
       p4Lep = *p4fin;
       fLepPdg = pdgc;
     }
     idp++;
  }
  
  if( doPol ){
    // boost this to NHL rest frame.
    TVector3 boostNHL = p4NHL.BoostVector();
    p4Lep_NHLRest = p4Lep; p4Lep_NHLRest.Boost( boostNHL );
    // write this
    fLPx = ( fixPol ) ? fFixedPolarisation.at(0) : p4Lep.Px() / p4Lep.P();
    fLPy = ( fixPol ) ? fFixedPolarisation.at(1) : p4Lep.Py() / p4Lep.P();
    fLPz = ( fixPol ) ? fFixedPolarisation.at(2) : p4Lep.Pz() / p4Lep.P();
  } else {
    fLPx = 0.0;
    fLPy = 0.0;
    fLPz = 0.0;
  }
  
  
  return p4NHL; // rest frame momentum!
}
//----------------------------------------------------------------------------
TVector3 NHLFluxCreator::PointToRandomPointInBBox( TVector3 detO_beam ) const
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
double NHLFluxCreator::GetAngDeviation( TLorentzVector p4par, TVector3 detO, bool seekingMax ) const
{
  TVector3 ppar = p4par.Vect(); assert( ppar.Mag() > 0.0 );
  TVector3 pparUnit = ppar.Unit();
  // let face be planar and perpendicular to vector Q
  // RETHERE: assuming Q = ( 0, 0, 1 ) == face perpendicular to z
  double Qx = 0.0, Qy = 0.0, Qz = 1.0;
  // plane: Qx . (x-xC) + Qy . (y-yC) + Qz . (z-zC) = 0
  // line: r(t) - r(D) = t * ppar
  // V0 \in plane && line
  // ==> Qx * ( x(t) - x(C) ) + Qy * ( y(t) - y(C) ) + Qz * ( z(t) - z(C) ) = 0
  // ==> t = ( \vec{Q} \cdot \vec{detO} ) / ( \vec{Q} \cdot \vec{ppar} )
  double nterm = Qx * detO.X() + Qy * detO.Y() + Qz * detO.Z();
  double dterm = Qx * ppar.X() + Qy * ppar.Y() + Qz * ppar.Z();
  double t = nterm / dterm;
  double x_incp = t * pparUnit.X(), y_incp = t * pparUnit.Y(), z_incp = t * pparUnit.Z();

  // sweep over plane perp to ppar, passing through centre, and calc intersection point
  // special case: parent is perfectly on axis so hits detector centre
  TVector3 IPdev( detO.X() - x_incp, detO.Y() - y_incp, detO.Z() - z_incp );
  bool parentHitsCentre = ( IPdev.Mag() < controls::kASmallNum );

  // RETHERE: see assumption about Q
  // to fix probably with a rotation of fLx, fLy by Euler angles onto Q-plane?
  // line: r(t) - r(incp) = t * IPdev
  // RETHERE: assume square face
  double ttx = ( IPdev.X() != 0.0 ) ? fLx / std::abs( IPdev.X() ) : 99999.9;
  double tty = ( IPdev.Y() != 0.0 ) ? fLy / std::abs( IPdev.Y() ) : 99999.9;
  double tt = std::max( ttx, tty ); // this defines how much the least sweep goes
  TVector3 atilde( tt * IPdev.X(), tt * IPdev.Y(), tt * IPdev.Z() );
  
  // fLT = d0 + dtilde = IPdev.Mag()
  // dist = dtilde

  double fLT = IPdev.Mag();
  double dist = atilde.Mag();
  
  assert( fLT > 0.0 );
  double detRadius = std::max( fLx, fLy ) / 2.0;

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
    TVector3 rVec = IPdev.Unit();
    // two IP with det. Closer(farther) has distance detRadius -(+) d( IP, detO )
    // actually, RETHERE: if IPdev endpoint lies inside detector have to go other way
    double dh = fLT + dist, dl = fLT - dist;
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
void NHLFluxCreator::GetAngDeviation( TLorentzVector p4par, TVector3 detO, TGeoManager * gm, double &zm, double &zp ) const
{
  // implementation of GetAngDeviation that uses ROOT geometry. More robust than analytical geom
  // (fewer assumptions about detector position)
  
  LOG( "NHL", pDEBUG )
    << "Entering GetAngDeviation( TLorentzVector, TVector3, TGeoManager * )...";

  assert( gm );
  TVector3 ppar = p4par.Vect(); assert( ppar.Mag() > 0.0 );
  TVector3 pparUnit = ppar.Unit();

  const double sx1 = TMath::Sin(fBx1), cx1 = TMath::Cos(fBx1);
  const double sz = TMath::Sin(fBz), cz = TMath::Cos(fBz);
  const double sx2 = TMath::Sin(fBx2), cx2 = TMath::Cos(fBx2);

  const double xun[3] = { cz, -cx1*sz, sx1*sz };
  const double yun[3] = { sz*cx2, cx1*cz*cx2 - sx1*sx2, -sx1*cz*cx2 - cx1*sx2 };
  const double zun[3] = { sz*sx2, cx1*cz*sx2 + sx1*cx2, -sx1*cz*sx2 + cx1*cx2 };
  
  //--- each unit vector{x,y,z}un has a theta and a phi. Extract these.
  
  double thx = TMath::ACos( xun[2] );
  double phx = 0.0;
  if( TMath::Sin( thx ) != 0.0 ){
    double xtmp = TMath::ACos( xun[0] / TMath::Sin( thx ) );
    phx = ( xun[1] > 0.0 ) ? xtmp : -1.0 * xtmp;
  }
  thx *= TMath::RadToDeg(); phx *= TMath::RadToDeg();
  
  double thy = TMath::ACos( yun[2] );
  double phy = 0.0;
  if( TMath::Sin( thy ) != 0.0 ){
    double ytmp = TMath::ACos( yun[0] / TMath::Sin( thy ) );
    phy = ( yun[1] > 0.0 ) ? ytmp : -1.0 * ytmp;
  }
  thy *= TMath::RadToDeg(); phy *= TMath::RadToDeg();

  double thz = TMath::ACos( zun[2] );
  double phz = 0.0;
  if( TMath::Sin( thz ) != 0.0 ){
    double ztmp = TMath::ACos( zun[0] / TMath::Sin( thz ) );
    phz = ( zun[1] > 0.0 ) ? ztmp : -1.0 * ztmp;
  }
  thz *= TMath::RadToDeg(); phz *= TMath::RadToDeg();

  /* plane is perpendicular to unit vector zun. Find where the momentum intersects this.
   * Formula for plane passing through point (x0, y0, z0) and perpendicular to vector (a,b,c) is
   * a.(x-x0) + b.(y-y0) + c.(z-z0) = 0 (1)
   * Formula for line passing through point (i0, j0, k0) and parallel to vector (u,v,w) is
   * (i,j,k)(t) = (i0, j0, k0) + t*(u,v,w) (2)
   * For x0, y0, z0, i0, j0, k0, a, b, c, u, v, w known, subbing in (2) in (1) gives
   * a.(u*t - detO.X()) + b.(v*t - detO.Y()) + c.(w*t - detO.Z()) = 0 (3)
   * ==> (\vec{ppar}.\vec{zun})*t = \vec{zun}.\vec{detO}
   * ==> t = \vec{zun}.\vec{detO} / \vec{ppar}.\vec{zun} [cm]
   */

  TVector3 detO_cm( (detO.X() + fDetOffset.at(0)) * units::m / units::cm, 
		    (detO.Y() + fDetOffset.at(1)) * units::m / units::cm,
		    (detO.Z() + fDetOffset.at(2)) * units::m / units::cm );
  double inProd = zun[0] * detO_cm.X() + zun[1] * detO_cm.Y() + zun[2] * detO_cm.Z(); // cm
  assert( pparUnit.X() * zun[0] + pparUnit.Y() * zun[1] + pparUnit.Z() * zun[2] != 0.0 );
  inProd /= ( pparUnit.X() * zun[0] + pparUnit.Y() * zun[1] + pparUnit.Z() * zun[2] );
  const double startPoint[3] = { fDx * units::m / units::cm + inProd * pparUnit.X(),
				 fDy * units::m / units::cm + inProd * pparUnit.Y(),
				 fDz * units::m / units::cm + inProd * pparUnit.Z() }; // cm
  const double sweepVect[3] = { (fCx + fDetOffset.at(0)) * units::m / units::cm - startPoint[0],
				(fCy + fDetOffset.at(1)) * units::m / units::cm - startPoint[1],
				(fCz + fDetOffset.at(2)) * units::m / units::cm - startPoint[2] }; // cm
  const double swvMag = std::sqrt( sweepVect[0]*sweepVect[0] + sweepVect[1]*sweepVect[1] + sweepVect[2]*sweepVect[2] ); assert( swvMag > 0.0 );

  // Note the geometry manager works in the *detector frame*. Transform to that.
  TVector3 detStartPoint( startPoint[0] - fCx * units::m / units::cm,
			  startPoint[1] - fCy * units::m / units::cm,
			  startPoint[2] - fCz * units::m / units::cm );
  TVector3 detSweepVect( sweepVect[0] / swvMag, sweepVect[1] / swvMag, sweepVect[2] / swvMag );

  TVector3 detPpar = ppar;

  TVector3 dumori(0.0, 0.0, 0.0); // tgt-hall frame origin is 0
  TVector3 detori( (fDetOffset.at(0)) * units::m / units::cm,
		   (fDetOffset.at(1)) * units::m / units::cm,
		   (fDetOffset.at(2)) * units::m / units::cm ); // for rotations of the detector

  LOG( "NHL", pDEBUG )
    << "\nStartPoint = " << utils::print::Vec3AsString( &detStartPoint )
    << "\nSweepVect  = " << utils::print::Vec3AsString( &detSweepVect )
    << "\nparent p3  = " << utils::print::Vec3AsString( &detPpar );
  
  detStartPoint = this->ApplyUserRotation( detStartPoint, detori, fDetRotation, false ); // passive transformation
  detSweepVect = this->ApplyUserRotation( detSweepVect, dumori, fDetRotation, true );
  detPpar = this->ApplyUserRotation( detPpar, dumori, fDetRotation, true );
  
  // now sweep along sweepVect until we hit either side of the detector. 
  // This will give us two points in space

  LOG( "NHL", pDEBUG )
    << "\npparUnit = ( " << pparUnit.X() << ", " << pparUnit.Y() << ", " << pparUnit.Z() << " )"
    << "\ndetO = ( " << detO_cm.X() << ", " << detO_cm.Y() << ", " << detO_cm.Z() << " ) [cm]"
    << "\nDecay point = ( " << fDx * units::m / units::cm << ", " << fDy * units::m / units::cm << ", " << fDz * units::m / units::cm << " ) [cm]"
    << "\nDetector centre at ( " << fCx * units::m / units::cm << ", " << fCy * units::m / units::cm << ", " << fCz * units::m / units::cm << " ) [cm]"
    << "\nDetector offset is ( " << fDetOffset.at(0) * units::m / units::cm << ", " << fDetOffset.at(1) * units::m / units::cm << ", " << fDetOffset.at(2) * units::m / units::cm << " ) [cm]"
    << "\ninProd = " << inProd << " [cm]"
    << "\nDetector rotation in extrinsic x-z-x = ( " << fBx1 << ", " << fBz << ", " << fBx2 << " ) [rad]"
    << "\nCurrent point (in user coords) = ( " << detStartPoint.X() << ", " << detStartPoint.Y() << ", " << detStartPoint.Z() << " ) [cm]"
    << "\nCurrent direction = ( " << detSweepVect.X() << ", " << detSweepVect.Y() << ", " << detSweepVect.Z() << " )";

  //gm->SetCurrentPoint( detStartPoint.X(), detStartPoint.Y(), detStartPoint.Z() );
  gGeoManager->SetCurrentPoint( detStartPoint.X(), detStartPoint.Y(), detStartPoint.Z() );
  //gm->SetCurrentDirection( detSweepVect.X(), detSweepVect.Y(), detSweepVect.Z() );
  gGeoManager->SetCurrentDirection( detSweepVect.X(), detSweepVect.Y(), detSweepVect.Z() );

  //TGeoNode * nextNode = gm->FindNextBoundaryAndStep( );
  TGeoNode * nextNode = gGeoManager->FindNextBoundaryAndStep( );
  
  if( nextNode == NULL ){
    LOG( "NHL", pWARN )
      << "Calculation failed. Returning nonsense.";
    zm = -999.9; zp = 999.9; return;
  }

  // sometimes the TGeoManager likes to hit the BBox and call this an entry point. Step forward again.
  //const double * tmpPoint = gm->GetCurrentPoint();
  const double * tmpPoint = gGeoManager->GetCurrentPoint();

  LOG( "NHL", pDEBUG )
    << "\ntmpPoint = ( " << tmpPoint[0] << ", " << tmpPoint[1] << ", " << tmpPoint[2] << " ) [ local, cm ]";

  /*
  if( std::abs(tmpPoint[0]) == fLxR/2.0 * units::m / units::cm ||
      std::abs(tmpPoint[1]) == fLyR/2.0 * units::m / units::cm ||
      std::abs(tmpPoint[2]) == fLzR/2.0 * units::m / units::cm )
    nextNode = gm->FindNextBoundaryAndStep();
  */
  
  assert( nextNode != NULL );

  //double minusPoint[3] = { (gm->GetCurrentPoint())[0], (gm->GetCurrentPoint())[1], (gm->GetCurrentPoint())[2] };
  double minusPoint[3] = { (gGeoManager->GetCurrentPoint())[0], (gGeoManager->GetCurrentPoint())[1], (gGeoManager->GetCurrentPoint())[2] };
  
  int bdIdx = 0; const int bdIdxMax = 1e+4;
  //double plusPoint[3] = { (gm->GetCurrentPoint())[0], (gm->GetCurrentPoint())[1], (gm->GetCurrentPoint())[2] };
  double plusPoint[3] = { (gGeoManager->GetCurrentPoint())[0], (gGeoManager->GetCurrentPoint())[1], (gGeoManager->GetCurrentPoint())[2] };
  //while( gm->FindNextBoundaryAndStep() && bdIdx < bdIdxMax ){
  while( gGeoManager->FindNextBoundaryAndStep() && bdIdx < bdIdxMax ){
    bdIdx++;
    if( bdIdx % 100  == 0 )
      LOG( "NHL", pDEBUG ) << "bdIdx = " << bdIdx;
    // explicit check against runaway exit
    //if( std::abs( (gm->GetCurrentPoint())[0] ) != fLxR/2.0 * units::m / units::cm &&
    //std::abs( (gm->GetCurrentPoint())[1] ) != fLyR/2.0 * units::m / units::cm &&
    //	std::abs( (gm->GetCurrentPoint())[2] ) != fLzR/2.0 * units::m / units::cm ){
    // plusPoint[0] = (gm->GetCurrentPoint())[0];
    // plusPoint[1] = (gm->GetCurrentPoint())[1];
    // plusPoint[2] = (gm->GetCurrentPoint())[2] ;
    if( std::abs( (gGeoManager->GetCurrentPoint())[0] ) != fLxR/2.0 * units::m / units::cm &&
	std::abs( (gGeoManager->GetCurrentPoint())[1] ) != fLyR/2.0 * units::m / units::cm &&
    	std::abs( (gGeoManager->GetCurrentPoint())[2] ) != fLzR/2.0 * units::m / units::cm ){
      plusPoint[0] = (gGeoManager->GetCurrentPoint())[0];
      plusPoint[1] = (gGeoManager->GetCurrentPoint())[1];
      plusPoint[2] = (gGeoManager->GetCurrentPoint())[2]; 
    }
  }

  // with 3 points and 1 vector we calculate the angles.
  // Points: D(decay), E(entry), X(exit) [all in local, cm]. Vector: detPpar [local, GeV/GeV]
  // angles are <DE, detPpar> and <DX, detPpar>

  // now obtain the angles themselves and return in deg.
  TVector3 decVec( (fDx - fCx) * units::m / units::cm, 
		   (fDy - fCy) * units::m / units::cm, 
		   (fDz - fCz) * units::m / units::cm );
  decVec = this->ApplyUserRotation( decVec, detori, fDetRotation, false ); // passive transformation
  TVector3 minusVec( minusPoint[0] - decVec.X(), minusPoint[1] - decVec.Y(), minusPoint[2] - decVec[2] );
  TVector3 plusVec( plusPoint[0] - decVec.X(), plusPoint[1] - decVec.Y(), plusPoint[2] - decVec[2] );

  double minusNum = detPpar.X() * minusVec.X() + detPpar.Y() * minusVec.Y() + detPpar.Z() * minusVec.Z();
  double minusDen = detPpar.Mag() * minusVec.Mag(); assert( minusDen > 0.0 );

  zm = TMath::ACos( minusNum / minusDen ) * TMath::RadToDeg();

  double plusNum = detPpar.X() * plusVec.X() + detPpar.Y() * plusVec.Y() + detPpar.Z() * plusVec.Z();
  double plusDen = detPpar.Mag() * plusVec.Mag(); assert( plusDen > 0.0 );

  zp = TMath::ACos( plusNum / plusDen ) * TMath::RadToDeg();

  if( zm > zp ){
    double tmpzp = zp;
    zp = zm;
    zm = zp;
  }

  LOG( "NHL", pDEBUG )
    << "\nIn DETECTOR coordinates:"
    << "\nentered at ( " << minusPoint[0] << ", " << minusPoint[1] << ", " << minusPoint[2] << " ) [cm]"
    << "\nexited  at ( " << plusPoint[0] << ", " << plusPoint[1] << ", " << plusPoint[2] << " ) [cm]"
    << "\nstarted at ( " << decVec.X() << ", " << decVec.Y() << ", " << decVec.Z() << " ) [cm]"
    << "\nmomentum   ( " << detPpar.X() << ", " << detPpar.Y() << ", " << detPpar.Z() << " )"
    << "\nmeaning zm = " << zm << ", zp = " << zp << " [deg]";
}
//----------------------------------------------------------------------------
double NHLFluxCreator::CalculateAcceptanceCorrection( TLorentzVector p4par, TLorentzVector p4NHL,
						      double SMECM, double zm, double zp ) const
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
  if( zp == zm ) return 1.0;

  double M = p4NHL.M();
  if( M == 0.0 ) return 1.0;

  TF1 * fNHL = ( TF1* ) gROOT->GetListOfFunctions()->FindObject( "fNHL" );
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

  TF1 * fSMv = ( TF1* ) gROOT->GetListOfFunctions()->FindObject( "fSMv" );
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
    }
    if( range2 < 1.0e-6 || range1 / range2 > 10.0 ) return 1.0; // sometimes this happens, gotta bail!
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

  LOG( "NHL", pINFO )
    << "\n range1 = " << range1 << " / range2 = " << range2;

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
int NHLFluxCreator::SelectMass( const double mN ) const
{
  const massHyp_t massesHyp[] = {
    kLight0Hyp,  kLight1Hyp,  kLight2Hyp,  kLight3Hyp,  kLight4Hyp,
    kLight5Hyp,  kLight6Hyp,  kLight7Hyp,  kLight8Hyp,  kLight9Hyp,
    kLightAHyp,  kLightBHyp,  kLightCHyp,  kLightDHyp,  kLightEHyp,
    kLightFHyp,  kLightGHyp,  kLightHHyp,  kLightIHyp,  kLightJHyp,
    kMedium0Hyp, kMedium1Hyp, kMedium2Hyp, kMedium3Hyp, kMedium4Hyp,
    kMedium5Hyp, kMedium6Hyp, kMedium7Hyp, kMedium8Hyp, kMedium9Hyp,
    kHeavy0Hyp,  kHeavy1Hyp,  kHeavy2Hyp,  kHeavy3Hyp,  kHeavy4Hyp,
    kHeavy5Hyp,  kHeavy6Hyp,  kHeavy7Hyp,  kHeavy8Hyp,  kHeavy9Hyp,
    kHeavyAHyp,  kHeavyBHyp,  kHeavyCHyp,  kHeavyDHyp,  kHeavyEHyp,
    kHeavyFHyp,  kHeavyGHyp,  kHeavyHHyp,  kHeavyIHyp,  kHeavyJHyp };
  const int nMasses = sizeof(massesHyp)/sizeof(massesHyp[0]) - 1;
  
  // because masses are kept in a *map*, gotta build array of the second elements!
  double masses[ nMasses + 1 ];
    for( int i = 0; i <= nMasses; i++ ){
      massHyp_t thisHyp = massesHyp[i];
      auto pos = massHypMap.find( thisHyp );
      masses[i] = pos->second;
      LOG("NHL", pINFO) 
	<< "At position " << i << " the mass hypothesis is " << masses[i];
    }
    
    assert( mN >= 0.0 );
    
    int mp = -1; double fmN = 0.0; // mass of the selected masspoint
    if( mN >= masses[ nMasses ] ){ mp = nMasses; fmN = masses[ nMasses ]; }
    while( masses[ mp + 1 ] < mN && mp < nMasses ){ mp++; } // decide interval
    
    // generally decide mass + point by closest endpoint in interval
    const double dLeft  = std::abs( mN - masses[ mp ] );
    const double dRight = masses[ mp + 1 ] - mN;
    
    LOG("NHL", pDEBUG) <<
      "Stats:" <<
      "\n Input mass: " << mN <<
      "\n Choice interval: [ " << masses[mp] << ", " << masses[mp+1] << " ] " <<
      "\n Left and right distance: " << dLeft << ", " << dRight <<
      "\n Chosen point: " << ( ( dLeft < dRight ) ? mp : mp + 1 ) << " ( " <<
      ( ( dLeft < dRight ) ? "LEFT )" : "RIGHT )" );
    
    fmN = ( dLeft < dRight ) ? masses[ mp ] : masses[ mp + 1 ];
    mp  = ( dLeft < dRight ) ? mp : mp + 1;
    
    return mp;
}
//----------------------------------------------------------------------------
std::string NHLFluxCreator::SelectFile( std::string fin, const double mN ) const
{
  std::string filePath = fin;
  filePath.append( "/FHC/" ); // RETHERE make this configurable!!!
  filePath.append( "EqualCouplings/" );
  const int mp = SelectMass( mN );
  filePath.append( Form( "mp%02d", mp ) );
  filePath.append( ".root" );
  return filePath;
}
//----------------------------------------------------------------------------
TH1F * NHLFluxCreator::GetFluxHist1F( std::string fin, int masspoint, bool isParticle ) const
{
  TFile * f = TFile::Open( fin.c_str() );
  assert( f );

  std::string histName( Form("hHNLFluxCenterAcc_%d", masspoint) );
  if( !isParticle ) histName.append("_bar");
  TObject * histObj = f->Get( histName.c_str() );
  
  LOG( "NHL", pDEBUG )
    << "Getting flux from histo with name " << histObj->GetName()
    << " and title " << histObj->GetTitle();
  
  TH1F * histPtr = dynamic_cast< TH1F* >( histObj );

  return histPtr;
}
//----------------------------------------------------------------------------
TH3D * NHLFluxCreator::GetFluxHist3D( std::string fin, std::string dirName, std::string hName ) const
{
  TFile * f = TFile::Open( fin.c_str() );
  
  TDirectory * baseDir = f->GetDirectory( "" );
  TDirectory * deepDir = baseDir->GetDirectory( dirName.c_str() );
  assert( deepDir );
  assert( deepDir->GetListOfKeys()->Contains( hName.c_str() ) );
  
  TH3D * histPtr = dynamic_cast< TH3D* >( deepDir->Get( hName.c_str() ) );
  return histPtr;
}
//----------------------------------------------------------------------------
std::vector< double > * NHLFluxCreator::GenerateVtx3X( TH3D * prodVtxHist ) const
{
  if( prodVtxHist->GetEntries() == 0 ){ // no production vtx? ok, set NHL to be 1km upstream
    std::vector< double > * vtxDir = new std::vector< double >();
    vtxDir->emplace_back( 0.0 );
    vtxDir->emplace_back( 0.0 );
    vtxDir->emplace_back( -100000.0 );
    vtxDir->emplace_back( 0.0 );
    return vtxDir;
  }

  double ux = 0.0, uy = 0.0, uz = 0.0;
  prodVtxHist->GetRandom3( ux, uy, uz );
  
  std::vector< double > * vtxDir = new std::vector< double >();
  vtxDir->emplace_back( ux );
  vtxDir->emplace_back( uy );
  vtxDir->emplace_back( uz );
  vtxDir->emplace_back( 0.0 ); // RETHERE: add some timing information?
  
  return vtxDir;
}
//----------------------------------------------------------------------------
void NHLFluxCreator::MakeBBox() const
{
  LOG( "NHL", pWARN )
    << "WARNING: This is a dummy (==unit-side) bounding box centred at config-given point";

  LOG( "NHL", pDEBUG )
    << "\nfCx, fCy, fCz = " << fCx << ", " << fCy << ", " << fCz
    << "\nfAx1, fAz, fAx2 = " << fAx1 << ", " << fAz << ", " << fAx2;


  fLx = 1.0; fLy = 1.0; fLz = 1.0;
}
//----------------------------------------------------------------------------
TVector3 NHLFluxCreator::ApplyUserRotation( TVector3 vec, bool doBackwards ) const
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
TVector3 NHLFluxCreator::ApplyUserRotation( TVector3 vec, TVector3 oriVec, std::vector<double> rotVec, bool doBackwards ) const
{
  double vx = vec.X(), vy = vec.Y(), vz = vec.Z();
  double ox = oriVec.X(), oy = oriVec.Y(), oz = oriVec.Z();

  LOG( "NHL", pDEBUG )
    << "\n original vec : " << utils::print::Vec3AsString( &vec )
    << "\n origin vec : " << utils::print::Vec3AsString( &oriVec );
  
  vx -= ox; vy -= oy; vz -= oz; // make this rotation about detector origin

  assert( rotVec.size() == 3 ); // want 3 Euler angles, otherwise this is unphysical.
  double Ax2 = ( doBackwards ) ? -rotVec.at(2) : rotVec.at(2);
  double Az  = ( doBackwards ) ? -rotVec.at(1) : rotVec.at(1);
  double Ax1 = ( doBackwards ) ? -rotVec.at(0) : rotVec.at(0);

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

  // back to beam frame
  vx += ox; vy += oy; vz += oz;
  TVector3 nvec( vx, vy, vz );
  LOG( "NHL", pDEBUG ) << "\nFinal vector: " << utils::print::Vec3AsString( &nvec );
  return nvec;
}
//____________________________________________________________________________
double NHLFluxCreator::CalculateDetectorAcceptanceSAA( TVector3 detO ) const
{
  // sang is solid-angle / 4pi
  double rad = std::sqrt( detO.X() * detO.X() + detO.Y() * detO.Y() + detO.Z() * detO.Z() );
  double sang = 1.0 - TMath::Cos( TMath::ATan( kRDET / rad ) ); sang *= 0.5;
  return sang;
}
//----------------------------------------------------------------------------
double NHLFluxCreator::CalculateAreaNormalisation()
{
  // for now this is just a square of length kRDET
  // returns 1 / area
  return 1.0 / ( kRDET * kRDET );
}
//----------------------------------------------------------------------------
void NHLFluxCreator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//----------------------------------------------------------------------------
void NHLFluxCreator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//----------------------------------------------------------------------------
void NHLFluxCreator::LoadConfig(void)
{
  if( fIsConfigLoaded ) return;

  LOG("NHL", pDEBUG)
    << "Loading flux-creation parameters from file...";

  this->GetParam( "NHL-Mass", fMass );
  this->GetParamVect( "NHL-LeptonMixing", fU4l2s );
  this->GetParam( "NHL-Majorana", fIsMajorana );
  
  this->GetParamVect( "Beam2User_T", fB2UTranslation );
  this->GetParamVect( "Beam2User_R", fB2URotation );
  this->GetParamVect( "DetCentre_User", fDetOffset );
  this->GetParamVect( "Beam2Det_R", fDetRotation );

  this->GetParamVect( "ParentPOTScalings", fScales );
  this->GetParam( "IncludePolarisation", doPol );
  this->GetParam( "FixPolarisationDirection", fixPol );
  this->GetParamVect( "NHL-PolDir", fFixedPolarisation );

  this->GetParam( "IsParentOnAxis", isParentOnAxis );
  this->GetParam( "UseBeamMomentum", fUseBeamMomentum );

  fCx = fB2UTranslation.at(0);
  fCy = fB2UTranslation.at(1);
  fCz = fB2UTranslation.at(2);
  
  fAx1 = fB2URotation.at(0);
  fAz  = fB2URotation.at(1);
  fAx2 = fB2URotation.at(2);

  fBx1 = fDetRotation.at(0);
  fBz  = fDetRotation.at(1);
  fBx2 = fDetRotation.at(2);

  POTScaleWeight = 1.0;
  if( utils::nhl::IsProdKinematicallyAllowed( kNHLProdMuon3Nue ) ) 
    POTScaleWeight = fScales[0]; // all POT contribute
  else if( utils::nhl::IsProdKinematicallyAllowed( kNHLProdPion2Muon ) ||
	   utils::nhl::IsProdKinematicallyAllowed( kNHLProdPion2Electron ) ) 
    POTScaleWeight = fScales[1]; // muons do not contribute
  else if( utils::nhl::IsProdKinematicallyAllowed( kNHLProdNeuk3Muon ) ||
	   utils::nhl::IsProdKinematicallyAllowed( kNHLProdNeuk3Electron ) )
    POTScaleWeight = fScales[2]; // only kaons contribute
  else if( utils::nhl::IsProdKinematicallyAllowed( kNHLProdKaon2Muon ) || 
	   utils::nhl::IsProdKinematicallyAllowed( kNHLProdKaon2Electron ) ||
	   utils::nhl::IsProdKinematicallyAllowed( kNHLProdKaon3Muon ) ||
	   utils::nhl::IsProdKinematicallyAllowed( kNHLProdKaon3Electron ) )
    POTScaleWeight = fScales[3]; // only charged kaons contribute

  LOG( "NHL", pDEBUG )
    << "Read the following parameters :"
    << "\n Mass = " << fMass
    << "\n couplings = " << fU4l2s.at(0) << " : " << fU4l2s.at(1) << " : " << fU4l2s.at(2)
    << "\n translation = " << fB2UTranslation.at(0) << ", " << fB2UTranslation.at(1) << ", " << fB2UTranslation.at(2)
    << "\n rotation = " << fB2URotation.at(0) << ", " << fB2URotation.at(1) << ", " << fB2URotation.at(2)
    << "\n isParentOnAxis = " << isParentOnAxis
    << "\n POTScaleWeight = " << POTScaleWeight;

  fIsConfigLoaded = true;
}
//____________________________________________________________________________
void NHLFluxCreator::SetGeomFile( string geomfile ) const
{
  fGeomFile = geomfile;
}
//____________________________________________________________________________
void NHLFluxCreator::ImportBoundingBox( TGeoBBox * box ) const
{
  fLxR = 2.0 * box->GetDX() * units::cm / units::m;
  fLyR = 2.0 * box->GetDY() * units::cm / units::m;
  fLzR = 2.0 * box->GetDZ() * units::cm / units::m;
}
//____________________________________________________________________________
