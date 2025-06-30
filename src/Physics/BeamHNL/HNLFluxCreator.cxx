//____________________________________________________________________________
/*
  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  
  Author: John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford
*/
//____________________________________________________________________________

#include "Physics/BeamHNL/HNLFluxCreator.h"

using namespace genie;
using namespace genie::hnl;
using namespace genie::hnl::enums;

//----------------------------------------------------------------------------
FluxCreator::FluxCreator() :
  FluxRecordVisitorI("genie::hnl::FluxCreator")
{

}
//----------------------------------------------------------------------------
FluxCreator::FluxCreator(string name) :
  FluxRecordVisitorI(name)
{

}
//----------------------------------------------------------------------------
FluxCreator::FluxCreator(string name, string config) :
  FluxRecordVisitorI(name, config)
{

}
//----------------------------------------------------------------------------
FluxCreator::~FluxCreator()
{

}
//----------------------------------------------------------------------------
void FluxCreator::ProcessEventRecord(GHepRecord * evrec) const
{
  // Adds the inital state HNL at the event record.
  // Also assigns the production vertex to evrec (this will be overwritten by subsequent modules)
  // Also adds (acceptance*nimpwt)^(-1) component of weight
  
  this->SetCurrentEntry( evrec->XSec() );
  
  if( !fIsUsingRootGeom ){
    this->SetUsingRootGeom(true); // must always be true
    
    gGeoManager = TGeoManager::Import( fGeomFile.c_str() );
    
    TGeoVolume * top_volume = gGeoManager->GetTopVolume();
    assert( top_volume );
    TGeoShape * ts = top_volume->GetShape();
    TGeoBBox * box = (TGeoBBox *) ts;
    
    this->ImportBoundingBox( box );  
  }

  if( fUsingDk2nu ){

    if( iCurrEntry < fFirstEntry ) iCurrEntry = fFirstEntry;
    
    if( iCurrEntry >= fFirstEntry ) {
      
      if( iCurrEntry == fFirstEntry ){
	FluxContainer * pfGnmf = new FluxContainer();
	fGnmf = *pfGnmf;
	delete pfGnmf;
      }
      
      fGnmf = this->MakeTupleFluxEntry( iCurrEntry, fCurrPath );
      
      if( std::abs(fGnmf.pdg) == genie::kPdgHNL ){ // only add particle if parent is valid

	LOG( "HNL", pDEBUG ) << fGnmf;
	
	double invAccWeight = fGnmf.nimpwt * fGnmf.acceptance;
	evrec->SetWeight( evrec->Weight() / invAccWeight );
	
	// scale by how many POT it takes to make the appropriate parent
	/*
	 * To incorporate populations of parents, we take the cumulative multiplicity
	 * i.e. HNL light enough to be made by every parent get scaled by 
	 * n1 = \sigma(p + target) / \sigma(p + target ; parent-producing)
	 * For HNL that are heavier than a muon, we don't take muons into account. So
	 * we up the scaling to incorporate their dropping out as
	 * n2 = \sigma(p + target) / \sigma(p + target ; parent-producing ; no muon) - etc.
	 */
	evrec->SetWeight( evrec->Weight() * POTScaleWeight );
	
	// set prod-vertex in cm, ns, NEAR coords
	TVector3 xVtxNear = fGnmf.startPoint; // in m, NEAR
	double tVtx = fGnmf.delay; // ns
	TLorentzVector pVtx( xVtxNear.X() * units::m / units::cm,
			     xVtxNear.Y() * units::m / units::cm,
			     xVtxNear.Z() * units::m / units::cm, tVtx ); // cm ns NEAR
	evrec->SetVertex( pVtx ); // HNL production vertex. NOT where HNL decays to visible FS.
	
	// construct Particle(0). Don't worry about daughter links at this stage.
	// this must be in USER coords
	TLorentzVector probeP4 = fGnmf.p4User; // USER
	GHepParticle ptHNL( fGnmf.pdg, kIStInitialState, -1, -1, -1, -1, probeP4, pVtx );
	evrec->AddParticle( ptHNL );
      }

      if( fGnmf.acceptance >= 0.0 ){
	// set some event information where subsequent events can see it.
	// Use the x4 position of the HNL. First, ensure the Vertex() is correctly set.
	TLorentzVector * vx4 = evrec->Particle(0)->GetX4();
	evrec->SetVertex( *vx4 );
	TLorentzVector tmpx4( fLPx, fLPy, fGnmf.parPdg, fGnmf.lepPdg );
	if( fLPz >= 0.0 ) evrec->SetXSec( 1.0 );
	else evrec->SetXSec( -1.0 );
	evrec->Particle(0)->SetPosition( tmpx4 );
	
      } // if( fGnmf.fgXYWgt >= 0 )
    } // if( iCurrEntry > fFirstEntry )
  } else {
    LOG( "HNL", pFATAL )
      << "No input dk2nu flux detected. Cannot proceed.";
    exit(1);
  }
}
//----------------------------------------------------------------------------
void FluxCreator::SetInputFluxPath(std::string finpath) const
{
  LOG( "HNL", pDEBUG ) << "Setting input path to " << finpath;
  fCurrPath = finpath;
}
//----------------------------------------------------------------------------
int FluxCreator::GetNFluxEntries() const
{
  if( fNEntries <= 0 ){
    this->OpenFluxInput( fCurrPath );
  }
  return fNEntries;
}
//----------------------------------------------------------------------------
void FluxCreator::SetCurrentEntry( int iCurr ) const
{
  iCurrEntry = iCurr;
}
//----------------------------------------------------------------------------
FluxContainer FluxCreator::RetrieveFluxInfo() const
{
  return fGnmf;
}
//----------------------------------------------------------------------------
void FluxCreator::SetUsingRootGeom( bool IsUsingRootGeom ) const
{
  fIsUsingRootGeom = IsUsingRootGeom;
}
//____________________________________________________________________________
FluxContainer FluxCreator::MakeTupleFluxEntry( int iEntry, std::string finpath ) const
{
  // This method creates 1 HNL from the flux info and saves the information
  // Essentially, it replaces a SMv with an HNL
  FluxContainer * tmpGnmf = new FluxContainer();
  FluxContainer gnmf = *tmpGnmf;
  delete tmpGnmf;

  // Open flux input and initialise trees
  if( iEntry == fFirstEntry ){
    this->OpenFluxInput( finpath );
    this->InitialiseTree();
    this->InitialiseMeta();
    if(fDoingOldFluxCalc) this->MakeBBox();
  } else if( iEntry < fFirstEntry ){
    this->FillNonsense( iEntry, gnmf );
    return gnmf;
  }

  TVector3 dumori( 0.0, 0.0, 0.0 ); // use to rotate VECTORS
  TVector3 originPoint( -(fCx + fDetOffset.at(0)), 
			-(fCy + fDetOffset.at(1)), 
			-(fCz + fDetOffset.at(2)) ); // use to rotate POINTS. m

  // All these in m
  TVector3 fCvec_beam( fCx, fCy, fCz );
  TVector3 fCvec = this->ApplyUserRotation( fCvec_beam ); // in NEAR coords
  fLepPdg = 0;
  
  ctree->GetEntry(iEntry);

  // explicitly check if there are any allowed decays for this parent
  bool canGoForward = true;
  switch( std::abs( decay_ptype ) ){
  case kPdgPiP:
    canGoForward = 
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdPion2Muon ) ||
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdPion2Electron ); break;
  case kPdgKP:
    canGoForward =
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdKaon2Muon ) || 
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdKaon2Electron ) ||
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdKaon3Muon ) ||
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdKaon3Electron ); break;
  case kPdgMuon:
    canGoForward =
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdMuon3Nue ); break; // SM nus are massless so doesn't matter
  case kPdgK0L:
    canGoForward =
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdNeuk3Muon ) ||
      utils::hnl::IsProdKinematicallyAllowed( kHNLProdNeuk3Electron ); break;
  }

  if( !canGoForward ){
    this->FillNonsense( iEntry, gnmf ); return gnmf; 
  }
    
  // turn cm to m and make origin wrt detector 
  fDx = decay_vx * units::cm / units::m; // BEAM, m
  fDy = decay_vy * units::cm / units::m;
  fDz = decay_vz * units::cm / units::m;

  if( !fSupplyingBEAM ){
    TVector3 tmpVec( fDx, fDy, fDz ); // NEAR
    tmpVec = this->ApplyUserRotation( tmpVec, false ); // BEAM
    fDx = tmpVec.X(); fDy = tmpVec.Y(); fDz = tmpVec.Z();
  }

  TVector3 fDvec( fDx, fDy, fDz ); // in BEAM coords
  TVector3 fDvec_beam = this->ApplyUserRotation( fDvec, true ); // in NEAR coords

  TVector3 fDvec_user( fDvec_beam.X() - fCx, fDvec_beam.Y() - fCy, fDvec_beam.Z() - fCz ); // in USER coords
  fDvec_user = this->ApplyUserRotation( fDvec_user, originPoint, fDetRotation, false );

  LOG( "HNL", pDEBUG )
    << "\nIn BEAM coords, fDvec = " << utils::print::Vec3AsString( &fDvec )
    << "\nIn NEAR coords, fDvec = " << utils::print::Vec3AsString( &fDvec_beam );

  TVector3 detO_beam( fCvec_beam.X() - fDvec_beam.X(),
		      fCvec_beam.Y() - fDvec_beam.Y(),
		      fCvec_beam.Z() - fDvec_beam.Z() ); // separation in NEAR coords
  TVector3 detO( fCvec.X() - fDvec.X(),
		 fCvec.Y() - fDvec.Y(),
		 fCvec.Z() - fDvec.Z() ); // separation in BEAM coords
  TVector3 detO_user( detO_beam.X(), detO_beam.Y(), detO_beam.Z() );

  detO_user = this->ApplyUserRotation( detO_user, dumori, fDetRotation, false ); // tgt-hall --> det
  
  double acc_saa = this->CalculateDetectorAcceptanceSAA( detO_user );
  
  // set parent mass

  double dpdpx = decay_pdpx, dpdpy = decay_pdpy, dpdpz = decay_pdpz; // BEAM GeV
  if( !fSupplyingBEAM ){
    TVector3 tmpVec( dpdpx, dpdpy, dpdpz ); // NEAR
    tmpVec = this->ApplyUserRotation( tmpVec, false ); // BEAM
    dpdpx = tmpVec.X(); dpdpy = tmpVec.Y(); dpdpz = tmpVec.Z();
  }

  switch( std::abs( decay_ptype ) ){
  case kPdgPiP: case kPdgKP: case kPdgMuon: case kPdgK0L:
    parentMass = PDGLibrary::Instance()->Find(decay_ptype)->Mass(); break;
  default:
    LOG( "HNL", pERROR ) << "Parent with PDG code " << decay_ptype << " not handled!"
			 << "\n\tProceeding, but results are possibly unphysical.";
    parentMass = PDGLibrary::Instance()->Find(decay_ptype)->Mass(); break;
  }
  parentMomentum = std::sqrt( dpdpx*dpdpx + dpdpy*dpdpy + dpdpz*dpdpz );
  parentEnergy = std::sqrt( parentMass*parentMass + parentMomentum*parentMomentum );

  TLorentzVector p4par = ( isParentOnAxis ) ? 
    TLorentzVector( parentMomentum * (detO.Unit()).X(), 
		    parentMomentum * (detO.Unit()).Y(),
		    parentMomentum * (detO.Unit()).Z(),
		    parentEnergy ) :
    TLorentzVector( dpdpx, dpdpy, dpdpz, parentEnergy ); // in BEAM coords

  TLorentzVector p4par_beam( dpdpx, dpdpy, dpdpz, parentEnergy ); // in BEAM coords
  TVector3 p3par_beam = p4par_beam.Vect();
  TVector3 p3par_near = this->ApplyUserRotation( p3par_beam, true );
  TLorentzVector p4par_near( p3par_near.X(), p3par_near.Y(), p3par_near.Z(), parentEnergy ); // in NEAR coords
  LOG( "HNL", pDEBUG )
    << "\nIn BEAM coords: p3par_beam = " << utils::print::Vec3AsString( &p3par_beam )
    << "\nIn NEAR coords: p3par_near = " << utils::print::Vec3AsString( &p3par_near );
  if( !isParentOnAxis ){
    // rotate p4par to NEAR coordinates
    TVector3 tmpv3 = ApplyUserRotation( p4par.Vect(), true );
    p4par.SetPxPyPzE( tmpv3.Px(), tmpv3.Py(), tmpv3.Pz(), p4par.E() );
  }

  TLorentzVector p4par_user = p4par_near;
  // rotate it to user coords
  TVector3 ppar_user = p4par_user.Vect();
  ppar_user = this->ApplyUserRotation( ppar_user, dumori, fDetRotation, false ); // tgt-hall --> det
  p4par_user.SetPxPyPzE( ppar_user.Px(), ppar_user.Py(), ppar_user.Pz(), p4par_user.E() );

  TVector3 boost_beta = p4par_beam.BoostVector(); // in BEAM coords

  // now calculate which decay channel produces the HNL.
  dynamicScores = this->GetProductionProbs( decay_ptype );
  assert( dynamicScores.size() > 0 );
  
  if( dynamicScores.find( kHNLProdNull ) != dynamicScores.end() ){ // exists kin allowed channel but 0 coupling
    this->FillNonsense( iEntry, gnmf ); return gnmf;
  }
  
  RandomGen * rnd = RandomGen::Instance();
  double score = rnd->RndGen().Uniform( 0.0, 1.0 );
  HNLProd_t prodChan;
  // compare with cumulative prob. If < 1st in map, pick 1st chan. If >= 1st and < (1st+2nd), pick 2nd, etc
  
  unsigned int imap = 0; double s1 = 0.0;
  std::map< HNLProd_t, double >::iterator pdit = dynamicScores.begin();
  while( score >= s1 && pdit != dynamicScores.end() ){
    s1 += (*pdit).second;
    if( parentMass > 0.495 ){
      LOG( "HNL", pDEBUG )
	<< "(*pdit).first = " << utils::hnl::ProdAsString( (*pdit).first )
	<< " : (*pdit).second = " << (*pdit).second;
    }
    if( score >= s1 ){
      imap++; pdit++;
    }
  }
  assert( imap < dynamicScores.size() ); // should have decayed to *some* HNL
  prodChan = (*pdit).first;

  // bookkeep this
  fProdChan = static_cast<int>(prodChan);
  switch( prodChan ){
  case kHNLProdNeuk3Electron: fNuProdChan = 1; fNuPdg = kPdgNuE; break;
  case kHNLProdNeuk3Muon: fNuProdChan = 2; fNuPdg = kPdgNuMu; break;
  case kHNLProdKaon2Electron: fNuProdChan = 4; fNuPdg = kPdgNuE; break;
  case kHNLProdKaon2Muon: fNuProdChan = 3; fNuPdg = kPdgNuMu; break;
  case kHNLProdKaon3Electron: fNuProdChan = 6; fNuPdg = kPdgNuE; break;
  case kHNLProdKaon3Muon: fNuProdChan = 5; fNuPdg = kPdgNuMu; break;
  case kHNLProdMuon3Nue:
  case kHNLProdMuon3Numu:
  case kHNLProdMuon3Nutau: fNuProdChan = 9; fNuPdg = kPdgNuE; break;
  case kHNLProdPion2Electron: fNuProdChan = 8; fNuPdg = kPdgNuE; break;
  case kHNLProdPion2Muon: fNuProdChan = 7; fNuPdg = kPdgNuMu; break;
  default: fNuProdChan = -999; fNuPdg = -999; break;
  }

  LOG( "HNL", pDEBUG )
    << "Selected channel: " << utils::hnl::ProdAsString( prodChan );
  
  // decay channel specified, now time to make kinematics
  TLorentzVector p4HNL_rest = HNLEnergy( prodChan, p4par ); 
  // this is a random direction rest-frame HNL. 
  fECM = p4HNL_rest.E();

  // we will now boost detO into rest frame, force rest to point to the new direction, boost the result, and compare the boost corrections
  double boost_correction_two = 0.0;
  
  // 17-Jun-22: Notice the time component needs to be nonzero to get this to work!

  // first guess: betaHNL ~= 1 . Do the Lorentz boosts knocking betaHNL downwards until we hit det centre
  double betaMag = boost_beta.Mag();
  double gamma   = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
  double betaLab = 1.0; // first guess

  // now make a TLorentzVector in lab frame to boost back to rest. 
  double timeBit = detO.Mag(); // / units::kSpeedOfLight ; // s
  TLorentzVector detO_4v( detO.X(), detO.Y(), detO.Z(), timeBit ); detO_4v.Boost( -boost_beta ); // BEAM with BEAM
  TVector3 detO_rest_unit = (detO_4v.Vect()).Unit();
  TLorentzVector p4HNL_rest_good( p4HNL_rest.P() * detO_rest_unit.X(),
				  p4HNL_rest.P() * detO_rest_unit.Y(),
				  p4HNL_rest.P() * detO_rest_unit.Z(),
				  p4HNL_rest.E() );

  double pLep_rest = std::sqrt( fLPx*fLPx + fLPy*fLPy + fLPz*fLPz );
  TLorentzVector p4Lep_rest_good( -1.0 * pLep_rest * detO_rest_unit.X(),
				  -1.0 * pLep_rest * detO_rest_unit.Y(),
				  -1.0 * pLep_rest * detO_rest_unit.Z(),
				  fLPE );

  // boost HNL into lab frame!

  TLorentzVector p4HNL_good = p4HNL_rest_good;
  p4HNL_good.Boost( boost_beta );
  boost_correction_two = p4HNL_good.E() / p4HNL_rest.E();

  TVector3 detO_unit = detO.Unit(); // BEAM

  TVector3 p4HNL_good_vect = p4HNL_good.Vect();
  TVector3 p4HNL_good_unit = p4HNL_good_vect.Unit();

  // now calculate how far away from target point we are.
  // dist = || detO x p4HNL || / || p4HNL_good || where x == cross product
  TVector3 distNum = detO.Cross( p4HNL_good_unit );
  double dist = distNum.Mag(); // m

  double prevDist = 2.0 * dist;

  while( betaLab > 0.0 && ( dist > 1.0e-3 && dist < prevDist  &&
	 std::abs(dist - prevDist) > 1.0e-3 * prevDist ) ){ // 1mm tolerance

    // that didn't work. Knock betaLab down a little bit and try again.
    prevDist = dist;
    betaLab -= 1.0e-4;
    timeBit = detO.Mag() / ( betaLab );
    detO_4v.SetXYZT( detO.X(), detO.Y(), detO.Z(), timeBit );
    detO_4v.Boost( -boost_beta );
    detO_rest_unit = (detO_4v.Vect()).Unit();
    p4HNL_rest_good.SetPxPyPzE( p4HNL_rest.P() * detO_rest_unit.X(),
				p4HNL_rest.P() * detO_rest_unit.Y(),
				p4HNL_rest.P() * detO_rest_unit.Z(),
				p4HNL_rest.E() );

    // boost into lab frame
    p4HNL_good = p4HNL_rest_good;
    p4HNL_good.Boost( boost_beta );
    
    detO_unit = detO.Unit();
    p4HNL_good_vect = p4HNL_good.Vect();
    p4HNL_good_unit = p4HNL_good_vect.Unit();

    distNum = detO.Cross( p4HNL_good_unit );
    dist = distNum.Mag(); // m
  }

  // but we don't care about that. We just want to obtain a proxy for betaHNL in lab frame.
  // Then we can use the dk2nu-style formula modified for betaHNL!

  /* 
   * it is NOT sufficient to boost this into lab frame! 
   * Only a small portion of the CM decays can possibly reach the detector, 
   * imposing a constraint on the allowed directions of p4HNL_rest. 
   * You will miscalculate the HNL energy if you just Boost here. 
   */
  // explicitly calculate the boost correction to lab-frame energy
  // in a dk2nu-like fashion. See bsim::CalcEnuWgt()
  //double betaHNL = p4HNL_rest.P() / p4HNL_rest.E();
  double betaHNL = p4HNL_good.P() / p4HNL_good.E();
  double costh_pardet = 0.0;
  double boost_correction = 0.0;
  if( parentMomentum > 0.0 ){
    costh_pardet = ( p4par_beam.X() * detO.X() +
		     p4par_beam.Y() * detO.Y() +
		     p4par_beam.Z() * detO.Z() ) / ( parentMomentum * detO.Mag() );
    if( costh_pardet < -1.0 ) costh_pardet = -1.0;
    if( costh_pardet > 1.0 ) costh_pardet = 1.0;
    boost_correction = 1.0 / ( gamma * ( 1.0 - betaMag * betaHNL * costh_pardet ) );
    // assume boost is on z' direction where z' = parent momentum direction, subbing betaMag ==> betaMag * costh_pardet
    if( true && boost_correction * p4HNL_rest.E() > p4HNL_rest.M() ) {
      boost_correction = 1.0 / ( gamma * ( 1.0 - betaMag * betaHNL * costh_pardet ) );
    } else {
      boost_correction = p4HNL_good.E() / p4HNL_rest_good.E();
    }
  }

  assert( boost_correction > 0.0 && boost_correction_two > 0.0 );

  // so now we have the random decay. Direction = parent direction, energy = what we calculated
  double EHNL = p4HNL_rest.E() * boost_correction;
  double MHNL = p4HNL_rest.M();
  double PHNL = std::sqrt( EHNL * EHNL - MHNL * MHNL );
  TVector3 pdu = ( p4par.Vect() ).Unit(); // in BEAM coords
  TLorentzVector p4HNL_rand( PHNL * pdu.X(), PHNL * pdu.Y(), PHNL * pdu.Z(), EHNL );

  // find random point in BBox and force momentum to point to that point
  
  double FDx = fDvec_beam.X();
  double FDy = fDvec_beam.Y();
  double FDz = fDvec_beam.Z();

  TVector3 absolutePoint = this->PointToRandomPointInBBox( ); // in NEAR coords, m
  TVector3 fRVec_beam( absolutePoint.X() - FDx, absolutePoint.Y() - FDy, absolutePoint.Z() - FDz ); // NEAR, m
  // rotate it and get unit
  TVector3 fRVec_unit = (this->ApplyUserRotation( fRVec_beam )).Unit(); // BEAM, m/m
  TVector3 fRVec_actualBeam = this->ApplyUserRotation( fRVec_beam ); // BEAM, m
  TVector3 fRVec_user = this->ApplyUserRotation( fRVec_beam, dumori, fDetRotation, false ); // USER m
  

  TVector3 rRVec_near( fRVec_beam.X() + FDx, fRVec_beam.Y() + FDy, fRVec_beam.Z() + FDz );
  TVector3 rRVec_beam = this->ApplyUserRotation( rRVec_near );
  TVector3 rRVec_user = this->ApplyUserRotation( rRVec_near, dumori, fDetRotation, false );
  /*
  rRVec_user.SetXYZ( rRVec_user.X() + (fCx + fDetOffset.at(0)),
		     rRVec_user.Y() + (fCy + fDetOffset.at(1)),
		     rRVec_user.Z() + (fCz + fDetOffset.at(2)) );
  */
  rRVec_user.SetXYZ( rRVec_user.X() + (fCx),
		     rRVec_user.Y() + (fCy),
		     rRVec_user.Z() + (fCz) );

  // force HNL to point along this direction
  TLorentzVector p4HNL( p4HNL_rand.P() * fRVec_unit.X(),
			p4HNL_rand.P() * fRVec_unit.Y(),
			p4HNL_rand.P() * fRVec_unit.Z(),
			p4HNL_rand.E() ); // in BEAM coords

  TVector3 pHNL_near = this->ApplyUserRotation( p4HNL.Vect(), true ); // BEAM --> NEAR coords
  TLorentzVector p4HNL_near( pHNL_near.X(), pHNL_near.Y(), pHNL_near.Z(), p4HNL.E() );

  LOG( "HNL", pDEBUG )
    << "\nRandom:  " << utils::print::P4AsString( &p4HNL_rand )
    << "\nPointed [NEAR]: " << utils::print::P4AsString( &p4HNL_near )
    << "\nRest:    " << utils::print::P4AsString( &p4HNL_rest );

  // update polarisation
  TLorentzVector p4Lep_good = p4Lep_rest_good; // in parent rest frame
  p4Lep_good.Boost( boost_beta ); // in lab frame
  TVector3 boost_beta_HNL = p4HNL_near.BoostVector();
  p4Lep_good.Boost( -boost_beta_HNL ); // in HNL rest frame

  fLPx = ( fixPol ) ? fFixedPolarisation.at(0) : p4Lep_good.Px() / p4Lep_good.P();
  fLPy = ( fixPol ) ? fFixedPolarisation.at(1) : p4Lep_good.Py() / p4Lep_good.P();
  fLPz = ( fixPol ) ? fFixedPolarisation.at(2) : p4Lep_good.Pz() / p4Lep_good.P();
  
  // calculate acceptance correction
  // first, get minimum and maximum deviation from parent momentum to hit detector in degrees
  double zm = 0.0, zp = 0.0;
  if( fIsUsingRootGeom ){
    this->GetAngDeviation( p4par_near, detO_beam, zm, zp ); // using NEAR and NEAR
  } else { // !fIsUsingRootGeom
    zm = ( isParentOnAxis ) ? 0.0 : this->GetAngDeviation( p4par_near, detO_beam, false );
    zp = this->GetAngDeviation( p4par_near, detO_beam, true );
  }

  if( zm == -999.9 && zp == 999.9 ){
    this->FillNonsense( iEntry, gnmf ); return gnmf;
  }

  if( isParentOnAxis ){ 
    double tzm = zm, tzp = zp;
    zm = 0.0;
    zp = (tzp - tzm)/2.0; // 1/2 * angular opening
  }

  fSMECM = decay_necm;
  fZm = zm; fZp = zp;
  double accCorr = this->CalculateAcceptanceCorrection( p4par, p4HNL_rest, decay_necm, zm, zp );
  //if( !fDoingOldFluxCalc ){
  if( fRerollPoints ){
    // if accCorr == 0 then we must ~bail and find the next event.~ 
    //                 We must reroll the point a bunch of times. Then we can skip.
    //                 Ideally we'd be able to tell how much detector lives within the reachable region [zm, zp]
    int iAccFail = 0; const int iAccFailBail = 10;
    while( iAccFail < iAccFailBail && accCorr == 0.0 ){
      LOG( "HNL", pNOTICE )
	<< "Point with separation " << utils::print::Vec3AsString( &fRVec_beam ) << " is unreachable "
	<< "by HNL from parent with momentum " << utils::print::P4AsString( &p4par ) << " !"
	<< "\nRerolling point. This is the " << iAccFail << "th try out of " << iAccFailBail;
      
      // find random point in BBox and force momentum to point to that point
      // first, separation in beam frame
      absolutePoint = this->PointToRandomPointInBBox( ); // always NEAR, m
      /*
      fRVec_beam.SetXYZ( absolutePoint.X() - (fCx + fDetOffset.at(0)),
			 absolutePoint.Y() - (fCy + fDetOffset.at(1)),
			 absolutePoint.Z() - (fCz + fDetOffset.at(2)) );
      */
      fRVec_beam.SetXYZ( absolutePoint.X() - (fCx),
			 absolutePoint.Y() - (fCy),
			 absolutePoint.Z() - (fCz) );
      // rotate it and get unit
      fRVec_unit = (this->ApplyUserRotation( fRVec_beam )).Unit(); // BEAM
      // force HNL to point along this direction
      p4HNL.SetPxPyPzE( p4HNL_rand.P() * fRVec_unit.X(),
			p4HNL_rand.P() * fRVec_unit.Y(),
			p4HNL_rand.P() * fRVec_unit.Z(),
			p4HNL_rand.E() ); // BEAM
      
      pHNL_near = this->ApplyUserRotation( p4HNL.Vect(), true ); // NEAR
      p4HNL_near.SetPxPyPzE( pHNL_near.X(), pHNL_near.Y(), pHNL_near.Z(), p4HNL.E() );
      
      LOG( "HNL", pDEBUG )
	<< "\nRandom:  " << utils::print::P4AsString( &p4HNL_rand )
	<< "\nPointed: " << utils::print::P4AsString( &p4HNL )
	<< "\nRest:    " << utils::print::P4AsString( &p4HNL_rest );
      
      // update polarisation
      p4Lep_good = p4Lep_rest_good; // in parent rest frame
      p4Lep_good.Boost( boost_beta ); // in lab frame
      boost_beta_HNL = p4HNL_near.BoostVector(); // NEAR coords
      p4Lep_good.Boost( -boost_beta_HNL ); // in HNL rest frame
      
      fLPx = ( fixPol ) ? fFixedPolarisation.at(0) : p4Lep_good.Px() / p4Lep_good.P();
      fLPy = ( fixPol ) ? fFixedPolarisation.at(1) : p4Lep_good.Py() / p4Lep_good.P();
      fLPz = ( fixPol ) ? fFixedPolarisation.at(2) : p4Lep_good.Pz() / p4Lep_good.P();
      
      // calculate acceptance correction
      // first, get minimum and maximum deviation from parent momentum to hit detector in degrees
      zm = 0.0; zp = 0.0;
      if( fIsUsingRootGeom ){
	this->GetAngDeviation( p4par_near, detO_beam, zm, zp ); // NEAR and NEAR
      } else { // !fIsUsingRootGeom
	zm = ( isParentOnAxis ) ? 0.0 : this->GetAngDeviation( p4par_near, detO_beam, false );
	zp = this->GetAngDeviation( p4par_near, detO_beam, true );
      }
      
      if( zm == -999.9 && zp == 999.9 ){
	this->FillNonsense( iEntry, gnmf ); return gnmf;
      }
      
      if( isParentOnAxis ){ 
	double tzm = zm, tzp = zp;
	zm = 0.0;
	zp = (tzp - tzm)/2.0; // 1/2 * angular opening
      }
      
      accCorr = this->CalculateAcceptanceCorrection( p4par_near, p4HNL_rest, decay_necm, zm, zp );
      iAccFail++;
    }
  }
  if( accCorr == 0.0 ){ // NOW we can give up and return.
    this->FillNonsense( iEntry, gnmf ); return gnmf;
  }
  
  // also have to factor in boost correction itself... that's same as energy boost correction squared
  // which means a true acceptance of...
  double acceptance = acc_saa * boost_correction * boost_correction * accCorr;

  // finally, a delay calculation
  // if SMv arrives at t=0, then HNL arrives at t = c * ( 1 - beta_HNL ) / L

  double detDist = std::sqrt( detO.X() * detO.X() +
			      detO.Y() * detO.Y() +
			      detO.Z() * detO.Z() ); // m
  const double kSpeedOfLightNs = units::kSpeedOfLight * units::ns / units::s; // m / ns
  double delay = detDist / kSpeedOfLightNs * ( 1.0 / betaHNL - 1.0 );
  delay *= units::ns / units::s;

  /*
  LOG( "HNL", pDEBUG )
    << "\ndetDist = " << detDist << " [m]"
    << "\nbetaHNL = " << betaHNL
    << "\ndelay = " << delay << " [ns]";
  */
  
  // write 4-position all this happens at

  double dxvx = decay_vx, dxvy = decay_vy, dxvz = decay_vz;

  if( !fSupplyingBEAM ){
    TVector3 tmpVec( dxvx, dxvy, dxvz ); // NEAR
    tmpVec = this->ApplyUserRotation( tmpVec, false ); // BEAM
    dxvx = tmpVec.X(); dxvy = tmpVec.Y(); dxvz = tmpVec.Z();
  }

  TLorentzVector x4HNL_beam( dxvx, dxvy, dxvz, delay ); // in cm, ns, BEAM coords
  TVector3 x3HNL_beam = x4HNL_beam.Vect();
  TVector3 x3HNL_near = this->ApplyUserRotation( x3HNL_beam, true );
  TLorentzVector x4HNL_near( x3HNL_near.X(), x3HNL_near.Y(), x3HNL_near.Z(), delay );

  TLorentzVector x4HNL( fDvec_user.X(), fDvec_user.Y(), fDvec_user.Z(), delay ); // in m, ns, USER coords
  TLorentzVector x4HNL_cm( units::m / units::cm * x4HNL.X(),
			   units::m / units::cm * x4HNL.Y(),
			   units::m / units::cm * x4HNL.Z(), delay ); // in cm, ns, USER

  LOG( "HNL", pDEBUG ) << "Filling gnmf...";

  // fill all the flux stuff now!
  
  int typeMod = 1;
  if( !fIsMajorana ) typeMod = ( decay_ptype > 0 ) ? 1 : -1;
  // fix for muons which are backwards...
  int parPdg = decay_ptype;
  if( std::abs(parPdg) == kPdgMuon ) typeMod *= -1;

  gnmf.evtno = iEntry;

  gnmf.pdg = typeMod * kPdgHNL;
  gnmf.parPdg = parPdg;
  gnmf.lepPdg = fLepPdg;
  gnmf.nuPdg = typeMod * fNuPdg;

  gnmf.prodChan = fProdChan;
  gnmf.nuProdChan = fNuProdChan;

  gnmf.startPoint.SetXYZ( fDvec_beam.X(), fDvec_beam.Y(), fDvec_beam.Z() ); // NEAR m
  gnmf.targetPoint.SetXYZ( fTargetPoint.X(), fTargetPoint.Y(), fTargetPoint.Z() ); // NEAR m
  gnmf.startPointUser.SetXYZ( fDvec_user.X() - fCx, fDvec_user.Y() - fCy, fDvec_user.Z() - fCz ); // USER m
  gnmf.targetPointUser.SetXYZ( fTargetPoint.X() - fCx, fTargetPoint.Y() - fCy, fTargetPoint.Z() - fCz ); // USER m
  gnmf.delay = delay; // ns

  gnmf.polz.SetXYZ( fLPx, fLPy, fLPz );

  TLorentzVector p4HNL_user = p4HNL_near;
  // rotate it to user coords
  TVector3 pHNL_user = p4HNL_user.Vect();
  pHNL_user = this->ApplyUserRotation( pHNL_user, dumori, fDetRotation, false ); // tgt-hall --> det
  p4HNL_user.SetPxPyPzE( pHNL_user.Px(), pHNL_user.Py(), pHNL_user.Pz(), p4HNL_user.E() );

  gnmf.p4 = p4HNL_near;
  gnmf.parp4 = p4par_near;
  gnmf.p4User = p4HNL_user;
  gnmf.parp4User = p4par_user;

  gnmf.Ecm = fECM;
  gnmf.nuEcm = fSMECM;

  gnmf.XYWgt = acc_saa;               
  gnmf.boostCorr = p4HNL.E() / fECM;
  gnmf.accCorr = accCorr;
  gnmf.zetaMinus = fZm;
  gnmf.zetaPlus = fZp;
  gnmf.acceptance = acceptance;
  gnmf.nimpwt = decay_nimpwt;
  
  return gnmf;
  
}
//----------------------------------------------------------------------------
void FluxCreator::FillNonsense( int iEntry, FluxContainer & gnmf ) const
{
  gnmf.evtno = iEntry;

  gnmf.pdg = -9999;
  gnmf.parPdg = -9999;
  gnmf.lepPdg = -9999;
  gnmf.nuPdg = -9999;

  gnmf.prodChan = -9999;
  gnmf.nuProdChan = -9999;

  gnmf.startPoint.SetXYZ(-9999.9, -9999.9, -9999.9);
  gnmf.targetPoint.SetXYZ(-9999.9, -9999.9, -9999.9);
  gnmf.startPointUser.SetXYZ(-9999.9, -9999.9, -9999.9);
  gnmf.targetPointUser.SetXYZ(-9999.9, -9999.9, -9999.9);

  gnmf.polz.SetXYZ(-9999.9, -9999.9, -9999.9);

  gnmf.p4.SetPxPyPzE(-9999.9, -9999.9, -9999.9, -9999.9);
  gnmf.parp4.SetPxPyPzE(-9999.9, -9999.9, -9999.9, -9999.9);
  gnmf.p4User.SetPxPyPzE(-9999.9, -9999.9, -9999.9, -9999.9);
  gnmf.parp4User.SetPxPyPzE(-9999.9, -9999.9, -9999.9, -9999.9);

  gnmf.Ecm = -9999.9;
  gnmf.nuEcm = -9999.9;
  gnmf.XYWgt = -9999.9;
  gnmf.boostCorr = -9999.9;
  gnmf.accCorr = -9999.9;
  gnmf.zetaMinus = -9999.9;
  gnmf.zetaPlus = -9999.9;
  gnmf.acceptance = -9999.9;

  gnmf.nimpwt = -9999.9;

  return;
}
//----------------------------------------------------------------------------
void FluxCreator::OpenFluxInput( std::string finpath ) const
{
  //if( std::strcmp( finpath.c_str(), fCurrPath.c_str() ) == 0 ) return;

  iCurrEntry = fFirstEntry;
  fCurrPath = finpath;
  finpath.append("/");

  LOG( "HNL", pDEBUG )
    << "Getting flux input from finpath = " << finpath.c_str();

  // recurse over files in this directory and add to chain
  if(!ctree){
    ctree = new TChain( "dkTree" );
    cmeta = new TChain( "dkMeta" );
  }

  if( fPathLoaded ) return;

  TSystemDirectory dir( finpath.c_str(), finpath.c_str() );
  TList * files = dir.GetListOfFiles(); int nFiles = 0;
  assert( files );
  files->Sort();

  TSystemFile * file;
  TString fname;
  TIter next(files);
  
  while( (file=( TSystemFile * ) next()) && !fPathLoaded ){
    fname = file->GetName();
    if( !file->IsDirectory() ){
      TString fullpath = TString( finpath.c_str() ) + fname;
      nFiles++;
      ctree->Add( fullpath );
      cmeta->Add( fullpath );
    }
  }

  if( !ctree ){ LOG( "HNL", pFATAL ) << "Could not open flux tree!"; }
  if( !cmeta ){ LOG( "HNL", pFATAL ) << "Could not open meta tree!"; }
  assert( ctree && cmeta );

  const int nEntriesInMeta = cmeta->GetEntries();
  int nEntries = ctree->GetEntries();

  fNEntries = nEntries;

  LOG( "HNL", pDEBUG )
    << "\nThere were " << nEntriesInMeta << " entries in meta with " << nEntries << " total nus"
    << "\n got from " << nFiles << " files";

  fPathLoaded = true;

  delete file;
  delete files;
}
//----------------------------------------------------------------------------
void FluxCreator::InitialiseTree() const
{
  potnum = 0.0;
  decay_ptype = 0;
  decay_vx = 0.0; decay_vy = 0.0; decay_vz = 0.0;
  decay_pdpx = 0.0; decay_pdpy = 0.0; decay_pdpz = 0.0;
  decay_nimpwt = 0.0;

  arSize = 0, anArSize = 0, trArSize = 0;
  djob = -9999;
  ppvx = -9999.9, ppvy = -9999.9, ppvz = -9999.9;
  decay_norig = -9999, decay_ndecay = -9999, decay_ntype = -9999;
  decay_ppdxdz = -9999.9, decay_ppdydz = -9999.9, decay_pppz = -9999.9, decay_ppenergy = -9999.9;
  decay_ppmedium = -9999;
  decay_muparpx = -9999.9, decay_muparpy = -9999.9, decay_muparpz = -9999.9, decay_mupare = -9999.9;
  
  tgtexit_tvx = -9999.9, tgtexit_tvy = -9999.9, tgtexit_tvz = -9999.9;
  tgtexit_tpx = -9999.9, tgtexit_tpy = -9999.9, tgtexit_tpz = -9999.9;
  tgtexit_tptype = -9999, tgtexit_tgen = -9999;

  for( int i = 0; i < maxArray; i++ ){
    nuray_px[i]  = -9999.9;
    nuray_py[i]  = -9999.9;
    nuray_pz[i]  = -9999.9;
    nuray_E[i]   = -9999.9;
    nuray_wgt[i] = -9999.9;

    ancestor_pdg[i]     = -9999;
    ancestor_startx[i]  = -9999.9;
    ancestor_starty[i]  = -9999.9;
    ancestor_startz[i]  = -9999.9;
    ancestor_startpx[i] = -9999.9;
    ancestor_startpy[i] = -9999.9;
    ancestor_startpz[i] = -9999.9;
    ancestor_stoppx[i]  = -9999.9;
    ancestor_stoppy[i]  = -9999.9;
    ancestor_stoppz[i]  = -9999.9;
    ancestor_polx[i]    = -9999.9;
    ancestor_poly[i]    = -9999.9;
    ancestor_polz[i]    = -9999.9;
    ancestor_pprodpx[i] = -9999.9;
    ancestor_pprodpy[i] = -9999.9;
    ancestor_pprodpz[i] = -9999.9;
    ancestor_nucleus[i] = -9999;
    
    traj_trkx[i]  = -9999.9;
    traj_trky[i]  = -9999.9;
    traj_trkz[i]  = -9999.9;
    traj_trkpx[i] = -9999.9;
    traj_trkpy[i] = -9999.9;
    traj_trkpz[i] = -9999.9;

    for( int j = 0; j < maxC; j++ ){
      ancestor_proc[i*maxC+j] = 0;
      ancestor_ivol[i*maxC+j] = 0;
      ancestor_imat[i*maxC+j] = 0;
    }
  }
  
  // necessary branches
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
  
  // extra branches
  if( ctree->GetBranch( "job" ) ) ctree->SetBranchAddress( "job", &djob );
  if( ctree->GetBranch( "decay_norig" ) ) ctree->SetBranchAddress( "decay_norig", &decay_norig );
  if( ctree->GetBranch( "decay_ndecay" ) ) ctree->SetBranchAddress( "decay_ndecay", &decay_ndecay );
  if( ctree->GetBranch( "decay_ntype" ) ) ctree->SetBranchAddress( "decay_ntype", &decay_ntype );
  if( ctree->GetBranch( "decay_ppdxdz" ) ) ctree->SetBranchAddress( "decay_ppdxdz", &decay_ppdxdz );
  if( ctree->GetBranch( "decay_ppdydz" ) ) ctree->SetBranchAddress( "decay_ppdydz", &decay_ppdydz );
  if( ctree->GetBranch( "decay_pppz" ) ) ctree->SetBranchAddress( "decay_pppz", &decay_pppz );
  if( ctree->GetBranch( "decay_ppenergy" ) ) ctree->SetBranchAddress( "decay_ppenergy", &decay_ppenergy );
  if( ctree->GetBranch( "decay_ppmedium" ) ) ctree->SetBranchAddress( "decay_ppmedium", &decay_ppmedium );
  if( ctree->GetBranch( "decay_ptype" ) ) ctree->SetBranchAddress( "decay_ptype", &decay_ptype );
  if( ctree->GetBranch( "decay_muparpx" ) ) ctree->SetBranchAddress( "decay_muparpx", &decay_muparpx );
  if( ctree->GetBranch( "decay_muparpy" ) ) ctree->SetBranchAddress( "decay_muparpy", &decay_muparpy );
  if( ctree->GetBranch( "decay_muparpz" ) ) ctree->SetBranchAddress( "decay_muparpz", &decay_muparpz );
  if( ctree->GetBranch( "decay_mupare" ) ) ctree->SetBranchAddress( "decay_mupare", &decay_mupare );
  
  if( ctree->GetBranch( "nuray_size" ) ) ctree->SetBranchAddress( "nuray_size", &arSize );
  if( ctree->GetBranch( "nuray_px" ) ) ctree->SetBranchAddress( "nuray_px", nuray_px );
  if( ctree->GetBranch( "nuray_py" ) ) ctree->SetBranchAddress( "nuray_py", nuray_py );
  if( ctree->GetBranch( "nuray_pz" ) ) ctree->SetBranchAddress( "nuray_pz", nuray_pz );
  if( ctree->GetBranch( "nuray_E" ) ) ctree->SetBranchAddress( "nuray_E", nuray_E );
  if( ctree->GetBranch( "nuray_wgt" ) ) ctree->SetBranchAddress( "nuray_wgt", nuray_wgt );

  if( ctree->GetBranch( "ancestor_size" ) ) ctree->SetBranchAddress( "ancestor_size", &anArSize );
  if( ctree->GetBranch( "ancestor_pdg" ) ) ctree->SetBranchAddress( "ancestor_pdg", ancestor_pdg );
  if( ctree->GetBranch( "ancestor_startx" ) ) ctree->SetBranchAddress( "ancestor_startx", ancestor_startx );
  if( ctree->GetBranch( "ancestor_starty" ) ) ctree->SetBranchAddress( "ancestor_starty", ancestor_starty );
  if( ctree->GetBranch( "ancestor_startz" ) ) ctree->SetBranchAddress( "ancestor_startz", ancestor_startz );
  if( ctree->GetBranch( "ancestor_startpx" ) ) ctree->SetBranchAddress( "ancestor_startpx", ancestor_startpx );
  if( ctree->GetBranch( "ancestor_startpy" ) ) ctree->SetBranchAddress( "ancestor_startpy", ancestor_starty );
  if( ctree->GetBranch( "ancestor_startpz" ) ) ctree->SetBranchAddress( "ancestor_startpz", ancestor_startpz );
  if( ctree->GetBranch( "ancestor_stoppx" ) ) ctree->SetBranchAddress( "ancestor_stoppx", ancestor_stoppx );
  if( ctree->GetBranch( "ancestor_stoppy" ) ) ctree->SetBranchAddress( "ancestor_stoppy", ancestor_stoppy );
  if( ctree->GetBranch( "ancestor_stoppz" ) ) ctree->SetBranchAddress( "ancestor_stoppz", ancestor_stoppz );
  if( ctree->GetBranch( "ancestor_polx" ) ) ctree->SetBranchAddress( "ancestor_polx", ancestor_polx );
  if( ctree->GetBranch( "ancestor_poly" ) ) ctree->SetBranchAddress( "ancestor_poly", ancestor_poly );
  if( ctree->GetBranch( "ancestor_polz" ) ) ctree->SetBranchAddress( "ancestor_polz", ancestor_polz );
  if( ctree->GetBranch( "ancestor_pprodpx" ) ) ctree->SetBranchAddress( "ancestor_pprodpx", ancestor_pprodpx );
  if( ctree->GetBranch( "ancestor_pprodpy" ) ) ctree->SetBranchAddress( "ancestor_pprodpy", ancestor_pprodpy );
  if( ctree->GetBranch( "ancestor_pprodpz" ) ) ctree->SetBranchAddress( "ancestor_pprodpz", ancestor_pprodpz );
  if( ctree->GetBranch( "ancestor_proc" ) ) ctree->SetBranchAddress( "ancestor_proc", ancestor_proc );
  if( ctree->GetBranch( "ancestor_ivol" ) ) ctree->SetBranchAddress( "ancestor_ivol", ancestor_ivol );
  if( ctree->GetBranch( "ancestor_imat" ) ) ctree->SetBranchAddress( "ancestor_imat", ancestor_imat );

  if( ctree->GetBranch( "ppvx" ) ) ctree->SetBranchAddress( "ppvx", &ppvx );
  if( ctree->GetBranch( "ppvy" ) ) ctree->SetBranchAddress( "ppvy", &ppvy );
  if( ctree->GetBranch( "ppvz" ) ) ctree->SetBranchAddress( "ppvz", &ppvz );

  if( ctree->GetBranch( "tgtexit_tvx" ) ) ctree->SetBranchAddress( "tgtexit_tvx", &tgtexit_tvx );
  if( ctree->GetBranch( "tgtexit_tvy" ) ) ctree->SetBranchAddress( "tgtexit_tvy", &tgtexit_tvy );
  if( ctree->GetBranch( "tgtexit_tvz" ) ) ctree->SetBranchAddress( "tgtexit_tvz", &tgtexit_tvz );

  if( ctree->GetBranch( "tgtexit_tpx" ) ) ctree->SetBranchAddress( "tgtexit_tpx", &tgtexit_tpx );
  if( ctree->GetBranch( "tgtexit_tpy" ) ) ctree->SetBranchAddress( "tgtexit_tpy", &tgtexit_tpy );
  if( ctree->GetBranch( "tgtexit_tpz" ) ) ctree->SetBranchAddress( "tgtexit_tpz", &tgtexit_tpz );

  if( ctree->GetBranch( "tgtexit_tptype" ) ) ctree->SetBranchAddress( "tgtexit_tptype", &tgtexit_tptype );
  if( ctree->GetBranch( "tgtexit_tgen" ) ) ctree->SetBranchAddress( "tgtexit_tgen", &tgtexit_tgen );

  if( ctree->GetBranch( "traj_size" ) ) ctree->SetBranchAddress( "traj_size", &trArSize );
  if( ctree->GetBranch( "traj_trkx" ) ) ctree->SetBranchAddress( "traj_trkx", traj_trkx );
  if( ctree->GetBranch( "traj_trky" ) ) ctree->SetBranchAddress( "traj_trky", traj_trky );
  if( ctree->GetBranch( "traj_trkz" ) ) ctree->SetBranchAddress( "traj_trkz", traj_trkz );
  if( ctree->GetBranch( "traj_trkpx" ) ) ctree->SetBranchAddress( "traj_trkpx", traj_trkpx );
  if( ctree->GetBranch( "traj_trkpy" ) ) ctree->SetBranchAddress( "traj_trkpy", traj_trkpy );
  if( ctree->GetBranch( "traj_trkpz" ) ) ctree->SetBranchAddress( "traj_trkpz", traj_trkpz );
  
}
//----------------------------------------------------------------------------
void FluxCreator::InitialiseMeta() const
{ 
  job = 0;
  pots = 0.0;
  
  beam0x = -9999.9;
  beam0y = -9999.9;
  beam0z = -9999.9;
  beamhwidth = -9999.9;
  beamvwidth = -9999.9;
  beamdxdz = -9999.9;
  beamdydz = -9999.9;
  mArSize = 0;

  for( int i = 0; i < maxC; i++ ){
    beamsim[i] = 0;
    physics[i] = 0;
    physcuts[i] = 0;
    tgtcfg[i] = 0;
    horncfg[i] = 0;
    dkvolcfg[i] = 0;
  }

  for( int i = 0; i < maxArray; i++ ){
    location_x[i] = -9999.9;
    location_y[i] = -9999.9;
    location_z[i] = -9999.9;

    for( int j = 0; j < maxC; j++ ){ location_name[i*maxC+j] = 0; }
  }

  // necessary branches
  cmeta->SetBranchAddress( "job",  &job  );
  cmeta->SetBranchAddress( "pots", &pots );

  // extra branches
  if( cmeta->GetBranch( "beamsim" ) ) cmeta->SetBranchAddress( "beamsim", beamsim );
  if( cmeta->GetBranch( "physics" ) ) cmeta->SetBranchAddress( "physics", physics );
  if( cmeta->GetBranch( "physcuts" ) ) cmeta->SetBranchAddress( "physcuts", physcuts );
  if( cmeta->GetBranch( "tgtcfg" ) ) cmeta->SetBranchAddress( "tgtcfg", tgtcfg );
  if( cmeta->GetBranch( "horncfg" ) ) cmeta->SetBranchAddress( "horncfg", horncfg );
  if( cmeta->GetBranch( "dkvolcfg" ) ) cmeta->SetBranchAddress( "dkvolcfg", dkvolcfg );
  if( cmeta->GetBranch( "beam0x" ) ) cmeta->SetBranchAddress( "beam0x", &beam0x );
  if( cmeta->GetBranch( "beam0y" ) ) cmeta->SetBranchAddress( "beam0y", &beam0y );
  if( cmeta->GetBranch( "beam0z" ) ) cmeta->SetBranchAddress( "beam0z", &beam0z );
  if( cmeta->GetBranch( "beamhwidth" ) ) cmeta->SetBranchAddress( "beamhwidth", &beamhwidth );
  if( cmeta->GetBranch( "beamvwidth" ) ) cmeta->SetBranchAddress( "beamvwidth", &beamvwidth );
  if( cmeta->GetBranch( "beamdxdz" ) ) cmeta->SetBranchAddress( "beamdxdz", &beamdxdz );
  if( cmeta->GetBranch( "beamdydz" ) ) cmeta->SetBranchAddress( "beamdydz", &beamdydz );
  if( cmeta->GetBranch( "arSize" ) ) cmeta->SetBranchAddress( "arSize", &mArSize );
  if( cmeta->GetBranch( "location_x" ) ) cmeta->SetBranchAddress( "location_x", location_x );
  if( cmeta->GetBranch( "location_y" ) ) cmeta->SetBranchAddress( "location_y", location_y );
  if( cmeta->GetBranch( "location_z" ) ) cmeta->SetBranchAddress( "location_z", location_z );
  if( cmeta->GetBranch( "location_name" ) ) cmeta->SetBranchAddress( "location_name", location_name );
}
//----------------------------------------------------------------------------
void FluxCreator::ReadBRs() const
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
  BR_K2e   = 1.6e-5; // From PDG 2021
  BR_K3mu  = kaon3muChannel->BranchingRatio();
  BR_K3e   = kaon3elChannel->BranchingRatio();

  BR_K03mu = 2.0 * neuk3muChannel->BranchingRatio(); // one from K0L-->mu+ and one from -->mu-
  BR_K03e  = 2.0 * neuk3elChannel->BranchingRatio();
}
//----------------------------------------------------------------------------
std::map< HNLProd_t, double > FluxCreator::GetProductionProbs( int parPDG ) const
{
  // check if we've calculated scores before
  switch( std::abs( parPDG ) ){
  case kPdgPiP : if( dynamicScores_pion.size() > 0 ) return dynamicScores_pion; break;
  case kPdgKP  : if( dynamicScores_kaon.size() > 0 ) return dynamicScores_kaon; break;
  case kPdgMuon: if( dynamicScores_muon.size() > 0 ) return dynamicScores_muon; break;
  case kPdgK0L : if( dynamicScores_neuk.size() > 0 ) return dynamicScores_neuk; break;
  default: LOG( "HNL", pWARN ) << "Unknown parent. Proceeding, but results may be unphysical"; break;
  }

  std::map< HNLProd_t, double > dynScores;

  // first get branching ratios to SM
  ReadBRs();
  // then get HNL parameter space
  
  double Ue42 = fU4l2s.at(0);
  double Um42 = fU4l2s.at(1);
  double Ut42 = fU4l2s.at(2);

  // also, construct an BRCalculator * object to handle the scalings.
  const Algorithm * algBRCalc = AlgFactory::Instance()->GetAlgorithm("genie::hnl::BRCalculator", "Default");
  const BRCalculator * BRCalc = dynamic_cast< const BRCalculator * >( algBRCalc );
  
  // first get pure kinematic part of the BRs
  double KScale[4] = { -1.0, -1.0, -1.0, -1.0 }, mixScale[4] = { -1.0, -1.0, -1.0, -1.0 };
  double totalMix = 0.0;
  switch( std::abs( parPDG ) ){
  case genie::kPdgMuon:
    KScale[0] = BRCalc->KinematicScaling( kHNLProdMuon3Numu );
    KScale[1] = BRCalc->KinematicScaling( kHNLProdMuon3Nue ); // same, convenience for later
    KScale[2] = BRCalc->KinematicScaling( kHNLProdMuon3Nutau ); // same, convenience for later
    mixScale[0] = 1.0 * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = 1.0 * Ue42 * KScale[1]; totalMix += mixScale[1];
    mixScale[2] = 1.0 * Ut42 * KScale[2]; totalMix += mixScale[2];

    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdMuon3Numu,  mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdMuon3Nue,   mixScale[1] / totalMix } ) );
    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdMuon3Nutau, mixScale[2] / totalMix } ) );

    // it can happen that HNL is not coupled to the only kinematically available channel.
    // Return bogus map if that happens
    if( totalMix <= 0.0 ){
      dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdNull, -999.9 } ) );
      dynamicScores_muon = dynScores;
      return dynScores;
    }

    dynamicScores_muon = dynScores;
    break;
  case genie::kPdgKP:
    KScale[0] = BRCalc->KinematicScaling( kHNLProdKaon2Muon );
    KScale[1] = BRCalc->KinematicScaling( kHNLProdKaon2Electron );
    KScale[2] = BRCalc->KinematicScaling( kHNLProdKaon3Muon );
    KScale[3] = BRCalc->KinematicScaling( kHNLProdKaon3Electron );
    mixScale[0] = BR_K2mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_K2e  * Ue42 * KScale[1]; totalMix += mixScale[1];
    mixScale[2] = BR_K3mu * Um42 * KScale[2]; totalMix += mixScale[2];
    mixScale[3] = BR_K3e  * Ue42 * KScale[3]; totalMix += mixScale[3];

    if( totalMix <= 0.0 ){
      dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdNull, -999.9 } ) );
      dynamicScores_pion = dynScores;
      return dynScores;
    }

    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdKaon2Muon,     mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdKaon2Electron, mixScale[1] / totalMix } ) );
    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdKaon3Muon,     mixScale[2] / totalMix } ) );
    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdKaon3Electron, mixScale[3] / totalMix } ) );

    dynamicScores_kaon = dynScores;
    break;
  case genie::kPdgPiP:

    KScale[0] = BRCalc->KinematicScaling( kHNLProdPion2Muon );
    KScale[1] = BRCalc->KinematicScaling( kHNLProdPion2Electron );
    mixScale[0] = BR_pi2mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_pi2e  * Ue42 * KScale[1]; totalMix += mixScale[1];

    if( totalMix <= 0.0 ){
      dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdNull, -999.9 } ) );
      dynamicScores_pion = dynScores;
      return dynScores;
    }

    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdPion2Muon,     mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdPion2Electron, mixScale[1] / totalMix } ) );

    dynamicScores_pion = dynScores;
    break;
  case genie::kPdgK0L:

    KScale[0] = BRCalc->KinematicScaling( kHNLProdNeuk3Muon );
    KScale[1] = BRCalc->KinematicScaling( kHNLProdNeuk3Electron );
    mixScale[0] = BR_K03mu * Um42 * KScale[0]; totalMix += mixScale[0];
    mixScale[1] = BR_K03e  * Ue42 * KScale[1]; totalMix += mixScale[1];
    
    if( totalMix <= 0.0 ){
      dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdNull, -999.9 } ) );
      dynamicScores_neuk = dynScores;
      return dynScores;
    }

    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdNeuk3Muon,     mixScale[0] / totalMix } ) );
    dynScores.insert( std::pair< HNLProd_t, double >( { kHNLProdNeuk3Electron, mixScale[1] / totalMix } ) );

    dynamicScores_neuk = dynScores;
    break;
  default:
    LOG( "HNL", pERROR )
      << "Unknown parent particle. Cannot make scales, exiting."; exit(1);
  }

  LOG( "HNL", pDEBUG )
    << "Score map now has " << dynScores.size() << " elements. Returning.";
  return dynScores;

}
//----------------------------------------------------------------------------
TLorentzVector FluxCreator::HNLEnergy( HNLProd_t hnldm, TLorentzVector p4par ) const
{
  // first boost to parent rest frame
  TLorentzVector p4par_rest = p4par;
  TVector3 boost_beta = p4par.BoostVector();
  p4par_rest.Boost( -boost_beta );

  LOG( "HNL", pDEBUG )
    << "Attempting to decay rest-system p4 = " << utils::print::P4AsString(&p4par_rest)
    << " as " << utils::hnl::ProdAsString( hnldm );
  
  // get PDGCodeList and truncate 1st member
  PDGCodeList fullList  = utils::hnl::ProductionProductList( hnldm );
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
    LOG("HNL", pERROR)
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
  double wmax = -1;
  for(int idec=0; idec<200; idec++) {
     double w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  wmax *= 2;

  LOG("HNL", pNOTICE)
     << "Max phase space gen. weight @ current HNL system: " << wmax;

  // Generate an unweighted decay
  RandomGen * rnd = RandomGen::Instance();
  
  bool accept_decay=false;
  unsigned int itry=0;
  while(!accept_decay)
    {
      itry++;
      
      if(itry > controls::kMaxUnweightDecayIterations) {
	// report, clean-up and return
	LOG("HNL", pWARN)
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
        LOG("HNL", pWARN)
	  << "Decay weight = " << w << " > max decay weight = " << wmax;
      }
      double gw = wmax * rnd->RndHadro().Rndm();
      accept_decay = (gw<=w);
      
      LOG("HNL", pINFO)
        << "Decay weight = " << w << " / R = " << gw
        << " - accepted: " << accept_decay;
      
    } //!accept_decay
  
  // Grab 0th entry energy and return that
  int idp = 0; TLorentzVector p4HNL, p4HNL_rest;
  // search for the charged lepton in this stack, get the 4-vector in parent rest frame
  TLorentzVector p4Lep, p4Lep_HNLRest;

  p4HNL.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );
  p4HNL_rest.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );
  p4Lep.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );
  p4Lep_HNLRest.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );

  for(std::vector<int>::const_iterator pdg_iter = decayList.begin(); pdg_iter != decayList.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);

     if( std::abs( pdgc ) == kPdgHNL ) p4HNL = *p4fin;
     if( std::abs( pdgc ) == kPdgElectron || 
	 std::abs( pdgc ) == kPdgMuon ||
	 std::abs( pdgc ) == kPdgTau ){
       p4Lep = *p4fin;
       fLepPdg = pdgc;
     }
     idp++;
  }
  
  if( doPol ){
    // boost this to HNL rest frame.
    TVector3 boostHNL = p4HNL.BoostVector();
    p4Lep_HNLRest = p4Lep; fLPE = p4Lep_HNLRest.E(); // still in parent rest frame here
    p4Lep_HNLRest.Boost( boostHNL );
    // save this
    fLPx = ( fixPol ) ? fFixedPolarisation.at(0) : p4Lep.Px() / p4Lep.P(); // note that this is for a true random decay.
    fLPy = ( fixPol ) ? fFixedPolarisation.at(1) : p4Lep.Py() / p4Lep.P(); // We still need to take the geometrical
    fLPz = ( fixPol ) ? fFixedPolarisation.at(2) : p4Lep.Pz() / p4Lep.P(); // constraint into account.
  } else {
    fLPx = 0.0;
    fLPy = 0.0;
    fLPz = 0.0;
  }
  
  delete [] mass;
  return p4HNL; // rest frame momentum!
}
//----------------------------------------------------------------------------
TVector3 FluxCreator::PointToRandomPointInBBox( ) const
{
  RandomGen * rnd = RandomGen::Instance();
  /*
  double ox = fCx + fDetOffset.at(0), oy = fCy + fDetOffset.at(1), oz = fCz + fDetOffset.at(2); // NEAR, m
  */
  double ox = fCx, oy = fCy, oz = fCz; // NEAR, m
  
  double rx = (rnd->RndGen()).Uniform( -fLx/2.0, fLx/2.0 ), 
    ry = (rnd->RndGen()).Uniform( -fLy/2.0, fLy/2.0 ),
    rz = (rnd->RndGen()).Uniform( -fLz/2.0, fLz/2.0 ); // USER, m

  double ux = (rx + fDetOffset.at(0)) * units::m / units::cm;
  double uy = (ry + fDetOffset.at(1)) * units::m / units::cm;
  double uz = (rz + fDetOffset.at(2)) * units::m / units::cm;
  TVector3 checkPoint( ux, uy, uz ); // USER, cm
  
  TVector3 originPoint( -ox, -oy, -oz );
  TVector3 dumori( 0.0, 0.0, 0.0 );
  if( !fDoingOldFluxCalc ){
    // user-coordinates of this point. [cm]
    //double ux = rx - ox, uy = ry - oy, uz = rz - oz;

    LOG( "HNL", pDEBUG )
      << "\nChecking point " << utils::print::Vec3AsString(&checkPoint) << " [m, user]";

    // check if the point is inside the geometry, otherwise do it again
    std::string pathString = this->CheckGeomPoint( ux, uy, uz ); int iNode = 1; // 1 past beginning
    int iBad = 0;
    while( pathString.find( "/", iNode ) == string::npos && iBad < 10 ){
      rx = (rnd->RndGen()).Uniform( -fLx/2.0, fLx/2.0 ); ux = (rx + fDetOffset.at(0)) * units::m / units::cm;
      ry = (rnd->RndGen()).Uniform( -fLy/2.0, fLy/2.0 ); uy = (ry + fDetOffset.at(1)) * units::m / units::cm;
      rz = (rnd->RndGen()).Uniform( -fLz/2.0, fLz/2.0 ); uz = (rz + fDetOffset.at(2)) * units::m / units::cm;
      checkPoint.SetXYZ( ux, uy, uz );
      pathString = this->CheckGeomPoint( ux, uy, uz ); iNode = 1;
      iBad++;
    }
    assert( pathString.find( "/", iNode ) != string::npos );
  }

  // turn u back into [m] from [cm]
  ux *= units::cm / units::m; uy *= units::cm / units::m; uz *= units::cm / units::m;
  // return the absolute point in space [NEAR, m] that we're pointing to!
  checkPoint.SetXYZ( ux, uy, uz );
  checkPoint = this->ApplyUserRotation( checkPoint, dumori, fDetRotation, true ); // det --> tgt-hall
  checkPoint.SetXYZ( checkPoint.X() + ox, checkPoint.Y() + oy, checkPoint.Z() + oz );

  TVector3 vec( ux, uy, uz ); // USER m
  LOG( "HNL", pDEBUG )
    << "\nPointing to this point in BBox (USER coords): " << utils::print::Vec3AsString( &vec ) << "[m]"
    << "\nIn NEAR coords this is " << utils::print::Vec3AsString( &checkPoint ) << "[m]";

  // update bookkeeping
  fTargetPoint = checkPoint;
  
  return checkPoint;
}
//----------------------------------------------------------------------------
double FluxCreator::GetAngDeviation( TLorentzVector p4par, TVector3 detO, bool seekingMax ) const
{
  TVector3 ppar = p4par.Vect(); assert( ppar.Mag() > 0.0 );
  TVector3 pparUnit = ppar.Unit();
  // let face be planar and perpendicular to vector Q
  // assuming Q = ( 0, 0, 1 ) == face perpendicular to z
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

  // see assumption about Q
  // to fix probably with a rotation of fLx, fLy by Euler angles onto Q-plane?
  // line: r(t) - r(incp) = t * IPdev
  // assume square face
  double ttx = ( IPdev.X() != 0.0 ) ? fLx / std::abs( IPdev.X() ) : 99999.9;
  double tty = ( IPdev.Y() != 0.0 ) ? fLy / std::abs( IPdev.Y() ) : 99999.9;
  double tt = std::max( ttx, tty ); // this defines how much the least sweep goes
  TVector3 atilde( tt * IPdev.X(), tt * IPdev.Y(), tt * IPdev.Z() );
  
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
    // actually, if IPdev endpoint lies inside detector have to go other way
    double dh = fLT + dist, dl = fLT - dist;
    // get those vectors and do inner product magic
    TVector3 ph( x_incp + dh * rVec.X(), y_incp + dh * rVec.Y(), z_incp + dh * rVec.Z() );
    TVector3 pl( x_incp - dl * rVec.X(), y_incp - dl * rVec.Y(), z_incp - dl * rVec.Z() );
    double thh = TMath::ACos( ( ph.X()*pparUnit.X() + ph.Y()*pparUnit.Y() + ph.Z()*pparUnit.Z() ) / ph.Mag() );
    double thl = TMath::ACos( ( pl.X()*pparUnit.X() + pl.Y()*pparUnit.Y() + pl.Z()*pparUnit.Z() ) / pl.Mag() );

    return ( seekingMax ) ? thh * 180.0 / constants::kPi : thl * 180.0 / constants::kPi;
  }

  LOG( "HNL", pERROR )
    << "Could not calculate the angle range for detector at ( " << detO.X() << ", " 
    << detO.Y() << ", " << detO.Z() << " ) [m] from HNL production point with parent momentum = ( "
    << ppar.X() << ", " << ppar.Y() << ", " << ppar.Z() << " ) [GeV]. Returning zero.";
  return 0.0;
}
//----------------------------------------------------------------------------
void FluxCreator::GetAngDeviation( TLorentzVector p4par, TVector3 detO, double &zm, double &zp ) const
{
  // implementation of GetAngDeviation that uses ROOT geometry. More robust than analytical geom
  // (fewer assumptions about detector position)

  TVector3 ppar = p4par.Vect(); assert( ppar.Mag() > 0.0 );
  TVector3 pparUnit = ppar.Unit();

  const double sx1 = TMath::Sin(fBx1), cx1 = TMath::Cos(fBx1);
  const double sz = TMath::Sin(fBz), cz = TMath::Cos(fBz);
  const double sx2 = TMath::Sin(fBx2), cx2 = TMath::Cos(fBx2);

  const double xun[3] = { cz, -cx1*sz, sx1*sz };
  const double yun[3] = { sz*cx2, cx1*cz*cx2 - sx1*sx2, -sx1*cz*cx2 - cx1*sx2 };
  const double zun[3] = { sz*sx2, cx1*cz*sx2 + sx1*cx2, -sx1*cz*sx2 + cx1*cx2 };

  LOG("HNL", pDEBUG)
    << "\nxun = ( " << xun[0] << ", " << xun[1] << ", " << xun[2] << " )"
    << "\nyun = ( " << yun[0] << ", " << yun[1] << ", " << yun[2] << " )"
    << "\nzun = ( " << zun[0] << ", " << zun[1] << ", " << zun[2] << " )";

  /*
  TVector3 detO_cm( (detO.X() + fDetOffset.at(0)) * units::m / units::cm, 
		    (detO.Y() + fDetOffset.at(1)) * units::m / units::cm,
		    (detO.Z() + fDetOffset.at(2)) * units::m / units::cm );
  */
  TVector3 detO_cm( (detO.X()) * units::m / units::cm, 
		    (detO.Y()) * units::m / units::cm,
		    (detO.Z()) * units::m / units::cm );
  double inProd = zun[0] * detO_cm.X() + zun[1] * detO_cm.Y() + zun[2] * detO_cm.Z(); // cm
  assert( pparUnit.X() * zun[0] + pparUnit.Y() * zun[1] + pparUnit.Z() * zun[2] != 0.0 );
  inProd /= ( pparUnit.X() * zun[0] + pparUnit.Y() * zun[1] + pparUnit.Z() * zun[2] );

  // using vector formulation, find point of closest approach between parent momentum from
  // parent decay point, and detector centre.

  TVector3 dumori(0.0, 0.0, 0.0); // tgt-hall frame origin is 0
  TVector3 detori( (fDetOffset.at(0)) * units::m / units::cm,
		   (fDetOffset.at(1)) * units::m / units::cm,
		   (fDetOffset.at(2)) * units::m / units::cm ); // for rotations of the detector

  // Do all this in NEAR coords and m to avoid ambiguity.

  TVector3 fCvec_near( fCx, fCy, fCz );  // NEAR m
  TVector3 fDvec( fDx, fDy, fDz ); // in BEAM coords, m
  TVector3 fDvec_beam = this->ApplyUserRotation( fDvec, true ); // in NEAR coords, m
  TVector3 detO_near( fCvec_near.X() - fDvec_beam.X(),
		      fCvec_near.Y() - fDvec_beam.Y(),
		      fCvec_near.Z() - fDvec_beam.Z() );
  TVector3 detO_user = this->ApplyUserRotation( detO_near, dumori, fDetRotation, false ); // tgt-hall --> det

  const double aConst[3] = { fDvec_beam.X(), fDvec_beam.Y(), fDvec_beam.Z() };
  const double aConstUser[3] = { -detO_user.X(), -detO_user.Y(), -detO_user.Z() };
  /*
  const double dConst[3] = { fCx + fDetOffset.at(0), fCy + fDetOffset.at(1), fCz + fDetOffset.at(2) };
  */
  const double dConst[3] = { fCx, fCy, fCz };
  const double nConst[3] = { pparUnit.X(), pparUnit.Y(), pparUnit.Z() };

  // formula for POCA is \vec{a} + \vec{n} * ( (\vec{d} - \vec{a}) \cdot \vec{n} )
  const double nMult = nConst[0] * (dConst[0] - aConst[0]) +
    nConst[1] * (dConst[1] - aConst[1]) +
    nConst[2] * (dConst[2] - aConst[2]);

  // don't use the actual POCA, force calculation from that point V0 such that z(V0) = z(C)
  // and V0 lies on line joining decay point and POCA
  
  const double POCA_m[3] = { aConst[0] + nMult * nConst[0],
			     aConst[1] + nMult * nConst[1],
			     aConst[2] + nMult * nConst[2] }; // NEAR
  /*
  const double zConstMult = ((fCz + fDetOffset.at(2)) - aConst[2]) / (POCA_m[2] - aConst[2]);
  */
  const double zConstMult = ((fCz) - aConst[2]) / (POCA_m[2] - aConst[2]);
  const double startPoint_m[3] = { aConst[0] + nMult * zConstMult * nConst[0],
				   aConst[1] + nMult * zConstMult * nConst[1],
				   aConst[2] + nMult * zConstMult * nConst[2] }; // NEAR

  LOG( "HNL", pDEBUG )
    << "\ndetO_cm = " << utils::print::Vec3AsString( &detO_cm )
    << "\npparUnit = " << utils::print::Vec3AsString( &pparUnit );

  const double startPoint[3] = { startPoint_m[0] * units::m / units::cm,
				 startPoint_m[1] * units::m / units::cm,
				 startPoint_m[2] * units::m / units::cm }; // NEAR, cm

  /*
  const double sweepVect[3] = { (fCx + fDetOffset.at(0)) * units::m / units::cm - startPoint[0],
				(fCy + fDetOffset.at(1)) * units::m / units::cm - startPoint[1],
				(fCz + fDetOffset.at(2)) * units::m / units::cm - startPoint[2] }; // NEAR, cm
  */
  const double sweepVect[3] = { (fCx) * units::m / units::cm - startPoint[0],
				(fCy) * units::m / units::cm - startPoint[1],
				(fCz) * units::m / units::cm - startPoint[2] }; // NEAR, cm
  const double swvMag = std::sqrt( sweepVect[0]*sweepVect[0] + sweepVect[1]*sweepVect[1] + sweepVect[2]*sweepVect[2] ); assert( swvMag > 0.0 );

  // Note the geometry manager works in the *detector frame*. Transform to that.
  /*
  TVector3 detStartPoint( startPoint[0] - (fCx + fDetOffset.at(0)) * units::m / units::cm,
			  startPoint[1] - (fCy + fDetOffset.at(1)) * units::m / units::cm,
			  startPoint[2] - (fCz + fDetOffset.at(2)) * units::m / units::cm );
  */
  TVector3 detStartPoint( startPoint[0] - (fCx) * units::m / units::cm,
			  startPoint[1] - (fCy) * units::m / units::cm,
			  startPoint[2] - (fCz) * units::m / units::cm );
  TVector3 detSweepVect( sweepVect[0], sweepVect[1], sweepVect[2] );

  TVector3 detPpar = ppar;

  LOG( "HNL", pDEBUG )
    << "\nStartPoint = " << utils::print::Vec3AsString( &detStartPoint )
    << "\nSweepVect  = " << utils::print::Vec3AsString( &detSweepVect )
    << "\nparent p3  = " << utils::print::Vec3AsString( &detPpar );
  
  detStartPoint = this->ApplyUserRotation( detStartPoint, detori, fDetRotation, false ); // passive transformation
  detSweepVect = this->ApplyUserRotation( detSweepVect, dumori, fDetRotation, true );
  detPpar = this->ApplyUserRotation( detPpar, dumori, fDetRotation, true );

  // first check that detStartPoint is not already in the detector! If it is, we should flag this now.
  std::string detPathString = this->CheckGeomPoint( detStartPoint.X(), detStartPoint.Y(), detStartPoint.Z() ); int iDNode = 1; // 1 past beginning
  bool startsInsideDet = ( detPathString.find("/", iDNode) != string::npos );

  TLorentzVector detPpar_4v( detPpar.X(), detPpar.Y(), detPpar.Z(), p4par.E() );
  
  // now sweep along sweepVect until we hit either side of the detector. 
  // This will give us two points in space

  double minusPoint[3] = { 0.0, 0.0, 0.0 }; double plusPoint[3] = { 0.0, 0.0, 0.0 }; // USER, cm
  if( startsInsideDet ){
    minusPoint[0] = (gGeoManager->GetCurrentPoint())[0];
    minusPoint[1] = (gGeoManager->GetCurrentPoint())[1];
    minusPoint[2] = (gGeoManager->GetCurrentPoint())[2];
  }

  gGeoManager->SetCurrentPoint( detStartPoint.X(), detStartPoint.Y(), detStartPoint.Z() );
  gGeoManager->SetCurrentDirection( detSweepVect.X() / swvMag, detSweepVect.Y() / swvMag, detSweepVect.Z() / swvMag );
  const double sStepSize = 0.05 * std::min( std::min( fLx, fLy ), fLz );

  // start stepping. Let's do this manually cause FindNextBoundaryAndStep() can be finicky.
  // if start inside detector, then only find exit point, otherwise find entry point.
  if( startsInsideDet ){ // only find exit
    
    // to avoid large time inefficiencies just go to the edge of the box and step back by 10% of the size
    // this is *roughly correct* if not hugely correct

    double currx = (gGeoManager->GetCurrentPoint())[0];
    double curry = (gGeoManager->GetCurrentPoint())[1];
    double currz = (gGeoManager->GetCurrentPoint())[2];
    double currDist = std::sqrt( currx*currx + curry*curry + currz*currz ); // cm

    double curdx = (gGeoManager->GetCurrentDirection())[0];
    double curdy = (gGeoManager->GetCurrentDirection())[1];
    double curdz = (gGeoManager->GetCurrentDirection())[2];

    double stepMod = 0.99;
    double boxSize = std::sqrt( fLxR*fLxR + fLyR*fLyR + fLzR*fLzR )/2.0; // m
    double halfBoxSize = boxSize/2.0;
    double desiredDist = stepMod * halfBoxSize * units::m / units::cm; // cm

    // now we're at distance d on the one side of our detector centre
    // we want to get to distance D on the other side of our detector centre
    // meaning we should step by d + D
    double largeStep = currDist + desiredDist;

    gGeoManager->SetCurrentPoint( currx + largeStep * curdx,
				  curry + largeStep * curdy,
				  currz + largeStep * curdz );

    /* // this is the very correct, very slow way of doing it
    while( detPathString.find("/", iDNode) != string::npos ){
      gGeoManager->SetCurrentPoint( (gGeoManager->GetCurrentPoint())[0] + (gGeoManager->GetCurrentDirection())[0] * sStepSize,
				    (gGeoManager->GetCurrentPoint())[1] + (gGeoManager->GetCurrentDirection())[1] * sStepSize,
				    (gGeoManager->GetCurrentPoint())[2] + (gGeoManager->GetCurrentDirection())[2] * sStepSize );
      detPathString = this->CheckGeomPoint( (gGeoManager->GetCurrentPoint())[0], (gGeoManager->GetCurrentPoint())[1], (gGeoManager->GetCurrentPoint())[2] );
    }
    */

    plusPoint[0] = (gGeoManager->GetCurrentPoint())[0];
    plusPoint[1] = (gGeoManager->GetCurrentPoint())[1];
    plusPoint[2] = (gGeoManager->GetCurrentPoint())[2];
  } else { // find entry and exit
    // first, entry

    // to avoid large time inefficiencies just go to the edge of the box and step back by 10% of the size
    // this is *roughly correct* if not hugely correct

    // find that point that lies on the line along direction from current point that is 90% of box size from box centre
    double currx = (gGeoManager->GetCurrentPoint())[0];
    double curry = (gGeoManager->GetCurrentPoint())[1];
    double currz = (gGeoManager->GetCurrentPoint())[2];
    double currDist = std::sqrt( currx*currx + curry*curry + currz*currz );

    double stepMod = 0.99;
    double boxSize = std::sqrt( fLxR*fLxR + fLyR*fLyR + fLzR*fLzR )/2.0; // m
    double halfBoxSize = boxSize / 2.0;
    double desiredDist = stepMod * halfBoxSize * units::m / units::cm;

    // we're at some distance D and want to get to distance d, which means we step forward by
    // R = D-d = D(1-d/D)
    double largeStep = currDist - desiredDist;

    double curdx = (gGeoManager->GetCurrentDirection())[0];
    double curdy = (gGeoManager->GetCurrentDirection())[1];
    double curdz = (gGeoManager->GetCurrentDirection())[2];

    gGeoManager->SetCurrentPoint( currx + largeStep * curdx,
				  curry + largeStep * curdy,
				  currz + largeStep * curdz );

    /* // slow but correct
    while( detPathString.find("/", iDNode) == string::npos ){
      gGeoManager->SetCurrentPoint( (gGeoManager->GetCurrentPoint())[0] + (gGeoManager->GetCurrentDirection())[0] * sStepSize,
				    (gGeoManager->GetCurrentPoint())[1] + (gGeoManager->GetCurrentDirection())[1] * sStepSize,
				    (gGeoManager->GetCurrentPoint())[2] + (gGeoManager->GetCurrentDirection())[2] * sStepSize );
      detPathString = this->CheckGeomPoint( (gGeoManager->GetCurrentPoint())[0], (gGeoManager->GetCurrentPoint())[1], (gGeoManager->GetCurrentPoint())[2] );
    }
    */

    minusPoint[0] = (gGeoManager->GetCurrentPoint())[0];
    minusPoint[1] = (gGeoManager->GetCurrentPoint())[1];
    minusPoint[2] = (gGeoManager->GetCurrentPoint())[2];

    // then, exit

    // quick and dirty way: reflect about the origin

    currx = (gGeoManager->GetCurrentPoint())[0];
    curry = (gGeoManager->GetCurrentPoint())[1];
    currz = (gGeoManager->GetCurrentPoint())[2];
    /*
    double cdx = fDetOffset.at(0) * units::m / units::cm - currx; // cm
    double cdy = fDetOffset.at(1) * units::m / units::cm - curry; // cm
    double cdz = fDetOffset.at(2) * units::m / units::cm - currz; // cm

    gGeoManager->SetCurrentPoint( fDetOffset.at(0) * units::m / units::cm + cdx,
				  fDetOffset.at(1) * units::m / units::cm + cdy,
				  fDetOffset.at(2) * units::m / units::cm + cdz );
    */

    gGeoManager->SetCurrentPoint( -currx, -curry, -currz );

    /* // slow but correct
    while( detPathString.find("/", iDNode) != string::npos ){
      gGeoManager->SetCurrentPoint( (gGeoManager->GetCurrentPoint())[0] + (gGeoManager->GetCurrentDirection())[0] * sStepSize,
				    (gGeoManager->GetCurrentPoint())[1] + (gGeoManager->GetCurrentDirection())[1] * sStepSize,
				    (gGeoManager->GetCurrentPoint())[2] + (gGeoManager->GetCurrentDirection())[2] * sStepSize );
      detPathString = this->CheckGeomPoint( (gGeoManager->GetCurrentPoint())[0], (gGeoManager->GetCurrentPoint())[1], (gGeoManager->GetCurrentPoint())[2] );
    }
    */

    plusPoint[0] = (gGeoManager->GetCurrentPoint())[0];
    plusPoint[1] = (gGeoManager->GetCurrentPoint())[1];
    plusPoint[2] = (gGeoManager->GetCurrentPoint())[2];
  }

  /*
  int bdIdx = 0; const int bdIdxMax = 1e+4;
  if(!startsInsideDet){
    while( gGeoManager->FindNextBoundaryAndStep() && bdIdx < bdIdxMax ){
      bdIdx++;
      if( bdIdx % 100  == 0 )
	LOG( "HNL", pDEBUG ) << "bdIdx = " << bdIdx;
      if( std::abs( (gGeoManager->GetCurrentPoint())[0] ) != fLxR/2.0 * units::m / units::cm &&
	  std::abs( (gGeoManager->GetCurrentPoint())[1] ) != fLyR/2.0 * units::m / units::cm &&
	  std::abs( (gGeoManager->GetCurrentPoint())[2] ) != fLzR/2.0 * units::m / units::cm ){
	plusPoint[0] = (gGeoManager->GetCurrentPoint())[0];
	plusPoint[1] = (gGeoManager->GetCurrentPoint())[1];
	plusPoint[2] = (gGeoManager->GetCurrentPoint())[2]; 
      }
    }
  }
  */

  TVector3 originPoint( -(fCx + fDetOffset.at(0)), -(fCy + fDetOffset.at(1)), -(fCz + fDetOffset.at(2)) );

  TVector3 vMUser( minusPoint[0] * units::cm / units::m, 
		   minusPoint[1] * units::cm / units::m,
		   minusPoint[2] * units::cm / units::m ); // USER, m
  TVector3 vPUser( plusPoint[0] * units::cm / units::m,
		   plusPoint[1] * units::cm / units::m,
		   plusPoint[2] * units::cm / units::m ); // USER, m

  // rotate to NEAR frame
  TVector3 vMNear = this->ApplyUserRotation( vMUser, originPoint, fDetRotation, true ); // det --> tgt-hall
  /*
  vMNear.SetXYZ( vMNear.X() + (fCx + fDetOffset.at(0)), 
		 vMNear.Y() + (fCy + fDetOffset.at(1)),
		 vMNear.Z() + (fCz + fDetOffset.at(2)) );
  */
  vMNear.SetXYZ( vMNear.X() + (fCx),
		 vMNear.Y() + (fCy),
		 vMNear.Z() + (fCz) );
  TVector3 vPNear = this->ApplyUserRotation( vPUser, originPoint, fDetRotation, true ); // det --> tgt-hall
  /*
  vPNear.SetXYZ( vPNear.X() + (fCx + fDetOffset.at(0)), 
		 vPNear.Y() + (fCy + fDetOffset.at(1)),
		 vPNear.Z() + (fCz + fDetOffset.at(2)) );
  */
  vPNear.SetXYZ( vPNear.X() + (fCx), 
		 vPNear.Y() + (fCy),
		 vPNear.Z() + (fCz) );

  // with 3 points and 1 vector we calculate the angles.
  // Points: D(decay), E(entry), X(exit) [all in local, cm]. Vector: detPpar [local, GeV/GeV]
  // angles are <DE, detPpar> and <DX, detPpar>

  // now obtain the angles themselves and return in deg.
  /*
  TVector3 decayPoint_user( (fDx - (fCx + fDetOffset.at(0))) * units::m / units::cm, 
			    (fDy - (fCy + fDetOffset.at(1))) * units::m / units::cm, 
			    (fDz - (fCz + fDetOffset.at(2))) * units::m / units::cm ); // USER, cm
  */
  TVector3 decayPoint_user( aConstUser[0] * units::m / units::cm, 
			    aConstUser[1] * units::m / units::cm,
			    aConstUser[2] * units::m / units::cm );
  TVector3 decayPoint_near = this->ApplyUserRotation( decayPoint_user, detori, fDetRotation, false ); // NEAR, cm
  /*
  decayPoint_near.SetXYZ( decayPoint_near.X() + (fCx + fDetOffset.at(0)) * units::m / units::cm,
			  decayPoint_near.Y() + (fCy + fDetOffset.at(1)) * units::m / units::cm,
			  decayPoint_near.Z() + (fCz + fDetOffset.at(2)) * units::m / units::cm );
  */
  decayPoint_near.SetXYZ( decayPoint_near.X() + (fCx) * units::m / units::cm,
			  decayPoint_near.Y() + (fCy) * units::m / units::cm,
			  decayPoint_near.Z() + (fCz) * units::m / units::cm );

  TVector3 minusVec( minusPoint[0] - decayPoint_user.X(), minusPoint[1] - decayPoint_user.Y(), minusPoint[2] - decayPoint_user[2] ); // USER, cm
  TVector3 plusVec( plusPoint[0] - decayPoint_user.X(), plusPoint[1] - decayPoint_user.Y(), plusPoint[2] - decayPoint_user[2] ); // USER, cm
  TVector3 startVec( detStartPoint[0] - decayPoint_user.X(), detStartPoint[1] - decayPoint_user.Y(), detStartPoint[2] - decayPoint_user.Z() ); // USER, cm

  TVector3 minusVec_near = this->ApplyUserRotation( minusVec, detori, fDetRotation, false ); // NEAR, cm
  TVector3 plusVec_near = this->ApplyUserRotation( plusVec, detori, fDetRotation, false ); // NEAR, cm
  TVector3 startVec_near = this->ApplyUserRotation( startVec, detori, fDetRotation, false ); // NEAR, cm

  double minusNum = startVec.X() * minusVec.X() + startVec.Y() * minusVec.Y() + startVec.Z() * minusVec.Z(); // USER AND USER
  double minusDen = startVec.Mag() * minusVec.Mag(); assert( minusDen > 0.0 ); // USER AND USER

  zm = TMath::ACos( minusNum / minusDen ) * TMath::RadToDeg();

  double plusNum = startVec.X() * plusVec.X() + startVec.Y() * plusVec.Y() + startVec.Z() * plusVec.Z(); // USER AND USER
  double plusDen = startVec.Mag() * plusVec.Mag(); assert( plusDen > 0.0 ); // USER AND USER

  zp = TMath::ACos( plusNum / plusDen ) * TMath::RadToDeg();

  if( zm > zp ){
    double tmpzp = zp;
    zp = zm;
    zm = tmpzp;
  }

  /*
  LOG( "HNL", pDEBUG )
    << "\nIn DETECTOR coordinates:"
    << "\nentered at ( " << minusPoint[0] << ", " << minusPoint[1] << ", " << minusPoint[2] << " ) [cm]"
    << "\nexited  at ( " << plusPoint[0] << ", " << plusPoint[1] << ", " << plusPoint[2] << " ) [cm]"
    << "\nstarted at ( " << decVec.X() << ", " << decVec.Y() << ", " << decVec.Z() << " ) [cm]"
    << "\nmomentum   ( " << detPpar.X() << ", " << detPpar.Y() << ", " << detPpar.Z() << " )"
    << "\nmeaning zm = " << zm << ", zp = " << zp << " [deg]";
  */
}
//----------------------------------------------------------------------------
double FluxCreator::CalculateAcceptanceCorrection( TLorentzVector p4par, TLorentzVector p4HNL,
						      double SMECM, double zm, double zp ) const
{
  /*
   * This method calculates HNL acceptance by taking into account the collimation effect
   * HNL are massive so Lorentz boost from parent CM ==> lab is more effective
   * This means that, given a desired range of lab-frame emission angles, there are
   * more rest-frame emission angles that map into this range. 
   * Find the measure of the rest-frame that maps onto the allowed lab-frame angles
   * and return the ratio over the relevant measure for a SM neutrino
   */

  assert( zm >= 0.0 && zp >= zm );
  if( zp == zm ) return 1.0;

  double M = p4HNL.M();
  if( M < 1.0e-3 ) return 1.0;

  TF1 * fHNL = ( TF1* ) gROOT->GetListOfFunctions()->FindObject( "fHNL" );
  if( !fHNL ){ fHNL = new TF1( "fHNL", labangle, 0.0, 180.0, 6 ); }
  fHNL->SetParameter( 0, p4par.E()  );
  fHNL->SetParameter( 1, p4par.Px() );
  fHNL->SetParameter( 2, p4par.Py() );
  fHNL->SetParameter( 3, p4par.Pz() );
  fHNL->SetParameter( 4, p4HNL.P()  );
  fHNL->SetParameter( 5, p4HNL.E()  );

  double ymax = fHNL->GetMaximum(), xmax = fHNL->GetMaximumX();
  double range1 = 0.0;

  if( fHNL->GetMinimum() >= fHNL->GetMaximum() ) return 1.0; // bail on constant function

  if( zm < fHNL->GetMinimum() ){ // really good collimation, ignore checks on zm
    double z0 = fHNL->GetMinimum();
    if( ymax > zp && xmax < 180.0 ){ // >=2 pre-images, add them together
      int nPreim = 0;

      // RETHERE: Make this more sophisticated! Assumes 2 preimages, 1 before and 1 after max
      double xl1 = fHNL->GetX( z0, 0.0, xmax );
      double xh1 = fHNL->GetX( zp, 0.0, xmax );
      double xl2 = fHNL->GetX( z0, xmax, 180.0 );
      double xh2 = fHNL->GetX( zp, xmax, 180.0 );

      range1 += std::abs( xl1 - xh1 ) + std::abs( xh2 - xl2 ); nPreim = 2;
    } else if( ymax > zp && xmax == 180.0 ){ // 1 pre-image, SMv-like case
      double xl = fHNL->GetX( z0 ), xh = fHNL->GetX( zp );
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
      double xl1 = fHNL->GetX( zm, 0.0, xmax );
      double xh1 = fHNL->GetX( zp, 0.0, xmax );
      double xl2 = fHNL->GetX( zm, xmax, 180.0 );
      double xh2 = fHNL->GetX( zp, xmax, 180.0 );

      range1 += std::abs( xl1 - xh1 ) + std::abs( xh2 - xl2 ); nPreim = 2;
    } else if( ymax > zp && xmax == 180.0 ){ // 1 pre-image, SMv-like case
      double xl = fHNL->GetX( zm ), xh = fHNL->GetX( zp );
      range1 = std::abs( xh - xl );

    } else if( zm < ymax && ymax <= zp ){ // 1 pre-image
      double xl = fHNL->GetX( zm, 0., xmax ), xh = fHNL->GetX( zm, xmax, 180.0 );
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

  // SMv deviates more from parent than HNL due to masslessness. This means a larger minimum of labangle
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
    
    TLorentzVector vcx( p4HNL.P(), 0.0, 0.0, p4HNL.E() ), 
      vcy( 0.0, p4HNL.P(), 0.0, p4HNL.E() ), 
      vcz( 0.0, 0.0, p4HNL.P(), p4HNL.E() );
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
double FluxCreator::labangle( double * x, double * par )
{
  double xrad = x[0] * TMath::DegToRad();
  double Ehad = par[0], pxhad = par[1], pyhad = par[2], pzhad = par[3];
  double phnl = par[4], Ehnl = par[5];

  TLorentzVector p4had( pxhad, pyhad, pzhad, Ehad );
  TVector3 boost_vec = p4had.BoostVector(); // beta of parent in lab frame

  // assume phi invariance so create HNL rest-frame momentum along y'z' plane
  TLorentzVector pncm( 0.0, phnl * TMath::Sin( xrad ), phnl * TMath::Cos( xrad ), Ehnl );

  // boost into lab frame
  pncm.Boost( boost_vec );
  
  // return lab frame theta wrt parent momentum in deg
  double num = pxhad * pncm.X() + pyhad * pncm.Y() + pzhad * pncm.Z();
  double den = p4had.P() * pncm.P();
  double theta = TMath::ACos( num / den ) * 180.0 / constants::kPi;
  return theta;
}
//----------------------------------------------------------------------------
void FluxCreator::MakeBBox() const
{
  if( fRadius < 0.0 ) fRadius = 1.0;
  LOG( "HNL", pFATAL )
    << "WARNING: This is a bounding box centred at config-given point and radius " << fRadius << " m";

  fLx = fRadius; fLy = fRadius; fLz = fRadius;
}
//----------------------------------------------------------------------------
TVector3 FluxCreator::ApplyUserRotation( TVector3 vec, bool doBackwards ) const
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
TVector3 FluxCreator::ApplyUserRotation( TVector3 vec, TVector3 oriVec, std::vector<double> rotVec, bool doBackwards ) const
{
  double vx = vec.X(), vy = vec.Y(), vz = vec.Z();
  double ox = oriVec.X(), oy = oriVec.Y(), oz = oriVec.Z();
  
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
  return nvec;
}
//____________________________________________________________________________
double FluxCreator::CalculateDetectorAcceptanceSAA( TVector3 detO ) const
{
  // sang is solid-angle / 4pi
  double rad = std::sqrt( detO.X() * detO.X() + detO.Y() * detO.Y() + detO.Z() * detO.Z() );
  double sang = 1.0 - TMath::Cos( TMath::ATan( kRDET / rad ) ); sang *= 0.5;
  return sang;
}
//----------------------------------------------------------------------------
double FluxCreator::CalculateAreaNormalisation()
{
  // for now this is just a square of length kRDET
  // returns 1 / area
  return 1.0 / ( kRDET * kRDET );
}
//----------------------------------------------------------------------------
void FluxCreator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//----------------------------------------------------------------------------
void FluxCreator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//----------------------------------------------------------------------------
void FluxCreator::LoadConfig(void)
{
  if( fIsConfigLoaded ) return;

  LOG("HNL", pDEBUG)
    << "Loading flux-creation parameters from file...";

  this->GetParam( "HNL-Mass", fMass );
  this->GetParamVect( "HNL-LeptonMixing", fU4l2s );
  this->GetParam( "HNL-Majorana", fIsMajorana );
  
  this->GetParamVect( "Near2User_T", fB2UTranslation );
  this->GetParamVect( "Near2User_R", fDetRotation );
  this->GetParamVect( "Near2Beam_R", fB2URotation );
  this->GetParamVect( "DetCentre_User", fDetOffset );

  this->GetParamVect( "ParentPOTScalings", fScales );
  this->GetParam( "DoOldFluxCalculation", fDoingOldFluxCalc );
  this->GetParam( "RerollPoints", fRerollPoints );
  this->GetParam( "CollectionRadius", fRadius );
  this->GetParam( "InputFluxesInBEAM", fSupplyingBEAM );
  this->GetParam( "IncludePolarisation", doPol );
  this->GetParam( "FixPolarisationDirection", fixPol );
  this->GetParamVect( "HNL-PolDir", fFixedPolarisation );

  this->GetParam( "IsParentOnAxis", isParentOnAxis );

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
  if( utils::hnl::IsProdKinematicallyAllowed( kHNLProdMuon3Nue ) ) 
    POTScaleWeight = fScales[0]; // all POT contribute
  else if( utils::hnl::IsProdKinematicallyAllowed( kHNLProdPion2Muon ) ||
	   utils::hnl::IsProdKinematicallyAllowed( kHNLProdPion2Electron ) ) 
    POTScaleWeight = fScales[1]; // muons do not contribute
  else if( utils::hnl::IsProdKinematicallyAllowed( kHNLProdNeuk3Muon ) ||
	   utils::hnl::IsProdKinematicallyAllowed( kHNLProdNeuk3Electron ) )
    POTScaleWeight = fScales[2]; // only kaons contribute
  else if( utils::hnl::IsProdKinematicallyAllowed( kHNLProdKaon2Muon ) || 
	   utils::hnl::IsProdKinematicallyAllowed( kHNLProdKaon2Electron ) ||
	   utils::hnl::IsProdKinematicallyAllowed( kHNLProdKaon3Muon ) ||
	   utils::hnl::IsProdKinematicallyAllowed( kHNLProdKaon3Electron ) )
    POTScaleWeight = fScales[3]; // only charged kaons contribute

  /*
  LOG( "HNL", pDEBUG )
    << "Read the following parameters :"
    << "\n Mass = " << fMass
    << "\n couplings = " << fU4l2s.at(0) << " : " << fU4l2s.at(1) << " : " << fU4l2s.at(2)
    << "\n translation = " << fB2UTranslation.at(0) << ", " << fB2UTranslation.at(1) << ", " << fB2UTranslation.at(2)
    << "\n rotation = " << fB2URotation.at(0) << ", " << fB2URotation.at(1) << ", " << fB2URotation.at(2)
    << "\n isParentOnAxis = " << isParentOnAxis
    << "\n POTScaleWeight = " << POTScaleWeight;
  */

  fIsConfigLoaded = true;
}
//____________________________________________________________________________
void FluxCreator::SetGeomFile( string geomfile ) const
{
  LOG( "HNL", pDEBUG ) << "Setting geometry file to " << geomfile;
  fGeomFile = geomfile;
}
//____________________________________________________________________________
void FluxCreator::SetFirstFluxEntry( int iFirst ) const
{
  fFirstEntry = iFirst;
}
//____________________________________________________________________________
void FluxCreator::ImportBoundingBox( TGeoBBox * box ) const
{
  LOG( "HNL", pDEBUG ) << "Importing bounding box...";
  fLxR = 2.0 * box->GetDX() * units::cm / units::m;
  fLyR = 2.0 * box->GetDY() * units::cm / units::m;
  fLzR = 2.0 * box->GetDZ() * units::cm / units::m;

  double testRadius = fRadius; // m
  if( !fDoingOldFluxCalc ){
    fLx = std::min( testRadius, fLxR );
    fLy = std::min( testRadius, fLyR );
    fLz = std::min( testRadius, fLzR );
  }
}
//____________________________________________________________________________
std::string FluxCreator::CheckGeomPoint( Double_t x, Double_t y, Double_t z ) const
{
  Double_t point[3];
  Double_t local[3];
  point[0] = x;
  point[1] = y;
  point[2] = z;
  TGeoVolume *vol = gGeoManager->GetTopVolume();
  TGeoNode *node = gGeoManager->FindNode(point[0], point[1], point[2]);
  gGeoManager->MasterToLocal(point, local);
  return gGeoManager->GetPath();
}
