/*!
 * A simple ROOT macro to make a ROOT geometry cylinder.
 * Based on $ROOTSYS/tutorials/geom/rootgeom.C
 *
 * \author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
 *          University of Oxford
 *
 * \cpright Copyright (c) 2003-2025, The GENIE Collaboration
 *          For the full text of the license visit http://copyright.genie-mc.org
 */

void makeCylinder()
{
    //--- define your favourite length units. Default is m
    const std::string lunits = "m";
    double uMult = 1.0;
    
    if( strcmp( lunits.c_str(), "m" ) == 0 ) {
	uMult = 100.0; // m to cm
    } else if( strcmp( lunits.c_str(), "mm" ) == 0 ) {
	uMult = 0.1; // mm to cm
    } else if( strcmp( lunits.c_str(), "cm" ) == 0 ) {
	uMult = 1.0; // cm to cm
    } else {
	std::cerr << "Unknown length units " << lunits.c_str() << ", please add a switch with the proper length conversion. Exiting now." << std::endl;
    }

    TGeoManager * geom = new TGeoManager( "cyl1", "A simple cylindrical detector" );

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);

    //--- make the top container volume
    const double cylRad = 2.6; // m
    const double cylHeight = 2.5; // m
    const double bigBoxSide = 2.0 * std::max( TMath::Sqrt(2.0) * cylRad, cylHeight );
    const double worldLen = 1.01 * bigBoxSide; // m

    TGeoVolume * topvol = geom->MakeBox( "TOP", Vacuum,
					 uMult * worldLen, uMult * worldLen, uMult * worldLen );
    geom->SetTopVolume( topvol );

    //--- do you want to rotate the detector?
    //--- default: no rotation - unit vectors are { (1,0,0), (0,1,0), (0,0,1) }
    //--- Specify 3 extrinsic Euler angles (x-z-x) to rotate these by. Origin is at box centre
    const double rx1DEG = 0.0, rzDEG = 90.0, rx2DEG = 90.0;
    const double rx1 = rx1DEG * TMath::DegToRad(),
	rz = rzDEG * TMath::DegToRad(),
	rx2 = rx2DEG * TMath::DegToRad();

    const double sx1 = TMath::Sin(rx1), cx1 = TMath::Cos(rx1);
    const double sz = TMath::Sin(rz), cz = TMath::Cos(rz);
    const double sx2 = TMath::Sin(rx2), cx2 = TMath::Cos(rx2);

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

    std::cout << xun[0] << " " << xun[1] << " " << xun[2] << std::endl;
    std::cout << yun[0] << " " << yun[1] << " " << yun[2] << std::endl;
    std::cout << zun[0] << " " << zun[1] << " " << zun[2] << std::endl;

    //--- make the actual cylinder and rotate it as desired
    //--- origin is at 1/2 cylinder height, endcap centre
    TGeoVolume * cyl = geom->MakeTube( "CYL", Al,
				       0.0, uMult * cylRad, uMult * cylHeight );
    cyl->SetLineColor(kGreen+2);
    TGeoTranslation * tr0 = new TGeoTranslation( 0.0, 0.0, 0.0 );
    TGeoRotation * rot0 = new TGeoRotation( "rot0", thx, phx, thy, phy, thz, phz );

    //--- add directly to top volume
    topvol->AddNode( cyl, 1, rot0 );

    //--- export this to a file
    geom->Export("./cylinder.root");

    //--- close the geometry
    geom->CloseGeometry();
    
}
