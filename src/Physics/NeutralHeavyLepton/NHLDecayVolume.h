//____________________________________________________________________________
/*!

\namespace genie::NHL::NHLDecayVolume

\brief     Utilities for extracting the position of the NHL decay vertex
           *** given production vertex and momentum ***

\author    John Plows <komninos-john.plows \at physics.ox.ac.uk>
          
\created   March 31st, 2022

\cpright   Copyright (c) 2003-2022, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _NHL_DECAY_VOLUME_H_
#define _NHL_DECAY_VOLUME_H_

#include <cmath>
#include <cassert>

#include <TVector3.h>
//#ifdef __GENIE_GEOM_DRIVERS_ENABLED__ // why do we crash with this guard on?
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
//#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__

//#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/UnitUtils.h"

namespace genie {
namespace NHL {

  namespace NHLDecayVolume {

    // --------------------------------------------------
    // Utilities
    // --------------------------------------------------

    // enforce chosen units
    void EnforceUnits( std::string length_units, std::string angle_units, std::string time_units );

    // calculate travel length in detector
    double CalcTravelLength( double betaMag, double CoMLifetime, double maxLength );

    // assign decay point given length
    TVector3 GetDecayPoint( double travelLength, TVector3 & entryPoint, TVector3 & momentum );

    // get max length in detector
    double GetMaxLength( TVector3 & entryPoint, TVector3 & exitPoint );

    // --------------------------------------------------
    // Simple Decay Volume - fallback if no geometry
    // --------------------------------------------------

    // It's a 1x1x1 m3 box around (0,0,0).
    void MakeSDV();

    // in case of SDV, calculate entry and exit points (if any) of some trajectory
    // if no entry/exit points return false
    bool SDVEntryAndExitPoints( TVector3 & startPoint, TVector3 momentum,
				TVector3 & entryPoint,
				TVector3 & exitPoint );

    // --------------------------------------------------
    // ROOT geometry stuff
    // --------------------------------------------------

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
    // use bounding box origin & sides
    void ImportBoundingBox( TGeoBBox * box );
    
    // get entry & exit points directly from volume
    bool VolumeEntryAndExitPoints( TVector3 & startPoint, TVector3 & momentum,
      TVector3 & entryPoint, TVector3 & exitPoint,
      TGeoManager * gm, TGeoVolume * vol );
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__

    // --------------------------------------------------

    extern double lunits; extern std::string lunitString;
    extern double aunits;
    extern double tunits; extern std::string tunitString;

    extern double fSx, fSy, fSz; //start point
    extern double fPx, fPy, fPz; //momentum
    extern double fEx, fEy, fEz; //entry point
    extern double fXx, fXy, fXz; //exit  point

    extern double fSxROOT, fSyROOT, fSzROOT; // start point in cm
    extern double fExROOT, fEyROOT, fEzROOT; // entry point in cm
    extern double fXxROOT, fXyROOT, fXzROOT; // exit  point in cm

    extern double fDx, fDy, fDz; //decay point
    extern double fOx, fOy, fOz; //origin
    extern double fLx, fLy, fLz; //dimensions

    extern double fDxROOT, fDyROOT, fDzROOT; // decay point in cm
    extern double fOxROOT, fOyROOT, fOzROOT; // origin in cm
    extern double fLxROOT, fLyROOT, fLzROOT; // dimensions in cm

    extern double fAx, fAy, fAz; //
    extern double fAlpha;
    extern double ft;
    
    extern double kNewSpeedOfLight;

  } // namespace NHLDecayVolume

} // namespace NHL
} // namespace genie

#endif // #ifndef _NHL_DECAY_VOLUME_H_
