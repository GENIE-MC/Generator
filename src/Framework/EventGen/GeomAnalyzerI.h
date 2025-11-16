//____________________________________________________________________________
/*!

\class    genie::GeomAnalyzerI

\brief    Defines the GENIE Geometry Analyzer Interface

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  July 13, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _GEOMETRY_ANALYZER_I_H_
#define _GEOMETRY_ANALYZER_I_H_

class TLorentzVector;
class TVector3;

namespace genie {

class PDGCodeList;
class PathLengthList;

class GeomAnalyzerI {

public :

  virtual ~GeomAnalyzerI();

  // define the GeomAnalyzerI interface

  virtual const PDGCodeList &
            ListOfTargetNuclei (void) = 0;

  virtual const PathLengthList &
            ComputeMaxPathLengths (void) = 0;
  virtual const PathLengthList &
            ComputePathLengths (
              const TLorentzVector & x, const TLorentzVector & p) = 0;
  virtual const TVector3 &
            GenerateVertex (
              const TLorentzVector & x, const TLorentzVector & p, int tgtpdg) = 0;

protected:

  GeomAnalyzerI();
};

}      // genie namespace

#endif // _GEOMETRY_ANALYZER_I_H_
