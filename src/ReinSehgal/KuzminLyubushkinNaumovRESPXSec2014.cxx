//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Oct 05, 2009 - CA
   Modified code to handle charged lepton scattering too.
   Also, the helicity amplitude code now returns a `const RSHelicityAmpl &'.
 @ July 23, 2010 - CA
   BaryonResParams, and BreitWignerI, BaryonResDataSetI implementations are
   now redundant. Get resonance parameters from BaryonResUtils and use the
   Breit-Weigner functions from utils::bwfunc.
 @ December 15, 2014 - JN, SD
   Add new version due to Jarek Nowak.  Based on parameters set in
   UserPhysicsOptions.xml.  Default is Berger-Sehgal (Phys. Rev. D76, 
   113004 (2007)) with new GA and new GV.  Additional model due to Kuzmin, 
   Lyubushkin, and Naumov (Mod. Phys. Lett. A19 (2004) 2815 and Phys. Part. 
   Nucl. 35 (2004) S133) also available.  Each of these models
   includes effect of lepton mass.  BS adds a new pole diagram.  
   New GA and GV form factors are based on studies with MiniBooNE data.
 @ Dec 22, 2014 - GP
   Incorporating changes from J. Nowak into a new class (was 
   ReinSehgalRESPXSec, now BergerSehgalRESPXSec2014).
*/
//____________________________________________________________________________

#include "ReinSehgal/KuzminLyubushkinNaumovRESPXSec2014.h"

using namespace genie;

//____________________________________________________________________________
KuzminLyubushkinNaumovRESPXSec2014::KuzminLyubushkinNaumovRESPXSec2014() :
BSKLNBaseRESPXSec2014("genie::KuzminLyubushkinNaumovRESPXSec2014")
{
  this->fKLN = true;
  this->fBRS = false;
}
//____________________________________________________________________________
KuzminLyubushkinNaumovRESPXSec2014::KuzminLyubushkinNaumovRESPXSec2014(string config) :
BSKLNBaseRESPXSec2014("genie::KuzminLyubushkinNaumovRESPXSec2014", config)
{
  this->fKLN = true;
  this->fBRS = false;
}
//____________________________________________________________________________
KuzminLyubushkinNaumovRESPXSec2014::~KuzminLyubushkinNaumovRESPXSec2014()
{

}
//____________________________________________________________________________
