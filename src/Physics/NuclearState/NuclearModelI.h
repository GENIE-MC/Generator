//____________________________________________________________________________
/*!

\class    genie::NuclearModelI

\brief    Pure abstract base class.
          Defines the NuclearModelI interface to be implemented by any physics
          model describing the distribution of nucleons within a nuclei

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  October 09, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          

 Important revisions after version 2.0.0 :
 @ Mar 18, 2016 - JJ (SD)
   Added option for GenerateNucleon() to be called with a target and a radius
   as the arguments. Currently used by LocalFGM. Calls
   GenerateNucleon() with the radius set to 0 for all other NuclearModelI
   implementations.

 @ Jul 2020 - Marco Roda
   Added fooks for FermiMomentum and LocalFermiMomentum

*/
//____________________________________________________________________________

#ifndef _NUCLEAR_MODEL_I_H_
#define _NUCLEAR_MODEL_I_H_

#include <string>

#include <TVector3.h>

#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Target.h"

namespace genie {

class NuclearModelI : public Algorithm {

public:
  virtual ~NuclearModelI() {};

  virtual bool           GenerateNucleon (const Target &) const = 0;
  virtual bool           GenerateNucleon (const Target & tgt,
                                          double hitNucleonRadius) const;

  virtual double         Prob            (double p, double w, const Target &) const = 0;
  virtual double         Prob            (double p, double w, const Target & tgt,
                                          double hitNucleonRadius) const;

  virtual NuclearModel_t ModelType       (const Target &) const = 0;

  virtual double         FermiMomentum( const Target &, int nucleon_pdg ) const ;
  virtual double         LocalFermiMomentum( const Target &, int nucleon_pdg, double radius ) const ; 
 

  inline double          RemovalEnergy   (void)           const
  {
    return fCurrRemovalEnergy;
  }

  inline double         Momentum        (void)           const
  {
    return fCurrMomentum.Mag();
  }

  inline const TVector3& Momentum3      (void)           const
  {
    return fCurrMomentum;
  }

  inline FermiMoverInteractionType_t GetFermiMoverInteractionType(void) const
  {
    return fFermiMoverInteractionType;
  }

  // These setters have to be const. I hate it. We should really update this class interface
  inline void SetMomentum3(const TVector3 & mom) const
  {
    fCurrMomentum = mom;
  };
  inline void SetRemovalEnergy(double E) const
  {
    fCurrRemovalEnergy = E;
  }

protected:
  NuclearModelI()
    : Algorithm()
    , fCurrRemovalEnergy(0)
    , fCurrMomentum(0,0,0)
    , fFermiMoverInteractionType(kFermiMoveDefault)
    , fKFTable(nullptr)
    , fKFTableName("Unspecified")
    {};
  NuclearModelI(std::string name)
    : Algorithm(name)
    , fCurrRemovalEnergy(0)
    , fCurrMomentum(0,0,0)
    , fFermiMoverInteractionType(kFermiMoveDefault)
    , fKFTable(nullptr)
    , fKFTableName("Unspecified")
    {};
  NuclearModelI(std::string name, std::string config)
    : Algorithm(name, config)
    , fCurrRemovalEnergy(0)
    , fCurrMomentum(0,0,0)
    , fFermiMoverInteractionType(kFermiMoveDefault)
    , fKFTable(nullptr)
    , fKFTableName("Unspecified")
    {};

  virtual void LoadConfig() ;

  const string & FermiMomentumTableName() const { return fKFTableName; }
  const genie::FermiMomentumTable & FermiMomentumTable() const { return *fKFTable ; }

  mutable double   fCurrRemovalEnergy;
  mutable TVector3 fCurrMomentum;
  mutable FermiMoverInteractionType_t fFermiMoverInteractionType;


 private:

  const genie::FermiMomentumTable * fKFTable; 
  string fKFTableName;

};

}         // genie namespace
#endif    // _NUCLEAR_MODEL_I_H_
