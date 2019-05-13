//____________________________________________________________________________
/*!

\class    genie::IBDXSecMap

\brief    Maps specific nuclei to appropriate cross section models.

\author   Corey Reed <cjreed \at nikhef.nl>
          Nikhef

\created  February 4, 2010

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _IBD_CROSSSECTION_MAP_H_
#define _IBD_CROSSSECTION_MAP_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class IBDXSecMap : public XSecAlgorithmI {
      
   public:
      IBDXSecMap();
      IBDXSecMap(string config);
      virtual ~IBDXSecMap();
      
      //-- XSecAlgorithmI interface implementation
      double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral        (const Interaction * i) const;
      bool   ValidProcess    (const Interaction * i) const;
      bool   ValidKinematics (const Interaction * i) const;

      //-- overload the Algorithm::Configure() methods to load private data
      //   members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

   protected:
      const XSecAlgorithmI* SelectModel(const Target & t) const;
      
   private:
      void LoadConfig (void);
      
      bool                                 fIsotopesUseSameModel; //  if true, force all nuclei with same Z to use the same cross section model
      const XSecAlgorithmI*                fDefaultModel;         //  the default xsec model (should "work" with any nucleus)
      std::map<int, const XSecAlgorithmI*> fRefinedModels;        //  specific models for the given nucleus (int being the PDG code of the nucleus)
      
   public:
      ClassDef(IBDXSecMap, 1) // maps a request for the xsec on a given nucleus to the most appropriate model for that nucleus

};

};

#endif // _IBD_CROSSSECTION_MAP_H_

