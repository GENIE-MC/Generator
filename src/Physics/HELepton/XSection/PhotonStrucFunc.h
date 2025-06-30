//____________________________________________________________________________
/*!

\class    genie::PhotonStrucFunc

\brief    Structure function using photon PDFs of nucleons

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PHOTON_STRUC_FUNC_H_
#define _PHOTON_STRUC_FUNC_H_

#include "Framework/Numerical/Spline.h"

#include <map>

using std::map;

namespace genie {

  struct SF_x {
    double Fp1,Fm1;
    double Fp2,Fm2;
    double Fp3,Fm3;
  };

  class PhotonStrucFunc
  {
    public:

      // ................................................................
      // GLRES form factor type
      //

      class PhotonStrucFuncTable 
      {
        public:
          PhotonStrucFuncTable() { }
          ~PhotonStrucFuncTable() { /* note: should delete the grids! */ }
          map< int, genie::Spline * > Table;
      };

      // ................................................................

      static PhotonStrucFunc * Instance(void);

      double EvalSF  ( int hitnuc, int hitlep, double x ) { return fSFTables[hitnuc].Table[hitlep]->Evaluate(x); }

    private:

      // Ctors & dtor
      PhotonStrucFunc();
      PhotonStrucFunc(const PhotonStrucFunc &);
     ~PhotonStrucFunc();

      // Self
      static PhotonStrucFunc * fgInstance;

      // These map holds all SF tables (interaction channel is the key)
      map<int, PhotonStrucFuncTable> fSFTables;


      // singleton cleaner
      struct Cleaner {
        void DummyMethodAndSilentCompiler(){}
          ~Cleaner(){
          if (PhotonStrucFunc::fgInstance !=0){
            delete PhotonStrucFunc::fgInstance;
            PhotonStrucFunc::fgInstance = 0;
          }
        }
      };
      friend struct Cleaner;
  };

} // genie namespace

#endif // _PHOTON_STRUC_FUNC_H_