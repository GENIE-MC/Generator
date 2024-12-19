//____________________________________________________________________________
/*!

\class    genie::HEDISStrcuFunc

\brief    Singleton class to load Structure Functions used in HEDIS.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_STRUC_FUNC_H_
#define _HEDIS_STRUC_FUNC_H_

#include "Framework/Numerical/BLI2D.h"
#include "Framework/Interaction/Interaction.h"

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>


using std::map;
using std::vector;
using std::string;
using std::to_string;

namespace genie {

  struct SF_info {

    string LHAPDFset;
    int    LHAPDFmember;
    bool   IsNLO;
    string Scheme;
    int    QrkThrs;
    int    NGridX;
    int    NGridQ2;
    double XGridMin;
    double Q2GridMin, Q2GridMax;
    double MassW, MassZ;
    double Rho, Sin2ThW;
    double Vud,Vus,Vub;
    double Vcd,Vcs,Vcb;
    double Vtd,Vts,Vtb;

  };

  inline bool operator== (const SF_info& a, const SF_info& b) {
    
    if ( a.LHAPDFmember != b.LHAPDFmember   ) return false;
    if ( a.LHAPDFset    != b.LHAPDFset      ) return false;
    if ( a.IsNLO        != b.IsNLO          ) return false;
    if ( a.Scheme       != b.Scheme         ) return false;
    if ( a.QrkThrs      != b.QrkThrs        ) return false;
    if ( a.NGridX       != b.NGridX         ) return false;
    if ( a.NGridQ2      != b.NGridQ2        ) return false;
    if ( abs(a.XGridMin-b.XGridMin)>1e-10   ) return false;
    if ( abs(a.Q2GridMin-b.Q2GridMin)>1e-10 ) return false;
    if ( abs(a.Q2GridMax-b.Q2GridMax)>1e-10 ) return false;
    if ( abs(a.MassW-b.MassW)>1e-10         ) return false;
    if ( abs(a.MassZ-b.MassZ)>1e-10         ) return false;
    if ( abs(a.Rho-b.Rho)>1e-10             ) return false;
    if ( abs(a.Sin2ThW-b.Sin2ThW)>1e-10     ) return false;
    if ( abs(a.Vud-b.Vud)>1e-10             ) return false;
    if ( abs(a.Vus-b.Vus)>1e-10             ) return false;
    if ( abs(a.Vub-b.Vub)>1e-10             ) return false;
    if ( abs(a.Vcd-b.Vcd)>1e-10             ) return false;
    if ( abs(a.Vcs-b.Vcs)>1e-10             ) return false;
    if ( abs(a.Vcb-b.Vcb)>1e-10             ) return false;
    if ( abs(a.Vtd-b.Vtd)>1e-10             ) return false;
    if ( abs(a.Vts-b.Vts)>1e-10             ) return false;
    if ( abs(a.Vtb-b.Vtb)>1e-10             ) return false;
    return true;

  }

  inline std::istream & operator>> (std::istream& is, SF_info& a) {

    string saux;
    std::getline (is,saux); //# Header
    std::getline (is,saux); //# Header
    std::getline (is,saux); //# Header
    std::getline (is,saux); //# Header
    std::getline (is,saux); //# LHAPDF set
    std::getline (is,saux); a.LHAPDFset=saux.c_str();
    std::getline (is,saux); //# LHAPDF member
    std::getline (is,saux); a.LHAPDFmember=atoi(saux.c_str());
    std::getline (is,saux); //#NLO
    std::getline (is,saux); a.IsNLO=atoi(saux.c_str());
    std::getline (is,saux); //# Mass scheme
    std::getline (is,saux); a.Scheme=saux.c_str();
    std::getline (is,saux); //# Quark threshold
    std::getline (is,saux); a.QrkThrs=atof(saux.c_str());
    std::getline (is,saux); //# NGridX
    std::getline (is,saux); a.NGridX=atoi(saux.c_str());
    std::getline (is,saux); //# NGridQ2
    std::getline (is,saux); a.NGridQ2=atoi(saux.c_str());
    std::getline (is,saux); //# XGridMin
    std::getline (is,saux); a.XGridMin=atof(saux.c_str());
    std::getline (is,saux); //# Q2min
    std::getline (is,saux); a.Q2GridMin=atof(saux.c_str());
    std::getline (is,saux); //# Q2max
    std::getline (is,saux); a.Q2GridMax=atof(saux.c_str());
    std::getline (is,saux); //# Mass W
    std::getline (is,saux); a.MassW=atof(saux.c_str());
    std::getline (is,saux); //# Mass Z
    std::getline (is,saux); a.MassZ=atof(saux.c_str());
    std::getline (is,saux); //# Rho
    std::getline (is,saux); a.Rho=atof(saux.c_str());
    std::getline (is,saux); //# Sin2ThW
    std::getline (is,saux); a.Sin2ThW=atof(saux.c_str());
    std::getline (is,saux); //# CKM
    std::getline (is,saux); a.Vud=atof(saux.c_str());
    std::getline (is,saux); a.Vus=atof(saux.c_str());
    std::getline (is,saux); a.Vub=atof(saux.c_str());
    std::getline (is,saux); a.Vcd=atof(saux.c_str());
    std::getline (is,saux); a.Vcs=atof(saux.c_str());
    std::getline (is,saux); a.Vcb=atof(saux.c_str());
    std::getline (is,saux); a.Vtd=atof(saux.c_str());
    std::getline (is,saux); a.Vts=atof(saux.c_str());
    std::getline (is,saux); a.Vtb=atof(saux.c_str());
    return is;

  }

  inline std::ostream & operator<< (std::ostream& os, const SF_info& a) {

    return os << '\n'
              << "#--------------------------------------------------------------------------------"    << '\n'
              << "# Metafile that stores information used to generate Structure Functions for HEDIS"   << '\n'
              << "#--------------------------------------------------------------------------------"    << '\n'
              << "# LHAPDF set"    << '\n'
              << a.LHAPDFset     << '\n'
              << "# LHAPDF member" << '\n'
              << a.LHAPDFmember  << '\n'
              << "# NLO"        << '\n' 
              << a.IsNLO         << '\n' 
              << "# Mass Scheme"       << '\n'
              << a.Scheme        << '\n'
              << "# Quark threshold"      << '\n'
              << a.QrkThrs       << '\n'
              << "# NX"       << '\n'
              << a.NGridX        << '\n'
              << "# NQ2"      << '\n'
              << a.NGridQ2       << '\n'
              << "# Xmin"     << '\n'
              << a.XGridMin      << '\n'
              << "# Q2min"    << '\n'
              << a.Q2GridMin     << '\n'
              << "# Q2max"    << '\n'
              << a.Q2GridMax     << '\n'
              << "# Mass W"        << '\n'
              << a.MassW         << '\n'
              << "# Mass Z"        << '\n'
              << a.MassZ         << '\n'
              << "# Rho"          << '\n'
              << a.Rho           << '\n'
              << "# Sin2ThW"      << '\n'
              << a.Sin2ThW       << '\n'
              << "# CKM"          << '\n'
              << a.Vud           << '\n'
              << a.Vus           << '\n'
              << a.Vub           << '\n'
              << a.Vcd           << '\n'
              << a.Vcs           << '\n'
              << a.Vcb           << '\n'
              << a.Vtd           << '\n'
              << a.Vts           << '\n'
              << a.Vtb           << '\n';

  }

  struct SF_xQ2 {

    double F1;
    double F2;
    double F3;

  };

  class HEDISStrucFunc
  {
    public:

      // ................................................................
      // HEDIS structure functions type
      //

      typedef enum StrucFuncType {  
        kMHTUndefined = 0,
        kSFT1, 
        kSFT2, 
        kSFT3, 
        kSFnumber, 
      } HEDISStrucFuncType_t;

      // ................................................................
      // HEDIS form factor type
      //

      class HEDISStrucFuncTable 
      {
        public:
          HEDISStrucFuncTable() { }
          ~HEDISStrucFuncTable() { /* note: should delete the grids! */ }
          map< HEDISStrucFunc::HEDISStrucFuncType_t, genie::BLI2DNonUnifGrid * > Table;
      };

      // ................................................................

      static HEDISStrucFunc * Instance(SF_info sfinfo);

      // method to return values of the SF for a particular channel in x and Q2
      SF_xQ2 EvalQrkSFLO  ( const Interaction * in, double x, double Q2 );
      SF_xQ2 EvalNucSFLO  ( const Interaction * in, double x, double Q2 ); 
      SF_xQ2 EvalNucSFNLO ( const Interaction * in, double x, double Q2 );

    private:

      // Ctors & dtor
      HEDISStrucFunc(SF_info sfinfo);
      HEDISStrucFunc(const HEDISStrucFunc &);
     ~HEDISStrucFunc();

      void CreateQrkSF    ( const Interaction * in, string sfFile );
      void CreateNucSF    ( const Interaction * in, string sfFile );

      string  QrkSFName ( const Interaction * in ); 
      string  NucSFName ( const Interaction * in ) ;
      int     QrkSFCode ( const Interaction * in ); 
      int     NucSFCode ( const Interaction * in ) ;

      // Self
      static HEDISStrucFunc * fgInstance;

      // These map holds all SF tables (interaction channel is the key)
      map<int, HEDISStrucFuncTable> fQrkSFLOTables;
      map<int, HEDISStrucFuncTable> fNucSFLOTables;
      map<int, HEDISStrucFuncTable> fNucSFNLOTables;

      SF_info fSF;
      vector<double> sf_x_array;
      vector<double> sf_q2_array;

      // singleton cleaner
      struct Cleaner {
        void DummyMethodAndSilentCompiler(){}
          ~Cleaner(){
          if (HEDISStrucFunc::fgInstance !=0){
            delete HEDISStrucFunc::fgInstance;
            HEDISStrucFunc::fgInstance = 0;
          }
        }
      };
      friend struct Cleaner;
  };

} // genie namespace

#endif // _HEDIS_STRUC_FUNC_H_
