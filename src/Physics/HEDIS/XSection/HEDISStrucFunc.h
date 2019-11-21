//____________________________________________________________________________
/*!

\class    genie::HEDISStrcuFunc

\brief    Singleton class to load Structure Functions used in HEDIS.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_STRUC_FUNC_H_
#define _HEDIS_STRUC_FUNC_H_

#include "Framework/Numerical/BLI2D.h"
#include "Framework/Interaction/HEDISChannel.h"

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
    
    bool ok = true;
    if ( a.LHAPDFmember != b.LHAPDFmember   ) { ok = false; std::cout << "LHAPDFmember" << std::endl; }
    if ( a.LHAPDFset    != b.LHAPDFset      ) { ok = false; std::cout << "LHAPDFset"    << std::endl; }
    if ( a.IsNLO        != b.IsNLO          ) { ok = false; std::cout << "IsNLO"        << std::endl; }
    if ( a.Scheme       != b.Scheme         ) { ok = false; std::cout << "Scheme"       << std::endl; }
    if ( a.QrkThrs      != b.QrkThrs        ) { ok = false; std::cout << "QrkThrs"      << std::endl; }
    if ( a.NGridX       != b.NGridX         ) { ok = false; std::cout << "NGridX"       << std::endl; }
    if ( a.NGridQ2      != b.NGridQ2        ) { ok = false; std::cout << "NGridQ2"      << std::endl; }
    if ( abs(a.XGridMin-b.XGridMin)>1e-10   ) { ok = false; std::cout << "XGridMin"     << std::endl; }
    if ( abs(a.Q2GridMin-b.Q2GridMin)>1e-10 ) { ok = false; std::cout << "Q2GridMin"    << std::endl; }
    if ( abs(a.Q2GridMax-b.Q2GridMax)>1e-10 ) { ok = false; std::cout << "Q2GridMax"    << std::endl; }
    if ( abs(a.MassW-b.MassW)>1e-10         ) { ok = false; std::cout << "MassW"        << std::endl; }
    if ( abs(a.MassZ-b.MassZ)>1e-10         ) { ok = false; std::cout << "MassZ"        << std::endl; }
    if ( abs(a.Rho-b.Rho)>1e-10             ) { ok = false; std::cout << "Rho"          << std::endl; }
    if ( abs(a.Sin2ThW-b.Sin2ThW)>1e-10     ) { ok = false; std::cout << "Sin2ThW"      << std::endl; }
    if ( abs(a.Vud-b.Vud)>1e-10             ) { ok = false; std::cout << "Vud"          << std::endl; }
    if ( abs(a.Vus-b.Vus)>1e-10             ) { ok = false; std::cout << "Vus"          << std::endl; }
    if ( abs(a.Vub-b.Vub)>1e-10             ) { ok = false; std::cout << "Vub"          << std::endl; }
    if ( abs(a.Vcd-b.Vcd)>1e-10             ) { ok = false; std::cout << "Vcd"          << std::endl; }
    if ( abs(a.Vcs-b.Vcs)>1e-10             ) { ok = false; std::cout << "Vcs"          << std::endl; }
    if ( abs(a.Vcb-b.Vcb)>1e-10             ) { ok = false; std::cout << "Vcb"          << std::endl; }
    if ( abs(a.Vtd-b.Vtd)>1e-10             ) { ok = false; std::cout << "Vtd"          << std::endl; }
    if ( abs(a.Vts-b.Vts)>1e-10             ) { ok = false; std::cout << "Vts"          << std::endl; }
    if ( abs(a.Vtb-b.Vtb)>1e-10             ) { ok = false; std::cout << "Vtb"          << std::endl; }
    return ok;

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

      typedef enum HEDISStrucFuncType {  
        kMHTUndefined = 0,
        kSFT1, 
        kSFT2, 
        kSFT3, 
        kSFnumber, 
      } 
      HEDISStrucFuncType_t;

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

      static HEDISStrucFunc * Instance(string basedir, SF_info sfinfo);

      void CreateQrkSF    ( HEDISQrkChannel_t ch, string sfFile );
      void CreateNucSF    ( HEDISNucChannel_t ch, string sfFile );

      // method to return values of the SF for a particular channel in x and Q2
      SF_xQ2 EvalQrkSFLO  ( HEDISQrkChannel_t ch, double x, double Q2 );
      SF_xQ2 EvalNucSFLO  ( HEDISNucChannel_t ch, double x, double Q2 ); 
      SF_xQ2 EvalNucSFNLO ( HEDISNucChannel_t ch, double x, double Q2 );

    private:

      // Ctors & dtor
      HEDISStrucFunc(string basedir, SF_info sfinfo);
      HEDISStrucFunc(const HEDISStrucFunc &);
     ~HEDISStrucFunc();

      // Self
      static HEDISStrucFunc * fgInstance;

      // These map holds all SF tables (interaction channel is the key)
      map<HEDISQrkChannel_t, HEDISStrucFuncTable> fQrkSFLOTables;
      map<HEDISNucChannel_t, HEDISStrucFuncTable> fNucSFLOTables;
      map<HEDISNucChannel_t, HEDISStrucFuncTable> fNucSFNLOTables;

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
