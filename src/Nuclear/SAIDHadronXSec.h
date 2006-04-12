//____________________________________________________________________________
/*!

\class    genie::SAIDHadronXSec

\brief    Singleton class to load & serve SAID cross section splines.
          Data provided by S.Dytmnan.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 12, 2006

*/
//____________________________________________________________________________

#ifndef _SAID_HADRON_CROSS_SECTIONS_H_
#define _SAID_HADRON_CROSS_SECTIONS_H_

namespace genie {

class Spline;

class SAIDHadronXSec
{
public:
  static SAIDHadronXSec * Instance (void);

  //! cross section splines

  const Spline * PiplusPElasXSecSpl  (void) const { return fSplPiplusPElas; }
  const Spline * PiplusPInelXSecSpl  (void) const { return fSplPiplusPInel; }
  const Spline * PiminusPElasXSecSpl (void) const { return fSplPiminusPElas;}
  const Spline * PiminusPInelXSecSpl (void) const { return fSplPiminusPInel;}
  const Spline * PiminusPCExXSecSpl  (void) const { return fSplPiminusPCEx; }
  const Spline * PiplusAbsXSecSpl    (void) const { return fSplPiplusAbs;   }
  const Spline * NPElasXSecSpl       (void) const { return fSplNPElas;      }
  const Spline * NPInelXSecSpl       (void) const { return fSplNPInel;      }
  const Spline * PPElasXSecSpl       (void) const { return fSplPPElas;      }
  const Spline * PPInelXSecSpl       (void) const { return fSplPPInel;      }

private:
  SAIDHadronXSec();
  SAIDHadronXSec(const SAIDHadronXSec & shx);
  ~SAIDHadronXSec();

  void LoadTables(void);

  static SAIDHadronXSec * fInstance;

  //! cross section splines

  Spline * fSplPiplusPElas;
  Spline * fSplPiplusPInel;
  Spline * fSplPiminusPElas;
  Spline * fSplPiminusPInel;
  Spline * fSplPiminusPCEx; 
  Spline * fSplPiplusAbs;    
  Spline * fSplNPElas;      
  Spline * fSplNPInel;      
  Spline * fSplPPElas;      
  Spline * fSplPPInel;      

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (SAIDHadronXSec::fInstance !=0) {
            delete SAIDHadronXSec::fInstance;
            SAIDHadronXSec::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif //_SAID_HADRON_CROSS_SECTIONS_H_

