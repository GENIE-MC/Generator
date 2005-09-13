//____________________________________________________________________________
/*!

\class   genie::NtpMCJobEnv

\brief   MINOS-style Ntuple Class to hold the MC Job enviroment (all GENIE
         environmental variables) in the output tree header

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 10, 2006

*/
//____________________________________________________________________________

#ifndef _NTP_MC_JOB_ENV_H_
#define _NTP_MC_JOB_ENV_H_

#include <ostream>

#include <TObject.h>
#include <TObjString.h>

using std::ostream;

namespace genie {

class NtpMCJobEnv : public TObject {

public :

  NtpMCJobEnv();
  NtpMCJobEnv(const NtpMCJobEnv & hdr);
  virtual ~NtpMCJobEnv();

  void Init (void);
  void Copy (const NtpMCJobEnv & hdr);

  void PrintToStream(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const NtpMCJobEnv & hdr);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breaking field data members not prefaced by "f" and mostly lowercase.
  TObjString gevgl;     ///< value of the GEVGL    env. var.
  TObjString gspload;   ///< value of the GSPLOAD  env. var.
  TObjString gspsave;   ///< value of the GSPSAVE  env. var.
  TObjString gmsgconf;  ///< value of the GMSGCONF env. var.

  ClassDef(NtpMCJobEnv, 1)
};

}      // genie namespace

#endif // _NTP_MC_JOB_ENV_H_
