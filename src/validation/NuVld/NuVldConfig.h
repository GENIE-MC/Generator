//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldConfig

\brief    NuValidator GUI configuration options

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  September 17, 2005
*/
//_____________________________________________________________________________

#ifndef _NUVLD_CONFIG_H_
#define _NUVLD_CONFIG_H_

namespace genie {
namespace nuvld {

class NuVldConfig {

public:

  NuVldConfig();
  NuVldConfig(const NuVldConfig & config);
  ~NuVldConfig();

  bool UseNeuGEN        (void) const { return fUseNeuGEN;        }
  bool UseCompactLayout (void) const { return fUseCompactLayout; }

  void AutoDetect (void);
  void Copy       (const NuVldConfig & config);

  void SetUseNeuGEN        (bool use) { fUseNeuGEN        = use; }
  void SetUseCompactLayout (bool use) { fUseCompactLayout = use; }

private:

  void Init();

  bool fUseNeuGEN;
  bool fUseCompactLayout;
};

} // nuvld namespace
} // genie namespace

#endif

