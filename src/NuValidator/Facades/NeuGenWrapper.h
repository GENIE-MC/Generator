//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenWrapper

\brief    NeuGEN's C++ Wrapper class

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_WRAPPER_H_
#define _NEUGEN_WRAPPER_H_

#include <string>
#include <ostream>

#include "Facades/NeuGenConfig.h"
#include "Facades/NeuGenCuts.h"
#include "Facades/NGInitState.h"
#include "Facades/NGFlavor.h"
#include "Facades/NGNucleus.h"
#include "Facades/NGCcNc.h"
#include "Facades/NGInteraction.h"
#include "Facades/NGKineVar.h"
#include "Facades/XSecVsEnergy.h"

using std::string;
using std::ostream;

extern "C" {
  void set_parameters_         (char*, int *,float * val);
  void set_pdfset_             (int *, int * , float *);
  void makestate_              (int*,int*,int*,int*,int*,int*,int*,int*,int*);
  void writestate_             (int*);
  void set_default_parameters_ (void);
  void sig_value_              (float *, int *, int *, float *) ;
  void dsig_value_             (float *, int *, float *, int *, int *, float *) ;
  void set_kvcuts_             (int*,float*,float*);
  void set_masks_              (int *, bool *, bool *, bool *);
}

namespace genie   {
namespace nuvld   {
namespace facades {

class NeuGenWrapper 
{

public:

  NeuGenWrapper();
  NeuGenWrapper(const char * name);
  NeuGenWrapper(const NeuGenConfig * config);
  ~NeuGenWrapper();

  const char * Name (void) const { return _name.c_str(); }

  void SetDefaultConfig (void);
  void Reconfigure      (const NeuGenConfig * config);
  void SetCuts          (NeuGenCuts * cuts);
  void NoCuts           (void);
  
  float XSec (float e, NGInteraction * ni, NeuGenCuts * cuts=0);

  float ExclusiveXSec(float e, NGInteraction * ni, NGFinalState * final, NeuGenCuts * cuts=0);

  float DiffXSec(float e, NGKineVar_t  kid, float kval, NGInteraction *ni, NeuGenCuts * cuts=0);

  float ExclusiveDiffXSec(float e, NGKineVar_t  kid, float kval, NGInteraction *ni, NGFinalState * final,
                             NeuGenCuts * cuts=0);
                  
  XSecVsEnergy * XSecSpline (float emin,float emax, int nbins, NGInteraction * ni, NeuGenCuts * cuts=0);

  XSecVsEnergy * ExclusiveXSecSpline(float emin,float emax, int nbins, 
                          NGInteraction * ni, NGFinalState * final, NeuGenCuts * cuts=0);

  void Print(ostream & stream) const;
  
  friend ostream & operator << (ostream & stream, const NeuGenWrapper & conf);
  
private:

  void Init(void);
  
  string         _name;
  NeuGenConfig * _config;

};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif
