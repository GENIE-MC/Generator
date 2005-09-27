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
#include "Facades/NGSF.h"
#include "Facades/NGInteraction.h"
#include "Facades/NGKineVar.h"
#include "Numerical/Spline.h"

using std::string;
using std::ostream;
using namespace genie;

extern "C" {

//-- the original NeuGEN function calls

  void print_configuration_      (void);
  void initialize_configuration_ (char*, int*, int*, bool*);
  void set_parameters_           (char*, int *,float *);
  void set_pdfset_               (int *, int * , float *);
  void makestate_                (int*,int*,int*,int*,int*,int*,int*,int*,int*);
  void writestate_               (int*);
  void set_default_parameters_   (void);
  void set_kvcuts_               (int*,float*,float*);
  void set_masks_                (int *, bool *, bool *, bool *);
  void sig_value_                (float *, int *, int *, float *) ;
  void dsig_value_               (float *, int *, float *, int *, int *, float *) ;
  void nu_structurefunctions_    (int *, int *, int *, float *, float *, int *, float *, float *, 
                                    float *, float *, float *, float *);
  void e_structurefunctions_     (int *, int *, float *, float *, int *, float *, float *);
  void ddsig_e_value_            (float *, int *, float *, int *, float *, int *, int *, float *);

  void gen_control_                   (char *, int *, bool *);
  void heplst_                        (int *);
  void generate_nu_event_             (int *, float *, int *, int *);
  void get_n_stdhep_entries_          (int *);
  void get_stdhep_particle_info_      (int *, int *, int *);
  void get_stdhep_particle_mothers_   (int *, int *, int *);
  void get_stdhep_particle_daughters_ (int *, int *, int *);
  void get_stdhep_particle_p4_        (int *, double *, double *, double *, double *);
  void get_stdhep_particle_v4_        (int *, double *, double *, double *, double *);
}

namespace genie   {

class EventRecord;
class InitialState;

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

  //-- methods for configuring and setting cuts
  
  void SetDefaultConfig (void);
  void Reconfigure      (const NeuGenConfig * config);
  void SetCuts          (NeuGenCuts * cuts);
  void NoCuts           (void);  

  //-- methods for getting NeuGEN v and e cross sections & structure functions
    
  float XSec
            (float e, NGInteraction * ni, NeuGenCuts * c=0);
  float ExclusiveXSec
            (float e, NGInteraction * ni, NGFinalState * f, NeuGenCuts * c=0);
  float DiffXSec
            (float e, NGKineVar_t  kv, float kval, NGInteraction *ni, NeuGenCuts * c=0);
  float ExclusiveDiffXSec
            (float e, NGKineVar_t  kv, float kval, NGInteraction *ni, NGFinalState * f, NeuGenCuts * c=0);
  float eDiff2Xsec
             (float e, NGKineVar_t  kv1, float kval1, NGKineVar_t  kv2, float kval2, NGInteraction * ni, NeuGenCuts * c=0);
  float StrucFunc
             (float x, float Q2, int A, NGInitState_t init, NGCcNc_t ccnc, int isf, int raw_dis);
  
  //-- methods for getting cross sections & structure functions predictions as splines

  Spline * XSecSpline
             (float emin,float emax, int nbins, NGInteraction * ni, NeuGenCuts * cuts=0);  
  Spline * ExclusiveXSecSpline
             (float emin,float emax, int nbins, NGInteraction * ni, NGFinalState * f, NeuGenCuts * c=0);                          
  Spline * StrucFuncSpline
            (NGKineVar_t xvar, float varmin, float varmax, int nbins,
             float fixedvar, int A, NGInitState_t init, NGCcNc_t ccnc,NGSF_t sf, int raw_dis = 2);
             
  //-- methods for generating an event and getting it as a GENIE event record

  void          GenControl    (char * name, int var);
  EventRecord * GenerateEvent (int nupdgc, float E, int A, int Z);
  void          PrintEvent    (void);

  //-- print state
                                  
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
