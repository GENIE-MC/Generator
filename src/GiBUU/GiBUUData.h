//____________________________________________________________________________
/*!

\class    genie::GiBUUData

\brief    Singleton to load and serve data tables provided by the GiBUU group

\ref      http://gibuu.physik.uni-giessen.de/GiBUU
          Specific references for each piece of data included in given below.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 30, 2009

*/
//____________________________________________________________________________

#ifndef _GIBUU_DATA_H_
#define _GIBUU_DATA_H_

#include "BaryonResonance/BaryonResonance.h"
#include "Interaction/InteractionType.h"

namespace genie {

class Spline;

class GiBUUData
{
public:

  class FormFactors;

  // access GiBUUData singleton instance
  static GiBUUData * Instance (void);

  // access form factor data
  const FormFactors & FF(void) const;

  //
  // Resonance form factors. 
  // Details given in Phys. Rev. C 79, 034601 (2009).
  //
  class FormFactors {

   public:

     FormFactors();
    ~FormFactors();

     // the following is non-zero for I=1/2 (N) resonances
     double C3V (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const;
     double C4V (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 
     double C5V (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 
     double C6V (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 
     double C3A (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 
     double C4A (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 
     double C5A (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 
     double C6A (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const;
  
     // the following is non-zero for I=3/2 (Delta) resonances
     double F1V (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 
     double F2V (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 
     double FA  (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const;  
     double FP  (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it) const; 

     double Q2min (void) const { return fMinQ2; }
     double Q2max (void) const { return fMaxQ2; }

   private:

     static double fMinQ2; ///< min Q2 for which resonance f/f data are given
     static double fMaxQ2; ///< max Q2 for which resonance f/f data are given

     // The first array index is the resonance id. 
     // Tina provided GiBUU form factor data for 13 resonances given below along with
     // the corresponding GENIE resonance ids
     //                 GENIE Resonance_t   as integer
     //   P33(1232) ->  kP33_1232            0
     //   S11(1535) ->  kS11_1535            1
     //   D13(1520) ->  kD13_1520            2
     //   S11(1650) ->  kS11_1650            3
     //   D15(1675) ->  kD15_1675            5
     //   S31(1620) ->  kS31_1620            6
     //   D33(1700) ->  kD33_1700            7
     //   P11(1440) ->  kP11_1440            8
     //   P13(1720) ->  kP13_1720           10
     //   F15(1680) ->  kF15_1680           11
     //   P31(1910) ->  kP31_1910           12
     //   F35(1905) ->  kF35_1905           14
     //   F37(1950) ->  kF37_1950           15
     // The remaining 3 array indices are:
     //     0  1  2  0 1    0   1  2  3   4   5   6   7   8   9  10  11
     //   [CC,NC,EM][n,p][F1V,F2V,FA,FP,C3V,C4V,C5V,C6V,C3A,C4A,C5A,C6A]
     //                  |-------------|                                  for I=1/2 resonances
     //                                |-------------------------------|  for I=3/2 resonances

     static const int kNRes    = 18;
     static const int kNCurr   =  3;
     static const int kNHitNuc =  2;
     static const int kNFFRes  = 12;

     //! actual form factor data = f(Q2)
     Spline * fFFRes [kNRes][kNCurr][kNHitNuc][kNFFRes];   

     //! func to retrieve interpolated form factor values
     double FFRes (double Q2, Resonance_t res, int nucleon_pdg, InteractionType_t it, int ffid) const;

     //! load all form factor data tables
     void LoadTables(void); 

     friend class GiBUUData;

  }; // FormFactors nested class

private:
  GiBUUData();
  GiBUUData(const GiBUUData & gibuu_data);
 ~GiBUUData();

  // load all data tables
  void LoadTables(void); 

  // singleton 'self'
  static GiBUUData * fInstance;

  // form factor data
  FormFactors * fFormFactors;

  // singleton cleaner
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (GiBUUData::fInstance !=0) {
            delete GiBUUData::fInstance;
            GiBUUData::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;

}; // GiBUUData class

}      // genie namespace

#endif // _GIBUU_DATA_H_

