//____________________________________________________________________________
/*!

\class   genie::rew::GReWeightINukeParams

\brief   Helper class for cross section model reweighting

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Jim Dobson <J.Dobson07 \at imperial.ac.uk>
         Imperial College London

\created Sep 10, 2009

*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INTRANUKE_PARAMS_H_
#define _G_REWEIGHT_INTRANUKE_PARAMS_H_

#include <map>

#include "PDG/PDGUtils.h"
#include "ReWeight/GSyst.h"

using std::map;

class TLorentzVector;

namespace genie {
namespace rew   {

 class GReWeightINukeParams {

 public:

   typedef enum EHadronType {
      kRwINukeUndefined = 0,
      kRwINukePion,
      kRwINukeNucl,
   } HadronType_t;

   static HadronType_t HadronTypeFromPdg(int pdgc) {
     if(pdg::IsPion   (pdgc)) return kRwINukePion;
     if(pdg::IsNucleon(pdgc)) return kRwINukeNucl;
     return kRwINukeUndefined;
   }

   GReWeightINukeParams();
  ~GReWeightINukeParams();

   class Fates;
   class MFP;

   Fates * FateParams         (int pdgc) const;
   MFP *   MeanFreePathParams (int pdgc) const;
   void    Reset              (void);
   void    Reconfigure        (void);
   void    SetCurTwkDial      (GSyst_t syst, double val);

   //.........................................................................
   //
   // nested class: Fates
   //
   //.........................................................................

   class Fates {
   public :
     Fates(HadronType_t hadtype = kRwINukeUndefined);
    ~Fates();

     double CurTwkDial    (GSyst_t s, double KE=-1.) const; ///< 
     bool   IsIncluded    (GSyst_t s) const;                ///< is included?
     bool   IsCushionTerm (GSyst_t s) const;                ///< is it a cushion term?
     bool   IsTweaked     (GSyst_t s) const;                ///< is included & tweaked to non-def value?
     bool   IsTweaked     (void) const;                     ///< is any param tweaked
     void   Reset         (void);                           ///<
     void   Reconfigure   (void);                           ///<
     double ChisqPenalty  (void) const;                     ///<
     void   SetCurTwkDial (GSyst_t s, double val);          ///< 

   private:    

     bool   IsHandled       (GSyst_t s) const;
     void   AddCushionTerms (void);
     bool   CheckUnitarity  (int n_points = 400);

     HadronType_t         fHadType;       ///<
     int                  fNSysts;        ///<
     int                  fNCushionTerms; ///<
     double               fCushAvgChisq;  ///<
     map<GSyst_t, double> fSystListMap;   ///< List of systematics included. 
     map<GSyst_t, bool>   fIsCushionMap;  ///< cushion term flag

     //JIMTODO - need to add bool fChi2AvgUpToDate which is set to false if 
     //SetTwkDial is called and is set to true if CheckFatesUnity is called. 
     //Also add conditional to GetAverageChi2 

   }; // Fates nested class


   //.........................................................................
   //
   // nested class: MFP
   //
   //.........................................................................

   class MFP {
   public :
     MFP(HadronType_t hadtype = kRwINukeUndefined);
    ~MFP();

     double ScaleFactor   (void) const;  ///< mean free path scale factor = 1 + twk_dial * frac_err
     double CurTwkDial    (void) const;  ///< current value of mfp tweak dial
     bool   IsIncluded    (void) const;  
     bool   IsTweaked     (void) const; 
     double ChisqPenalty  (void) const;  ///<
     void   Reset         (void);
     void   SetCurTwkDial (double val);

   private:    
     HadronType_t fHadType;     ///<
     GSyst_t      fSyst;        ///<
     double       fTwkDial;     ///<
     bool         fIsIncluded;  ///<

   }; // MFP nested class


 private:

    void Init(void);

    Fates * fParmPionFates;
    Fates * fParmNuclFates; 
    MFP *   fParmPionMFP;   
    MFP *   fParmNuclMFP;

 }; //GReWeightINukeParams

}      // rew   namespace
}      // genie namespace

#endif // _G_REWEIGHT_INTRANUKE_PARAMS_H_


