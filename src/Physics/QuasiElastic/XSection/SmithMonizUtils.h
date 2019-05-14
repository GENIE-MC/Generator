//____________________________________________________________________________
/*!

\class    genie::SmithMonizUtils

\brief    Contains auxiliary functions for Smith-Moniz model. \n

\ref      [1] R.A.Smith and E.J.Moniz, Nuclear Physics  B43, (1972) 605-622 \n
          [2] K.S. Kuzmin, V.V. Lyubushkin, V.A.Naumov Eur. Phys. J. C54, (2008) 517-538

\author   Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          adapted from  fortran code provided by \n
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>, Joint Institute for Nuclear Research \n
          Vladimir Lyubushkin, Joint Institute for Nuclear Research \n
          Vadim Naumov <vnaumov@theor.jinr.ru>, Joint Institute for Nuclear Research  \n
          based on code of \n
          Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>, University of Liverpool & STFC Rutherford Appleton Lab

\created  May 05, 2017

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SMITH_MONIZ_UTILS_H_
#define _SMITH_MONIZ_UTILS_H_

#include <TLorentzVector.h>

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"

namespace genie {



class SmithMonizUtils : public Algorithm {

public:
        SmithMonizUtils();
        SmithMonizUtils(string config);
        virtual ~SmithMonizUtils();
        void SetInteraction(const Interaction * i);
        double GetBindingEnergy(void) const;
        double GetFermiMomentum(void) const;
        double GetTheta_k(double v, double qv) const;
        double GetTheta_p(double pv, double v, double qv, double &E_p) const;
        double E_nu_thr_SM(void) const;
        Range1D_t Q2QES_SM_lim(void) const;
        Range1D_t vQES_SM_lim(double Q2) const;
        Range1D_t kFQES_SM_lim(double nu, double Q2) const;
        static double rho(double P_Fermi, double T_Fermi, double p);
        double PhaseSpaceVolume(KinePhaseSpace_t ps) const;
        
        //! methods overloading the Algorithm() interface implementation
        //! to build the fragmentation function from configuration data
        void Configure(const Registry & config);
        void Configure(string config);

private:
        template <class C>
        class Func1D
        {
            public: 
                Func1D(const C &obj, double (C::*f)(double) const):obj_(obj), f_(f){}
                ~Func1D(){}
                double operator()(double d) {return (obj_.*f_)( d);}
            private:
                const C &obj_;
                double (C::*f_)(double) const;
        };
        
        void   LoadConfig (void);
        double QEL_EnuMin_SM(double E_nu) const;
        double Q2lim1_SM(double Q2) const;
        double Q2lim2_SM(double Q2) const;
        double LambdaFUNCTION(double a, double b, double c) const;
        void DMINFC(Func1D<SmithMonizUtils> &F, double A,double B, double EPS, double DELTA, double &X, double &Y, bool &LLM) const;
        double vQES_SM_low_bound  (double Q2) const;
        double vQES_SM_upper_bound(double Q2) const;
        
        map<int, double> fNucRmvE;
        string fKFTable;
        bool fUseParametrization;
        
        const Interaction *  fInteraction;
        
        // Some often used variables of class.
        // To not calculate them again and again and for speed increase
        // they are initialized at once for multiple use
        double  E_nu;       ///<  Neutrino energy (GeV)
        double  m_lep;      ///<  Mass of final charged lepton (GeV)
        double  mm_lep;     ///<  Squared mass of final charged lepton (GeV)
        double  m_ini;      ///<  Mass of initial hadron or hadron system (GeV)
        double  mm_ini;     ///<  Sqared mass of initial hadron or hadron system (GeV)
        double  m_fin;      ///<  Mass of final hadron or hadron system (GeV)
        double  mm_fin;     ///<  Squared mass of final hadron or hadron system (GeV)
        double  m_tar;      ///<  Mass of target nucleus (GeV)
        double  mm_tar;     ///<  Squared mass of target nucleus (GeV)
        double  m_rnu;      ///<  Mass of residual nucleus (GeV)
        double  mm_rnu;     ///<  Squared mass of residual nucleus (GeV)
        double  P_Fermi;    ///<  Maximum value of Fermi momentum of target nucleon (GeV)
        double  E_BIN;      ///<  Binding energy (GeV)
mutable	double  Enu_in;     ///<  Running neutrino energy (GeV)

     
};
  

}       // genie namespace
#endif  // _SMITH_MONIZ_UTILS_H_

