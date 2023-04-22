//____________________________________________________________________________
/*!

\class    genie::DCCSPPPXSec

\brief
Class calculate differental cross-sections
\f[
\frac{d^4\sigma}{dQ^2dWd\cos\theta_\pi d\phi_\pi}
\f]
or
\f[
\frac{d^2\sigma}{dQ^2dW}
\f]
for the specific neutrino energy (in the lab frame), where:


Variable             | Description
---------------------|-----------------------------------------------------
\f$W\f$              | The invariant mass
\f$Q^2\f$            | Minus of the square of the 4-momentum transfer
\f$\cos\theta_\pi\f$ | Cosine of the polar angle of the final pion in the \f$\pi\f$N rest frame
\f$\phi_\pi\f$       | The azimuthal angle of the final pion in the \f$\pi\f$N rest frame
for the following channels:
-#  \f$\ell           + p \to \ell            + p + \pi^0\f$
-#  \f$\ell           + p \to \ell            + n + \pi^+\f$
-#  \f$\ell           + n \to \ell            + n + \pi^0\f$
-#  \f$\ell           + n \to \ell            + p + \pi^-\f$
-#  \f$\nu            + p \to \ell^-          + p + \pi^+\f$
-#  \f$\nu            + n \to \ell^-          + p + \pi^0\f$
-#  \f$\nu            + n \to \ell^-          + n + \pi^+\f$
-#  \f$\overline{\nu} + n \to \ell^+          + n + \pi^-\f$
-#  \f$\overline{\nu} + p \to \ell^+          + n + \pi^0\f$
-#  \f$\overline{\nu} + p \to \ell^+          + p + \pi^-\f$
-#  \f$\nu            + p \to \nu             + p + \pi^0\f$
-#  \f$\nu            + p \to \nu             + n + \pi^+\f$
-#  \f$\nu            + n \to \nu             + n + \pi^0\f$
-#  \f$\nu            + n \to \nu             + p + \pi^-\f$
-#  \f$\overline{\nu} + p \to \overline{\nu}  + p + \pi^0\f$
-#  \f$\overline{\nu} + p \to \overline{\nu}  + n + \pi^+\f$
-#  \f$\overline{\nu} + n \to \overline{\nu}  + n + \pi^0\f$
-#  \f$\overline{\nu} + n \to \overline{\nu}  + p + \pi^-\f$



\ref
          1. https://www.phy.anl.gov/theory/research/anl-osaka-pwa/ (see also references on the website)
          2. https://www.rcnp.osaka-u.ac.jp/~anl-osk/               (see also references on the website)


\authors  Igor Kakorin <kakorin@jinr.ru>,
          Joint Institute for Nuclear Research,
          Bogoliubov Laboratory of Theoretical Physics \n
          adapted from fortran code provided by \n
          Toru Sato <tsato@rcnp.osaka-u.ac.jp>
          Osaka University, Department of Physics \n
          based on code of \n
          Costas Andreopoulos <constantinos.andreopoulos@cern.ch>
          University of Liverpool, Department of Physics

\created  Apr 12, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DCC_SPP_PXSEC_H_
#define _DCC_SPP_PXSEC_H_

#include <vector>
#include <complex>
//#include <memory>
#include <string>

#include "Framework/Interaction/SppChannel.h"
#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {


  class XSecIntegratorI;

  class DCCSPPPXSec: public XSecAlgorithmI {


    public:
        DCCSPPPXSec();
        DCCSPPPXSec(std::string config);
        virtual ~DCCSPPPXSec();

        // implement the XSecAlgorithmI interface
        double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
        double Integral     (const Interaction * i) const;
        bool   ValidProcess (const Interaction * i) const;
        bool   ValidKinematics(const Interaction * i) const;

        // overload the Algorithm::Configure() methods to load private data
        // members from configuration options
        void Configure(const Registry & config);
        void Configure(std::string config);

    private:
        /// Maximal number of \f$W\f$-nodes in the table with ANL-Osaka multipole amplitudes
        static const int maxW    = 110;
        /// Maximal number of \f$Q^2\f$-nodes in the table with ANL-Osaka multipole amplitudes
        static const int maxQ2   = 28;
        /// Maximal number of the multipoles \f$E_{l+}, E_{l-}, M_{l+}, M_{l-}, S_{l+}, S_{l-}, L_{l+}, L_{l-}\f$
        static const int maxmu   = 8;
        /// Maximal number of the orbital momentum in the table with ANL-Osaka multipole amplitudes, \f$l\in[0,\textrm{maxl}-1]\f$
        static const int maxl    = 6;
        /// Maximal number of isospin indices \f$a^V_{\frac{3}{2}}, a^V_{\frac{1}{2}}, a^V_{0}, a^A_{\frac{3}{2}}, a^A_{\frac{1}{2}}\f$
        static const int maxiso  = 5;
        /// Real and imaginary parts of complex number
        static const int maxzp   = 2;
        /// Maximal number of coefficients of cubic spline \f$S_i(x) = y_i+b_i(x-x_i)+{c_i}(x-x_i)^2+d_i(x-x_i)^3\f$
        static const int maxcf   = 4;


        void LoadConfig (void);
        /// Return vector part multipole
        std::complex<double> MultipoleV(const Interaction * interaction, int mult, int l) const;
        /// Return axial-vector part multipole
        std::complex<double> MultipoleA(const Interaction * interaction, int multipole, int l) const;
        /// Return the value of isospin amplitude
        double IsospinAmplitude(const Interaction * interaction, int iso) const;
        /// Return file with multupole amplitude table
        std::string FindDataTableFile(const std::string &basename, bool &ok) const;
        /// Read ANL-Osaka multipole amplitudes from file
        void ReadMultipoleTable(void);
        /// Precalculate cubic spline coefficients using W-nodes and table with ANL-Osaka multipole amplitudes to speed up computations
        void InitializeSplineCoefficients(void);
        /// Calculate position in table with ANL-Osaka multipole amplitudes
        int MultipoleTblIndx (int imu, int il, int iso, int izpart, int iq2, int icf, int iw) const;
        /// Calculate coefficients \f$b_i, c_i, d_i\f$ for cubic spline \f$S_i(x) = y_i + b_i(x - x_i) + {c_i}(x-x_i)^2 + {d_i}(x - x_i)^3\f$ for a function given in points \f$x_i, y_i\f$
        void Spline (int n, const double * x, const double * y, double * b, double * c, double * d) const;
        /// Return the nearest node to the value \f$W\f$ in the table DCCSPPPXSec::Wnodes
        int Wnode (double W) const;
        /// Return the nearest node to the value \f$Q^2\f$ in the table DCCSPPPXSec::Q2nodes
        int Q2node (double Q2) const;
        /// Calculate legendre polynomials and their first two derivatives
        void CalculateLegendre(double x, double l[3][maxl+1]) const;


        /// Table with complex values of ANL-Osaka multipole amplitudes: \f$E_{l\pm}(W, Q^2), M_{l\pm}(W, Q^2), S_{l\pm}(W, Q^2), L_{l\pm}(W, Q^2)\f$, the layout of table is determined by DCCSPPPXSec::MultipoleTblIndx
        std::vector<double> MultipoleTbl = std::vector<double>(maxmu*maxl*maxiso*maxzp*maxQ2*maxcf*maxW);
        /// Table with \f$W\f$-nodes
        std::vector<double> Wnodes = std::vector<double>(maxW);
        /// Table with \f$Q^2\f$-nodes
        std::vector<double> Q2nodes = std::vector<double>(maxQ2);
        /// Paths to check when searching for file with ANL-Osaka multipole amplitudes
        std::vector<std::string> fDataPaths;
        /// Table of Fermi momentum (kF) constants for various nuclei
        std::string fKFTable;
        /// Use parametrization for Fermi momentum instead of table?
        bool fUseRFGParametrization;
        /// Account for Pauli blocking?
        bool fUsePauliBlocking;
        /// External EM cross section scaling factor
        double fXSecScaleEM;
        /// External CC cross section scaling factor
        double fXSecScaleCC;
        /// External NC cross section scaling factor
        double fXSecScaleNC;
        /// sin^2(Weingberg angle)
        double   fSin2Wein;
        /// |Vud| (magnitude ud-element of CKM-matrix)
        double   fVud;

        const XSecIntegratorI * fXSecIntegrator;

  };


}       // genie namespace

#endif  // _DCC_SPP_PXSEC_H_
