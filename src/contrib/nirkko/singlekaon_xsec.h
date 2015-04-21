#include <iostream>
#include <iomanip>
#include <math.h>

// Class initialisation
class singlekaon_xsec {
  
  // Physics parameters
  double pi, amLam, am, amEta, aml, amSig, amk, ampi, Vus;
  
  // Input by interaction list generator
  int ilep, ik;
  
  // Input by the kinematics generator
  double Enu, Ekaon, pkvec;
  
  // SU(3) parameters, maybe in UserPhysicsOptions.xml
  double GeVtocm, fpi, d, f, g, amup, amun, Fm1, Fm2;
  
  // Output calculated by cross-section function
  double Elep, alepvec, aqvec, angkq, aq0;

public:

  // Threshold for given reaction
  double threshold;

  // Initialise cross-section calculation
  void init(double Etot, int type, int reac);
  
  // Calculate cross-section
  double diffxsec(double Tlep, double Tkaon, double theta, double phikq);
  
  // Calculate matrix elements
  double Amatrix_NN(double theta, double phikq);
  double Amatrix_NP(double theta, double phikq);
  double Amatrix_PP(double theta, double phikq);

};

