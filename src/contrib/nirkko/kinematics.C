// *********************************************************************
//  Kinematics of single kaon production - Martti Nirkko (28th Nov 2014)
//  Compile and run in terminal:     root -l -q kinematics.C+g
// *********************************************************************

#include "code/singlekaon_xsec.cxx"
#include <TMath.h>
#include <TVector3.h>

void kinematics() {
  // Validate kinematics of single kaon production
  // Have: T_l, T_K, cos(theta_l), phi_l, phi_Kq (?)
  // Want: p_l, p_K, p_N (3-vectors)
  
  // Specify amount of output
  const int VERBOSE = 0;    // verbosity (0/1)
  const int DEBUG   = 0;    // debugging mode
  
  // Use random input variables (0: use fixed inputs to verify GENIE results)
  const int RANDOM  = 1;
  
  // Initialise random number seed
  srand (time(NULL));
  const double pi = TMath::Pi();
  
  const int type = 2;       // lepton:   1=electron, 2=muon, 3=tau
  const int reac = 3;       // reaction: 1=NN, 2=NP, 3=PP
  
  std::string strl;
  double ml = 0.;                     // lepton mass [GeV]
  if      (type==1) { ml = 0.510998928e-3; strl = "e"; }
  else if (type==2) { ml = 0.1056583715;   strl = "mu"; }
  else if (type==3) { ml = 1.77682;        strl = "tau"; }
  else {std::cout<<"ERROR: Invalid lepton type!"<<std::endl; return;}
  
  std::string strK;
  double mK = 0.;                    // kaon mass [GeV]
  if      (reac==1) { mK = 0.493677; strK = "K+"; }
  else if (reac==2) { mK = 0.497614; strK = "K0"; }
  else if (reac==3) { mK = 0.493677; strK = "K+"; }
  else {std::cout<<"ERROR: Invalid reaction!"<<std::endl; return;}
  
  std::string strN0, strN1;
  double mN = 0.;                    // nucleon mass [GeV]
  if      (reac==1) { mN = 0.939565379; strN0 = "n"; strN1 = "n"; }
  else if (reac==2) { mN = 0.939565379; strN0 = "n"; strN1 = "p"; }
  else if (reac==3) { mN = 0.938272046; strN0 = "p"; strN1 = "p"; }
  else {std::cout<<"ERROR: Invalid reaction!"<<std::endl; return;}
  
  std::string strReac = Form( "nu_%s + %s -> %s + %s + %s", strl.c_str(), strN0.c_str(),
                               strl.c_str(), strN1.c_str(), strK.c_str() );
  
  // Reaction threshold
  const double threshold = ((ml+mN+mK)*(ml+mN+mK)-mN*mN)/(2.0*mN);
  const double Emax = 3.1;
  
  // Neutrino energy [GeV] - chosen at random between threshold and Emax
  double Enu;
  if (RANDOM)   Enu = threshold + (Emax-threshold)*(rand()/(double)RAND_MAX);
  else          Enu = 1.0;
  
  // Initialize reaction
  singlekaon_xsec *inst = new singlekaon_xsec();
  inst->init(Enu, type, reac);
  
  // INPUT PARAMETERS
  double TK_max, TK, Tl_max, Tl, costh_l, phi_l, phi_Kq, xsec;
  int evt;
  if (RANDOM) {
    xsec = -999.;
    while(xsec <= 0) {
      TK_max  = Enu - mK - ml;                            // maximal allowed kaon kinetic energy
      TK      = TK_max*(rand()/(double)RAND_MAX);         // kaon kinetic energy [0, TK_max]
      Tl_max  = Enu - mK - ml - TK;                       // maximal allowed lepton kinetic energy
      Tl      = Tl_max*(rand()/(double)RAND_MAX);         // lepton kinetic energy [0, Tl_max]
      costh_l = 2.*(rand()/(double)RAND_MAX)-1.;          // lepton polar angle [-1, 1]
      phi_l   = pi*(2.*(rand()/(double)RAND_MAX)-1.);     // lepton azimuthal angle [-pi, pi]
      phi_Kq  = pi*(2.*(rand()/(double)RAND_MAX)-1.);     // kaon azimuthal angle [-pi, pi]
      xsec = inst->diffxsec(Tl,TK,acos(costh_l),phi_Kq);  // 4D-differential cross-section
    }
  } else {
    xsec = -999.;
    evt = 0;                // choose event from input table
    double input[9][5] = { 
      {0.188312,   0.108214,   0.706911,   -0.66872,   5.75773},
      {0.234579,   0.0893942,  0.74982,     3.035,     6.10204},
      {0.196701,   0.0362663,  0.663934,    2.64581,   4.82931},
      {0.0191452,  0.0575003,  0.67804,    -2.96401,   4.11426},
      {0.0529443,  0.0119084,  0.127732,   -1.57085,   1.10203},
      {0.139563,   0.143596,   0.946546,    0.80543,   2.34818},
      {0.200567,   0.13117,    0.851108,   -0.806353,  1.98515},
      {0.0372796,  0.0501203,  0.148006,    0.81699,   4.35815},
      {0.0806993,  0.172774,   0.903585,    2.42196,   3.6006}
    };
    
    TK      = input[evt][0];
    Tl      = input[evt][1];
    costh_l = input[evt][2];
    phi_l   = input[evt][3];
    phi_Kq  = input[evt][4];
    xsec = inst->diffxsec(Tl,TK,acos(costh_l),phi_Kq);
    if (xsec <= 0) {
      printf("ERROR: INVALID XSEC! (%lf)\n", xsec);
      return;
    }
  }
  
  // Print input parameters to screen
  printf("\n\n");
  printf("INITIAL INPUT PARAMETERS\n");
  printf("------------------------\n");
  printf("Neutrino energy:\t\tE_nu\t= %6.3lf GeV\n",Enu);
  printf("Lepton kinetic energy:\t\tT_l\t= %6.3lf GeV\n",Tl);
  printf("Kaon kinetic energy:\t\tT_K\t= %6.3lf GeV\n",TK);
  printf("Lepton polar angle:\t\tcosth_l\t= %6.3lf\n",costh_l);
  printf("Lepton azimuthal angle:\t\tphi_l\t= %6.3lf pi\n",phi_l/pi);
  printf("Kaon azimuthal angle:\t\tphi_Kq\t= %6.3lf pi\n",phi_Kq/pi);
  printf("Resulting cross-section:\td4xsec\t= %10.3e [units]\n",xsec);
  printf("\n\n");
  
  // Some directly obtainable parameters
  const double theta_l = acos(costh_l);
  const double TN = Enu - Tl - TK - ml - mK;
  const double El = Tl + ml;
  const double EK = TK + mK;
  const double EN = TN + mN;
  
  // Momentum 3-vectors
  TVector3 p_l, p_Q, p_K, p_N, p_Kq, p_Nq;
  
  // Unit 3-vectors
  TVector3 e1(1,0,0);
  TVector3 e2(0,1,0);
  TVector3 e3(0,0,1);
  
  double pl, pK, pN;    // Total momentum
  double bl, bK, bN;    // Lorentz beta
  double gl, gK, gN;    // Lorentz gamma
  
  // Lepton momentum
  gl = Tl/ml + 1;
  bl = sqrt(1.-1./gl/gl);
  pl = bl*gl*ml;
  
  // Kaon momentum
  gK = TK/mK + 1;
  bK = sqrt(1.-1./gK/gK);
  pK = bK*gK*mK;
  
  // Nucleon momentum
  gN = TN/mN + 1;
  bN = sqrt(1.-1./gN/gN);
  pN = bN*gN*mN;
  
  // Calculate 3-momentum of lepton
  p_l[0] = pl*sin(theta_l)*cos(phi_l);
  p_l[1] = pl*sin(theta_l)*sin(phi_l);
  p_l[2] = pl*costh_l;
  if (VERBOSE) {
    printf("Lepton 3-mom in lab:\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_l[0], p_l[1], p_l[2], pl);
  }
  
  // Calculate 3-momentum of momentum transfer to hadronic system
  p_Q = -p_l;
  p_Q[2] += Enu;
  double pQ = p_Q.Mag();
  if (VERBOSE) {
    printf("Q-squared 3-mom in lab:\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_Q[0], p_Q[1], p_Q[2], pQ);
    printf("\nQ-squared 3-mom in Q2:\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           0., 0., pQ, pQ);
  }
  TVector3 dirQ = p_Q.Unit();
  
  // In the momentum transfer plane, new angles... (see eq.17 of notes)
  double costh_Kq = (pQ*pQ + pK*pK + mN*mN - (Enu-El-EK+mN)*(Enu-El-EK+mN))/(2.*pQ*pK);
  double theta_Kq = acos(costh_Kq);
  
  // Get rotation angle from z-axis to Q^2-direction
  double rot_Z = acos(e3*(p_Q.Unit()));         // rotate counter-clockwise, towards x-axis by [0,pi]
  
  // Get rotation angle from x-axis to Q^2-direction (projected to transverse plane)
  int sign = 0;
  if (p_Q[1] > 0) sign = +1;
  else            sign = -1;
  TVector3 p_Qt(p_Q[0],p_Q[1],0);
  double rot_X = sign*acos(e1*(p_Qt.Unit()));   // rotate around z-axis by [-pi,pi]
  
  // Make sure that x' is in the nu-q plane
  TVector3 p_nu = Enu*e3;
  TVector3 yrot = (p_Q.Cross(p_nu)).Unit();     // must be parallel to y' axis
  p_nu.RotateZ(-rot_X);
  p_nu.RotateY(-rot_Z);
  if (VERBOSE) {
    printf("Neutrino 3-mom in Q2:\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c",
           p_nu[0], p_nu[1], p_nu[2], p_nu.Mag());
    if ( abs(p_nu[1])>1.e-6 )
      printf("--> WARNING: NOT IN X'Z'-PLANE!\n");
    else
      printf("\n");
    printf("\n");
  }
  
  // Make sure that q points the same way as the y'-axis
  TVector3 unu = p_nu.Unit();
  TVector3 uQ  = p_Q.Unit();
  uQ.RotateZ(-rot_X);
  uQ.RotateY(-rot_Z);
  if (DEBUG) {
    printf("Unit vector nu in Q2:\t(%6.3lf | %6.3lf | %6.3lf )\n",unu[0],unu[1],unu[2]);
    printf("Unit vector q  in Q2:\t(%6.3lf | %6.3lf | %6.3lf )",uQ[0],uQ[1],uQ[2]);
    if ( abs(uQ[0])>1.e-6 || abs(uQ[1])>1.e-6 || abs(uQ[2])<(1-1.e-6) )
      printf(" --> WARNING: NOT EQUAL TO Z' UNIT VECTOR!\n");
    else
      printf("\n");
    printf("\n");
  }
  
  // Check what happens to the unit vectors when rotated (for debugging purposes)
  if (DEBUG) {
    TVector3 u1=e1, u2=e2, u3=e3;     // unit vectors
    u1.RotateZ(-rot_X); u2.RotateZ(-rot_X); u3.RotateZ(-rot_X);
    u1.RotateY(-rot_Z); u2.RotateY(-rot_Z); u3.RotateY(-rot_Z);
    printf("Unit vector X in Q2:\t(%6.3lf | %6.3lf | %6.3lf )\n",u1[0],u1[1],u1[2]);
    printf("Unit vector Y in Q2:\t(%6.3lf | %6.3lf | %6.3lf )\n",u2[0],u2[1],u2[2]);
    printf("Unit vector Z in Q2:\t(%6.3lf | %6.3lf | %6.3lf )\n",u3[0],u3[1],u3[2]);
    printf("\n");
    u1=e1; u2=e2; u3=e3;
    
    u1.RotateY(rot_Z); u2.RotateY(rot_Z); u3.RotateY(rot_Z);
    u1.RotateZ(rot_X); u2.RotateZ(rot_X); u3.RotateZ(rot_X);
    printf("Unit vector X' in lab:\t(%6.3lf | %6.3lf | %6.3lf )\n",u1[0],u1[1],u1[2]);
    printf("Unit vector Y' in lab:\t(%6.3lf | %6.3lf | %6.3lf )\n",u2[0],u2[1],u2[2]);
    printf("Unit vector Z' in lab:\t(%6.3lf | %6.3lf | %6.3lf )\n",u3[0],u3[1],u3[2]);
    printf("\n");
    printf("CHECK yrot in lab:\t(%6.3lf | %6.3lf | %6.3lf )\n",yrot[0],yrot[1],yrot[2]);
    printf("CHECK q1*x2' = q2*x1':\t%6.3lf = %6.3lf\n", p_Q[0]*u1[1], p_Q[1]*u1[0]);
    printf("\n");
  }
  
  // Calculate 3-momentum of kaon (in Q-frame)
  p_Kq[0] = pK*sin(theta_Kq)*cos(phi_Kq);
  p_Kq[1] = pK*sin(theta_Kq)*sin(phi_Kq);
  p_Kq[2] = pK*costh_Kq;
  
  // Calculate 3-momentum of nucleon (in Q-frame)
  p_Nq = -p_Kq;
  p_Nq[2] += pQ;
  if (VERBOSE) {
    printf("Kaon 3-mom in Q2:\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_Kq[0], p_Kq[1], p_Kq[2], p_Kq.Mag());
    printf("Nucleon 3-mom in Q2:\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n\n",
           p_Nq[0], p_Nq[1], p_Nq[2], p_Nq.Mag());
  }
  
  // Rotate particles into correct frame
  p_K = p_Kq;
  p_K.RotateY(rot_Z);
  p_K.RotateZ(rot_X);
  p_N = p_Nq;
  p_N.RotateY(rot_Z);
  p_N.RotateZ(rot_X);
  if (VERBOSE) {
    printf("Kaon 3-mom in lab:\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_K[0], p_K[1], p_K[2], p_K.Mag());
    printf("Nucleon 3-mom in lab:\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_N[0], p_N[1], p_N[2], p_N.Mag());
    printf("\n\n");
  }
  
  // Compare alternative rotation by Chris
  if (DEBUG) {
    TVector3 p_K2 = p_Kq;
    p_K2.RotateUz( p_Q.Unit() );
    TVector3 p_N2 = p_Nq;
    p_N2.RotateUz( p_Q.Unit() );
    printf("Kaon 3-mom in Q2:\t\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_Kq[0], p_Kq[1], p_Kq[2], p_Kq.Mag());
    printf("Kaon 3-mom in lab (Martti):\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_K[0], p_K[1], p_K[2], p_K.Mag());
    printf("Kaon 3-mom in lab (Chris):\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_K2[0], p_K2[1], p_K2[2], p_K2.Mag());
    printf("\n");
    printf("Nucleon 3-mom in Q2:\t\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_Nq[0], p_Nq[1], p_Nq[2], p_Nq.Mag());
    printf("Nucleon 3-mom in lab (Martti):\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_N[0], p_N[1], p_N[2], p_N.Mag());
    printf("Nucleon 3-mom in lab (Chris):\t(%6.3lf | %6.3lf | %6.3lf ) --> %6.3lf GeV/c\n",
           p_N2[0], p_N2[1], p_N2[2], p_N2.Mag());
    printf("\n\n");
  }
  
  // Check total momenta
  if (DEBUG) {
    printf("Check lepton momentum:\t%5.3lf = %5.3lf\n",pl,p_l.Mag());
    printf("Check nucleon momentum:\t%5.3lf = %5.3lf\n",pN,p_N.Mag());
    printf("Check kaon momentum:\t%5.3lf = %5.3lf\n",pK,p_K.Mag());
    printf("\n\n");
  }
  
  if (!RANDOM) {
    printf("COMPARE THIS TO EVENT #%d FROM CHRIS! (email 30.10.2014)\n\n",evt+1);
    return;
  }
  
  // Check that [numu -> mu+Q] conserves energy+momentum
  // -----------------------------------------------------
  if (VERBOSE) {
    double dX1 = p_l[0]+p_Q[0];
    double dY1 = p_l[1]+p_Q[1];
    double dZ1 = p_l[2]+p_Q[2] - Enu;
    double dE1 = (Enu-El + El) - (Enu);
    
    printf("CONSERVATION OF REACTION (nu -> l + q)\n");
    printf("---------------------------------------------------------\n");
    printf(" Particle |   px       py        pz       E        m \n");
    printf("---------------------------------------------------------\n");
    printf(" neutrino | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",0.,0.,Enu,Enu,0.);
    printf("---------------------------------------------------------\n");
    printf(" lepton   | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",p_l[0],p_l[1],p_l[2],El,ml);
    printf(" Q^2      | %6.3lf   %6.3lf   %6.3lf   %6.3lf     ---  \n",p_Q[0],p_Q[1],p_Q[2],Enu-El);
    printf("---------------------------------------------------------\n");
    printf(" FIN-INIT | %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",dX1,dY1,dZ1,dE1);
    printf("---------------------------------------------------------\n");
    printf("\n\n");
  }
  
  // Check that [Q+N -> K+N] conserves energy+momentum
  // -------------------------------------------------
  if (VERBOSE) {
    double dX2 = p_Kq[0]+p_Nq[0] - pQ*e3[0];
    double dY2 = p_Kq[1]+p_Nq[1] - pQ*e3[1];
    double dZ2 = p_Kq[2]+p_Nq[2] - pQ*e3[2];
    double dE2 = (EN+EK) - (Enu-El+mN);
    
    printf("CONSERVATION OF REACTION (q + N -> K + N) - ROTATED FRAME\n");
    printf("---------------------------------------------------------\n");
    printf(" Particle |   px       py        pz       E        m \n");
    printf("---------------------------------------------------------\n");
    printf(" Q^2      | %6.3lf   %6.3lf   %6.3lf   %6.3lf     ---  \n",pQ*e3[0],pQ*e3[1],pQ*e3[2],Enu-El);
    printf(" nucleon  | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",0.,0.,0.,mN,mN);
    printf("---------------------------------------------------------\n");
    printf(" kaon     | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",p_Kq[0],p_Kq[1],p_Kq[2],EK,mK);
    printf(" nucleon  | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",p_Nq[0],p_Nq[1],p_Nq[2],EN,mN);
    printf("---------------------------------------------------------\n");
    printf(" FIN-INIT | %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",dX2,dY2,dZ2,dE2);
    printf("---------------------------------------------------------\n");
    printf("\n\n");
  }
  
  // Check that [Q+N -> K+N] conserves energy+momentum
  // -------------------------------------------------
  if (VERBOSE) {
    double dX3 = p_K[0]+p_N[0] - p_Q[0];
    double dY3 = p_K[1]+p_N[1] - p_Q[1];
    double dZ3 = p_K[2]+p_N[2] - p_Q[2];
    double dE3 = (EN+EK) - (Enu-El+mN);
    
    printf("CONSERVATION OF REACTION (q + N -> K + N) - NORMAL FRAME\n");
    printf("---------------------------------------------------------\n");
    printf(" Particle |   px       py        pz       E        m \n");
    printf("---------------------------------------------------------\n");
    printf(" Q^2      | %6.3lf   %6.3lf   %6.3lf   %6.3lf     ---  \n",p_Q[0],p_Q[1],p_Q[2],Enu-El);
    printf(" nucleon  | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",0.,0.,0.,mN,mN);
    printf("---------------------------------------------------------\n");
    printf(" kaon     | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",p_K[0],p_K[1],p_K[2],EK,mK);
    printf(" nucleon  | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",p_N[0],p_N[1],p_N[2],EN,mN);
    printf("---------------------------------------------------------\n");
    printf(" FIN-INIT | %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",dX3,dY3,dZ3,dE3);
    printf("---------------------------------------------------------\n");
    printf("\n\n");
  }
  
  // Check that entire reaction conserves energy+momentum
  // ----------------------------------------------------
  double dX = p_l[0]+p_N[0]+p_K[0];
  double dY = p_l[1]+p_N[1]+p_K[1];
  double dZ = p_l[2]+p_N[2]+p_K[2] - Enu;
  double dE = (El+EN+EK) - (Enu+mN);
  
  printf("KINEMATICS OF PARTICLES IN REACTION (%s)\n",strReac.c_str());
  printf("---------------------------------------------------------\n");
  printf(" Particle |   px       py        pz       E        m \n");
  printf("---------------------------------------------------------\n");
  printf(" neutrino | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",0.,0.,Enu,Enu,0.);
  printf(" nucleon  | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",0.,0.,0.,mN,mN);
  printf("---------------------------------------------------------\n");
  printf(" lepton   | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",p_l[0],p_l[1],p_l[2],El,ml);
  printf(" nucleon  | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",p_N[0],p_N[1],p_N[2],EN,mN);
  printf(" kaon     | %6.3lf   %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",p_K[0],p_K[1],p_K[2],EK,mK);
  printf("---------------------------------------------------------\n");
  printf(" FIN-INIT | %6.3lf   %6.3lf   %6.3lf   %6.3lf \n",dX,dY,dZ,dE);
  printf("---------------------------------------------------------\n");
  printf("\n\n");
  
  const double EPS = 1e-6;
  if (abs(dX)<EPS && abs(dY)<EPS && abs(dZ)<EPS && abs(dE)<EPS)
    printf("INFO: ALL OK - energy & momentum are conserved.\n\n");
  else
    printf("WARNING: something is wrong, check E/p conservation!\n\n");
  return;
}

