{
//________________________________________________________________________________
/*
  makeCovMatrix.C
  
  A simple example ROOT script which turns a hard-coded matrix in the script (tDat)
  into a small ROOT file with just the matrix. The covariance matrix ROOT file
  is used as an input to grwghtcov, the covariance reweighting routine.

  The error matrix provided is for reweighting the four z-expansion parameters
  for the axial form factor. The numerical values are from e-Print arXiv: 1603.03048 [hep-ph],
  which is calculated from combining Eqs. (32) and (33) to form an error matrix.

 */
#include <TMatrixD.h>

  int nP = 4; // matrix size
  TFile *fMat = new TFile("tmat.out.root","recreate");

  // covariance matrix
  double tDat[nP][nP] =
    {{0.0154582, 0.0451836, -0.215641, 0.20647},
     {0.0451836, 1.08091, -2.38702, 1.0386},
     {-0.215641, -2.38702, 6.53568, -4.76577},
     {0.20647, 1.0386, -4.76577, 7.39832} };

  // flatten matrix into 1-D array
  double tFlt[nP* nP] = {0.};
  for (int i=0;i<nP;i++) {
  for (int j=0;j<nP;j++) {
    if (tFlt[j*nP+i] != tDat[i][j] && i > j)
    {
      std::cout<< "Matrix needs to be symmetric!" << std::endl;
      std::cout<< i<<","<<j<<": "<<tFlt[j*nP+i] <<","<< tDat[i][j]<<std::endl;
      exit(1);
    }
    tFlt[i*nP+j] = tDat[i][j];
  }}
  
  // make into TMatrixD object and write to file
  TMatrixD tMat(nP,nP);
  tMat.SetMatrixArray(tFlt);
  tMat->Write("tMat",TObject::kOverwrite);

  fMat->Close();
  exit(0);
}
