//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenConfig

\brief    Encapsulation of NeuGEN's Configuration Card

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#include "Facades/NeuGenConfig.h"

using std::endl;
using namespace genie::nuvld::facades;

ClassImp(NeuGenConfig)

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
  namespace facades {
    ostream & operator << (ostream & stream, const NeuGenConfig & conf)
    {
      conf.Print(stream);
      return stream;
    }
  } 
 }
}    
//____________________________________________________________________________
NeuGenConfig::NeuGenConfig() :
fName("default")
{
  this->SetBestParameters();
}
//____________________________________________________________________________
NeuGenConfig::NeuGenConfig(const char * name) :
fName(name)
{
  this->SetBestParameters();
}
//____________________________________________________________________________
NeuGenConfig::NeuGenConfig(const NeuGenConfig * config) 
{
  this->Copy(config);
}
//____________________________________________________________________________
NeuGenConfig::~NeuGenConfig()
{

}
//____________________________________________________________________________
void NeuGenConfig::SetBestParameters(void)
{
  fPdfGrp   =  kNGDefPdfGrp;
  fPdfSet   =  kNGDefPdfSet;
  fMaQel    =  kNGDefMaQel;
  fMaRes    =  kNGDefMaRes;
  fMaCoh    =  kNGDefMaCoh;
  fFa0Qel   =  kNGDefFa0Qel;
  fEtaQel   =  kNGDefEtaQel;
  fOmegaRes =  kNGDefOmegaRes;
  fZRes     =  kNGDefZRes;
  fR0Coh    =  kNGDefR0Coh;
  fREICoh   =  kNGDefREICoh;
  fNuDisFF  =  kNGDefNuDisFF;
    
  for(unsigned int istate = 0; istate < kNInSt; istate++)   {
                                          fKnoA[istate] = kNGDefKnoA[istate];
                                          fKnoB[istate] = kNGDefKnoB[istate];
                                          fKnoC[istate] = kNGDefKnoC[istate];
  }

  for(unsigned int imulti = 0; imulti < kNMlt; imulti++)
          for(unsigned int istate = 0; istate < kNInSt; istate++)
                      fDisRes[imulti][istate] = kNGDefDisRes[imulti][istate];
}
//____________________________________________________________________________
float NeuGenConfig::KnoA(NGInitState_t initial_state) const
{
 bool is_valid = this->IsValidInitState(initial_state);
 
 if( is_valid ) {

    int istate = this->InitState2IPos(initial_state);

    return fKnoA[istate];
    
 } else return -1;
}
//____________________________________________________________________________
float NeuGenConfig::KnoB(NGInitState_t initial_state) const
{
 bool is_valid = this->IsValidInitState(initial_state);
 
 if( is_valid ) {

    int istate = this->InitState2IPos(initial_state);

    return fKnoB[istate];
    
 } else return -1;
}
//____________________________________________________________________________
float NeuGenConfig::KnoC(NGInitState_t initial_state) const
{
 bool is_valid = this->IsValidInitState(initial_state);
 
 if( is_valid ) {

    int istate = this->InitState2IPos(initial_state);

    return fKnoC[istate];
    
 } else return -1;
}
//____________________________________________________________________________
float NeuGenConfig::DisRes(unsigned int multiplicity,
                                             NGInitState_t initial_state) const
{
 bool is_valid = this->IsValidMultiplicity(multiplicity) &&
                                   this->IsValidInitState(initial_state);
 if( is_valid ) {

    int imulti = this->Multiplicity2IPos (multiplicity);
    int istate = this->InitState2IPos    (initial_state);

    return fDisRes[imulti][istate];

 } else return -1;
}
//____________________________________________________________________________
void NeuGenConfig::SetKnoA(NGInitState_t initial_state, float kno)
{
 bool is_valid = this->IsValidInitState(initial_state);

 if( is_valid ) {
    int istate = this->InitState2IPos(initial_state);
    
    fKnoA[istate] = kno;
 }
}
//____________________________________________________________________________
void NeuGenConfig::SetKnoB(NGInitState_t initial_state, float kno)
{
 bool is_valid = this->IsValidInitState(initial_state);

 if( is_valid ) {
    int istate = this->InitState2IPos(initial_state);
    
    fKnoB[istate] = kno;
 }
}
//____________________________________________________________________________
void NeuGenConfig::SetKnoC(NGInitState_t initial_state, float kno)
{
 bool is_valid = this->IsValidInitState(initial_state);

 if( is_valid ) {
    int istate = this->InitState2IPos(initial_state);
    
    fKnoC[istate] = kno;
 }
}
//____________________________________________________________________________
void NeuGenConfig::SetDisRes(unsigned int multiplicity,
                                    NGInitState_t initial_state, float dis_res)
{
 bool is_valid = this->IsValidMultiplicity(multiplicity) &&
                                    this->IsValidInitState(initial_state);
 if( is_valid ) {

    int imulti = this->Multiplicity2IPos (multiplicity);
    int istate = this->InitState2IPos    (initial_state);

    fDisRes[imulti][istate] = dis_res;
 }
}
//____________________________________________________________________________
void NeuGenConfig::Copy(const NeuGenConfig * config) 
{
  fName     = config->fName;
  fPdfGrp   = config->fPdfGrp;
  fPdfSet   = config->fPdfSet;
  fMaQel    = config->fMaQel;
  fMaRes    = config->fMaRes;
  fMaCoh    = config->fMaCoh;    
  fFa0Qel   = config->fFa0Qel;   
  fEtaQel   = config->fEtaQel;   
  fOmegaRes = config->fOmegaRes; 
  fZRes     = config->fZRes;     
  fR0Coh    = config->fR0Coh;    
  fREICoh   = config->fREICoh;   
  fNuDisFF  = config->fNuDisFF;
 
  for(unsigned int istate = 0; istate < kNInSt; istate++)   {
                                      fKnoA[istate] = config->fKnoA[istate];
                                      fKnoB[istate] = config->fKnoB[istate];
                                      fKnoC[istate] = config->fKnoC[istate];
  }

  for(unsigned int imult = 0; imult < kNMlt; imult++)
             for(unsigned int istate = 0; istate < kNInSt; istate++)
                    fDisRes[imult][istate] = config->fDisRes[imult][istate];
}
//____________________________________________________________________________
void NeuGenConfig::Print(ostream & stream) const
{
  stream << "Axial Mass (Quasi Elastic):.............." << fMaQel    << endl;
  stream << "Axial Mass (Resonance):.................." << fMaRes    << endl;
  stream << "Axial mass (Coherent):..................." << fMaCoh    << endl;
  stream << "Fa(Q^2 = 0) (Quasi Elastic):............." << fFa0Qel   << endl;
  stream << "Elastic scattering parameter:............" << fEtaQel   << endl;
  stream << "R-S Model parameter Omega:..............." << fOmegaRes << endl;
  stream << "R-S Model parameter Z:..................." << fZRes     << endl;
  stream << "Nuclear size scale param:................" << fR0Coh    << endl;
  stream << "Re/Im for pion scattering amplitude:....." << fREICoh   << endl;
  stream << "Neutrino DIS Scale Factor:..............." << fNuDisFF  << endl;
  stream << "KNO Hadronization Param A / vp .........."
                               << fKnoA[this->InitState2IPos(e_vp)]  << endl;
  stream << "KNO Hadronization Param A / vn .........."
                               << fKnoA[this->InitState2IPos(e_vn)]  << endl;
  stream << "KNO Hadronization Param A / vbp ........."
                               << fKnoA[this->InitState2IPos(e_vbp)] << endl;
  stream << "KNO Hadronization Param A / vbn ........."
                               << fKnoA[this->InitState2IPos(e_vbn)] << endl;
  stream << "KNO Hadronization Param B / vp .........."
                               << fKnoB[this->InitState2IPos(e_vp)]  << endl;
  stream << "KNO Hadronization Param B / vn .........."
                               << fKnoB[this->InitState2IPos(e_vn)]  << endl;
  stream << "KNO Hadronization Param B / vbp ........."
                               << fKnoB[this->InitState2IPos(e_vbp)] << endl;
  stream << "KNO Hadronization Param B / vbn ........."
                               << fKnoB[this->InitState2IPos(e_vbn)] << endl;
  stream << "KNO Hadronization Param C / vp .........."
                               << fKnoC[this->InitState2IPos(e_vp)]  << endl;
  stream << "KNO Hadronization Param C / vn .........."
                               << fKnoC[this->InitState2IPos(e_vn)]  << endl;
  stream << "KNO Hadronization Param C / vbp ........."
                               << fKnoC[this->InitState2IPos(e_vbp)] << endl;
  stream << "KNO Hadronization Param C / vbn ........."
                               << fKnoC[this->InitState2IPos(e_vbn)] << endl;
  stream << "DIS/RES Tuning Param - 2 / vp ......" 
         << fDisRes[this->Multiplicity2IPos(2)][this->InitState2IPos(e_vp)]
         << endl;
  stream << "DIS/RES Tuning Param - 3 / vp ......" 
         << fDisRes[this->Multiplicity2IPos(3)][this->InitState2IPos(e_vp)]
         << endl;
  stream << "DIS/RES Tuning Param - 2 / vn ......" 
         << fDisRes[this->Multiplicity2IPos(2)][this->InitState2IPos(e_vn)]
         << endl;
  stream << "DIS/RES Tuning Param - 3 / vn ......" 
         << fDisRes[this->Multiplicity2IPos(3)][this->InitState2IPos(e_vn)]
         << endl;
  stream << "DIS/RES Tuning Param - 2 / vbp ....." 
         << fDisRes[this->Multiplicity2IPos(2)][this->InitState2IPos(e_vbp)]
         << endl;
  stream << "DIS/RES Tuning Param - 3 / vbp ....." 
         << fDisRes[this->Multiplicity2IPos(3)][this->InitState2IPos(e_vbp)]
         << endl;
  stream << "DIS/RES Tuning Param - 2 / vbn ....." 
         << fDisRes[this->Multiplicity2IPos(2)][this->InitState2IPos(e_vbn)]
         << endl;
  stream << "DIS/RES Tuning Param - 3 / vbn ....." 
         << fDisRes[this->Multiplicity2IPos(3)][this->InitState2IPos(e_vbn)]
         << endl;
}
//____________________________________________________________________________
int NeuGenConfig::Multiplicity2IPos(unsigned int multiplicity) const
{
  switch(multiplicity) {
    case (2): return  0; break;
    case (3): return  1; break;
    default : return -1; break;
  }  
}                                              
//____________________________________________________________________________
int NeuGenConfig::InitState2IPos(NGInitState_t state) const
{
  switch(state) {
    case (e_vp):   return  0; break;
    case (e_vn):   return  1; break;
    case (e_vbp):  return  2; break;
    case (e_vbn):  return  3; break;
    default :      return -1; break;
  }
}                                         
//____________________________________________________________________________
bool NeuGenConfig::IsValidMultiplicity(unsigned int multiplicity) const
{
  if(multiplicity == 2 || multiplicity == 3) return true;
  else return false;
}
//____________________________________________________________________________
bool NeuGenConfig::IsValidInitState(NGInitState_t state) const
{
  bool is_valid_state = false;
  
  switch(state) {
    case (e_vp):
    case (e_vn):
    case (e_vbp):
    case (e_vbn):   is_valid_state = true;   break;
    default:        is_valid_state = false;  break;
  }

  return is_valid_state;
}
//____________________________________________________________________________

