//_________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/XSection/NievesSimoVacasMECPXSec2016.h"
#include "Physics/Multinucleon/XSection/MECHadronTensor.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//_________________________________________________________________________
NievesSimoVacasMECPXSec2016::NievesSimoVacasMECPXSec2016() :
XSecAlgorithmI("genie::NievesSimoVacasMECPXSec2016")
{

}
//_________________________________________________________________________
NievesSimoVacasMECPXSec2016::NievesSimoVacasMECPXSec2016(string config) :
XSecAlgorithmI("genie::NievesSimoVacasMECPXSec2016", config)
{

}
//_________________________________________________________________________
NievesSimoVacasMECPXSec2016::~NievesSimoVacasMECPXSec2016()
{

}
//_________________________________________________________________________
double NievesSimoVacasMECPXSec2016::XSec(
  const Interaction * interaction, KinePhaseSpace_t kps) const
{
// The basic quantity calculated is d2sigma/(dTmu dcos_mu)

    // Get hadron tensor
    MECHadronTensor * hadtensor = MECHadronTensor::Instance();

    int targetpdg = interaction->InitState().Tgt().Pdg();
    int Arequest = pdg::IonPdgCodeToA(targetpdg);
    int Zrequest = pdg::IonPdgCodeToZ(targetpdg);

    // To generate cross-sections for nuclei other than those with hadron 
    // tensors we need to pull both the full cross-section and 
    // the pn initial state fraction.
    // Non-isoscalar nuclei are beyond the original published Valencia model
    // and scale with A according to the number of pp, pn, or nn pairs
    // the probe is expected to find.  
    // There is some by-hand optimization here, skipping the delta part when
    // only the total cross-section is requested.
    // Possible future models without a Delta had tensor would also use that
    // flag to call this without computing the Delta part.

    int tensorpdg = targetpdg;

    if( ! hadtensor->KnownTensor(tensorpdg) ) {

        // Decide which hadron tensor to use for different ranges of A.

        if( Arequest == 4 && Zrequest == 2 ){
            tensorpdg = kPdgTgtC12;
            // This is for helium 4, but use carbon tensor
            // the use of nuclear density parameterization is suspicious
            // but some users (MINERvA) need something not nothing.
            // The pn will be exactly 1/3, but pp and nn will be ~1/4
            // Because the combinatorics are different.
            // Could do lithium beryllium boron which you don't need
        }
        else if (Arequest < 9){
            // refuse to do D, T, He3, Li, and some Be, B
            // actually it would work technically, maybe except D, T
            MAXLOG("NievesSimoVacasMEC", pWARN, 10)
                << "Asked to scale to deuterium through boron "
                << targetpdg << " nope, lets not do that.";
            return 0;
        }
        else if( Arequest >= 9 && Arequest < 15){
            tensorpdg = kPdgTgtC12;
            //}
            // could explicitly put in nitrogen for air
            //else if ( Arequest >= 14 && A < 15) { // AND CHANGE <=14 to <14.
            //  tensorpdg = kPdgTgtN14;
        } 
        else if( Arequest >= 15 && Arequest < 22){
            tensorpdg = kPdgTgtO16;
        } 
        else if( Arequest >= 22 && Arequest < 33){
            // of special interest, this gets Al27 and Si28
            tensorpdg = 1000140280;
        } 
        else if(Arequest >= 33 && Arequest < 50){
            // of special interest, this gets Ar40 and Ti48   
            tensorpdg = kPdgTgtCa40;
        } 
        else if( Arequest >= 50 && Arequest < 90){
            // pseudoFe56, also covers many other ferrometals and Ge
            tensorpdg = 1000280560;
        } 
        else if( Arequest >= 90 && Arequest < 160){
            // use Ba112 = PseudoCd.  Row5 of Periodic table useless. Ag, Xe?
            tensorpdg = 1000561120;
        } 
        else if( Arequest >= 160 ){
            // use Rf208 = pseudoPb
            tensorpdg = 1001042080;   
        } 
        else {
            MAXLOG("NievesSimoVacasMEC", pWARN, 10) 
                << "Asked to scale to a nucleus " 
                << targetpdg << " which we don't know yet.";
            return 0;
        }  
    }

    // Check that the input kinematical point is within the range
    // in which hadron tensors are known (for chosen target)
    double Ev    = interaction->InitState().ProbeE(kRfLab);
    double Tl    = interaction->Kine().GetKV(kKVTl);
    double costl = interaction->Kine().GetKV(kKVctl);
    double ml    = interaction->FSPrimLepton()->Mass();
    double Q0    = 0;
    double Q3    = 0;
    genie::utils::mec::Getq0q3FromTlCostl(Tl, costl, Ev, ml, Q0, Q3);
    const vector <genie::BLI2DNonUnifGrid *> &
        tensor_table = hadtensor->TensorTable(
                tensorpdg, MECHadronTensor::kMHTValenciaFullAll);
    double Q0min = tensor_table[0]->XMin();
    double Q0max = tensor_table[0]->XMax();
    double Q3min = tensor_table[0]->YMin();
    double Q3max = tensor_table[0]->YMax();
    if(Q0 < Q0min || Q0 > Q0max || Q3 < Q3min || Q3 > Q3max) {
        return 0.0;
    }

    // By default, compute will compute the full cross-section.
    // If a resonance is set, will will calculate the part of the cross-section
    // with an internal Delta line without a final state pion (usually called PPD
    // for pioness Delta decay).
    // If a {p,n} hit dinucleon was set will calculate the cross-section
    // for that component only (either full or PDD cross-section)
    bool delta = interaction->ExclTag().KnownResonance(); 
    bool pn    = (interaction->InitState().Tgt().HitNucPdg() == kPdgClusterNP);
    //LOG("NievesSimoVacasMEC", pDEBUG) << "delta: " << delta << ", pn: " << pn;

    double xsec_all = 0;
    double xsec_pn  = 0;
    if(delta) {
        xsec_all = genie::utils::mec::TensorContraction(interaction,
                tensorpdg, MECHadronTensor::kMHTValenciaDeltaAll );
        xsec_pn  = genie::utils::mec::TensorContraction(interaction,
                tensorpdg, MECHadronTensor::kMHTValenciaDeltapn  );
    } else {
        xsec_all = genie::utils::mec::TensorContraction(interaction,
                tensorpdg, MECHadronTensor::kMHTValenciaFullAll  );
        xsec_pn  = genie::utils::mec::TensorContraction(interaction,
                tensorpdg, MECHadronTensor::kMHTValenciaFullpn   );
    }

    // now again, run this only if targetpdg != tensorpdg.
    // would need to trap and treat He3, T, D special here.
    if( ! hadtensor->KnownTensor(targetpdg) ) {
        // if we need to scale, figure it out here.

        double PP = Zrequest;
        double NN = Arequest - PP;
        double P  = pdg::IonPdgCodeToZ(tensorpdg);
        double N  = pdg::IonPdgCodeToA(tensorpdg) - P;

        double scale_pn = TMath::Sqrt( (PP*NN)/(P*N) );
        double scale_pp = TMath::Sqrt( (PP * (PP - 1.)) / (P * (P - 1.)) );
        double scale_nn = TMath::Sqrt( (NN * (NN - 1.)) / (N * (N - 1.)) );

        LOG("NievesSimoVacasMEC", pDEBUG) 
            << "Scale pn pp nn for (" << targetpdg << ", " << tensorpdg << ")"
            << " : " << scale_pn << " " << scale_pp << " " << scale_nn;

        // this is an approximation in at least three senses.
        // we are scaling from an isoscalar nucleus using p and n counting
        // we are not using the right qvalue in the had tensor
        // we are not scaling the Delta faster than the non-Delta.
        // the guess is that these are good approximations.
        // a test we could document is to scale from O16 to N14 or C12 using this algorithm
        // and see how many percent deviation we see from the full calculation.

        double temp_all = xsec_all;
        double temp_pn  = xsec_pn * scale_pn;
        int nupdg = interaction->InitState().ProbePdg();
        if(nupdg > 0){
            temp_all = xsec_pn * scale_pn + (xsec_all - xsec_pn) * scale_nn;
        } else {
            temp_all = xsec_pn * scale_pn + (xsec_all - xsec_pn) * scale_pp;
        }
        xsec_all = temp_all;
        xsec_pn  = temp_pn;
    } 

    double xsec = (pn) ? xsec_pn : xsec_all;

    // Apply given scaling factor
    xsec *= fXSecScale;

    if(kps!=kPSTlctl) {
        LOG("NievesSimoVacasMEC", pWARN)
            << "Doesn't support transformation from "
            << KinePhaseSpace::AsString(kPSTlctl) << " to "
            << KinePhaseSpace::AsString(kps);
        xsec = 0;
    }

    return xsec;
}
//_________________________________________________________________________
double NievesSimoVacasMECPXSec2016::Integral(
        const Interaction * interaction) const 
{
    double xsec = fXSecIntegrator->Integrate(this,interaction);
    return xsec;
}
//_________________________________________________________________________
bool NievesSimoVacasMECPXSec2016::ValidProcess(
        const Interaction * interaction) const 
{
    if (interaction->TestBit(kISkipProcessChk)) return true;

    const ProcessInfo & proc_info = interaction->ProcInfo();
    if (!proc_info.IsMEC()) {
        return false;
    }
    return true;
}
//_________________________________________________________________________
void NievesSimoVacasMECPXSec2016::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void NievesSimoVacasMECPXSec2016::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//_________________________________________________________________________
void NievesSimoVacasMECPXSec2016::LoadConfig(void)
{
	// Cross section scaling factor
	GetParam( "MEC-CC-XSecScale", fXSecScale ) ;

	fXSecIntegrator =
        dynamic_cast<const XSecIntegratorI *> (
                this->SubAlg("NumericalIntegrationAlg"));
    assert(fXSecIntegrator);

}
//_________________________________________________________________________
