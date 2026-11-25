#include "Physics/NuclearState/NoFGM.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Conventions/Controls.h"

#include <TVector3.h>

using namespace genie;

//____________________________________________________________________________
NoFGM::NoFGM() :
NuclearModelI("genie::NoFGM")
{
    return;
}
//____________________________________________________________________________
NoFGM::NoFGM(string config) :
NuclearModelI("genie::NoFGM", config)
{
return;
}
//____________________________________________________________________________
NoFGM::~NoFGM()
{
    return;
}
//____________________________________________________________________________
bool   NoFGM::GenerateNucleon (const Target & t, double) const {
    return GenerateNucleon(t);
};
//____________________________________________________________________________
bool NoFGM::GenerateNucleon (const Target & t) const { 
    fCurrRemovalEnergy = utils::nuclear::BindEnergyPerNucleon(t);
    fCurrMomentum = TVector3(0,0,0);
    return true;
};
double NoFGM::Prob (double p, double , const Target & ) const {
    return p > controls::kASmallNum ? 0 : 1 ;
};
//____________________________________________________________________________
double NoFGM::Prob (double p, double w, const Target & t, double) const {
    return Prob(p, w, t);
};
//____________________________________________________________________________
double NoFGM::LocalFermiMomentum( const Target & , int , double  ) const {
    return 0;
};
//____________________________________________________________________________
void NoFGM::Configure (const Registry & /*config*/) {
    return;
};
//____________________________________________________________________________
void NoFGM::Configure (string /*param_set*/) {
    return;
};
//____________________________________________________________________________
void   NoFGM::LoadConfig (void) {
    return;
};
