//____________________________________________________________________________
/*!

\class   genie::MCModel

\brief   A collection of cross section models

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#include <sstream>

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "ReWeight/MCModel.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;

//___________________________________________________________________________
MCModel::MCModel()
{
  fName = "unamed mc model";
  this->Initialize();
}
//___________________________________________________________________________
MCModel::MCModel(string name)
{
  fName = name;
  this->Initialize();
}
//___________________________________________________________________________
MCModel::MCModel(const MCModel & model)
{
  this->Copy(model);
}
//___________________________________________________________________________
MCModel::~MCModel()
{
  fXSecModelList.clear();
}
//___________________________________________________________________________
void MCModel::Copy(const MCModel & model)
{
  this->Initialize();

  map<string, const XSecAlgorithmI *>::const_iterator iter;
  for(iter = model.fXSecModelList.begin();
                           iter != model.fXSecModelList.end(); ++iter) {
     string key = iter->first;
     const XSecAlgorithmI * alg = iter->second;

     fXSecModelList.insert(
             map<string, const XSecAlgorithmI *>::value_type(key, alg));
  }
  fName = model.fName;
}
//___________________________________________________________________________
void MCModel::UseXSecAlg(const ProcessInfo & proc, const AlgId & algid)
{
  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * alg =
              dynamic_cast<const XSecAlgorithmI *>
                         (algf->GetAlgorithm(algid.Name(), algid.Config()));

  string key = this->BuildKey(proc);

  fXSecModelList.insert(
                 map<string, const XSecAlgorithmI *>::value_type(key, alg));
}
//___________________________________________________________________________
void MCModel::UseXSecAlg(
    const ProcessInfo & proc, const InitialState & init, const AlgId & algid)
{
  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * alg =
              dynamic_cast<const XSecAlgorithmI *>
                         (algf->GetAlgorithm(algid.Name(), algid.Config()));

  string key = this->BuildKey(proc, init);

  fXSecModelList.insert(
                 map<string, const XSecAlgorithmI *>::value_type(key, alg));
}
//___________________________________________________________________________
const XSecAlgorithmI * MCModel::XSecAlg (
                                       const Interaction * interaction) const
{
  if(!interaction) {
      LOG("ReWeight", pWARN) << "Null interaction!!";
      return 0;
  }

  LOG("ReWeight", pDEBUG)
     << "Finding cross section algorithm for: \n" << interaction->AsString();

  const ProcessInfo &  proc = interaction->ProcInfo();
  const InitialState & init = interaction->InitState();

  const XSecAlgorithmI * alg = 0;

  string key = this->BuildKey(proc, init);
  if(fXSecModelList.count(key) == 1)
  {
    map<string, const XSecAlgorithmI *>::const_iterator iter;
    iter = fXSecModelList.find(key);
    alg  = iter->second;

    LOG("ReWeight", pDEBUG)
           << "Key = " << key << " -> AlgId = " << alg->Id();
    return alg;
  }

  key = this->BuildKey(proc);
  if(fXSecModelList.count(key) == 1)
  {
    map<string, const XSecAlgorithmI *>::const_iterator iter;
    iter = fXSecModelList.find(key);
    alg  = iter->second;

    LOG("ReWeight", pDEBUG)
           << "Key = " << key << " -> AlgId = " << alg->Id();
    return alg;
  }

  LOG("ReWeight", pWARN)
       << "No cross section model for the input interaction";
  return 0;
}
//___________________________________________________________________________
void MCModel::Initialize(void)
{
  fXSecModelList.clear();
}
//___________________________________________________________________________
string MCModel::BuildKey(const ProcessInfo & proc) const
{
  ostringstream key;
  key << "PROC:" << proc.AsString();
  return key.str();
}
//___________________________________________________________________________
string MCModel::BuildKey(
                   const ProcessInfo & proc, const InitialState & init) const
{
  ostringstream key;
  key << "PROC:" << proc.AsString() << ";INIT:"  <<  init.AsString();
  return key.str();
}
//___________________________________________________________________________

