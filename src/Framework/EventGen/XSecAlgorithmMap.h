//____________________________________________________________________________
/*!

\class   genie::XSecAlgorithmMap

\brief   An Interaction -> XSecAlgorithmI associative container. The container
         is being built for the loaded EventGeneratorList and for the input
         InitialState object.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created January 23, 2006

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XSEC_ALGORITHM_MAP_H_
#define _XSEC_ALGORITHM_MAP_H_

#include <map>
#include <string>
#include <ostream>

using std::map;
using std::string;
using std::ostream;

namespace genie {

class XSecAlgorithmMap;
class XSecAlgorithmI;
class Interaction;
class InteractionList;
class InitialState;
class EventGeneratorList;

ostream & operator << (ostream & stream, const XSecAlgorithmMap & xsmap);

class XSecAlgorithmMap : public map<string, const XSecAlgorithmI *> {

public :

  XSecAlgorithmMap();
  XSecAlgorithmMap(const XSecAlgorithmMap & xsmap);
  ~XSecAlgorithmMap();

  void UseGeneratorList (const EventGeneratorList * list);
  void BuildMap         (const InitialState & init_state);

  const XSecAlgorithmI *  FindXSecAlgorithm  (const Interaction * in) const;
  const InteractionList & GetInteractionList (void) const;

  void Reset (void);
  void Copy  (const XSecAlgorithmMap & xsmap);
  void Print (ostream & stream) const;

  XSecAlgorithmMap & operator =  (const XSecAlgorithmMap & xsmap);
  friend ostream &   operator << (ostream & stream, const XSecAlgorithmMap & xsmap);

private:

  void Init    (void);
  void CleanUp (void);

  const EventGeneratorList * fEventGeneratorList;

  InitialState *    fInitState;
  InteractionList * fInteractionList;
};

}      // genie namespace

#endif // _XSEC_ALGORITHM_MAP_H_
