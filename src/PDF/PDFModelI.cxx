//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "PDF/PDFModelI.h"

using namespace genie;

//____________________________________________________________________________
PDFModelI::PDFModelI() :
Algorithm()
{

}
//____________________________________________________________________________
PDFModelI::PDFModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
PDFModelI::PDFModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
PDFModelI::~PDFModelI() 
{ 

}
//____________________________________________________________________________



