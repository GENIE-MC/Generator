//____________________________________________________________________________
/*!

\program testRegistry

\brief   test program used for testing / debugging GENIE's Registry

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created May 4, 2004

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TH1F.h>

#include "Registry/Registry.h"
#include "Messenger/Messenger.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
 // manual override of Registry mesg level
 //Messenger * msg = Messenger::Instance();
 //msg->SetPriorityLevel("Registry", pDEBUG);

 // some objects stored at the test configuration registries
 RgAlg some_algorithm("alg-name","alg-config");

 TH1F * h1 = new TH1F("h1","a TH1F to be passed as config param", 100,-5,5);
 h1->Fill(3,100);
 h1->Fill(2,50);



 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Basic Tests:
 // Build a Registry, unlock it, add some vars, lock it, print it.
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 LOG("Main", pINFO) << "***** Basic Registry tests *****";

 LOG("Main", pINFO) << "Building, Unlocking, Setting, Locking, Printing...";

 Registry registry("example-registry");

 registry.UnLock();

 registry.Set("var-bool-1",        true);
 registry.Set("var-double-1",     2.289);
 registry.Set("var-double-2",     4.190);
 registry.Set("var-int-1",            0);
 registry.Set("var-int-2",           11);
 registry.Set(string("var-th1f-1"),  h1);
 registry.Set("myalg",   some_algorithm);

 registry.Lock();

 LOG("Main", pINFO) << "registry:\n" << registry;

 // try to override var-int-2 and add a new var-int-3 to the locked registry

 LOG("Main", pINFO)
              << "Trying to set variables in a locked registry - should fail";

 registry.Set("var-int-2",    12);
 registry.Set("var-int-3",    89);

 LOG("Main", pINFO) << "registry:\n" << registry;

 // do the same and add a couple of string vars, but now unlock it first

 LOG("Main", pINFO)
   << "Unlock the registry first and then set the variables - should succeed";

 registry.UnLock();

 registry.Set("var-int-2",    12);
 registry.Set("var-int-3",    89);
 registry.Set("var-string-1", string("this is a string"));
 registry.Set("var-string-2", string("and this is another string"));

 registry.Lock();

 LOG("Main", pINFO) << "registry\n" << registry;

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Testing copy constructor
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 LOG("Main", pINFO)
    << "***** Testing copy constructor: Registry registry2(registry) *****";

 Registry registry2(registry);

 LOG("Main", pINFO) << "Registry clone:    \n" << registry2;
 LOG("Main", pINFO) << "Original registry: \n" << registry;

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Testing operator () overloading
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 LOG("Main", pINFO) << "***** Testing operator () *****";

 registry2.UnLock();

 registry2("added an integer with ()", 7    );
 registry2("added a boolean with ()",  true );
 registry2("added a double with ()",   3.14 );
 registry2("added a string with ()",   "ok" );

 registry2.Lock();

 LOG("Main", pINFO) << "registry:\n" << registry2;

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Testing data retrieval
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 LOG("Main", pINFO) << "***** Testing data retrieval *****";

 bool   retrieved_bool   = false;
 int    retrieved_int    = 0;
 double retrieved_double = 3.14;
 string retrieved_string = "random init";

 registry.Get("var-bool-1",   retrieved_bool  );
 registry.Get("var-int-1",    retrieved_int   );
 registry.Get("var-double-1", retrieved_double);
 registry.Get("var-string-1", retrieved_string);

 LOG("Main", pINFO) << "retrieved-bool   = " << retrieved_bool;
 LOG("Main", pINFO) << "retrieved-int    = " << retrieved_int;
 LOG("Main", pINFO) << "retrieved-double = " << retrieved_double;
 LOG("Main", pINFO) << "retrieved-string = " << retrieved_string;

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Test Copy()
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 LOG("Main", pINFO) << "***** Testing Copy(const Registry & reg) *****";

 Registry registry3;        // create an empty registry
 registry3.Copy(registry2); // copy registry2 contents to registry3

 LOG("Main", pINFO) << "Registry clone:    \n" << registry3;
 LOG("Main", pINFO) << "Original registry: \n" << registry2;

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Test Individual Item Locking
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 LOG("Main", pINFO) << "***** Testing individual item locking *****";

 Registry registry4("example-registry-with-locked-items");

 registry4.UnLock();

 registry4("an int variable",    19     );
 registry4("a bool variable",    true   );
 registry4("a double variable",  2.71   );
 registry4("a string variable", "hello" );

 LOG("Main", pINFO)
    << "Initial registry: \n" << registry4;

 // lock two of the variables

 registry4.LockItem("an int variable");
 registry4.LockItem("a double variable");

 LOG("Main", pINFO)
    << "Registry with locked keys: \n" << registry4;

 LOG("Main", pINFO)
    << "Attempting to change locked items in unlocked registry";

 // try to change the locked variables - this should fail even though
 // the registry itself is unclocked

 registry4.Set("an int variable",   25);
 registry4.Set("a double variable", 129320.21);

 LOG("Main", pINFO)
    << "Should have failed to change locked entries: \n" << registry4;

 // inhibit individual item locking

 LOG("Main", pINFO) << "Inhibit item locking";
 registry4.InhibitItemLocks();

 LOG("Main", pINFO) 
   << "registry with item locks inhibited: \n" << registry4;

 // re-try to change the locked variables
 LOG("Main", pINFO)
    << "Retrying to change locked items in unlocked registry with inhibited item locking";

 registry4.Set("an int variable",   25);
 registry4.Set("a double variable", 9.21);

 LOG("Main", pINFO) << "registry: \n" << registry4;

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Test Assignment operator =
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 LOG("Main", pINFO) << "***** Testing operator = *****";

 Registry & registry5 = registry4;

 LOG("Main", pINFO) << "Printing registry set with the assignment operator = ";
 LOG("Main", pINFO) << "registry: \n" << registry5;
 
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Test operator +=
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 LOG("Main", pINFO) << "***** Testing operator += *****";

 registry5 += registry;

 LOG("Main", pINFO) << "Printing registry after adding values with the += operator ";
 LOG("Main", pINFO) << "registry: \n" << registry5;
}

