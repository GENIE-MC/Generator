//____________________________________________________________________________
/*!

\program testRegistry

\brief   test program used for testing / debugging GENIE's Registry

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 4, 2004
*/
//____________________________________________________________________________

#include "Config/Registry.h"
#include "Messenger/Messenger.h"

using namespace genie;

int main(int argc, char ** argv)
{
 //-- Build a Registry, unlock it, add some vars, lock it, print it.

 Registry registry("example-registry");

 registry.UnLock();

 registry.Set("var-bool-1",   true);
 registry.Set("var-double-1", 2.289);
 registry.Set("var-double-2", 4.190);
 registry.Set("var-int-1",    0);
 registry.Set("var-int-2",    11);

 registry.Lock();

 LOG("Main", pINFO) << ENDL << registry;

 //-- Try to override var-int-2 and add a new var-int-3 to the locked registry

 registry.Set("var-int-2",    12);
 registry.Set("var-int-3",    89);

 LOG("Main", pINFO) << ENDL << registry;

 //-- Do the same and add a couple of string vars, but now unlock it first

 registry.UnLock();

 registry.Set("var-int-2",    12);
 registry.Set("var-int-3",    89);
 registry.Set("var-string-1", "this is a string");
 registry.Set("var-string-2", "and this is another string");

 registry.Lock();

 LOG("Main", pINFO) << ENDL << registry;

 //-- Testing copy constructor

 Registry registry2(registry);

 LOG("Main", pINFO) << "Testing copy constructor: Registry registry2(registry)";
 
 LOG("Main", pINFO) << ENDL << registry2;

 //-- Testing () operator overloading

 LOG("Main", pINFO) << "Testing () operator overloading";

 registry2.UnLock();

 registry2("int .....added with ()", 7    );
 registry2("bool ....added with ()", true );
 registry2("double ..added with ()", 3.14 );
 registry2("string ..added with ()", "ok" );

 registry2.Lock();

 LOG("Main", pINFO) << ENDL << registry2;

 //-- Testing data retrieval

 LOG("Main", pINFO) << "Testing data retrieval";

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
}

