//
// Test GENIE event server client
//
// C.Andreopoulos
//

gROOT->Reset();

TSocket * sock = 0;

bool handshake      (void);
void configure      (string particle_list);
void generate_event (string opt);
void shutdown       (void);

//..........................................................................
void client_test(int port = 9090)
{
  bool is_alive = handshake();
  assert(is_alive);

  configure("load-splines neutrino-list=14,-14,12,-12 target-list=1000260560");

  generate_event("14011011  1  14  1.849474  0.740705  53.184013 14 1000260560  0.012129 -0.941614 16.443062");
  generate_event("14011011  2  14  2.187474  0.362736  22.218212 14 1000260560  0.120932  0.239200  3.212121");
  generate_event("14011011  3  14  1.094340  0.127128  18.210291 14 1000260560  0.239001 -0.129101  8.029121");

  shutdown();
}
//..........................................................................
bool handshake(void)
{
// Check that the GENIE event server is alive
//
   sock = new TSocket("localhost", port);
   if(!sock) return false;

   cout << "Pinging the GENIE event server" << endl;

   sock->Send("RUB GENIE LAMP");
   cout << "Sent: RUB GENIE LAMP" << endl;

   char mesg[2048];
   sock->Recv(mesg,2048);
   cout << "Received: " << mesg << endl;

   bool is_alive = (strcmp(mesg,"YOU HAVE 3 WISHES!") == 0);
   return is_alive;
}
//..........................................................................
void configure(string particle_list)
{
// Configure the GENIE event server
// mesg sytnax:
//   CONFIG: load-splines neutrino-list=comma_separated_list_of_pdg_codes target_list=comma_separated_list_of_pdg_codes
//
   cout << "Configuring the GENIE event server" << endl;

   string cmd = "CONFIG: " + particle_list;

   sock->Send(cmd.c_str());
   cout << "Sent: " << cmd << endl;

   char mesg[2048];
   sock->Recv(mesg,2048);
   cout << "Received: " << mesg << endl;

   bool is_configured = (strcmp(mesg,"CONFIGURATION COMPLETED") == 0);
   assert(is_configured);
}
//..........................................................................
void generate_event(string opt)
{
// Generate an event
// mesg sytnax:
//   EVTVTX: irun ievt ipdgnunosoc vtx_x vtx_y vtx_z ipdgnu ipdgtgt px_nu py_nu pz_nu
//
   string cmd = "EVTVTX: " + opt;

   sock->Send(cmd.c_str());
   cout << "Sent: " << cmd << endl;

   while(1) {
      char mesg[2048];
      sock->Recv(mesg,2048);
      cout << "Received: " << mesg << endl;

      bool exit_loop = (strcmp(mesg,"FAILED")==0) || 
                       (strcmp(mesg,"EVENT GENERATED")==0);

      if(exit_loop) break;
  }
}
//..........................................................................
void shutdown(void)
{
// Shutdown server
//
   sock->Send("SHUTDOWN");
   cout << "Sent: SHUTDOWN" << endl;

   char mesg[2048];
   sock->Recv(mesg,2048);
   cout << "Received: " << mesg << endl;

   sock->Close();
}
//..........................................................................

