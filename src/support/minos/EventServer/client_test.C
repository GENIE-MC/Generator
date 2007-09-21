//
// Test GENIE event server client
//
// C.Andreopoulos
//

gROOT->Reset();

int port = 9090;

TSocket * sock = 0;

bool handshake      (void);
void configure      (string particle_list);
void request_xsec   (string opt);
void request_event  (string opt);
void shutdown       (void);

//..........................................................................
void client_test()
{
  // handhake with genie event server
  //
  bool is_alive = handshake();
  assert(is_alive);

  // configure genie event server
  //
  configure("load-splines neutrino-list=14,-14,12,-12 target-list=1000260560");

  // request total cross section data (sum for all enabled channels) for
  // some initial state
  //
  request_xsec("14 1000260560");

  // request the generation of some events
  //
  request_event("14011011  1  14  1.849474  0.740705  53.184013 14 1000260560  0.012129 -0.941614 16.443062");
  request_event("14011011  2  14  2.187474  0.362736  22.218212 14 1000260560  0.120932  0.239200  3.212121");
  request_event("14011011  3  14  1.094340  0.127128  18.210291 14 1000260560  0.239001 -0.129101  8.029121");

  // shutdown the genie event server
  //
  shutdown();
}
//..........................................................................
bool handshake(void)
{
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
// Sytnax:
//   mesg sent:
//      CONFIG: load-splines neutrino-list=comma_separated_list_of_pdg_codes target_list=comma_separated_list_of_pdg_codes
//
   cout << "Configuring the GENIE event server" << endl;

   string cmd = "CONFIG: " + particle_list;

   sock->Send(cmd.c_str());
   cout << "Sent: " << cmd << endl;

   char mesg[2048];
   sock->Recv(mesg,2048);
   cout << "Received: " << mesg << endl;

   bool is_configured = (strcmp(mesg,"CONFIG COMPLETED") == 0);
   assert(is_configured);
}
//..........................................................................
void request_xsec(string opt)
{
// Syntax:
//   mesg sent:
//      XSEC: ipdgnu ipdgtgt
//   mesg recv (E in GeV, Sig in 1E-38 cm^2):
//      XSECSPL: npoints
//      knot(0) E(0) Sig(0) 
//      knot(1) E(1) Sig(1)
//      knot(2) E(2) Sig(2)
//      ....
//      

   string cmd = "XSEC: " + opt;

   sock->Send(cmd.c_str());
   cout << "Sent: " << cmd << endl;

   while(1) {
      char mesg[2048];
      sock->Recv(mesg,2048);
      cout << "Received: " << mesg << endl;

      bool exit_loop = (strcmp(mesg,"FAILED")==0) || 
                       (strcmp(mesg,"XSEC SENT")==0);

      if(exit_loop) break;
  }
}
//..........................................................................
void request_event(string opt)
{
// Syntax:
//   mesg sent:
//       EVTVTX: irun ievt ipdgnunosoc vtx_x vtx_y vtx_z ipdgnu ipdgtgt px_nu py_nu pz_nu
//   mesg recv:
//       EVTREC: irun ievt ipdgnunosoc vtx_x vtx_y vtx_z
//       ipdgnu ipdgtgt px_nu py_nu pz_nu
//       int_type iaction nucleon struck-quark xbj_int ybj_int w2_int q2_int
//       total-xsection diff-xsection ihadmod
//       STDHEP: nlines
//       indx ist ipdg jmo1 jmo2 jda1 jda2 px py pz E mass vx vy vz time
//        ...
//       indx ist ipdg jmo1 jmo2 jda1 jda2 px py pz E mass vx vy vz time
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
   sock->Send("SHUTDOWN");
   cout << "Sent: SHUTDOWN" << endl;

   char mesg[2048];
   sock->Recv(mesg,2048);
   cout << "Received: " << mesg << endl;

   sock->Close();
}
//..........................................................................

