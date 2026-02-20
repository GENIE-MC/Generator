#include <string.h>
// declaration of PYTHIA6 common clocks
#ifndef WIN32
# define pyjets pyjets_
# define pydat1 pydat1_
# define pydat2 pydat2_
# define pydat3 pydat3_
# define pydat4 pydat4_
# define pydatr pydatr_
# define pysubs pysubs_
# define pypars pypars_
# define pyint1 pyint1_
# define pyint2 pyint2_
# define pyint3 pyint3_
# define pyint4 pyint4_
# define pyint5 pyint5_
# define pyint6 pyint6_
# define pyint7 pyint7_
# define pyint8 pyint8_
# define pyint9 pyint9_
# define pyuppr pyuppr_
# define pymssm pymssm_
# define pyssmt pyssmt_
# define pyints pyints_
# define pybins pybins_
#else
# define pyjets PYJETS
# define pydat1 PYDAT1
# define pydat2 PYDAT2
# define pydat3 PYDAT3
# define pydat4 PYDAT4
# define pydatr PYDATR
# define pysubs PYSUBS
# define pypars PYPARS
# define pyint1 PYINT1
# define pyint2 PYINT2
# define pyint3 PYINT3
# define pyint4 PYINT4
# define pyint5 PYINT5
# define pyint6 PYINT6
# define pyint7 PYINT7
# define pyint8 PYINT8
# define pyint9 PYINT9
# define pyuppr PYUPPR
# define pymssm PYMSSM
# define pyssmt PYSSMT
# define pyints PYINTS
# define pybins PYBINS
#endif

int pyjets[2+5*4000+2*2*5*4000];
int pydat1[200+2*200+200+2*200];
int pydat2[4*500+2*4*500+2*2000+2*4*4];
int pydat3[3*500+2*8000+2*8000+5*8000];  /* KNDCAY=8000 */
char pydat4[2*500*16];
int pydatr[6+2*100];
int pysubs[2+500+81*2+2*200];
int pypars[200+2*200+200+2*200];
int pyint1[400+2*400];
int pyint2[500+2*500+2*20*500+2*4*40];
int pyint3[2*81*2+3*1000+2*1000];
int pyint4[500+2*5*500];
int pyint5[1+3*501+2*3*501];
char pyint6[501*28];
int pyint7[2*6*7*7];
int pyint8[2*5*13];
int pyint9[2*4*13];
int pyuppr[1+7*20+1+2*10+2*5*20+2*11]; /* PYUPPR DOES NOT EXIST IN PYTHIA6 AT ALL!!! */
int pymssm[100+2*100];
int pyssmt[2*4*4+2*2*2+2*2*2+2*4+2*2+2*4*16+2*4*4+2*2*2+2*2*2];
int pyints[2*20];
int pybins[4+1000+2*20000];

void *pythia6_common_address(const char* name) {
   if      (!strcmp(name,"PYJETS")) return pyjets;
   else if (!strcmp(name,"PYDAT1")) return pydat1;
   else if (!strcmp(name,"PYDAT2")) return pydat2;
   else if (!strcmp(name,"PYDAT3")) return pydat3;
   else if (!strcmp(name,"PYDAT4")) return pydat4;
   else if (!strcmp(name,"PYDATR")) return pydatr;
   else if (!strcmp(name,"PYSUBS")) return pysubs;
   else if (!strcmp(name,"PYPARS")) return pypars;
   else if (!strcmp(name,"PYINT1")) return pyint1;
   else if (!strcmp(name,"PYINT2")) return pyint2;
   else if (!strcmp(name,"PYINT3")) return pyint3;
   else if (!strcmp(name,"PYINT4")) return pyint4;
   else if (!strcmp(name,"PYINT5")) return pyint5;
   else if (!strcmp(name,"PYINT6")) return pyint6;
   else if (!strcmp(name,"PYINT7")) return pyint7;
   else if (!strcmp(name,"PYINT8")) return pyint8;
   else if (!strcmp(name,"PYINT9")) return pyint9;
   else if (!strcmp(name,"PYUPPR")) return pyuppr;
   else if (!strcmp(name,"PYMSSM")) return pymssm;
   else if (!strcmp(name,"PYSSMT")) return pyssmt;
   else if (!strcmp(name,"PYINTS")) return pyints;
   else if (!strcmp(name,"PYBINS")) return pybins;
   return 0;
}
   
