
The first spline files for pion-nucleus total cross sections were made at Univ. of Pittsburgh in 2018.  
Those results had issues wihch made certain cases (e.g. pi charge exchange for Ar40) give odd results 
that aren't consistent with the undelying data.  In addition, only CEM03 was available for the first
effort and it has proven unable to fit data.  Finally, the code was still using pion-nucleus elastic
cross sections and this has been discontinued.  The new effort (2025) was undertaken by
Mohamed Ismail (mis90@pitt.edu) and Steve Dytman (dytman@pitt.edu) at Univ., of Pittsburgh
to fix these and other issues.  Dytman should be contacted for all questions.

The normalizing cross sections are now total reaction cross sections instead of the total cross
sections used for the 2018 versions.  Ashery (A>4) and Lehmann (A<=4) are used for 50<Tpi<400 MeV for 
targets where data exists.  INCL calculations are used at Tpi<=50 MeV to give good threshold effects.  
hN2018 calculations are used for Tpi>400 MeV where there is presently no data and calculations are mostly
reliable (minimal nuclear effects).

Splines are built at A=3, 7, 12, 27, 56, 93, and 209 for charge exchange (cex), inelastic (inel),
absorption, and pion production (pipro).  Files labeled e,g. c_pip_combined_abs.txt have the 
amalgam of the 3 sources of input data.  Files labeled e.g. pip12_abs.txt have the
results of smoothing via ROOT.  These files are all read into GENIE via INukeHadroData2025.
If probe is pi0, cross sections are modified inside HAIntranuke2025 according to results from hN.

Only hA pionA splines have been modified.  Similar updates of pA and nA splines are anticipated
in the future. 

references:
Ashery, et al., Phys Rev C 23, 2173 (1981)
Lehmann, et al., Phys Rev C 55, 2391 (1997)
SAID: downloaded from web site httphttps://gwdac.phys.gwu.edu/ 2009-2018.
