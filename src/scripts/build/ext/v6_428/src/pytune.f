C*********************************************************************
 
C...PYTUNE
C...Presets for a few specific underlying-event and min-bias tunes
C...Note some tunes require external pdfs to be linked (e.g. 105:QW),
C...others require particular versions of pythia (e.g. the SCI and GAL
C...models). See below for details.
      SUBROUTINE PYTUNE(MYTUNE)
C
C ITUNE    NAME (detailed descriptions below)
C     0 Default : No settings changed => defaults.
C
C ====== Old UE, Q2-ordered showers ====================================
C   100       A : Rick Field's CDF Tune A                     (Oct 2002)
C   101      AW : Rick Field's CDF Tune AW                    (Apr 2006)
C   102      BW : Rick Field's CDF Tune BW                    (Apr 2006)
C   103      DW : Rick Field's CDF Tune DW                    (Apr 2006)
C   104     DWT : As DW but with slower UE ECM-scaling        (Apr 2006)
C   105      QW : Rick Field's CDF Tune QW using CTEQ6.1M            (?)
C   106 ATLAS-DC2: Arthur Moraes' (old) ATLAS tune ("Rome")          (?)
C   107     ACR : Tune A modified with new CR model           (Mar 2007)
C   108      D6 : Rick Field's CDF Tune D6 using CTEQ6L1             (?)
C   109     D6T : Rick Field's CDF Tune D6T using CTEQ6L1            (?)
C ---- Professor Tunes : 110+ (= 100+ with Professor's tune to LEP) ----
C   110   A-Pro : Tune A, with LEP tune from Professor        (Oct 2008)
C   111  AW-Pro : Tune AW, -"-                                (Oct 2008)
C   112  BW-Pro : Tune BW, -"-                                (Oct 2008)
C   113  DW-Pro : Tune DW, -"-                                (Oct 2008)
C   114 DWT-Pro : Tune DWT, -"-                               (Oct 2008)
C   115  QW-Pro : Tune QW, -"-                                (Oct 2008)
C   116 ATLAS-DC2-Pro: ATLAS-DC2 / Rome, -"-                  (Oct 2008)
C   117 ACR-Pro : Tune ACR, -"-                               (Oct 2008)
C   118  D6-Pro : Tune D6, -"-                                (Oct 2008)
C   119 D6T-Pro : Tune D6T, -"-                               (Oct 2008)
C ---- Professor's Q2-ordered Perugia Tune : 129 -----------------------
C   129 Pro-Q2O : Professor Q2-ordered tune                   (Feb 2009)
C ---- LHC tune variations on Pro-Q2O 
C   136 Q12-F1  : Variation with wide fragmentation function (Mar 2012)
C   137 Q12-F2  : Variation with narrow fragmentation function (Mar 2012)
C
C ====== Intermediate and Hybrid Models ================================
C   200    IM 1 : Intermediate model: new UE, Q2-ord. showers, new CR
C   201     APT : Tune A w. pT-ordered FSR                    (Mar 2007)
C   211 APT-Pro : Tune APT, with LEP tune from Professor      (Oct 2008)
C   221 Perugia APT  : "Perugia" update of APT-Pro            (Feb 2009)
C   226 Perugia APT6 : "Perugia" update of APT-Pro w. CTEQ6L1 (Feb 2009)
C
C ====== New UE, interleaved pT-ordered showers, annealing CR ==========
C   300      S0 : Sandhoff-Skands Tune using the S0 CR model  (Apr 2006)
C   301      S1 : Sandhoff-Skands Tune using the S1 CR model  (Apr 2006)
C   302      S2 : Sandhoff-Skands Tune using the S2 CR model  (Apr 2006)
C   303     S0A : S0 with "Tune A" UE energy scaling          (Apr 2006)
C   304    NOCR : New UE "best try" without col. rec.         (Apr 2006)
C   305     Old : New UE, original (primitive) col. rec.      (Aug 2004)
C   306 ATLAS-CSC: Arthur Moraes' (new) ATLAS tune w. CTEQ6L1 (?)
C ---- Professor Tunes : 310+ (= 300+ with Professor's tune to LEP)
C   310   S0-Pro : S0 with updated LEP pars from Professor    (Oct 2008)
C   311   S1-Pro : S1 -"-                                     (Oct 2008)
C   312   S2-Pro : S2 -"-                                     (Oct 2008)
C   313  S0A-Pro : S0A -"-                                    (Oct 2008)
C   314 NOCR-Pro : NOCR -"-                                   (Oct 2008)
C   315  Old-Pro : Old -"-                                    (Oct 2008)
C   316  ATLAS MC08 : pT-ordered showers, CTEQ6L1             (2008)
C ---- Peter's Perugia Tunes : 320+ ------------------------------------
C   320 Perugia 0 : "Perugia" update of S0-Pro                (Feb 2009)
C   321 Perugia HARD : More ISR, More FSR, Less MPI, Less BR, Less HAD
C   322 Perugia SOFT : Less ISR, Less FSR, More MPI, More BR, More HAD
C   323 Perugia 3 : Alternative to Perugia 0, with different ISR/MPI
C                   balance & different scaling to LHC & RHIC (Feb 2009)
C   324 Perugia NOCR : "Perugia" update of NOCR-Pro           (Feb 2009)
C   325 Perugia * : "Perugia" Tune w. (external) MRSTLO* PDFs (Feb 2009)
C   326 Perugia 6 : "Perugia" Tune w. (external) CTEQ6L1 PDFs (Feb 2009)
C   327 Perugia 10: Alternative to Perugia 0, with more FSR   (May 2010)
C                   off ISR, more BR breakup, more strangeness
C   328 Perugia K : Alternative to Perugia 2010, with a       (May 2010)   
C                   K-factor applied to MPI cross sections
C ---- Professor's pT-ordered Perugia Tune : 329 -----------------------
C   329 Pro-pTO   : Professor pT-ordered tune w. S0 CR model  (Feb 2009)
C ---- Tunes introduced in 6.4.23:
C   330 ATLAS MC09 : pT-ordered showers, LO* PDFs             (2009)
C   331 ATLAS MC09c : pT-ordered showers, LO* PDFs, better CR (2009)
C   334 Perugia 10 NOCR : Perugia 2010 with no CR, less MPI   (Oct 2010)
C   335 Pro-pT*   : Professor Tune with LO*                   (Mar 2009)
C   336 Pro-pT6   : Professor Tune with CTEQ6LL               (Mar 2009)
C   339 Pro-pT**  : Professor Tune with LO**                  (Mar 2009)
C   340 AMBT1   : First ATLAS tune including 7 TeV data       (May 2010)
C   341 Z1      : First CMS tune including 7 TeV data         (Aug 2010)
C   342 Z1-LEP  : CMS tune Z1, with improved LEP parameters   (Oct 2010)
C   343 Z2      : Retune of Z1 by Field w CTEQ6L1 PDFs            (2010)
C   344 Z2-LEP  : Retune of Z1 by Skands w CTEQ6L1 PDFs       (Feb 2011)
C   345 AMBT2B-CT6L : 2nd ATLAS MB tune, vers 'B', w CTEQ6L1  (Jul 2011)
C   346 AUET2B-CT6L : UE tune accompanying AMBT2B             (Jul 2011)
C   347 AUET2B-CT66 : AUET2 with CTEQ 6.6 NLO PDFs            (Nov 2011)
C   348 AUET2B-CT10 : AUET2 with CTEQ 10 NLO PDFs             (Nov 2011)
C   349 AUET2B-NN21 : AUET2 with NNPDF 2.1 NLO PDFs           (Nov 2011)
C   350 Perugia 2011 : Retune of Perugia 2010 incl 7-TeV data (Mar 2011)
C   351 P2011 radHi : Variation with alphaS(pT/2) 
C   352 P2011 radLo : Variation with alphaS(2pT)
C   353 P2011 mpiHi : Variation with more semi-hard MPI
C   354 P2011 noCR  : Variation without color reconnections
C   355 P2011 LO**  : Perugia 2011 using MSTW LO** PDFs       (Mar 2011)
C   356 P2011 C6    : Perugia 2011 using CTEQ6L1 PDFs         (Mar 2011)
C   357 P2011 T16   : Variation with PARP(90)=0.32 away from 7 TeV
C   358 P2011 T32   : Variation with PARP(90)=0.16 awat from 7 TeV
C   359 P2011 TeV   : Perugia 2011 optimized for Tevatron     (Mar 2011)
C   360 S Global    : Schulz-Skands Global fit                (Mar 2011)
C   361 S 7000      : Schulz-Skands at 7000 GeV               (Mar 2011)
C   362 S 1960      : Schulz-Skands at 1960 GeV               (Mar 2011)
C   363 S 1800      : Schulz-Skands at 1800 GeV               (Mar 2011)
C   364 S 900       : Schulz-Skands at 900 GeV                (Mar 2011)
C   365 S 630       : Schulz-Skands at 630 GeV                (Mar 2011)
C
C   370 P12       : Retune of Perugia 2011 w CTEQ6L1          (Oct 2012)
C   371 P12-radHi : Variation with alphaS(pT/2) 
C   372 P12-radLo : Variation with alphaS(2pT)
C   373 P12-mpiHi : Variation with more semi-hard MPI 
C   374 P12-loCR  : Variation using lower CR strength -> more Nch
C   375 P12-noCR  : Variation without any color reconnections
C   376 P12-FL    : Variation with more longitudinal fragmentation
C   377 P12-FT    : Variation with more transverse fragmentation
C   378 P12-M8LO  : Variation using MSTW 2008 LO PDFs     
C   379 P12-LO**  : Variation using MRST LO** PDFs     
C   380 P12-val0  : Variation with PARP(87)=0D0               (Jul 2013)
C   381 P12-ueHi  : Variation with lower pT0 (more soft UE activity)
C   382 P12-ueLo  : Variation with higher pT0 (less soft UE activity)
C   383 P12-IBK   : Perugia 2012 with Innsbruck ee fragmentation parameters

C   390 IBK-CTEQ5L    : Innsbruck pp tune with CTEQ5 LO PDFs  (Jul 2013)
C   391 IBK-CTEQ6LL   :    with CTEQ6LL LO PDFs
C   392 IBK-MSTW08LO  :    with MSTW08 LO PDFS
C   393 IBK-CTEQ66NLO :    with CTEQ6 NLO PDFs
C   394 IBK-CT10NLO   :    with CT10 NLO PDFs
C   395 IBK-MSTW08NLO :    with MSTW08 NLO PDFs
C   396 IBK-MSTW08LO* :    with MSTW07 LO* PDFs
C   397 IBK-MRSTLO**  :    with MRSTMCal (LO**) PDFs
C   398 IBK-CT09MC2   :    with CTEQ09MC2 PDFs  

C ======= The Uppsala models ===========================================
C  1201   SCI 0 : Soft-Colour-Interaction model. Org pars     (Dec 1998)
C  1202   SCI 1 : SCI 0. Tevatron MB retuned (Skands)         (Oct 2006)
C  1401   GAL 0 : Generalized area-law model. Org pars        (Dec 1998)
C  1402   GAL 1 : GAL 0. Tevatron MB retuned (Skands)         (Oct 2006)
C
C More details;
C
C Quick Dictionary:
C      BE : Bose-Einstein
C      BR : Beam Remnants
C      CR : Colour Reconnections
C      HAD: Hadronization
C      ISR/FSR: Initial-State Radiation / Final-State Radiation
C      FSI: Final-State Interactions (=CR+BE)
C      MB : Minimum-bias
C      MI : Multiple Interactions
C      UE : Underlying Event
C
C=======================================================================
C TUNES OF OLD FRAMEWORK (Q2-ORDERED ISR AND FSR, NON-INTERLEAVED UE)
C=======================================================================
C
C   A (100) and AW (101). CTEQ5L parton distributions
C...*** NB : SHOULD BE RUN WITH PYTHIA 6.2 (e.g. 6.228) ***
C...***      CAN ALSO BE RUN WITH PYTHIA 6.406+
C...Key feature: extensively compared to CDF data (R.D. Field).
C...* Large starting scale for ISR (PARP(67)=4)
C...* AW has even more radiation due to smaller mu_R choice in alpha_s.
C...* See: http://www.phys.ufl.edu/~rfield/cdf/
C
C   BW (102). CTEQ5L parton distributions
C...*** NB : SHOULD BE RUN WITH PYTHIA 6.2 (e.g. 6.228) ***
C...***      CAN ALSO BE RUN WITH PYTHIA 6.406+
C...Key feature: extensively compared to CDF data (R.D. Field).
C...NB: Can also be run with Pythia 6.2 or 6.312+
C...* Small starting scale for ISR (PARP(67)=1)
C...* BW has more radiation due to smaller mu_R choice in alpha_s.
C...* See: http://www.phys.ufl.edu/~rfield/cdf/
C
C   DW (103) and DWT (104). CTEQ5L parton distributions
C...*** NB : SHOULD BE RUN WITH PYTHIA 6.2 (e.g. 6.228) ***
C...***      CAN ALSO BE RUN WITH PYTHIA 6.406+
C...Key feature: extensively compared to CDF data (R.D. Field).
C...NB: Can also be run with Pythia 6.2 or 6.312+
C...* Intermediate starting scale for ISR (PARP(67)=2.5)
C...* DWT has a different reference energy, the same as the "S" models
C...  below, leading to more UE activity at the LHC, but less at RHIC.
C...* See: http://www.phys.ufl.edu/~rfield/cdf/
C
C   QW (105). CTEQ61 parton distributions
C...*** NB : SHOULD BE RUN WITH PYTHIA 6.2 (e.g. 6.228) ***
C...***      CAN ALSO BE RUN WITH PYTHIA 6.406+
C...Key feature: uses CTEQ61 (external pdf library must be linked)
C
C   ATLAS-DC2 (106). CTEQ5L parton distributions
C...*** NB : SHOULD BE RUN WITH PYTHIA 6.2 (e.g. 6.228) ***
C...***      CAN ALSO BE RUN WITH PYTHIA 6.406+
C...Key feature: tune used by the ATLAS collaboration.
C
C   ACR (107). CTEQ5L parton distributions
C...*** NB : SHOULD BE RUN WITH PYTHIA 6.412+    ***
C...Key feature: Tune A modified to use annealing CR.
C...NB: PARP(85)=0D0 and amount of CR is regulated by PARP(78).
C
C   D6 (108) and D6T (109). CTEQ6L parton distributions
C...Key feature: Like DW and DWT but retuned to use CTEQ6L PDFs.
C
C   A-Pro, BW-Pro, etc (111, 112, etc). CTEQ5L parton distributions
C   Old UE model, Q2-ordered showers.
C...Key feature: Rick Field's family of tunes revamped with the
C...Professor Q2-ordered final-state shower and fragmentation tunes
C...presented by Hendrik Hoeth at the Perugia MPI workshop in Oct 2008.
C...Key feature: improved descriptions of LEP data.
C
C   Pro-Q2O (129). CTEQ5L parton distributions
C   Old UE model, Q2-ordered showers.
C...Key feature: Complete retune of old model by Professor, including
C...large amounts of both LEP and Tevatron data.
C...Note that PARP(64) (ISR renormalization scale pre-factor) is quite
C...extreme in this tune, corresponding to using mu_R = pT/3 .
C
C=======================================================================
C INTERMEDIATE/HYBRID TUNES (MIX OF NEW AND OLD SHOWER AND UE MODELS)
C=======================================================================
C
C   IM1 (200). Intermediate model, Q2-ordered showers,
C   CTEQ5L parton distributions
C...Key feature: new UE model w Q2-ordered showers and no interleaving.
C...* "Rap" tune of hep-ph/0402078, modified with new annealing CR.
C...* See: Sjostrand & Skands: JHEP 03(2004)053, hep-ph/0402078.
C
C   APT (201). Old UE model, pT-ordered final-state showers,
C   CTEQ5L parton distributions
C...Key feature: Rick Field's Tune A, but with new final-state showers
C
C   APT-Pro (211). Old UE model, pT-ordered final-state showers,
C   CTEQ5L parton distributions
C...Key feature: APT revamped with the Professor pT-ordered final-state
C...shower and fragmentation tunes presented by Hendrik Hoeth at the
C...Perugia MPI workshop in October 2008.
C
C   Perugia-APT (221). Old UE model, pT-ordered final-state showers,
C   CTEQ5L parton distributions
C...Key feature: APT-Pro with final-state showers off the MPI,
C...lower ISR renormalization scale to improve agreement with the
C...Tevatron Drell-Yan pT measurements and with improved energy scaling
C...to min-bias at 630 GeV.
C
C   Perugia-APT6 (226). Old UE model, pT-ordered final-state showers,
C   CTEQ6L1 parton distributions.
C...Key feature: uses CTEQ6L1 (external pdf library must be linked),
C...with a slightly lower pT0 (2.0 instead of 2.05) due to the smaller
C...UE activity obtained with CTEQ6L1 relative to CTEQ5L.
C
C=======================================================================
C TUNES OF NEW FRAMEWORK (PT-ORDERED ISR AND FSR, INTERLEAVED UE)
C=======================================================================
C
C   S0 (300) and S0A (303). CTEQ5L parton distributions
C...Key feature: large amount of multiple interactions
C...* Somewhat faster than the other colour annealing scenarios.
C...* S0A has a faster energy scaling of the UE IR cutoff, borrowed
C...  from Tune A, leading to less UE at the LHC, but more at RHIC.
C...* Small amount of radiation.
C...* Large amount of low-pT MI
C...* Low degree of proton lumpiness (broad matter dist.)
C...* CR Type S (driven by free triplets), of medium strength.
C...* See: Pythia6402 update notes or later.
C
C   S1 (301). CTEQ5L parton distributions
C...Key feature: large amount of radiation.
C...* Large amount of low-pT perturbative ISR
C...* Large amount of FSR off ISR partons
C...* Small amount of low-pT multiple interactions
C...* Moderate degree of proton lumpiness
C...* Least aggressive CR type (S+S Type I), but with large strength
C...* See: Sandhoff & Skands: FERMILAB-CONF-05-518-T, in hep-ph/0604120.
C
C   S2 (302). CTEQ5L parton distributions
C...Key feature: very lumpy proton + gg string cluster formation allowed
C...* Small amount of radiation
C...* Moderate amount of low-pT MI
C...* High degree of proton lumpiness (more spiky matter distribution)
C...* Most aggressive CR type (S+S Type II), but with small strength
C...* See: Sandhoff & Skands: FERMILAB-CONF-05-518-T, in hep-ph/0604120.
C
C   NOCR (304). CTEQ5L parton distributions
C...Key feature: no colour reconnections (NB: "Best fit" only).
C...* NB: <pT>(Nch) problematic in this tune.
C...* Small amount of radiation
C...* Small amount of low-pT MI
C...* Low degree of proton lumpiness
C...* Large BR composite x enhancement factor
C...* Most clever colour flow without CR ("Lambda ordering")
C
C   ATLAS-CSC (306). CTEQ6L parton distributions
C...Key feature: 11-parameter ATLAS tune of the new framework.
C...* Old (pre-annealing) colour reconnections a la 305.
C...* Uses CTEQ6 Leading Order PDFs (must be interfaced externally)
C
C   S0-Pro, S1-Pro, etc (310, 311, etc). CTEQ5L parton distributions.
C...Key feature: the S0 family of tunes revamped with the Professor
C...pT-ordered final-state shower and fragmentation tunes presented by
C...Hendrik Hoeth at the Perugia MPI workshop in October 2008.
C...Key feature: improved descriptions of LEP data.
C
C   ATLAS MC08 (316). CTEQ6L1 parton distributions
C...Key feature: ATLAS tune of the new framework using CTEQ6L1 PDFs
C...* Warning: uses Peterson fragmentation function for heavy quarks
C...* Uses CTEQ6 Leading Order PDFs (must be interfaced externally)
C
C   Perugia-0 (320). CTEQ5L parton distributions.
C...Key feature: S0-Pro retuned to more Tevatron data. Better Drell-Yan
C...pT spectrum, better <pT>(Nch) in min-bias, and better scaling to
C...630 GeV than S0-Pro. Also has a slightly smoother mass profile, more
C...beam-remnant breakup (more baryon number transport), and suppression
C...of CR in high-pT string pieces.
C
C   Perugia-HARD (321). CTEQ5L parton distributions.
C...Key feature: More ISR, More FSR, Less MPI, Less BR
C...Uses pT/2 as argument of alpha_s for ISR, and a higher Lambda_FSR.
C...Has higher pT0, less intrinsic kT, less beam remnant breakup (less
C...baryon number transport), and more fragmentation pT.
C...Multiplicity in min-bias is LOW, <pT>(Nch) is HIGH,
C...DY pT spectrum is HARD.
C
C   Perugia-SOFT (322). CTEQ5L parton distributions.
C...Key feature: Less ISR, Less FSR, More MPI, More BR
C...Uses sqrt(2)*pT as argument of alpha_s for ISR, and a lower
C...Lambda_FSR. Has lower pT0, more beam remnant breakup (more baryon
C...number transport), and less fragmentation pT.
C...Multiplicity in min-bias is HIGH, <pT>(Nch) is LOW,
C...DY pT spectrum is SOFT
C
C   Perugia-3 (323). CTEQ5L parton distributions.
C...Key feature: variant of Perugia-0 with more extreme energy scaling
C...properties while still agreeing with Tevatron data from 630 to 1960.
C...More ISR and less MPI than Perugia-0 at the Tevatron and above and
C...allows FSR off the active end of dipoles stretched to the remnant.
C
C   Perugia-NOCR (324). CTEQ5L parton distributions.
C...Key feature: Retune of NOCR-Pro with better scaling properties to
C...lower energies and somewhat better agreement with Tevatron data
C...at 1800/1960.
C
C   Perugia-* (325). MRST LO* parton distributions for generators
C...Key feature: first attempt at using the LO* distributions
C...(external pdf library must be linked).
C
C   Perugia-6 (326). CTEQ6L1 parton distributions
C...Key feature: uses CTEQ6L1 (external pdf library must be linked).
C
C   Perugia-2010 (327). CTEQ5L parton distributions
C...Key feature: Retune of Perugia 0 to attempt to better describe 
C...strangeness yields at RHIC and at LEP. Also increased the amount 
C...of FSR off ISR following the conclusions in arXiv:1001.4082. 
C...Increased the amount of beam blowup, causing more baryon transport
C...into the detector, to further explore this possibility. Using 
C...a new color-reconnection model that relies on determining a thrust
C...axis for the events and then computing reconnection probabilities for
C...the individual string pieces based on the actual string densities
C...per rapidity interval along that thrust direction.
C
C   Perugia-K (328). CTEQ5L parton distributions 
C...Key feature: uses a ``K'' factor on the MPI cross sections
C...This gives a larger rate of minijets and pushes the underlying-event 
C...activity towards higher pT. To compensate for the increased activity 
C...at higher pT, the infared regularization scale is larger for this tune.
C
C   Pro-pTO (329). CTEQ5L parton distributions
C...Key feature: Complete retune of new model by Professor, including
C...large amounts of both LEP and Tevatron data. Similar to S0A-Pro.
C
C   ATLAS MC09 (330). LO* parton distributions
C...Key feature: Good overall agreement with Tevatron and early LHC data.
C...Similar to Perugia *.
C
C   ATLAS MC09c (331). LO* parton distributions
C...Key feature: Good overall agreement with Tevatron and 900-GeV LHC data.
C...Similar to Perugia *. Retuned CR model with respect to MC09.
C
C   Pro-pT* (335) LO* parton distributions
C...Key feature: Retune of Pro-PTO with MRST LO* PDFs.
C
C   Pro-pT6 (336). CTEQ6L1 parton distributions
C...Key feature: Retune of Pro-PTO with CTEQ6L1 PDFs.
C
C   Pro-pT** (339). LO** parton distributions
C...Key feature: Retune of Pro-PTO with MRST LO** PDFs.
C
C   AMBT1 (340). LO* parton distributions
C...Key feature: First ATLAS tune including 7-TeV LHC data.
C...Mainly retuned CR and mass distribution with respect to MC09c.
C...Note: cannot be run standalone since it uses external PDFs.
C
C   CMSZ1 (341). CTEQ5L parton distributions
C...Key feature: First CMS tune including 7-TeV LHC data.
C...Uses many of the features of AMBT1, but uses CTEQ5L PDFs, 
C...has a lower pT0 at the Tevatron, which scales faster with energy. 
C
C   Z1-LEP (342). CTEQ5L parton distributions
C...Key feature: CMS tune Z1 with improved LEP parameters, mostly 
C...taken from the Professor/Perugia tunes, with a few minor updates.
C
C...More recent Perugia tunes: see arXiv:1005.3457
C
C...Schulz-Skands tunes: see arXiv:1103.3649 

 
C...Global statements
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
 
C...SAVE statements
      SAVE /PYDAT1/,/PYPARS/

C...Internal parameters
      PARAMETER(MXTUNS=500)
      CHARACTER*8 CHDOC
      PARAMETER (CHDOC='Aug 2013')
      CHARACTER*16 CHNAMS(0:MXTUNS), CHNAME
      CHARACTER*42 CHMSTJ(50), CHMSTP(100), CHPARP(100),
     &    CHPARJ(100), CHMSTU(101:121), CHPARU(101:121)
      CHARACTER*60 CH60
      CHARACTER*70 CH70
      DATA (CHNAMS(I),I=0,1)/'Default',' '/
      DATA (CHNAMS(I),I=100,119)/
     &    'Tune A','Tune AW','Tune BW','Tune DW','Tune DWT','Tune QW',
     &    'ATLAS DC2','Tune ACR','Tune D6','Tune D6T',
     1    'Tune A-Pro','Tune AW-Pro','Tune BW-Pro','Tune DW-Pro',
     1    'Tune DWT-Pro','Tune QW-Pro','ATLAS DC2-Pro','Tune ACR-Pro',
     1    'Tune D6-Pro','Tune D6T-Pro'/
      DATA (CHNAMS(I),I=120,129)/
     &     9*' ','Pro-Q2O'/
      DATA (CHNAMS(I),I=130,139)/
     &     'Q12','Q12-radHi','Q12-radLo','Q12-mpiHi','Q12-noCR',
     &     'Q12-M','Q12-F1','Q12-F2','Q12-LE','Q12-TeV'/
      DATA (CHNAMS(I),I=300,309)/
     &    'Tune S0','Tune S1','Tune S2','Tune S0A','NOCR','Old',
     5    'ATLAS-CSC Tune','Yale Tune','Yale-K Tune',' '/
      DATA (CHNAMS(I),I=310,316)/
     &    'Tune S0-Pro','Tune S1-Pro','Tune S2-Pro','Tune S0A-Pro',
     &    'NOCR-Pro','Old-Pro','ATLAS MC08'/
      DATA (CHNAMS(I),I=320,329)/
     &    'Perugia 0','Perugia HARD','Perugia SOFT',
     &    'Perugia 3','Perugia NOCR','Perugia LO*',
     &    'Perugia 6','Perugia 10','Perugia K','Pro-pTO'/
      DATA (CHNAMS(I),I=330,349)/
     &     'ATLAS MC09','ATLAS MC09c',2*' ','Perugia 10 NOCR','Pro-PT*',
     &     'Pro-PT6',' ',' ','Pro-PT**',
     4     'Tune AMBT1','Tune Z1','Tune Z1-LEP','Tune Z2','Tune Z2-LEP',
     4     'AMBT2B-CT6L1','AUET2B-CT6L1','AUET2B-CT66','AUET2B-CT10',
     4     'AUET2B-NN21'/
      DATA (CHNAMS(I),I=350,359)/
     &     'Perugia 2011','P2011 radHi','P2011 radLo','P2011 mpiHi',
     &     'P2011 noCR','P2011 M(LO**)', 'P2011 CTEQ6L1',
     &     'P2011 T16','P2011 T32','P2011 Tevatron'/
      DATA (CHNAMS(I),I=360,369)/
     &     'S Global','S 7000','S 1960','S 1800',
     &     'S 900','S 630', 4*' '/
      DATA (CHNAMS(I),I=370,379)/
     &     'P12','P12-radHi','P12-radLo','P12-mpiHi','P12-loCR',
     &     'P12-noCR','P12-FL','P12-FT','P12-M8LO','P12-LO**'/
      DATA (CHNAMS(I),I=380,399)/
     &     'P12-val0','P12-ueHi','P12-ueLo','P12-IBK',6*' ',
     9     'Innsbruck C5LO','Innsbruck C6LO','Innsbruck M8LO',
     &     'Innsbruck C66NLO','Innsbruck C10NLO',
     &     'Innsbruck M8NLO','Innsbruck LO*','Innsbruck LO**',
     &     'Innsbruck C9MC2',
     &     ' '/
      DATA (CHNAMS(I),I=200,229)/
     &    'IM Tune 1','Tune APT',8*' ',
     &    ' ','Tune APT-Pro',8*' ',
     &    ' ','Perugia APT',4*' ','Perugia APT6',3*' '/
      DATA (CHNAMS(I),I=400,409)/
     &    'GAL Tune 0','SCI Tune 0','GAL Tune 1','SCI Tune 1',6*' '/
      DATA (CHMSTJ(I),I=11,20)/
     &    'HAD choice of fragmentation function(s)',4*' ',
     &    'HAD treatment of small-mass systems',4*' '/
      DATA (CHMSTJ(I),I=41,50)/
     &    'FSR type (Q2 or pT) for old framework',9*' '/
      DATA (CHMSTP(I),I=1,10)/
     &    2*' ','INT switch for choice of LambdaQCD',7*' '/
      DATA (CHMSTP(I),I=31,40)/
     &    2*' ','"K" switch for K-factor on/off & type',7*' '/
      DATA (CHMSTP(I),I=51,100)/
     5    'PDF set','PDF set internal (=1) or pdflib (=2)',8*' ',
     6    'ISR master switch',2*' ','ISR alphaS type',2*' ',
     6    'ISR coherence option for 1st emission',
     6    'ISR phase space choice & ME corrections',' ',
     7    'ISR IR regularization scheme',' ',
     7    'IFSR scheme for non-decay FSR',8*' ',
     8    'UE model',
     8    'UE hadron transverse mass distribution',5*' ',
     8    'BR composite scheme','BR color scheme',
     9    'BR primordial kT compensation',
     9    'BR primordial kT distribution',
     9    'BR energy partitioning scheme',2*' ',
     9    'FSI color (re-)connection model',5*' '/
      DATA (CHPARP(I),I=1,10)/
     &    'ME/UE LambdaQCD',9*' '/
      DATA (CHPARP(I),I=31,40)/
     &    ' ','"K" K-factor',8*' '/
      DATA (CHPARP(I),I=61,100)/
     6     'ISR LambdaQCD','ISR IR cutoff',' ',
     6     'ISR renormalization scale prefactor',
     6     2*' ','ISR Q2max factor',3*' ',
     7     'IFSR Q2max factor in non-s-channel procs',
     7     'IFSR LambdaQCD (outside resonance decays)',4*' ',
     7     'FSI color reco high-pT damping strength',
     7     'FSI color reconnection strength',
     7     'BR composite x enhancement','BR breakup suppression',
     8     2*'UE IR cutoff at reference ecm',
     8     2*'UE mass distribution parameter',
     8     'UE gg color correlated fraction','UE total gg fraction',
     8     'UE qq enhancement at low pT','UE qq enh scale / pT0',
     8     'UE IR cutoff reference ecm',
     8     'UE IR cutoff ecm scaling power',
     9     'BR primordial kT width <|kT|>',' ',
     9     'BR primordial kT UV cutoff',7*' '/
      DATA (CHPARJ(I),I=1,30)/
     &     'HAD diquark suppression','HAD strangeness suppression',
     &     'HAD strange diquark suppression',
     &     'HAD vector diquark suppression','HAD P(popcorn)',
     &     'HAD extra popcorn B(s)-M-B(s) supp',
     &     'HAD extra popcorn B-M(s)-B supp',
     &     3*' ',
     1     'HAD P(vector meson), u and d only',
     1     'HAD P(vector meson), contains s',
     1     'HAD P(vector meson), heavy quarks',
     1     'HAD P(L=1;S=0,J=1)','HAD P(L=1;S=1,J=0)',
     1     'HAD P(L=1;S=1,J=1)','HAD P(L=1;S=1,J=2)', 
     1     'HAD extra spin-3/2 baryon supp',
     1     'HAD extra leading-baryon supp',' ',
     2     'HAD fragmentation pT',' ',' ',' ',
     2     'HAD eta0 suppression',"HAD eta0' suppression",4*' '/
      DATA (CHPARJ(I),I=41,90)/
     4     'HAD string parameter a(Meson)','HAD string parameter b',
     4     2*' ','HAD string a(Baryon)-a(Meson)',
     4     'HAD Lund(=0)-Bowler(=1) rQ (rc)',
     4     'HAD Lund(=0)-Bowler(=1) rb',3*' ',
     5     3*' ', 'HAD charm parameter','HAD bottom parameter',5*' ',
     6     10*' ',10*' ',
     8     'FSR LambdaQCD (inside resonance decays)',
     &     'FSR IR cutoff',8*' '/
      DATA (CHMSTU(I),I=111,120)/
     1     ' ','INT n(flavors) for LambdaQCD',8*' '/
      DATA (CHPARU(I),I=111,120)/
     1     ' ','INT LambdaQCD',8*' '/
      
C...1) Shorthand notation
      M13=MSTU(13)
      M11=MSTU(11)
      IF (MYTUNE.LE.MXTUNS.AND.MYTUNE.GE.0) THEN
        CHNAME=CHNAMS(MYTUNE)
        IF (MYTUNE.EQ.0) GOTO 9999
      ELSE
        CALL PYERRM(9,'(PYTUNE:) Tune number > max. Using defaults.')
        GOTO 9999
      ENDIF
      
C...  2) Hello World
      IF (M13.GE.1) WRITE(M11,5000) CHDOC
      
C...  Hardcode some defaults
C...  Get Lambda from PDF
      MSTP(3)  =  2
C...  CTEQ5L1 PDFs
      MSTP(52) =  1
      MSTP(51) =  7
C...  No K-factor 
      MSTP(33) =  0
C...  Low-pT qq enhancement factor and pT/pT0 ratio
      PARP(87) = 0.7D0
      PARP(88) = 0.5D0
C...  Hard-initialize L=1 meson rates to old default: 0.0
      PARJ(14) = 0D0
      PARJ(15) = 0D0
      PARJ(16) = 0D0
      PARJ(17) = 0D0

C...  3) Tune parameters
      ITUNE = MYTUNE
 
C=======================================================================
C...ATLAS MC08

      IF (ITUNE.EQ.316) THEN
        
        IF (M13.GE.1) WRITE(M11,5010) ITUNE, CHNAME
        IF (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.405))THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ENDIF

C...First set some explicit defaults from 6.4.20
C...# Old defaults
        MSTJ(11) = 4
C...# Old default flavour parameters
        PARJ(1)  =   0.1
        PARJ(2)  =   0.3  
        PARJ(3)  =   0.40 
        PARJ(4)  =   0.05 
        PARJ(11) =   0.5  
        PARJ(12) =   0.6 
        PARJ(21) = 0.36
        PARJ(41) = 0.30
        PARJ(42) = 0.58
        PARJ(46) = 1.0
        PARJ(82) = 1.0

C...PDFs: CTEQ6L1 for 326
        MSTP(52)=2
        MSTP(51)=10042

C...UE and ISR switches
        MSTP(81)=21
        MSTP(82)=4
        MSTP(70)=0
        MSTP(72)=1

C...CR:
        MSTP(95)=2
        PARP(78)=0.3
        PARP(77)=0.0
        PARP(80)=0.1

C...Primordial kT
        PARP(91)=2.0D0
        PARP(93)=5.0D0

C...MPI:
        PARP(82)=2.1
        PARP(83)=0.8
        PARP(84)=0.7
        PARP(89)=1800.0
        PARP(90)=0.16

C...FSR inside resonance decays
        PARJ(81)=0.29

C...Fragmentation (warning: uses Peterson)
        MSTJ(11)=3   
        PARJ(54)=-0.07
        PARJ(55)=-0.006
        
        IF (M13.GE.1) THEN
          CH60='Tuned by ATLAS, ATL-PHYS-PUB-2010-002'
          WRITE(M11,5030) CH60
          CH60='Physics model: '//
     &         'T. Sjostrand & P. Skands, hep-ph/0408302'
          WRITE(M11,5030) CH60
          CH60='CR by P. Skands & D. Wicke, hep-ph/0703081'
          WRITE(M11,5030) CH60
          
C...Output
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          IF (MSTP(70).EQ.0) THEN
            WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          ENDIF
          WRITE(M11,5040) 64, MSTP(64), CHMSTP(64)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5040) 67, MSTP(67), CHMSTP(67)
          WRITE(M11,5050) 67, PARP(67), CHPARP(67)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5040) 70, MSTP(70), CHMSTP(70)
          WRITE(M11,5040) 72, MSTP(72), CHMSTP(72)          
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          WRITE(M11,5060) 82, PARJ(82), CHPARJ(82)
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5040) 88, MSTP(88), CHMSTP(88)
          WRITE(M11,5040) 89, MSTP(89), CHMSTP(89)
          WRITE(M11,5050) 79, PARP(79), CHPARP(79)
          WRITE(M11,5050) 80, PARP(80), CHPARP(80)
          WRITE(M11,5040) 91, MSTP(91), CHMSTP(91)
          WRITE(M11,5050) 91, PARP(91), CHPARP(91)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          IF (MSTP(95).GE.1) THEN
            WRITE(M11,5050) 78, PARP(78), CHPARP(78)
            IF (MSTP(95).GE.2) WRITE(M11,5050) 77, PARP(77), CHPARP(77)
          ENDIF

        ENDIF
 
C=======================================================================
C...ATLAS MC09, MC09c, AMBT1, AMBT2B, AUET2B + NLO PDF vars
C...CMS Z1 (R. Field), Z1-LEP

      ELSEIF (ITUNE.EQ.330.OR.ITUNE.EQ.331.OR.ITUNE.EQ.340.OR.
     &       ITUNE.GE.341.AND.ITUNE.LE.349) THEN
        
        IF (M13.GE.1) WRITE(M11,5010) ITUNE, CHNAME
        IF (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.405))THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ENDIF

C...pT-ordered shower default for everything
        MSTJ(41) = 12

C...FSR inside resonance decays, base value (modified by individual tunes)
        PARJ(81) = 0.29

C...First set some explicit defaults from 6.4.20
        IF (ITUNE.LE.341.OR.ITUNE.EQ.343) THEN
C...  # Old defaults
          MSTJ(11) = 4
C...# Old default flavour parameters
          PARJ(1)  =   0.1
          PARJ(2)  =   0.3  
          PARJ(3)  =   0.40 
          PARJ(4)  =   0.05 
          PARJ(11) =   0.5  
          PARJ(12) =   0.6 
          PARJ(21) = 0.36
          PARJ(41) = 0.30
          PARJ(42) = 0.58
          PARJ(46) = 1.0
          PARJ(82) = 1.0
        ELSE IF (ITUNE.LE.344) THEN
C...# For Zn-LEP tunes, use tuned flavour parameters from Professor/Perugia
          PARJ( 1) = 0.08D0
          PARJ( 2) = 0.21D0
          PARJ( 3) = 0.94
          PARJ( 4) = 0.04D0
          PARJ(11) = 0.35D0
          PARJ(12) = 0.35D0
          PARJ(13) = 0.54
          PARJ(25) = 0.63
          PARJ(26) = 0.12
C...# Switch on Bowler:
          MSTJ(11) = 5
C...# Fragmentation
          PARJ(21) = 0.34D0
          PARJ(41) = 0.35D0
          PARJ(42) = 0.80D0
          PARJ(47) = 1.0
          PARJ(81) = 0.26D0
          PARJ(82) = 1.0D0
        ELSE 
C... A*T2 tunes, from ATL-PHYS-PUB-2011-008
          PARJ( 1) = 0.073
          PARJ( 2) = 0.202
          PARJ( 3) = 0.950
          PARJ( 4) = 0.033
          PARJ(11) = 0.309
          PARJ(12) = 0.402
          PARJ(13) = 0.544
          PARJ(25) = 0.628
          PARJ(26) = 0.129
C...# Switch on Bowler:
          MSTJ(11) = 5
C...  # Fragmentation
          PARJ(21) = 0.30
          PARJ(41) = 0.368
          PARJ(42) = 1.004
          PARJ(47) = 0.873
          PARJ(81) = 0.256
          PARJ(82) = 0.830
        ENDIF

C...Default scales and alphaS choices
        IF (ITUNE.GE.345) THEN
          MSTP(3) = 1
          PARU(112) = 0.192
          PARP(1)   = 0.192
          PARP(61)  = 0.192
        ENDIF

C...PDFs: MRST LO* 
        MSTP(52) = 2
        MSTP(51) = 20650
        IF (ITUNE.EQ.341.OR.ITUNE.EQ.342) THEN
C...Z1 uses CTEQ5L
          MSTP(52) = 1
          MSTP(51) = 7
        ELSEIF (ITUNE.EQ.343.OR.ITUNE.EQ.344) THEN
C...Z2 uses CTEQ6L
          MSTP(52) = 2
          MSTP(51) = 10042
        ELSEIF (ITUNE.EQ.345.OR.ITUNE.EQ.346) THEN 
C...AMBT2B, AUET2B use CTEQ6L1 
          MSTP(52) = 2
          MSTP(51) = 10042          
        ELSEIF (ITUNE.EQ.347) THEN 
C...AUET2B-CT66 uses CTEQ66 NLO PDFs
          MSTP(52) = 2
          MSTP(51) = 10550
        ELSEIF (ITUNE.EQ.348) THEN 
C...AUET2B-CT10 uses CTEQ10 NLO PDFs
          MSTP(52) = 2
          MSTP(51) = 10800
        ELSEIF (ITUNE.EQ.349) THEN 
C...AUET2B-NN21 uses NNPDF 2.1 NLO PDF
          MSTP(52) = 2
          MSTP(51) = 192800
        ENDIF

C...UE and ISR switches
        MSTP(81) = 21
        MSTP(82) = 4
        MSTP(70) = 0
        MSTP(72) = 1

C...CR:
        MSTP(95) = 6
        PARP(78) = 0.3
        PARP(77) = 0.0
        PARP(80) = 0.1
        IF (ITUNE.EQ.331) THEN
          PARP(78) = 0.224          
        ELSEIF (ITUNE.EQ.340) THEN
C...AMBT1
          PARP(77) = 1.016D0
          PARP(78) = 0.538D0
        ELSEIF (ITUNE.GE.341.AND.ITUNE.LE.344) THEN
C...Z1 and Z2 use the AMBT1 CR values
          PARP(77) = 1.016D0
          PARP(78) = 0.538D0
        ELSEIF (ITUNE.EQ.345) THEN
C...AMBT2B
          PARP(77) = 0.357D0
          PARP(78) = 0.235D0
        ELSEIF (ITUNE.EQ.346) THEN
C...AUET2B
          PARP(77) = 0.491D0
          PARP(78) = 0.311D0
        ELSEIF (ITUNE.EQ.347) THEN
C...AUET2B-CT66
          PARP(77) = 0.505D0
          PARP(78) = 0.385D0
        ELSEIF (ITUNE.EQ.348) THEN
C...AUET2B-CT10
          PARP(77) = 0.125D0
          PARP(78) = 0.309D0
        ELSEIF (ITUNE.EQ.349) THEN
C...AUET2B-NN21
          PARP(77) = 0.498D0
          PARP(78) = 0.354D0
        ENDIF

C...MPI:
        PARP(82) = 2.3
        PARP(83) = 0.8
        PARP(84) = 0.7
        PARP(89) = 1800.0
        PARP(90) = 0.25
        IF (ITUNE.EQ.331) THEN
          PARP(82) = 2.315
          PARP(90) = 0.2487
        ELSEIF (ITUNE.EQ.340) THEN
          PARP(82) = 2.292D0
          PARP(83) = 0.356D0
          PARP(84) = 0.651
          PARP(90) = 0.25D0
        ELSEIF (ITUNE.EQ.341.OR.ITUNE.EQ.342) THEN
          PARP(82) = 1.932D0
          PARP(83) = 0.356D0
          PARP(84) = 0.651
          PARP(90) = 0.275D0
        ELSEIF (ITUNE.EQ.343.OR.ITUNE.EQ.344) THEN
          PARP(82) = 1.832D0
          PARP(83) = 0.356D0
          PARP(84) = 0.651
          PARP(90) = 0.275D0
        ELSEIF (ITUNE.EQ.345) THEN
          PARP(82) = 2.34
          PARP(83) = 0.356
          PARP(84) = 0.605
          PARP(90) = 0.246
        ELSEIF (ITUNE.EQ.346) THEN
          PARP(82) = 2.26
          PARP(83) = 0.356
          PARP(84) = 0.443
          PARP(90) = 0.249
        ELSEIF (ITUNE.EQ.347) THEN
          PARP(82) = 1.87
          PARP(83) = 0.356
          PARP(84) = 0.561
          PARP(90) = 0.189
        ELSEIF (ITUNE.EQ.348) THEN
          PARP(82) = 1.89
          PARP(83) = 0.356
          PARP(84) = 0.415
          PARP(90) = 0.182
        ELSEIF (ITUNE.EQ.349) THEN
          PARP(82) = 1.86
          PARP(83) = 0.356
          PARP(84) = 0.588
          PARP(90) = 0.177
        ENDIF
        
C...Primordial kT
        PARP(91) = 2.0D0
        PARP(93) = 5D0
        IF (ITUNE.GE.340) THEN
          PARP(93) = 10D0
        ENDIF
        IF (ITUNE.GE.345) THEN
          PARP(91) = 2.0
        ENDIF

C...ISR
        IF (ITUNE.EQ.345.OR.ITUNE.EQ.346) THEN
          MSTP(64) = 2
          PARP(62) = 1.13
          PARP(64) = 0.68
          PARP(67) = 1.0
        ELSE IF (ITUNE.EQ.347) THEN
          MSTP(64) = 2
          PARP(62) = 0.946
          PARP(64) = 1.032
          PARP(67) = 1.0
        ELSE IF (ITUNE.EQ.348) THEN
          MSTP(64) = 2
          PARP(62) = 0.312
          PARP(64) = 0.939
          PARP(67) = 1.0
        ELSE IF (ITUNE.EQ.349) THEN
          MSTP(64) = 2
          PARP(62) = 1.246
          PARP(64) = 0.771
          PARP(67) = 1.0
        ELSE IF (ITUNE.GE.340) THEN
          PARP(62) = 1.025
        ENDIF

C...FSR off ISR (LambdaQCD) for A*ET2B tunes
        IF (ITUNE.GE.345) THEN
          MSTP(72) = 2
          PARP(72) = 0.527
          IF (ITUNE.EQ.348) THEN
            PARP(72) = 0.537
          ENDIF
        ENDIF

        IF (M13.GE.1) THEN
          IF (ITUNE.LT.340) THEN
            CH60='Tuned by ATLAS, ATL-PHYS-PUB-2010-002'
          ELSEIF (ITUNE.EQ.340) THEN
            CH60='Tuned by ATLAS, ATLAS-CONF-2010-031'
          ELSEIF (ITUNE.EQ.341) THEN
            CH60='AMBT1 Tuned by ATLAS, ATLAS-CONF-2010-031'
            WRITE(M11,5030) CH60
            CH60='Z1 variation tuned by R. D. Field (CMS)'
          ELSEIF (ITUNE.EQ.342) THEN
            CH60='AMBT1 Tuned by ATLAS, ATLAS-CONF-2010-031'
            WRITE(M11,5030) CH60
            CH60='Z1 variation retuned by R. D. Field (CMS)'
            WRITE(M11,5030) CH60
            CH60='Z1-LEP variation retuned by Professor / P. Skands'
          ELSEIF (ITUNE.EQ.343) THEN
            CH60='AMBT1 Tuned by ATLAS, ATLAS-CONF-2010-031'
            WRITE(M11,5030) CH60
            CH60='Z2 variation retuned by R. D. Field (CMS)'
          ELSEIF (ITUNE.EQ.344) THEN
            CH60='AMBT1 Tuned by ATLAS, ATLAS-CONF-2010-031'
            WRITE(M11,5030) CH60
            CH60='Z2 variation retuned by R. D. Field (CMS)'
            WRITE(M11,5030) CH60
            CH60='Z2-LEP variation retuned by Professor / P. Skands'
          ELSEIF (ITUNE.EQ.345.OR.ITUNE.EQ.346) THEN
            CH60='A*T2B tunes by ATLAS, ATL-PHYS-PUB-2011-009'
          ELSEIF (ITUNE.GE.347) THEN
            CH60='A*T2B-NLO tunes by ATLAS, ATL-PHYS-PUB-2011-014'
            WRITE(M11,5030) CH60
            CH60='Warning: NLO PDFs are NOT recommended!'
          ENDIF
          WRITE(M11,5030) CH60
          CH60='Physics Model: '//
     &         'T. Sjostrand & P. Skands, hep-ph/0408302'
          WRITE(M11,5030) CH60
          CH60='CR by P. Skands & D. Wicke, hep-ph/0703081'
          WRITE(M11,5030) CH60

C...Output
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          IF (MSTP(3).EQ.1) THEN
            WRITE(M11,6100) 112, MSTU(112), CHMSTU(112)
            WRITE(M11,6110) 112, PARU(112), CHPARU(112)
            WRITE(M11,5050)   1, PARP(1)  , CHPARP(  1)
          ENDIF
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          IF (MSTP(3).EQ.1) THEN
            WRITE(M11,5050)  72, PARP(72) , CHPARP( 72)
            WRITE(M11,5050)  61, PARP(61) , CHPARP( 61)
          ENDIF
          WRITE(M11,5040) 64, MSTP(64), CHMSTP(64)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5040) 67, MSTP(67), CHMSTP(67)
          WRITE(M11,5050) 67, PARP(67), CHPARP(67)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5040) 70, MSTP(70), CHMSTP(70)
          IF (MSTP(70).EQ.0) THEN
            WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          ENDIF
          WRITE(M11,5040) 72, MSTP(72), CHMSTP(72)
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5050) 72, PARP(72), CHPARP(72)
          WRITE(M11,5060) 82, PARJ(82), CHPARJ(82)
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5040) 88, MSTP(88), CHMSTP(88)
          WRITE(M11,5040) 89, MSTP(89), CHMSTP(89)
          WRITE(M11,5050) 79, PARP(79), CHPARP(79)
          WRITE(M11,5050) 80, PARP(80), CHPARP(80)
          WRITE(M11,5040) 91, MSTP(91), CHMSTP(91)
          WRITE(M11,5050) 91, PARP(91), CHPARP(91)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          IF (MSTP(95).GE.1) THEN
            WRITE(M11,5050) 78, PARP(78), CHPARP(78)
            IF (MSTP(95).GE.2) WRITE(M11,5050) 77, PARP(77), CHPARP(77)
          ENDIF

        ENDIF

C=======================================================================
C...S0, S1, S2, S0A, NOCR, Rap,
C...S0-Pro, S1-Pro, S2-Pro, S0A-Pro, NOCR-Pro, Rap-Pro
C...Perugia 0, HARD, SOFT, 3, LO*, 6, 2010, K
C...Pro-pTO, Pro-PT*, Pro-PT6, Pro-PT**
C...Perugia 2011 (incl variations)
C...Schulz-Skands tunes
      ELSEIF ((ITUNE.GE.300.AND.ITUNE.LE.305)
     &    .OR.(ITUNE.GE.310.AND.ITUNE.LE.315)
     &    .OR.(ITUNE.GE.320.AND.ITUNE.LE.329)
     &    .OR.(ITUNE.GE.334.AND.ITUNE.LE.336).OR.ITUNE.EQ.339
     &    .OR.(ITUNE.GE.350.AND.ITUNE.LE.389)) THEN
        IF (M13.GE.1) WRITE(M11,5010) ITUNE, CHNAME
        IF (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.405))THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ELSEIF(ITUNE.GE.320.AND.ITUNE.LE.339.AND.ITUNE.NE.324.AND.
     &         ITUNE.NE.334.AND.
     &        (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.419)))
     &        THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ELSEIF((ITUNE.EQ.327.OR.ITUNE.EQ.328.OR.ITUNE.GE.350).AND.
     &         (MSTP(181).LE.5.OR.
     &         (MSTP(181).EQ.6.AND.MSTP(182).LE.422)))
     &        THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ENDIF
 
C...Use 327 as base tune for 350-359 and 370-379 (Perugia 2011 and 2012)
        ITUNSV = ITUNE
        IF (ITUNE.GE.350.AND.ITUNE.LE.359) ITUNE = 327
        IF (ITUNE.GE.370.AND.ITUNE.LE.389) ITUNE = 327
C...Use 320 as base tune for 360+ (Schulz-Skands)
        IF (ITUNE.GE.360) ITUNE = 320

C...HAD: Use Professor's LEP pars if ITUNE >= 310
C...(i.e., for S0-Pro, S1-Pro etc, and for Perugia tunes)
        IF (ITUNE.LT.310) THEN
C...# Old defaults
          MSTJ(11) = 4
C...# Old default flavour parameters
          PARJ(1)  =   0.1
          PARJ(2)  =   0.3  
          PARJ(3)  =   0.40 
          PARJ(4)  =   0.05 
          PARJ(11) =   0.5  
          PARJ(12) =   0.6 
          PARJ(21) = 0.36
          PARJ(41) = 0.30
          PARJ(42) = 0.58
          PARJ(46) = 1.0
          PARJ(82) = 1.0
          
        ELSEIF (ITUNE.GE.310) THEN
C...# Tuned flavour parameters:
          PARJ(1)  = 0.073
          PARJ(2)  = 0.2
          PARJ(3)  = 0.94
          PARJ(4)  = 0.032
          PARJ(11) = 0.31
          PARJ(12) = 0.4
          PARJ(13) = 0.54
          PARJ(25) = 0.63
          PARJ(26) = 0.12
C...# Always use pT-ordered shower:
          MSTJ(41) = 12
C...# Switch on Bowler:
          MSTJ(11) = 5
C...# Fragmentation
          PARJ(21) = 0.313
          PARJ(41) = 0.49
          PARJ(42) = 1.2
          PARJ(47) = 1.0
          PARJ(81) = 0.257
          PARJ(82) = 0.8

C...HAD: fragmentation pT (only if not using professor) - HARD and SOFT
          IF (ITUNE.EQ.321) PARJ(21) = 0.34D0
          IF (ITUNE.EQ.322) PARJ(21) = 0.28D0

C...HAD: P-2010 and P-K use different strangeness parameters 
C...     indicated by LEP and RHIC yields.
C...(only 5% different from Professor values, so should be within acceptable
C...theoretical uncertainty range)
C...(No attempt made to retune other flavor parameters post facto)
          IF (ITUNE.EQ.327.OR.ITUNE.EQ.328.OR.ITUNE.EQ.334) THEN
            PARJ( 1) = 0.08D0
            PARJ( 2) = 0.21D0
            PARJ( 4) = 0.04D0
            PARJ(11) = 0.35D0
            PARJ(12) = 0.35D0
            PARJ(21) = 0.36D0
            PARJ(41) = 0.35D0
            PARJ(42) = 0.90D0
            PARJ(81) = 0.26D0
            PARJ(82) = 1.0D0
          ENDIF 
        ENDIF
 
C...Remove middle digit now for Professor variants, since identical pars
        ITUNEB=ITUNE
        IF (ITUNE.GE.310.AND.ITUNE.LE.319) THEN
          ITUNEB=(ITUNE/100)*100+MOD(ITUNE,10)
        ENDIF
 
C...PDFs: all use CTEQ5L as starting point
        MSTP(52) = 1
        MSTP(51) = 7
        IF (ITUNE.EQ.325.OR.ITUNE.EQ.335) THEN
C...MRST LO* for 325 and 335
          MSTP(52) = 2
          MSTP(51) = 20650
        ELSEIF (ITUNE.EQ.326.OR.ITUNE.EQ.336) THEN
C...CTEQ6L1 for 326 and 336
          MSTP(52) = 2
          MSTP(51) = 10042
        ELSEIF (ITUNE.EQ.339) THEN
C...MRST LO** for 339
          MSTP(52) = 2
          MSTP(51) = 20651
        ENDIF
 
C...LambdaQCD choice: 327 and 328 use hardcoded, others get from PDF
        MSTP(3) = 2
        IF (ITUNE.EQ.327.OR.ITUNE.EQ.328.OR.ITUNE.EQ.334) THEN
          MSTP(3)   = 1
C...Hardcode CTEQ5L values for ME and ISR
          MSTU(112) = 4
          PARU(112) = 0.192D0
          PARP(61)  = 0.192D0
          PARP( 1)  = 0.192D0
C...but use LEP value also for non-res FSR
          PARP(72)  = 0.260D0
        ENDIF

C...ISR: use Lambda_MSbar with default scale for S0(A)
        MSTP(64) = 2
        PARP(64) = 1D0
        IF (ITUNE.EQ.320.OR.ITUNE.EQ.323.OR.ITUNE.EQ.324.OR.ITUNE.EQ.334
     &       .OR.ITUNE.EQ.326.OR.ITUNE.EQ.327.OR.ITUNE.EQ.328) THEN
C...Use Lambda_MC with muR^2=pT^2 for most central Perugia tunes
          MSTP(64) = 3
          PARP(64) = 1D0
        ELSEIF (ITUNE.EQ.321) THEN
C...Use Lambda_MC with muR^2=(1/2pT)^2 for Perugia HARD
          MSTP(64) = 3
          PARP(64) = 0.25D0
        ELSEIF (ITUNE.EQ.322) THEN
C...Use Lambda_MSbar with muR^2=2pT^2 for Perugia SOFT
          MSTP(64) = 2
          PARP(64) = 2D0
        ELSEIF (ITUNE.EQ.325) THEN
C...Use Lambda_MC with muR^2=2pT^2 for Perugia LO*
          MSTP(64) = 3
          PARP(64) = 2D0
        ELSEIF (ITUNE.EQ.329.OR.ITUNE.EQ.335.OR.ITUNE.EQ.336.OR.
     &         ITUNE.EQ.339) THEN
C...Use Lambda_MSbar with P64=1.3 for Pro-pT0
          MSTP(64) = 2
          PARP(64) = 1.3D0
          IF (ITUNE.EQ.335) PARP(64) = 0.92D0
          IF (ITUNE.EQ.336) PARP(64) = 0.89D0
          IF (ITUNE.EQ.339) PARP(64) = 0.97D0
        ENDIF
 
C...ISR : power-suppressed power showers above s_color (since 6.4.19)
        MSTP(67) = 2
        PARP(67) = 4D0
C...Perugia tunes have stronger suppression, except HARD
        IF ((ITUNE.GE.320.AND.ITUNE.LE.328).OR.ITUNE.EQ.334) THEN
          PARP(67) = 1D0
          IF (ITUNE.EQ.321) PARP(67) = 4D0
          IF (ITUNE.EQ.322) PARP(67) = 0.25D0
        ENDIF
 
C...ISR IR cutoff type and FSR off ISR setting:
C...Smooth ISR, low FSR-off-ISR
        MSTP(70) = 2
        MSTP(72) = 0
        IF (ITUNEB.EQ.301) THEN
C...S1, S1-Pro: sharp ISR, high FSR
          MSTP(70) = 0
          MSTP(72) = 1
        ELSEIF (ITUNE.EQ.320.OR.ITUNE.EQ.324.OR.ITUNE.EQ.326
     &        .OR.ITUNE.EQ.325) THEN
C...Perugia default is smooth ISR, high FSR-off-ISR
          MSTP(70) = 2
          MSTP(72) = 1
        ELSEIF (ITUNE.EQ.321) THEN
C...Perugia HARD: sharp ISR, high FSR-off-ISR (but no dip-to-BR rad)
          MSTP(70) = 0
          PARP(62) = 1.25D0
          MSTP(72) = 1
        ELSEIF (ITUNE.EQ.322) THEN
C...Perugia SOFT: scaling sharp ISR, low FSR-off-ISR
          MSTP(70) = 1
          PARP(81) = 1.5D0
          MSTP(72) = 0
        ELSEIF (ITUNE.EQ.323) THEN
C...Perugia 3: sharp ISR, high FSR-off-ISR (with dipole-to-BR radiating)
          MSTP(70) = 0
          PARP(62) = 1.25D0
          MSTP(72) = 2
        ELSEIF (ITUNE.EQ.327.OR.ITUNE.EQ.328.OR.ITUNE.EQ.334) THEN
C...Perugia 2010/K: smooth ISR, high FSR-off-ISR (with dipole-to-BR radiating)
          MSTP(70) = 2
          MSTP(72) = 2
        ENDIF
 
C...FSR activity: Perugia tunes use a lower PARP(71) as indicated 
C...by Professor tunes (with HARD and SOFT variations)
        PARP(71) = 4D0
        IF ((ITUNE.GE.320.AND.ITUNE.LE.328).OR.ITUNE.EQ.334) THEN 
          PARP(71) = 2D0
          IF (ITUNE.EQ.321) PARP(71) = 4D0
          IF (ITUNE.EQ.322) PARP(71) = 1D0
        ENDIF
        IF (ITUNE.EQ.329) PARP(71) = 2D0
        IF (ITUNE.EQ.335) PARP(71) = 1.29D0
        IF (ITUNE.EQ.336) PARP(71) = 1.72D0
        IF (ITUNE.EQ.339) PARP(71) = 1.20D0

C...FSR: Lambda_FSR scale (only if not using professor)
        IF (ITUNE.LT.310) PARJ(81) = 0.23D0
        IF (ITUNE.EQ.321) PARJ(81) = 0.30D0
        IF (ITUNE.EQ.322) PARJ(81) = 0.20D0

C...K-factor : only 328 uses a K-factor on the UE cross sections
        MSTP(33) = 0
        IF (ITUNE.EQ.328) THEN
          MSTP(33) = 10
          PARP(32) = 1.5
        ENDIF
C...UE on, new model
        MSTP(81) = 21
 
C...UE: hadron-hadron overlap profile (expOfPow for all)
        MSTP(82) = 5
C...UE: Overlap smoothness (1.0 = exponential; 2.0 = gaussian)
        PARP(83) = 1.6D0
        IF (ITUNEB.EQ.301) PARP(83) = 1.4D0
        IF (ITUNEB.EQ.302) PARP(83) = 1.2D0
C...NOCR variants have very smooth distributions
        IF (ITUNEB.EQ.304) PARP(83) = 1.8D0
        IF (ITUNEB.EQ.305) PARP(83) = 2.0D0
        IF ((ITUNE.GE.320.AND.ITUNE.LE.328).OR.ITUNE.EQ.334) THEN
C...Perugia variants have slightly smoother profiles by default
C...(to compensate for more tail by added radiation)
C...Perugia-SOFT has more peaked distribution, NOCR less peaked
          PARP(83) = 1.7D0
          IF (ITUNE.EQ.322) PARP(83) = 1.5D0
          IF (ITUNE.EQ.327) PARP(83) = 1.5D0
          IF (ITUNE.EQ.328) PARP(83) = 1.5D0
C...NOCR variants have smoother mass profiles
          IF (ITUNE.EQ.324) PARP(83) = 1.8D0
          IF (ITUNE.EQ.334) PARP(83) = 1.8D0
        ENDIF
C...Professor-pT0 also has very smooth distribution
        IF (ITUNE.EQ.329) PARP(83) = 1.8
        IF (ITUNE.EQ.335) PARP(83) = 1.68
        IF (ITUNE.EQ.336) PARP(83) = 1.72
        IF (ITUNE.EQ.339) PARP(83) = 1.67

C...UE: pT0 = 1.85 for S0, S0A, 2.0 for Perugia version
        PARP(82) = 1.85D0
        IF (ITUNEB.EQ.301) PARP(82) = 2.1D0
        IF (ITUNEB.EQ.302) PARP(82) = 1.9D0
        IF (ITUNEB.EQ.304) PARP(82) = 2.05D0
        IF (ITUNEB.EQ.305) PARP(82) = 1.9D0
        IF ((ITUNE.GE.320.AND.ITUNE.LE.328).OR.ITUNE.EQ.334) THEN
C...Perugia tunes (def is 2.0 GeV, HARD has higher, SOFT has lower,
C...Perugia-3 has more ISR, so higher pT0, NOCR can be slightly lower,
C...CTEQ6L1 slightly lower, due to less activity, and LO* needs to be
C...slightly higher, due to increased activity.
          PARP(82) = 2.0D0
          IF (ITUNE.EQ.321) PARP(82) = 2.3D0
          IF (ITUNE.EQ.322) PARP(82) = 1.9D0
          IF (ITUNE.EQ.323) PARP(82) = 2.2D0
          IF (ITUNE.EQ.324) PARP(82) = 1.95D0
          IF (ITUNE.EQ.325) PARP(82) = 2.2D0
          IF (ITUNE.EQ.326) PARP(82) = 1.95D0
          IF (ITUNE.EQ.327) PARP(82) = 2.05D0
          IF (ITUNE.EQ.328) PARP(82) = 2.45D0
          IF (ITUNE.EQ.334) PARP(82) = 2.15D0
        ENDIF
C...Professor-pT0 maintains low pT0 vaue
        IF (ITUNE.EQ.329) PARP(82) = 1.85D0
        IF (ITUNE.EQ.335) PARP(82) = 2.10D0
        IF (ITUNE.EQ.336) PARP(82) = 1.83D0
        IF (ITUNE.EQ.339) PARP(82) = 2.28D0

C...UE: IR cutoff reference energy and default energy scaling pace
        PARP(89) = 1800D0
        PARP(90) = 0.16D0
C...S0A, S0A-Pro have tune A energy scaling
        IF (ITUNEB.EQ.303) PARP(90) = 0.25D0
        IF ((ITUNE.GE.320.AND.ITUNE.LE.328).OR.ITUNE.EQ.334) THEN
C...Perugia tunes explicitly include MB at 630 to fix energy scaling
          PARP(90) = 0.26
          IF (ITUNE.EQ.321) PARP(90) = 0.30D0
          IF (ITUNE.EQ.322) PARP(90) = 0.24D0
          IF (ITUNE.EQ.323) PARP(90) = 0.32D0
          IF (ITUNE.EQ.324) PARP(90) = 0.24D0
C...LO* and CTEQ6L1 tunes have slower energy scaling
          IF (ITUNE.EQ.325) PARP(90) = 0.23D0
          IF (ITUNE.EQ.326) PARP(90) = 0.22D0
        ENDIF
C...Professor-pT0 has intermediate scaling
        IF (ITUNE.EQ.329) PARP(90) = 0.22D0
        IF (ITUNE.EQ.335) PARP(90) = 0.20D0
        IF (ITUNE.EQ.336) PARP(90) = 0.20D0
        IF (ITUNE.EQ.339) PARP(90) = 0.21D0

C...BR: MPI initiator color connections rap-ordered by default
C...NOCR variants are Lambda-ordered, Perugia SOFT & 2010 random-ordered
        MSTP(89) = 1
        IF (ITUNEB.EQ.304.OR.ITUNE.EQ.324) MSTP(89) = 2
        IF (ITUNE.EQ.322) MSTP(89) = 0
        IF (ITUNE.EQ.327) MSTP(89) = 0
        IF (ITUNE.EQ.328) MSTP(89) = 0
 
C...BR: BR-g-BR suppression factor (higher values -> more beam blowup)
        PARP(80) = 0.01D0
        IF (ITUNE.GE.320.AND.ITUNE.LE.328) THEN
C...Perugia tunes have more beam blowup by default
          PARP(80) = 0.05D0
          IF (ITUNE.EQ.321) PARP(80) = 0.01
          IF (ITUNE.EQ.323) PARP(80) = 0.03
          IF (ITUNE.EQ.324) PARP(80) = 0.01
          IF (ITUNE.EQ.327) PARP(80) = 0.1
          IF (ITUNE.EQ.328) PARP(80) = 0.1
        ENDIF
 
C...BR: diquarks (def = valence qq and moderate diquark x enhancement)
        MSTP(88) = 0
        PARP(79) = 2D0
        IF (ITUNEB.EQ.304) PARP(79) = 3D0
        IF (ITUNE.EQ.329) PARP(79) = 1.18
        IF (ITUNE.EQ.335) PARP(79) = 1.11
        IF (ITUNE.EQ.336) PARP(79) = 1.10
        IF (ITUNE.EQ.339) PARP(79) = 3.69

C...BR: Primordial kT, parametrization and cutoff, default is 2 GeV
        MSTP(91) = 1
        PARP(91) = 2D0
        PARP(93) = 10D0
C...Perugia-HARD only uses 1.0 GeV
        IF (ITUNE.EQ.321) PARP(91) = 1.0D0
C...Perugia-3 only uses 1.5 GeV
        IF (ITUNE.EQ.323) PARP(91) = 1.5D0
C...Professor-pT0 uses 7-GeV cutoff
        IF (ITUNE.EQ.329) PARP(93) = 7.0
        IF (ITUNE.EQ.335) THEN
          PARP(91) = 2.15
          PARP(93) = 6.79
        ELSEIF (ITUNE.EQ.336) THEN
          PARP(91) = 1.85
          PARP(93) = 6.86
        ELSEIF (ITUNE.EQ.339) THEN
          PARP(91) = 2.11
          PARP(93) = 5.08
        ENDIF

C...FSI: Colour Reconnections - Seattle algorithm is default (S0)
        MSTP(95) = 6
C...S1, S1-Pro: use S1
        IF (ITUNEB.EQ.301) MSTP(95) = 2
C...S2, S2-Pro: use S2
        IF (ITUNEB.EQ.302) MSTP(95) = 4
C...NOCR, NOCR-Pro, Perugia-NOCR: use no CR
        IF (ITUNE.EQ.304.OR.ITUNE.EQ.314.OR.ITUNE.EQ.324.OR.
     &       ITUNE.EQ.334) MSTP(95) = 0
C..."Old" and "Old"-Pro: use old CR
        IF (ITUNEB.EQ.305) MSTP(95) = 1
C...Perugia 2010 and K use Paquis model
        IF (ITUNE.EQ.327.OR.ITUNE.EQ.328) MSTP(95) = 8
 
C...FSI: CR strength and high-pT dampening, default is S0
        PARP(77) = 0D0
        IF (ITUNE.LT.320.OR.ITUNE.EQ.329.OR.ITUNE.GE.335) THEN
          PARP(78) = 0.2D0
          IF (ITUNEB.EQ.301) PARP(78) = 0.35D0
          IF (ITUNEB.EQ.302) PARP(78) = 0.15D0
          IF (ITUNEB.EQ.304) PARP(78) = 0.0D0
          IF (ITUNEB.EQ.305) PARP(78) = 1.0D0
          IF (ITUNE.EQ.329) PARP(78) = 0.17D0
          IF (ITUNE.EQ.335) PARP(78) = 0.14D0
          IF (ITUNE.EQ.336) PARP(78) = 0.17D0
          IF (ITUNE.EQ.339) PARP(78) = 0.13D0
        ELSE
C...Perugia tunes also use high-pT dampening : default is Perugia 0,*,6
          PARP(78) = 0.33
          PARP(77) = 0.9D0
          IF (ITUNE.EQ.321) THEN
C...HARD has HIGH amount of CR
            PARP(78) = 0.37D0
            PARP(77) = 0.4D0
          ELSEIF (ITUNE.EQ.322) THEN
C...SOFT has LOW amount of CR
            PARP(78) = 0.15D0
            PARP(77) = 0.5D0
          ELSEIF (ITUNE.EQ.323) THEN
C...Scaling variant appears to need slightly more than default
            PARP(78) = 0.35D0
            PARP(77) = 0.6D0
          ELSEIF (ITUNE.EQ.324.OR.ITUNE.EQ.334) THEN
C...NOCR has no CR
            PARP(78) = 0D0
            PARP(77) = 0D0
          ELSEIF (ITUNE.EQ.327) THEN
C...2010
            PARP(78) = 0.035D0
            PARP(77) = 1D0
          ELSEIF (ITUNE.EQ.328) THEN
C...K
            PARP(78) = 0.033D0
            PARP(77) = 1D0
          ENDIF
        ENDIF
 
C================
C...Perugia 2011 and 2012 tunes 
C...(written as modifications on top of Perugia 2010)
C================
        IF ( (ITUNSV.GE.350.AND.ITUNSV.LE.359) 
     &       .OR.(ITUNSV.GE.370.AND.ITUNSV.LE.389) ) THEN
          ITUNE = ITUNSV
C...  Scale setting for matching applications.
C...  Switch to 5-flavor CMW LambdaQCD = 0.26 for all shower activity
C...  (equivalent to a 5-flavor MSbar LambdaQCD = 0.26/1.6 = 0.16)
          MSTP(64) = 2
          MSTU(112) = 5
C...  This sets the Lambda scale for ISR, IFSR, and FSR
          PARP(61) = 0.26D0
          PARP(72) = 0.26D0
          PARJ(81) = 0.26D0
C...  This sets the Lambda scale for QCD hard interactions (important for the 
C...  UE dijet cross sections. Here we still use an MSbar value, rather than 
C...  a CMW one, in order not to hugely increase the UE jettiness. The CTEQ5L
C...  value corresponds to a Lambda5 of 0.146 for comparison, so quite close.)
          PARP(1) = 0.16D0
          PARU(112) = 0.16D0
C...  For matching applications, PARP(71) and PARP(67) = 1
          PARP(67) = 1D0
          PARP(71) = 1D0
C...  Primordial kT: only use 1 GeV
          MSTP(91) = 1
          PARP(91) = 1D0
C...  ADDITIONAL LESSONS WRT PERUGIA 2010
C...  ALICE taught us: need less baryon transport than SOFT
          MSTP(89) = 0
          PARP(80) = 0.015
C...  Small adjustments at LEP (slightly softer frag functions, esp for baryons)
          PARJ(21) = 0.33
          PARJ(41) = 0.35
          PARJ(42) = 0.8
          PARJ(45) = 0.55
C...  Increase Lambda/K ratio and other strange baryon yields 
          PARJ(1) = 0.087D0
          PARJ(3) = 0.95D0
          PARJ(4) = 0.043D0
          PARJ(6) = 1.0D0
          PARJ(7) = 1.0D0
C...  Also reduce total strangeness yield a bit, with higher K*/K
          PARJ(2) = 0.19D0
          PARJ(12) = 0.40D0
C...  Perugia 2011 default is sharp ISR, dipoles to BR radiating, pTmax individual
          MSTP(70) = 0
          MSTP(72) = 2
          PARP(62) = 1.5D0
C...  Holger taught us a smoother proton is preferred at high energies
C...  Just use a simple Gaussian 
          MSTP(82) = 3
C...  Scaling of pt0 cutoff
          PARP(90) = 0.265
C...  Now retune pT0 to give right UE activity.
C...  Low CR strength indicated by LHC tunes 
C...  (also keep low to get <pT>(Nch) a bit down for pT>100MeV samples)
          PARP(78) = 0.036D0
C...  Choose 7 TeV as new reference scale
          PARP(89) = 7000.0D0
          PARP(82) = 2.93D0          
C================
C...  P2011 Variations
C================
          IF (ITUNE.EQ.351) THEN
C...  radHi: high Lambda scale for ISR, IFSR, and FSR
C...  ( ca 10% more particles at LEP after retune )
            PARP(61) = 0.52D0
            PARP(72) = 0.52D0
            PARJ(81) = 0.52D0
C...  Retune cutoff scales to compensate partially
C...  (though higher cutoff causes faster multiplicity drop at low energies)
            PARP(62) = 1.75D0
            PARJ(82) = 1.75D0
            PARP(82) = 3.00D0
C...  Needs faster cutoff scaling than nominal variant for same <Nch> scaling
C...  (since more radiation otherwise generates faster mult growth)
            PARP(90) = 0.28  
          ELSEIF (ITUNE.EQ.352) THEN
C...  radLo: low Lambda scale for ISR, IFSR, and FSR
C...  ( ca 10% less particles at LEP after retune )
            PARP(61) = 0.13D0
            PARP(72) = 0.13D0
            PARJ(81) = 0.13D0
C...  Retune cutoff scales to compensate partially
            PARP(62) = 1.00D0
            PARJ(82) = 0.75D0
            PARP(82) = 2.95D0 
C...  Needs slower cutoff scaling than nominal variant for same <Nch> scaling
C...  (since less radiation otherwise generates slower mult growth)
            PARP(90) = 0.24
          ELSEIF (ITUNE.EQ.353) THEN
C...  mpiHi: high Lambda scale for MPI
            PARP(1) = 0.26D0
            PARU(112) = 0.26D0
            PARP(82) = 3.35D0
            PARP(90) = 0.26D0
          ELSEIF (ITUNE.EQ.354) THEN
            MSTP(95) = 0
            PARP(82) = 3.05D0
          ELSEIF (ITUNE.EQ.355) THEN
C...  LO**
            MSTP(52) = 2
            MSTP(51) = 20651
            PARP(62) = 1.5D0
C...  Compensate for higher <pT> with less CR
            PARP(78) = 0.034
            PARP(82) = 3.40D0 
C...  Need slower energy scaling than CTEQ5L
            PARP(90) = 0.23D0 
          ELSEIF (ITUNE.EQ.356) THEN
C...  CTEQ6L1
            MSTP(52) = 2
            MSTP(51) = 10042
            PARP(82) = 2.65D0
C...  Need slower cutoff scaling than CTEQ5L
            PARP(90) = 0.22D0 
          ELSEIF (ITUNE.EQ.357) THEN
C...  T16
            PARP(90) = 0.16
          ELSEIF (ITUNE.EQ.358) THEN
C...  T32
            PARP(90) = 0.32
          ELSEIF (ITUNE.EQ.359) THEN
C...  Tevatron
            PARP(89) = 1800D0
            PARP(90) = 0.28 
            PARP(82) = 2.10 
            PARP(78) = 0.05 
          ENDIF
          
C================
C...  Perugia 2012 Variations
C================
          IF (ITUNE.GE.370) THEN
C...  CTEQ6L1 Baseline
            MSTP(52) = 2
            MSTP(51) = 10042
            PARP(82) = 2.65D0
C...  Needs slower cutoff scaling than CTEQ5L
            PARP(90) = 0.24D0 
C...  Slightly lower CR strength than Perugia 2011
            PARP(78) = 0.035D0
C...  Adjusted fragmentation parameters wrt 2011
            PARJ(1)  = 0.085D0
            PARJ(2)  = 0.2
            PARJ(3)  = 0.92
            PARJ(25) = 0.70
            PARJ(26) = 0.135
            PARJ(41) = 0.45
            PARJ(42) = 1.0
            PARJ(45) = 0.86
          ENDIF
C... Variations
          IF (ITUNE.EQ.371) THEN
C...  radHi: high Lambda scale for ISR, IFSR, and FSR
C...  ( ca 10% more particles at LEP after retune )
            PARP(61) = 0.52D0
            PARP(72) = 0.52D0
            PARJ(81) = 0.52D0
C...  Retune cutoff scales to compensate partially
C...  (though higher cutoff causes faster multiplicity drop at low energies)
            PARP(62) = 1.75D0
            PARJ(82) = 1.75D0
            PARP(82) = 2.725D0
C...  Needs faster cutoff scaling than nominal variant for same <Nch> scaling
C...  (since more radiation otherwise generates faster mult growth)
            PARP(90) = 0.25
          ELSEIF (ITUNE.EQ.372) THEN
C...  radLo: low Lambda scale for ISR, IFSR, and FSR
C...  ( ca 10% less particles at LEP after retune )
            PARP(61) = 0.13D0
            PARP(72) = 0.13D0
            PARJ(81) = 0.13D0
C...  Retune cutoff scales to compensate partially
            PARP(62) = 1.00D0
            PARJ(82) = 0.75D0
            PARP(82) = 2.6D0 
C...  Needs slower cutoff scaling than nominal variant for same <Nch> scaling
C...  (since less radiation otherwise generates slower mult growth)
            PARP(90) = 0.23
          ELSEIF (ITUNE.EQ.373) THEN
C...  mpiHi: high Lambda scale for MPI
            PARP(1) = 0.26D0
            PARU(112) = 0.26D0
            PARP(82) = 3.0D0
            PARP(90) = 0.24D0
          ELSEIF (ITUNE.EQ.374) THEN
C... LOCR : uses global CR model. Less extreme alternative to noCR. 
            MSTP(95) = 6
            PARP(78) = 0.25D0
            PARP(82) = 2.7D0
            PARP(83) = 1.50D0
            PARP(90) = 0.24
          ELSEIF (ITUNE.EQ.375) THEN
C... NOCR : with higher pT0
            MSTP(95) = 0
            PARP(82) = 2.80D0
          ELSEIF (ITUNE.EQ.376) THEN
C... hadF1 (harder frag function, smaller n.p. pT)
            PARJ(21) = 0.30
            PARJ(41) = 0.36
            PARJ(42) = 1.0
            PARJ(45) = 0.75
          ELSEIF (ITUNE.EQ.377) THEN
C... hadF2 (softer frag function, larger n.p. pT)
            PARJ(21) = 0.36 
            PARJ(41) = 0.45
            PARJ(42) = 0.75
            PARJ(45) = 0.9
          ELSEIF (ITUNE.EQ.378) THEN
C... MSTW08LO
            MSTP(52) = 2
            MSTP(51) = 21000
            PARP(82) = 2.9D0 
C...Uses a large LambdaQCD MSbar value (close to CMW one)
C...(Nominally, MSTW 2008 alphaS(mZ) = 0.139)
            PARP(1) = 0.26D0
            PARU(112) = 0.26D0
C...Tentative (fast) energy scaling
            PARP(90) = 0.29
          ELSEIF (ITUNE.EQ.379) THEN
C... MSTW LO**
            MSTP(52) = 2
            MSTP(51) = 20651
            PARP(62) = 1.5D0
C... Use a smaller LambdaQCD MSbar than with CTEQ
            PARP(1) = 0.14D0
            PARU(112) = 0.14D0
C...  Compensate for higher <pT> with less CR
            PARP(78) = 0.034
            PARP(82) = 3.25D0 
C...Tentative scaling
            PARP(90) = 0.25
          ELSEIF (ITUNE.EQ.380) THEN
C...  val0: remove artificial valence-domination of low-pT scatterings
C...  slightly faster energy scaling of pT0 cutoff (slower mult growth)
            PARP(87)=0D0
            PARP(90)=0.245
          ELSEIF (ITUNE.EQ.381) THEN
C...  ueHi: lower pT0 value, slower pT0 scaling
            PARP(82)=2.46D0   
            PARP(90)=0.23
          ELSEIF (ITUNE.EQ.382) THEN
C...  ueLo: higher pT0 value, faster pT0 scaling
            PARP(82)=2.92D0
            PARP(90)=0.26
          ELSEIF (ITUNE.EQ.383) THEN
C...  IBK: same as Perugia 2012, but with Innsbruck ee fragm parameters 
C...  Different Lambdas
            MSTP(3)  =  1
C...  Lund+Bowler scheme for HQ fragment. 
            MSTJ(11) = 5
C...  old baryon model     
            MSTJ(12) = 2
C...  2=PYSHOW  12=PYPTFS for gluon and photon emiss.
            MSTJ(41) = 12
C...  Lambda_LLA  
            PARJ(81) = 0.261  
C...  p_tmin cutoff (set by hand)             
            PARJ(82) = 0.90    
C...  sigma_pt
            PARJ(21) = 0.329   
C...  A of LSFF
            PARJ(41) = 0.425   
C...  B of LSFF
            PARJ(42) = 1.65    
C...  r_c
            PARJ(46) = 1.42    
C...  r_b
            PARJ(47) = 0.975   
C...  reset popcorn parameters 
            PARJ( 6) = 0.5
            PARJ( 7) = 0.5
C...  V_u,d
            PARJ(11) = 0.549   
C...  V_s
            PARJ(12) = 0.450   
C...  V_c,b
            PARJ(13) = 0.500   
C...  L=1 mesons rates
            PARJ(17) = 0.20    
            PARJ(14) = 0.12   
            PARJ(15) = 0.04   
            PARJ(16) = 0.12   
C...  eta suppr.
            PARJ(25) = 1.000
C...  eta-prime suppr.
            PARJ(26) = 0.245   
C...  s/u
            PARJ( 2) = 0.268   
C...  qq/q
            PARJ( 1) = 0.128   
C...  su/du
            PARJ( 3) = 0.772   
C...  (qq)_1
            PARJ( 4) = 0.05    
C...  end-point baryon suppress.           
            PARJ(19) = 0.402   
C...  reset a(Baryon)-a(Meson) parameter to default value
            PARJ(45) = 0.50
          ENDIF 
C================
C...Schulz-Skands 2011 tunes 
C...(written as modifications on top of Perugia 0)
C================
        ELSEIF (ITUNSV.GE.360.AND.ITUNSV.LE.365) THEN
          ITUNE = ITUNSV

          IF (ITUNE.EQ.360) THEN
            PARP(78) = 0.40D0
            PARP(82) = 2.19D0
            PARP(83) = 1.45D0
            PARP(89) = 1800.0D0
            PARP(90) = 0.27D0
          ELSEIF (ITUNE.EQ.361) THEN
            PARP(78) = 0.20D0
            PARP(82) = 2.75D0
            PARP(83) = 1.73D0
            PARP(89) = 7000.0D0
          ELSEIF (ITUNE.EQ.362) THEN
            PARP(78) = 0.31D0
            PARP(82) = 1.97D0
            PARP(83) = 1.98D0
            PARP(89) = 1960.0D0
          ELSEIF (ITUNE.EQ.363) THEN
            PARP(78) = 0.35D0
            PARP(82) = 1.91D0
            PARP(83) = 2.02D0
            PARP(89) = 1800.0D0
          ELSEIF (ITUNE.EQ.364) THEN
            PARP(78) = 0.33D0
            PARP(82) = 1.69D0
            PARP(83) = 1.92D0
            PARP(89) = 900.0D0
          ELSEIF (ITUNE.EQ.365) THEN
            PARP(78) = 0.47D0
            PARP(82) = 1.61D0
            PARP(83) = 1.50D0
            PARP(89) = 630.0D0
          ENDIF

        ENDIF
        
C...Switch off trial joinings
        MSTP(96) = 0
 
C...S0 (300), S0A (303)
        IF (ITUNEB.EQ.300.OR.ITUNEB.EQ.303) THEN
          IF (M13.GE.1) THEN
            CH60='see P. Skands & D. Wicke, hep-ph/0703081'
            WRITE(M11,5030) CH60
            CH60='M. Sandhoff & P. Skands, in hep-ph/0604120'
            WRITE(M11,5030) CH60
            CH60='and T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            IF (ITUNE.GE.310) THEN
              CH60='LEP parameters tuned by Professor,'//
     &             ' hep-ph/0907.2973'
              WRITE(M11,5030) CH60
            ENDIF
          ENDIF
 
C...S1 (301)
        ELSEIF(ITUNEB.EQ.301) THEN
          IF (M13.GE.1) THEN
            CH60='see M. Sandhoff & P. Skands, in hep-ph/0604120'
            WRITE(M11,5030) CH60
            CH60='and T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            IF (ITUNE.GE.310) THEN
              CH60='LEP parameters tuned by Professor,'//
     &             ' hep-ph/0907.2973'
              WRITE(M11,5030) CH60
            ENDIF
          ENDIF
 
C...S2 (302)
        ELSEIF(ITUNEB.EQ.302) THEN
          IF (M13.GE.1) THEN
            CH60='see M. Sandhoff & P. Skands, in hep-ph/0604120'
            WRITE(M11,5030) CH60
            CH60='and T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            IF (ITUNE.GE.310) THEN
              CH60='LEP parameters tuned by Professor,'//
     &             ' hep-ph/0907.2973'
              WRITE(M11,5030) CH60
            ENDIF
          ENDIF
 
C...NOCR (304)
        ELSEIF(ITUNEB.EQ.304) THEN
          IF (M13.GE.1) THEN
            CH60='"best try" without colour reconnections'
            WRITE(M11,5030) CH60
            CH60='see P. Skands & D. Wicke, hep-ph/0703081'
            WRITE(M11,5030) CH60
            CH60='and T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            IF (ITUNE.GE.310) THEN
              CH60='LEP parameters tuned by Professor,'//
     &             ' hep-ph/0907.2973'
              WRITE(M11,5030) CH60
            ENDIF
          ENDIF
 
C..."Lo FSR" retune (305)
        ELSEIF(ITUNEB.EQ.305) THEN
          IF (M13.GE.1) THEN
            CH60='"Lo FSR retune" with primitive colour reconnections'
            WRITE(M11,5030) CH60
            CH60='see T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            IF (ITUNE.GE.310) THEN
              CH60='LEP parameters tuned by Professor,'//
     &             ' hep-ph/0907.2973'
              WRITE(M11,5030) CH60
            ENDIF
          ENDIF
 
C...Perugia Tunes (320-328 and 334)
        ELSEIF((ITUNE.GE.320.AND.ITUNE.LE.328).OR.ITUNE.EQ.334) THEN
          IF (M13.GE.1) THEN
            CH60='Tuned by P. Skands, hep-ph/1005.3457'
            WRITE(M11,5030) CH60
            CH60='Physics Model: '//
     &           'T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            IF (ITUNE.LE.326) THEN
              CH60='CR by P. Skands & D. Wicke, hep-ph/0703081'
              WRITE(M11,5030) CH60
              CH60='LEP parameters tuned by Professor, hep-ph/0907.2973'
              WRITE(M11,5030) CH60
            ENDIF
            IF (ITUNE.EQ.325) THEN
              CH70='NB! This tune requires MRST LO* pdfs to be '//
     &            'externally linked'
              WRITE(M11,5035) CH70
            ELSEIF (ITUNE.EQ.326) THEN
              CH70='NB! This tune requires CTEQ6L1 pdfs to be '//
     &            'externally linked'
              WRITE(M11,5035) CH70
            ELSEIF (ITUNE.EQ.321) THEN
              CH60='NB! This tune has MORE ISR & FSR / LESS UE & BR'
              WRITE(M11,5030) CH60
            ELSEIF (ITUNE.EQ.322) THEN
              CH60='NB! This tune has LESS ISR & FSR / MORE UE & BR'
              WRITE(M11,5030) CH60
            ENDIF
          ENDIF
 
C...Professor-pTO (329)
        ELSEIF(ITUNE.EQ.329.OR.ITUNE.EQ.335.OR.ITUNE.EQ.336.OR.
     &         ITUNE.EQ.339) THEN
          IF (M13.GE.1) THEN
            CH60='Tuned by Professor, hep-ph/0907.2973'
            WRITE(M11,5030) CH60 
            CH60='Physics Model: '//
     &           'T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            CH60='CR by P. Skands & D. Wicke, hep-ph/0703081'
            WRITE(M11,5030) CH60
          ENDIF
 
C...Perugia 2011 Tunes (350-359)
        ELSEIF(ITUNE.GE.350.AND.ITUNE.LE.359) THEN
          IF (M13.GE.1) THEN
            CH60='Tuned by P. Skands, hep-ph/1005.3457'
            WRITE(M11,5030) CH60
            CH60='Physics Model: '//
     &           'T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            CH60='CR by P. Skands & D. Wicke, hep-ph/0703081'
            WRITE(M11,5030) CH60
            IF (ITUNE.EQ.355) THEN
              CH70='NB! This tune requires MRST LO** pdfs to be '//
     &            'externally linked'
              WRITE(M11,5035) CH70
            ELSEIF (ITUNE.EQ.356) THEN
              CH70='NB! This tune requires CTEQ6L1 pdfs to be '//
     &            'externally linked'
              WRITE(M11,5035) CH70
            ENDIF
          ENDIF

C...Schulz-Skands Tunes (360-365)
        ELSEIF(ITUNE.GE.360.AND.ITUNE.LE.365) THEN
          IF (M13.GE.1) THEN
            CH60='Tuned by H. Schulz & P. Skands, MCNET-11-07'
            WRITE(M11,5030) CH60
            CH60='Based on Perugia 0, hep-ph/1005.3457'
            WRITE(M11,5030) CH60
            CH60='Physics Model: '//
     &           'T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            CH60='CR by P. Skands & D. Wicke, hep-ph/0703081'
            WRITE(M11,5030) CH60
          ENDIF
 
C...Perugia 2012 Tunes (370-389)
        ELSEIF(ITUNE.GE.370.AND.ITUNE.LE.389) THEN
          IF (M13.GE.1) THEN
            CH60='Tuned by P. Skands, hep-ph/1005.3457'
            WRITE(M11,5030) CH60
            IF (ITUNE.EQ.383) THEN
              CH60='with Innsbruck (IBK) ee fragmentation parameters'
              WRITE(M11,5030) CH60
            ENDIF
            CH60='Physics Model: '//
     &           'T. Sjostrand & P. Skands, hep-ph/0408302'
            WRITE(M11,5030) CH60
            CH60='CR by P. Skands & D. Wicke, hep-ph/0703081'
            WRITE(M11,5030) CH60
            IF (ITUNE.EQ.378) THEN
            ELSEIF (ITUNE.EQ.379) THEN
              CH70='NB! This tune requires MRST 2008 LO** pdfs to be '//
     &            'externally linked'
              WRITE(M11,5035) CH70
            ELSE 
              CH70='NB! This tune requires CTEQ6L1 pdfs to be '//
     &            'externally linked'
              WRITE(M11,5035) CH70
            ENDIF
          ENDIF

        ENDIF
 
C...Output
        IF (M13.GE.1) THEN
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          IF (MSTP(33).GE.10) THEN
            WRITE(M11,5050) 32, PARP(32), CHPARP(32)
          ENDIF
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          IF (MSTP(3).EQ.1) THEN
            WRITE(M11,6100) 112, MSTU(112), CHMSTU(112)
            WRITE(M11,6110) 112, PARU(112), CHPARU(112)
            WRITE(M11,5050)   1, PARP(1)  , CHPARP(  1)
          ENDIF
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          IF (MSTP(3).EQ.1) THEN 
            WRITE(M11,5050)  72, PARP(72) , CHPARP( 72)
            WRITE(M11,5050)  61, PARP(61) , CHPARP( 61)
          ENDIF
          WRITE(M11,5040) 64, MSTP(64), CHMSTP(64)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5040) 67, MSTP(67), CHMSTP(67)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5050) 67, PARP(67), CHPARP(67)
          WRITE(M11,5040) 72, MSTP(72), CHMSTP(72)
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5040) 70, MSTP(70), CHMSTP(70)
          IF (MSTP(70).EQ.0) THEN
            WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          ELSEIF (MSTP(70).EQ.1) THEN
            WRITE(M11,5050) 81, PARP(81), CHPARP(62)
            CH60='(Note: PARP(81) replaces PARP(62).)'
            WRITE(M11,5030) CH60
          ENDIF
          WRITE(M11,5060) 82, PARJ(82), CHPARJ(82)
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          IF (MSTP(70).EQ.2) THEN
            CH60='(Note: PARP(82) replaces PARP(62).)'
            WRITE(M11,5030) CH60
          ENDIF
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          IF (MSTP(82).EQ.5) THEN
            WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          ELSEIF (MSTP(82).EQ.4) THEN
            WRITE(M11,5050) 83, PARP(83), CHPARP(83)
            WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          ENDIF
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5040) 88, MSTP(88), CHMSTP(88)
          WRITE(M11,5040) 89, MSTP(89), CHMSTP(89)
          WRITE(M11,5050) 79, PARP(79), CHPARP(79)
          WRITE(M11,5050) 80, PARP(80), CHPARP(80)
          WRITE(M11,5040) 91, MSTP(91), CHMSTP(91)
          WRITE(M11,5050) 91, PARP(91), CHPARP(91)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          IF (MSTP(95).GE.1) THEN
            WRITE(M11,5050) 78, PARP(78), CHPARP(78)
            IF (MSTP(95).GE.2) WRITE(M11,5050) 77, PARP(77), CHPARP(77)
          ENDIF

        ENDIF
 
C=======================================================================
C...Innsbruck tunes (provided by N. Firdous and G. Rudolph, Innsbruck)
C...390-395
      ELSEIF (ITUNE.GE.390.AND.ITUNE.LE.395) THEN 
        IF (M13.GE.1) WRITE(M11,5010) ITUNE, CHNAME
        IF (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.419))THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ENDIF

C...  1) Set the IBK ee fragmentation parameters (March 2012)
C...  Lund+Bowler scheme for HQ fragment. 
        MSTJ(11) = 5
C...  old baryon model     
        MSTJ(12) = 2
C...  2=PYSHOW  12=PYPTFS for gluon and photon emiss.
        MSTJ(41) = 12
C...  Lambda_LLA  
        PARJ(81) = 0.261  
C...  p_tmin cutoff (set by hand)             
        PARJ(82) = 0.90    
C...  sigma_pt
        PARJ(21) = 0.329   
C...  A of LSFF
        PARJ(41) = 0.425   
C...  B of LSFF
        PARJ(42) = 1.65    
C...  r_c
        PARJ(46) = 1.42    
C...  r_b
        PARJ(47) = 0.975   
C...  V_u,d
        PARJ(11) = 0.549   
C...  V_s
        PARJ(12) = 0.450   
C...  V_c,b
        PARJ(13) = 0.500   
C...  L=1 mesons rates
        PARJ(17) = 0.20    
        PARJ(14) = 0.12   
        PARJ(15) = 0.04   
        PARJ(16) = 0.12   
C...  eta suppr.
        PARJ(25) = 1.000
C...  eta-prime suppr.
        PARJ(26) = 0.245   
C...  s/u
        PARJ( 2) = 0.268   
C...  qq/q
        PARJ( 1) = 0.128   
C...  su/du
        PARJ( 3) = 0.772   
C...  (qq)_1
        PARJ( 4) = 0.05    
C...  end-point baryon suppress.           
        PARJ(19) = 0.402   
C...  reset a(Baryon)-a(Meson) parameter to default value
        PARJ(45) = 0.50

C...  2) Set the global IBK pp tune parameters        
C...  Different Lambda_QCD    
        MSTP(  3) = 1
C...  N_flavors = 5
        MSTU(112) = 5
C...  MPI & BR master switch   
        MSTP( 81) = 21
C...  alpha_s(Q**2) choice in ISR  (def=2)                  
        MSTP( 64) = 2
C...  ISR regularisation (def=1)
        MSTP( 70) = 2
C...  ptmax scale for rad betw ISR partons (def=1)
        MSTP( 72) = 2
C...  MPI structure: matter overlap (def=4) 
        MSTP( 82) = 5
C...  collapse of junction configur. (def=1) 
        MSTP( 88) = 0
C...  CR: annealing model (def=1) 
        MSTP( 95) = 6 
C...  Lam_QCD for ISR 
        PARP( 61) = 0.190
C...  K-factor in alpha_s for ISR (def=1.) 
        PARP( 64) = 1.0
C...  max.virt. scale factor for ISR  (def=4.)                  
        PARP( 67) = 1.0
C...  max.virt. scale factor for FSR (def=4.) 
        PARP( 71) = 1.0
C...  CR suppression for fast moving strings (def=0.)
        PARP( 77) = 0.90
C...  PT0 reference Ecm (def=1800 GeV)
        PARP( 89) = 7000.0
C...  beam remnant x enhancement  (def=2.)              
        PARP( 79) = 1.50
C...  beam remnant breakup suppression (def=0.1)        
        PARP( 80) = 0.06
C...  intrinsic kT width (def=2.0) 
        PARP( 91) = 2.0
C...  intrinsic kT cutoff(def=5.0) 
        PARP( 93) = 10.0

C...  3) Set the tune-specific IBK pp tune parameters
        IF (ITUNE.EQ.390) THEN
C...  CTEQ5L          
          MSTP(51)=7  
          MSTP(52)=1
          PARP(82)=2.942
          PARP(90)=0.2450
          PARP(83)=1.817
          PARP(78)=0.433 
          PARP( 1)=0.163
          PARU(112)=0.163
          PARP(72)=0.531
        ELSEIF (ITUNE.EQ.391) THEN
C...  CTEQ6LL
          MSTP(51)=10042
          MSTP(52)=2
          PARP(82)=2.625
          PARP(90)=0.2178
          PARP(83)=1.863
          PARP(78)=0.461 
          PARP( 1)=0.141
          PARU(112)=0.141 
          PARP(72)=0.475
        ELSEIF (ITUNE.EQ.392) THEN
C...  MSTW08LO
          MSTP(51)=21000
          MSTP(52)=2
          PARP(82)=2.889
          PARP(90)=0.2832
          PARP(83)=1.785
          PARP(78)=0.478 
          PARP( 1)=0.199
          PARU(112)=0.199
          PARP(72)=0.657            
        ELSEIF (ITUNE.EQ.393) THEN
C...  CTEQ66 NLO
          MSTP(51)=10550
          MSTP(52)=2
          PARP(82)=2.172
          PARP(90)=0.1818
          PARP(83)=1.939
          PARP(78)=0.513 
          PARP( 1)=0.173
          PARU(112)=0.173 
          PARP(72)=0.456
        ELSEIF (ITUNE.EQ.394) THEN
C...  CT10 NLO
          MSTP(51)=10800
          MSTP(52)=2
          PARP(82)=2.090
          PARP(90)=0.1687
          PARP(83)=1.939
          PARP(78)=0.517 
          PARP( 1)=0.177 
          PARU(112)=0.177
          PARP(72)=0.463 
        ELSEIF (ITUNE.EQ.395) THEN
C...  MSTW08NLO
          MSTP(51)=21100
          MSTP(52)=2
          PARP(82)=1.773
          PARP(90)=0.1780
          PARP(83)=1.882
          PARP(78)=0.590 
          PARP( 1)=0.161 
          PARU(112)=0.161
          PARP(72)=0.367 
        ELSEIF (ITUNE.EQ.396) THEN
C...  MRST07LO* 
          MSTP(51)=20650
          MSTP(52)=2
          PARP(82)=2.619
          PARP(90)=0.2286
          PARP(83)=1.812
          PARP(78)=0.471 
          PARP( 1)=0.082 
          PARU(112)=0.082 
          PARP(72)=0.500
        ELSEIF (ITUNE.EQ.397) THEN
C...  MRSTMCal (LO**)
          MSTP(51)=20651
          MSTP(52)=2
          PARP(82)=2.802
          PARP(90)=0.2220
          PARP(83)=1.821
          PARP(78)=0.441 
          PARP( 1)=0.080
          PARU(112)=0.080 
          PARP(72)=0.519
        ELSEIF (ITUNE.EQ.398) THEN
C...CT09MC2 
          MSTP(51)=10772
          MSTP(52)=2
          PARP(82)=2.355
          PARP(90)=0.2062
          PARP(83)=1.893
          PARP(78)=0.509 
          PARP( 1)=0.058 
          PARU(112)=0.058 
          PARP(72)=0.401
        ENDIF

C...Output
        IF (M13.GE.1) THEN
          CH60='Tune provided by N. Firdous & G. Rudolph (Innsbruck)'
          WRITE(M11,5030) CH60
          CH60='Physics Model: '//
     &         'T. Sjostrand & P. Skands, hep-ph/0408302'
          WRITE(M11,5030) CH60
          CH60='CR by P. Skands & D. Wicke, hep-ph/0703081'
          WRITE(M11,5030) CH60
          IF (ITUNE.GE.391) THEN
            CH70='NB ! This tune requires LHAPDF to be '//
     &           'externally linked'
            WRITE(M11,5035) CH70
          ENDIF
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          IF (MSTP(33).GE.10) THEN
            WRITE(M11,5050) 32, PARP(32), CHPARP(32)
          ENDIF
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          IF (MSTP(3).EQ.1) THEN
            WRITE(M11,6100) 112, MSTU(112), CHMSTU(112)
            WRITE(M11,6110) 112, PARU(112), CHPARU(112)
            WRITE(M11,5050)   1, PARP(1)  , CHPARP(  1)
          ENDIF
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          IF (MSTP(3).EQ.1) THEN 
            WRITE(M11,5050)  72, PARP(72) , CHPARP( 72)
            WRITE(M11,5050)  61, PARP(61) , CHPARP( 61)
          ENDIF
          WRITE(M11,5040) 64, MSTP(64), CHMSTP(64)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5040) 67, MSTP(67), CHMSTP(67)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5050) 67, PARP(67), CHPARP(67)
          WRITE(M11,5040) 72, MSTP(72), CHMSTP(72)
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5040) 70, MSTP(70), CHMSTP(70)
          IF (MSTP(70).EQ.0) THEN
            WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          ELSEIF (MSTP(70).EQ.1) THEN
            WRITE(M11,5050) 81, PARP(81), CHPARP(62)
            CH60='(Note: PARP(81) replaces PARP(62).)'
            WRITE(M11,5030) CH60
          ENDIF
          WRITE(M11,5060) 82, PARJ(82), CHPARJ(82)
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          IF (MSTP(70).EQ.2) THEN
            CH60='(Note: PARP(82) replaces PARP(62).)'
            WRITE(M11,5030) CH60
          ENDIF
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          IF (MSTP(82).EQ.5) THEN
            WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          ELSEIF (MSTP(82).EQ.4) THEN
            WRITE(M11,5050) 83, PARP(83), CHPARP(83)
            WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          ENDIF
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5040) 88, MSTP(88), CHMSTP(88)
          WRITE(M11,5040) 89, MSTP(89), CHMSTP(89)
          WRITE(M11,5050) 79, PARP(79), CHPARP(79)
          WRITE(M11,5050) 80, PARP(80), CHPARP(80)
          WRITE(M11,5040) 91, MSTP(91), CHMSTP(91)
          WRITE(M11,5050) 91, PARP(91), CHPARP(91)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          IF (MSTP(95).GE.1) THEN
            WRITE(M11,5050) 78, PARP(78), CHPARP(78)
            IF (MSTP(95).GE.2) WRITE(M11,5050) 77, PARP(77), CHPARP(77)
          ENDIF

        ENDIF
C=======================================================================
C...ATLAS-CSC 11-parameter tune (By A. Moraes)
      ELSEIF (ITUNE.EQ.306) THEN
        IF (M13.GE.1) WRITE(M11,5010) ITUNE, CHNAME
        IF (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.405))THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ENDIF
 
C...PDFs
        MSTP(52) = 2
        MSTP(54) = 2
        MSTP(51) = 10042
        MSTP(53) = 10042
C...ISR
C        PARP(64) = 1D0
C...UE on, new model.
        MSTP(81) = 21
C...Energy scaling
        PARP(89) = 1800D0
        PARP(90) = 0.22D0
C...Switch off trial joinings
        MSTP(96) = 0
C...Primordial kT cutoff
 
        IF (M13.GE.1) THEN
          CH60='see presentations by A. Moraes (ATLAS),'
          WRITE(M11,5030) CH60
          CH60='and T. Sjostrand & P. Skands, hep-ph/0408302'
          WRITE(M11,5030) CH60
          WRITE(M11,5030) ' '
          CH70='NB! This tune requires CTEQ6.1 pdfs to be '//
     &        'externally linked'
          WRITE(M11,5035) CH70
        ENDIF
C...Smooth ISR, low FSR
        MSTP(70) = 2
        MSTP(72) = 0
C...pT0
        PARP(82) = 1.9D0
C...Transverse density profile.
        MSTP(82) = 4
        PARP(83) = 0.3D0
        PARP(84) = 0.5D0
C...ISR & FSR in interactions after the first (default)
        MSTP(84) = 1
        MSTP(85) = 1
C...No double-counting (default)
        MSTP(86) = 2
C...Companion quark parent gluon (1-x) power
        MSTP(87) = 4
C...Primordial kT compensation along chaings (default = 0 : uniform)
        MSTP(90) = 1
C...Colour Reconnections
        MSTP(95) = 1
        PARP(78) = 0.2D0
C...Lambda_FSR scale.
        PARJ(81) = 0.23D0
C...Rap order, Valence qq, qq x enhc, BR-g-BR supp
        MSTP(89) = 1
        MSTP(88) = 0
C   PARP(79) = 2D0
        PARP(80) = 0.01D0
C...Peterson charm frag, and c and b hadr parameters
        MSTJ(11) = 3
        PARJ(54) = -0.07
        PARJ(55) = -0.006
C...  Output
        IF (M13.GE.1) THEN
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5040) 70, MSTP(70), CHMSTP(70)
          WRITE(M11,5040) 72, MSTP(72), CHMSTP(72)
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          CH60='(Note: PARJ(81) changed from 0.14! See update notes)'
          WRITE(M11,5030) CH60
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5040) 88, MSTP(88), CHMSTP(88)
          WRITE(M11,5040) 89, MSTP(89), CHMSTP(89)
          WRITE(M11,5040) 90, MSTP(90), CHMSTP(90)
          WRITE(M11,5050) 79, PARP(79), CHPARP(79)
          WRITE(M11,5050) 80, PARP(80), CHPARP(80)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          WRITE(M11,5050) 78, PARP(78), CHPARP(78)

        ENDIF
 
C=======================================================================
C...Tunes A, AW, BW, DW, DWT, QW, D6, D6T (by R.D. Field, CDF)
C...(100-105,108-109), ATLAS-DC2 Tune (by A. Moraes, ATLAS) (106)
C...A-Pro, DW-Pro, etc (100-119), and Pro-Q2O (129)
      ELSEIF ((ITUNE.GE.100.AND.ITUNE.LE.106).OR.ITUNE.EQ.108.OR.
     &      ITUNE.EQ.109.OR.(ITUNE.GE.110.AND.ITUNE.LE.116).OR.
     &      ITUNE.EQ.118.OR.ITUNE.EQ.119.OR.ITUNE.EQ.129) THEN
        IF (M13.GE.1.AND.ITUNE.NE.106.AND.ITUNE.NE.129) THEN
          WRITE(M11,5010) ITUNE, CHNAME
          CH60='see R.D. Field, in hep-ph/0610012'
          WRITE(M11,5030) CH60
          CH60='and T. Sjostrand & M. v. Zijl, PRD36(1987)2019'
          WRITE(M11,5030) CH60
          IF (ITUNE.GE.110.AND.ITUNE.LE.119) THEN
            CH60='LEP parameters tuned by Professor, hep-ph/0907.2973'
            WRITE(M11,5030) CH60
          ENDIF
        ELSEIF (M13.GE.1.AND.ITUNE.EQ.129) THEN
          WRITE(M11,5010) ITUNE, CHNAME
          CH60='Tuned by Professor, hep-ph/0907.2973'
          WRITE(M11,5030) CH60
          CH60='Physics Model: '//
     &         'T. Sjostrand & M. v. Zijl, PRD36(1987)2019'
          WRITE(M11,5030) CH60
        ENDIF
 
C...Make sure we start from old default fragmentation parameters
        PARJ(81) = 0.29
        PARJ(82) = 1.0
 
C...Use Professor's LEP pars if ITUNE >= 110
C...(i.e., for A-Pro, DW-Pro etc)
        IF (ITUNE.LT.110) THEN
C...# Old defaults
          MSTJ(11) = 4
          PARJ(1)  =   0.1
          PARJ(2)  =   0.3  
          PARJ(3)  =   0.40 
          PARJ(4)  =   0.05 
          PARJ(11) =   0.5  
          PARJ(12) =   0.6 
          PARJ(21) = 0.36
          PARJ(41) = 0.30
          PARJ(42) = 0.58
          PARJ(46) = 1.0
          PARJ(81) = 0.29
          PARJ(82) = 1.0
        ELSE
C...# Tuned flavour parameters:
          PARJ(1)  = 0.073
          PARJ(2)  = 0.2
          PARJ(3)  = 0.94
          PARJ(4)  = 0.032
          PARJ(11) = 0.31
          PARJ(12) = 0.4
          PARJ(13) = 0.54
          PARJ(25) = 0.63
          PARJ(26) = 0.12
C...# Switch on Bowler:
          MSTJ(11) = 5
C...# Fragmentation
          PARJ(21) = 0.325
          PARJ(41) = 0.5
          PARJ(42) = 0.6
          PARJ(47) = 0.67
          PARJ(81) = 0.29
          PARJ(82) = 1.65
        ENDIF
 
C...Remove middle digit now for Professor variants, since identical pars
        ITUNEB=ITUNE
        IF (ITUNE.GE.110.AND.ITUNE.LE.119) THEN
          ITUNEB=(ITUNE/100)*100+MOD(ITUNE,10)
        ENDIF
 
C...Multiple interactions on, old framework
        MSTP(81) = 1
C...Fast IR cutoff energy scaling by default
        PARP(89) = 1800D0
        PARP(90) = 0.25D0
C...Default CTEQ5L (internal), except for QW: CTEQ61 (external)
        MSTP(51) = 7
        MSTP(52) = 1
        IF (ITUNEB.EQ.105) THEN
          MSTP(51) = 10150
          MSTP(52) = 2
        ELSEIF(ITUNEB.EQ.108.OR.ITUNEB.EQ.109) THEN
          MSTP(52) = 2
          MSTP(54) = 2
          MSTP(51) = 10042
          MSTP(53) = 10042
        ENDIF
C...Double Gaussian matter distribution.
        MSTP(82) = 4
        PARP(83) = 0.5D0
        PARP(84) = 0.4D0
C...FSR activity.
        PARP(71) = 4D0
C...Fragmentation functions and c and b parameters
C...(only if not using Professor)
        IF (ITUNE.LE.109) THEN
          MSTJ(11) = 4
          PARJ(54) = -0.05
          PARJ(55) = -0.005
        ENDIF
 
C...Tune A and AW
        IF(ITUNEB.EQ.100.OR.ITUNEB.EQ.101) THEN
C...pT0.
          PARP(82) = 2.0D0
c...String drawing almost completely minimizes string length.
          PARP(85) = 0.9D0
          PARP(86) = 0.95D0
C...ISR cutoff, muR scale factor, and phase space size
          PARP(62) = 1D0
          PARP(64) = 1D0
          PARP(67) = 4D0
C...Intrinsic kT, size, and max
          MSTP(91) = 1
          PARP(91) = 1D0
          PARP(93) = 5D0
C...AW : higher ISR IR cutoff, but also larger alphaS, more intrinsic kT
          IF (ITUNEB.EQ.101) THEN
            PARP(62) = 1.25D0
            PARP(64) = 0.2D0
            PARP(91) = 2.1D0
            PARP(92) = 15.0D0
          ENDIF
 
C...Tune BW (larger alphaS, more intrinsic kT. Smaller ISR phase space)
        ELSEIF (ITUNEB.EQ.102) THEN
C...pT0.
          PARP(82) = 1.9D0
c...String drawing completely minimizes string length.
          PARP(85) = 1.0D0
          PARP(86) = 1.0D0
C...ISR cutoff, muR scale factor, and phase space size
          PARP(62) = 1.25D0
          PARP(64) = 0.2D0
          PARP(67) = 1D0
C...Intrinsic kT, size, and max
          MSTP(91) = 1
          PARP(91) = 2.1D0
          PARP(93) = 15D0
 
C...Tune DW
        ELSEIF (ITUNEB.EQ.103) THEN
C...pT0.
          PARP(82) = 1.9D0
c...String drawing completely minimizes string length.
          PARP(85) = 1.0D0
          PARP(86) = 1.0D0
C...ISR cutoff, muR scale factor, and phase space size
          PARP(62) = 1.25D0
          PARP(64) = 0.2D0
          PARP(67) = 2.5D0
C...Intrinsic kT, size, and max
          MSTP(91) = 1
          PARP(91) = 2.1D0
          PARP(93) = 15D0
 
C...Tune DWT
        ELSEIF (ITUNEB.EQ.104) THEN
C...pT0.
          PARP(82) = 1.9409D0
C...Run II ref scale and slow scaling
          PARP(89) = 1960D0
          PARP(90) = 0.16D0
c...String drawing completely minimizes string length.
          PARP(85) = 1.0D0
          PARP(86) = 1.0D0
C...ISR cutoff, muR scale factor, and phase space size
          PARP(62) = 1.25D0
          PARP(64) = 0.2D0
          PARP(67) = 2.5D0
C...Intrinsic kT, size, and max
          MSTP(91) = 1
          PARP(91) = 2.1D0
          PARP(93) = 15D0
 
C...Tune QW
        ELSEIF(ITUNEB.EQ.105) THEN
          IF (M13.GE.1) THEN
            WRITE(M11,5030) ' '
            CH70='NB! This tune requires CTEQ6.1 pdfs to be '//
     &           'externally linked'
            WRITE(M11,5035) CH70
          ENDIF
C...pT0.
          PARP(82) = 1.1D0
c...String drawing completely minimizes string length.
          PARP(85) = 1.0D0
          PARP(86) = 1.0D0
C...ISR cutoff, muR scale factor, and phase space size
          PARP(62) = 1.25D0
          PARP(64) = 0.2D0
          PARP(67) = 2.5D0
C...Intrinsic kT, size, and max
          MSTP(91) = 1
          PARP(91) = 2.1D0
          PARP(93) = 15D0
 
C...Tune D6 and D6T
        ELSEIF(ITUNEB.EQ.108.OR.ITUNEB.EQ.109) THEN
          IF (M13.GE.1) THEN
            WRITE(M11,5030) ' '
            CH70='NB! This tune requires CTEQ6L pdfs to be '//
     &           'externally linked'
            WRITE(M11,5035) CH70
          ENDIF
C...The "Rick" proton, double gauss with 0.5/0.4
          MSTP(82) = 4
          PARP(83) = 0.5D0
          PARP(84) = 0.4D0
c...String drawing completely minimizes string length.
          PARP(85) = 1.0D0
          PARP(86) = 1.0D0
          IF (ITUNEB.EQ.108) THEN
C...D6: pT0, Run I ref scale, and fast energy scaling
            PARP(82) = 1.8D0
            PARP(89) = 1800D0
            PARP(90) = 0.25D0
          ELSE
C...D6T: pT0, Run II ref scale, and slow energy scaling
            PARP(82) = 1.8387D0
            PARP(89) = 1960D0
            PARP(90) = 0.16D0
          ENDIF
C...ISR cutoff, muR scale factor, and phase space size
          PARP(62) = 1.25D0
          PARP(64) = 0.2D0
          PARP(67) = 2.5D0
C...Intrinsic kT, size, and max
          MSTP(91) = 1
          PARP(91) = 2.1D0
          PARP(93) = 15D0
 
C...Old ATLAS-DC2 5-parameter tune
        ELSEIF(ITUNEB.EQ.106) THEN
          IF (M13.GE.1) THEN
            WRITE(M11,5010) ITUNE, CHNAME
            CH60='see A. Moraes et al., SN-ATLAS-2006-057,'
            WRITE(M11,5030) CH60
            CH60='    R. Field in hep-ph/0610012,'
            WRITE(M11,5030) CH60
            CH60='and T. Sjostrand & M. v. Zijl, PRD36(1987)2019'
            WRITE(M11,5030) CH60
          ENDIF
C...  pT0.
          PARP(82) = 1.8D0
C...  Different ref and rescaling pacee
          PARP(89) = 1000D0
          PARP(90) = 0.16D0
C...  Parameters of mass distribution
          PARP(83) = 0.5D0
          PARP(84) = 0.5D0
C...  Old default string drawing
          PARP(85) = 0.33D0
          PARP(86) = 0.66D0
C...  ISR, phase space equivalent to Tune B
          PARP(62) = 1D0
          PARP(64) = 1D0
          PARP(67) = 1D0
C...  FSR
          PARP(71) = 4D0
C...  Intrinsic kT
          MSTP(91) = 1
          PARP(91) = 1D0
          PARP(93) = 5D0
 
C...Professor's Pro-Q2O Tune
        ELSEIF(ITUNE.EQ.129) THEN
          PARP(62) = 2.9
          PARP(64) = 0.14
          PARP(67) = 2.65
          PARP(82) = 1.9
          PARP(83) = 0.83
          PARP(84) = 0.6
          PARP(85) = 0.86
          PARP(86) = 0.93
          PARP(89) = 1800D0
          PARP(90) = 0.22
          MSTP(91) = 1
          PARP(91) = 2.1
          PARP(93) = 5.0
 
        ENDIF
 
C...  Output
        IF (M13.GE.1) THEN
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5050) 67, PARP(67), CHPARP(67)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          WRITE(M11,5060) 82, PARJ(82), CHPARJ(82)
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5050) 85, PARP(85), CHPARP(85)
          WRITE(M11,5050) 86, PARP(86), CHPARP(86)
          WRITE(M11,5040) 91, MSTP(91), CHMSTP(91)
          WRITE(M11,5050) 91, PARP(91), CHPARP(91)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)

        ENDIF
 
C=======================================================================
C... ACR, tune A with new CR (107)
      ELSEIF(ITUNE.EQ.107.OR.ITUNE.EQ.117) THEN
        IF (M13.GE.1) THEN
          WRITE(M11,5010) ITUNE, CHNAME
          CH60='Tune A modified with new colour reconnections'
          WRITE(M11,5030) CH60
          CH60='PARP(85)=0D0 and amount of CR is regulated by PARP(78)'
          WRITE(M11,5030) CH60
          CH60='see P. Skands & D. Wicke, hep-ph/0703081,'
          WRITE(M11,5030) CH60
          CH60='    R. Field, in hep-ph/0610012 (Tune A),'
          WRITE(M11,5030) CH60
          CH60='and T. Sjostrand & M. v. Zijl, PRD36(1987)2019'
          WRITE(M11,5030) CH60
          IF (ITUNE.EQ.117) THEN
            CH60='LEP parameters tuned by Professor, hep-ph/0907.2973'
            WRITE(M11,5030) CH60
          ENDIF
        ENDIF
        IF (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.406))THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune. Using defaults.')
          GOTO 100
        ENDIF
 
C...Make sure we start from old default fragmentation parameters
        PARJ(81) = 0.29
        PARJ(82) = 1.0
 
C...Use Professor's LEP pars if ITUNE >= 110
C...(i.e., for A-Pro, DW-Pro etc)
        IF (ITUNE.LT.110) THEN
C...# Old defaults
          MSTJ(11) = 4
C...# Old default flavour parameters
          PARJ(21) = 0.36
          PARJ(41) = 0.30
          PARJ(42) = 0.58
          PARJ(46) = 1.0
          PARJ(82) = 1.0
        ELSE
C...# Tuned flavour parameters:
          PARJ(1)  = 0.073
          PARJ(2)  = 0.2
          PARJ(3)  = 0.94
          PARJ(4)  = 0.032
          PARJ(11) = 0.31
          PARJ(12) = 0.4
          PARJ(13) = 0.54
          PARJ(25) = 0.63
          PARJ(26) = 0.12
C...# Switch on Bowler:
          MSTJ(11) = 5
C...# Fragmentation
          PARJ(21) = 0.325
          PARJ(41) = 0.5
          PARJ(42) = 0.6
          PARJ(47) = 0.67
          PARJ(81) = 0.29
          PARJ(82) = 1.65
        ENDIF
 
        MSTP(81) = 1
        PARP(89) = 1800D0
        PARP(90) = 0.25D0
        MSTP(82) = 4
        PARP(83) = 0.5D0
        PARP(84) = 0.4D0
        MSTP(51) = 7
        MSTP(52) = 1
        PARP(71) = 4D0
        PARP(82) = 2.0D0
        PARP(85) = 0.0D0
        PARP(86) = 0.66D0
        PARP(62) = 1D0
        PARP(64) = 1D0
        PARP(67) = 4D0
        MSTP(91) = 1
        PARP(91) = 1D0
        PARP(93) = 5D0
        MSTP(95) = 6
C...P78 changed from 0.12 to 0.09 in 6.4.19 to improve <pT>(Nch)
        PARP(78) = 0.09D0
C...Frag functions (only if not using Professor)
        IF (ITUNE.LE.109) THEN
          MSTJ(11) = 4
          PARJ(54) = -0.05
          PARJ(55) = -0.005
        ENDIF
 
C...Output
        IF (M13.GE.1) THEN
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5050) 67, PARP(67), CHPARP(67)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          WRITE(M11,5060) 82, PARJ(82), CHPARJ(82)
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5050) 85, PARP(85), CHPARP(85)
          WRITE(M11,5050) 86, PARP(86), CHPARP(86)
          WRITE(M11,5040) 91, MSTP(91), CHMSTP(91)
          WRITE(M11,5050) 91, PARP(91), CHPARP(91)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          WRITE(M11,5050) 78, PARP(78), CHPARP(78)

        ENDIF
 
C=======================================================================
C...Intermediate model. Rap tune
C...(retuned to post-6.406 IR factorization)
      ELSEIF(ITUNE.EQ.200) THEN
        IF (M13.GE.1) THEN
          WRITE(M11,5010) ITUNE, CHNAME
          CH60='see T. Sjostrand & P. Skands, JHEP03(2004)053'
          WRITE(M11,5030) CH60
        ENDIF
        IF (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.405))THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ENDIF
C...PDF
        MSTP(51) = 7
        MSTP(52) = 1
C...ISR
        PARP(62) = 1D0
        PARP(64) = 1D0
        PARP(67) = 4D0
C...FSR
        PARP(71) = 4D0
        PARJ(81) = 0.29D0
C...UE
        MSTP(81) = 11
        PARP(82) = 2.25D0
        PARP(89) = 1800D0
        PARP(90) = 0.25D0
C...  ExpOfPow(1.8) overlap profile
        MSTP(82) = 5
        PARP(83) = 1.8D0
C...  Valence qq
        MSTP(88) = 0
C...  Rap Tune
        MSTP(89) = 1
C...  Default diquark, BR-g-BR supp
        PARP(79) = 2D0
        PARP(80) = 0.01D0
C...  Final state reconnect.
        MSTP(95) = 1
        PARP(78) = 0.55D0
C...Fragmentation functions and c and b parameters
        MSTJ(11) = 4
        PARJ(54) = -0.05
        PARJ(55) = -0.005
C...  Output
        IF (M13.GE.1) THEN
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5050) 67, PARP(67), CHPARP(67)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5040) 88, MSTP(88), CHMSTP(88)
          WRITE(M11,5040) 89, MSTP(89), CHMSTP(89)
          WRITE(M11,5050) 79, PARP(79), CHPARP(79)
          WRITE(M11,5050) 80, PARP(80), CHPARP(80)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          WRITE(M11,5050) 78, PARP(78), CHPARP(78)

        ENDIF
 
C...APT(201), APT-Pro (211), Perugia-APT (221), Perugia-APT6 (226).
C...Old model for ISR and UE, new pT-ordered model for FSR
      ELSEIF(ITUNE.EQ.201.OR.ITUNE.EQ.211.OR.ITUNE.EQ.221.OR
     &       .ITUNE.EQ.226) THEN
        IF (M13.GE.1) THEN
          WRITE(M11,5010) ITUNE, CHNAME
          CH60='see P. Skands & D. Wicke, hep-ph/0703081 (Tune APT),'
          WRITE(M11,5030) CH60
          CH60='    R.D. Field, in hep-ph/0610012 (Tune A)'
          WRITE(M11,5030) CH60
          CH60='    T. Sjostrand & M. v. Zijl, PRD36(1987)2019'
          WRITE(M11,5030) CH60
          CH60='and T. Sjostrand & P. Skands, hep-ph/0408302'
          WRITE(M11,5030) CH60
          IF (ITUNE.EQ.211.OR.ITUNE.GE.221) THEN
            CH60='LEP parameters tuned by Professor, hep-ph/0907.2973'
            WRITE(M11,5030) CH60
          ENDIF
        ENDIF
        IF (MSTP(181).LE.5.OR.(MSTP(181).EQ.6.AND.MSTP(182).LE.411))THEN
          CALL PYERRM(9,'(PYTUNE:) linked PYTHIA version incompatible'//
     &        ' with tune.')
        ENDIF
C...First set as if Pythia tune A
C...Multiple interactions on, old framework
        MSTP(81) = 1
C...Fast IR cutoff energy scaling by default
        PARP(89) = 1800D0
        PARP(90) = 0.25D0
C...Default CTEQ5L (internal)
        MSTP(51) = 7
        MSTP(52) = 1
C...Double Gaussian matter distribution.
        MSTP(82) = 4
        PARP(83) = 0.5D0
        PARP(84) = 0.4D0
C...FSR activity.
        PARP(71) = 4D0
c...String drawing almost completely minimizes string length.
        PARP(85) = 0.9D0
        PARP(86) = 0.95D0
C...ISR cutoff, muR scale factor, and phase space size
        PARP(62) = 1D0
        PARP(64) = 1D0
        PARP(67) = 4D0
C...Intrinsic kT, size, and max
        MSTP(91) = 1
        PARP(91) = 1D0
        PARP(93) = 5D0
C...Use 2 GeV of primordial kT for "Perugia" version
        IF (ITUNE.EQ.221) THEN
          PARP(91) = 2D0
          PARP(93) = 10D0
        ENDIF
C...Use pT-ordered FSR
        MSTJ(41) = 12
C...Lambda_FSR scale for pT-ordering
        PARJ(81) = 0.23D0
C...Retune pT0 (changed from 2.1 to 2.05 in 6.4.20)
        PARP(82) = 2.05D0
C...Fragmentation functions and c and b parameters
C...(overwritten for 211, i.e., if using Professor pars)
        PARJ(54) = -0.05
        PARJ(55) = -0.005
 
C...Use Professor's LEP pars if ITUNE == 211, 221, 226
        IF (ITUNE.LT.210) THEN
C...# Old defaults
          MSTJ(11) = 4
C...# Old default flavour parameters
          PARJ(21) = 0.36
          PARJ(41) = 0.30
          PARJ(42) = 0.58
          PARJ(46) = 1.0
          PARJ(82) = 1.0
        ELSE
C...# Tuned flavour parameters:
          PARJ(1)  = 0.073
          PARJ(2)  = 0.2
          PARJ(3)  = 0.94
          PARJ(4)  = 0.032
          PARJ(11) = 0.31
          PARJ(12) = 0.4
          PARJ(13) = 0.54
          PARJ(25) = 0.63
          PARJ(26) = 0.12
C...# Always use pT-ordered shower:
          MSTJ(41) = 12
C...# Switch on Bowler:
          MSTJ(11) = 5
C...# Fragmentation
          PARJ(21) = 3.1327e-01
          PARJ(41) = 4.8989e-01
          PARJ(42) = 1.2018e+00
          PARJ(47) = 1.0000e+00
          PARJ(81) = 2.5696e-01
          PARJ(82) = 8.0000e-01
        ENDIF
 
C...221, 226 : Perugia-APT and Perugia-APT6
        IF (ITUNE.EQ.221.OR.ITUNE.EQ.226) THEN
 
          PARP(64) = 0.5D0
          PARP(82) = 2.05D0
          PARP(90) = 0.26D0
          PARP(91) = 2.0D0
C...The Perugia variants use Steve's showers off the old MPI
          MSTP(152) = 1
C...And use a lower PARP(71) as suggested by Professor tunings
C...(although not certain that applies to Q2-pT2 hybrid)
          PARP(71) = 2.5D0
 
C...Perugia-APT6 uses CTEQ6L1 and a slightly lower pT0
          IF (ITUNE.EQ.226) THEN
            CH70='NB! This tune requires CTEQ6L1 pdfs to be '//
     &           'externally linked'
            WRITE(M11,5035) CH70
            MSTP(52) = 2
            MSTP(51) = 10042
            PARP(82) = 1.95D0
          ENDIF
 
        ENDIF
 
C...  Output
        IF (M13.GE.1) THEN
          WRITE(M11,5030) ' '
          WRITE(M11,5040) 51, MSTP(51), CHMSTP(51)
          WRITE(M11,5040) 52, MSTP(52), CHMSTP(52)
          WRITE(M11,5040)  3, MSTP( 3), CHMSTP( 3)
          WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          WRITE(M11,5050) 64, PARP(64), CHPARP(64)
          WRITE(M11,5050) 67, PARP(67), CHPARP(67)
          WRITE(M11,5040) 68, MSTP(68), CHMSTP(68)
          CH60='(Note: MSTP(68) is not explicitly (re-)set by PYTUNE)'
          WRITE(M11,5030) CH60
          WRITE(M11,5070) 41, MSTJ(41), CHMSTJ(41)
          WRITE(M11,5050) 71, PARP(71), CHPARP(71)
          WRITE(M11,5060) 81, PARJ(81), CHPARJ(81)
          WRITE(M11,5040) 33, MSTP(33), CHMSTP(33)
          WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
          WRITE(M11,5050) 82, PARP(82), CHPARP(82)
          WRITE(M11,5050) 89, PARP(89), CHPARP(89)
          WRITE(M11,5050) 90, PARP(90), CHPARP(90)
          WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          WRITE(M11,5050) 83, PARP(83), CHPARP(83)
          WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          IF (MSTP(82).GE.2) THEN
            WRITE(M11,5050) 87, PARP(87), CHPARP(87)
            IF (PARP(87).GE.0D0) 
     &           WRITE(M11,5050) 88, PARP(88), CHPARP(88)            
          ENDIF
          WRITE(M11,5050) 85, PARP(85), CHPARP(85)
          WRITE(M11,5050) 86, PARP(86), CHPARP(86)
          WRITE(M11,5040) 91, MSTP(91), CHMSTP(91)
          WRITE(M11,5050) 91, PARP(91), CHPARP(91)
          WRITE(M11,5050) 93, PARP(93), CHPARP(93)

        ENDIF
 
C======================================================================
C...Uppsala models: Generalized Area Law and Soft Colour Interactions
      ELSEIF(CHNAME.EQ.'GAL Tune 0'.OR.CHNAME.EQ.'GAL Tune 1') THEN
        IF (M13.GE.1) THEN
          WRITE(M11,5010) ITUNE, CHNAME
          CH60='see J. Rathsman, PLB452(1999)364'
          WRITE(M11,5030) CH60
          CH60='and T. Sjostrand & M. v. Zijl, PRD36(1987)2019'
          WRITE(M11,5030) CH60
        ENDIF
C...GAL Recommended settings from Uppsala web page 
        MSTP(95) = 13
        PARP(78) = 0.10
        MSTJ(16) = 0
        PARJ(42) = 0.45
        PARJ(82) = 2.0
        PARP(62) = 2.0
        MSTP(81) = 1
        MSTP(82) = 1
        PARP(81) = 1.9
        MSTP(92) = 1
        IF(CHNAME.EQ.'GAL Tune 1') THEN
C...GAL retune (P. Skands) to get better min-bias <Nch> at Tevatron
          MSTP(82) = 4
          PARP(83) = 0.25D0
          PARP(84) = 0.5D0
          PARP(82) = 1.75
          IF (M13.GE.1) THEN
            WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
            WRITE(M11,5050) 82, PARP(82), CHPARP(82)
            WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
            WRITE(M11,5050) 83, PARP(83), CHPARP(83)
            WRITE(M11,5050) 84, PARP(84), CHPARP(84)
          ENDIF
        ELSE
          IF (M13.GE.1) THEN
            WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
            WRITE(M11,5050) 81, PARP(81), CHPARP(81)
            WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          ENDIF
        ENDIF
C...Output
        IF (M13.GE.1) THEN
          WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          WRITE(M11,5060) 82, PARJ(82), CHPARJ(82)
          WRITE(M11,5040) 92, MSTP(92), CHMSTP(92)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          WRITE(M11,5050) 78, PARP(78), CHPARP(78)
          WRITE(M11,5060) 42, PARJ(42), CHPARJ(42)
          WRITE(M11,5070) 16, MSTJ(16), CHMSTJ(16)
        ENDIF
      ELSEIF(CHNAME.EQ.'SCI Tune 0'.OR.CHNAME.EQ.'SCI Tune 1') THEN
        IF (M13.GE.1) THEN
          WRITE(M11,5010) ITUNE, CHNAME
          CH60='see A.Edin et al, PLB366(1996)371, Z.Phys.C75(1997)57,'
          WRITE(M11,5030) CH60
          CH60='and T. Sjostrand & M. v. Zijl, PRD36(1987)2019'
          WRITE(M11,5030) CH60
          WRITE(M11,5030) ' '
          CH70='NB! The SCI model must be run with modified '//
     &        'Pythia v6.215:'
          WRITE(M11,5035) CH70
          CH70='available from http://www.isv.uu.se/thep/MC/scigal/'
          WRITE(M11,5035) CH70
          WRITE(M11,5030) ' '
        ENDIF
C...SCI Recommended settings from Uppsala web page (as per 22/08 2006)
        MSTP(81) = 1
        MSTP(82) = 1
        PARP(81) = 2.2
        MSTP(92) = 1
        MSTP(95) = 11
        PARP(78) = 0.50
        MSTJ(16) = 0
        IF (CHNAME.EQ.'SCI Tune 1') THEN
C...SCI retune (P. Skands) to get better min-bias <Nch> at Tevatron
          MSTP(81) = 1
          MSTP(82) = 3
          PARP(82) = 2.4
          PARP(83) = 0.5D0
          PARP(62) = 1.5
          PARP(84) = 0.25D0
          IF (M13.GE.1) THEN
            WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
            WRITE(M11,5050) 82, PARP(82), CHPARP(82)
            WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
            WRITE(M11,5050) 83, PARP(83), CHPARP(83)
            WRITE(M11,5050) 62, PARP(62), CHPARP(62)
          ENDIF
        ELSE
          IF (M13.GE.1) THEN
            WRITE(M11,5040) 81, MSTP(81), CHMSTP(81)
            WRITE(M11,5050) 81, PARP(81), CHPARP(81)
            WRITE(M11,5040) 82, MSTP(82), CHMSTP(82)
          ENDIF
        ENDIF
C...Output
        IF (M13.GE.1) THEN
          WRITE(M11,5040) 92, MSTP(92), CHMSTP(92)
          WRITE(M11,5040) 95, MSTP(95), CHMSTP(95)
          WRITE(M11,5050) 78, PARP(78), CHPARP(78)
          WRITE(M11,5070) 16, MSTJ(16), CHMSTJ(16)
        ENDIF
 
      ELSE
        IF (MSTU(13).GE.1) WRITE(M11,5020) ITUNE
 
      ENDIF
 
C...Output of LEP parameters, common to all models
      IF (M13.GE.1) THEN
        WRITE(M11,5080) 
        WRITE(M11,5070) 11, MSTJ(11), CHMSTJ(11)
        IF (MSTJ(11).EQ.3) THEN
          CH60='Warning: using Peterson fragmentation function'
          WRITE(M11,5030) CH60 
        ENDIF
        
        WRITE(M11,5060)  1, PARJ( 1), CHPARJ( 1)
        WRITE(M11,5060)  2, PARJ( 2), CHPARJ( 2)
        WRITE(M11,5060)  3, PARJ( 3), CHPARJ( 3)
        WRITE(M11,5060)  4, PARJ( 4), CHPARJ( 4)
        WRITE(M11,5060)  5, PARJ( 5), CHPARJ( 5)
        WRITE(M11,5060)  6, PARJ( 6), CHPARJ( 6)
        WRITE(M11,5060)  7, PARJ( 7), CHPARJ( 7)
        
        WRITE(M11,5060) 11, PARJ(11), CHPARJ(11)
        WRITE(M11,5060) 12, PARJ(12), CHPARJ(12)
        WRITE(M11,5060) 13, PARJ(13), CHPARJ(13)
        
        WRITE(M11,5060) 14, PARJ(14), CHPARJ(14)
        WRITE(M11,5060) 15, PARJ(15), CHPARJ(15)
        WRITE(M11,5060) 16, PARJ(16), CHPARJ(16)
        WRITE(M11,5060) 17, PARJ(17), CHPARJ(17)
        WRITE(M11,5060) 18, PARJ(18), CHPARJ(18)
        WRITE(M11,5060) 19, PARJ(19), CHPARJ(19)
        
        WRITE(M11,5060) 21, PARJ(21), CHPARJ(21)
        
        WRITE(M11,5060) 25, PARJ(25), CHPARJ(25)
        WRITE(M11,5060) 26, PARJ(26), CHPARJ(26)
        
        WRITE(M11,5060) 41, PARJ(41), CHPARJ(41)
        WRITE(M11,5060) 42, PARJ(42), CHPARJ(42)
        WRITE(M11,5060) 45, PARJ(45), CHPARJ(45)
        
        IF (MSTJ(11).LE.3) THEN
          WRITE(M11,5060) 54, PARJ(54), CHPARJ(54)
          WRITE(M11,5060) 55, PARJ(55), CHPARJ(55)
        ELSE
          WRITE(M11,5060) 46, PARJ(46), CHPARJ(46)
        ENDIF
        IF (MSTJ(11).EQ.5) WRITE(M11,5060) 47, PARJ(47), CHPARJ(47)
      ENDIF
        
 100  IF (MSTU(13).GE.1) WRITE(M11,6000)
 
 9999 RETURN
 
 5000 FORMAT(1x,78('*')/' *',76x,'*'/' *',3x,'PYTUNE : ',
     &    'Presets for underlying-event (and min-bias)',21x,'*'/' *',
     &    12x,'Last Change : ',A8,' - P. Skands',30x,'*'/' *',76x,'*')
 5010 FORMAT(' *',3x,I4,1x,A16,52x,'*')
 5020 FORMAT(' *',3x,'Tune ',I4, ' not recognized. Using defaults.')
 5030 FORMAT(' *',3x,10x,A60,3x,'*')
 5035 FORMAT(' *',3x,A70,3x,'*')
 5040 FORMAT(' *',5x,'MSTP(',I2,') = ',I12,3x,A42,3x,'*')
 5050 FORMAT(' *',5x,'PARP(',I2,') = ',F12.4,3x,A40,5x,'*')
 5060 FORMAT(' *',5x,'PARJ(',I2,') = ',F12.4,3x,A40,5x,'*')
 5070 FORMAT(' *',5x,'MSTJ(',I2,') = ',I12,3x,A40,5x,'*')
 5080 FORMAT(' *',3x,'----------------------------',42('-'),3x,'*')
 6100 FORMAT(' *',5x,'MSTU(',I3,')= ',I12,3x,A42,3x,'*')
 6110 FORMAT(' *',5x,'PARU(',I3,')= ',F12.4,3x,A42,3x,'*')
C 5140 FORMAT(' *',5x,'MSTP(',I3,')= ',I12,3x,A40,5x,'*')
C 5150 FORMAT(' *',5x,'PARP(',I3,')= ',F12.4,3x,A40,5x,'*')
 6000 FORMAT(' *',76x,'*'/1x,32('*'),1x,'END OF PYTUNE',1x,31('*'))
C 6040 FORMAT(' *',5x,'MSWI(',I1,')  = ',I12,3x,A40,5x,'*')
C 6050 FORMAT(' *',5x,'PARSCI(',I1,')= ',F12.4,3x,A40,5x,'*')
 
      END
