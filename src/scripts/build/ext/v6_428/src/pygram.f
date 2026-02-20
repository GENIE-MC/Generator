C*********************************************************************
 
C...PYGRAM
C...Universal Extra Dimensions Model (UED)
C...Computation of the Graviton mass.

      SUBROUTINE PYGRAM(IN)

C...Double precision and integer declarations
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)

C...Pythia commonblocks
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)      
C...UED Pythia common
      COMMON/PYPUED/IUED(0:99),RUED(0:99)

C...Local variables
      INTEGER KCFLA,NMAX
      PARAMETER(KCFLA=450,NMAX=5000)
      DIMENSION YVEC(5000),RESVEC(5000)
      COMMON/INTSAV/YSAV,YMAX,RESMAX
      COMMON/UEDGRA/XMPLNK,XMD,RINV,NDIM
      COMMON/KAPPA/XKAPPA

C...External function (used in call to PYGAUS)
      EXTERNAL PYGRAW

C...SAVE statements
      SAVE /PYDAT1/,/PYDAT2/,/PYPUED/,/INTSAV/

C...Initialization
      NDIM=IUED(4)
      RINV=RUED(1)
      XMD=RUED(2)
      PI=PARU(1)

C...Initialize for numerical integration
      XMPLNK=2.4D+18
      XKAPPA=DSQRT(2.D0)/XMPLNK      

C...For NDIM=2, compute graviton mass distribution numerically
      IF(NDIM.EQ.2)THEN
        
C...  For first event: tabulate distribution of stepwise integrals:
C...  int_y1^y2 dy dGamma/dy , with y = MG*/MgammaKK
        IF(IN.EQ.0)THEN
          RESMAX = 0D0
          YMAX   = 0D0
          DO 100 I=1,NMAX
            YSAV = (I-0.5)/DBLE(NMAX)
            TOL       = 1D-6
C...Integral of PYGRAW from 0 to 1, with precision TOL, for given YSAV
            RESINT    = PYGAUS(PYGRAW,0D0,1D0,TOL)
            YVEC(I)   = YSAV
            RESVEC(I) = RESINT
C...  Save max of distribution (for accept/reject below)
            IF(RESINT.GT.RESMAX)THEN
              RESMAX = RESINT
              YMAX   = YVEC(I)
            ENDIF
 100      CONTINUE
        ENDIF
        
C...  Generate Mg for each graviton (1D0 ensures a minimal open phase space)
        PCUJET=1D0
        KCGAKK=KCFLA+23
        XMGAMK=PMAS(KCGAKK,1)
        
C...  Pick random graviton mass, accept according to stored integrals
        AMMAX=DSQRT(XMGAMK**2-2D0*XMGAMK*PCUJET)
 110    RMG=AMMAX*PYR(0)
        X=RMG/XMGAMK        

C...  Bin enumeration starts at 1, but make sure always in range
        IBIN=INT(NMAX*X)+1
        IBIN=MIN(IBIN,NMAX)        
        IF(RESVEC(IBIN)/RESMAX.LT.PYR(0)) GOTO 110
        
C...  For NDIM=4 and 6, the analytical expression for the
C...  graviton mass distribution integral is used.
      ELSEIF(NDIM.EQ.4.OR.NDIM.EQ.6)THEN
        
C...  Ensure minimal open phase space (max(mG*) < m(gamma*))
        PCUJET=1D0
        
C...  KK photon (?) compressed code and mass
        KCGAKK=KCFLA+23
        XMGAMK=PMAS(KCGAKK,1)
        
C...  Find maximum of (dGamma/dMg)
        IF(IN.EQ.0)THEN
          RESMAX=0D0
          YMAX=0D0
          DO 120 I=1,NMAX-1 
            Y=I/DBLE(NMAX)
            RESINT=Y**(NDIM-3)*(1D0/(1D0-Y**2))*(1D0+DCOS(PI*Y))
            IF(RESINT.GE.RESMAX)THEN
              RESMAX=RESINT
              YMAX=Y
            ENDIF
 120      CONTINUE
        ENDIF
        
C...  Pick random graviton mass, accept/reject
        AMMAX=DSQRT(XMGAMK**2-2D0*XMGAMK*PCUJET)
 130    RMG=AMMAX*PYR(0)
        X=RMG/XMGAMK
        DGADMG=X**(NDIM-3)*(1./(1.-X**2))*(1.+DCOS(PI*X))
        IF(DGADMG/RESMAX.LT.PYR(0)) GOTO 130
        
C...  If the user has not chosen N=2,4 or 6, STOP
      ELSE
        WRITE(MSTU(11),*) '(PYGRAM:) BAD VALUE N(LARGE XD) =',NDIM,
     &       ' (MUST BE 2, 4, OR 6) '
        CALL PYSTOP(6002)
      ENDIF
      
C...  Now store the sampled Mg
      PMAS(39,1)=RMG
      
      RETURN
      END
