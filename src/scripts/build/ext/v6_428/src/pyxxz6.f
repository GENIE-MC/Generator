 
C*********************************************************************
 
C...PYXXZ6
C...Used in the calculation of  inoi -> inoj + f + ~f.
 
      FUNCTION PYXXZ6(X)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      COMMON/PYINTS/XXM(20)
      COMPLEX*16 CXC
      COMMON/PYINTC/XXC(10),CXC(8)
      SAVE /PYDAT1/,/PYINTC/
 
C...Local variables.
      COMPLEX*16 QLLS,QRRS,QRLS,QLRS,QLLU,QRRU,QLRT,QRLT
      DOUBLE PRECISION PYXXZ6,X
      DOUBLE PRECISION XM12,XM22,XM32,S,S13,WPROP2
      DOUBLE PRECISION WW,WF1,WF2,WFL1,WFL2
      DOUBLE PRECISION SIJ
      DOUBLE PRECISION XMV,XMG,XMSU1,XMSU2,XMSD1,XMSD2
      DOUBLE PRECISION OL2
      DOUBLE PRECISION S23MIN,S23MAX,S23AVE,S23DEL
      INTEGER I
 
C...Statement functions.
C...Integral from x to y of (t-a)(b-t) dt.
      TINT(X,Y,A,B)=(X-Y)*(-(X**2+X*Y+Y**2)/3D0+(B+A)*(X+Y)/2D0-A*B)
C...Integral from x to y of (t-a)(b-t)/(t-c) dt.
      TINT2(X,Y,A,B,C)=(X-Y)*(-0.5D0*(X+Y)+(B+A-C))-
     &LOG(ABS((X-C)/(Y-C)))*(C-B)*(C-A)
C...Integral from x to y of (t-a)(b-t)/(t-c)**2 dt.
      TINT3(X,Y,A,B,C)=-(X-Y)+(C-A)*(C-B)*(Y-X)/(X-C)/(Y-C)+
     &(B+A-2D0*C)*LOG(ABS((X-C)/(Y-C)))
C...Integral from x to y of (t-a)/(b-t) dt.
      UTINT(X,Y,A,B)=LOG(ABS((X-A)/(B-X)*(B-Y)/(Y-A)))/(B-A)
C...Integral from x to y of 1/(t-a) dt.
      TPROP(X,Y,A)=LOG(ABS((X-A)/(Y-A)))
 
      XM12=XXC(1)**2
      XM22=XXC(2)**2
      XM32=XXC(3)**2
      S=XXC(4)**2
      S13=X
 
      S23AVE=XM22+XM32-0.5D0/X*(X+XM32-XM12)*(X+XM22-S)
      S23DEL=0.5D0/X*SQRT( ( (X-XM12-XM32)**2-4D0*XM12*XM32)*
     &( (X-XM22-S)**2  -4D0*XM22*S  ) )
 
      S23MIN=(S23AVE-S23DEL)
      S23MAX=(S23AVE+S23DEL)
 
      XMSD1=XXC(5)**2
      XMSD2=XXC(7)**2
      XMSU1=XXC(6)**2
      XMSU2=XXC(8)**2
 
      XMV=XXC(9)
      XMG=XXC(10)
      QLLS=CXC(1)
      QLLU=CXC(2)
      QLRS=CXC(3)
      QLRT=CXC(4)
      QRLS=CXC(5)
      QRLT=CXC(6)
      QRRS=CXC(7)
      QRRU=CXC(8)
      WPROP2=(S13-XMV**2)**2+(XMV*XMG)**2
      SIJ=2D0*XXC(2)*XXC(4)*S13
      IF(XMV.LE.1000D0) THEN
        OL2=ABS(QLLS)**2+ABS(QRRS)**2+ABS(QLRS)**2+ABS(QRLS)**2
        OLR=-2D0*DBLE(QLRS*DCONJG(QLLS)+QRLS*DCONJG(QRRS))
        WW=(OL2*2D0*TINT(S23MAX,S23MIN,XM22,S)
     &  +OLR*SIJ*(S23MAX-S23MIN))/WPROP2
        IF(XXC(5).LE.10000D0) THEN
          WFL1=4D0*(DBLE(QLLS*DCONJG(QLLU))*
     &    TINT2(S23MAX,S23MIN,XM22,S,XMSD1)-
     &    .5D0*DBLE(QLLS*DCONJG(QLRT))*SIJ*TPROP(S23MAX,S23MIN,XMSD2)+
     &    DBLE(QLRS*DCONJG(QLRT))*TINT2(S23MAX,S23MIN,XM22,S,XMSD2)-
     &    .5D0*DBLE(QLRS*DCONJG(QLLU))*SIJ*TPROP(S23MAX,S23MIN,XMSD1))
     &    *(S13-XMV**2)/WPROP2
        ELSE
          WFL1=0D0
        ENDIF
 
        IF(XXC(6).LE.10000D0) THEN
          WFL2=4D0*(DBLE(QRRS*DCONJG(QRRU))*
     &    TINT2(S23MAX,S23MIN,XM22,S,XMSU1)-
     &    .5D0*DBLE(QRRS*DCONJG(QRLT))*SIJ*TPROP(S23MAX,S23MIN,XMSU2)+
     &    DBLE(QRLS*DCONJG(QRLT))*TINT2(S23MAX,S23MIN,XM22,S,XMSU2)-
     &    .5D0*DBLE(QRLS*DCONJG(QRRU))*SIJ*TPROP(S23MAX,S23MIN,XMSU1))
     &    *(S13-XMV**2)/WPROP2
        ELSE
          WFL2=0D0
        ENDIF
      ELSE
        WW=0D0
        WFL1=0D0
        WFL2=0D0
      ENDIF
      IF(XXC(5).LE.10000D0) THEN
        WF1=2D0*ABS(QLLU)**2*TINT3(S23MAX,S23MIN,XM22,S,XMSD1)
     &  +2D0*ABS(QLRT)**2*TINT3(S23MAX,S23MIN,XM22,S,XMSD2)
     &  - 2D0*DBLE(QLRT*DCONJG(QLLU))*
     &  SIJ*UTINT(S23MAX,S23MIN,XMSD1,XM22+S-S13-XMSD2)
      ELSE
        WF1=0D0
      ENDIF
      IF(XXC(6).LE.10000D0) THEN
        WF2=2D0*ABS(QRRU)**2*TINT3(S23MAX,S23MIN,XM22,S,XMSU1)
     &  +2D0*ABS(QRLT)**2*TINT3(S23MAX,S23MIN,XM22,S,XMSU2)
     &  - 2D0*DBLE(QRLT*DCONJG(QRRU))*
     &  SIJ*UTINT(S23MAX,S23MIN,XMSU1,XM22+S-S13-XMSU2)
      ELSE
        WF2=0D0
      ENDIF
 
      PYXXZ6=(WW+WF1+WF2+WFL1+WFL2)
 
      IF(PYXXZ6.LT.0D0) THEN
        WRITE(MSTU(11),*) ' NEGATIVE WT IN PYXXZ6 '
        WRITE(MSTU(11),*) (XXC(I),I=1,5)
        WRITE(MSTU(11),*) (XXC(I),I=6,10)
        WRITE(MSTU(11),*) WW,WF1,WF2,WFL1,WFL2
        WRITE(MSTU(11),*) S23MIN,S23MAX
        PYXXZ6=0D0
      ENDIF
 
      RETURN
      END
