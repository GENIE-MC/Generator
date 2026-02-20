 
C*********************************************************************
 
C...PYSGTC
C...Subprocess cross sections for Technicolor processes.
C...Auxiliary to PYSIGH.
 
      SUBROUTINE PYSGTC(NCHN,SIGS)
 
C...Double precision and integer declarations
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYTCSM/ITCM(0:99),RTCM(0:99)
      COMMON/PYSGCM/ISUB,ISUBSV,MMIN1,MMAX1,MMIN2,MMAX2,MMINA,MMAXA,
     &KFAC(2,-40:40),COMFAC,FACK,FACA,SH,TH,UH,SH2,TH2,UH2,SQM3,SQM4,
     &SHR,SQPTH,TAUP,BE34,CTH,X(2),SQMZ,SQMW,GMMZ,GMMW,
     &AEM,AS,XW,XW1,XWC,XWV,POLL,POLR,POLLL,POLRR
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYPARS/,/PYINT1/,/PYINT2/,
     &/PYINT3/,/PYINT4/,/PYTCSM/,/PYSGCM/
C...Local arrays and complex variables
      DIMENSION WDTP(0:400),WDTE(0:400,0:5)
      COMPLEX*16 SSMZ,SSMR,SSMO,DETD,F2L,F2R,DARHO,DZRHO,DAOME,DZOME
      COMPLEX*16 SSMX,DAAST,DZAST,DWAST
      COMPLEX*16 DAA,DZZ,DAZ,DWW,DWRHO
      COMPLEX*16 ZTC(6,6),YTC(6,6),DGGS,DGGT,DGGU,DGVS,DGVT,DGVU
      COMPLEX*16 DQQS,DQQT,DQQU,DQTS,DQGS,DTGS
      COMPLEX*16 DVVS,DVVT,DVVU
      INTEGER INDX(6)
 
C...Combinations of weak mixing angle.
      TANW=SQRT(XW/XW1)
      CT2W=(1D0-2D0*XW)/(2D0*XW/TANW)
 
C...Convert almost equivalent technicolor processes into
C...a few basic processes, and set distinguishing parameters.
      IF(ISUB.GE.361.AND.ISUB.LE.380) THEN
        SQTV=RTCM(12)**2
        SQTA=RTCM(13)**2
        SN2W=2D0*SQRT(XW*XW1)
        CS2W=1D0-2D0*XW
        CT2W=CS2W/SN2W
        CSXI=COS(ASIN(RTCM(3)))
        CSXIP=COS(ASIN(RTCM(4)))
        QUPD=2D0*RTCM(2)-1D0
        Q2UD=RTCM(2)**2+(RTCM(2)-1D0)**2
        CAB2=0D0
        VOGP=0D0
        VRGP=0D0
        AOGP=0D0
        ARGP=0D0
        VXGP=0D0
        AXGP=0D0
        VAGP=0D0
        VZGP=0D0
        VWGP=0D0
C... rho_tc0, etc. -> W_L W_L, W_L W_T
        IF(ISUB.EQ.361) THEN
           KFA=24
           KFB=24
           CAB2=RTCM(3)**4
           AXGP=-RTCM(3)/(2D0*SQRT(XW))/RTCM(49)
           ARGP=RTCM(3)/(2D0*SQRT(XW))/RTCM(13)
           VOGP=RTCM(3)/(2D0*SQRT(XW))/RTCM(12)
C...Multiply by sqrt(2) to account for W^+_T W^-_L + W^+_L W^-_T.
           AXGP = SQRT(2D0)*AXGP
           ARGP = SQRT(2D0)*ARGP
           VOGP = SQRT(2D0)*VOGP
C... rho_tc0 -> W_L pi_tc-
        ELSEIF(ISUB.EQ.362) THEN
           KFA=24
           KFB=KTECHN+211
           ISUB=361
           CAB2=RTCM(3)**2*(1D0-RTCM(3)**2)
C... pi_tc pi_tc
        ELSEIF(ISUB.EQ.363) THEN
           KFA=KTECHN+211
           KFB=KTECHN+211
           ISUB=361
           CAB2=(1D0-RTCM(3)**2)**2
C... rho_tc0/omega_tc -> gamma pi_tc
        ELSEIF(ISUB.EQ.364) THEN
           KFA=22
           KFB=KTECHN+111
           ISUB=361
           VOGP=CSXI/RTCM(12)
           VRGP=VOGP*QUPD
           VAGP=2D0*QUPD*CSXI
           VZGP=QUPD*CSXI*(1D0-4D0*XW)/SN2W
C... gamma pi_tc'
        ELSEIF(ISUB.EQ.365) THEN
           KFA=22
           KFB=KTECHN+221
           ISUB=361
           VRGP=CSXIP/RTCM(12)
           VOGP=VRGP*QUPD
           VAGP=2D0*Q2UD*CSXIP
           VZGP=CSXIP/SN2W*(1D0-4D0*XW*Q2UD)
C... Z pi_tc
        ELSEIF(ISUB.EQ.366) THEN
           KFA=23
           KFB=KTECHN+111
           ISUB=361
           VOGP=CSXI*CT2W/RTCM(12)
           VRGP=-QUPD*CSXI*TANW/RTCM(12)
           VAGP=QUPD*CSXI*(1D0-4D0*XW)/SN2W
           VZGP=-QUPD*CSXI*CS2W/XW1
C... Z pi_tc'
        ELSEIF(ISUB.EQ.367) THEN
           KFA=23
           KFB=KTECHN+221
           ISUB=361
C...RTCM(48) is the M_V for the techni-a
           VXGP=-CSXIP/SN2W/RTCM(48)
           VRGP=CSXIP*CT2W/RTCM(12)
           VOGP=-QUPD*CSXIP*TANW/RTCM(12)
           VAGP=CSXIP*(1D0-4D0*Q2UD*XW)/SN2W
           VZGP=2D0*CSXIP*(CS2W+4D0*Q2UD*XW**2)/SN2W**2
C... W_T pi_tc
        ELSEIF(ISUB.EQ.368) THEN
           KFA=24
           KFB=KTECHN+211
           ISUB=361
C...RTCM(49) is the M_A for the techni-a
           AXGP=-CSXI/(2D0*SQRT(XW))/RTCM(49)
           VOGP=CSXI/(2D0*SQRT(XW))/RTCM(12)
           ARGP=CSXI/(2D0*SQRT(XW))/RTCM(13)
           VAGP=QUPD*CSXI/(2D0*SQRT(XW))
           VZGP=-QUPD*CSXI/(2D0*SQRT(XW1))
C... rho_tc+, a_T+ -> W_L Z_L, W_T Z_L
        ELSEIF(ISUB.EQ.370) THEN
           KFA=24
           KFB=23
           CAB2=RTCM(3)**4
           ARGP=-RTCM(3)/(2D0*SQRT(XW))/RTCM(13)
           AXGP=RTCM(3)/(2D0*SQRT(XW))/RTCM(49)
C... W_L pi_tc0
        ELSEIF(ISUB.EQ.371) THEN
           KFA=24
           KFB=KTECHN+111
           ISUB=370
           CAB2=RTCM(3)**2*(1D0-RTCM(3)**2)
C... Z_L pi_tc+
        ELSEIF(ISUB.EQ.372) THEN
           KFA=KTECHN+211
           KFB=23
           ISUB=370
           CAB2=RTCM(3)**2*(1D0-RTCM(3)**2)
C... pi_tc+ pi_tc0
        ELSEIF(ISUB.EQ.373) THEN
           KFA=KTECHN+211
           KFB=KTECHN+111
           ISUB=370
           CAB2=(1D0-RTCM(3)**2)**2
C... gamma pi_tc+
        ELSEIF(ISUB.EQ.374) THEN
           KFA=KTECHN+211
           KFB=22
           ISUB=370
           VRGP=QUPD*CSXI/RTCM(12)
           VWGP=QUPD*CSXI/(2D0*SQRT(XW))
           AXGP=-CSXI/RTCM(49)
C... Z_T pi_tc+
        ELSEIF(ISUB.EQ.375) THEN
           KFA=KTECHN+211
           KFB=23
           ISUB=370
           VRGP=-QUPD*CSXI*TANW/RTCM(12)
           ARGP=CSXI/(2D0*SQRT(XW*XW1))/RTCM(13)
           VWGP=-QUPD*CSXI/(2D0*SQRT(XW1))
           AXGP=-CSXI*CT2W/RTCM(49)
C... W_T pi_tc0
        ELSEIF(ISUB.EQ.376) THEN
           KFA=24
           KFB=KTECHN+111
           ISUB=370
           VRGP=0D0
           ARGP=-CSXI/(2D0*SQRT(XW))/RTCM(13)
           AXGP=CSXI/(2D0*SQRT(XW))/RTCM(49)
C... W_T pi_tc0'
        ELSEIF(ISUB.EQ.377) THEN
           KFA=24
           KFB=KTECHN+221
           ISUB=370
           VRGP=CSXIP/(2D0*SQRT(XW))/RTCM(12)
           VWGP=CSXIP/(2D0*XW)
           VXGP=-CSXIP/(2D0*SQRT(XW))/RTCM(48)
C... gamma W+
        ELSEIF(ISUB.EQ.378) THEN
           KFA=24
           KFB=22
           ISUB=370
           VRGP=QUPD*RTCM(3)/RTCM(12)
           AXGP=-RTCM(3)/RTCM(49)
C... gamma Z
        ELSEIF(ISUB.EQ.379) THEN
           KFA=23
           KFB=22
           ISUB=361
           VOGP=RTCM(3)/RTCM(12)
           VRGP=QUPD*RTCM(3)/RTCM(12)
        ELSEIF(ISUB.EQ.380) THEN
           KFA=23
           KFB=23
           ISUB=361
           VOGP=RTCM(3)*CT2W/RTCM(12)
           VRGP=-QUPD*RTCM(3)*TANW/RTCM(12)
        ENDIF
      ENDIF
 
C...QCD 2 -> 2 processes: corrections from virtual technicolor exchange.
      IF(ISUB.GE.381.AND.ISUB.LE.388) THEN
        IF(ITCM(5).LE.4) THEN
          SQDQQS=1D0/SH2
          SQDQQT=1D0/TH2
          SQDQQU=1D0/UH2
          SQDGGS=SQDQQS
          SQDGGT=SQDQQT
          SQDGGU=SQDQQU
          REDGGS=1D0/SH
          REDGGT=1D0/TH
          REDGGU=1D0/UH
          REDGTU=1D0/UH/TH
          REDGSU=1D0/SH/UH
          REDGST=1D0/SH/TH
          REDQST=1D0/SH/TH
          REDQTU=1D0/UH/TH
          SQDLGS=0D0
          SQDLGT=0D0
          SQDQTS=SQDQQS
        ELSEIF(ITCM(5).EQ.5) THEN
          TANT3=RTCM(21)
          IF(ITCM(2).EQ.0) THEN
            IMDL=1
          ELSE
            IMDL=2
          ENDIF
          ALPRHT=2.16D0*(3D0/ITCM(1))
          SIN2T=2D0*TANT3/(TANT3**2+1D0)
          SINT3=TANT3/SQRT(TANT3**2+1D0)
          XIG=SQRT(PYALPS(SH)/ALPRHT)
          X12=(RTCM(29)*SQRT(1D0-RTCM(29)**2)*COS(RTCM(30))+
     &    RTCM(31)*SQRT(1D0-RTCM(31)**2)*COS(RTCM(32)))/SQRT(2D0)/SIN2T
          X21=(RTCM(29)*SQRT(1D0-RTCM(29)**2)*SIN(RTCM(30))+
     &    RTCM(31)*SQRT(1D0-RTCM(31)**2)*SIN(RTCM(32)))/SQRT(2D0)/SIN2T
          X11=(.25D0*(RTCM(29)**2+RTCM(31)**2+2D0)-
     &    SINT3**2)*2D0/SIN2T
          X22=(.25D0*(2D0-RTCM(29)**2-RTCM(31)**2)-
     &    SINT3**2)*2D0/SIN2T
 
          SM1122=.5D0*(2D0-RTCM(29)**2-RTCM(31)**2)*RTCM(28)**2
          SM1112=X12*RTCM(28)**2*SIN2T
          SM1121=-X21*RTCM(28)**2*SIN2T
          SM2212=-SM1112
          SM2221=-SM1121
          SM1221=-.5D0*((1D0-RTCM(29)**2)*SIN(2D0*RTCM(30))+
     &    (1D0-RTCM(31)**2)*SIN(2D0*RTCM(32)))*RTCM(28)**2
 
C.........SH LOOP
          ZTC(1,1)=DCMPLX(SH,0D0)
          CALL PYWIDT(3100021,SH,WDTP,WDTE)
          IF(WDTP(0).GT.RTCM(33)*SHR) WDTP(0)=RTCM(33)*SHR
          ZTC(2,2)=DCMPLX(SH-PMAS(PYCOMP(3100021),1)**2,-SHR*WDTP(0))
          CALL PYWIDT(3100113,SH,WDTP,WDTE)
          ZTC(3,3)=DCMPLX(SH-PMAS(PYCOMP(3100113),1)**2,-SHR*WDTP(0))
          CALL PYWIDT(3400113,SH,WDTP,WDTE)
          ZTC(4,4)=DCMPLX(SH-PMAS(PYCOMP(3400113),1)**2,-SHR*WDTP(0))
          CALL PYWIDT(3200113,SH,WDTP,WDTE)
          ZTC(5,5)=DCMPLX(SH-PMAS(PYCOMP(3200113),1)**2,-SHR*WDTP(0))
          CALL PYWIDT(3300113,SH,WDTP,WDTE)
          ZTC(6,6)=DCMPLX(SH-PMAS(PYCOMP(3300113),1)**2,-SHR*WDTP(0))
          ZTC(1,2)=(0D0,0D0)
          ZTC(1,3)=DCMPLX(SH*XIG,0D0)
          ZTC(1,4)=ZTC(1,3)
          ZTC(1,5)=ZTC(1,2)
          ZTC(1,6)=ZTC(1,2)
          ZTC(2,3)=DCMPLX(SH*XIG*X11,0D0)
          ZTC(2,4)=DCMPLX(SH*XIG*X22,0D0)
          ZTC(2,5)=DCMPLX(SH*XIG*X12,0D0)
          ZTC(2,6)=DCMPLX(SH*XIG*X21,0D0)
          ZTC(3,4)=-SM1122
          ZTC(3,5)=-SM1112
          ZTC(3,6)=-SM1121
          ZTC(4,5)=-SM2212
          ZTC(4,6)=-SM2221
          ZTC(5,6)=-SM1221
 
          DO 110 I=1,5
            DO 100 J=I+1,6
               ZTC(J,I)=ZTC(I,J)
  100       CONTINUE
  110     CONTINUE
          CALL PYLDCM(ZTC,6,6,INDX,D)
          DO 130 I=1,6
            DO 120 J=1,6
             YTC(I,J)=(0D0,0D0)
              IF(I.EQ.J) YTC(I,J)=(1D0,0D0)
  120       CONTINUE
  130     CONTINUE
 
          DO 140 I=1,6
            CALL PYBKSB(ZTC,6,6,INDX,YTC(1,I))
  140     CONTINUE
          DGGS=YTC(1,1)
          DVVS=YTC(2,2)
          DGVS=YTC(1,2)
 
          XIG=SQRT(PYALPS(-TH)/ALPRHT)
C.........TH LOOP
          ZTC(1,1)=DCMPLX(TH)
          ZTC(2,2)=DCMPLX(TH-PMAS(PYCOMP(3100021),1)**2)
          ZTC(3,3)=DCMPLX(TH-PMAS(PYCOMP(3100113),1)**2)
          ZTC(4,4)=DCMPLX(TH-PMAS(PYCOMP(3400113),1)**2)
          ZTC(5,5)=DCMPLX(TH-PMAS(PYCOMP(3200113),1)**2)
          ZTC(6,6)=DCMPLX(TH-PMAS(PYCOMP(3300113),1)**2)
          ZTC(1,2)=(0D0,0D0)
          ZTC(1,3)=DCMPLX(TH*XIG,0D0)
          ZTC(1,4)=ZTC(1,3)
          ZTC(1,5)=ZTC(1,2)
          ZTC(1,6)=ZTC(1,2)
          ZTC(2,3)=DCMPLX(TH*XIG*X11,0D0)
          ZTC(2,4)=DCMPLX(TH*XIG*X22,0D0)
          ZTC(2,5)=DCMPLX(TH*XIG*X12,0D0)
          ZTC(2,6)=DCMPLX(TH*XIG*X21,0D0)
          ZTC(3,4)=-SM1122
          ZTC(3,5)=-SM1112
          ZTC(3,6)=-SM1121
          ZTC(4,5)=-SM2212
          ZTC(4,6)=-SM2221
          ZTC(5,6)=-SM1221
          DO 160 I=1,5
            DO 150 J=I+1,6
               ZTC(J,I)=ZTC(I,J)
  150       CONTINUE
  160     CONTINUE
          CALL PYLDCM(ZTC,6,6,INDX,D)
          DO 180 I=1,6
            DO 170 J=1,6
              YTC(I,J)=(0D0,0D0)
              IF(I.EQ.J) YTC(I,J)=(1D0,0D0)
  170       CONTINUE
  180     CONTINUE
          DO 190 I=1,6
            CALL PYBKSB(ZTC,6,6,INDX,YTC(1,I))
  190     CONTINUE
          DGGT=YTC(1,1)
          DVVT=YTC(2,2)
          DGVT=YTC(1,2)
 
          XIG=SQRT(PYALPS(-UH)/ALPRHT)
C.........UH LOOP
          ZTC(1,1)=DCMPLX(UH,0D0)
          ZTC(2,2)=DCMPLX(UH-PMAS(PYCOMP(3100021),1)**2)
          ZTC(3,3)=DCMPLX(UH-PMAS(PYCOMP(3100113),1)**2)
          ZTC(4,4)=DCMPLX(UH-PMAS(PYCOMP(3400113),1)**2)
          ZTC(5,5)=DCMPLX(UH-PMAS(PYCOMP(3200113),1)**2)
          ZTC(6,6)=DCMPLX(UH-PMAS(PYCOMP(3300113),1)**2)
          ZTC(1,2)=(0D0,0D0)
          ZTC(1,3)=DCMPLX(UH*XIG,0D0)
          ZTC(1,4)=ZTC(1,3)
          ZTC(1,5)=ZTC(1,2)
          ZTC(1,6)=ZTC(1,2)
          ZTC(2,3)=DCMPLX(UH*XIG*X11,0D0)
          ZTC(2,4)=DCMPLX(UH*XIG*X22,0D0)
          ZTC(2,5)=DCMPLX(UH*XIG*X12,0D0)
          ZTC(2,6)=DCMPLX(UH*XIG*X21,0D0)
          ZTC(3,4)=-SM1122
          ZTC(3,5)=-SM1112
          ZTC(3,6)=-SM1121
          ZTC(4,5)=-SM2212
          ZTC(4,6)=-SM2221
          ZTC(5,6)=-SM1221
          DO 210 I=1,5
            DO 200 J=I+1,6
               ZTC(J,I)=ZTC(I,J)
  200       CONTINUE
  210     CONTINUE
          CALL PYLDCM(ZTC,6,6,INDX,D)
          DO 230 I=1,6
            DO 220 J=1,6
              YTC(I,J)=(0D0,0D0)
              IF(I.EQ.J) YTC(I,J)=(1D0,0D0)
  220       CONTINUE
  230     CONTINUE
          DO 240 I=1,6
            CALL PYBKSB(ZTC,6,6,INDX,YTC(1,I))
  240     CONTINUE
          DGGU=YTC(1,1)
          DVVU=YTC(2,2)
          DGVU=YTC(1,2)
 
          IF(IMDL.EQ.1) THEN
            DQQS=DGGS+DVVS*DCMPLX(TANT3**2)-DGVS*DCMPLX(2D0*TANT3)
            DQQT=DGGT+DVVT*DCMPLX(TANT3**2)-DGVT*DCMPLX(2D0*TANT3)
            DQQU=DGGU+DVVU*DCMPLX(TANT3**2)-DGVU*DCMPLX(2D0*TANT3)
            DQTS=DGGS-DVVS-DGVS*DCMPLX(TANT3-1D0/TANT3)
            DQGS=DGGS-DGVS*DCMPLX(TANT3)
            DTGS=DGGS+DGVS*DCMPLX(1D0/TANT3)
          ELSE
            DQQS=DGGS+DVVS*DCMPLX(1D0/TANT3**2)+DGVS*DCMPLX(2D0/TANT3)
            DQQT=DGGT+DVVT*DCMPLX(1D0/TANT3**2)+DGVT*DCMPLX(2D0/TANT3)
            DQQU=DGGU+DVVU*DCMPLX(1D0/TANT3**2)+DGVU*DCMPLX(2D0/TANT3)
            DQTS=DGGS+DVVS*DCMPLX(1D0/TANT3**2)+DGVS*DCMPLX(2D0/TANT3)
            DQGS=DGGS+DGVS*DCMPLX(1D0/TANT3)
            DTGS=DGGS+DGVS*DCMPLX(1D0/TANT3)
          ENDIF
 
          SQDQTS=ABS(DQTS)**2
          SQDQQS=ABS(DQQS)**2
          SQDQQT=ABS(DQQT)**2
          SQDQQU=ABS(DQQU)**2
          SQDLGS=ABS(DCMPLX(SH)*DQGS-DCMPLX(1D0))**2
          REDLGS=DBLE(DQGS)
          SQDHGS=ABS(DCMPLX(SH)*DTGS-DCMPLX(1D0))**2
          REDHGS=DBLE(DTGS)
          SQDLGT=ABS(DCMPLX(TH)*DGGT-DCMPLX(1D0))**2
 
          SQDGGS=ABS(DGGS)**2
          SQDGGT=ABS(DGGT)**2
          SQDGGU=ABS(DGGU)**2
          REDGGS=DBLE(DGGS)
          REDGGT=DBLE(DGGT)
          REDGGU=DBLE(DGGU)
          REDGTU=DBLE(DGGU*DCONJG(DGGT))
          REDGSU=DBLE(DGGU*DCONJG(DGGS))
          REDGST=DBLE(DGGS*DCONJG(DGGT))
          REDQST=DBLE(DQQS*DCONJG(DQQT))
          REDQTU=DBLE(DQQT*DCONJG(DQQU))
        ENDIF
      ENDIF
 
 
C...Differential cross section expressions.
 
      IF(ISUB.LE.190) THEN
        IF(ISUB.EQ.149) THEN
C...g + g -> eta_tc
          KCTC=PYCOMP(KTECHN+331)
          CALL PYWIDT(KTECHN+331,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=COMFAC*0.5D0/((SH-PMAS(KCTC,1)**2)**2+HS**2)
          IF(ABS(SHR-PMAS(KCTC,1)).GT.PARP(48)*PMAS(KCTC,2)) FACBW=0D0
          HP=SH
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 250
          HI=HP*WDTP(3)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=HI*FACBW*HF
  250     CONTINUE
 
        ELSEIF(ISUB.EQ.165) THEN
C...q + qbar -> l+ + l- (including contact term for compositeness)
          ZRATR=XWC*SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)
          ZRATI=XWC*SH*GMMZ/((SH-SQMZ)**2+GMMZ**2)
          KFF=IABS(KFPR(ISUB,1))
          EF=KCHG(KFF,1)/3D0
          AF=SIGN(1D0,EF+0.1D0)
          VF=AF-4D0*EF*XWV
          VALF=VF+AF
          VARF=VF-AF
          FCOF=1D0
          IF(KFF.LE.10) FCOF=3D0
          WID2=1D0
          IF(KFF.EQ.6) WID2=WIDS(6,1)
          IF(KFF.EQ.7.OR.KFF.EQ.8) WID2=WIDS(KFF,1)
          IF(KFF.EQ.17.OR.KFF.EQ.18) WID2=WIDS(KFF,1)
          DO 260 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 260
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI+0.1D0)
            VI=AI-4D0*EI*XWV
            VALI=VI+AI
            VARI=VI-AI
            FCOI=1D0
            IF(IABS(I).LE.10) FCOI=FACA/3D0
            IF((ITCM(5).EQ.1.AND.IABS(I).LE.2).OR.ITCM(5).EQ.2) THEN
              FGZA=(EI*EF+VALI*VALF*ZRATR+RTCM(42)*SH/
     &        (AEM*RTCM(41)**2))**2+(VALI*VALF*ZRATI)**2+
     &        (EI*EF+VARI*VARF*ZRATR)**2+(VARI*VARF*ZRATI)**2
            ELSE
              FGZA=(EI*EF+VALI*VALF*ZRATR)**2+(VALI*VALF*ZRATI)**2+
     &        (EI*EF+VARI*VARF*ZRATR)**2+(VARI*VARF*ZRATI)**2
            ENDIF
            FGZB=(EI*EF+VALI*VARF*ZRATR)**2+(VALI*VARF*ZRATI)**2+
     &      (EI*EF+VARI*VALF*ZRATR)**2+(VARI*VALF*ZRATI)**2
            FGZAB=AEM**2*(FGZA*UH2/SH2+FGZB*TH2/SH2)
            IF((ITCM(5).EQ.3.AND.IABS(I).EQ.2).OR.(ITCM(5).EQ.4.AND.
     &      MOD(IABS(I),2).EQ.0)) FGZAB=FGZAB+SH2/(2D0*RTCM(41)**4)
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=COMFAC*FCOI*FCOF*FGZAB*WID2
  260     CONTINUE
 
        ELSEIF(ISUB.EQ.166) THEN
C...q + q'bar -> l + nu_l (including contact term for compositeness)
          WFAC=(1D0/4D0)*(AEM/XW)**2*UH2/((SH-SQMW)**2+GMMW**2)
          WCIFAC=WFAC+SH2/(4D0*RTCM(41)**4)
          KFF=IABS(KFPR(ISUB,1))
          FCOF=1D0
          IF(KFF.LE.10) FCOF=3D0
          DO 280 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 280
            IA=IABS(I)
            DO 270 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 270
              JA=IABS(J)
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 270
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 270
              FCOI=1D0
              IF(IA.LE.10) FCOI=VCKM((IA+1)/2,(JA+1)/2)*FACA/3D0
              WID2=1D0
              IF((I.GT.0.AND.MOD(I,2).EQ.0).OR.(J.GT.0.AND.
     &        MOD(J,2).EQ.0)) THEN
                IF(KFF.EQ.5) WID2=WIDS(6,2)
                IF(KFF.EQ.7) WID2=WIDS(8,2)*WIDS(7,3)
                IF(KFF.EQ.17) WID2=WIDS(18,2)*WIDS(17,3)
              ELSE
                IF(KFF.EQ.5) WID2=WIDS(6,3)
                IF(KFF.EQ.7) WID2=WIDS(8,3)*WIDS(7,2)
                IF(KFF.EQ.17) WID2=WIDS(18,3)*WIDS(17,2)
              ENDIF
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=COMFAC*FCOI*FCOF*WFAC*WID2
              IF((ITCM(5).EQ.3.AND.IA.LE.2.AND.JA.LE.2).OR.ITCM(5).EQ.4)
     &        SIGH(NCHN)=COMFAC*FCOI*FCOF*WCIFAC*WID2
  270       CONTINUE
  280     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.200) THEN
        IF(ISUB.EQ.191) THEN
C...q + qbar -> rho_tc0.
          KCTC=PYCOMP(KTECHN+113)
          SQMRHT=PMAS(KCTC,1)**2
          CALL PYWIDT(KTECHN+113,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=12D0*COMFAC/((SH-SQMRHT)**2+HS**2)
          IF(ABS(SHR-PMAS(KCTC,1)).GT.PARP(48)*PMAS(KCTC,2)) FACBW=0D0
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          ALPRHT=2.16D0*(3D0/ITCM(1))
          HP=(1D0/6D0)*(AEM**2/ALPRHT)*(SQMRHT**2/SH)
          XWRHT=(1D0-2D0*XW)/(4D0*XW*(1D0-XW))
          BWZR=XWRHT*SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)
          BWZI=XWRHT*SH*GMMZ/((SH-SQMZ)**2+GMMZ**2)
          DO 290 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 290
            IA=IABS(I)
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI+0.1D0)
            VI=AI-4D0*EI*XWV
            VALI=0.5D0*(VI+AI)
            VARI=0.5D0*(VI-AI)
            HI=HP*((EI+VALI*BWZR)**2+(VALI*BWZI)**2+
     &      (EI+VARI*BWZR)**2+(VARI*BWZI)**2)
            IF(IA.LE.10) HI=HI*FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=HI*FACBW*HF
  290     CONTINUE
 
        ELSEIF(ISUB.EQ.192) THEN
C...q + qbar' -> rho_tc+/-.
          KCTC=PYCOMP(KTECHN+213)
          SQMRHT=PMAS(KCTC,1)**2
          CALL PYWIDT(KTECHN+213,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=12D0*COMFAC/((SH-SQMRHT)**2+HS**2)
          IF(ABS(SHR-PMAS(KCTC,1)).GT.PARP(48)*PMAS(KCTC,2)) FACBW=0D0
          ALPRHT=2.16D0*(3D0/ITCM(1))
          HP=(1D0/6D0)*(AEM**2/ALPRHT)*(SQMRHT**2/SH)*
     &    (0.25D0/XW**2)*SH**2/((SH-SQMW)**2+GMMW**2)
          DO 310 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 310
            IA=IABS(I)
            DO 300 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 300
              JA=IABS(J)
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 300
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 300
              KCHR=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              HF=SHR*(WDTE(0,1)+WDTE(0,(5-KCHR)/2)+WDTE(0,4))
              HI=HP
              IF(IA.LE.10) HI=HI*VCKM((IA+1)/2,(JA+1)/2)*FACA/3D0
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=HI*FACBW*HF
  300       CONTINUE
  310     CONTINUE
 
        ELSEIF(ISUB.EQ.193) THEN
C...q + qbar -> omega_tc0.
          KCTC=PYCOMP(KTECHN+223)
          SQMOMT=PMAS(KCTC,1)**2
          CALL PYWIDT(KTECHN+223,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=12D0*COMFAC/((SH-SQMOMT)**2+HS**2)
          IF(ABS(SHR-PMAS(KCTC,1)).GT.PARP(48)*PMAS(KCTC,2)) FACBW=0D0
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          ALPRHT=2.16D0*(3D0/ITCM(1))
          HP=(1D0/6D0)*(AEM**2/ALPRHT)*(SQMOMT**2/SH)*
     &    (2D0*RTCM(2)-1D0)**2
          BWZR=(0.5D0/(1D0-XW))*SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)
          BWZI=(0.5D0/(1D0-XW))*SH*GMMZ/((SH-SQMZ)**2+GMMZ**2)
          DO 320 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 320
            IA=IABS(I)
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI+0.1D0)
            VI=AI-4D0*EI*XWV
            VALI=0.5D0*(VI+AI)
            VARI=0.5D0*(VI-AI)
            HI=HP*((EI-VALI*BWZR)**2+(VALI*BWZI)**2+
     &      (EI-VARI*BWZR)**2+(VARI*BWZI)**2)
            IF(IA.LE.10) HI=HI*FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=HI*FACBW*HF
  320     CONTINUE
 
        ELSEIF(ISUB.EQ.194) THEN
C...f + fbar -> f' + fbar' via s-channel rho_tc, omega_tc a_T0.
C...Default final state is e+e-
          KFA=KFPR(ISUBSV,1)
          ALPRHT=2.16D0*(3D0/ITCM(1))
          HP=AEM**2*COMFAC

          SN2W=2D0*SQRT(XW*XW1)
C          TANW=SQRT(PARU(102)/(1D0-PARU(102)))
C          CT2W=(1D0-2D0*PARU(102))/(2D0*PARU(102)/TANW)
 
          QUPD=2D0*RTCM(2)-1D0
          FAR=SQRT(AEM/ALPRHT)
          FAO=FAR*QUPD
          FZR=FAR*CT2W
          FZO=-FAO*TANW
C...RTCM(47) is the ratio g_{rho_T}/g_{a_T}
          FZX=-FAR/SN2W*RTCM(47)
          SFAR=FAR**2
          SFAO=FAO**2
          SFZR=FZR**2
          SFZO=FZO**2
          SFZX=FZX**2
          CALL PYWIDT(23,SH,WDTP,WDTE)
          SSMZ=DCMPLX(1D0-PMAS(23,1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+113,SH,WDTP,WDTE)
          SSMR=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+113),1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+223,SH,WDTP,WDTE)
          SSMO=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+223),1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+115,SH,WDTP,WDTE)
          SSMX=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+115),1)**2/SH,WDTP(0)/SHR)
C...Propagator including a_T^0
          DETD=(FAR*FZO-FAO*FZR)**2+SSMZ*SSMR*SSMO-SFZR*SSMO-
     $    SFZO*SSMR-SFAR*SSMO*SSMZ-SFAO*SSMR*SSMZ
C...Add in techni-a contribution
          DETD=SSMX*DETD-SFZX*(SSMR*SSMO-SFAO*SSMR-SFAR*SSMO)
          DAA=(-SSMX*(SFZO*SSMR+SFZR*SSMO-SSMO*SSMR*SSMZ)-
     $     SFZX*SSMR*SSMO)/DETD/SH
          DZZ=-(SFAO*SSMR+SFAR*SSMO-SSMO*SSMR)/DETD/SH*SSMX
          DAZ=(FAR*FZR*SSMO+FAO*FZO*SSMR)/DETD/SH*SSMX
 
          XWRHT=1D0/(4D0*XW*(1D0-XW))
          KFF=IABS(KFPR(ISUB,1))
          EF=KCHG(KFF,1)/3D0
          AF=SIGN(1D0,EF+0.1D0)
          VF=AF-4D0*EF*XWV
          VALF=0.5D0*(VF+AF)
          VARF=0.5D0*(VF-AF)
          FCOF=1D0
          IF(KFF.LE.10) FCOF=3D0
 
          WID2=1D0
          IF(KFF.GE.6.AND.KFF.LE.8) WID2=WIDS(KFF,1)
          IF(KFF.EQ.17.OR.KFF.EQ.18) WID2=WIDS(KFF,1)
          DZZ=DZZ*DCMPLX(XWRHT,0D0)
          DAZ=DAZ*DCMPLX(SQRT(XWRHT),0D0)
 
          DO 330 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 330
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI+0.1D0)
            VI=AI-4D0*EI*XWV
            VALI=0.5D0*(VI+AI)
            VARI=0.5D0*(VI-AI)
            FCOI=FCOF
            IF(IABS(I).LE.10) FCOI=FCOI/3D0
            DIFLL=ABS(EI*EF*DAA+VALI*VALF*DZZ+DAZ*(EI*VALF+EF*VALI))**2
            DIFRR=ABS(EI*EF*DAA+VARI*VARF*DZZ+DAZ*(EI*VARF+EF*VARI))**2
            DIFLR=ABS(EI*EF*DAA+VALI*VARF*DZZ+DAZ*(EI*VARF+EF*VALI))**2
            DIFRL=ABS(EI*EF*DAA+VARI*VALF*DZZ+DAZ*(EI*VALF+EF*VARI))**2
            FACSIG=(DIFLL+DIFRR)*((UH-SQM4)**2+SH*SQM4)+
     &      (DIFLR+DIFRL)*((TH-SQM3)**2+SH*SQM3)
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=HP*FCOI*FACSIG*WID2
  330     CONTINUE
 
        ELSEIF(ISUB.EQ.195) THEN
C...f + fbar' -> f'' + fbar''' via s-channel rho_tc+, a_T+
          KFA=KFPR(ISUBSV,1)
          KFB=KFA+1
          ALPRHT=2.16D0*(3D0/ITCM(1))
          FACTC=COMFAC*(AEM**2/12D0/XW**2)*(UH-SQM3)*(UH-SQM4)*3D0
 
          FWR=SQRT(AEM/ALPRHT)/(2D0*SQRT(XW))
C...RTCM(47) is the ratio g_{rho_T}/g_{a_T}
C
C...Propagator including a_T^+
          FWX=-FWR*RTCM(47)
          CALL PYWIDT(24,SH,WDTP,WDTE)
          SSMZ=DCMPLX(1D0-PMAS(24,1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+213,SH,WDTP,WDTE)
          SSMR=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+213),1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+215,SH,WDTP,WDTE)
          SSMX=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+215),1)**2/SH,WDTP(0)/SHR)
          DETD=SSMX*(SSMZ*SSMR-DCMPLX(FWR**2,0D0))-
     &     DCMPLX(FWX**2,0D0)*SSMR
          DWW=SSMR*SSMX/DETD/SH
          FCOF=1D0
          IF(KFA.LE.8) FCOF=3D0
          HP=FACTC*ABS(DWW)**2*FCOF
 
          DO 350 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 350
            IA=IABS(I)
            DO 340 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 340
              JA=IABS(J)
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 340
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 340
              KCHR=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              HI=HP
              IF(IA.LE.10) HI=HI*VCKM((IA+1)/2,(JA+1)/2)/3D0
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=HI*WIDS(KFA,(5-KCHR)/2)*WIDS(KFB,(5+KCHR)/2)
  340       CONTINUE
  350     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.380) THEN
        ALPRHT=2.16D0*(3D0/ITCM(1))
        IF(ISUB.EQ.361) THEN
          FAR=SQRT(AEM/ALPRHT)
          FAO=FAR*QUPD
          FZR=FAR*CT2W
          FZO=-FAO*TANW
C...RTCM(47) is the ratio g_{rho_T}/g_{a_T}
          FZX=-FAR/SN2W*RTCM(47)
          SFAR=FAR**2
          SFAO=FAO**2
          SFZR=FZR**2
          SFZO=FZO**2
          SFZX=FZX**2
          CALL PYWIDT(23,SH,WDTP,WDTE)
          SSMZ=DCMPLX(1D0-PMAS(23,1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+113,SH,WDTP,WDTE)
          SSMR=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+113),1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+223,SH,WDTP,WDTE)
          SSMO=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+223),1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+115,SH,WDTP,WDTE)
          SSMX=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+115),1)**2/SH,WDTP(0)/SHR)
          DETD=(FAR*FZO-FAO*FZR)**2+SSMZ*SSMR*SSMO-SFZR*SSMO-
     $    SFZO*SSMR-SFAR*SSMO*SSMZ-SFAO*SSMR*SSMZ
C...Add in techni-a contribution
          DETD=SSMX*DETD-SFZX*(SSMR*SSMO-SFAO*SSMR-SFAR*SSMO)
          DARHO=-(SSMX*(-FAR*SFZO+FAO*FZO*FZR+FAR*SSMO*SSMZ)-
     $     SFZX*FAR*SSMO)/DETD/SH
          DZRHO=-(-FZR*SFAO+FAO*FZO*FAR+FZR*SSMO)/DETD/SH*SSMX
          DAOME=-(SSMX*(-FAO*SFZR+FAR*FZO*FZR+FAO*SSMR*SSMZ)-
     $     SFZX*FAO*SSMR)/DETD/SH
          DZOME=-(-FZO*SFAR+FAR*FAO*FZR+FZO*SSMR)/DETD/SH*SSMX
          DAAST=-FZX*(FAO*FZO*SSMR+FAR*FZR*SSMO)/DETD/SH
          DZAST=-FZX*(SSMR*SSMO-SFAO*SSMR-SFAR*SSMO)/DETD/SH
          DAA=(-SSMX*(SFZO*SSMR+SFZR*SSMO-SSMO*SSMR*SSMZ)-
     $     SFZX*SSMR*SSMO)/DETD/SH
          DZZ=-(SFAO*SSMR+SFAR*SSMO-SSMO*SSMR)/DETD/SH*SSMX
          DAZ=(FAR*FZR*SSMO+FAO*FZO*SSMR)/DETD/SH*SSMX
 
C...f + fbar -> gamma pi_tc, gamma pi_tc', Z pi_tc, Z pi_tc',
C...W+W-, W pi_tc, pi_T pi_T, etc.
          FACA=(SH**2*BE34**2-(TH-UH)**2)
          VFAC=(TH**2+UH**2-2D0*SQM3*SQM4)
          AFAC=(TH**2+UH**2-2D0*SQM3*SQM4+4D0*SH*SQM3)
          FANOM=SQRT(PARU(1)*AEM)*ITCM(1)/PARU(2)**2/RTCM(1)
          HP=(1D0/24D0)*AEM**2*COMFAC*3D0*SH 
          DO 370 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 370
            IA=IABS(I)
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI+0.1D0)
            VI=AI-4D0*EI*XWV
            VALI=0.25D0*(VI+AI) ! = \zeta_{iL} in PRD67-115011
            VARI=0.25D0*(VI-AI) ! = \zeta_{iR} in PRD67-115011
C...........Eqs. (5) and (6) in LSTC-rates.pdf
            F2L=(EI*DARHO+VALI*DZRHO/SQRT(XW*XW1))*VRGP
            F2L=F2L+(EI*DAOME+VALI*DZOME/SQRT(XW*XW1))*VOGP
            F2L=F2L+(EI*DAAST+VALI*DZAST/SQRT(XW*XW1))*VXGP
            F2L=F2L+FANOM*(VAGP*(EI*DAA+VALI*DAZ/SQRT(XW*XW1))+
     $                    VZGP*(EI*DAZ+VALI*DZZ/SQRT(XW*XW1)))
            F2R=(EI*DARHO+VARI*DZRHO/SQRT(XW*XW1))*VRGP
            F2R=F2R+(EI*DAOME+VARI*DZOME/SQRT(XW*XW1))*VOGP
            F2R=F2R+(EI*DAAST+VARI*DZAST/SQRT(XW*XW1))*VXGP
            F2R=F2R+FANOM*(VAGP*(EI*DAA+VARI*DAZ/SQRT(XW*XW1))+
     $                    VZGP*(EI*DAZ+VARI*DZZ/SQRT(XW*XW1)))
            HI=(ABS(F2L)**2+ABS(F2R)**2)*VFAC
C...........Eqs. (5) and (7) in LSTC-rates.pdf
            F2L=(EI*DARHO+VALI*DZRHO/SQRT(XW*XW1))*ARGP
            F2L=F2L+(EI*DAOME+VALI*DZOME/SQRT(XW*XW1))*AOGP
            F2L=F2L+(EI*DAAST+VALI*DZAST/SQRT(XW*XW1))*AXGP
            F2R=(EI*DARHO+VARI*DZRHO/SQRT(XW*XW1))*ARGP
            F2R=F2R+(EI*DAOME+VARI*DZOME/SQRT(XW*XW1))*AOGP
            F2R=F2R+(EI*DAAST+VARI*DZAST/SQRT(XW*XW1))*AXGP
            HJ=(ABS(F2L)**2+ABS(F2R)**2)*AFAC
C
C...........Eqs. (24) in PRD67-115011 with DAA, etc.terms dropped.
C
c$$$            F2L=EI*(DARHO/FAR+(DAA+CT2W*DAZ))+
c$$$     $      VALI*(CT2W*DZRHO/FZR+(CT2W*DZZ+DAZ))/SQRT(XW*XW1)
c$$$            F2R=EI*(DARHO/FAR+(DAA+CT2W*DAZ))+
c$$$     $      VARI*(CT2W*DZRHO/FZR+(CT2W*DZZ+DAZ))/SQRT(XW*XW1)
            F2L=EI*DARHO/FAR + VALI*CT2W*DZRHO/FZR/SQRT(XW*XW1)
            F2R=EI*DARHO/FAR + VARI*CT2W*DZRHO/FZR/SQRT(XW*XW1)
            HK=(ABS(F2L)**2+ABS(F2R)**2)*2D0*FACA*CAB2/SH
            HI=HI+HJ+HK
            IF(IA.LE.10) HI=HI/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            IF(KFA.EQ.KFB) THEN
               SIGH(NCHN)=HI*HP*WIDS(PYCOMP(KFA),1)
            ELSEIF(ISUBSV.EQ.362.OR.ISUBSV.EQ.368) THEN
               SIGH(NCHN)=HI*HP*WIDS(PYCOMP(KFA),2)*WIDS(PYCOMP(KFB),3)
               NCHN=NCHN+1
               ISIG(NCHN,1)=I
               ISIG(NCHN,2)=-I
               ISIG(NCHN,3)=2
               SIGH(NCHN)=HI*HP*WIDS(PYCOMP(KFA),3)*WIDS(PYCOMP(KFB),2)
            ELSE 
               SIGH(NCHN)=HI*HP*WIDS(PYCOMP(KFA),2)*WIDS(PYCOMP(KFB),2)
            ENDIF
  370     CONTINUE
 
        ELSEIF(ISUB.EQ.370) THEN
C...f + fbar' -> W_L Z_L, W_L Z_T, W_T, Z_L, W_L pi_tc, Z_L pi_tc, pi_tc pi_tc
C...f + fbar' -> gamma pi_tc, etc.
          FACA=(SH**2*BE34**2-(TH-UH)**2)
          FANOM=SQRT(PARU(1)*AEM)*ITCM(1)/PARU(2)**2/RTCM(1)
          VFAC=(TH**2+UH**2-2D0*SQM3*SQM4)
          AFAC=(TH**2+UH**2-2D0*SQM3*SQM4+4D0*SH*SQM3)
          ALPRHT=2.16D0*(3D0/ITCM(1))
          FACHP=(1D0/48D0)*AEM**2/XW*COMFAC*3D0*SH
          FWR=SQRT(AEM/ALPRHT)/(2D0*SQRT(XW))
C...RTCM(47) is the ratio g_{rho_T}/g_{a_T}
          FWX=-FWR*RTCM(47)
          CALL PYWIDT(24,SH,WDTP,WDTE)
          SSMZ=DCMPLX(1D0-PMAS(24,1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+213,SH,WDTP,WDTE)
          SSMR=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+213),1)**2/SH,WDTP(0)/SHR)
          CALL PYWIDT(KTECHN+215,SH,WDTP,WDTE)
          SSMX=DCMPLX(1D0-PMAS(PYCOMP(KTECHN+215),1)**2/SH,WDTP(0)/SHR)
          DETD=SSMX*(SSMZ*SSMR-DCMPLX(FWR**2,0D0))-
     &     DCMPLX(FWX**2,0D0)*SSMR
          DWW=SSMR*SSMX/DETD/SH
          DWRHO=-DCMPLX(FWR,0D0)*SSMX/DETD/SH
          DWAST=-DCMPLX(FWX,0D0)*SSMR/DETD/SH
          HP=FACHP*(AFAC*ABS(DWRHO*ARGP+DWAST*AXGP)**2+
     $    VFAC*ABS(FANOM*DWW*VWGP+DWRHO*VRGP+DWAST*VXGP)**2)
C
C...........Eq. (25) in PRD67-115011 with DWW term dropped.
C
c$$$          HP=HP+.5D0*FACHP*CAB2*FACA/XW/SH*ABS(DWW + DWRHO/FWR)**2
          HP=HP+.5D0*FACHP*CAB2*FACA/XW/SH*ABS(DWRHO/FWR)**2
C...Add in W_L Z_T axial and vector contributions.
          IF(ISUBSV.EQ.370) HP=HP+FACHP*RTCM(3)**2*(
     $    (TH**2+UH**2-2D0*SQM3*SQM4+4D0*SH*SQM4)*     !AFAC w/ switched masses.
     $    ABS(DWRHO/RTCM(13)-DWAST/RTCM(49)*CS2W)**2/SN2W**2+
     $    VFAC*QUPD**2*XW/XW1*ABS(DWRHO)**2/RTCM(12)**2)
          DO 410 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 410
            IA=IABS(I)
            DO 400 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 400
              JA=IABS(J)
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 400
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 400
              KCHR=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              HI=HP
              IF(IA.LE.10) HI=HI*VCKM((IA+1)/2,(JA+1)/2)/3D0
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              IF(ISUBSV.EQ.374.OR.ISUBSV.EQ.378) THEN
                SIGH(NCHN)=HI*WIDS(PYCOMP(KFA),(5-KCHR)/2)
              ELSE
                SIGH(NCHN)=HI*WIDS(PYCOMP(KFA),(5-KCHR)/2)*
     &          WIDS(PYCOMP(KFB),2)
              ENDIF
  400       CONTINUE
  410     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.390) THEN
        IF(ISUB.EQ.381) THEN
C...f + f' -> f + f' (g exchange)
          FACQQ1=COMFAC*AS**2*4D0/9D0*(SH2+UH2)*SQDQQT
          FACQQB=COMFAC*AS**2*4D0/9D0*((SH2+UH2)*SQDQQT*FACA-
     &    MSTP(34)*2D0/3D0*UH2*REDQST)
          FACQQ2=COMFAC*AS**2*4D0/9D0*(SH2+TH2)*SQDQQU
          FACQQI=-COMFAC*AS**2*4D0/9D0*MSTP(34)*2D0/3D0*SH2/(TH*UH)
          RATQQI=(FACQQ1+FACQQ2+FACQQI)/(FACQQ1+FACQQ2)
          IF(ITCM(5).GE.1.AND.ITCM(5).LE.4) THEN
C...Modifications from contact interactions (compositeness)
            FACCI1=FACQQ1+COMFAC*(SH2/RTCM(41)**4)
            FACCIB=FACQQB+COMFAC*(8D0/9D0)*(AS*RTCM(42)/RTCM(41)**2)*
     &      (UH2/TH+UH2/SH)+COMFAC*(5D0/3D0)*(UH2/RTCM(41)**4)
            FACCI2=FACQQ2+COMFAC*(8D0/9D0)*(AS*RTCM(42)/RTCM(41)**2)*
     &      (SH2/TH+SH2/UH)+COMFAC*(5D0/3D0)*(SH2/RTCM(41)**4)
            FACCI3=FACQQ1+COMFAC*(UH2/RTCM(41)**4)
            RATCII=(FACCI1+FACCI2+FACQQI)/(FACCI1+FACCI2)
          ELSEIF(ITCM(5).EQ.5) THEN
            FACCI1=FACQQ1
            FACCIB=FACQQB
            FACCI2=FACQQ2
            FACCI3=FACQQ1
CSM.......Check this change from
CSM            RATCII=1D0
            RATCII=RATQQI
          ENDIF
          DO 430 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 430
            DO 420 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).EQ.0) GOTO 420
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              IF(ITCM(5).LE.0.OR.(ITCM(5).EQ.1.AND.(IA.GE.3.OR.
     &        JA.GE.3))) THEN
                SIGH(NCHN)=FACQQ1
                IF(I.EQ.-J) SIGH(NCHN)=FACQQB
              ELSE
                SIGH(NCHN)=FACCI1
                IF(I*J.LT.0) SIGH(NCHN)=FACCI3
                IF(I.EQ.-J) SIGH(NCHN)=FACCIB
              ENDIF
              IF(I.EQ.J) THEN
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=2
                IF(ITCM(5).LE.0.OR.(ITCM(5).EQ.1.AND.IA.GE.3)) THEN
                  SIGH(NCHN-1)=0.5D0*FACQQ1*RATQQI
                  SIGH(NCHN)=0.5D0*FACQQ2*RATQQI
                ELSE
                  SIGH(NCHN-1)=0.5D0*FACCI1*RATCII
                  SIGH(NCHN)=0.5D0*FACCI2*RATCII
                ENDIF
              ENDIF
  420       CONTINUE
  430     CONTINUE
 
        ELSEIF(ISUB.EQ.382) THEN
C...f + fbar -> f' + fbar' (q + qbar -> q' + qbar' only)
          CALL PYWIDT(21,SH,WDTP,WDTE)
          FACQQF=COMFAC*AS**2*4D0/9D0*(TH2+UH2)
          FACQQB=FACQQF*SQDQQS*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          IF(ITCM(5).EQ.1) THEN
C...Modifications from contact interactions (compositeness)
            FACCIB=FACQQB
            DO 440 I=1,2
              FACCIB=FACCIB+COMFAC*(UH2/RTCM(41)**4)*(WDTE(I,1)+
     &        WDTE(I,2)+WDTE(I,4))
  440       CONTINUE
          ELSEIF(ITCM(5).GE.2.AND.ITCM(5).LE.4) THEN
            FACCIB=FACQQB+COMFAC*(UH2/RTCM(41)**4)*
     &      (WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          ELSEIF(ITCM(5).EQ.5) THEN
            FACQQB=FACQQF*SQDQQS*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4)-
     &      WDTE(5,1)-WDTE(5,2)-WDTE(5,4))
            FACCIB=FACQQF*SQDQTS*(WDTE(5,1)+WDTE(5,2)+WDTE(5,4))
          ENDIF
          DO 450 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 450
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            IF(ITCM(5).LE.0.OR.(ITCM(5).EQ.1.AND.IABS(I).GE.3)) THEN
              SIGH(NCHN)=FACQQB
            ELSEIF(ITCM(5).EQ.5) THEN
              SIGH(NCHN)=FACQQB
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=-I
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACCIB
            ELSE
              SIGH(NCHN)=FACCIB
            ENDIF
  450     CONTINUE
 
        ELSEIF(ISUB.EQ.383) THEN
C...f + fbar -> g + g (q + qbar -> g + g only)
          FACGG1=COMFAC*AS**2*32D0/27D0*(UH/TH-(2D0+MSTP(34)*1D0/4D0)*
     &    UH2/SH2+9D0/4D0*TH*UH/SH2*SQDLGS)
          FACGG2=COMFAC*AS**2*32D0/27D0*(TH/UH-(2D0+MSTP(34)*1D0/4D0)*
     &    TH2/SH2+9D0/4D0*TH*UH/SH2*SQDLGS)
          IF(ITCM(5).EQ.5) THEN
            FACGG3=COMFAC*AS**2*32D0/27D0*(UH/TH-(2D0+MSTP(34)*1D0/4D0)*
     &      UH2/SH2+9D0/4D0*TH*UH/SH2*SQDHGS)
            FACGG4=COMFAC*AS**2*32D0/27D0*(TH/UH-(2D0+MSTP(34)*1D0/4D0)*
     &      TH2/SH2+9D0/4D0*TH*UH/SH2*SQDHGS)
          ENDIF
          DO 460 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 460
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=0.5D0*FACGG1
            IF(ITCM(5).EQ.5.AND.IABS(I).EQ.5) SIGH(NCHN)=0.5D0*FACGG3
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=2
            SIGH(NCHN)=0.5D0*FACGG2
            IF(ITCM(5).EQ.5.AND.IABS(I).EQ.5) SIGH(NCHN)=0.5D0*FACGG4
  460     CONTINUE
 
        ELSEIF(ISUB.EQ.384) THEN
C...f + g -> f + g (q + g -> q + g only)
          FACQG1=COMFAC*AS**2*4D0/9D0*((2D0+MSTP(34)*1D0/4D0)*UH2/TH2-
     &    UH/SH-9D0/4D0*SH*UH/TH2*SQDLGT)*FACA
          FACQG2=COMFAC*AS**2*4D0/9D0*((2D0+MSTP(34)*1D0/4D0)*SH2/TH2-
     &    SH/UH-9D0/4D0*SH*UH/TH2*SQDLGT)
          DO 480 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.10) GOTO 480
            DO 470 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 470
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 470
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQG1
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACQG2
  470       CONTINUE
  480     CONTINUE
 
        ELSEIF(ISUB.EQ.385) THEN
C...g + g -> f + fbar (g + g -> q + qbar only)
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 500
          IDC0=MDCY(21,2)-1
C...Begin by d, u, s flavours.
          FLAVWT=0D0
          IF(MDME(IDC0+1,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(1,1)**2/SH))
          IF(MDME(IDC0+2,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(2,1)**2/SH))
          IF(MDME(IDC0+3,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(3,1)**2/SH))
          FACQQ1=COMFAC*AS**2*1D0/6D0*(UH/TH-(2D0+MSTP(34)*1D0/4D0)*
     &    UH2/SH2+9D0/4D0*TH*UH/SH2*SQDLGS)*FLAVWT*FACA
          FACQQ2=COMFAC*AS**2*1D0/6D0*(TH/UH-(2D0+MSTP(34)*1D0/4D0)*
     &    TH2/SH2+9D0/4D0*TH*UH/SH2*SQDLGS)*FLAVWT*FACA
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ1
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=2
          SIGH(NCHN)=FACQQ2
C...Next c and b flavours: modified that and uhat for fixed
C...cos(theta-hat).
          DO 490 IFL=4,5
          SQMAVG=PMAS(IFL,1)**2
          IF(MDME(IDC0+IFL,1).GE.1.AND.SH.GT.4.04D0*SQMAVG) THEN
            BE34=SQRT(1D0-4D0*SQMAVG/SH)
            THQ=-0.5D0*SH*(1D0-BE34*CTH)
            UHQ=-0.5D0*SH*(1D0+BE34*CTH)
            THUHQ=THQ*UHQ-SQMAVG*SH
            IF(MSTP(34).EQ.0) THEN
              FACQQ1=UHQ/THQ-2D0*UHQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/THQ**2
              FACQQ2=THQ/UHQ-2D0*THQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/UHQ**2
            ELSE
              FACQQ1=UHQ/THQ-2.25D0*UHQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &        THQ**2+0.5D0*SQMAVG*(THQ+SQMAVG)/THQ**2-SQMAVG**2/(SH*THQ)
              FACQQ2=THQ/UHQ-2.25D0*THQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &        UHQ**2+0.5D0*SQMAVG*(UHQ+SQMAVG)/UHQ**2-SQMAVG**2/(SH*UHQ)
            ENDIF
            IF(ITCM(5).GE.5) THEN
              IF(IFL.EQ.4) THEN
                FACQQ1=FACQQ1+2.25D0*SQMAVG*(THQ-UHQ)/(SH*THQ)*REDLGS+
     &          2.25D0*THQ*UHQ/SH2*SQDLGS
                FACQQ2=FACQQ2+2.25D0*SQMAVG*(UHQ-THQ)/(SH*UHQ)*REDLGS+
     &          2.25D0*THQ*UHQ/SH2*SQDLGS
              ELSE
                FACQQ1=FACQQ1+2.25D0*SQMAVG*(THQ-UHQ)/(SH*THQ)*REDHGS+
     &          2.25D0*THQ*UHQ/SH2*SQDHGS
                FACQQ2=FACQQ2+2.25D0*SQMAVG*(UHQ-THQ)/(SH*UHQ)*REDHGS+
     &          2.25D0*THQ*UHQ/SH2*SQDHGS
              ENDIF
            ENDIF
            FACQQ1=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ1*BE34
            FACQQ2=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ2*BE34
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1+2*(IFL-3)
            SIGH(NCHN)=FACQQ1
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=2+2*(IFL-3)
            SIGH(NCHN)=FACQQ2
          ENDIF
  490     CONTINUE
  500     CONTINUE
 
        ELSEIF(ISUB.EQ.386) THEN
C...g + g -> g + g
          IF(ITCM(5).LE.4) THEN
            FACGG1=COMFAC*AS**2*9D0/4D0*(SH2/TH2+2D0*SH/TH+3D0+
     &      2D0*TH/SH+TH2/SH2)*FACA
            FACGG2=COMFAC*AS**2*9D0/4D0*(UH2/SH2+2D0*UH/SH+3D0+
     &      2D0*SH/UH+SH2/UH2)*FACA
            FACGG3=COMFAC*AS**2*9D0/4D0*(TH2/UH2+2D0*TH/UH+3D0+
     &      2D0*UH/TH+UH2/TH2)
          ELSE
            GST=  (12D0 + 40D0*TH/SH + 56D0*TH2/SH2 + 32D0*TH**3/SH**3 +
     &      16D0*TH**4/SH**4 + SQDGGS*(4D0*SH2 + 16D0*SH*TH + 16D0*TH2)+
     &      4D0*REDGST*(SH + 2D0*TH)*
     &      (2D0*SH**3 - 3D0*SH2*TH - 2D0*SH*TH2 + 2D0*TH**3)/SH2 +
     &      2D0*REDGGS*(2D0*SH - 12D0*TH2/SH - 8D0*TH**3/SH2) +
     &      2D0*REDGGT*(4D0*SH - 22D0*TH - 68D0*TH2/SH - 60D0*TH**3/SH2-
     &      32D0*TH**4/SH**3 - 16D0*TH**5/SH**4) +
     &      SQDGGT*(16D0*SH2 + 16D0*SH*TH + 68D0*TH2 + 144D0*TH**3/SH +
     &      96D0*TH**4/SH2 + 32D0*TH**5/SH**3 + 16D0*TH**6/SH**4))/16D0
            GSU=  (12D0 + 40D0*UH/SH + 56D0*UH2/SH2 + 32D0*UH**3/SH**3 +
     &      16D0*UH**4/SH**4 + SQDGGS*(4D0*SH2 + 16D0*SH*UH + 16D0*UH2)+
     &      4D0*REDGSU*(SH + 2D0*UH)*
     &      (2D0*SH**3 - 3D0*SH2*UH - 2D0*SH*UH2 + 2D0*UH**3)/SH2 +
     &      2D0*REDGGS*(2D0*SH - 12D0*UH2/SH - 8D0*UH**3/SH2) +
     &      2D0*REDGGU*(4D0*SH - 22D0*UH - 68D0*UH2/SH - 60D0*UH**3/SH2-
     &      32D0*UH**4/SH**3 - 16D0*UH**5/SH**4) +
     &      SQDGGU*(16D0*SH2 + 16D0*SH*UH + 68D0*UH2 + 144D0*UH**3/SH +
     &      96D0*UH**4/SH2 + 32D0*UH**5/SH**3 + 16D0*UH**6/SH**4))/16D0
            GUT=  (12D0 - 16D0*TH*(TH - UH)**2*UH/SH**4 +
     &      4D0*REDGGU*(2D0*TH**5 - 15D0*TH**4*UH - 48D0*TH**3*UH2 -
     &      58D0*TH2*UH**3 - 10D0*TH*UH**4 + UH**5)/SH**4 +
     &      4D0*REDGGT*(TH**5 - 10D0*TH**4*UH - 58D0*TH**3*UH2 -
     &      48D0*TH2*UH**3 - 15D0*TH*UH**4 + 2D0*UH**5)/SH**4 +
     &      4D0*SQDGGU*(4D0*TH**6 + 20D0*TH**5*UH + 57D0*TH**4*UH2 +
     &      72D0*TH**3*UH**3+ 38D0*TH2*UH**4+4D0*TH*UH**5 +UH**6)/SH**4+
     &      4D0*SQDGGT*(4D0*UH**6 + 4D0*TH**5*UH + 38D0*TH**4*UH2 +
     &      72D0*TH**3*UH**3 +57D0*TH2*UH**4+20D0*TH*UH**5+TH**6)/SH**4+
     &      2D0*REDGTU*((TH - UH)**2* (TH**4 + 20D0*TH**3*UH +
     &      30D0*TH2*UH2 + 20D0*TH*UH**3 + UH**4) +
     &      SH2*(7D0*TH**4 + 52D0*TH**3*UH + 274D0*TH2*UH2 +
     &      52D0*TH*UH**3 + 7D0*UH**4))/(2D0*SH**4))/16D0
            FACGG1=COMFAC*AS**2*9D0/4D0*GST*FACA
            FACGG2=COMFAC*AS**2*9D0/4D0*GSU*FACA
            FACGG3=COMFAC*AS**2*9D0/4D0*GUT
          ENDIF
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 510
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=0.5D0*FACGG1
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=2
          SIGH(NCHN)=0.5D0*FACGG2
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=3
          SIGH(NCHN)=0.5D0*FACGG3
  510     CONTINUE
 
        ELSEIF(ISUB.EQ.387) THEN
C...q + qbar -> Q + Qbar
          SQMAVG=0.5D0*(SQM3+SQM4)-0.25D0*(SQM3-SQM4)**2/SH
          THQ=-0.5D0*SH*(1D0-BE34*CTH)
          UHQ=-0.5D0*SH*(1D0+BE34*CTH)
          FACQQB=COMFAC*AS**2*4D0/9D0*((THQ**2+UHQ**2)/SH2+
     &    2D0*SQMAVG/SH)
          IF(ITCM(5).GE.5) THEN
            IF(MINT(55).EQ.5.OR.MINT(55).EQ.6) THEN
              FACQQB=FACQQB*SH2*SQDQTS
            ELSE
              FACQQB=FACQQB*SH2*SQDQQS
            ENDIF
          ENDIF
          IF(MSTP(35).GE.1) FACQQB=FACQQB*PYHFTH(SH,SQMAVG,0D0)
          WID2=1D0
          IF(MINT(55).EQ.6) WID2=WIDS(6,1)
          IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=WIDS(MINT(55),1)
          FACQQB=FACQQB*WID2
          DO 520 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 520
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQB
  520     CONTINUE
 
        ELSEIF(ISUB.EQ.388) THEN
C...g + g -> Q + Qbar
          SQMAVG=0.5D0*(SQM3+SQM4)-0.25D0*(SQM3-SQM4)**2/SH
          THQ=-0.5D0*SH*(1D0-BE34*CTH)
          UHQ=-0.5D0*SH*(1D0+BE34*CTH)
          THUHQ=THQ*UHQ-SQMAVG*SH
          IF(MSTP(34).EQ.0) THEN
            FACQQ1=UHQ/THQ-2D0*UHQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/THQ**2
            FACQQ2=THQ/UHQ-2D0*THQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/UHQ**2
          ELSE
            FACQQ1=UHQ/THQ-2.25D0*UHQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &      THQ**2+0.5D0*SQMAVG*(THQ+SQMAVG)/THQ**2-SQMAVG**2/(SH*THQ)
            FACQQ2=THQ/UHQ-2.25D0*THQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &      UHQ**2+0.5D0*SQMAVG*(UHQ+SQMAVG)/UHQ**2-SQMAVG**2/(SH*UHQ)
          ENDIF
          IF(ITCM(5).GE.5) THEN
            IF(MINT(55).EQ.5.OR.MINT(55).EQ.6) THEN
              FACQQ1=FACQQ1+2.25D0*SQMAVG*(THQ-UHQ)/(SH*THQ)*REDHGS+
     &        2.25D0*THQ*UHQ/SH2*SQDHGS
              FACQQ2=FACQQ2+2.25D0*SQMAVG*(UHQ-THQ)/(SH*UHQ)*REDHGS+
     &        2.25D0*THQ*UHQ/SH2*SQDHGS
            ELSE
              FACQQ1=FACQQ1+2.25D0*SQMAVG*(THQ-UHQ)/(SH*THQ)*REDLGS+
     &        2.25D0*THQ*UHQ/SH2*SQDLGS
              FACQQ2=FACQQ2+2.25D0*SQMAVG*(UHQ-THQ)/(SH*UHQ)*REDLGS+
     &        2.25D0*THQ*UHQ/SH2*SQDLGS
            ENDIF
          ENDIF
          FACQQ1=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ1
          FACQQ2=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ2
          IF(MSTP(35).GE.1) THEN
            FATRE=PYHFTH(SH,SQMAVG,2D0/7D0)
            FACQQ1=FACQQ1*FATRE
            FACQQ2=FACQQ2*FATRE
          ENDIF
          WID2=1D0
          IF(MINT(55).EQ.6) WID2=WIDS(6,1)
          IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=WIDS(MINT(55),1)
          FACQQ1=FACQQ1*WID2
          FACQQ2=FACQQ2*WID2
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 530
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ1
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=2
          SIGH(NCHN)=FACQQ2
  530     CONTINUE
        ENDIF
      ENDIF
 
CMRENNA--
 
      RETURN
      END
