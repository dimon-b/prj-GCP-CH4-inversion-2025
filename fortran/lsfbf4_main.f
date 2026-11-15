C--------------------------------------------------------------------
C
C  Program : lsfbf4-4.f
C
C          KD must be sufficiently larger than the value of ((910
C       +2920)/(mean data interval expressed in day)+data number of 
C       the record). KO must be sufficiently larger than the value 
C       of (910+2920+length of the data record expressed in day).  
C       If you want to change KD, you must change KD of subroutine 
C       or function Fourier, Frc1, Reinsch, Setup, Chol1d, Smooth,
C       Extrad1 and Judge1,too.  If you want to change KO, you must 
C       change N or KO in subroutines Bwf, Scycle, Yac and Judge1,
C       too.
C
C                                 By T. Nakazawa (Mar. 8, 1991)
C
C  ------------------------------------------------------------------
      subroutine lsfbf4(ndata0,yearin,monthin,dayin,yin,tim,
     1           yout,ytrn,ndata3,OUTVALUE)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(KD=30000,KO=60000,KL=30000,
     1          NCFS=26,MCFS=16,ITESS=10)
C
      INTEGER YEARIN,MONTHIN,DAYIN
      DIMENSION A(20),T(KD),Y(KD),TR(KD),SD(KD),atemp(20)
      DIMENSION A1(20),B1(20),C1(20)
      DIMENSION BS(KO),X2(KO),S(KO),SS(KO)
      DIMENSION YT(KO),YI(KO),YS(KO),CDEV(KO),YM(KO),yit(KO)
      DIMENSION tim(KO),yout(KO),ytrn(KO)
      DIMENSION BS1(KD),SIG(KD),SPQF(KD),SPQC(KD,3)
      DIMENSION TR1(KD),TF1(KD),X1(KO)
      DIMENSION W1(KO),W2(KO),BS2(KD)
      DIMENSION SMAX(50),SMIN(50),MAXD(50),MIND(50)
      DIMENSION TRB(KD),MDI(KD),MDF(KD)
      DIMENSION IT(2),A2(20),B2(20),C2(20),W3(KO),TYD(KD)
      DIMENSION WVAR(KD)
      DIMENSION TBL(KD),YBL(KD),TYDBL(KD),TNBL(KD),YNBL(KD),TYDNBL(KD)
      DIMENSION T0(KD),Y0(KD),TYD0(KD)
      DIMENSION STYDOUT(KL),SYOUT(KL),DETR(KD),DESE(KD),VARI(KD)
      DIMENSION TYDIN(KD),TIN(KD),YIN(KD)
      DIMENSION YEARIN(KD),MONTHIN(KD),DAYIN(KD),JULIN(KD)
      DIMENSION OUTTYD(KL),OUTY(KL)
C 
      CHARACTER INFILE*80,OUTFILE(10)*80,LABEL(14)*10,YACLABEL*80
      CHARACTER INDATAM*1,IBWF*1,IFILE*1
      CHARACTER STATION*8,CLASS*5,INUM*1
!     logical INUM
C
      DATA LABEL /'sin(t)','con(t)','sin(t/2)','con(t/2)',
     1     'sin(t/3)','con(t/3)','sin(t/4)','con(t/4)','sin(t/5)',
     2     'con(t/5)','sin(t/6)','con(t/6)','const.','t'/
C
      COMMON /SYEAR/IYS
      COMMON /DAYYEAR/DYEAR,DYEAR0,DYEAR1,DYEAR2
!
!     DATA NCF/26/
!     DATA MCF/16/
!     DATA ITES/10/
      NCF = NCFS
      MCF = MCFS
      ITES = ITESS
C
      PI=3.141592653589793D0
C
      DYEAR0 = 365.25D0
      DYEAR1 = 365.0D0
      DYEAR2 = 366.0D0
C 
      ITYPE = 3    !data format  : year, mon, day, data
C
      INDATA = KD  !for unlimited rows
C
      INUM = 'y'   ! Y if there are missing values
!     INUM = .T.   ! Y if there are missing values
C
      STATION = 'test2222' !for writing to a data file
C
      CUTL = 24    !cutoff length for digital filter
      NH  = 6      !no . of harmonics  ( permitted value : 0-6 )
      M = 1        !every how many days output is requested
C
C ---- naming of output files --------
C
      OUTFILE(1) = STATION//'.lst'
      OUTFILE(2) = STATION//'.fit'
      OUTFILE(3) = STATION//'.bal'
      OUTFILE(4) = STATION//'.out'
      OUTFILE(5) = STATION//'.dev'
      OUTFILE(6) = STATION//'.rej'
      OUTFILE(7) = STATION//'.bwf'
C
      IFILE = 'n'  !we output only first two files; type yes for four more
C
      IBWF = 'n'   !i am quite sure of this
C
C 
!     OPEN(UNIT=1,FILE=INFILE,STATUS='OLD')
      OPEN(UNIT=2,FILE=OUTFILE(1))
C
      CALL MDINPUT1(YEARIN,MONTHIN,DAYIN,TYDIN,JULIN,NDATA0)
C
      IYS = INT(TYDIN(1))
      AIYS = 0.0
      IF(MOD(IYS,4).EQ.0) AIYS=1.0
      DO 35 I=1,NDATA0
        DYEAR = DYEAR1
        IF(MOD(INT(TYDIN(I)),4).EQ.0) DYEAR = DYEAR2
        JDAY = JULIN(I)-1
        NYEAR= INT(TYDIN(I))-IYS
        LYEAR= INT(TYDIN(I))/4-IYS/4
        IF(MOD(INT(TYDIN(I)),4).EQ.0) LYEAR=LYEAR-1
        YDAY = NYEAR*DYEAR1 + LYEAR*1.0
        IF(NYEAR.NE.0) YDAY = YDAY + AIYS
        TIN(I) = JDAY + YDAY
!       print *, YEARIN(i),MONTHIN(i),DAYIN(i),TYDIN(i),JULIN(i)
   35 CONTINUE
      print*, 'Conversion completed YEARIN,MONTHIN,DAYIN,TYDIN,JULIN'
C
      IYS = INT(TYDIN(1))
C
      I2 = 0
      J  = 0
      IF(INUM.EQ.'y') THEN               !this works well, i have commented it 
!     IF (INUM) THEN               !this works well, i have commented it 
        DO 36 I = 1,NDATA0               !for the purpose of simplicity; i 
          IF(YIN(I).EQ.OUTVALUE) THEN    !always use INUM = 'n'
            J = J+1
            OUTTYD(J) = TYDIN(I)
            OUTY(J)   = YIN(I)      
          ELSE
            I2 = I2+1
            TYD(I2) = TYDIN(I)
            T(I2)   = TIN(I)
            Y(I2)   = YIN(I)
          ENDIF
   36   CONTINUE
      ELSE
        DO 37 I = 1,NDATA0
          I2 = I2+1
          TYD(I2) = TYDIN(I)
          T(I2)   = TIN(I)
          Y(I2)   = YIN(I)
   37   CONTINUE
      ENDIF
C
      NDATA    = I2
      NDATAOUT = J
C
      CLOSE(1)
C 
      WRITE(2,*) ' Number of data ',NDATA+NDATAOUT
      WRITE(*,*) ' Number of data ',NDATA+NDATAOUT
      WRITE(2,*) ' Missing data   ',NDATAOUT
      WRITE(*,*) ' Missing data   ',NDATAOUT
      IF(NDATAOUT.GT.0) THEN                       !commented for simplicity
        DO 38 I=1,NDATAOUT
          WRITE(2,*) OUTTYD(I),OUTY(I)
!         WRITE(*,*) OUTTYD(I),OUTY(I)
   38   CONTINUE
      ENDIF
C
  137 WRITE(2,*) 'Data for fitting ',NDATA
      WRITE(*,*) 'Data for fitting ',NDATA
      WRITE(2,610) (TYD(I),I=1,NDATA)        !actual time 1997.****
!     WRITE(*,610) (TYD(I),I=1,NDATA)
  610 FORMAT(7F10.3)
C
      WRITE(2,611) (INT(T(I)),I=1,NDATA)     !some sort of day number 0-ndata
!     WRITE(*,611) (INT(T(I)),I=1,NDATA)
  611 FORMAT(7I10)
C 
      DO 40 I=1,NDATA-1
        IF(T(I).GE.T(I+1))THEN
          WRITE(*,*) I,T(I),TYD(I)
          WRITE(*,*) I,T(I+1),TYD(I+1)
          !PAUSE 'Check the time record in input data'
        ENDIF
   40 CONTINUE
C   
      WRITE(2,615) (Y(I),I=1,NDATA)         !CO2 data to be fitted
!     WRITE(*,615) (Y(I),I=1,NDATA)
  615 FORMAT(7F10.3)
C 
      EXD1 = 910.D0
      EXD2 = 2920.D0
      DO 4500 I=1,KD
        SIG(I) = 1.D0
 4500 CONTINUE
C 
C Extension of the data record
C
      I=0
   50 I=I+1
      X2(I) = T(1)-EXD1+1.D0*DFLOAT(I-1)     !x2 is set here for the first time - pkp
      X1(I) = X2(I)                          ! = day# - 910 + i
      IF(X2(I).GT.T(NDATA)+EXD2) GO TO 55
      GO TO 50
   55 NN=I-1
      DO 230 I=1,NN
        IF(DABS(X2(I)-T(1)).LE.1.) MI=I
        IF(DABS(X2(I)-T(NDATA)).LE.1.) MF=I
        X2(I)=X2(I)/DYEAR0                   ! DYEAR0 = 365.25
  230 CONTINUE
      DO 90 I=1,NDATA
        T(I) = T(I)/DYEAR0
   90 CONTINUE
C
C
C   SAVE THE ORIGINAL DATA SET
C
      NDATA0 = NDATA
      DO 91 I=1,NDATA
        T0(I)   = T(I)
        Y0(I)   = Y(I)
        TYD0(I) = TYD(I)
   91 CONTINUE
C
C
C ====  Reject outlier ===========================================
C
      WRITE(*,701) 
      WRITE(2,701)
  701 FORMAT(/1H 72(1H*),
     &       /' Reject data with more than 3 sigma variance')

C
C Fit of linear eq. + Fourier to the data
C
      WRITE(*,*) 'Fit of linear eq. + Fourier to the data'
      WRITE(2,*) 'Fit of linear eq. + Fourier to the data'
      NH       = 2*NH
      NDATAOUT = 0
      ITEOUT   = 0
  174 ITEOUT   = ITEOUT+1               !ITERATION and data rejection starts....
      NTERM    = NH+2
!     WRITE(*,702) ITEOUT
!     WRITE(2,702) ITEOUT
! 702 FORMAT(/1H ,' ITERATION = ',I4)
      CALL FRC1(NTERM,NH,NDATA,T,Y,A)      ! A is probably calculated here - pkp
      !print *, 'A 1= ', A
      do ii = 1, 20
         atemp(ii) = a(ii)
      end do
      SUM=0.D0
      DO 160 I=1,NDATA
        YX    = SEASON1(A,PI,T(I),NH)
        TR(I) = Y(I)-YX
        SUM   = SUM+(TR(I)-A(NH+1)-A(NH+2)*T(I))**2
  160 CONTINUE
      SIGM = DSQRT(SUM/DFLOAT(NDATA))
C      
      WRITE(2,700)
  700 FORMAT(/' COEFFICIENTS A(K)= ')
      WRITE(2,710) LABEL(13),A(NH+1)
      WRITE(2,710) LABEL(14),A(NH+2)
      DO 161 K=1,NTERM
        WRITE(2,710) LABEL(K),A(K)
  161 CONTINUE
  710 FORMAT(1H ,A10,E14.6)
C
      if ( SIGM .le. 0.001 ) SIGM = 0.01
      WRITE(2,640) SIGM
      WRITE(*,640) SIGM
  640 FORMAT(1H ,7HST-DEV=,E14.6)   !FORMAT(/1H ,7HST-DEV=,E14.6/)  - pkp
 
      DO 171 I=1,NN
        YI(I) = A(NH+1)+A(NH+2)*X2(I)      ! YI is calculated here - pkp
        yit(i) = yi(i)
        !print*,'YI(I),A(NH+1),X2(I)',i,YI(I),A(NH+1),X2(I),A(NH)
  171 CONTINUE
      YACLABEL= 'From Linear Trend '
      CALL YAC(X1,YI,NN,MI,MF,YACLABEL)
C
C
C Fit of Reinsch + Fourier to the data
C
      T5    = DYEAR0*5.0/DYEAR0
      DD    = (T(NDATA)-T(1))/DFLOAT(NDATA)
      RM    = (T5/(2.D0*PI))**4/DD
      NTERM = NH
      !print *, 'A 2 (before calling FRC1= ', A
C
      ITE = 0
  175 ITE = ITE+1
      SM  = SMOOTH(T,TR,SIG,NDATA,RM,BS1)
      CALL REINSCH(T,BS1,SIG,SPQC,SPQF,NDATA,0.D0)
      LI=1
      DO 176 I=1,NDATA
        TR1(I) = SPLINT(NDATA,T,SPQF,SPQC,KD,LI,T(I))
        SD(I)  = Y(I)-TR1(I)
  176 CONTINUE
      CALL FRC1(NTERM,NH,NDATA,T,SD,A)
      SUM=0.D0
      DO 177 I=1,NDATA
        TR(I) = Y(I)-SEASON1(A,PI,T(I),NH)
        SUM   = SUM+(TR(I)-TR1(I))**2
  177 CONTINUE
      SIGM = DSQRT(SUM/DFLOAT(NDATA))
      IF(ITE.EQ.1) THEN
        SIGMW = SIGM
      ELSE
        IF(1.D0-SIGM/SIGMW.LT.0.001D0) GO TO 179
        SIGMW = SIGM
      ENDIF
      GO TO 175
C
  179 SIGM3 = SIGM*3.0 
      !print *, 'A 2 (after calling FRC1= ', A
      J     = 0
      J2    = 0
      DO 120 I=1,NDATA
        WVAR(I) = TR1(I) + SEASON1(A,PI,T(I),NH)
        VAR     = Y(I)-WVAR(I)
        IF(DABS(VAR).LE.SIGM3) THEN
          J = J+1
          TBL(J)   = T(I)
          YBL(J)   = Y(I)
          TYDBL(J) = TYD(I)
        ELSE
          J2 = J2+1
          TNBL(J2)   = T(I)
          YNBL(J2)   = Y(I)
          TYDNBL(J2) = TYD(I)
        ENDIF
  120 CONTINUE
      IF(J2.EQ.0) GOTO 178
      NDATA  = J
      NDATA2 = J2
      DO 121 I=1,NDATA
        T(I)   = TBL(I)
        Y(I)   = YBL(I)
        TYD(I) = TYDBL(I)
  121 CONTINUE
!     WRITE(*,703)
!     WRITE(2,703) 
      DO 122 I=1,NDATA2
!       WRITE(*,704) TYDNBL(I),TNBL(I),YNBL(I)
!       WRITE(2,704) TYDNBL(I),TNBL(I),YNBL(I)
        NDATAOUT = NDATAOUT + 1
        STYDOUT(NDATAOUT) = TYDNBL(I)
        SYOUT(NDATAOUT)   = YNBL(I)
  122 CONTINUE
  703 FORMAT(/1H ,' Rejected data')
  704 FORMAT(1H ,3(F12.5))
      GOTO 174                        !ITERATION and data rejection ends ....
C
C  
C 
  178 WRITE(*,705) NDATAOUT
      WRITE(2,705) NDATAOUT
  705 FORMAT(/1H ,'Number of Rejected data = ',I5)
      IF(IFILE.EQ.'y') THEN                        !we set ifile = n
        IF(NDATAOUT.NE.0) THEN
          CALL SORT1(STYDOUT,SYOUT,NDATAOUT,OUTFILE(6))
        ENDIF
      ENDIF
C         
C   
C         
C ** Interpolation of missing data and extrapolation            **
C **  the data to be added to both ends of the record            **
C
      SM=SMOOTH(T,TR,SIG,NDATA,RM,BS1)
      CALL REINSCH(T,BS1,SIG,SPQC,SPQF,NDATA,0.D0)
      SPT=30.D0/DYEAR0
      NDATA1=0
C
C --  Left side of the record  -----------------------------------
C
      MD=1
      CALL EXTRAD1(TF1,X2(MI)-SPT,X2(1),DD,IP)
      LI=1
      DO 220 I=1,IP
        TR1(I)=SPLINT(NDATA,T,SPQF,SPQC,KD,LI,TF1(I))
  220 CONTINUE
C
      MDI(MD) = NDATA1+1
      MDF(MD) = MDI(MD)+IP-1
      NDATA1  = NDATA1+IP
      DM      = DYEAR0/2.0
      DO 221 I=1,NDATA-1
C
        NDATA1=NDATA1+1
        TF1(NDATA1) = T(I)
        TR1(NDATA1) = TR(I)
        IF(T(I+1)-T(I).LT.DM) GO TO 221
C
C -- Interpolation of missing data   -----------------------------
C
        MD = MD+1
        CALL EXTRAD1(W1,T(I),T(I+1),DD,IP)
        LI = 1
        DO 222 J=1,IP
          W2(J) = SPLINT(NDATA,T,SPQF,SPQC,KD,LI,W1(J))
  222   CONTINUE
C
        MDI(MD) = NDATA1+1
        MDF(MD) = MDI(MD)+IP-1
        DO 223 J=1,IP
          NDATA1      = NDATA1+1
          TF1(NDATA1) = W1(J)
          TR1(NDATA1) = W2(J)
  223   CONTINUE
C
  221 CONTINUE
C
      NDATA1      = NDATA1+1
      TF1(NDATA1) = T(NDATA)
      TR1(NDATA1) = TR(NDATA)
C
C --  Right side of the record   ---------------------------------
C
      MD = MD+1
      CALL EXTRAD1(W1,X2(MF)+SPT,X2(NN),DD,IP)
      LI = 1
      DO 224 I=1,IP
        W2(I) = SPLINT(NDATA,T,SPQF,SPQC,KD,LI,W1(I))
  224 CONTINUE
C
      MDI(MD) = NDATA1+1
      MDF(MD) = MDI(MD)+IP-1
      DO 225 J=1,IP
        NDATA1      = NDATA1+1
        TF1(NDATA1) = W1(J)
        TR1(NDATA1) = W2(J)
  225 CONTINUE
C 
C -- Remove of the approximate trend from the data  ---------------
C
      LI = 1
      DO 312 I=1,NDATA
        BS1(I) = SPLINT(NDATA,T,SPQF,SPQC,KD,LI,T(I))
        Y(I)   = Y(I)-BS1(I)
  312 CONTINUE
      LI = 1
      DO 313 I=1,NDATA1
        BS2(I) = SPLINT(NDATA,T,SPQF,SPQC,KD,LI,TF1(I))
        TR1(I) = TR1(I)-BS2(I)
  313 CONTINUE
      LI = 1
      DO 314 I=1,NN
        BS(I) = SPLINT(NDATA,T,SPQF,SPQC,KD,LI,X2(I))
  314 CONTINUE
C
C 
C =====   Start of the data fitting  ================================
C
      WRITE(*,669) 
      WRITE(2,669)
  669 FORMAT(/1H 72(1H*),/' Start of the data fitting ')  !fitting '/) - pkp
C
C
C ---  coefficients of digital filters  ---------------------------
C
      NCF = NCF/2
C
      IF(IBWF.EQ.'y') THEN
        OPEN(UNIT=10,FILE=OUTFILE(7))
        WRITE(10,60) 
   60   FORMAT('  Butterworth filter for long-term variation'/)
      ENDIF
C      
      CALL BWFCOEF(A1,B1,C1,NCF,CUTL,IBWF)
C
      CUT1= 0
      IF(NH.NE.0) CUT1= DBLE(NH/2)
      CUTS= 12.0D0/(CUT1+1.0D0)
      MCF = MCF/2
      IF(IBWF.EQ.'y') WRITE(10,65)
   65 FORMAT(///'  Butterworth filter for short-term variation'/)  
      CALL BWFCOEF(A2,B2,C2,MCF,CUTS,IBWF)
      IF(IBWF.EQ.'y') CLOSE(10)
C
      IJF   = 2
      NH    = NH+2
      NTERM = NH
C
C -------------------------------------------------------------------
C
      DO 256 IJ=1,IJF
C 
        ISTP  = 0
        SIGM1 = SIGM
C 
        ITE = 0
  210   ITE = ITE+1
C 
!       WRITE(*,670) IJ,ITE
        WRITE(2,670) IJ,ITE
  670   FORMAT('  ITERATION=',I2,1H-,I2) !FORMAT(/'  ITERATION=',I2,1H-,I2) - pkp
 
        LI = 1
        DO 180 I=2,NDATA1
          PP = (TR1(I)-TR1(I-1))/(TF1(I)-TF1(I-1))
          QQ = TR1(I)-PP*TF1(I)
          DO 181 J=LI,NN
            SS(J) = PP*X2(J)+QQ
            IF(DABS(X2(J)-TF1(I)).LT.0.002D0) GO TO 182
  181     CONTINUE
  182     LI=J+1
  180   CONTINUE
C
        IF(IJ.EQ.1) THEN
          CALL BWF(SS,YM,NN,A1,B1,C1,NCF)
          LI = 1
          DO 195 I=1,NDATA
            DO 196 J=LI,NN
              IF(DABS(X2(J)-T(I)).GT.0.002D0) GO TO 196
              SD(I)  = Y(I)-YM(J)
              TRB(I) = YM(J)
  196       CONTINUE
  195     CONTINUE
          CALL FRC1(NTERM,NH,NDATA,T,SD,A)
      !print *, 'A 3 (after calling FRC1= ', A
C
C ----    smoothing of seasonal cycle with Butterworth filter ----
C
          DO 281 I=1,NN
            W2(I) = SEASON1(A,PI,X2(I),NH)
  281     CONTINUE
      !print '(15f8.2)',(W2(i), i = 1, n)
      !print '(15f8.2)',(A2(i), i = 1, 20)
      !print '(15f8.2)',(b2(i), i = 1, 20)
      !print '(15f8.2)',(c2(i), i = 1, 20)
          CALL BWF(W2,W3,NN,A2,B2,C2,MCF)
          N=0
          DO 280 I=MI,MF,3       !mi = 912, mf=2006, nterm=nh=8, n=365
            N = N+1
            W1(N) = X2(I)
            W2(N) = W3(I)                       ! A is becoming zero here!!!
  280     CONTINUE
          CALL FRC1(NTERM,NH,N,W1,W2,A)
      !print*,'A 4(after calling FRC1=',NN,MCF,Mi,MF,NTERM,N,NH,A
      !print '(15f8.2)',(W1(i), i = 1, n)
      !print '(15f8.2)',(W2(i), i = 1, n)
      !print '(15f8.2)',(W3(i), i = 1, n)
C
C ---   end of smoothing -----
C
          SUM = 0.D0
          LI  = 1
          DO 190 I=1,NDATA
            TR(I) = Y(I)-SEASON1(A,PI,T(I),NH)
            SUM   = SUM+(TR(I)-TRB(I))**2
            DO 191 J=LI,NDATA1
              IF(DABS(TF1(J)-T(I)).GT.0.002D0) GO TO 191
              TR1(J)=TR(I)
              LI = J
  191       CONTINUE
  190     CONTINUE
          SIGM = DSQRT(SUM/DFLOAT(NDATA))
          CALL JUDGE(ISTP,SIGM1,SIGM,YM,MDI,MDF,MD,X2,TF1,TR1,NN)
          IF(ISTP.EQ.1.OR.ITE.EQ.ITES) GO TO 258
        ELSE
C 
          CALL BWF(SS,YI,NN,A2,B2,C2,MCF)
          SUM = 0.D0
          LI  = 1
          DO 250 I=1,NDATA
            DO 251 J=LI,NN
              IF(DABS(X2(J)-T(I)).GT.0.002D0) GO TO 251
              SUM = SUM+(TR(I)-YI(J))**2
              LI  = J
  251       CONTINUE
  250     CONTINUE
          SIGM = DSQRT(SUM/DFLOAT(NDATA))
          CALL JUDGE(ISTP,SIGM1,SIGM,YI,MDI,MDF,MD,X2,TF1,TR1,NN)
          IF(ISTP.EQ.1.OR.ITE.EQ.ITES) GO TO 257
        END IF
 
!       WRITE(*,600) SIGM
        WRITE(2,600) SIGM
  600   FORMAT('  SIGM=',E14.6)  !FORMAT(/'  SIGM=',E14.6) - pkp
        GO TO 210
 
  258   LI = 1
        DO 270 J=1,NDATA1
          DO 271 K=LI,NN
            IF(DABS(X2(K)-TF1(J)).GT.0.002D0) GO TO 271
            TR1(J) = TR1(J)-YM(K)
            LI     = K
  271     CONTINUE
  270   CONTINUE
        LI = 1
        DO 272 J=1,NDATA
          DO 273 K=LI,NDATA1
            IF(DABS(TF1(K)-T(J)).GT.0.002D0) GO TO 273
            TR(J) = TR1(K)
            LI    = K
  273     CONTINUE
  272   CONTINUE
  257   WRITE(*,600) SIGM
        WRITE(2,600) SIGM
        IT(IJ) = ITE
 
  256 CONTINUE
C 
C  ----  Recover of the approximate trend removed above  -----
C
      DO 200 I=1,NN
        YM(I) = YM(I)+BS(I)
        YS(I) = SEASON1(A,PI,X2(I),NH)
        YT(I) = YI(I)+YS(I)+YM(I)     !problem here - ys and yi are zero, ym okay - pkp
        S(I)  = YT(I)-YM(I)           !Basically A and YI are getting set to 0 somewhere
        !print '(2i7,24f8.2)', i,NH,Atemp,X2(I),YS(I),yit(i)
  200 CONTINUE
 
      DO 709 I=2,NN
        CDEV(I) = (YM(I)-YM(I-1))/(1.D0/DYEAR0)
  709 CONTINUE
C
C
C ===========   Output of the results  =========================
C
      WRITE(*,668)
      WRITE(2,668)
  668 FORMAT(/1H 72(1H*),/' Results'/)
C
C
      WRITE(2,671) IT(1),IT(2)
  671 FORMAT(1H ,'No. of ITERATION(Long- and short-term)',I2,1H-,I2)  
      WRITE(2,700)
      DO 162 K=1,NTERM
        WRITE(2,710) LABEL(K),A(K)
  162 CONTINUE
C
      WRITE(2,600) SIGM
C
      YACLABEL = 'From the Fitted Curve'
      CALL YAC(X1,YT,NN,MI,MF,YACLABEL)
      YACLABEL = 'From the Trend Curve '
      CALL YAC(X1,YM,NN,MI,MF,YACLABEL)
C
C
C  ---  Seasonal Amplitude  ------------------------------------
C
      CALL SCYCLE(YS,X1,NN,MI,MI+364,KI,KF,SMAX,SMIN,MAXD,MIND)
      WRITE(2,650)
  650 FORMAT(/' * AVERAGE SEASONAL CYCLE')
      WRITE(2,651)
  651 FORMAT(/5X,15X,'MAX',14X,'MIN',14X,'AMP')
      WRITE(2,652) SMAX(KI),MAXD(KI),SMIN(KI),MIND(KI),
     &             SMAX(KI)+DABS(SMIN(KI))
  652 FORMAT(17X,F7.3,1H(,I4,1H),5X,F7.3,1H(,I4,1H),5X,F7.3)
C 
      CALL SCYCLE2(S,X1,NN,MI,MF,KI,KF,SMAX,SMIN,MAXD,MIND)
      WRITE(2,653)
  653 FORMAT(/' * YEAR-TO-YEAR SEASONAL CYCLE')
      WRITE(2,651)
      WRITE(2,654)(IYS+I-KI,SMAX(I),MAXD(I),SMIN(I),MIND(I),
     &                      SMAX(KI)+DABS(SMIN(I)),I=KI,KF)
  654 FORMAT(6X,I6,5X,F7.3,1H(,I4,1H),5X,F7.3,1H(,I4,1H),5X,F7.3)
C
C
C ---   deviation from the curve  -------------------------------
C
      LI=MI-5
      DO 261 I=1,NDATA
        DO 262 J=LI,MF+5
          IF(ABS(X2(J)-T(I)).GT.0.002D0) GO TO 262
          TRB(I) = Y(I)+BS1(I)-YT(J)
          SS(I)  = YT(J)
          W2(I)  = YS(J)+YM(J)
          W1(I)  = Y(I)+BS1(I)-W2(I)
          LI=J
  262   CONTINUE
  261 CONTINUE
C 
C      WRITE(2,751)
C  751 FORMAT(/' * DATA DEVIATIONS FROM THE FITTED CURVE(=T+S+I)'/)
C      WRITE(2,752)(I,TRB(I),INT(T(I)*DYEAR0+.5),SS(I),I=1,NDATA)
C  752 FORMAT(2(I5,F8.3,1H(,I5,F8.3,1H)))
C      WRITE(2,753)
C  753 FORMAT(/' * DATA DEVIATIONS FROM THE FITTED CURVE(=T+S)'/)
C      WRITE(2,752)(I,W1(I),INT(T(I)*DYEAR0+.5),W2(I),I=1,NDATA)
C
      IF(IFILE.EQ.'y') THEN
        OPEN(UNIT=12,FILE=OUTFILE(3))
        OPEN(UNIT=13,FILE=OUTFILE(4))
        OPEN(UNIT=14,FILE=OUTFILE(5))
      ENDIF
C
      LI = MI-5
      SUM = 0.0D0
      DO 126 I=1,NDATA0
        DO 127 J=LI,MF+5
          IF(DABS(X2(J)-T0(I)).GT.0.002D0) GOTO 127
          DESE(I) = Y0(I)-(YS(J)+YI(J))
          DETR(I) = Y0(I)-YM(J)
          VARI(I) = Y0(I)-YT(J)
          SUM     = SUM + VARI(I)**2
          LI      = J
  127   CONTINUE
  126 CONTINUE
      SIGM  = DSQRT(SUM/NDATA0)
      SIGM3 = SIGM*3.0
      NVARI1 = 0
      NVARI2 = 0
      DO 128 I=1,NDATA0
        IF(DABS(VARI(I)).LT.SIGM3) THEN
          NVARI1 = NVARI1 + 1
          IF(IFILE.EQ.'y') THEN
            WRITE(12,100) I,TYD0(I),Y0(I),VARI(I)
          ENDIF
        ELSE
          NVARI2 = NVARI2 + 1
          IF(IFILE.EQ.'y') THEN
            WRITE(13,100) I,TYD0(I),Y0(I),VARI(I)
          ENDIF
        ENDIF
        IF(IFILE.EQ.'y') THEN
          WRITE(14,110) I,TYD0(I),INT(T(I)*DYEAR0+.5),
     &                  Y0(I),DESE(I),DETR(I),VARI(I)
        ENDIF
  128 CONTINUE
  100 FORMAT(1H ,I5,2X,3(F12.5))
  110 FORMAT(1H ,I5,2X,F12.5,I7,4(F12.5))
      IF(IFILE.EQ.'y') THEN
        CLOSE(12)
        CLOSE(13)
        CLOSE(14)
      ENDIF
C
C
C -----------------------------------------------------------------
      MI2 = MI-30 
      MF2 = MF+30
C
      OPEN(UNIT=1,FILE=OUTFILE(2))
      STYEAR = DBLE(IYS)
      ndata3 = 0
      yout = 0.0d0
      ytrn = 0.0d0
      tim  = 0.0d0
      DO 124 I=MI2,MF2,M
        TIME     = X1(I)/DYEAR0+STYEAR
        CURVE    = YT(I)
        TREND    = YM(I)
        SEASN    = YS(I)+YI(I)
        ONHM     = YS(I)
        OUTHM    = YI(I)
        RATE     = CDEV(I) 
        ndata3 = ndata3 + 1
        yout (ndata3) = YT(I)
        ytrn (ndata3) = YM(I)
        tim  (ndata3) = time
        WRITE(1,125) TIME,CURVE,TREND,SEASN,ONHM,OUTHM,RATE
!       WRITE(*,125) TIME,CURVE,TREND,SEASN,ONHM,OUTHM,RATE
  124 CONTINUE
      CLOSE(1)
  125 FORMAT(1H ,7(F11.5))
C
C  ------------------------------------------------------------------
C
  129 WRITE(2,750) INT(DD*DYEAR0)
      WRITE(*,750) INT(DD*DYEAR0)
  750 FORMAT(/'  Mean Date Interval (days)   :',I6)
      WRITE(2,754) NVARI1
      WRITE(*,754) NVARI1
      WRITE(2,755) NVARI2
      WRITE(*,755) NVARI2
  754 FORMAT(1H ,' Number of data (<  3 sigma) :',I6)
  755 FORMAT(1H ,' Number of data (>= 3 sigma) :',I6)
!     WRITE(2,756) (OUTFILE(K),K=1,2)
!     WRITE(*,756) (OUTFILE(K),K=1,2)
! 756 FORMAT(/'  The output files from this curve fitting are :'/
!    &4X,A15/4X,A15)
      IF(IFILE.EQ.'y') THEN
        WRITE(2,757) (OUTFILE(K),K=3,6)
        WRITE(*,757) (OUTFILE(K),K=3,6)
      ENDIF
      IF(IBWF.EQ.'y') THEN
        WRITE(2,758) OUTFILE(7)
        WRITE(*,758) OUTFILE(7)
      ENDIF
  757 FORMAT(3(4X,A15/),4X,A15)
  758 FORMAT(4X,A15)
C
      CLOSE(2)
C 
      return !STOP
      END
