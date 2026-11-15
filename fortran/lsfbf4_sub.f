C
C***********   Subroutines &  Functions *****************************
C
C
      SUBROUTINE FRC1(NTERM,NH,NDATA,T,Y,A)
C--------------------------------------------------------------------
C
C This subroutine program is used to estimate the coefficients of
C Fourier harmonics and is called in LSFCO2 or LSFFR6.
C Subroutine: Fourier, Invmat
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KD=30000)
      DIMENSION T(KD),A(20),F(20,KD),Y(KD)
C
      PI=3.141592653589793D0
C 
      DO 90 I=1,20
        A(I) = 0.D0
   90 CONTINUE
 
      IF(NH.EQ.0) GO TO 80
      DO 100 I=1,NH
        DO 100 K=1,NDATA
          IF(MOD(I,2).EQ.1) F(I,K)=DSIN(2.D0*DBLE((I-1)/2+1)*PI*T(K))
          IF(MOD(I,2).EQ.0) F(I,K)=DCOS(2.D0*DBLE((I-1)/2+1)*PI*T(K))
  100 CONTINUE
C
   80 IF(NH.EQ.0.AND.NH.EQ.NTERM) GO TO 120
      IF(NTERM.EQ.NH) GO TO 200
      DO 110 I=NH+1,NH+2
        DO 110 K=1,NDATA
          IF(I.EQ.NH+1) F(I,K) = 1.
          IF(I.EQ.NH+2) F(I,K) = T(K)
  110 CONTINUE
  200 CALL FOURIER(F,Y,NTERM,NDATA,A)
 
  120 RETURN
      END
C
C
C
      SUBROUTINE FOURIER(F,Y,NOT,NOD,A)
C----------------------------------------------------------------
C
C LEAST SQUARE FITTING WITH FOURIER HARMONICS
C                                             
C----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(KD=30000) 
      DIMENSION F(20,KD),A(20),B(20),C(20,20),D(20,20),Y(KD)
C
C GENERATE NORMAL MATRIX C
C
      DO 110 I=1,NOT
        DO 110 J=1,NOT
          SUM = 0.D0
          DO 120 K=1,NOD
            SUM = SUM+F(I,K)*F(J,K)
  120     CONTINUE
          C(I,J) = SUM
  110 CONTINUE
C
C GENERATE RIGHT-SIDE MATRIX B
C
      DO 130 I=1,NOT
        SUM = 0.D0
        DO 140 K=1,NOD
          SUM = SUM+F(I,K)*Y(K)
  140   CONTINUE
        B(I) = SUM
  130 CONTINUE
C
C INVERT NORMAL MATRIX C
C
      CALL INVMAT(C,D,NOT)
C
C MATRIX D IS INVERSE OF NORMAL MATRIX
C FORM MATRIX PRODUCT OF INVERSE AND RIGHT-SIDE MATRIX
C
      DO 150 I=1,NOT
        A(I) = 0.D0
        DO 150 J=1,NOT
          A(I) = A(I)+D(I,J)*B(J)
  150 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE INVMAT(A,B,N)
C----------------------------------------------------------------
C
C MATRIX INVERSION BY ELIMINATION WITH PARTIAL PIVITING
C ORIGINAL MATRIX=A, INVERSE MATRIX=B                                          
C
C----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(20,20),B(20,20)
C
      ESP=0.0000001D0
C
C CONSTRUCT IDENTITY MATRIX B(I,J)=I
C
      DO 6 I=1,N
        DO 5 J=1,N
          IF(I-J) 4,3,4
    3     B(I,J)=1.D0
          GO TO 5
    4     B(I,J)=0.D0
    5   CONTINUE
    6 CONTINUE
C
C LOCATE MAXIMUM MAGNITUDE A(I,K) ON OR BELOW MAIN DIAGONAL
C
      DEL = 1.D0
      DO 45 K=1,N
        IF(K-N) 12,30,30
   12   IMAX = K
        AMAX = DABS(A(K,K))
        KP1  = K+1
        DO 20 I=KP1,N
          IF(AMAX-DABS(A(I,K))) 15,20,20
   15     IMAX = I
          AMAX = DABS(A(I,K))
   20   CONTINUE
C
C INTERCHANGE ROWS IMAX AND K IF IMAX NOT EQUAL TO K
C
        IF(IMAX-K) 25,30,25
   25   DO 29 J=1,N
          ATMP      = A(IMAX,J)
          A(IMAX,J) = A(K,J)
          A(K,J)    = ATMP
          BTMP      = B(IMAX,J)
          B(IMAX,J) = B(K,J)
          B(K,J)    = BTMP
   29   CONTINUE
        DEL=-DEL
   30   CONTINUE
C
C TEST FOR SINGULAR MATRIX
C
        IF(DABS(A(K,K))-ESP) 93,93,35
   35   DEL = A(K,K)*DEL
C
C DIVIDE PIVOT ROW BY ITS MAIN DIAGONAL ELEMENT
C
        DIV = A(K,K)
        DO 38 J=1,N
          A(K,J) = A(K,J)/DIV
          B(K,J) = B(K,J)/DIV
   38   CONTINUE
C
C REPLACE EACH ROW BY LINEAR COMBINATION WITH PIVOT ROW
C
        DO 43 I=1,N
          AMULT = A(I,K)
          IF(I-K) 39,43,39
   39     DO 42 J=1,N
            A(I,J) = A(I,J)-AMULT*A(K,J)
            B(I,J)=B(I,J)-AMULT*B(K,J)
   42    CONTINUE
   43   CONTINUE
C
   45 CONTINUE
      GOTO 99
C
C
   93 WRITE(*,113) K; !pause
  113 FORMAT(25H SINGULAR MATRIX FOR K=, I2)
C
   99 RETURN
      END
C
C
C
      SUBROUTINE REINSCH(X1,Y1,DY1,C,A,ND,S)
C--------------------------------------------------------------------
C
C Smoothing of data points by Reinsch-type cubic spline
C
C       Input
C            X,Y,DY ; arrays with x(i),y(i),dy(i), as abscissa,
C                     ordinate and relative weight of i-th data
C                     point, respectively.  The components of the
C                     X must be strictly increasing.
C
C            S      ; a non-negative parameter which controls the
C                     extent of smoothing; the spline function f
C                     is determined such that
C
C                     SUM((f(i)-y(i))/dy(i))**2<=S,
C
C                     where equality holds unless f describes a
C                     straight line.
C
C       *************************************************************
C
C       Output
C            C,A    ; arrays, collecting the coefficients of the
C                     cubic spline f such that with
C                     h=xx-x(i) we have f(xx)=((c(i,3)*h+c(i,2))*h
C                     +c(i,1))*h+a(i)
C                     if x(i)<=xx<x(i+1), i=n1,......n2-1 where n1
C                     n2 are number of first and last data point,
C                     n2>n1.  Furthermore, a(n2)=f(x(n2)) and 
C                     c(n2,2)=0 while c(n2,1) and c(n2,3) are left
C                     undefined.
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(KD=30000) 
      DIMENSION X1(KD),Y1(KD),DY1(KD),X(KD),Y(KD),DY(KD),C(KD,3),
     &          A(KD),R1(KD),R2(KD),T(KD),T1(KD),U(KD),V(KD),R(KD)
C 
C This loop is used to shift data no by +1 to avoid zero or
C       negative subscription in arrays.
C
      DO 100 I=ND+1,2,-1
        X(I)  = X1(I-1)
        Y(I)  = Y1(I-1)
        DY(I) = DY1(I-1)
  100 CONTINUE
C 
C The meaning of N1=2 and N2=ND+1 is the same as above.
C
      N1 = 2
      N2 = ND+1
      M1 = N1-1
      M2 = N2+1
      R(M1)  = 0.D0
      R(N1)  = 0.D0
      R1(N2) = 0.D0
      R2(N2) = 0.D0
      R2(M2) = 0.D0
      U(M1)  = 0.D0
      U(N1)  = 0.D0
      U(N2)  = 0.D0
      U(M2)  = 0.D0
      P      = 0.D0
C 
      M1 = N1+1
      M2 = N2-1
      H  = X(M1)-X(N1)
      F  = (Y(M1)-Y(N1))/H
      DO 110 I=M1,M2
        G     = H
        H     = X(I+1)-X(I)
        E     = F
        F     = (Y(I+1)-Y(I))/H
        A(I)  =  F-E
        T(I)  =  2.D0*(G+H)/3.D0
        T1(I) =  H/3.D0
        R2(I) =  DY(I-1)/G
        R(I)  =  DY(I+1)/H
        R1(I) = -DY(I)/G-DY(I)/H
  110 CONTINUE
C 
      DO 120 I=M1,M2
        C(I,1) = R(I)*R(I)+R1(I)*R1(I)+R2(I)*R2(I)
        C(I,2) = R(I)*R1(I+1)+R1(I)*R2(I+1)
        C(I,3) = R(I)*R2(I+2)
  120 CONTINUE
      F2 = -S
C 
  200 DO 130 I=M1,M2
        R1(I-1) = F*R(I-1)
        R2(I-2) = G*R(I-2)
        R(I)    = 1.D0/(P*C(I,1)+T(I)-F*R1(I-1)-G*R2(I-2))
        U(I)    = A(I)-R1(I-1)*U(I-1)-R2(I-2)*U(I-2)
        F       = P*C(I,2)+T1(I)-H*R1(I-1)
        G       = H
        H       = C(I,3)*P
  130 CONTINUE
 
      DO 140 I=M2,M1,-1
        U(I) = R(I)*U(I)-R1(I)*U(I+1)-R2(I)*U(I+2)
  140 CONTINUE
C
      E = 0.D0
      H = 0.D0
      DO 150 I=N1,M2
        G    = H
        H    = (U(I+1)-U(I))/(X(I+1)-X(I))
        V(I) = (H-G)*DY(I)*DY(I)
        E    = E+V(I)*(H-G)
  150 CONTINUE
C
      V(N2)  = -H*DY(N2)*DY(N2)
      G      = V(N2)
      E      = E-G*H
      G      = F2
      F2     = E*P*P
      IF(F2.GE.S.OR.F2.LE.G) GO TO 210
C
      F      = 0.D0
      H      = (V(M1)-V(N1))/(X(M1)-X(N1))
      DO 160 I=M1,M2
        G    = H
        H    = (V(I+1)-V(I))/(X(I+1)-X(I))
        G    = H-G-R1(I-1)*R(I-1)-R2(I-2)*R(I-2)
        F    = F+G*R(I)*G
        R(I) =G
  160 CONTINUE
      H = E-P*F
      IF(H.LE.0.D0) GO TO 210
      P = P+(S-F2)/((DSQRT(S/E)+P)*H)
      GO TO 200
C 
C use negative branch of square root, if the sequence of
C  abscissae x(i) is strictly decreasing.
C 
  210 DO 170 I=N1,N2
        A(I)   = Y(I)-P*V(I)
        C(I,2) = U(I)
  170 CONTINUE
      DO 180 I=N1,M2
        H      = X(I+1)-X(I)
        C(I,3) = (C(I+1,2)-C(I,2))/(3.D0*H)
        C(I,1) = (A(I+1)-A(I))/H-(H*C(I,3)+C(I,2))*H
  180 CONTINUE
C 
C This loop is used to shift coefficients
C
      DO 190 I=N1,N2
        A(I-1)   = A(I)
        C(I-1,1) = C(I,1)
        C(I-1,2) = C(I,2)
        C(I-1,3) = C(I,3)
  190 CONTINUE
C
      RETURN
      END
C
C
C
      FUNCTION SPLINT(N,X,YS,C,IND,LI,T)
C--------------------------------------------------------------------
C
C Calculation of Y values at arbitrary X values using the output of
C Reinsch-type spline.                                           
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(N),YS(IND),C(IND,3)
C 
      IF(LI.LT.1.OR.LI.GT.N) LI=1
      IF(T.GE.X(LI)) GO TO 30
      DO 10 I=2,LI
        IF(T.LT.X(I)) GO TO 50
   10 CONTINUE
   30 DO 40 I=LI+1,N-1
        IF(T.LT.X(I)) GO TO 50
   40 CONTINUE
      I  = N
   50 LI = I-1
      D  = T-X(LI)
      SPLINT=((C(LI,3)*D+C(LI,2))*D+C(LI,1))*D+YS(LI)
      RETURN
      END
C
C
C
      FUNCTION SMOOTH(X,Y,DY,NPOINT,RM,F)
C--------------------------------------------------------------------
C
C  Constructs the cubic smoothing spline F to given data(X(I),Y(I),
C  I=1,NPOINT) and RM on the basis of COX's idea.
C
C   -------------------------- INPUT --------------------------------
C    X(I)   : data abscissae, assumed to be strictly increasing.
C    Y(I)   : corresponding data ordinates.
C    DY(I)  : estimate of uncertainty in data, assumed to be positive
C    NPOINT : number of data points, assumed greater than 1.
C    RM     : P=1.0/(1.0+RM)
C
C   --------------------------- OUTPUT ------------------------------
C    F(I)    : the sequence of smoothed ordinates.
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(KD=30000)
      DIMENSION A(KD,4),DY(KD),V(KD,7),X(KD),Y(KD),F(KD)
C 
      P   = 1.D0/(1.D0+RM)
      CALL SETUP(X,DY,Y,NPOINT,V,A(1,4))
      CALL CHOL1D(P,V,A(1,4),NPOINT,A(1,3),A(1,1))
      SFP = 0.D0
      DO 10 I=1,NPOINT
        SFP = SFP+(A(I,1)*DY(I))**2
   10 CONTINUE
        SFP = SFP*36.D0*(1.D0-P)**2
C 
C Compute pol. coefficients from Q*U(in A(.,1)).
C
      SMOOTH = SFP
      SIX1MP = 6.D0*(1.D0-P)
      DO 20 I=1,NPOINT
        A(I,1) = Y(I)-SIX1MP*DY(I)**2*A(I,1)
        F(I)=A(I,1)
   20 CONTINUE
      SIXP = 6.D0*P
      DO 30 I=1,NPOINT
        A(I,3) = A(I,3)*SIXP
   30 CONTINUE
      NPM1 = NPOINT-1
      DO 40 I=1,NPM1
        A(I,4)=(A(I+1,3)-A(I,3))/V(I,4)
        A(I,2)=(A(I+1,1)-A(I,1))/V(I,4)
     &        -(A(I,3)+A(I,4)/3.D0*V(I,4))/2.D0*V(I,4)
   40 CONTINUE
C 
      RETURN
      END
C
C
C
      SUBROUTINE SETUP(X,DX,Y,NPOINT,V,QTY)
C--------------------------------------------------------------------
C
C  To be called in SMOOTH
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(KD=30000)
      DIMENSION DX(KD),QTY(KD),V(KD,7),X(KD),Y(KD)
C 
      NPM1   = NPOINT-1
      V(1,4) = X(2)-X(1)
      DO 11 I=2,NPM1
        V(I,4) = X(I+1)-X(I)
        V(I,1) = DX(I-1)/V(I-1,4)
        V(I,2) = -DX(I)/V(I,4)-DX(I)/V(I-1,4)
        V(I,3)=DX(I+1)/V(I,4)
   11 CONTINUE
C
      V(NPOINT,1) = 0.D0
      DO 12 I=2,NPM1
        V(I,5)=V(I,1)**2+V(I,2)**2+V(I,3)**2
   12 CONTINUE
      IF(NPM1.LT.3) GO TO 14
      DO 13 I=3,NPM1
        V(I-1,6)=V(I-1,2)*V(I,1)+V(I-1,3)*V(I,2)
   13 CONTINUE
   14 V(NPM1,6)=0.D0
      IF(NPM1.LT.4) GO TO 16
      DO 15 I=4,NPM1
        V(I-2,7)=V(I-2,3)*V(I,1)
   15 CONTINUE
   16 V(NPM1-1,7)=0.D0
      V(NPM1,7)=0.D0
C 
C Construct Q-transp.*Y in QTY.
C
      PREV = (Y(2)-Y(1))/V(1,4)
      DO 21 I=2,NPM1
        DIFF   = (Y(I+1)-Y(I))/V(I,4)
        QTY(I) = DIFF-PREV
   21 PREV=DIFF
C 
      RETURN
      END
C
C
C
      SUBROUTINE CHOL1D(P,V,QTY,NPOINT,U,QU)
C--------------------------------------------------------------------
C
C To be called in SMOOTH
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KD=30000)
      DIMENSION QTY(KD),QU(KD),U(KD),V(KD,7)
C 
      NPM1=NPOINT-1
C 
C Construct 6*(1-P)*Q-transp.*(D**2)*Q+P*R
C
      SIX1MP = 6.D0*(1.D0-P)
      TWOP   = 2.D0*P
      DO 2 I=2,NPM1
        V(I,1) = SIX1MP*V(I,5)+TWOP*(V(I-1,4)+V(I,4))
        V(I,2) = SIX1MP*V(I,6)+P*V(I,4)
        V(I,3)=SIX1MP*V(I,7)
    2 CONTINUE
      NPM2 = NPOINT-2
      IF(NPM2.GE.2) GO TO 10
      U(1) = 0.D0
      U(2) = QTY(2)/V(2,1)
      U(3) = 0.D0
      GO TO 41
C 
C Factorization
C
   10 DO 20 I=2,NPM2
        RATIO    = V(I,2)/V(I,1)
        V(I+1,1) = V(I+1,1)-RATIO*V(I,2)
        V(I+1,2) = V(I+1,2)-RATIO*V(I,3)
        V(I,2)   = RATIO
        RATIO    = V(I,3)/V(I,1)
        V(I+2,1) = V(I+2,1)-RATIO*V(I,3)
        V(I,3)=RATIO
   20 CONTINUE
C
C Forward substitution
C
      U(1)   = 0.D0
      V(1,3) = 0.D0
      U(2)   = QTY(2)
      DO 30 I=2,NPM2
        U(I+1)=QTY(I+1)-V(I,2)*U(I)-V(I-1,3)*U(I-1)
   30 CONTINUE
C
C Back substitution
C
      U(NPOINT) = 0.D0
      U(NPM1)   = U(NPM1)/V(NPM1,1)
      DO 40 I=NPM2,2,-1
        U(I)=U(I)/V(I,1)-U(I+1)*V(I,2)-U(I+2)*V(I,3)
   40 CONTINUE
C
C Construct Q*U
C
   41 PREV=0.D0
      DO 50 I=2,NPOINT
        QU(I)=(U(I)-U(I-1))/V(I-1,4)
        QU(I-1)=QU(I)-PREV
        PREV=QU(I)
   50 CONTINUE
      QU(NPOINT) = -QU(NPOINT)
C 
      RETURN
      END
C
C
C
      SUBROUTINE BWF(X,XX,NN,A,B1,B2,M)
C--------------------------------------------------------------------
C
C Butterworth filter
C X ; input data        XX ; output data     NN ; data no
C A, B1, B2 ; filter coefficients    M ; coefficients no
C       Y ; working space
C
C output data are passed through the filter in reverse order
C       to compensate phase delay.
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KO=30000)
C 
      DIMENSION X(KO),Y(KO,20),XX(KO),A(20),B1(20),B2(20)
C 
      DO 100 J=1,2
C
        DO 110 K=1,M
          IF(J.EQ.2) GO TO 200
          Y(1,K) = X(1)
          Y(2,K) = X(1)
          GO TO 110
  200     Y(1,K) = XX(1)
          Y(2,K) = XX(1)
  110   CONTINUE
C
        DO 120 I=3,NN
          DO 125 K=1,M
            IF(K.GE.2) GO TO 130
            IF(J.EQ.1) Y(I,K)=A(K)*(X(I)+2.D0*X(I-1)+X(I-2))
     &                       +B1(K)*Y(I-1,K)
     &                       +B2(K)*Y(I-2,K)
            IF(J.EQ.2) Y(I,K)=A(K)*(XX(I)+2.D0*XX(I-1)+XX(I-2))
     &                       +B1(K)*Y(I-1,K)
     &                       +B2(K)*Y(I-2,K)
            GO TO 125
  130       Y(I,K)=A(K)*(Y(I,K-1)+2.D0*Y(I-1,K-1)+Y(I-2,K-1))
     &            +B1(K)*Y(I-1,K) 
     &            +B2(K)*Y(I-2,K)
  125     CONTINUE
  120   CONTINUE
C 
        N2 = NN/2
        DO 140 I=1,N2
          Y1          = Y(I,M)
          Y2          = Y(NN-I+1,M)
          Y(I,M)      = Y2
          Y(NN-I+1,M) =Y1
  140   CONTINUE
        DO 150 I=1,NN
          XX(I) = Y(I,M)
  150   CONTINUE
C
  100 CONTINUE
C 
      RETURN
      END
C
C
C
      SUBROUTINE YAC(X1,YI,NN,MI,MF,YACLABEL)
C--------------------------------------------------------------------
C
C Calculation of yearly averaged values from time series YI(X1).
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KO=30000) 
      DIMENSION X1(KO),YI(KO),AMC(50),ICD(50)
      CHARACTER YACLABEL*80
      COMMON /SYEAR/IYS
C
C     ITYD1 = INT(TYD1)
C 
      DO 100 I=1,50
        AMC(I) = 0.D0
        ICD(I) = 0
  100 CONTINUE
 
      I0 = 0
      KI = 0
      KF = 0
      IX1= INT(X1(1))
      DO 170 I=1,NN
        IXI = INT(X1(I))
        IF(IXI.LT.0) I1 = 1+IABS((IX1+1)/365)+(IXI+1)/365
        IF(IXI.LT.0) I0 = I1
        IF(IXI.GE.0) I1 = I0+1+IXI/365
        IF(KI.EQ.0.AND.I.GE.MI) KI = I1
        IF(KF.EQ.0.AND.I.GE.MF) KF = I1
        AMC(I1) = AMC(I1)+YI(I)
        ICD(I1) = ICD(I1)+1
  170 CONTINUE
      DO 240 MC=KI,KF
        AMC(MC) = AMC(MC)/DFLOAT(ICD(MC))
  240 CONTINUE
      WRITE(2,670)
      WRITE(2,675) YACLABEL
  670 FORMAT(/' * ANNUAL MEAN VALUES')
  675 FORMAT(1H ,A40/)
      DO 180 MC=KI,KF
        WRITE(2,680) IYS+MC-KI,AMC(MC),ICD(MC)
  180 CONTINUE  
  680 FORMAT(1H ,I6,2X,F8.3,1X,1H(,I3,1H))
      WRITE(2,660)
  660 FORMAT(1H 72(1H-))
      RETURN
      END
C
C
C
      FUNCTION SEASON1(A,PI,T,NH)
C--------------------------------------------------------------------
C
C    Calculation of the seasonal cycle
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(20)
C
      N  = NH/2
      YX = 0.D0
      IF(NH.EQ.0) GO TO 110
      DO 100 I=1,N
        YX = YX+A(2*I-1)*SIN(DFLOAT(2*I)*PI*T)
     &         +A(2*I)  *COS(DFLOAT(2*I)*PI*T)
  100 CONTINUE
C
  110 SEASON1=YX
C 
      RETURN
      END
C
C
C
      SUBROUTINE SCYCLE(S,X,NN,MI,MF,KI,KF,SMAX,SMIN,MAXD,MIND)
C--------------------------------------------------------------------
C
C    Calculation of the amplitude, maximum and minimum concentrations
C    of the seasonal cycle in each year.
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KO=30000)
      DIMENSION S(KO),X(KO),SMIN(50),SMAX(50),MIND(50),MAXD(50)
C 
      DO 1 I=1,50
        SMIN(I) = 0.D0
        SMAX(I) = 0.D0
        MIND(I) = 0
        MAXD(I) = 0
    1 CONTINUE
C
      I0 = 0
      KI = 0
      KF = 0
      IX = INT(X(1))
      DO 100 I=1,NN
        IXI = INT(X(I))
        IF(IXI.LT.0) I1 = 1+IABS((IX+1)/365)+(IXI+1)/365
        IF(IXI.LT.0) I0 = I1
        IF(IXI.GE.0) I1 = I0+1+IXI/365
        IF(KI.EQ.0.AND.I.GE.MI) KI = I1
        IF(KF.EQ.0.AND.I.GE.MF) KF = I1
        IF(S(I).GT.SMIN(I1)) GO TO 110
        SMIN(I1) = S(I)
        MIND(I1) = INT(X(I)+0.5D0)
        GO TO 100
  110   IF(S(I).LT.SMAX(I1)) GO TO 100
        SMAX(I1) = S(I)
        MAXD(I1) = INT(X(I)+0.5D0)
  100 CONTINUE
 
      DO 120 I=KI,KF
       IF(MAXD(I).LT.0) MAXD(I)=(IABS(MAXD(I)+1)/365+1)*365+MAXD(I)+1
       IF(MAXD(I).GE.0) MAXD(I)= MAXD(I)-FLOAT(MAXD(I)/365)*365+1
       IF(MIND(I).LT.0) MIND(I)=(IABS(MIND(I)+1)/365+1)*365+MIND(I)+1
       IF(MIND(I).GE.0) MIND(I)= MIND(I)-FLOAT(MIND(I)/365)*365+1
  120 CONTINUE
C 
      RETURN
      END
C
C
C
      SUBROUTINE SCYCLE2(S,X,NN,MI,MF,KI,KF,SMAX,SMIN,MAXD,MIND)
C--------------------------------------------------------------------
C
C    Calculation of the maximum and minimum concentrations
C    and their timing of the seasonal cycle in each year.
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KO=30000)
      DIMENSION S(KO),X(KO),SMIN(50),SMAX(50),MIND(50),MAXD(50),
     &          ITYEAR(0:50)
C
      COMMON /SYEAR/IYS
      write(*,*) 'start seasaonl'
C 
      DO 1 I=1,50
        SMIN(I) = 0.D0
        SMAX(I) = 0.D0
        MIND(I) = 0
        MAXD(I) = 0
    1 CONTINUE
C
      I0  = 0
      KI  = 0
      KF  = 0
      IAY = 0
      I1  = 0
      ITYEAR(I1) = 0
      I  = 1
C
   50 IF(X(I).GE.0) THEN
  100   I1 = I1 + 1
        IDYEAR = 365
        IF(MOD(IYS+IAY,4).EQ.0) IDYEAR = 366
        DO 10 J = I,I+IDYEAR-1
          IF((KI.EQ.0).AND.(J.GE.MI)) KI = I1
          IF((KF.EQ.0).AND.(J.GE.MF)) KF = I1
          IF(S(J).GT.SMIN(I1)) GO TO 20
          SMIN(I1) = S(J)
          MIND(I1) = INT(X(J))
   20     IF(S(J).LT.SMAX(I1)) GO TO 10
          SMAX(I1) = S(J)
          MAXD(I1) = INT(X(J))
   10   CONTINUE
        IAY = IAY + 1
        ITYEAR(I1) =  ITYEAR(I1-1) + IDYEAR
        I = (I+IDYEAR-1)+1
        IF(I.GT.MF+366) GO TO 200
        GO TO 100
      ELSE
        I = I+1
        GO TO 50
      ENDIF
C
  200 DO 120 I=KI,KF
        MAXD(I) = MAXD(I) + 1 - ITYEAR(I-1)
        MIND(I) = MIND(I) + 1 - ITYEAR(I-1)
  120 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE EXTRAD1(W,XMI,XMF,DD,IP)
C--------------------------------------------------------------------
C
C    Time of the data that are added to both ends of the record or 
C  That are interpolated within the record.
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(KD=30000) 
      DIMENSION W(KD)
C 
      IB = INT(DABS(XMF-XMI)/DD)
      DS = DABS(XMF-XMI)/DFLOAT(IB)
      IP = IB-1
      IF(XMF.GT.XMI) GO TO 100
C 
C ---- left side of the record  -------------------------
C
      XMF = XMI-DS*DFLOAT(IB)
      DO 110 I=1,IP
        W(I) = XMF+DS*DFLOAT(I)
  110 CONTINUE
      GO TO 120
C 
C ---- right side of the record or interpolation  -------
C
  100 DO 130 I=1,IP
        W(I) = XMI+DS*DFLOAT(I)
  130 CONTINUE
C 
  120 RETURN
      END
C
C
C
      SUBROUTINE JUDGE(ISTP,SIGM1,SIGM,Y,MDI,MDF,MD,X2,TF1,TR1,NN)
C--------------------------------------------------------------------
C
C This subroutine is used in the program "LSFCO2F" to judge the
C convergence of an iterative procedure.
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KD=30000,KO=60000)
      DIMENSION Y(KD),X2(KO),TF1(KO),TR1(KO),MDI(KD),MDF(KD)
 
      IF(1.D0-SIGM/SIGM1.GT.0.0001D0) THEN
        SIGM1 = SIGM
      ELSE
        ISTP = 1
        GO TO 200
      END IF
C 
      LI = 1
      DO 100 I=1,MD
        DO 110 J=MDI(I),MDF(I)
          DO 120 K=LI,NN
            IF(DABS(X2(K)-TF1(J)).GT.0.002D0) GO TO 120
            TR1(J) = Y(K)
            LI     = K
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
C 
  200 RETURN
      END
C
C
C
      SUBROUTINE AIC(T,SD,A,PI,NH,NDATA,AI)
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KD=30000) 
      DIMENSION T(KD),SD(KD),A(20),R(KD)
C 
      S = 0.D0
      DO 36 I=1,NDATA
        R(I) = SD(I)-SEASON1(A,PI,T(I),NH)
        S    = S+R(I)*R(I)
   36 CONTINUE
C 
      S  = S/0.0009D0
      AI = NDATA*DLOG(S)+2*NH
C 
      WRITE(*,1000) AI
 1000 FORMAT(/'   AIC = ',F15.5)
C 
      RETURN
      END
C
C
C
      SUBROUTINE SORT1(T,Y,NDATA,OUTFILE)
C--------------------------------------------------------------------
C
C    SORT OF OUTLIERS
C
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (KL=100)
      DIMENSION T(KL),Y(KL)
      CHARACTER*80 OUTFILE
C
      DO 10 K=1,NDATA-1
        M=K
        DO 20 I=K+1,NDATA
          IF(T(I).LT.T(M)) M=I
   20   CONTINUE
        TEMPT = T(K)
        T(K)  = T(M)
        T(M)  = TEMPT
        TEMPY = Y(K)
        Y(K)  = Y(M)
        Y(M)  = TEMPY
   10 CONTINUE
C
      OPEN(UNIT=21,FILE=OUTFILE)
      DO 30 I=1,NDATA
        WRITE(*,100) T(I),Y(I)
        WRITE(21,100) T(I),Y(I)
   30 CONTINUE
      CLOSE(21)
  100 FORMAT(1H ,2(F12.5))
      RETURN
      END
C
C
C
      SUBROUTINE BWFCOEF(SA,SB,SC,NUM,CUT,IBWF)
C--------------------------------------------------------------------
C
C   Calculation of coefficients of butterworth filter
C
C   FS : Sampling periods (1 day)
C   FC : Cut off periods (day)
C   N  : Order of Butterworth filter
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX Z1,Z2,F1,F2,TF
      PARAMETER (KI=20,KF=2000)
      DIMENSION SI(0:KI),SA(KI),SB(KI),SC(KI)
      DIMENSION FQ(KF),FQ2(KF),TF(KF),TFT(KF)
      CHARACTER IBWF*1
C
      COMMON /DAYYEAR/DYEAR,DYEAR0,DYEAR1,DYEAR2
      PI=3.141592653589793D0
C
      FS = 1.0D0
      FC = CUT/12.0D0*DYEAR0
      N  = NUM*2
C
      FS=1.D0/FS
      FC=1.D0/FC
      FC=FC/FS
      WC=DTAN(PI*FC)/(1.D0/0.5D0-1.D0)**(1.D0/DBLE(2*N))
C
      DO 10 I=0,N/2-1
        SI(I)=DSIN(PI*(2*DBLE(I)+1.D0)/(2*DBLE(N)))
   10 CONTINUE
C
      DO 20 I=1,N/2
        SD=WC*WC+2*WC*SI(I-1)+1.D0
        SA(I)=WC*WC/SD
        SB(I)=-2*(WC*WC-1.D0)/SD
        SC(I)=-(WC*WC-2*WC*SI(I-1)+1.D0)/SD
   20 CONTINUE
C
C------------------------------------------------------
C
      IF(IBWF.NE.'y') GOTO 80
C
      WRITE(10,*) N/2
      DO 30 I=1,N/2
        WRITE(10,*) SA(I),SB(I),SC(I)
   30 CONTINUE
      WRITE(10,*) 
C
C--- transfer function  -------------------------------
C
      TL =   1.0D0
      TY =   DYEAR0
      DO 40 I=2,1000
        T = DBLE(I)
        FQ(I)  = TL/T
        FQ2(I) = TY/T
        Z1    = CMPLX(DCOS(2.0D0*PI*FQ(I)),DSIN(2.0D0*PI*FQ(I)))
        Z2    = CMPLX(DCOS(2.0D0*PI*FQ(I)),-DSIN(2.0D0*PI*FQ(I)))
        TF(I) = CMPLX(1.0D0,1.0D0)
        DO 50 J=1,N/2
          F1    = SA(J)*(Z1**2+2.0D0*Z1+1.0D0)/(Z1**2-SB(J)*Z1-SC(J))
          F2    = SA(J)*(Z2**2+2.0D0*Z2+1.0D0)/(Z2**2-SB(J)*Z2-SC(J))
          TF(I) = TF(I)*F1*F2
   50   CONTINUE
        TFT(I) = DBLE(TF(I))
   40 CONTINUE
C
C
      DO 60 I = 2,1000
        T = DBLE(I)
        WRITE(10,110) FQ2(I),T,TFT(I)
   60 CONTINUE
      WRITE(10,*)
  110 FORMAT(1H ,3(F15.6))
   80 RETURN
      END
C
C
C
      SUBROUTINE MDINPUT1(YEAR,MONTH,DAY,TYDIN,JULIN,NDATA)
C---------------------------------------------------------------
C
C     This subroutine calculates the Time and Julian day
C     corresponding the set of Year, Month and Day.
C
C---------------------------------------------------------------
C
      PARAMETER (KD=30000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER MON1(12),MON2(12),TMON1(20),TMON2(20)
      INTEGER YEAR(KD),MONTH(KD),DAY(KD)
      DIMENSION TYDIN(KD),JULIN(KD)
C
      data MON1 /31,28,31,30,31,30,31,31,30,31,30,31/
      data MON2 /31,29,31,30,31,30,31,31,30,31,30,31/
C
      TMON1(1) = 0
      TMON2(1) = 0
      DO 10 I=2,12
        TMON1(I) = TMON1(I-1) + MON1(I-1)
        TMON2(I) = TMON2(I-1) + MON2(I-1)
   10 CONTINUE
C
      DO 30 I=1,NDATA
        IM     = INT(MONTH(I))
        IF(MOD(YEAR(I),4).NE.0) THEN
          TMONTH = TMON1(IM)
          DYEAR  = 365.0D0
        ELSE
          TMONTH = TMON2(IM)
          DYEAR  = 366.0D0
        ENDIF
        JULIN(I) = DBLE(TMONTH) + DAY(I)           ! PKP
        TYDIN(I) = YEAR(I) + (JULIN(I)-1+0.5)/DYEAR !365.0
   30 CONTINUE   
C
      RETURN
      END
C
C
C
      SUBROUTINE MDINPUT2(YEARIN1,MONTHIN1,DAYIN1,HOUR1,YIN1,
     &                    NDATA1,INUM,OUTVALUE,
     &                    TYDIN,JULIN,YIN,NDATA0)
      PARAMETER (KL=60000,KD=30000,KM=30000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER YEARIN1(KL),MONTHIN1(KL),DAYIN1(KL),HOUR1(KL),
     &          OUTYEAR(KM),OUTMONTH(KM),OUTDAY(KM),OUTHOUR(KM),
     &          YEARIN2(KL),MONTHIN2(KL),DAYIN2(KL),
     &         YEARIN(KD),MONTHIN(KD),DAYIN(KD)
      DIMENSION YIN1(KL),OUTY(KM),YIN2(KL),YIN(KD)
      DIMENSION TYDIN(KD),JULIN(KD),TYIN(30)
C
      CHARACTER INUM*1
C
      I2 = 0
      J  = 0
      IF(INUM.EQ.'y') THEN
        DO 10 I=1,NDATA1
          IF(YIN1(I).EQ.OUTVALUE) THEN
            J = J+1
            OUTYEAR(J) = YEARIN1(I)
            OUTMONTH(J)= MONTHIN1(I)
            OUTDAY(J)  = DAYIN1(I)
            OUTHOUR(J) = HOUR1(I)
            OUTY(J)    = YIN1(I)
          ELSE
            I2 = I2 + 1
            YEARIN2(I2) = YEARIN1(I)
            MONTHIN2(I2)= MONTHIN1(I)
            DAYIN2(I2)  = DAYIN1(I)
            YIN2(I2)    = YIN1(I)
          ENDIF
   10   CONTINUE
      ELSE
        DO 20 I=1,NDATA1
          I2 = I2+1
          YEARIN2(I2)  = YEARIN1(I)
          MONTHIN2(I2) = MONTHIN1(I)
          DAYIN2(I2)   = DAYIN1(I)
          YIN2(I2)     = YIN1(I)
   20   CONTINUE
      ENDIF
C
      NDATA2   = I2
      NDATAOUT = J
C
      I = 1
      K = 1
      TYIN(K) = YIN2(I)
      J = 0
      DO 30 I=2,NDATA2+1
        K = K+1
        TYIN(K) = TYIN(K-1) + YIN2(I)
        IF((DAYIN2(I).NE.DAYIN2(I-1))
     &      .OR.(MONTHIN2(I).NE.MONTHIN2(I-1))) THEN
          J = J+1
          YIN(J)     = TYIN(K-1)/(K-1)
          YEARIN(J)  = YEARIN2(I-1)
          MONTHIN(J) = MONTHIN2(I-1)
          DAYIN(J)   = DAYIN2(I-1)
          K = 1
          TYIN(K) = YIN2(I)
        ENDIF
   30 CONTINUE
      NDATA0 = J
C
      CALL MDINPUT1(YEARIN,MONTHIN,DAYIN,TYDIN,JULIN,NDATA0)
C
      WRITE(2,*) 'Number of data ',NDATA2+NDATAOUT
      WRITE(*,*) 'Number of data ',NDATA2+NDATAOUT
      WRITE(2,*) 'Missing data   ',NDATAOUT
      WRITE(*,*) 'Missing data   ',NDATAOUT
      IF(NDATAOUT.GT.0) THEN
        DO 40 I=1,NDATAOUT
          WRITE(2,*) OUTYEAR(I),OUTMONTH(I),OUTDAY(I),
     &               OUTHOUR(I),OUTY(I)
          WRITE(*,*) OUTYEAR(I),OUTMONTH(I),OUTDAY(I),
     &               OUTHOUR(I),OUTY(I)
   40   CONTINUE
      ENDIF
C
      RETURN
      END       
C
C
C
      SUBROUTINE MDINPUT3(TYDIN1,JULIN1,HOUR1,YIN1,
     &                    NDATA1,INUM,OUTVALUE,
     &                    TYDIN,JULIN,YIN,NDATA0)
      PARAMETER (KL=60000,KD=30000,KM=30000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER HOUR1(KL),OUTHOUR(KM)
      DIMENSION TYDIN1(KL),JULIN1(KL),YIN1(KL),
     &          OUTTYD(KM),OUTJUL(KM),OUTY(KM)
      DIMENSION TYDIN2(KL),JULIN2(KL),YIN2(KL),
     &          TYDIN(KD),JULIN(KD),YIN(KD)
      DIMENSION TYIN(30)
C
      CHARACTER INUM*1
C
      I2 = 0
      J  = 0
      IF(INUM.EQ.'y') THEN
        DO 10 I=1,NDATA1
          IF(YIN1(I).EQ.OUTVALUE) THEN
            J = J+1
            OUTTYD(J)  = TYDIN1(I)
            OUTJUL(J)  = JULIN1(I)
            OUTHOUR(J) = HOUR1(I)
            OUTY(J)    = YIN1(I)
          ELSE
            I2 = I2 + 1
            TYDIN2(I2)   = TYDIN1(I)
            JULIN2(I2) = JULIN1(I)
            YIN2(I2)   = YIN1(I)
          ENDIF
   10   CONTINUE
      ELSE
        DO 20 I=1,NDATA1
          I2 = I2+1
          TYDIN2(I2) = TYDIN1(I)
          JULIN2(I2) = JULIN1(I)
          YIN2(I2)   = YIN1(I)
   20   CONTINUE
      ENDIF
C
      NDATA2   = I2
      NDATAOUT = J
C
      I = 1
      K = 1
      TYIN(K) = YIN2(I)
      J = 0
      DO 30 I=2,NDATA2+1
        K = K+1
        TYIN(K) = TYIN(K-1) + YIN2(I)
        IF(JULIN2(I).NE.JULIN2(I-1)) THEN
          J = J+1
          YIN(J)   = TYIN(K-1)/(K-1)
          TYDIN(J) = TYDIN2(I-1)
          JULIN(J) = JULIN2(I-1) 
          K = 1
          TYIN(K) = YIN2(I)
        ENDIF
   30 CONTINUE
      NDATA0 = J
C
C
      WRITE(2,*) 'Number of data ',NDATA2+NDATAOUT
      WRITE(*,*) 'Number of data ',NDATA2+NDATAOUT
      WRITE(2,*) 'Missing data   ',NDATAOUT
      WRITE(*,*) 'Missing data   ',NDATAOUT
      IF(NDATAOUT.GT.0) THEN
        DO 40 I=1,NDATAOUT
          WRITE(2,*) OUTTYD(I),OUTJUL(I),OUTHOUR(I),OUTY(I)
          WRITE(*,*) OUTTYD(I),OUTJUL(I),OUTHOUR(I),OUTY(I)
   40   CONTINUE
      ENDIF
C
      RETURN
      END       
      

