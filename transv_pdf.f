C
*********************************************************************
*                                                                   *
*                    TRANSVERSITY DENSITIES                         *
*                                                                   *
*   INPUT:   ISET = number of the parton set :                      *
*              ISET = 1  NEXT-TO-LEADING ORDER  (MS-bar)            * 
*                        (DATA FILE 'transmaxnlo.grid' UNIT=11, TO  *
*                         BE DEFINED BY THE USER )                  * 
*              ISET = 2  LEADING ORDER                              * 
*                        (DATA FILE 'transmaxlo_new.grid' UNIT=22   *
*                                                                   *
*            X  = Bjorken-x       (between  1.E-4  and  1)          *
*            Q2 = scale in GeV**2 (between  0.8  and   1.E6)        *
*             (for values outside the allowed range the program     *
*              writes a warning and extrapolates to the x and       *
*              Q2 values requested)                                 *
*                                                                   *
*   OUTPUT:  U = x * DELTA u                                        * 
*            D = x * DELTA d                                        *
*            UB = x * UBAR                                          *
*            DB = x * DBAR                                          *   
*            ST = x * DELTA STRANGE = x * DELTA STRANGE(BAR)        *     
*                                                                   *
*          (  For the parton distributions always x times           *
*                   the distribution is returned   )                *
*                                                                   *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINI / IINI , and  IINI     *
*            has always to be zero when PARPOL is called for the    *
*            first time or when 'ISET' has been changed.            *
*                                                                   *
* ALEXEI PROKUDIN 13/12/2006 IF YOU WANT TO USE THE PROGRAMM IN C++ *
* PLEASE CHANGE INTINI -> INTINI2, IINI -> IINI2, DFINT -> DFINT2   *
* OVERWISE INITILISATION OF PARPOLT WILL COINCISE WITH INITIALISATION *
* OF PARPOL...                                                      *
*********************************************************************
*
      SUBROUTINE PARPOLT(ISET, X, Q2, U, D, UB, DB, ST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPART=6, NX=42, NQ=30, NARG=2)
      DIMENSION XUF(NX,NQ), XDF(NX,NQ), XUBF(NX,NQ), 
     1          XDBF(NX,NQ), XSF(NX,NQ), PARTON (NPART,NQ,NX-1), 
     2          QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      COMMON / INTINI2 / IINI2
      SAVE XUF, XDF, XUBF, XDBF, XSF, NA, ARRF
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.8D0, 1.0D0, 1.25d0, 1.5D0, 2.d0, 2.5D0, 
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 
     4           1.0D5, 1.8D5, 3.2D5, 5.8D5, 1.0D6  /
       DATA XB / 
     1           1.D-4, 1.5D-4, 2.2D-4, 3.2D-4, 4.8D-4, 7.D-4,
     2           1.D-3, 1.5D-3, 2.2D-3, 3.2D-3, 4.8D-3, 7.D-3,
     3           1.D-2, 1.5D-2, 2.2D-2, 3.2D-2, 5.0D-2, 7.5D-2,
     4           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5           0.3, 0.325, 0.35, 0.375, 0.4, 0.45,  0.5, 0.55,
     6           0.6, 0.65,  0.7,  0.75,  0.8, 0.85,  0.9, 1.0 /
*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.1.0D-4) .OR. (X.GT.1.0D0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
          GOTO 60
       ENDIF
       IF ( (Q2.LT.0.8D0) .OR. (Q2.GT.1.D6) ) THEN
           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
          GOTO 60
       ENDIF
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
*    FILE - NO. = 11 FOR NLO  ( FIRST NUMBER IN THE GRID: 1.102E-03 )     
*    FILE - NO. = 22 FOR LO   ( FIRST NUMBER IN THE GRID: 9.5433E-02 )
      IF (IINI2.NE.0) GOTO 16
      IF (ISET.EQ.1) THEN
       IIREAD=11       
       OPEN(UNIT=11,FILE='transmaxnlo_new.grid',STATUS='OLD')
      ELSE IF (ISET.EQ.2) THEN
       IIREAD=22
       OPEN(UNIT=22,FILE='transmaxlo_new.grid',STATUS='OLD')
      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'PARTON INTERPOLATION: ISET OUT OF RANGE')
        GOTO 60
      END IF
C
       DO 15 M = 1, NX-1 
       DO 15 N = 1, NQ
       READ(IIREAD,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1                 PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M)
  90   FORMAT (6(1PE12.4))
  15   CONTINUE
C
      IINI2 = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.D0-XB(IX)
        XUF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0)
        XDF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0)
        XUBF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**8 * XB0**0.5) 
        XDBF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**8 * XB0**0.5) 
        XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**8 * XB0**0.5)
  20  CONTINUE
        XUF(NX,IQ) = 0.D0
        XDF(NX,IQ) = 0.D0
        XUBF(NX,IQ) = 0.D0
        XDBF(NX,IQ) = 0.D0
        XSF(NX,IQ)  = 0.D0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      U = DFINT2(NARG,XT,NA,ARRF,XUF) * (1.D0-X)**3 * X
      D = DFINT2(NARG,XT,NA,ARRF,XDF) * (1.D0-X)**4 * X
      UB = DFINT2(NARG,XT,NA,ARRF,XUBF) * (1.D0-X)**8 * X**0.5
      DB = DFINT2(NARG,XT,NA,ARRF,XDBF)  * (1.D0-X)**8 * X**0.5
      ST = DFINT2(NARG,XT,NA,ARRF,XSF)  * (1.D0-X)**8 * X**0.5
 60   RETURN
      END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION DFINT2(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARG(5),NENT(5),ENT(73),TABLE(1200)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      DFINT2=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      DFINT2=DFINT2+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END

