      SUBROUTINE STOP(N,X,DX,F,FNOR,G,SX,SF,IRETCD,ITER,MAXTKN,         
     +       KMAXDU,TERMCD,EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT)

C   detects reasons to stop iterations
C     * TERMCD=0 - do not stop;

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER KPRSOL
      DOUBLE PRECISION             EPSSOL,EPSDU,EPSMIN,MAXDU
      INTEGER          TERMCD
      LOGICAL          MAXTKN
      DOUBLE PRECISION             X(1),DX(1),F(1),G(1)
      DOUBLE PRECISION             SX(1),SF(1)
      TERMCD=0
	  KPRSOL=2

      IF(IRETCD.NE.1)GOTO 10
C  in LSERCH was not done appropriate step
      TERMCD=3
      IF(KPRSOL.GT.0) WRITE(6, 103) ITER,EPSDU,FNOR
      FLUSH (6)
      RETURN
  10  CONTINUE

C  MAX scaled memeber of residual vector:
      CONV=0.D0
      DO 20 I=1,N
  20  CONV=DMAX1(CONV,SF(I)*DABS(F(I)))
      IF(CONV.GT.EPSSOL) GOTO 30
C  we have approximate solution (if EPSSOL not too large)
      TERMCD=1
      IF(KPRSOL.GT.0)  WRITE(6,101) ITER,CONV,EPSSOL,FNOR
      FLUSH(6)
      RETURN
  30  CONTINUE

C  MAX scaled member of corrections vector
      DUNOR=0.D0
      DO 40 I=1,N
      SINV=1/SX(I)
  40  DUNOR=DMAX1(DUNOR,DABS(DX(I)/DMAX1(DABS(X(I)),SINV)))
      IF(DUNOR.GT.EPSDU) GOTO 50

C  correction is less than allowed
      TERMCD=2
      IF(KPRSOL.GT.0) WRITE(6, 102) ITER,DUNOR,EPSDU,CONV,FNOR
      FLUSH(6)
      RETURN
  50  CONTINUE

      IF(ITER.LT.LIMIT)GOTO 60

C  iterations limit exhausted
      TERMCD=4
      IF(KPRSOL.GT.0) WRITE(6, 104) ITER,CONV,EPSSOL,DUNOR,EPSDU
      FLUSH(6)
      RETURN
  60  CONTINUE

      IF(.NOT.MAXTKN)GOTO 70

C  step was made with length=MAXDU
      KMAXDU=KMAXDU+1
      IF(KMAXDU.LT.5)GOTO 70
C  5 steps was made with length=MAXDU
      TERMCD=5
      IF(KPRSOL.GT.0) WRITE(6, 105) ITER,MAXDU,CONV,EPSSOL
      FLUSH(6)
      RETURN
  70  CONTINUE

      KMAXDU=0
C   define relative "speed" of decreasing of FNOR
      TERM=0.D0


      DO 80 I=1,N
      SINV=1/SX(I)
      TERM=DMAX1(TERM,DABS(G(I))*DMAX1(X(I)      ,SINV)/                
     + DMAX1(FNOR,DFLOAT(N)/2.D0))

C     PRINT 222 removed Koshmanova N.V.
  80  CONTINUE
      IF(TERM.GT.EPSMIN)GOTO 90
C  here wa are at lockal minimum of FNOR
      TERMCD=6
      IF(KPRSOL.GT.0) WRITE (6,106) ITER,TERM,EPSMIN,FNOR,CONV,EPSSOL
      FLUSH(6)
      RETURN
  90  CONTINUE
C   nothing happened
      IF(KPRSOL.GE.2) WRITE(6,107) CONV,DUNOR,TERM
      FLUSH(6)
      RETURN

 101  FORMAT('SUCCESS.'/
     *       '  AT ',I5,'-TH ITERATION CONVERGED TO SOLUTION    '/
     *       '  WITH RESIDUAL <= ',E12.6 ,' ( < ',E12.6,' )'/
     *       '  1/2 OF SQUARED L-2 NORM OF RESIDUAL =',E12.6 )

 102  FORMAT('FAILURE IS POSSIBLE. CHECK RESIDUAL.'/
     *       '  AT ',I5,' ITERATION MAX. STEP =',E12.6,' ( <',E12.6,')'/ 
     *       '  RESIDUAL =',E12.6/
     *       '  1/2 SQUARED NORM OF RESIDUAL =',E12.6)

 103  FORMAT('FAILURE IS POSSIBLE. CHECK RESIDUAL.'/
     *       '  AT ',I5,' ITERATION IMPOSSIBLE OBTAIN'/   
     *       '  ALLOWABLE STEP > ',E12.6 /       
     *       '  1/2 SQUARED NORM OF RESID. =',E12.6)

 104  FORMAT('FAILURE.'/
     *       '  EXHAUSTED LIMIT ',I4,' OF ITERATIONS'/
     *       '  RESIDUAL =',E12.6,'( > ',E12.6,' );'/
     *       '  STEP =',E12.6,'( > ',E12.6,' ).')

 105  FORMAT('FAILURE.'/
     *       '  AT ',I4,'-TH ITERATION'/
     *       '  5 STEPS WITH LENGTH ',E12.6,' WAS MADE.'/
     *       '  RESIDUAL =',E12.6,' ( > ',E12.6,' ).')

 106  FORMAT('FAILURE.'/
     *      '    LOCAL MINIMUM OF 1/2*|F|*|F|  :  '/ 
     *      '     ITERATION                     ',I4/
     *      '     MEASURE OF GRADIENT OF RESID. ',E12.6,
     *      '(<',E12.6 ,')'/
     *      '     1/2 SQUARED L-2 NORM OF RESID.',E12.6/
     *      '     RESIDUAL                      ',E12.6/
     *      '     REQUIRED                      ',E12.6)

 107  FORMAT('  RESIDUAL=',E12.6,', STEP=',E12.6,
     *',  RESIDUAL DECR.  =', E12.6)

      END
