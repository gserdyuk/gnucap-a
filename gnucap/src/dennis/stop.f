      SUBROUTINE STOP(N,X,DX,F,FNOR,G,SX,SF,IRETCD,ITER,MAXTKN,         
     +       KMAXDU,TERMCD,EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT)

C   �/� O�PE�E�EH�� �P���H� OCTAHOBA.
C     * TERMCD=0 - HET OCTAHOBA;
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
C  B LSERCH HAM HE��A�OC� C�E�AT� ��OB�ETBOP�TE��H�� �A�
C  in LSERCH was not done appropriate step
      TERMCD=3
      IF(KPRSOL.GT.0) WRITE(6, 103) ITER,EPSDU,FNOR

*                     PRINT    103, ITER,EPSDU,FNOR
      RETURN
  10  CONTINUE

C  MAX.MAC�TA��POBAHH�� ��EH HEB��K�:
C  MAX scaled memeber of residual vector:
      CONV=0.D0
      DO 20 I=1,N
  20  CONV=DMAX1(CONV,SF(I)*DABS(F(I)))
      IF(CONV.GT.EPSSOL) GOTO 30
C  M� �MEEM �P������TE��HOE PE�EH�E
C   ( EC�� EPSSOL HE O�EH� BE��KO)
C  we have approximate solution (if EPSSOL not too large)
      TERMCD=1
      IF(KPRSOL.GT.0)  PRINT    101, ITER,CONV,EPSSOL,FNOR
      RETURN
  30  CONTINUE

C  MAX. MAC�TA��POBAHH�� ��EH �O�PABK�
C  MAX scaled member of corrections vector
      DUNOR=0.D0
      DO 40 I=1,N
      SINV=1/SX(I)
  40  DUNOR=DMAX1(DUNOR,DABS(DX(I)/DMAX1(DABS(X(I)),SINV)))
      IF(DUNOR.GT.EPSDU) GOTO 50
C  �O�PABKA MEH��E M�H. �O��CT�MO�.
C  correction is less than allowed
      TERMCD=2
      IF(KPRSOL.GT.0) WRITE(6, 102) ITER,DUNOR,EPSDU,CONV,FNOR
*                      PRINT    102, ITER,DUNOR,EPSDU,CONV,FNOR
      RETURN
  50  CONTINUE

      IF(ITER.LT.LIMIT)GOTO 60
C  �C�EP�AH ��M�T �TEPA���
C  iterations limit exhausted
      TERMCD=4
      IF(KPRSOL.GT.0) WRITE(6, 104) ITER,CONV,EPSSOL,DUNOR,EPSDU
*                      PRINT    104, ITER,CONV,EPSSOL,DUNOR,EPSDU
      RETURN
  60  CONTINUE

      IF(.NOT.MAXTKN)GOTO 70
C  ��� C�E�AH �A� ���HO� MAXDU
C  step was made with length=MAXDU
      KMAXDU=KMAXDU+1
      IF(KMAXDU.LT.5)GOTO 70
C  ���O C�E�AHO 5 �A�OB ���HO� MAXDU
C  5 steps was made with length=MAXDU
      TERMCD=5
      IF(KPRSOL.GT.0) WRITE(6, 105) ITER,MAXDU,CONV,EPSSOL
*                      PRINT    105, ITER,MAXDU,CONV,EPSSOL
      RETURN
  70  CONTINUE

      KMAXDU=0
C   O�PE�E�EH�E OTHOC�TE��HO� CKOPOCT� �MEH��EH��  FNOR
C   define relative "speed" of decreasing of FNOR
      TERM=0.D0


      DO 80 I=1,N
      SINV=1/SX(I)
      TERM=DMAX1(TERM,DABS(G(I))*DMAX1(X(I)      ,SINV)/                
     + DMAX1(FNOR,DFLOAT(N)/2.D0))
C     ��A�EH PRINT 222 C COOTB. FORMATOM 8.12.89. KO�MAHOBA H.B.  *****
C     PRINT 222 removed Koshmanova N.V.
  80  CONTINUE
      IF(TERM.GT.EPSMIN)GOTO 90
C  �TAK,M� HAXO��MC� B �OKA��HOM M�H�M�ME FNOR.
C  here wa are at lockal minimum of FNOR
      TERMCD=6
      IF(KPRSOL.GT.0) PRINT    106, ITER,TERM,EPSMIN,FNOR,CONV,EPSSOL
      RETURN
  90  CONTINUE
C   H��E�O HE �PO��O��O
C   nothing happened
      IF(KPRSOL.GE.2) PRINT    107, CONV,DUNOR,TERM

      RETURN
C 101  FORMAT('     HA ',I5,' �TEPA��� �OCT��H�TA CXO��MOCT�'/
C     *       '  K PE�EH�� C HEB��KO� <= ',E12.6 ,' ( < ',E12.6,' )'/
C     *       '  �O�OB�HA KBA�PATA L-2 HOPM� HEB��K� =',E12.6 )
 101  FORMAT('SUCCESS.'/
     *       '  AT ',I5,'-TH ITERATION CONVERGED TO SOLUTION    '/
     *       '  WITH RESIDUAL <= ',E12.6 ,' ( < ',E12.6,' )'/
     *       '  1/2 OF SQUARED L-2 NORM OF RESIDUAL =',E12.6 )
C 102  FORMAT('     HA ',I5,' �TEPA��� MAX.�A�=',E12.6,' ( <',E12.6,')'/
C     *       '  HEB��KA=',E12.6/
C     *       '  1/2 KBA�PATA HOPM� HEB��K� =',E12.6)
 102  FORMAT('FAILURE IS POSSIBLE. CHECK RESIDUAL.'/
     *       '  AT ',I5,' ITERATION MAX. STEP =',E12.6,' ( <',E12.6,')'/ 
     *       '  RESIDUAL =',E12.6/
     *       '  1/2 SQUARED NORM OF RESIDUAL =',E12.6)
C 103  FORMAT('     HA ',I5,' �TEPA��� HE��A�OC� �O����T� �P�EM�EM��'/
C     *       '  �A� > ',E12.6/
C     *       '  1/2 KBA�PATA HOPM� HEB��K� =',E12.6)
 103  FORMAT('FAILURE IS POSSIBLE. CHECK RESIDUAL.'/
     *       '  AT ',I5,' ITERATION IMPOSSIBLE OBTAIN'/   
     *       '  ALLOWABLE STEP > ',E12.6 /       
     *       '  1/2 SQUARED NORM OF RESID. =',E12.6)
C 104  FORMAT('     �B�.  �C�EP�AH ��M�T:',I4,' �TEPA���.'/
C     *       '  HEB��KA =',E12.6,'( > ',E12.6,' );'/
C     *       '  �A� =',E12.6,'( > ',E12.6,' ).')
 104  FORMAT('FAILURE.'/
     *       '  EXHAUSTED LIMIT ',I4,' OF ITERATIONS'/
     *       '  RESIDUAL =',E12.6,'( > ',E12.6,' );'/
     *       '  STEP =',E12.6,'( > ',E12.6,' ).')
C 105  FORMAT('     CTPAHHO.. BCE�O ',I4,' �TEPA��� , A ��E C�E�AHO'/
C     *       '  5 �A�OB ���HO� ',E12.6/
C     *       '  �P��EM HEB��KA =',E12.6,' ( > ',E12.6,' ).')
 105  FORMAT('FAILURE.'/
     *       '  AT ',I4,'-TH ITERATION'/
     *       '  5 STEPS WITH LENGTH ',E12.6,' WAS MADE.'/
     *       '  RESIDUAL =',E12.6,' ( > ',E12.6,' ).')
C 106  FORMAT('  M� �O�A�� B �OKA��H�� M�H�M�M :  '/
C     *       '              �TEPA���                         ',I4/
C     *       '              MEPA �PA�. HOPM� HEB��K�      ',E12.6,
C     *                                               '(<',E12.6 ,')'/
C     *       '              1/2 KBA�PATA L-2 HOPM� HEB��K�',E12.6/
C     *       '               HEB��KA                      ',E12.6/
C     *       '               TPE�OBAH�E                   ',E12.6)
 106  FORMAT('FAILURE.'/
     *      '    LOCAL MINIMUM OF 1/2*|F|*|F|  :  '/ 
     *      '     ITERATION                     ',I4/
     *      '     MEASURE OF GRADIENT OF RESID. ',E12.6,
     *      '(<',E12.6 ,')'/
     *      '     1/2 SQUARED L-2 NORM OF RESID.',E12.6/
     *      '     RESIDUAL                      ',E12.6/
     *      '     REQUIRED                      ',E12.6)

C 107  FORMAT('   HEB��KA=',E12.6,',  �A�=',E12.6,',  CKOP.�MEH��.HEB.='
C     *       , E12.6)
 107  FORMAT('  RESIDUAL=',E12.6,', STEP=',E12.6,
     *',  RESIDUAL DECR.  =', E12.6)

C     DEBUG SUBTRACE,INIT
      END
