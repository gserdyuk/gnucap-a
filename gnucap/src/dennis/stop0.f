      SUBROUTINE STOP0(N,F,X,SF,SX,ICODE,KMAXDU,
     +EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT)
C     check for stop conditions at the beginning of 
C     the Newton's algorithm

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION F(1), X(1)
      DOUBLE PRECISION SF(1),SX(1)
      INTEGER ICODE,KMAXDU

      DOUBLE PRECISION  EPSSOL,EPSDU,EPSMIN,MAXDU
      integer KPRSOL

      DOUBLE PRECISION CONV


      KPRSOL=1
C      �CTAHOB�T� C�ET��K KO���ECTBA
C      MAKC�MA��H�X �A�OB HA 0
C      set maxstep counter to 0
      KMAXDU=0

C      O�PE�E�EH�E HEB��K� B HA�A��HO� TO�KE
C      calculate residual in initial point
      ICODE=0
      CONV=0.0D0
      DO 10 J=1,N,1
   10 CONV=DMAX1(CONV,SF(J)*DABS(F(J)))
      IF(CONV.LE.1.D-02*EPSSOL)ICODE=1

C  �E�AT� COO��EH��
C  print message
      IF(ICODE.EQ.1.AND.KPRSOL.GT.0) WRITE(6, 100) CONV
      IF(ICODE.EQ.1.AND.markup.eq.1) then
	write(6,111)
	WRITE(6,100) CONV
        write(6,112)
      endIF

*      IF(ICODE.EQ.1)                 PRINT    100, CONV

      RETURN
C  100 FORMAT(2X,'     �CXO�HA� TO�KA �B��ETC� PE�EH�EM'/
C     *       2X,' C�CTEM� �PABHEH��. HEB��KA :',E12.6   )
  100 FORMAT(1X,'SUCCESS (FORMALLY, BUT CHECK THE SYSTEM):',
     *       1X,'INITIAL POINT IS A SOLUTION OF THE SYSTEM.'/
     *       1x,'RESIDUAL=',E12.6   )
  111 format('.@w2f1s')
  112 format('.@w2f1e')

C     DEBUG SUBTRACE,INIT
      END


