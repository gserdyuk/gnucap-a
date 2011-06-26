      SUBROUTINE LSEARCH(N,X,FNOR,G,Y,SX,SF,IRETCD,MAXTKN,XN,FN,  
     +                FNOR1,LAMBDA, FUNCT,
     +                EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT, ADD_DATA)
C
C ��HE�H�� �O�CK.   MO������POBAHH�� A��OP�TM �O���CTE�HA - APM��O
C linear search in Newton's direction. 
C Modified Goldstain-Armijo algorithm
C$LARGE:Y
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  external FUNCT
      DOUBLE PRECISION X(1),F(1),G(1),Y(1),SX(1),SF(1),XN(1),FN(1)
      INTEGER           IRETCD
      LOGICAL           MAXTKN
      DOUBLE PRECISION              MAXDU,EPSSOL,EPSDU,EPSMIN
      DOUBLE PRECISION              NEWTLN,MINLBD,LAMBDA,LTEMP,LPREV
      INTEGER ADD_DATA
C �������� ����
*      print * ,' N= ',N,' FNOR= ',FNOR
*      print * ,' ','       X     ', '  ', '      G     ','  ',' Y   '
*      do i=1,N
*      print 933, X(I), G(I),Y(I)
*      enddo
*933   format (1x, 3(E12.6,2x))     

C �H���A��� A��OP�TMA
C initialization
      MAXTKN=.FALSE.
      IRETCD=2
      ALPHA=1.D-4

C  O�PE�E�EH�E L-2 HOPM� BEKTOPA DU
C L-2 norm of DU vector
      NEWTLN=0.D0
      DO 10 I=1,N
   10 NEWTLN=NEWTLN+(Y(I)*SX(I))**2
      NEWTLN=DSQRT(NEWTLN)
      IF(NEWTLN.LE.MAXDU) GO TO 20

C H��TOHOBCK�� �A� �O���E MAXDU
C newton step larger than MAXDU
      REL=MAXDU/NEWTLN
	  print *, 'maxdu= ',maxdu, 'newtln=', newtln
	  
	  
      DO 30 I=1,N
   30 Y(I)=REL*Y(I)
      NEWTLN=MAXDU
   20 CONTINUE

C  B���C�EH�E CKOPOCT� �MEH��EH��
C  "speed" of decreasing
      SLOPE=0.D0
      DO 40 I=1,N
   40 SLOPE=SLOPE-G(I)*Y(I)

C  OTHOC�TE��HA� ���HA �A�A
C  relative step length
      RELLEN=0.D0
      DO 50 I=1,N
      SINV=1/SX(I)
   50 RELLEN=DMAX1(RELLEN,DABS(Y(I))/DMAX1(DABS(X(I)),SINV))

C  M�H�MA��HO �O��CT�MA� ���HA �A�A
C  min available step length 
      MINLBD=EPSDU/RELLEN
C  �H���A���A��� LAMBDA
C  initialize LAMBDA
      LAMBDA=1.0D0

C_ B���C��M LAMBDA_____________________________________________________
C  calculate LAMBDA
 1000 CONTINUE
      DO 65 I=1,N
   65 XN(I)=X(I)-LAMBDA*Y(I)

      CALL NEF(N,XN,FN,IFLAG, SF,FNOR1, FUNCT, ADD_DATA)
C  B���C�EH�E L-2 HOPM� F+ �POBE�EHO B NEF
C  L-2 norm of F+ calculated in NEF
C  �POBEPKA FT<=FC+ALPHA*LAMBDA*SLOPE
C  check FT<=FC+ALPHA*LAMBDA*SLOPE
      IF((FNOR1-FNOR).GT.(ALPHA*LAMBDA*SLOPE))GOTO 80
C  ECT� XOPO�A� TO�KA
C  good point have found
      IRETCD=0
C MAXTKN=?
      IF(LAMBDA.EQ.1.0D0.AND.NEWTLN.GT.0.99D0*MAXDU) MAXTKN=.TRUE.
C_BO�BPAT______________________________________________________________
C returen
      RETURN

   80 IF(LAMBDA.GE.MINLBD) GO TO 90
C  XOPO�E� TO�K� HA�T� HE��A�OC
C  there are no good points at all.
      IRETCD=1
      RETURN


   90 CONTINUE
C  COO��EH�E O HA�A�E O�HOMEPHO�O �O�CKA
C   (HA �EPBOM �A�E,KO��A LAMBDA=1)
C  start of 1-dimentional search (at the first step, LAMBDA=1)
      IF(LAMBDA.EQ.1.D0)PRINT    501

C  �MEH��AEM LAMBDA
C  decrease LAMBDA
      IF(LAMBDA.LT.1) GO TO 100

C  �EPBOE �PO��EH�E. KBA�PAT��HA� �HTEP�O�����
C  first dividing. Quadratic interpolation
      LTEMP=-SLOPE/(2*(FNOR1-FNOR-SLOPE))
      GO TO 110

C  K����ECKA� �HTEP�O�����
C  cubic interpolation
  100 CONTINUE
      DIV=1/(LAMBDA-LPREV)
      V1=FNOR1-FNOR-LAMBDA*SLOPE
      V2=FPREV1-FNOR-LPREV*SLOPE
      A=DIV*(V1/(LAMBDA**2)-V2/(LPREV**2))
      B=DIV*(-V1*LPREV/(LAMBDA**2)+V2*LAMBDA/(LPREV**2))
C     DISC=B*B-3.*A*SLOPE
C  EC�� A=0 -K��. �HTEP�O����� B�PO��AETC� B KBA�PAT��H��.
C  if A=0, cubic interpolation degrades to quadratic
      IF(A.EQ.0.D0) LTEMP=-SLOPE/(2.D0*B)
C  HEB�PO��EHHA� �HTEP�O�����
C  non-degraded interpolation
      IF(A.NE.0.D0) LTEMP=-B/(3.D0*A)+
     &DSQRT((B/(3.D0*A))**2-SLOPE/(3.D0*A))
C  �POBEPKA LTEMP>0.5*LAMBDA
C  check LTEMP>0.5*LAMBDA
      IF(LTEMP.GT.LAMBDA/2)LTEMP=LAMBDA/2.D0
  110 CONTINUE

C �EPE�P�CBOEH�E . LTEMP B���C�EHO.
C reassign. LTEMP calculated
      LPREV=LAMBDA
      FPREV1=FNOR1

C  �POBEPKA LTEMP<=0.1*LAMBDA
C  check LTEMP<=0.1*LAMBDA

      IF(LTEMP.LE.0.1D0*LAMBDA) LTEMP=0.1D0*LAMBDA
      LAMBDA=LTEMP

C  COO��EH�E O XO�E O�HOMEPHO�O �O�CKA
C  information about 1-dim search
                      PRINT    502, FNOR1,LAMBDA
C  �MEH��EH�E LAMBDA OKOH�EHO
C  decreasing of LAMBDA is finished
      IF(IRETCD.LT.2) RETURN
      GO TO 1000
C     DEBUG SUBTRACE,INIT(NEWTLN,REL,MINLBD,RELLEN,FNOR1,FNOR,SLOPE,
C    *              DISC, LTEMP,LAMBDA,LPREV,DIV,A,B)

  501 FORMAT(15X,'1-DIMENSIONAL SEARCH: ')
  502 FORMAT(15X,'   FNOR=',E12.6,',  LAMBDA=',E12.6)



      END
