      SUBROUTINE solve(N, ICODE_RETURN, DFDX, U, F, 
     +                 EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT, 
     +                 FUNCT, kprsol, add_data)
****************************************************
*         Parameters 
*N			IN: int, size 
*ICODE_RETURN	OUT: int, return code
*DFDX		INTERNAL: double prec(N,N+1), jacobian
*U			INOUT: double prec(N), vector unknowns
*F			OUT: double prec(N), vector function
*EPSSOL		IN: double prec - solution prec
*EPSDU		IN: double prec - du (step) precision
*EPSMIN		IN: double prec - local minimum estimation
*MAXDU		IN: double prec - max step
*LIMIT		IN: integer - max iterations
*FUNCT		IN: function "pointer"	 
*kprsol     IN: debug level
*add_data   IN: pointer to additional data to be passed inside FUNCT - passed via (void *)
*
*		Internbal variables
*Y			double prec(N), vector correction
*GR			double prec(N), vector gradient
*UN			double prec(N), vector unknowns, temporary
*FN			double prec(N), vector func, temporary
*SF			double prec(N), vector function scales
*SX			double prec(N), vector unknowns scales
*INDEX		integer(N), vector perestanovok
	 
C************************************************
C  driver of nonlinear simultaneous equations solving
C         using Newton's method
C
C  see:     Dennis J.,jr, Schnabel R.
C       Numerical methods of unconstrained minimization
C       and solving of nonlinear equations.
C       Translation: Moscow, MIR,1988ç.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  integer ICODE_RETURN
	  external FUNCT
 	  DOUBLE PRECISION U(N)
      DOUBLE PRECISION F(N)
      DOUBLE PRECISION DFDX(N,N+1)      
      DOUBLE PRECISION EPSSOL,EPSDU,EPSMIN,MAXDU
	  INTEGER LIMIT
      INTEGER add_data
*  allocated in this module: 	  
	  DOUBLE PRECISION Y(N),GR(N),UN(N),FN(N)
	  double precision WA1(N) 
      DOUBLE PRECISION SF(N),SX(N)
	  
C     IRANG PABEH PAúMEPHOCTé úAäAþé.
C     IRANG is equal to size of task
C change IRANG. Gennady

      INTEGER IRANG

C  the following 3 lines - for LINEQ1
      INTEGER         		INDEX(N),NER 
      DOUBLE COMPLEX         DET
      DOUBLE PRECISION       ARGD

      INTEGER KPRSOL
	  	  
      INTEGER TERMCD,KMAXDU
      DOUBLE PRECISION TS,TLS,FNOR,FNORN
      LOGICAL MAXTKN
            
      ic_func=0
      ic_jac=0
	  ICODE_RETURN=0
c initialization
      IRANG=N
      IRANG1 =IRANG+1
      N1=IRANG1
      ITER=0
      TLS=1.D0
      TERMCD=0

C   calculate machine epsilon
      CALL MACHEP(EPSIM)
	  print *, 'epsim=', epsim
C  set and check input information
      CALL NEINCK(N,EPSIM,IRANG,SF,SX,U,TERMCD,
     +            EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT)
	  ICODE_RETURN=TERMCD
      IF(TERMCD.GE.0) GOTO 10
      IF(KPRSOL.GT.0) PRINT   1010, TERMCD
      call exit(73)
c #define errTermcd 73  	/*  bad input parameters for nonlin solving */
      STOP
   10 CONTINUE

      CALL NEF (N,U,FN,IFLAG, SF, FNOR, FUNCT, add_data)
C   FLAG- if error occures in subroutines of element's models

C  check - may be we are at thesolution point ???
      CALL STOP0(N,FN,U,SF,SX,TERMCD,KMAXDU, EPSSOL,
     +           EPSDU,EPSMIN,MAXDU,LIMIT)
	  ICODE_RETURN=TERMCD
C      IF(TERMCD.EQ.1)RETURN gs: 24-dec-2003 - moved down
      
C compute finite-differential jacobian
      ML=N
	  MU=N
	  epsfcn=0
	  iflag=0

      call fdjac1_d (FUNCT,N,U,FN,DFDX,N,IFLAG,ML,MU,epsfcn,wa1,wa1, 
     *               add_data)
	  if (kprsol. ge. 4) then 
        print *,'__jacobian at iter 1__'
  	    do ii=1,N
	      print 555, (DFDX(ii,j),j=1,N)
	    enddo
	  endif
      print *,'_______________________'
      
      IF(TERMCD.EQ.1)GOTO 1100

C  compute gradient of 1/2 squared norm of residual by variable
      CALL GRADIE(IRANG,N,DFDX,FN,SF,GR)
*      print *, 'gradient: ..'
*	  print *, GR

C F <- FN
      DO I=1,N
      F(I)=FN(I)
	  END DO

C iterations
 1000 IF(TERMCD.NE.0)GOTO 1100
      call flush
      ITER=ITER+1
C prints
      IF(KPRSOL.GE.2) PRINT 5, ITER
      IF(KPRSOL.GE.3) PRINT 2, (U(I),I=1,N)
      IF(KPRSOL.GE.3) PRINT 1, (F(I),I=1,N)


C ********************************************

C  solve affine model

*      print *, 'F before LINEQ1 ******', F
      DO I=1,N
      DFDX(I,N1)=F(I)
      ENDDO
	  
      CALL LINEQ1(DFDX,IRANG,N,IRANG1,1,INDEX,NER,DET,ARGD)
C  INTEGER INDEX(IRANG),COMPLEX DET,REAL ARGD,INTEGER NER  are
C  declared already. NER=1-OK,NER=0-zero column detected
   
      DO I=1,N
      Y(I)=DFDX(I,N1)
	  ENDDO
      
*      print *, 'Y after LINEQ1 *******', Y
      IF(KPRSOL.GE.3) PRINT 3, (Y(I),I=1,N)
 
C   calculate angle between gradient and correction-vector. Estim.
C   of condition number of Jacobian. L2 norm of gradient and corection

      SLOPE=0.D0
      GRNOR=0.D0
      YNOR=0.D0
      DO II=1,N
      SLOPE=SLOPE-GR(II)*Y(II)
      GRNOR=GRNOR+GR(II)*GR(II)
 	  YNOR=YNOR+Y(II)*Y(II)
 	  ENDDO
	  
      YNOR=DSQRT(YNOR)
      GRNOR=DSQRT(GRNOR)
C  COS(GAMMA)
      COSGA=SLOPE/(GRNOR*YNOR)
      IF(KPRSOL.GE.2)  PRINT 3001, SLOPE,GRNOR,YNOR, COSGA
 3001 FORMAT(2X,'ZNEWSOL : SLOPE=',E12.5,' GRNOR=',E12.5,' YNOR=',      
     + E12.5/12X,'COSGA=',E20.14)
	  if(cosga.gt.1.d0) cosga=1.d0
	  if(cosga.lt.-1.d0) cosga=-1.d0

      DGAMMA=DACOS(COSGA)
      GA=(DGAMMA/3.14D0)*180.D0
      IF(KPRSOL.GE.2) PRINT    2003, GA
 2003 FORMAT(' * ANGLE BTWN. GRAD. AND CORR.  =',F8.4)
      ESTIM=-1.D0/COSGA
      IF(KPRSOL.GE.2) PRINT    2004, ESTIM
 2004 FORMAT(' ** LOWER ESTIM. JACOBIAN CONDIT. NUMBER =',E14.7)

C .2. linear search
*      print*, 'LSEARCH ENTERED *******'
*	  print *,(Y(I),I=1,N)
      CALL LSEARCH(N,U,FNOR,GR,Y,SX,SF,IRETCD,MAXTKN,
     +            UN,FN,FNORN,TLS, FUNCT,
     +            EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT, add_data)

*      print*, 'LSEARCH EXITED *******'
*	  print *,(Y(I),I=1,N)
*	  print *,'TLS=',TLS


C  finish iteration ...
C    Jacobian:
      call fdjac1_d (FUNCT,N,UN,FN,DFDX,N,IFLAG,ML,MU,epsfcn,wa1,wa1,
     *               add_data)
	  if (kprsol. ge. 4) then 
	    print *, '__DFDX__'
	    do i=1,N
	      print 555, (DFDX(i,j),j=1,N)
555     format(24(2x,E8.2))		  
	    enddo
      endif  
      print *, '____________________'


C    gradient:
      CALL GRADIE(IRANG,N,DFDX,FN,SF,GR)
*      print *, 'gradient: ..'
*	  print *, GR

C  check stop conditions

      CALL STOP(N,U,Y,FN,FNORN,GR,SX,SF, IRETCD,ITER, 
     +          MAXTKN,KMAXDU,TERMCD,
     +          EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT)
	  ICODE_RETURN=TERMCD
      IF(KPRSOL.GE.3) PRINT 2, (U(I),I=1,N)

C  and more prints
      IF(KPRSOL.GE.3) PRINT 4, FNOR,TS,TLS
    5 FORMAT(2X,'ITERATION ',I4)
    1 FORMAT(2X,'  F= ..'/(3X,6(E12.5)))
    2 FORMAT(2X,'  U= ..'/(3X,6(E12.5)))
    3 FORMAT(2X,' DU= ..'/(3X,6(E12.5)))
    4 FORMAT(2X,'FNOR=',E12.5,2X,'TS=',E12.5,2X,'TLS=',E12.5)


C  U <- UN , F <- FN , FNOR <- FNORN
      DO 110 I=1,N
      U(I)=UN(I)
  110 F(I)=FN(I)
      FNOR=FNORN
      IF(KPRSOL.GE.3) PRINT 2, (U(I),I=1,N)
      
      GOTO 1000
C___K_O_H_E_ã__________________________________________________________
C___E N D _____________________________________________________________
1100  CONTINUE
      IF(KPRSOL.GE.2) PRINT 1200
1200  FORMAT(2X,'ZNEWSOL (EXIT)  : ')
      IF(KPRSOL.GE.3) PRINT 2, (U(I),I=1,N)
      
      print *, ' ##STATISTICS: ',' ITERATIONS= ', ITER,
     *        ' FUNC_CALLS= ',ic_func,
     *        ' JAC_CALLS= ',ic_jac,
     *        ' TERMCODE= ',TERMCD
      
      call flush()
      
*      if(store_jac.eq.1) then
* store jacobian in file - by columns  -- for inexact solver as preconditioner
*        OPEN(UNIT=99,FILE='JACOBIAN')
*        write (unit=99,fmt=*) N
*        do j=1,N+1
*            write (unit=99,fmt=*),'--'
*            write (unit=99,fmt=*)(DFDX(i,j),i=1,N)
*        enddo
*        CLOSE(99)
*      end if

      RETURN

 1010 FORMAT(2X,'  ##### ZNEWSOL:  ATTENTION !!!            #######'/   
     +    2X,'     FATAL ERROR IN INPUT DATA            .'/       2X,'  
     +          ERROR CODE     =',I2)
      END


* wrapper for nonlinear function with norm calculation)

      SUBROUTINE NEF (N,U,FN,IFLAG, SF, FNOR, FUNCT, ADD_DATA)
  
	  external FUNCT
	  INTEGER N
	  DOUBLE PRECISION FN(1), U(1), SF(1)
	  DOUBLE PRECISION FNOR
      INTEGER ADD_DATA
	  INTEGER IFLAG
	  integer i
	  
c      print *,'NEF entered'
c	  print *, 'N=',N
c      print *, 'U=',(U(i), i=1,N)
c      print *, 'FN=',(FN(i), i=1,N)
c      print *, 'IFLAG=',IFLAG
c      print *, 'SF',(SF(i), i=1,N)
c      print *, 'FNOR=', FNOR
C      print *, FUNCT
c      call flush()
	  
	  CALL FUNCT(N,U, FN, IFLAG, ADD_DATA)

c      print *, 'NEF before norm calculation'
c      print *, 'N=',N
c      print *, 'U=',(U(i), i=1,N)
c      print *, 'FN=',(FN(i), i=1,N)
c      print *, 'IFLAG=',IFLAG
c      print *, 'SF',(SF(i), i=1,N)
c      print *, 'FNOR=', FNOR
C      print *, FUNCT
	  
	  FNOR=0.0D0
      DO I=1,N
      FNOR=FNOR+SF(I)*SF(I)*FN(I)*FN(I)
c      print *, 'FNOR[i]: ',i,' ',FNOR
	  enddo
      FNOR=FNOR/2.D0

c      print *, 'FNOR: ', FNOR
c	  print *, 'NEF exited'
c      call flush()
	  
	  return 
	  end
