      SUBROUTINE MACHEP   (EPSIM)

C   ð/ð BùþéCìEHéñ MAûéHHOçO üðCéìOH.
C   subr. to compute machine epsilon

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION EPSIM, DS
C      print 5
C    5 format(2x,'subr MACHEP')

      EPSIM=1.D0
   10 CONTINUE
      EPSIM=EPSIM/2.D0

*      print 15, epsim
      
      DS=EPSIM+1.D0
      IF(DS.NE.1.D0) GOTO 10
           
      EPSIM=2.D0*EPSIM
C      EPSIM=.1084202172E-18 
   
*      print 15, epsim
*   15 format(2x,'EPSIM=',d16.10)

      RETURN
C     DEBUG INIT(EPSIM),SUBTRACE
      END
