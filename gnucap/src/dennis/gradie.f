      SUBROUTINE GRADIE(NTOT,N,DJ,F,SF,GR)

C  calculate gradient from Jacobian and vector of residuals F
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION F(1),SF(1)
      DOUBLE PRECISION GR(1)
C$LARGE: DJ
      DOUBLE PRECISION DJ(NTOT,1)

      
      DO 10 I=1,N
      GR(I)=0.D0
      DO 10 J=1,N
      GR(I)=GR(I)+DJ(J,I)*F(J)*SF(J)**2
   10 CONTINUE


      RETURN 

C     DEBUG SUBTRACE,INIT(GR)
      END
