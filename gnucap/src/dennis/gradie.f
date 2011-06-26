      SUBROUTINE GRADIE(NTOT,N,DJ,F,SF,GR)
C  calculate gradient from Jacobian and vector of residuals F
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION F(1),SF(1)
      DOUBLE PRECISION GR(1)
C$LARGE: DJ
      DOUBLE PRECISION DJ(NTOT,1)

C ðåþáôáôø ÷èïä
*      print * ,'$$$$$$$$$$$$$$$       GRADIE IN    $$$$$$$$$$$$$$'
*      print * ,' N= ',N,' NTOT=', NTOT
*
*      print * ,' ','       F     ', '  ', '      G     '
*      do i=1,N
*      print 933, F(I), GR(I)
*      enddo
*933   format (1x, 2(E12.6,2x))  
*
*      print * ,' ','       SF     '
*      do i=1,N
*      print 9331, SF((I+1)/2)
*      enddo
*9331   format (1x, (E12.6,2x))  

*      do i=1,2
*      do j=1,N
*      print * ,'DJ(', J,' , ',I, ')=', DJ(J,I)
*      enddo
*      enddo
      
      DO 10 I=1,N
      GR(I)=0.D0
      DO 10 J=1,N
      GR(I)=GR(I)+DJ(J,I)*F(J)*SF(J)**2
   10 CONTINUE
C  äìñ äBõX COCEäHéX üìEMEHTOB ðPéMEHñETCñ OäéH é TOT öE SCALE,
C  T.K.CõTø RE é IM OäHOçO þéCìA.
C  for two heighbour elements one SCALE is used, because RE and IM
C  are parts of the same value

*     ðåþáôø
*      print * ,'$$$$$$$$$$$$$$$       GRADIE OUT    $$$$$$$$$$$$$$'
*      print * ,' N= ',N,' NTOT=', NTOT
*      print * ,' ','       F     ', '  ', '      G     '
*      do i=1,N
*      print 933, F(I), GR(I)
*      enddo

      RETURN

C     DEBUG SUBTRACE,INIT(GR)
      END
