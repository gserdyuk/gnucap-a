c auxilliary library for 
c vector operations and other sspecial operations for nonlinear equations
c
c
c l2normsc - l2 norm scaled
      function l2normsc(N,a,s)
      implicit none
      integer N
      double precision l2normsc
      double precision a(N),s(N)
      integer i
      double precision tmp
      l2normsc=0
      do i=1,N
        tmp=s(i)*a(i)
        l2normsc=l2normsc+tmp*tmp
      enddo
      l2normsc=sqrt(l2normsc)
      return
      end
      

c linfnormsc - l-inf norm scaled
      function linfnormsc(N,a,s)
      implicit none
      integer N
      double precision linfnormsc
      double precision a(N),s(N)
      integer i
      double precision tmp
      linfnormsc=0
      do i=1,N
        tmp=abs(s(i)*a(i))
        linfnormsc=max(tmp*tmp,linfnormsc)
      enddo
      return
      end

c     subroutine clean            
      subroutine CLEAN(N,X)
      integer N
      double precision X(N)
      do i=1,N
        X(i)=0
      enddo
      return 
      end

c**********************************************************
c function to calculate machine precision (epsilon)
      function macheps()
      implicit none
      double precision macheps, tmp
      macheps=1d0
      tmp=1.d0+macheps
      do while (tmp.ne.1.d0)
        macheps=macheps/2d0
        tmp=1.d0+macheps
      end do
      macheps=2*macheps
      return
      end
c**********************************************************                        
      subroutine fdjac(N,x,f,NEF,pp,sx,teta,jac)      
c N     size
c x     variables
c f     vector-function
c pp    pointer to parameters (for C)
c pj    add parameters for jacobian
c npj   their number
c nef   name of program for nonlinear function
c sx    scalinf by x
c teta  precision
c jac   jcobian itself
      implicit none
      integer N,pp
      double precision x(N),f(N),sx(N),teta
      double precision fj(N),jac(N,N)
      external NEF
      double precision sqrteta, stepsizej, tempj
      integer i,j,icode
      sqrteta=sqrt(teta)
      do j=1,N 
        stepsizej=sqrteta*max(abs(x(j)),1./sx(j))*sign(1,x(j))
        tempj=x(j)
        x(j)=x(j)+stepsizej
        stepsizej=x(j)-tempj
        call NEF(N,x,fj,icode,pp)
            do i=1,N
                jac(i,j)=(fj(i)-f(i))/stepsizej
            enddo
        x(j)=tempj
      enddo
      return
      end
      
c**********************************************************      
      subroutine gradjac(N,jac,sf,fv,grad)
      implicit none
      integer N,i,j
      double precision jac(N,N),sf(N),fv(N), grad(N)
      do i=1,N
        grad(i)=0
        do j=1,N
            grad(i)=grad(i)+jac(j,i)*fv(j)*sf(j)**2
        enddo
      enddo
      return 
      end
      
c**********************************************************      
c  ########## to be completed
      subroutine FDGRAD(N,sf,sx,NEF,pF,grad)
      integer N,pF
      double precision sf(N), sx(N), grad(N)
      external NEF
      do i=1,N
        grad(i)=0
      enddo
      return 
      end
c**********************************************************                        
      subroutine fdjacdc(N,x,f,NEF,pp,sx,teta,jac,M)      
c computes finite difference jacobian for DC component only. To be used 
c as block preconditioner for inexact solver
c N     size (full size of system)
c x     variables
c f     vector-function
c pp    pointer to parameters (for C)
c pj    additional  parameters for jacobian
c npj   their number
c nef   name of program for nonlinear function
c sx    scalinf by x
c teta  precision
c jac   jcobian itself
c M     half of the number of dc components
c let nfc - number of frequencies = N/(2*M). Then both x and f carry interesting us values at position:
c 1,2*nfc+1,4*nfc+1,...2*i*nfc+1 - so RE components and at DC only :)
c
c  thus, this jacobian is df_j/dx_i, where i,j=1,2nfc+1,4nfc+1....
c
c
      implicit none
      integer N,pp,M
      double precision x(N),f(N),sx(N),teta
      double precision fj(N),jac(M,M)
      external NEF
      double precision sqrteta, stepsizej, tempj
      integer i,j,icode,j1,i1,nfc
      sqrteta=sqrt(teta)
      nfc=N/(2*M)
*       print *,'N=',N,' M=',M, ' nfc=',nfc
      do j=1,M 
c     so -  for real part only
        j1=2*(j-1)*nfc+1
*        print *,'j=',j,' j1=',j1
        stepsizej=sqrteta*max(abs(x(j1)),1./sx(j1))*sign(1,x(j1))
        tempj=x(j1)
        x(j1)=x(j1)+stepsizej
        stepsizej=x(j1)-tempj
        call NEF(N,x,fj,icode,pp)
            do i=1,M
*        print *,'I=',i,' I1=',i1
                i1=2*(i-1)*nfc+1
                jac(i,j)=(fj(i1)-f(i1))/stepsizej
            enddo
        x(j1)=tempj        
        
      enddo
      return
      end
      