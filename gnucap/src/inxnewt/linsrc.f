c linsrcder()
c linear search routine. using 2 point (plus derivative) for 
c quadratic approx or 3-point (plus derivative) for qubic
c see Dennis Schnabe Numerical Methods for Unconstrained Optimization
c and Nonlinear Equations - alg. A6.3.1
c searches from "xc" in direction "p" for point "xc+ld*p" where
c f(x+) <= f(xc)+alpha*ld*transp(grad)*p
c alpha = 1.e-4

      subroutine linsrcder(N,xc,fc,NEF,grad,p,sx,maxstep,steptol,
     +                     retcode,xplus,fplus,fvplus,maxtaken,
     +                     ptr,lambda,sf)
c
c input parameters 
c N     size
c xc(N)     current point
c fc        current scalar value of function
c NEF       function name
c grad(N)   gradient of fc      
c p(N)      step (diretion of search)
c sx(N)     x scale
c maxstep   max tep length
c steptol   step tolerance
c output parameters
c retcode   return code 0 - OK, 1 - search failed, 2- continue work 
c xplus(N)  next poijnt (searched)      
c fplus     scalar: l2norm(NEF(xplus))
c fvplus    NEF(xplus)
c maxtaken  is max step taken? 1- yes, 0 - no
c       pf        parameters for NEF
c       npf       their number
c - removed voth pf and npf. int ptr is used instead. points to special strcuture (C)
c ptr       structure pointer made as int
c lambda    damping factor (info)
c sf(N)     scale of F    
c
      implicit none
      integer N,ptr,retcode,maxtaken
      double precision xc(N),fc,grad(N),p(N),sx(N),maxstep,steptol
      double precision xplus(N),fplus,fvplus(N),lambda,sf(N)
      double precision l2normsc,ddot
      external NEF,l2normsc,ddot

      double precision alpha,newtlen,initslope,rellength,minlambda
      parameter (alpha=1.e-4)
      integer i, icode
      double precision ltemp,k,t11,t12,t21,t22,u1,u2,aa,bb,lprev
      double precision fplusprev,disc
      
      maxtaken=0
      retcode=2
            
      newtlen=l2normsc(N,p,sx)
      if(newtlen.ge.maxstep) then
        call dscal(N,maxstep/newtlen,p,1)
        newtlen=maxstep
      endif
      initslope=ddot(N,grad,1,p,1)
*####      print *,'grad ',grad
*####      print *,'p ',p
      
*####      print *,'INITSLOPE',initslope
      rellength=0
      do i=1,N
        rellength=max(rellength,p(i)/max(abs(xc(i)),1/sx(i))) 
      enddo  
      minlambda=steptol/rellength
      
      lambda=1
      do while (retcode.eq.2)
*####        print *, 'LINSRCDER retcode=',retcode,' lambda=',lambda
        call dcopy(N,xc,1,xplus,1)
        call daxpy(N,lambda,p,1,xplus,1)
        call NEF(N,xplus,fvplus,icode,ptr)
        fplus=l2normsc(N,fvplus,sf)
*####        print *, 'FPLUS=', fplus
*####        print *, 'fc =',fc
*####        print *, ' alpha =', alpha
*####        print *, ' lambda =', lambda
*####        print *,'alpha*lambda*initslope ',alpha*lambda*initslope
*        if(fplus.le.(fc+alpha*lambda*initslope)) then
        if(fplus.le.fc*(1-alpha*lambda)) then
*####        print *,' #*********** COMPARE',-initslope,' ',fc*fc,' ',
*####     $  ' ',initslope+fc*fc
     
            retcode=0
            if(lambda.eq.1.and.newtlen.ge.0.99*maxstep) then
                maxtaken=1
            endif
        elseif(lambda.lt.minlambda) then
c ... no appropriate lambda in direction p
            retcode=1
            call dcopy(N,xc,1,xplus,1)
        else
c ... divide step 
            if(lambda.eq.1) then
c quadratic approximation            
*####                print *,' FIRST divide ###'
                ltemp=initslope/(2*(fplus-fc-initslope))
*####                print *,'initslope,fplus,fc=',initslope,' ',fplus,' ',fc
*####                print *,' ltemp=',ltemp
            else
*####                print *,' SECOND etc. divide ###'

c qubic approximation
c 
c |aa|     | t11 t12 |   | u1 |
c |  | = k*|         | * |    |
c |bb|     | t21 t22 |   | u2 |
c
*####                print *,'lprev: ',lprev
                
                k   =   1./(lambda-lprev)
                t11 =   1./(lambda**2)
                t12 =   -1./(lprev**2)
                t21 =   -lprev/(lambda**2)
                t22 =   lambda/(lprev**2)
                u1  =   fplus-fc-lambda*initslope
                u2  =   fplusprev-fc-lprev*initslope
                aa  =   k*(t11*u1+t12*u2)
                bb  =   k*(t21*u1+t22*u2)
                
*###                print *,'vals1: ',k,' ',t11,' ',t12,' ',t21,' ',t22
*###                print *,'vals2: ',u1,' ',u2,' ',aa,' ',bb
                
                disc=bb**2-3*aa*initslope
                
                if(aa.eq.0) then
                    ltemp=-initslope/(2*bb)
                else
                    ltemp=(-bb+sqrt(disc))/(3*aa)
                endif
*###                print *,'ltemp= ',ltemp
            endif
            lprev=lambda
            fplusprev=fplus
            if(ltemp.le.0.1*lambda) then
                lambda=0.1*lambda
            else
                lambda=ltemp
            endif    
        endif
      enddo  
      
      return
      end
      
c********************************************************
c linsrcdiv2()
c linear search by division by 2
c searches from "x" in direction dx until
c f(x+) <= f(xc)*(1+alpha*ld
c see C.T.Kelley "Iterative methods for Linear and Nonlinear Equations"
c SIAM Frontiers in Applied Math. 16 - condition 8.1
      subroutine linsrcdiv2(N,x,xp,dx,NEF,ptr,fvp,fplus,sf,
     +                      fcurr,retcode)
      implicit none
      integer N,ptr,retcode, icode
      double precision x(N),xp(N),dx(N),fvp(N),fplus,sf(N)
      external NEF, l2normsc
      double precision fcurr,l2normsc,lambda,lambdamin, alpha
      parameter(lambdamin=1.e-12)
      parameter(alpha=1.e-4)
      retcode=0
      lambda=1
      
*####      print *,'  ---                fcurr=',fcurr
      
      call dcopy(N,x,1,xp,1)                        
      call daxpy(N,lambda,dx,1,xp,1)
      call NEF(N,xp,fvp,icode,ptr)      
      fplus=l2normsc(N,fvp,sf)                        
*      print *,'upd:   x:',x
*      print *,'upd:   dx:',dx
*      print *,'upd:   xp:',xp
*      print *,'upd:   fvp:',fvp
*      print *,'upd:   fplus:',fplus

c  |f_plus| < (1-a*l)*|f_curr|  -  cond 8.1 from Kelley, SIAM Frontiers 16
      do while (fplus.ge.(1-alpha*lambda)*fcurr.and.retcode.eq.0)
        lambda=lambda/2
        call dcopy(N,x,1,xp,1)                        
        call daxpy(N,lambda,dx,1,xp,1)
        call NEF(N,xp,fvp,icode,ptr)
        fplus=l2normsc(N,fvp,sf)

*      print *,'upd:   fplus:',fplus
*      print *,'upd:   fvp:',fvp
        
        if(lambda.le.lambdamin) then
            retcode=1
        endif
        print *,'  --- lambda=',lambda,' fplus=',fplus
      enddo
      return
      end

c********************************************************
c linsrc3p()
c linear search routine. using 3 point (without derivative) for 
c quadratic approx 
c see C.T.Kelley "Iterative methods for Linear and Nonlinear Equations"
c SIAM Frontiers in Applied Math. 16
c searches from "xc" in direction "p" for point "xc+ld*p" where
c f(x+) <= f(xc)*(1-alpha*ld)
c alpha = 1.e-4

      subroutine linsrc3p(N,xc,fc,NEF,p,sx,maxstep,steptol,
     +                     retcode,xplus,fplus,fvplus,maxtaken,
     +                     ptr,lambda,sf)
c
c input parameters 
c N     size
c xc(N)     current point
c fc        current scalar value of function
c NEF       function name
c             grad(N)   - is absent !! 
c p(N)      step (diretion of search)
c sx(N)     x scale
c maxstep   max tep length
c steptol   step tolerance
c output parameters
c retcode   return code 0 - OK, 1 - search failed, 2- continue work , 3 - unlimited loop in step division
c xplus(N)  next poijnt (searched)      
c fplus     scalar: l2norm(NEF(xplus))
c fvplus    NEF(xplus)
c maxtaken  is max step taken? 1- yes, 0 - no
c       pf        parameters for NEF
c       npf       their number
c ptr       is used instead of pf and npf, INTEGER
c lambda    damping factor (info)
c sf(N)     scale of F    
c
      implicit none
      integer N,ptr,retcode,maxtaken
      double precision xc(N),fc,p(N),sx(N),maxstep,steptol
      double precision xplus(N),fplus,fvplus(N),lambda,sf(N)
      double precision l2normsc,ddot
      external NEF,l2normsc,ddot

      double precision alpha,newtlen,rellength,minlambda
      parameter (alpha=1.e-5)
      integer i,icode
      double precision ltemp,fplusprev, lprev
      double precision sigma1
      parameter (sigma1=0.5)
      double precision kk,aa,bb
      integer count, countmax
      parameter (countmax = 15)
      maxtaken=0
      retcode=2
            
      newtlen=l2normsc(N,p,sx)
      if(newtlen.ge.maxstep) then
        call dscal(N,maxstep/newtlen,p,1)
        newtlen=maxstep
      endif

c      print *,"linsrc3p"
      
c      print *,"+++++++++++++ p= ",p
c      print *,"+++++++++++++ xc= ",xc
c      print *,"+++++++++++++ sx= ",sx      
            
      rellength=0
      do i=1,N
        rellength=max(rellength,abs(p(i))/max(abs(xc(i)),1/sx(i))) 
      enddo  
      minlambda=steptol/rellength
*###      print *,"steptol= ",steptol, " rellength=",rellength
      
      lambda=1
      count=0
      do while (retcode.eq.2)
*###        print *,'while entered, lambda=',lambda
        call dcopy(N,xc,1,xplus,1)
        call daxpy(N,lambda,p,1,xplus,1)
        call NEF(N,xplus,fvplus,icode,ptr)
        fplus=l2normsc(N,fvplus,sf)
        if(fplus.lt.fc*(1-alpha*lambda)) then
*###            print *,'retcode <= 0'
            retcode=0
            if(lambda.eq.1.and.newtlen.ge.0.99*maxstep) then
                maxtaken=1
            endif
        elseif(lambda.lt.minlambda) then
c ... no appropriate lambda in direction p
*###            print *,'retcode <=1'
            retcode=1
            call dcopy(N,xc,1,xplus,1)
        else
*###            print *,'  divide step:'
c ... divide step 
            if(lambda.eq.1) then
c just divide step
                ltemp=sigma1*lambda
            else
c qubic 3-point apprximtion
c lambda - current value of lambda (l_c), lprev - (l_minus) - previous lambda
c fplus=F(l_c); frevplus = F(lprev); fc=F(when lambda=0)
                kk=lambda*lprev*(lambda-lprev)
                aa=1./kk*(lprev*(fplus-fc)-lambda*(fplusprev-fc))   
*###            print *,'  --- lambda!=1, aa=',aa

                if(aa.gt.0) then
c ok, second derivative is positive, quadratic interpolation is ok,                 
c polynom(lambda) is: aa*lambda**2+bb*lambda+cc=0
*###                    print *,' aa >0 - usual step division'
                    bb=1./kk*
     +                  (lambda**2*(fplusprev-fc)-lprev**2*(fplus-fc))
                    ltemp=-bb/(2.*aa)
*###                print *,'ltemp=',ltemp
*###                print *,'aa=',aa ,'bb=',bb
                else
*###                    print *,'aa <0 - no interpolation; simply div. step'
c negative or zero second derivative - cant make quadratic interpolation. simply divide step
                    ltemp=sigma1*lambda
                endif
            count=count+1
            if (count. gt. countmax) retcode=3
            endif
c update points
            lprev=lambda
            fplusprev=fplus
            if(ltemp.le.0.1*lambda) then
                lambda=0.1*lambda
            elseif (ltemp.ge.0.5*lambda) then
                lambda=0.5*lambda
            else
                lambda=ltemp
            endif    
        endif
        print *,' -- lambda next=',lambda,' fplus=',fplus,
     +          ' retcode=',retcode
      enddo  
      
      return
      end
           