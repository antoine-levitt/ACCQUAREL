MODULE rootfinding_tools
CONTAINS

SUBROUTINE scrsho(func)
! Reference: W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling, Numerical Recipes in Fortran 90: the art of parallel scientific computing, Cambridge University Press, 1986-1996.
! For interactive "dumb terminal" use. Produce a crude graph of the function func over the prompted-for interval x1,x2. Query for another plot until the user signals satisfaction.
! Parameters: Number of horizontal and vertical positions in display.
  IMPLICIT NONE
  INTERFACE
    DOUBLE PRECISION FUNCTION func(x)
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN) :: x
    END FUNCTION func
  END INTERFACE
  INTEGER,PARAMETER :: ISCR=60,JSCR=21
  INTEGER :: i,j,jz
  DOUBLE PRECISION :: dx,dyj,x,x1,x2,ybig,ysml
  DOUBLE PRECISION,DIMENSION(ISCR) :: y
  CHARACTER(1),DIMENSION(ISCR,JSCR) :: scr
  CHARACTER(1) :: blank=' ',zero='-',yy='l',xx='-',ff='x'

  do
     write(*,*) ' Enter x1,x2 (= to stop)' ! Query for another plot; quit if x1=x2.
     read(*,*) x1,x2
     if (x1==x2) RETURN
     scr(1,1:JSCR)=yy ! Fill vertical sides with character 'l'.
     scr(ISCR,1:JSCR)=yy
     scr(2:ISCR-1,1)=xx ! Fill top, bottom with character '-'.
     scr(2:ISCR-1,JSCR)=xx
     scr(2:ISCR-1,2:JSCR-1)=blank ! Fill interior with blanks.
     dx=(x2-x1)/(ISCR-1)
     x=x1
     do i=1,ISCR ! Evaluate the function at equal intervals.
        y(i)=func(x)
        x=x+dx
     end do
     ysml=min(minval(y(:)),0.D0) ! Limits will include 0.
     ybig=max(maxval(y(:)),0.D0)
     if (ybig==ysml) ybig=ysml+1.D0 ! Be sure to separate top and bottom.
     dyj=(JSCR-1)/(ybig-ysml)
     jz=1-ysml*dyj  ! Note which row corresponds to 0.
     scr(1:ISCR,jz)=zero
     do i=1,ISCR ! Place an indicator at function height and 0.
        j=1+(y(i)-ysml)*dyj
        scr(i,j)=ff
     end do
     write (*,'(1x,1p,e10.3,1x,80a1)') ybig,(scr(i,JSCR),i=1,ISCR)
     do j=JSCR-1,2,-1 ! Display.
        write (*,'(12x,80a1)') (scr(i,j),i=1,ISCR)
     end do
     write (*,'(1x,1p,e10.3,1x,80a1)') ysml,(scr(i,1),i=1,ISCR)
     write (*,'(12x,1p,e10.3,40x,e10.3)') x1,x2
  end do
END SUBROUTINE scrsho

FUNCTION rtbis(func,x1,x2,xacc)
! Reference: W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling, Numerical Recipes in Fortran 90: the art of parallel scientific computing, Cambridge University Press, 1986-1996.
! Using bisection, find the root of a function func known to lie between x1 and x2. The root, returned as rtbis, will be refined until its accuracy is ±xacc.
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: x1,x2,xacc
  DOUBLE PRECISION :: rtbis
  INTERFACE
    DOUBLE PRECISION FUNCTION func(x)
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN) :: x
    END FUNCTION func
  END INTERFACE
  INTEGER,PARAMETER :: MAXIT=1000 ! Maximum allowed number of bisections.
  INTEGER :: j
  DOUBLE PRECISION :: dx,f,fmid,xmid

  fmid=func(x2)
  f=func(x1)
  if (f*fmid>=0.D0) STOP'Function RTBIS: root must be bracketed.'
  if (f<0.D0) then ! Orient the search so that f>0 lies at x+dx.
     rtbis=x1
     dx=x2-x1
  else
     rtbis=x2
     dx=x1-x2
  end if
  do j=1,MAXIT ! Bisection loop.
     dx=dx*0.5D0
     xmid=rtbis+dx
     fmid=func(xmid)
     if (fmid<=0.D0) rtbis=xmid
     if (abs(dx)<xacc.or.fmid==0.D0) RETURN
  end do
  STOP'Function RTBIS: too many bisections.'
END FUNCTION rtbis

FUNCTION rtsafe(funcd,x1,x2,xacc)
! Reference: W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling, Numerical Recipes in Fortran 90: the art of parallel scientific computing, Cambridge University Press, 1986-1996.
! Using a combination of Newton-Raphson and bisection, find the root of a function bracketed between x1 and x2. The root, returned as the function value rtsafe, will be refined until its accuracy is known within ±xacc. funcd is a user-supplied subroutine that returns both the function value and the first derivative of the function.
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: x1,x2,xacc
  DOUBLE PRECISION :: rtsafe
  INTERFACE
    SUBROUTINE funcd(x,fval,fderiv)
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN) :: x
      DOUBLE PRECISION,INTENT(OUT) :: fval,fderiv
    END SUBROUTINE funcd
  END INTERFACE

  INTEGER,PARAMETER :: MAXIT=100 ! Maximum allowed number of iterations.
  INTEGER :: j
  DOUBLE PRECISION :: df,dx,dxold,f,fh,fl,temp,xh,xl

  call funcd(x1,fl,df)
  call funcd(x2,fh,df)
  write(55,*)x1,fl,x2,fh
  if ((fl>0.D0.and.fh>0.D0).or.(fl<0.D0.and.fh<0.D0)) &
  & STOP'Function RTSAFE: the root must be bracketed.'
  if (fl==0.D0) then
     rtsafe=x1
     RETURN
  else if (fh==0.D0) then
     rtsafe=x2
     RETURN
  else if (fl<0.D0) then ! Orient the search so that f(xl)<0.
     xl=x1
     xh=x2
  else
     xh=x1
     xl=x2
  end if
  rtsafe=.5D0*(x1+x2) ! Initialize the guess for root,
  dxold=abs(x2-x1)    ! the “stepsize before last,”
  dx=dxold            ! and the last step.
  call funcd(rtsafe,f,df)
  do j=1,MAXIT ! Loop over allowed iterations.
     if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f)>0.D0.or.abs(2.D0*f)>abs(dxold*df)) then ! bisect if Newton out of range, or not decreasing fast enough.
        dxold=dx
        dx=.5D0*(xh-xl)
        rtsafe=xl+dx
        if (xl==rtsafe) RETURN ! Change in root is negligible.
     else                      ! Newton step acceptable. Take it.
        dxold=dx
        dx=f/df
        temp=rtsafe
        rtsafe=rtsafe-dx
        if (temp==rtsafe) RETURN
     end if
     if (abs(dx)<xacc) RETURN ! Convergence criterion.
     call funcd(rtsafe,f,df) ! One new function evaluation per iteration.
     write(55,*)rtsafe,f,df
     if (f<0.D0) then ! Maintain the bracket on the root.
        xl=rtsafe
     else
        xh=rtsafe
     end if
  end do
  STOP'Function RTSAFE: exceeded the maximum of iterations.'
END FUNCTION rtsafe
END MODULE
