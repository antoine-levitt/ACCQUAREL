MODULE optimization_tools
CONTAINS

  SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
    ! Subroutine for initially braketing a minimum
    ! Reference: W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling, Numerical Recipes in Fortran 90: the art of parallel scientific computing, Cambridge University Press, 1986-1996.
    ! Given a function func, and given distinct initial points ax and bx, this routine searches in the downhill direction (defined by the function as evaluated at the initial points) and returns new points ax, bx, cx that bracket a minimum of the function. Also returned are the function values at the three points, fa, fb, and fc.
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(INOUT) :: ax,bx
    DOUBLE PRECISION,INTENT(OUT) :: cx,fa,fb,fc
    INTERFACE
       DOUBLE PRECISION FUNCTION func(x)
         IMPLICIT NONE
         DOUBLE PRECISION,INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE

    DOUBLE PRECISION,PARAMETER :: GOLD=1.618034D0,GLIMIT=100.D0,TINY=1.D-20 ! GOLD is the default ration by which successive intervals are magnified; GLIMIT is the maximum magnification allowed for a parabolic-fit step, TINY is used to prevent any possible division by zero.
    DOUBLE PRECISION :: dum,fu,q,r,u,ulim

    fa=func(ax)
    fb=func(bx)
    if (fb>fa) THEN ! Switch roles of a and b so that we can go downhill in the direction from a to b
       dum=ax ; ax=bx ; bx=dum
       dum=fb ; fb=fa ; fa=dum
    end if
    cx=bx+GOLD*(bx-ax) ! First guess for c.
    fc=func(cx)
    do
       if (fb<fc) RETURN
       r=(bx-ax)*(fb-fc) ! Compute u by parabolic extrapolation from a, b, c.
       q=(bx-cx)*(fb-fa)
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
       ulim=bx+GLIMIT*(cx-bx) ! We won't go farther than this. Test various possibilities:
       if ((bx-u)*(u-cx)>0.D0) THEN ! Parabolic u is between b and c: try it.
          fu=func(U)
          if (fu<fc) then ! Got a minimum between b and c.
             ax=bx
             fa=fb
             bx=u
             fb=fu
             RETURN
          else if (fu>fb) then ! Got a minimum between A and u.
             cx=u
             fc=fu
             RETURN
          end if
          u=cx+GOLD*(cx-bx) ! Parabolic fit was no use. Use default magnification.
          fu=func(u)
       else if ((cx-u)*(u-ulim)>0.D0) then ! Parabolic fit is between c and its allowed limit.
          fu=func(u)
          if (fu<fc) then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
             fb=fc
             fc=fu
             fu=func(u)
          end if
       else if ((u-ulim)*(ulim-cx)>=0.D0) THEN ! Limit parabolic u to maximum allowed value.
          u=ulim
          fu=func(u)
       else ! Reject parabolic u, use default magnification.
          u=cx+GOLD*(cx-bx)
          fu=func(u)
       end if
       ax=bx ; bx=cx ; cx=u ! Eliminate oldest point and continue.
       fa=fb ; fb=fc ; fc=fu
    end do
  END SUBROUTINE mnbrak

  FUNCTION golden(ax,bx,cx,func,tol,xmin)
    ! Golden Section Search
    ! Reference: W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling, Numerical Recipes in Fortran 90: the art of parallel scientific computing, Cambridge University Press, 1986-1996.
    ! Given a function func, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is between ax and cx, and func(bx) is less than both func(ax) and func(cx)), this routine performs a golden section search for the minimum, isolating it to a fractional precision of about tol. The abscissa of the minimum is returned as xmin, and the minimum function value is returned as golden, the returned function value.
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN) :: ax,bx,cx,tol
    DOUBLE PRECISION,INTENT(OUT) :: xmin
    DOUBLE PRECISION :: golden
    INTERFACE
       DOUBLE PRECISION FUNCTION func(x)
         IMPLICIT NONE
         DOUBLE PRECISION,INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE

    DOUBLE PRECISION,PARAMETER :: R=0.61803399D0,C=1.D0-R ! The golden ratios.
    DOUBLE PRECISION :: f1,f2,x0,x1,x2,x3

    x0=ax ! At any given time we will keep trace of four points: x0,x1,x2,x3.
    x3=cx
    if (abs(cx-bx)>abs(bx-ax)) then ! Make X0 to X1 the smaller segment, and fill in the new point to be tried.
       x1=bx
       x2=bx+C*(cx-bx)
    else
       x2=bx
       x1=bx-C*(bx-ax)
    end if
    f1=func(x1)
    f2=func(x2) ! The initial function evaluations. Note that we never need to evaluate the function at the original endpoints.
    DO ! Do-while-loop: We keep returning here.
       if (abs(x3-x0) <= tol*(abs(x1)+abs(x2))) exit
       if (f2<f1) then                 ! One possible outcome,
          x0=x1 ; x1=x2 ; x2=R*x1+C*x3 ! its housekeeping,
          f1=f2 ; f2=func(x2)          ! and a new function evaluation.
       else                            ! The other outcome,
          x3=x2 ; x2=x1 ; x1=R*x2+C*x0 ! its housekeeping,
          f2=f1 ; f1=func(x1)          ! and its new function evaluation.
       end if
    end do ! Back to see if are done.
    if (f1<f2) then ! We are done. Output the best of the two current values.
       golden=f1
       xmin=x1
    else
       golden=f2
       xmin=x2
    end if
  END FUNCTION golden

  FUNCTION brent(ax,bx,cx,func,tol,xmin)
    ! Brent's method
    ! Reference: W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling, Numerical Recipes in Fortran 90: the art of parallel scientific computing, Cambridge University Press, 1986-1996.
    ! Given a function func, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is between ax and cx, and func(bx) is less than both func(ax) and func(cx)), this routine isolates the minimum to a fractional precision of about tol using Brent’s method. The abscissa of the minimum is returned as xmin, and the minimum function value is returned as brent, the returned function value.
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN) :: ax,bx,cx,tol
    DOUBLE PRECISION,INTENT(OUT) :: xmin
    DOUBLE PRECISION :: brent
    INTERFACE
       DOUBLE PRECISION FUNCTION func(x)
         IMPLICIT NONE
         DOUBLE PRECISION,INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE

    INTEGER,PARAMETER :: ITMAX=100 ! The maximum number of iterations.
    DOUBLE PRECISION,PARAMETER :: CGOLD=0.3819660D0,ZEPS=1.D-10 ! CGOLD is the golden ratio, ZEPS is a small number that protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero.
    INTEGER :: iter
    DOUBLE PRECISION :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

    a=min(ax,cx) ; b=max(ax,cx) ! a and b must be in ascending order, though the input abcissas need not be.
    v=bx ! Initializations...
    w=v
    x=v
    e=0.D0 ! This will be the distance moved on the step before last.
    fx=func(x)
    fv=fx
    fw=fx
    do iter=1,ITMAX
       xm=0.5D0*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.D0*tol1
       if (abs(x-xm)<=(tol2-0.5D0*(b-a))) then ! Test for done here.
          xmin=x ! Arrive here ready to exit with the best values.
          brent=fx
          RETURN
       end if
       if (abs(e)>tol1) then ! Construct a trial parabolic fit.
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.D0*(q-r)
          if (q>0.D0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if (abs(p)>=abs(0.5D0*q*etemp).or.p<=q*(a-x).or.p>=q*(b-x)) then ! These conditions determine the acceptability of the parabolic fit.
             e=merge(a-x,b-x,x>=xm) ! Here it is not o.k., so we take the golden section step into the larger of the two segments.
             d=CGOLD*e
          else
             d=p/q ! Take the parabolic step.
             u=x+d
             if (u-a<tol2.or.b-u<tol2) d=sign(tol1,xm-x)
          end if
       else ! Take the golden section step into the larger of the two segments.
          e=merge(a-x,b-x,x>=xm)
          d=CGOLD*e
       end if
       u=merge(x+d,x+sign(tol1,d),abs(d)>=tol1) ! Arrive here with d computed either from parabolic fit, or else from golden section.
       fu=func(u) ! This is the one function evaluation per iteration,
       if (fu<=fx) then ! Now we have to decide what to do with our function evaluation. Housekeeping follows:
          if (u>=x) then
             a=x
          else
             b=x
          end if
          v=w ; w=x ; x=u
          fv=fw ; fw=fx ; fx=fu
       else
          if (u<x) then
             a=u
          else
             b=u
          end if
          if (fu<=fw.or.w==x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if (fu<=fv.or.v==x.or.v==w) then
             v=u
             fv=fu
          end if
       end if
    end do ! Done with housekeeping. Back for another iteration.
    WRITE(*,*)'Subroutine BRENT: exceeded the maximum number of iterations.'
  END FUNCTION brent
END MODULE optimization_tools
