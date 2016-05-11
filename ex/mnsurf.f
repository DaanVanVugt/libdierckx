cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc        mnsurf : surfit test program                                cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(80),y(80),z(80),w(80),tx(15),ty(15),c(200),wrk1(12000),
     * wrk2(6000),xx(11),yy(11),zz(121)
      integer iwrk(300)
      integer i,ier,iopt,is,j,kwrk,kx,ky,lwrk1,lwrk2,m,mx,my,nc,
     * nmax,nx,nxest,ny,nyest
      real ai,delta,eps,fp,s,ww,xb,xe,yb,ye
c  we fetch the number of data points
      read(5,900) m
      write(6,905) m
c  we fetch the co-ordinate and function values of each data point.
      write(6,910)
      do 10 i=1,m
        read(5,915) x(i),y(i),z(i)
        if((i/2)*2.ne.i) go to 10
        j = i-1
        write(6,920) j,x(j),y(j),z(j),i,x(i),y(i),z(i)
  10  continue
c  we fetch an estimate of the standard deviation of the data values.
      read(5,925) delta
      write(6,930) delta
c  the weights are set equal to delta**(-1)
      ww = 1./delta
      do 20 i=1,m
        w(i) = ww
  20  continue
c  we set up the boundaries of the approximation domain.
      xb = -2.
      xe = 2.
      yb = -2.
      ye = 2.
c we generate a rectangular grid for evaluating the splines.
      mx = 11
      my = 11
      do 30 i=1,11
        ai = i-6
        xx(i) = ai*0.4
        yy(i) = xx(i)
  30  continue
c  we set up the dimension information
      nxest = 15
      nyest = 15
      nmax = 15
      kwrk = 300
      lwrk1 = 12000
      lwrk2 = 6000
c  we choose a value for eps
      eps=0.1e-05
c  main loop for the different spline approximations.
      do 300 is=1,6
        go to (110,120,130,140,150,160),is
c  we start computing the least-squares bicubic polynomial (large s)
 110    iopt = 0
        kx = 3
        ky = 3
        s = 900000.
        go to 200
c  iopt=1 from the second call on.
 120    iopt = 1
        s = 200.
        go to 200
c  a value for s within its confidence interval
 130    s = m
        go to 200
c  overfitting (s too small)
 140    s = 20.
        go to 200
c  we change the degrees of the spline
 150    iopt = 0
        kx = 5
        ky = 5
        s = m
        go to 200
c  finally, we also calculate a least-squares spline approximation
c  with specified knots.
 160    iopt = -1
        kx = 3
        ky = 3
        nx = 11
        ny = 11
        j = kx+2
        do 170 i=1,3
          ai = i-2
          tx(j) = ai
          ty(j) = ai
          j = j+1
 170    continue
c  determination of the spline approximation.
 200    call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     *   nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c  printing of the fitting results.
        if(iopt.ge.0) go to 210
        write(6,935) kx,ky
        go to 220
 210    write(6,940) kx,ky
        write(6,945) s
 220    write(6,950) fp,ier
        write(6,955) nx
        write(6,960)
        write(6,965) (tx(i),i=1,nx)
        write(6,970) ny
        write(6,960)
        write(6,965) (ty(i),i=1,ny)
        nc = (nx-kx-1)*(ny-ky-1)
        write(6,975)
        write(6,980) (c(i),i=1,nc)
c  evaluation of the spline approximation.
        call bispev(tx,nx,ty,ny,c,kx,ky,xx,mx,yy,my,zz,
     *   wrk2,lwrk2,iwrk,kwrk,ier)
        write(6,1000)
        write(6,985) (xx(i),i=1,mx)
        write(6,990)
        do 230 j=1,my
          write(6,995) yy(j),(zz(i),i=j,121,11)
 230    continue
 300  continue
      stop
c  format statements.
 900  format(i3)
 905  format(1h1,i3,12h data points)
 910  format(1h0,2(2x,1hi,5x,4hx(i),6x,4hy(i),6x,4hz(i),6x))
 915  format(3f10.4)
 920  format(1x,2(i3,3f10.4,5x))
 925  format(e20.6)
 930  format(1x,40hestimate of standard deviation of z(i) =,e15.6)
 935  format(32h0least-squares spline of degrees,2i3)
 940  format(28h0smoothing spline of degrees,2i3)
 945  format(20h smoothing factor s=,f9.0)
 950  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 955  format(1x,42htotal number of knots in the x-direction =,i3)
 960  format(1x,22hposition of the knots )
 965  format(5x,10f7.3)
 970  format(1x,42htotal number of knots in the y-direction =,i3)
 975  format(23h0b-spline coefficients )
 980  format(5x,8f9.4)
 985  format(1h0,1hx,2x,11f7.1)
 990  format(3x,1hy)
 995  format(1x,f4.1,11f7.3)
 1000 format(1h0,33hspline evaluation on a given grid)
      end

