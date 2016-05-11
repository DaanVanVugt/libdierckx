cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc               mnregr : regrid test program                         cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(11),y(11),z(121),tx(17),ty(17),c(300),wrk(850),f(121),
     * wk(132)
      integer iwrk(60),iw(22)
      real ai,fp,s,xb,xe,yb,ye
      integer kx,ky,kwrk,lwrk,m,mx,my,m1,m2,nc,nx,nxest,ny,nyest,
     * i,ier,is,iopt,j
c  we fetch the number of x-coordinates of the grid.
      read(5,900) mx
c  we fetch the x-coordinates of the grid.
      read(5,905) (x(i),i=1,mx)
c  we fetch the number of y-coordinates of the grid.
      read(5,900) my
c  we fetch the y-coordinates of the grid.
      read(5,905) (y(i),i=1,my)
c  we fetch the function values at the grid points.
      m = mx*my
      read(5,910) (z(i),i=1,m)
c  printing of the input data.
      write(6,915)
      write(6,920) (y(i),i=1,6)
      write(6,925)
      m1 = 1
      do 10 i=1,mx
        m2 = m1+5
        write(6,930) x(i),(z(j),j=m1,m2)
        m1 = m1+my
  10  continue
      write(6,920) (y(i),i=7,my)
      write(6,925)
      m1 = 7
      do 20 i=1,mx
        m2 = m1+4
        write(6,930) x(i),(z(j),j=m1,m2)
        m1 = m1+my
  20  continue
c  we set up the boundaries of the approximation domain.
      xb = x(1)
      yb = y(1)
      xe = x(mx)
      ye = y(my)
c  we set up the dimension information
      nxest = 17
      nyest = 17
      lwrk = 850
      kwrk = 60
c  main loop for the different spline approximations
      do 300 is=1,6
        go to (110,120,130,140,150,160),is
c  we start computing the least-squares bicubic polynomial
 110    iopt = 0
        kx = 3
        ky = 3
        s = 10.
        go to 200
c  iopt=1 from the second call on
 120    iopt = 1
        s = 0.22
        go to 200
c  overfitting (s too small)
 130    s = 0.1
        go to 200
c  an interpolating spline
 140    s = 0.
        go to 200
c  we change the degrees of the spline
 150    kx = 5
        ky = 5
        s = 0.2
        iopt = 0
        go to 200
c  finally we also calculate a least-squares spline approximation
c  with specified knots.
 160    iopt = -1
        kx = 3
        ky = 3
        nx = 11
        ny = 11
        j = kx+2
        do 170 i=1,3
          ai = i-2
          tx(j) = ai*0.5
          ty(j) = tx(j)
          j = j+1
 170    continue
c  determination of the spline approximation.
 200    call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,
     *   nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
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
        call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,f,
     *   wk,132,iw,22,ier)
        write(6,985)
        write(6,920) (y(i),i=1,my,2)
        write(6,925)
        m1 = 1
        do 230 i=1,mx,2
          m2 = m1+my-1
          write(6,930) x(i),(f(j),j=m1,m2,2)
          m1 = m1+2*my
 230    continue
 300  continue
      stop
c  format statements.
 900  format(i2)
 905  format(11f5.1)
 910  format(11f7.4)
 915  format(15h1the input data)
 920  format(1h0,8x,1hy,4x,6(4x,f4.1))
 925  format(1h ,7x,1hx)
 930  format(6x,f4.1,5x,6f8.4)
 935  format(32h0least-squares spline of degrees,2i3)
 940  format(28h0smoothing spline of degrees,2i3)
 945  format(20h smoothing factor s=,f7.2)
 950  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 955  format(1x,42htotal number of knots in the x-direction =,i3)
 960  format(1x,22hposition of the knots )
 965  format(5x,10f6.2)
 970  format(1x,42htotal number of knots in the y-direction =,i3)
 975  format(23h0b-spline coefficients )
 980  format(5x,8f9.4)
 985  format(1h0,37hspline values at selected grid points)
      end

