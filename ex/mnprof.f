cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc              mnprof : profil test program                          cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real fac,facx,u
      integer i,ier,iopt,j,kx,kx1,ky,ky1,m,mx,my,m0,m1,m2,m3,nc,
     * nkx1,nky1,nx,ny
      real tx(15),ty(15),c(100),x(6),y(6),z(36),cc(15)
c  we set up the grid points for evaluating the tensor product splines.
      mx = 6
      my = 6
      m = mx*my
      do 10 i=1,6
      x(i) = (i-1)*0.2
      y(i) = x(i)
  10  continue
c  loop for different spline degrees with respect to the x-variable
      do 300 kx=3,5,2
c  the knots in the x-direction
        tx(kx+2) = 0.4
        tx(kx+3) = 0.7
        tx(kx+4) = 0.9
        kx1 = kx+1
        nx = 3+2*kx1
        j = nx
        do 20 i=1,kx1
          tx(i) = 0.
          tx(j) = 1.
          j = j-1
  20    continue
c  loop for different spline degrees with respect to the y-variable
      do 200 ky=2,3
c  the knots in the y-direction
        ty(ky+2) = 0.3
        ty(ky+3) = 0.8
        ky1 = ky+1
        ny = 2+2*ky1
        j = ny
        do 30 i=1,ky1
          ty(i) = 0.
          ty(j) = 1.
          j = j-1
  30    continue
c  we generate the b-spline coefficients for the test function x*y
        nkx1 = nx-kx1
        nky1 = ny-ky1
        do 40 i=1,nky1
          c(i) = 0.
  40    continue
        do 50 i=2,nkx1
          c((i-1)*nky1+1) = 0.
  50    continue
        fac = kx*ky
        m0 = 1
        do 70 i=2,nkx1
          m1 = m0+nky1
          facx = (tx(i+kx)-tx(i))/fac
          do 60 j=2,nky1
            m2 = m0+1
            m3 = m1+1
            c(m3) = c(m1)+c(m2)-c(m0)+facx*(ty(j+ky)-ty(j))
            m0 = m0+1
            m1 = m1+1
  60      continue
          m0 = m0+1
  70    continue
c  printing of the spline information
        write(6,900) kx,ky
        write(6,910)
        write(6,920) (tx(i),i=1,nx)
        write(6,930)
        write(6,920) (ty(i),i=1,ny)
        nc = nkx1*nky1
        write(6,940)
        write(6,950) (c(i),i=1,nc)
c  we calculate a number of profiles f(y)=s(u,y)
        iopt = 0
        m0 = 1
        do 80 i=1,mx
          u = x(i)
          call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,15,cc,ier)
          write(6,955) u
          write(6,950) (cc(j),j=1,nky1)
c  evaluation of the one-dimensional spline f(y)
          call splev(ty,ny,cc,ky,y,z(m0),my,ier)
          m0 = m0+my
  80    continue
        write(6,960)
        write(6,970) (y(i),i=1,my)
        write(6,980)
        m2 = 0
        do 100 i=1,mx
          m1 = m2+1
          m2 = m2+my
          write(6,990) x(i),(z(j),j=m1,m2)
 100    continue
c  we calculate a number of profiles g(x)=s(x,u)
        iopt = 1
        m0 = 1
        do 120 i=1,my
          u = y(i)
          call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,15,cc,ier)
          write(6,995) u
          write(6,950) (cc(j),j=1,nkx1)
c  evaluation of the one-dimensional spline g(x)
          call splev(tx,nx,cc,kx,x,z(m0),mx,ier)
          m0 = m0+mx
 120    continue
        write(6,960)
        write(6,970) (y(i),i=1,my)
        write(6,980)
        do 140 i=1,mx
          write(6,990) x(i),(z(j),j=i,m,mx)
 140    continue
 200    continue
 300  continue
      stop
c  format statements.
 900  format(33h0tensor product spline of degrees,2i3)
 910  format(1x,40hposition of the knots in the x-direction)
 920  format(1x,15f5.1)
 930  format(1x,40hposition of the knots in the y-direction)
 940  format(23h b-spline coefficients )
 950  format(1x,8f9.4)
 955  format(45h0b-spline coefficients of the profile f(y)=s(,f3.1,
     * 3h,y))
 960  format(1h0,37hspline values at selected grid points)
 970  format(1h0,8x,1hy,4x,6(4x,f4.1))
 980  format(1h ,7x,1hx)
 990  format(6x,f4.1,5x,6f8.2)
 995  format(47h0b-spline coefficients of the profile g(x)=s(x,,f3.1,
     * 1h))
      end
