cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc              mnsuev : surev test program                           cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real fac
      integer i,idim,ier,j,m,mu,mv,m0,m1,m2,m3,nc,nu4,nv4,nu,nv,l
      real tu(11),tv(10),c(126),u(6),v(6),f(108),wrk(48)
      integer iwrk(12)
c  we set up the grid points for evaluating the spline surface.
      mu = 6
      mv = 6
      do 10 i=1,6
      u(i) = (i-1)*0.2
      v(i) = u(i)
  10  continue
c  the interior knots with respect to the u-variable.
      tu(5) = 0.4
      tu(6) = 0.7
      tu(7) = 0.9
      nu = 11
c  the interior knots with respect to the v-variable.
      tv(5) = 0.3
      tv(6) = 0.8
      nv = 10
c  the boundary knots
      do 20 i=1,4
        tu(i) = 0.
        tv(i) = 0.
        tu(i+7) = 1.
        tv(i+6) = 1.
  20  continue
c  we generate the b-spline coefficients for the test surface
c        x = u*v    y = v**2    z = u+v     0 <= u,v <= 1
c  the dimension of the surface
      idim = 3
c  the number of b-spline coefficients for each co-ordinate
      nu4 = nu-4
      nv4 = nv-4
      nc = nu4*nv4
c  the coefficients for x = u*v
      do 30 i=1,nv4
        c(i) = 0.
  30  continue
      do 40 i=2,nu4
        c((i-1)*nv4+1) = 0.
  40  continue
      m0 = 1
      do 60 i=2,nu4
        m1 = m0+nv4
        fac = (tu(i+3)-tu(i))/9.
        do 50 j=2,nv4
          m2 = m0+1
          m3 = m1+1
          c(m3) = c(m1)+c(m2)-c(m0)+fac*(tv(j+3)-tv(j))
          m0 = m0+1
          m1 = m1+1
  50    continue
        m0 = m0+1
  60  continue
c  the coefficients for y = v**2.
      l = nc
      m0 = l+1
      m1 = m0+1
      c(m0) = 0.
      c(m1) = 0.
      do 70 i=3,nv4
        c(m1+1) = c(m1)+(tv(i+3)-tv(i))*((c(m1)-c(m0))/(tv(i+2)-tv(i-1))
     *    +(tv(i+2)-tv(i))/3.)
        m0 = m1
        m1 = m0+1
  70  continue
      do 80 i=1,nv4
        m0 = l+i
        fac = c(m0)
        do 80 j=1,nu4
          m0 = m0+nv4
          c(m0) = fac
  80  continue
c  the coefficients for z = u+v
      l = l+nc
      m0 = l+1
      c(m0) = 0.
      do 90 i=2,nv4
        m1 = m0+1
        c(m1) = c(m0)+(tv(i+3)-tv(i))/3.
        m0 = m1
  90  continue
      do 100 i=1,nv4
        m0 = l+i
        do 100 j=2,nu4
          m1 = m0+nv4
          c(m1) = c(m0)+(tu(j+3)-tu(j))/3.
          m0 = m1
 100  continue
c  evaluation of the spline surface
      call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,108,wrk,48,iwrk,12,ier)
c  printing of the results
      write(6,900)
      write(6,910)
      write(6,920) (tu(i),i=1,nu)
      write(6,930)
      write(6,920) (tv(i),i=1,nv)
      write(6,940)
      m1 = 0
      do 110 l=1,idim
        m0 = m1+1
        m1 = m1+nc
        write(6,950) (c(j),j=m0,m1)
 110  continue
      write(6,960)
      write(6,970) (v(i),i=1,mv)
      write(6,980)
      m = mu*mv
      m0 = 0
      do 130 i=1,mu
        write(6,990) u(i)
        m1 = m0
        do 120 l=1,idim
          m2 = m1+1
          m3 = m1+mv
          write(6,995) (f(j),j=m2,m3)
          m1 = m1+m
 120    continue
        m0 = m0+mv
 130  continue
      stop
c  format statements.
 900  format(23h0bicubic spline surface)
 910  format(1x,40hposition of the knots in the u-direction)
 920  format(1x,15f5.1)
 930  format(1x,40hposition of the knots in the v-direction)
 940  format(23h b-spline coefficients )
 950  format(5x,8f9.4)
 960  format(1h0,37hspline values at selected grid points)
 970  format(1h0,2x,1hv,6(3x,f4.1))
 980  format(1h ,1x,1hu)
 990  format(1h ,f4.1)
 995  format(5x,6f7.3)
      end
