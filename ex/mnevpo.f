cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc              mnevpo : evapol test program                          cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real fac,r1,r2
      integer i,ir,j,m,mx,my,m0,m1,m2,nc,nu4,nv4,nu,nv
      real tu(11),tv(10),c(42),x(6),y(6),f(36)
      external r1,r2
c  we set up the grid points for evaluating the polar spline.
      mx = 6
      my = 6
      do 10 i=1,6
      x(i) = (2*i-7)*0.1
      y(i) = x(i)
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
c  the number of b-spline coefficients
      nu4 = nu-4
      nv4 = nv-4
      nc = nu4*nv4
c  we generate the b-spline coefficients for the function s(u,v)=u**2
      m0 = 1
      m1 = m0+nv4
      c(m0) = 0.
      c(m1) = 0.
      do 70 i=3,nu4
        m2 = m1+nv4
        c(m2) = c(m1)+(tu(i+3)-tu(i))*((c(m1)-c(m0))/(tu(i+2)-tu(i-1))
     *    +(tu(i+2)-tu(i))/3.)
        m0 = m1
        m1 = m2
  70  continue
      do 80 i=1,nu4
        m0 = (i-1)*nv4+1
        fac = c(m0)
        do 80 j=2,nv4
          m0 = m0+1
          c(m0) = fac
  80  continue
      write(6,900)
      write(6,910)
      write(6,920) (tu(i),i=1,nu)
      write(6,930)
      write(6,920) (tv(i),i=1,nv)
      write(6,940)
      write(6,950) (c(j),j=1,nc)
c the spline s(u,v) defines a function f(x,y) through the transformation
c    x = r(v)*u*cos(v)   y = r(v)*u*sin(v)
c we consider two different functions r(v)
      do 300 ir=1,2
        go to (110,130),ir
c if r(v) =1 and s(u,v) = u**2 then f(x,y) = x**2+y**2
c evaluation of f(x,y)
 110    m = 0
        do 120 i=1,my
          do 120 j=1,mx
            m = m+1
            f(m) = evapol(tu,nu,tv,nv,c,r1,x(j),y(i))
 120    continue
        write(6,960)
        go to 200
c if r(v) = (1+cos(v)**2)/2 and s(u,v) = u**2 then f(x,y) =
c    4*(x**2+y**2)**3/(4*x**4 +y**4 +4*x**2*y**2)
c evaluation of f(x,y)
 130    m = 0
        do 140 i=1,my
          do 140 j=1,mx
            m = m+1
            f(m) = evapol(tu,nu,tv,nv,c,r2,x(j),y(i))
 140    continue
        write(6,965)
 200    write(6,970) (x(i),i=1,mx)
        write(6,975)
        m1 = 0
        do 210 j=1,my
          m0 = m1+1
          m1 = m1+mx
          write(6,980) y(j),(f(m),m=m0,m1)
 210    continue
 300  continue
      stop
c  format statements.
 900  format(24h0polar spline evaluation)
 910  format(1x,40hposition of the knots in the u-direction)
 920  format(1x,15f5.1)
 930  format(1x,40hposition of the knots in the v-direction)
 940  format(23h b-spline coefficients )
 950  format(5x,8f9.4)
 960  format(1h0,30hf(x,y) corresponding to r(v)=1)
 965  format(1h0,44hf(x,y) corresponding to r(v)=(1+cos(v)**2)/2)
 970  format(1h0,1hx,2x,6f7.1)
 975  format(3x,1hy)
 980  format(1x,f4.1,6f7.3)
      end
      real function r1(v)
      real v
      r1 = 1.
      return
      end
      real function r2(v)
      real v,cos
      r2 = (1.+cos(v)**2)*0.5
      return
      end

