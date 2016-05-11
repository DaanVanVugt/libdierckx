cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mnfour : fourco test program                       cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real c(20),t(20),wrk1(20),wrk2(20),alfa(10),ress(10),resc(10)
      integer i,ier,j,k,k1,m,n,nk1
      real ak,rc,rs
c  as an example we calculate some integrals of the form
c          / 1                               / 1
c         !    x * sin(alfa*x) dx  and      !   x * cos(alfa*x) dx
c      0 /                               0 /
c
c  we will represent y = x as a cubic spline.
      k = 3
      k1 = k+1
c  we fetch the knots of the cubic spline
      n = 2*k1+4
c  the boundary knots
      j = n
      do 10 i=1,k1
         t(i) = 0.
         t(j) = 0.1e+01
         j = j-1
  10  continue
c  the interior knots
      t(5) = 0.1e+0
      t(6) = 0.3e+0
      t(7) = 0.4e+0
      t(8) = 0.8e+0
c  find the b-spline representation of y=x
      nk1 = n-k1
      ak = k
      c(1) = 0.
      do 20 i=2,nk1
         j = i+k
         c(i) = c(i-1)+(t(j)-t(i))/ak
  20  continue
c  print the data for the spline.
      write(6,900) k
      write(6,905)
      write(6,910) (t(i),i=1,n)
      write(6,915)
      write(6,920) (c(i),i=1,nk1)
c  fetch the different values for alfa
      m = 8
      alfa(1) = 0.
      alfa(2) = 0.1e-02
      do 30 i=3,m
        alfa(i) = -alfa(i-1)*0.1e+02
  30  continue
c  calculate the fourier integrals of the cubic spline
      call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
c  print the results
      write(6,925)
      do 40 i=1,m
c  fetch the exact values of the integrals
        call exfour(alfa(i),rs,rc)
        write(6,930) alfa(i),ress(i),rs,resc(i),rc
  40  continue
      stop
c  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,12f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(1h0,2x,4halfa,9x,4hress,9x,4hexas,9x,4hresc,9x,4hexac)
 930  format(1x,e8.1,4f13.5)
      end
      subroutine exfour(alfa,rs,rc)
c  subroutine exfour calculates the integrals
c                 / 1
c      rs =      !    x*sin(alfa*x) dx    and
c             0 /
c                 / 1
c      rc =      !    x*cos(alfa*x) dx
c             0 /
      integer k,k2
      real aa,ak,alfa,cc,c1,half,one,rc,rs,ss,s1,three
c  ..function references..
      real cos,sin,abs
c  set constants
      one = 0.1e+01
      three = 0.3e+01
      half = 0.5e0
      if(abs(alfa).lt.one) go to 10
c  integration by parts
      aa = one/alfa
      cc = cos(alfa)
      ss = sin(alfa)
      rs = (ss*aa-cc)*aa
      rc = ((cc-one)*aa+ss)*aa
      go to 50
c  using the series expansions of sin(alfa*x) and cos(alfa*x)
  10  rs = 0.
      rc = half
      if(alfa.eq.0.) go to 50
      rs = alfa/three
      ss = rs
      cc = rc
      aa = -alfa*alfa
      do 20 k=1,21
        k2 = 2*(k-1)
        ak = (k2+2)*(k2+5)
        ss = ss*aa/ak
        s1 = rs+ss
        if(s1.eq.rs)go to 30
        rs = s1
  20  continue
  30  do 40 k=1,21
        k2 = 2*(k-1)
        ak = (k2+1)*(k2+4)
        cc = cc*aa/ak
        c1 = rc+cc
        if(c1.eq.rc)go to 50
        rc = c1
  40  continue
  50  return
      end
