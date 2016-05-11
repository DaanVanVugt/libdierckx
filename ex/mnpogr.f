cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                  mnpogr : pogrid test program                      cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ..local scalars..
      real cv,del,ermax,er0,exz0,fp,pi,r,sum,sv,x,y,z0,one,ai,s
      integer i,ier,is,j,k,kwrk,lwrk,m,mu,mv,nc,nuest,nu,nvest,nv
c  ..local arrays..
      integer ider(2),iopt(3),iwrk(100),iw(29)
      real u(9),v(20),z(180),c(300),tu(50),tv(50),f(180),wk(116),
     * exact(180),err(9),sp(9),wrk(1600)
c  ..function references..
      real abs,atan2,cos,sin,tespog
c  ..
c  set constants
      one = 1
      pi = atan2(0.,-one)
c we set up the radius of the disc
      r = one
c we set up the number of u (radius)-values of the grid.
      mu = 9
c we set up the u-coordinates of the grid.
      do 10 i=1,mu
         ai = i
         u(i) = ai*0.1
  10  continue
c we set up the number of v (angle)-values of the grid
      mv = 20
c we set up the v-coordinates of the grid.
      del = pi*0.1
      do 20 j=1,mv
         ai = j-1
         v(j) = ai*del-pi
  20  continue
c we fetch the data values at the grid points.
      m = mu*mv
      read(5,900) (z(i),i=1,m)
c we fetch the data value at the origin.
      read(5,900) z0
c we print the data values at the grid points. we also compute and print
c the exact value of the test function underlying the data.
      write(6,905)
      write(6,910) (i,i=1,mu)
      write(6,915)
      exz0 = tespog(0.,0.)
      er0 = abs(exz0-z0)
      ermax = er0
      sum = er0
      do 40 j=1,mv
         cv = cos(v(j))
         sv = sin(v(j))
         k = j
         do 30 i=1,mu
            x = u(i)*cv
            y = u(i)*sv
            exact(k) = tespog(x,y)
            err(i) = abs(exact(k)-z(k))
            sum = sum+err(i)
            if(err(i).gt.ermax) ermax = err(i)
            k = k+mv
  30     continue
         write(6,920) j,(z(k),k=j,m,mv)
         write(6,925) (exact(k),k=j,m,mv)
  40  continue
      ai = m+1
      sum = sum/ai
      write(6,930) z0,exz0
      write(6,935) sum,ermax
c  we set up the dimension information
      nuest = 16
      nvest = 27
      kwrk = 100
      lwrk = 1600
c main loop for the different spline approximations
      do 300 is=1,6
        go to (110,120,130,140,150,160),is
c  we start computing a set of spline approximations with
c  only c0-continuity at the origin,
 110    iopt(2) = 0
        ider(2) = 0
c  non-vanishing at the boundary of the disc,
        iopt(3) = 0
c  with a data value at the origin.
        ider(1) = 0
c  initialisation
        iopt(1) = 0
c  a large value for s for computing the least-squares polynomial
        s = 5.
        go to 200
c  iopt(1) = 1 from the second call on
 120    s = 0.1
        iopt(1) = 1
        go to 200
c  an interpolating spline
 130    s = 0.
        go to 200
c  a second set of approximations with c1-continuity at the origin
 140    iopt(2) = 1
c  vanishing at the boundary of the disc.
        iopt(3) = 1
c  exact value at the origin.
        ider(1) = 1
        z0 = exz0
c reinitialization
        iopt(1) = 0
        s = 0.1
        go to 200
c  no data value at the origin
 150    ider(1) = -1
c  vanishing partial derivatives at the origin
        ider(2) = 1
c reinitialization
        iopt(1) = 0
        go to 200
c finally we calculate the least-squares spline according to the current
c  set of knots
 160    iopt(1) = -1
 200    call pogrid(iopt,ider,mu,u,mv,v,z,z0,r,s,nuest,nvest,
     *    nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c printing of the fitting results.
        if(iopt(1).ge.0) go to 210
        write(6,940)
        go to 220
 210    write(6,945) s
 220    write(6,950) iopt(2)
        if(ider(2).eq.1) write(6,955)
        if(iopt(3).eq.1) write(6,960)
        write(6,965) fp,ier
        write(6,970) nu
        write(6,975)
        write(6,980) (tu(i),i=1,nu)
        write(6,985) nv
        write(6,975)
        write(6,980) (tv(i),i=1,nv)
        nc = (nu-4)*(nv-4)
        write(6,990)
        write(6,980) (c(i),i=1,nc)
c  evaluation of the spline approximation
        call bispev(tu,nu,tv,nv,c,3,3,u,mu,v,mv,f,wk,116,iw,29,ier)
        write(6,995)
        write(6,910) (i,i=1,mu,2)
        write(6,915)
        er0 = abs(exz0-c(1))
        ermax = er0
        sum = er0
        do 240 j=1,mv
          k = j
          do 230 i=1,mu
            sp(i) = f(k)
            err(i) = abs(exact(k)-f(k))
            sum = sum+err(i)
            if(err(i).gt.ermax) ermax = err(i)
            k = k+mv
 230      continue
          if( (j/3)*3 .ne.j ) go to 240
          write(6,920) j,(sp(i),i=1,mu,2)
          write(6,925) (err(i),i=1,mu,2)
 240    continue
        sum = sum/ai
        write(6,1000) c(1),er0
        write(6,935) sum,ermax
 300  continue
      stop
 900  format(10f8.3)
 905  format(49h1data value (exact function value) at grid points)
 910  format(8h u(i),i=,3x,9(i1,7x))
 915  format(8h v(j),j=)
 920  format(1x,i5,9(2x,f6.3))
 925  format(7x,9(2h (,f5.3,1h)))
 930  format(23h0data value at (0,0) = ,f7.3,5x,14hexact value = ,f7.3)
 935  format(19h0mean abs. error = ,f9.3,5x,18hmax. abs. error = ,f9.3)
 940  format(21h0least-squares spline)
 945  format(25h0smoothing spline with s=,f7.2)
 950  format(1x,35horder of continuity at the origin =,i3)
 955  format(1x,43hvanishing partial derivatives at the origin)
 960  format(1x,37hvanishing at the boundary of the disc)
 965  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 970  format(1x,42htotal number of knots in the u-direction =,i3)
 975  format(1x,22hposition of the knots )
 980  format(5x,8f9.4)
 985  format(1x,42htotal number of knots in the v-direction =,i3)
 990  format(23h0b-spline coefficients )
 995  format(50h0spline value (approximation error) at grid points)
1000  format(25h0spline value at (0,0) = ,f7.3,5x,8herror = ,f7.3)
      end
c
      real function tespog(x,y)
c function program tespog calculates the value of the test function
c underlying the data.
c  ..
c  ..scalar arguments..
      real x,y,f
c  ..
      f = 1.-((3.*x-1.)**2+(3.*y-1.)**2)/(11.-6.*(x+y))
      tespog = f-(1.-x**2-y**2)*(x+y)*54./121.
      return
      end
