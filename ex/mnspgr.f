cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                  mnspgr : spgrid test program                      cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ..local scalars..
      real del,ermax,erf,exr0,exr1,fp,pi,sum,r0,r1,one,ai,s
      integer i,ier,is,j,k,kwrk,lwrk,m,mu,mv,nc,nuest,nu,nvest,nv
c  ..local arrays..
      integer ider(4),iopt(3),iwrk(70),iw(25)
      real u(11),v(14),r(154),c(300),tu(25),tv(25),f(154),wk(100),
     * exact(154),err(14),sp(14),wrk(1500)
c  ..function references..
      real abs,atan2,tesspg
c  ..
c  set constants
      one = 1
      pi = atan2(0.,-one)
      del = pi*0.05
c we set up the number of u (latitude)-values of the grid.
      mu = 11
c we set up the u-coordinates of the grid.
      read(5,905)(iw(i),i=1,mu)
      do 10 i=1,mu
         ai = iw(i)
         u(i) = ai*del
  10  continue
c we set up the number of v (longitude)-values of the grid
      mv = 14
c we set up the v-coordinates of the grid.
      read(5,905)(iw(i),i=1,mv)
      do 20 i=1,mv
         ai = iw(i)
         v(i) = ai*del
  20  continue
c we fetch the data values at the grid points.
      m = mu*mv
      read(5,900) (r(i),i=1,m)
c we print the data values at the grid points. we also compute and print
c the exact value of the test function underlying the data.
      write(6,910)
      write(6,915) (j,j=1,mv,2)
      write(6,920)
      exr0 = tesspg(0.,0.)
      exr1 = tesspg(pi,0.)
      ermax = 0.
      sum = 0.
      l = 0
      do 40 i=1,mu
        l = (i-1)*mv+1
        k = 1
        do 30 j=1,7
          exact(l) = tesspg(u(i),v(k))
          erf = abs(exact(l)-r(l))
          sum = sum+erf
          if(erf.gt.ermax) ermax = erf
          sp(j) = exact(l)
          err(j) = r(l)
          l = l+2
          k = k+2
  30    continue
        write(6,925) i,(err(j),j=1,7)
        write(6,930) (sp(j),j=1,7)
  40  continue
      write(6,915) (j,j=2,mv,2)
      write(6,920)
      do 60 i=1,mu
        l = (i-1)*mv+2
        k = 2
        do 50 j=1,7
          exact(l) = tesspg(u(i),v(k))
          erf = abs(exact(l)-r(l))
          sum = sum+erf
          if(erf.gt.ermax) ermax = erf
          sp(j) = exact(l)
          err(j) = r(l)
          l = l+2
          k = k+2
  50    continue
        write(6,925) i,(err(j),j=1,7)
        write(6,930) (sp(j),j=1,7)
  60  continue
      ai = m
      sum = sum/ai
      write(6,935) sum,ermax
      write(6,940) exr0,exr1
c  we set up the dimension information
      nuest = 19
      nvest = 21
      kwrk = 70
      lwrk = 1500
c main loop for the different spline approximations
      do 300 is=1,6
        go to (110,120,130,140,150,160),is
c  we start computing a set of spline approximations with
c  only c0-continuity at the poles,
 110    iopt(2) = 0
        ider(2) = 0
        iopt(3) = 0
        ider(4) = 0
c  with no data values at the poles.
        ider(1) = -1
        ider(3) = -1
c  initialisation
        iopt(1) = 0
c  a large value for s for computing the least-squares polynomial
        s = 60.
        go to 200
c  iopt(1) = 1 from the second call on
 120    s = 0.05
        iopt(1) = 1
        go to 200
c  an interpolating spline
 130    s = 0.
        go to 200
c  a second set of approximations with c1-continuity at the poles
 140    iopt(2) = 1
        iopt(3) = 1
c  exact values at the poles.
        ider(1) = 1
        ider(3) = 1
        r0 = exr0
        r1 = exr1
c reinitialization
        iopt(1) = 0
        s = 0.05
        go to 200
c  vanishing derivatives at the poles
 150    ider(2) = 1
        ider(4) = 1
c reinitialization
        iopt(1) = 0
        go to 200
c finally we calculate the least-squares spline according to the current
c  set of knots
 160    iopt(1) = -1
 200    call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,
     *    nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c printing of the fitting results.
        if(iopt(1).ge.0) go to 210
        write(6,945)
        go to 220
 210    write(6,950) s
 220    write(6,955) iopt(2),iopt(3)
        if(ider(2).eq.1) write(6,960)
        if(ider(4).eq.1) write(6,965)
        write(6,970) fp,ier
        write(6,975) nu
        write(6,980)
        write(6,985) (tu(i),i=1,nu)
        write(6,990) nv
        write(6,980)
        write(6,985) (tv(i),i=1,nv)
        nc = (nu-4)*(nv-4)
        write(6,995)
        write(6,985) (c(i),i=1,nc)
c  evaluation of the spline approximation
        call bispev(tu,nu,tv,nv,c,3,3,u,mu,v,mv,f,wk,100,iw,25,ier)
        write(6,1000)
        write(6,915) (j,j=1,mv,2)
        write(6,920)
        ermax = 0.
        sum = 0.
        do 240 i=1,mu
          k = i
          do 230 j=1,mv
            sp(j) = f(k)
            err(j) = abs(exact(k)-f(k))
            sum = sum+err(j)
            if(err(j).gt.ermax) ermax = err(j)
            k = k+1
 230      continue
          if( (i/2)*2 .ne.i ) go to 240
          write(6,925) i,(sp(j),j=1,mv,2)
          write(6,930) (err(j),j=1,mv,2)
 240    continue
        sum = sum/ai
        write(6,935) sum,ermax
        write(6,1005) c(1),c(nc)
 300  continue
      stop
 900  format(7f8.3)
 905  format(14i3)
 910  format(49h1data value (exact function value) at grid points)
 915  format(8h0v(j),j=,3x,7(i2,6x))
 920  format(8h u(i),i=)
 925  format(1x,i5,7(2x,f6.3))
 930  format(7x,7(2h (,f5.3,1h)))
 935  format(19h0mean abs. error = ,f9.3,5x,18hmax. abs. error = ,f9.3)
 940  format(30h function values at the poles ,f7.3,5x,f7.3)
 945  format(21h0least-squares spline)
 950  format(25h0smoothing spline with s=,f7.2)
 955  format(1x,34horder of continuity at the poles =,2i5)
 960  format(1x,37hvanishing derivatives at the pole u=0)
 965  format(1x,38hvanishing derivatives at the pole u=pi)
 970  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 975  format(1x,42htotal number of knots in the u-direction =,i3)
 980  format(1x,22hposition of the knots )
 985  format(5x,8f9.4)
 990  format(1x,42htotal number of knots in the v-direction =,i3)
 995  format(23h0b-spline coefficients )
1000  format(50h0spline value (approximation error) at grid points)
1005  format(28h spline values at the poles ,f7.3,5x,f7.3)
      end
c
      real function tesspg(u,v)
c function program tesspg calculates the value of the test function
c underlying the data.
      real u,v,sin,cos
      tesspg = 2./(4.1+cos(3.*u)+3.*cos(v+v+u*0.25)*sin(u)**2)
      return
      end

