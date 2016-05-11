cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc             mnsphe : sphere test program                           cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real teta(192),phi(192),r(192),w(192),tp(30),tt(30),c(300),
     * p(9),t(9),f(81),wrk1(12000),wrk2(72)
      integer iwrk(300)
      real eps,fp,pi,pi2,pi4,s,scale,scp,sct,ai
      integer i,ier,iopt,j,kwrk,lwrk1,lwrk2,l1,l2,l,m,np,npest,nt,ntest
     * ,is,i1,i2,nc,ntt,npp
      real atan,testsp
c  set constants
      pi4 = atan(0.1e+01)
      pi = pi4*4
      pi2 = pi+pi
      scale = pi4/0.45e+02
c  we fetch the number of data points.
      m = 192
c  we fetch and print the latitude - longitude coordinates of the data
c  points (in degrees).
      write(6,900)
      write(6,905)
      l2 = 0
      do 10 i=1,48
         l1 = l2+1
         l2 = l2+4
         read(5,910) (teta(l),phi(l),l=l1,l2)
         write(6,915)(teta(l),phi(l),l=l1,l2)
  10  continue
c  we set up the weights, scale into radians the latitude-longitude
c  coordinates and calculate the function values.
      do 20 i=1,m
         w(i) = 0.1e+01
         teta(i) = teta(i)*scale
         phi(i) = phi(i)*scale
         if(teta(i).gt.pi) teta(i) = pi
         if(phi(i).gt.pi2) phi(i) = pi2
         r(i) = testsp(teta(i),phi(i))
  20  continue
c  we set up the coordinates of the grid points for the evaluation of
c  the spline approximations.
      sct = pi/8
      scp = pi2/8
      do 30 i=1,8
         ai = i-1
         t(i) = ai*sct
         p(i) = ai*scp
  30  continue
      t(9) = pi
      p(9) = pi2
c we set up the dimension information
      ntest = 15
      npest = 19
      lwrk1 = 12000
      lwrk2 = 72
      kwrk = 300
c  we choose a value for eps
      eps = 0.1e-05
c  main loop for the different spline approximations
      do 300 is=1,4
        go to (110,120,130,140),is
c  we start computing the least-squares constrained polynomial (large s)
 110    iopt = 0
        s = 500.
        go to 200
c  iopt = 1 from the second call on.
 120    iopt = 1
        s = 135.
        go to 200
 130    s = 15.
        go to 200
c  a least-squares spherical spline with specified knots.
 140    iopt = -1
c  we set up the number of knots.
        nt = 11
        np = 15
c  we set up the position of the interior knots of the spline.
        ntt = nt-8
        do 150 i=1,ntt
          ai = i
          j = i+4
          tt(j) = ai*pi4
 150    continue
        npp = np-8
        do 160 i=1,npp
          ai = i
          j = i+4
          tp(j) = ai*pi4
 160    continue
c  determination of the spline approximation.
 200    call sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps,nt,tt,
     *   np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c  printing of the fitting results.
        if(iopt.ge.0) go to 210
        write(6,920)
        go to 220
 210    write(6,925)
        write(6,930) s
 220    write(6,935) fp,ier
        write(6,940) nt
        write(6,945)
        write(6,950) (tt(i),i=1,nt)
        write(6,955) np
        write(6,945)
        write(6,950) (tp(i),i=1,np)
        nc = (nt-4)*(np-4)
        write(6,960)
        write(6,965) (c(i),i=1,nc)
c  evaluation of the spline approximation.
        call bispev(tt,nt,tp,np,c,3,3,t,9,p,9,f,wrk2,lwrk2,
     *   iwrk,kwrk,ier)
        write(6,970) (p(i),i=1,9)
        write(6,975)
        i2 = 0
        do 230 i=1,9
          i1 = i2+1
          i2 = i2+9
          write(6,980) t(i),(f(j),j=i1,i2)
 230    continue
 300  continue
      stop
c  format statements.
 900  format(55h1latitude-longitude values of the data points (degrees))
 905  format(1h0,4(3x,10hteta   phi,3x))
 910  format(8f6.0)
 915  format(1h ,4(3x,f4.0,2x,f4.0,3x))
 920  format(50h0least-squares spline approximation on the sphere.)
 925  format(32h0smoothing spline on the sphere.)
 930  format(20h smoothing factor s=,f9.0)
 935  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 940  format(1x,45htotal number of knots in the teta-direction =,i3)
 945  format(1x,22hposition of the knots )
 950  format(5x,8f8.4)
 955  format(1x,44htotal number of knots in the phi-direction =,i3)
 960  format(23h0b-spline coefficients )
 965  format(5x,8f9.4)
 970  format(9h      phi,9f7.3)
 975  format(6h  teta)
 980  format(1h ,f6.3,2x,9f7.3)
      end
      real function testsp(v,u)
c function program testsp calculates the value of a test function for
c the sphere package.
c ..
      real cos,cu,cv,rad1,rad2,rad3,sin,sqrt,su,sv,u,v
      cu = cos(u)
      cv = cos(v)
      su = sin(u)
      sv = sin(v)
      rad1 = (cu*sv*0.2)**2+(su*sv)**2+(cv*0.5)**2
      rad2 = (cu*sv)**2+(su*sv*0.5)**2+(cv*0.2)**2
      rad3 = (cu*sv*0.5)**2+(su*sv*0.2)**2+cv**2
      testsp = 1./sqrt(rad1) + 1./sqrt(rad2) + 1./sqrt(rad3)
      return
      end
