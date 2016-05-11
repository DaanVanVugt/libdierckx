cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mnpola : polar test program                        cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(200),y(200),z(200),w(200),u(200),v(200),tu(30),tv(30),
     * c(300),exact(200),f(200),wrk1(15000),wrk2(5700)
      integer iopt(3),iwrk(500)
      integer i,is,ier,kwrk,l,lwrk1,lwrk2,l1,l2,m,m1,m2,nc,nu,nv,
     * nuest,nvest
      real s,fp,eps,sum,ermax,error,ai
      real abs,rad1,rad2,testpo,evapol
      external rad1,rad2
c  we fetch the number of data points.
      m1 = 200
      m2 =  90
c  we fetch and print the coordinates and function values of the data.
      write(6,900)
      write(6,905)
      l2 = 0
      do 10 i=1,50
         l1 = l2+1
         l2 = l2+4
         read(5,910) (x(l),y(l),z(l),l=l1,l2)
  10  continue
      write(6,915)(x(l),y(l),z(l),l=1,m1)
c  we calculate the exact function values and set up the weights w(i)=
c  (0.01)**(-1) (0.01 is an estimate for the standard deviation of the
c  error in z(i)). at the same time we calculate the mean and maximum
c  errors for the data values.
      sum = 0.
      ermax = 0.
      do 20 i=1,m1
         w(i) = 0.1e03
         exact(i) = testpo(x(i),y(i))
         error = abs(z(i)-exact(i))
         sum = sum+error
         if(error.gt.ermax) ermax = error
  20  continue
      ai = m1
      sum = sum/ai
      write(6,920) sum,ermax
c  we set up the dimension information
      nuest = 15
      nvest = 19
      lwrk1 = 15000
      lwrk2 = 5700
      kwrk = 500
c  we choose a value for eps
      eps = 0.1e-05
c  main loop for the different spline approximations
      do 400 is=1,5
        go to (110,120,130,140,160),is
c  we determine a number of smoothing spline approximations on the unit
c  disk x**2+y**2 <= 1.
c  all the data points are considered.
 110    m = m1
c  we set up the smoothing factor.
        s = 1500.
c  the approximations are not restricted at the boundaries of the disk
        iopt(3) = 0
c  we request c2-continuity at the origin.
        iopt(2) = 2
c  at the first call of polar iopt(1) must be zero.
        iopt(1) = 0
        go to 200
c  iopt(1) = 1 from the second call on
 120    iopt(1) = 1
        s = 200.
        go to 200
 130    s = 170.
        go to 200
c  we determine a smoothing spline approximation on the ellips
c  3*x**2+3*y**2-4*x*y<=1.
c  we only consider the data points inside this domain.
 140    m = m2
        ai = m
c  the given function has a constant value 0.4 at the boundary of the
c  ellips. we calculate new data values by substracting this constant
c  from the old ones.
        do 150 i=1,m
          z(i) = z(i)-0.4
 150    continue
c  given these data we will then determine approximations which are
c  identically zero at the boundary of the ellips.
        iopt(3) = 1
c  we still request c2-continuity at the origin.
        iopt(2) = 2
c  reinitialization for the knots.
        iopt(1) = 0
c  we set up the smoothing factor.
        s = 90.
        go to 250
c  at the last call we will determine the least-squares spline
c  approximation corresponding to the current set of knots
 160    iopt(1) = -1
        go to 250
c  determination of the spline approximation on the disk
 200    call polar(iopt,m,x,y,z,w,rad1,s,nuest,nvest,eps,nu,tu,
     *   nv,tv,u,v,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
        nc = (nu-4)*(nv-4)
c  we calculate the function values at the different points.
        do 220 i=1,m
            f(i) = evapol(tu,nu,tv,nv,c,rad1,x(i),y(i))
 220    continue
        write(6,925) s
        go to 300
c  determination of the spline approximation on the ellips.
 250    call polar(iopt,m,x,y,z,w,rad2,s,nuest,nvest,eps,nu,tu,
     *   nv,tv,u,v,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c  we determine the b-spline coefficients for the spline approximations
c  of the given function.
        nc = (nu-4)*(nv-4)
        do 260 i=1,nc
            c(i) = c(i)+0.4
 260    continue
c  we calculate the function values at the different points.
        do 270 i=1,m
            f(i) = evapol(tu,nu,tv,nv,c,rad2,x(i),y(i))
 270    continue
        if(iopt(1).lt.0) go to 280
        write(6,930) s
        go to 300
 280    write(6,935)
 300    write(6,940) fp,ier
        write(6,945) nu
        write(6,950)
        write(6,955) (tu(i),i=1,nu)
        write(6,960) nv
        write(6,950)
        write(6,955) (tv(i),i=1,nv)
        write(6,965)
        write(6,970) (c(i),i=1,nc)
c  we determine mean and maximum errors.
        sum = 0.
        ermax = 0.
        do 350 i=1,m
          error = abs(f(i)-exact(i))
          sum = sum+error
          if(error.gt.ermax) ermax = error
 350    continue
        sum = sum/ai
        write(6,975)
        write(6,980)
        write(6,915)(x(l),y(l),f(l),l=2,m,3)
        write(6,920) sum,ermax
 400  continue
      stop
c  format statements
 900  format(15h1the input data)
 905  format(1h0,3(3x,1hx,6x,1hy,6x,1hz,5x))
 910  format(12f6.3)
 915  format(1h ,3(3f7.3,2x))
 920  format(14h0mean error = ,f7.4,5x,13hmax. error = ,f7.4)
 925  format(38h0smoothing spline on the disk with s =,f5.0)
 930  format(40h0smoothing spline on the ellips with s =,f5.0)
 935  format(35h0least-squares spline on the ellips)
 940  format(27h0sum of squared residuals =,e15.6,5x,12herror flag =,i5)
 945  format(1x,42htotal number of knots in the u-direction =,i3)
 950  format(1x,22hposition of the knots )
 955  format(5x,8f8.4)
 960  format(1x,42htotal number of knots in the v-direction =,i3)
 965  format(23h0b-spline coefficients )
 970  format(5x,8f9.4)
 975  format(33h0spline values at selected points)
 980  format(1h0,3(3x,1hx,6x,1hy,6x,1hf,5x))
      end
      real function rad1(v)
c  function program rad1 defines in polar coordinates, the boundary of
c  the approximation domain  x**2+y**2<=1.
      real v
      rad1 = 1.
      return
      end
      real function rad2(v)
c  function program rad2 defines in polar coordinates, the boundary of
c  the approximation domain  3*x**2+3*y**2-4*x*y<=1.
      real sqrt,sin,v
      rad2 = 1./sqrt(3.-2.*sin(v+v))
      return
      end
      real function testpo(x,y)
c  function program testpo evaluates the test function for the polar
c  package.
      real x,y
      testpo=(x**2+y**2)/((x+y)**2+0.5)
      return
      end

