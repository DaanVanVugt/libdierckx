cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mncosp : cocosp test program                       cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(10),y(10),w(10),sx(10),s2(10),t(20),c(20),e(20),wrk(550)
      integer iwrk(450)
      logical bind(20)
      integer i,ier,is,j,j3,kwrk,lwrk,m,maxbin,maxtr,n,n4,n6
      real sq
c  the absciss values of the data points.
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)/0.25,
     * 0.5,0.75,1.25,1.75,2.25,2.75,3.25,6.25,12.25/
c  the ordinate values of the data points.
      data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10)/17.0,
     * 15.2,13.8,12.2,11.0,10.1,9.4,8.6,6.1,3.5/
c  m denotes the number of data points.
      m = 10
c  we set up the weights of the data points.
      do 10 i=1,m
         w(i) = 1.
  10  continue
c  we set up the dimension information.
      maxtr = 100
      maxbin = 10
      lwrk = 550
      kwrk = 450
c  we fetch the knots of the cubic spline
      n = 11
      n4 = n-4
      n6 = n-6
c  the interior knots
      t(5) = 1.6
      t(6) = 2.5
      t(7) = 6.0
c  the boundary knots
      do 20 i=1,4
        t(i) = x(1)
        t(i+7) = x(m)
  20  continue
c  loop for the different spline approximations
      do 500 is=1,3
         go to (110,130,150),is
c  a convex spline approximation
 110     write(6,900)
         do 120 j=1,n6
            e(j) = -1.
 120     continue
         go to 200
c  a concave spline approximation (a straight line)
 130     write(6,905)
         do 140 j=1,n6
            e(j) = 1.
 140     continue
         go to 200
c  no convexity/concavity constraints
 150     write(6,910)
         do 160 j=1,n6
            e(j) = 0.
 160     continue
c  we determine the spline approximation.
 200     call cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,
     *     bind,wrk,lwrk,iwrk,kwrk,ier)
c  printing of the results.
         write(6,915) ier
         write(6,920) sq
         write(6,925) n
         write(6,930)
         write(6,935) (t(i),i=1,n)
         write(6,940)
         do 300 j=1,n6
            j3 = j+3
            if(bind(j)) write(6,945) t(j3)
 300     continue
         write(6,950)
         write(6,955) (c(i),i=1,n4)
c  we evaluate the second order derivative of the spline.
         call splder(t,n,c,3,2,x,s2,m,wrk,ier)
         write(6,960)
         do 400 i=1,m
            write(6,965) i,x(i),y(i),sx(i),s2(i)
 400     continue
 500  continue
      stop
c  format statements
 900  format(28h0convex spline approximation)
 905  format(29h0concave spline approximation)
 910  format(35h0unconstrained spline approximation)
 915  format(16h error flag ier=,i2)
 920  format(1x,28hsum of squared residuals sq=,e10.3)
 925  format(1x,24htotal number of knots n=,i2)
 930  format(1x,21hposition of the knots)
 935  format(5x,8f7.2)
 940  format(1x,24hthe knots where s''(x)=0)
 945  format(5x,f7.2)
 950  format(1x,21hb-spline coefficients)
 955  format(5x,4f12.4)
 960  format(3h0 i,6x,4hx(i),5x,4hy(i),4x,7hs(x(i)),3x,9hs''(x(i)))
 965  format(1x,i2,5x,f5.2,5x,f4.1,5x,f5.2,5x,f5.2)
      end
