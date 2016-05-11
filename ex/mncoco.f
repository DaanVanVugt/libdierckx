cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mncoco : concon test program                       cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(16),y(16),w(16),v(16),sx(16),s2(16),t(20),c(20),wrk(550),
     * s(3)
      integer iwrk(450)
      logical bind(20)
      integer i,ier,iopt,is,j,j3,kwrk,lwrk,m,maxbin,maxtr,n,nest,n4,n6
      real sq
c  the absciss values of the data points.
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),
     * x(12),x(13),x(14),x(15),x(16)/0.1,0.3,0.5,0.7,0.9,1.25,1.75,
     * 2.25,2.75,3.5,4.5,5.5,6.5,7.5,8.5,9.5/
c  the ordinate values of the data points.
      data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),
     * y(12),y(13),y(14),y(15),y(16)/0.124,0.234,0.256,0.277,0.278,
     * 0.291,0.308,0.311,0.315,0.322,0.317,0.326,0.323,0.321,0.322,
     * 0.328/
c  m denotes the number of data points.
      m = 16
c  we set up the weights of the data points.
      do 10 i=1,m
         w(i) = 1.
  10  continue
      w(1) = 10.
      w(2) = 3.
      w(16) = 10.
c  we will determine concave approximations
      do 20 i=1,m
         v(i) = 1.
  20  continue
c  we set up the dimension information.
      nest = 20
      maxtr = 100
      maxbin = 10
      lwrk = 550
      kwrk = 450
c  we set up the different s-values.
      s(1) = 0.2
      s(2) = 0.04
      s(3) = 0.0002
c  initialization.
      iopt = 0
c  loop for the different spline approximations
      do 300 is=1,3
c  we determine the concave spline approximation.
         call concon(iopt,m,x,y,w,v,s(is),nest,maxtr,maxbin,n,t,c,sq,sx,
     *     bind,wrk,lwrk,iwrk,kwrk,ier)
c  printing of the results.
         write(6,900) s(is),ier
         write(6,905) sq
         write(6,910) n
         write(6,915)
         write(6,920) (t(i),i=1,n)
         write(6,925)
         n6 = n-6
         do 100 j=1,n6
            j3 = j+3
            if(bind(j)) write(6,930) t(j3)
 100     continue
         write(6,935)
         n4 = n-4
         write(6,940) (c(i),i=1,n4)
c  we evaluate the second order derivative of the spline.
         call splder(t,n,c,3,2,x,s2,m,wrk,ier)
         write(6,945)
         do 200 i=1,m
            write(6,950) i,x(i),y(i),sx(i),s2(i)
 200     continue
c  iopt=1 from the second call on.
         iopt = 1
 300  continue
      stop
c  format statements
 900  format(48h0upper limit for the sum of squared residuals s=,
     * e8.1,5x,15herror flag ier=,i2)
 905  format(1x,28hsum of squared residuals sq=,e10.3)
 910  format(1x,24htotal number of knots n=,i2)
 915  format(1x,21hposition of the knots)
 920  format(5x,8f7.2)
 925  format(1x,24hthe knots where s''(x)=0)
 930  format(5x,f7.2)
 935  format(1x,21hb-spline coefficients)
 940  format(5x,4f12.6)
 945  format(3h0 i,5x,4hx(i),6x,4hy(i),4x,7hs(x(i)),4x,9hs''(x(i)))
 950  format(1x,i2,5x,f4.2,5x,f5.3,5x,f5.3,5x,f8.4)
      end
