cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                mncurf : curfit test program                        cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(25),y(25),w(25),t(35),c(35),wrk(1000),sp(25)
      integer iwrk(35)
      real ai,fp,s,xb,xe
      integer i,ier,iopt,is,j,k,l,lwrk,l1,l2,m,n,nest,nk1
c  the ordinate values of the data points
      data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),
     * y(12),y(13),y(14),y(15),y(16),y(17),y(18),y(19),y(20),y(21),
     * y(22),y(23),y(24),y(25)/1.0,1.0,1.4,1.1,1.0,1.0,4.0,9.0,13.0,
     * 13.4,12.8,13.1,13.0,14.0,13.0,13.5,10.0,2.0,3.0,2.5,2.5,2.5,
     * 3.0,4.0,3.5/
c  m denotes the number of data points
      m = 25
c  we set up the abscissae and weights of the data points
      do 10 i=1,m
         ai = i-1
         x(i) = ai
         w(i) = 1.0
  10  continue
c  we set up the boundaries of the approximation interval
      xb = x(1)
      xe = x(m)
c  we set up the dimension information.
      nest = 35
      lwrk = 1000
c  loop for the different spline degrees.
      do 400 k=3,5,2
c  loop for the different spline approximations of degree k
         do 300 is=1,7
            go to (110,120,130,140,150,160,170),is
c  we start computing the least-squares polynomial (large value for s).
 110        iopt = 0
            s = 1000.
            go to 200
c  iopt=1 from the second call on
 120        iopt = 1
            s = 60.
            go to 200
c  a smaller value for s to get a closer approximation
 130        s = 10.
            go to 200
c  a larger value for s to get a smoother approximation
 140        s = 30.
            go to 200
c  if a satisfactory fit is obtained  we can calculate a spline of equal
c  quality of fit ( same value for s ) but possibly with fewer knots by
c  specifying iopt=0
 150        s = 30.
            iopt = 0
            go to 200
c  we calculate an interpolating spline
 160        s = 0.
            go to 200
c  finally, we also calculate a least-squares spline function with
c  specified knots
 170        iopt = -1
            j = k+2
            do 180 l=1,7
               ai =3*l
               t(j) = ai
               j = j+1
 180        continue
            n = 9+2*k
 200        call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,
     *       iwrk,ier)
c  printing of the results.
            if(iopt.ge.0) go to 210
            write(6,910) k
            go to 220
 210        write(6,915) k
            write(6,920) s
 220        write(6,925) fp,ier
            write(6,930) n
            write(6,935)
            write(6,940) (t(i),i=1,n)
            nk1 = n-k-1
            write(6,945)
            write(6,950) (c(i),i=1,nk1)
            write(6,955)
c  evaluation of the spline approximation
            call splev(t,n,c,k,x,sp,m,ier)
            do 230 i=1,5
               l1 = (i-1)*5+1
               l2 = l1+4
               write(6,960) (x(l),y(l),sp(l),l=l1,l2)
 230        continue
 300     continue
 400  continue
      stop
 910  format(32h0least-squares spline of degree ,i1)
 915  format(28h0smoothing spline of degree ,i1)
 920  format(20h smoothing factor s=,f5.0)
 925  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
 930  format(1x,24htotal number of knots n=,i3)
 935  format(1x,22hposition of the knots )
 940  format(5x,12f6.1)
 945  format(23h0b-spline coefficients )
 950  format(5x,8f9.4)
 955  format(1h0,5(1x,2hxi,3x,2hyi,2x,5hs(xi),1x))
 960  format(1h ,5(f4.1,1x,f4.1,1x,f4.1,2x))
      end
