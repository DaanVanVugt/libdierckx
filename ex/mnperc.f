cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                mnperc : percur test program                        cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(27),y(27),w(27),t(37),c(37),wrk(1400),sp(27)
      integer iwrk(37)
      real al,fp,s
      integer i,ier,iopt,is,j,k,l,lwrk,l1,l2,m,m1,n,nest,nk1
c  the data absciss values
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),
     * x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20),x(21),
     * x(22),x(23),x(24),x(25),x(26)/0.0,3.922,7.843,11.765,15.686,
     * 19.608,23.509,27.451,31.373,35.294,39.216,43.137,47.059,50.980,
     * 54.902,58.824,62.745,66.667,70.588,74.510,78.431,82.353,86.275,
     * 90.196,94.118,98.039/
c  the data ordinate values
      data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),
     * y(12),y(13),y(14),y(15),y(16),y(17),y(18),y(19),y(20),y(21),
     * y(22),y(23),y(24),y(25),y(26)/10.099,14.835,21.453,25.022,22.427,
     * 22.315,22.070,19.673,16.754,13.983,11.973,12.286,16.129,21.560,
     * 28.041,39.205,59.489,72.559,75.960,79.137,75.925,68.809,55.758,
     * 39.915,22.006,12.076/
c  m denotes the number of data points
      m = 27
c  the period of the spline is determined by x(m)
      x(m) = 100.
      y(m) = y(1)
c  we set up the weights of the data points
      m1 = m-1
      do 10 i=1,m1
         w(i) = 1.0
  10  continue
c  we set up the dimension information.
      nest = 37
      lwrk = 1400
c  loop for the different spline degrees.
      do 400 k=3,5,2
c  loop for the different spline approximations of degree k
         do 300 is=1,7
            go to (110,120,130,140,150,160,170),is
c  we start computing the least-squares constant (large value for s).
 110        iopt = 0
            s = 65000.
            go to 200
c  iopt=1 from the second call on
 120        iopt = 1
            s = 500.
            go to 200
c  a smaller value for s to get a closer approximation
 130        s = 5.
            go to 200
c  a larger value for s to get a smoother approximation
 140        s = 20.
            go to 200
c  if a satisfactory fit is obtained  we can calculate a spline of equal
c  quality of fit ( same value for s ) but possibly with fewer knots by
c  specifying iopt=0
 150        s = 20.
            iopt = 0
            go to 200
c  we calculate an interpolating periodic spline.
 160        s = 0.
            go to 200
c  finally, we also calculate a least-squares periodic spline function
c  with specified knots.
 170        iopt = -1
            n = 11+2*k
            j = k+2
            do 180 l=1,9
               al = l*10
               t(j) = al
               j = j+1
 180        continue
c  determine the periodic spline approximation
 200        call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,
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
            do 230 i=1,9
               l1 = (i-1)*3+1
               l2 = l1+2
               write(6,960) (x(l),y(l),sp(l),l=l1,l2)
 230        continue
 300     continue
 400  continue
      stop
 910  format(41h0least-squares periodic spline of degree ,i1)
 915  format(37h0smoothing periodic spline of degree ,i1)
 920  format(20h smoothing factor s=,f7.0)
 925  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
 930  format(1x,24htotal number of knots n=,i3)
 935  format(1x,22hposition of the knots )
 940  format(5x,8f8.3)
 945  format(23h0b-spline coefficients )
 950  format(5x,8f8.4)
 955  format(1h0,3(3x,2hxi,6x,2hyi,4x,5hs(xi),3x))
 960  format(1h ,3(f7.3,1x,f7.3,1x,f7.3,2x))
      end
