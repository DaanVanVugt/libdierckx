cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                mncloc : clocur test program                        cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(38),w(19),u(19),t(40),c(80),wrk(1500),sp(40)
      integer iwrk(40)
      real al,del,fp,s
      integer i,idim,ier,iopt,ipar,is,i1,i2,j,j1,k,l,lwrk,l1,m,mx,
     * n,nc,nest,nk1
c  the data absciss values
      data x(1),x(3),x(5),x(7),x(9),x(11),x(13),x(15),x(17),x(19),x(21),
     * x(23),x(25),x(27),x(29),x(31),x(33),x(35)/-4.7,-7.048,-6.894,
     * -3.75,-1.042,0.938,2.5,3.524,4.511,5.0,4.886,3.524,3.2,1.302,
     * -1.424,-3.0,-3.064,-3.665/
c  the data ordinate values
      data x(2),x(4),x(6),x(8),x(10),x(12),x(14),x(16),x(18),x(20),
     * x(22),x(24),x(26),x(28),x(30),x(32),x(34),x(36)/0.0,2.565,
     * 5.785,6.495,5.909,5.318,4.33,2.957,1.642,0.0,-1.779,-2.957,
     * -5.543,-7.386,-8.075,-5.196,-2.571,-1.334/
c  m denotes the number of data points
      m = 19
c  the first and last data point coincide
      x(2*m-1) = x(1)
      x(2*m) = x(2)
c  we set up the weights and parameter values of the data points
      do 10 i=1,m
         w(i) = 1.0
         al = (i-1)*20
         u(i) = al
  10  continue
c  we set up the dimension information.
      nest = 40
      lwrk = 1500
      nc = 80
      mx = 38
c  we will determine a planar closed curve   x=sx(u) , y=sy(u)
      idim = 2
c  for the first approximations we will use cubic splines
      k = 3
c  we will also supply the parameter values u(i)
      ipar = 1
c  loop for the different approximating spline curves
      do 400 is=1,9
         go to (110,120,130,140,150,160,170,180,190),is
c  we start computing the least-squares point ( s very large)
 110     iopt = 0
         s = 900.
         go to 300
c  iopt =  1 from the second call on
 120     iopt = 1
         s = 10.
         go to 300
c  a smaller value for s to get a closer approximation
 130     s = 0.1
         go to 300
c  a larger value for s to get a smoother approximation
 140     s = 0.5
         go to 300
c  if a satisfactory fit is obtained we can calculate a curve of equal
c  quality of fit (same value for s) but possibly with fewer knots by
c  specifying iopt=0
 150     iopt = 0
         s = 0.5
         go to 300
c  we determine a spline curve with respect to the same smoothing
c  factor s,  but now we let the program determine parameter values u(i)
 160     ipar = 0
         iopt = 0
         s = 0.5
         go to 300
c  we choose a different degree of spline approximation
 170     k = 5
         iopt = 0
         s = 0.5
         go to 300
c  we determine an interpolating curve
 180     s = 0.
         go to 300
c  finally we calculate a least-squares spline curve with specified
c  knots
 190     iopt =-1
         n = 9+2*k
         j = k+2
         del = (u(m)-u(1))*0.125
         do 200 l=1,7
            al = l
            t(j) = u(1)+al*del
            j = j+1
 200     continue
c  determine the approximating closed curve
 300     call clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,
     *    c,fp,wrk,lwrk,iwrk,ier)
c  printing of the results.
         if(iopt.ge.0) go to 310
         write(6,910) k,ipar
         go to 320
 310     write(6,915) k,ipar
         write(6,920) s
 320     write(6,925) fp,ier
         write(6,930) n
         write(6,935)
         if(ipar.eq.1) write(6,940) (t(i),i=1,n)
         if(ipar.eq.0) write(6,950) (t(i),i=1,n)
         nk1 = n-k-1
         write(6,945)
         write(6,950) (c(l),l=1,nk1)
         write(6,955)
         i1 = n+1
         i2 = n+nk1
         write(6,950) (c(l),l=i1,i2)
         write(6,960)
c  we evaluate the spline curve
         call curev(idim,t,n,c,nc,k,u,m,sp,mx,ier)
         do 330 i=1,9
           l = (i-1)*4+1
           l1 = l+1
           j = l+2
           j1 = j+1
           write(6,965) x(l),x(l1),sp(l),sp(l1),x(j),x(j1),sp(j),sp(j1)
 330     continue
 400  continue
      stop
 910  format(38h0least-squares closed curve of degree ,i1,7h  ipar=,i1)
 915  format(34h0smoothing closed curve of degree ,i1,7h  ipar=,i1)
 920  format(20h smoothing factor s=,f7.1)
 925  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
 930  format(1x,24htotal number of knots n=,i3)
 935  format(1x,22hposition of the knots )
 940  format(5x,10f6.0)
 945  format(1x,30hb-spline coefficients of sx(u))
 950  format(5x,8f9.4)
 955  format(1x,30hb-spline coefficients of sy(u))
 960  format(1h0,2(4x,2hxi,7x,2hyi,6x,6hsx(ui),3x,6hsy(ui)))
 965  format(1h ,8f9.4)
      end
