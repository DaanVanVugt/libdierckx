cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                mnconc : concur test program                        cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(62),w(31),u(31),t(50),c(100),wrk(1400),xx(62),db(6),de(6),
     * cp(24),dd(12),sp(62)
      integer iwrk(50)
      real ai,atan,del,fp,pi,s,sigma,ww
      integer i,ib,idim,ie,ier,iopt,is,i1,i2,j,j1,k,kk,k1,l,lwrk,l1,l2,
     * m,mx,n,nb,nc,ndd,ne,nest,nk1,np
c  the data absciss values
      data x(1),x(3),x(5),x(7),x(9),x(11),x(13),x(15),x(17),x(19),x(21),
     * x(23),x(25),x(27),x(29),x(31),x(33),x(35),x(37),x(39),x(41),
     * x(43),x(45),x(47),x(49),x(51),x(53),x(55),x(57),x(59),x(61)/
     * -3.109,-2.188,-1.351,-0.605,0.093,0.451,0.652,0.701,0.518,0.277,
     * 0.008,-0.291,-0.562,-0.679,-0.637,-0.425,-0.049,0.575,1.334,
     * 2.167,3.206,4.099,4.872,5.710,6.330,6.741,6.928,6.965,6.842,
     * 6.593,6.269/
c  the data ordinate values
      data x(2),x(4),x(6),x(8),x(10),x(12),x(14),x(16),x(18),x(20),
     * x(22),x(24),x(26),x(28),x(30),x(32),x(34),x(36),x(38),x(40),
     * x(42),x(44),x(46),x(48),x(50),x(52),x(54),x(56),x(58),x(60),
     * x(62)/3.040,2.876,2.634,2.183,1.586,1.010,0.382,-0.218,-0.632,
     * -0.879,-0.981,-0.886,-0.642,-0.195,0.373,1.070,1.607,2.165,2.618,
     * 2.905,2.991,2.897,2.615,2.164,1.617,0.977,0.383,-0.194,-0.665,
     * -0.901,-1.010/
c  m denotes the number of data points
      m = 31
c  set up the parameter values for the data points
      pi = atan(1.0)*4.0
      del =pi*0.1
      do 10 i=1,m
         ai = i-11
         u(i) = ai*del
  10  continue
c  the weights are taken as 1/sigma with sigma an estimate of the
c  standard deviation of the data points.
      sigma = 0.04
      ww = 1./sigma
      do 20 i=1,m
        w(i) = ww
  20  continue
c  accordingly, the smoothing factor is chosen as s = m
      s = m
c  we have a planar curve  x = sx(u) , y = sy(u)
      idim = 2
c  begin point derivatives of the curve
      db(1) = -pi
      db(2) = 3.0
      db(3) = 3.0
      db(4) = 0.
      db(5) = 0.
      db(6) = -2.0
c  end point derivatives of the curve
      de(1) = pi*2.0
      de(2) = -1.0
      de(3) = -1.0
      de(4) = 0.
      de(5) = 0.
      de(6) = 2.0
c  we set up the dimension information.
      ndd = 12
      np = 24
      nb = 6
      ne = 6
      nest = 50
      lwrk = 1400
      nc = 100
      mx = 62
c  for the first approximations we will use cubic splines.
      k = 3
c  loop for the different spline curves
      do 500 is=1,7
         go to (110,120,130,140,150,160,170), is
c  no derivative constraints
 110     iopt = 0
         ib = 0
         ie = 0
         go to 200
c  fixed end points.
 120     iopt = 0
         ib = 1
         ie = 1
         go to 200
c  first derivative constraint at the end point
 130     iopt = 0
         ib = 2
         ie = 1
         go to 200
c  first derivative constraints at begin and end point.
 140     iopt = 0
         ib = 2
         ie = 2
         go to 200
c  we choose quintic splines with second derivative constraints.
 150     iopt = 0
         k = 5
         ib = 3
         ie = 3
         go to 200
c  we choose another s-value and continue with the set of knots found at
c  the last call of concur.
 160     iopt = 1
         s = 26.
         go to 200
c  finally we also calculate a least-squares curve with specified knots
 170     iopt = -1
         j = k+2
         do 180 l=1,5
            ai = l-2
            t(j) = ai*pi*0.5
            j = j+1
 180     continue
         n = 7+2*k
c  determination of the spline curve.
 200     call concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,
     *    nest,n,t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
c  printing of the results.
         if(iopt.ge.0) go to 210
         write(6,910) k
         go to 220
 210     write(6,915) k
         write(6,920) s
 220     write(6,925) ib,ie
         write(6,930) fp,ier
         write(6,935) n
         write(6,940)
         write(6,945) (t(i),i=1,n)
         nk1 = n-k-1
         write(6,950)
         write(6,945) (c(l),l=1,nk1)
         write(6,955)
         i1 = n+1
         i2 = n+nk1
         write(6,945) (c(l),l=i1,i2)
c  calculate derivatives at the begin point.
         k1 = k+1
         kk = k1/2
         call cualde(idim,t,n,c,nc,k1,u(1),dd,ndd,ier)
         write(6,960)
         do 300 i=1,kk
            l = i-1
            l1 = l*idim+1
            l2 = l1+1
            write(6,970) l,dd(l1),dd(l2)
 300     continue
c  calculate derivatives at the end point.
         call cualde(idim,t,n,c,nc,k1,u(m),dd,ndd,ier)
         write(6,965)
         do 350 i=1,kk
            l = i-1
            l1 = l*idim+1
            l2 = l1+1
            write(6,970) l,dd(l1),dd(l2)
 350     continue
c  we evaluate the spline curve
         call curev(idim,t,n,c,nc,k,u,m,sp,mx,ier)
         write(6,975)
         do 400 i=1,5
           l = (i-1)*12+3
           l1 = l+1
           j = l+6
           j1 = j+1
           write(6,980) x(l),x(l1),sp(l),sp(l1),x(j),x(j1),sp(j),sp(j1)
 400     continue
 500  continue
      stop
 910  format(31h0least-squares curve of degree ,i1)
 915  format(27h0smoothing curve of degree ,i1)
 920  format(20h smoothing factor s=,f5.0)
 925  format(37h number of derivative constraints ib=,i2,5x,3hie=,i2)
 930  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
 935  format(1x,24htotal number of knots n=,i3)
 940  format(1x,22hposition of the knots )
 945  format(5x,8f8.4)
 950  format(1x,30hb-spline coefficients of sx(u))
 955  format(1x,30hb-spline coefficients of sy(u))
 960  format(1x,30hderivatives at the begin point)
 965  format(1x,28hderivatives at the end point)
 970  format(5x,6horder=,i2,2f9.4)
 975  format(1h0,2(4x,2hxi,7x,2hyi,6x,6hsx(ui),3x,6hsy(ui)))
 980  format(1h ,8f9.4)
      end
