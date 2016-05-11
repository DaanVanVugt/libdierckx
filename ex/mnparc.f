cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                mnparc : parcur test program                        cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(64),w(32),u(32),t(40),c(80),wrk(1200),sp(64)
      integer iwrk(40)
      real al,del,fp,s,ub,ue
      integer i,idim,ier,iopt,ipar,is,i1,i2,j,j1,k,l,lwrk,l1,m,mx,
     * n,nc,nest,nk1
c  the data parameter values
      data u(1),u(2),u(3),u(4),u(5),u(6),u(7),u(8),u(9),u(10),u(11),
     * u(12),u(13),u(14),u(15),u(16),u(17),u(18),u(19),u(20),u(21),
     * u(22),u(23),u(24),u(25),u(26),u(27),u(28),u(29),u(30),u(31),
     * u(32)/120.,128.,133.,136.,138.,141.,144.,146.,149.,151.,154.,
     * 161.,170.,180.,190.,200.,210.,220.,230.,240.,250.,262.,269.,
     * 273.,278.,282.,287.,291.,295.,299.,305.,315./
c  the data absciss values
      data x(1),x(3),x(5),x(7),x(9),x(11),x(13),x(15),x(17),x(19),x(21),
     * x(23),x(25),x(27),x(29),x(31),x(33),x(35),x(37),x(39),x(41),
     * x(43),x(45),x(47),x(49),x(51),x(53),x(55),x(57),x(59),x(61),
     * x(63)/-1.5141,-2.0906,-1.9253,-0.8724,-0.3074,-0.5534,0.0192,
     * 1.2298,2.5479,2.4710,1.7063,1.1183,0.5534,0.4727,0.3574,0.1998,
     * 0.2882,0.2613,0.2652,0.2805,0.4112,0.9377,1.3527,1.5564,1.6141,
     * 1.6333,1.1567,0.8109,0.2498,-0.2306,-0.7571,-1.1222/
c  the data ordinate values
      data x(2),x(4),x(6),x(8),x(10),x(12),x(14),x(16),x(18),x(20),
     * x(22),x(24),x(26),x(28),x(30),x(32),x(34),x(36),x(38),x(40),
     * x(42),x(44),x(46),x(48),x(50),x(52),x(54),x(56),x(58),x(60),
     * x(62),x(64)/0.5150,1.3412,2.6094,3.2358,2.7401,2.7823,3.5932,
     * 3.8353,2.5863,1.3105,0.6841,0.2575,0.2460,0.3689,0.2460,0.2998,
     * 0.3651,0.3343,0.3881,0.4573,0.5918,0.7110,0.4035,0.0769,-0.3920,
     * -0.8570,-1.3412,-1.5641,-1.7409,-1.7178,-1.2989,-0.5572/
c  m denotes the number of data points
      m = 32
c  we set up the weights of the data points
      do 10 i=1,m
         w(i) = 1.0
  10  continue
c  we set up the dimension information.
      nest = 40
      lwrk = 1200
      nc = 80
      mx = 64
c  we will determine a planar curve   x=sx(u) , y=sy(u)
      idim = 2
c  for the first approximations we will use cubic splines
      k = 3
c  we will also supply the parameter values u(i)
      ipar = 1
      ub = 120.
      ue = 320.
c  loop for the different approximating spline curves
      do 400 is=1,9
         go to (110,120,130,140,150,160,170,180,190),is
c  we start computing a polynomial curve ( s very large)
 110     iopt = 0
         s = 100.
         go to 300
c  iopt =  1 from the second call on
 120     iopt = 1
         s = 1.0
         go to 300
c  a smaller value for s to get a closer approximation
 130     s = 0.05
         go to 300
c  a larger value for s to get a smoother approximation
 140     s = 0.25
         go to 300
c  if a satisfactory fit is obtained we can calculate a curve of equal
c  quality of fit (same value for s) but possibly with fewer knots by
c  specifying iopt=0
 150     iopt = 0
         s = 0.25
         go to 300
c  we determine a spline curve with respect to the same smoothing
c  factor s,  but now we let the program determine parameter values u(i)
 160     ipar = 0
         iopt = 0
         s = 0.25
         go to 300
c  we choose a different degree of spline approximation
 170     k = 5
         iopt = 0
         s = 0.25
         go to 300
c  we determine an interpolating curve
 180     s = 0.
         go to 300
c  finally we calculate a least-squares spline curve with specified
c  knots
 190     iopt =-1
         n = 9+2*k
         j = k+2
         del = (ue-ub)*0.125
         do 200 l=1,7
            al = l
            t(j) = ub+al*del
            j = j+1
 200     continue
c  determine the approximating curve
 300     call parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,
     *    nc,c,fp,wrk,lwrk,iwrk,ier)
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
         do 330 i=1,8
           l = (i-1)*8+3
           l1 = l+1
           j = l+4
           j1 = j+1
           write(6,965) x(l),x(l1),sp(l),sp(l1),x(j),x(j1),sp(j),sp(j1)
 330     continue
 400  continue
      stop
 910  format(31h0least-squares curve of degree ,i1,7h  ipar=,i1)
 915  format(27h0smoothing curve of degree ,i1,7h  ipar=,i1)
 920  format(20h smoothing factor s=,f7.2)
 925  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
 930  format(1x,24htotal number of knots n=,i3)
 935  format(1x,22hposition of the knots )
 940  format(5x,10f6.0)
 945  format(1x,30hb-spline coefficients of sx(u))
 950  format(5x,8f8.4)
 955  format(1x,30hb-spline coefficients of sy(u))
 960  format(1h0,2(4x,2hxi,7x,2hyi,6x,6hsx(ui),3x,6hsy(ui)))
 965  format(1h ,8f9.4)
      end
