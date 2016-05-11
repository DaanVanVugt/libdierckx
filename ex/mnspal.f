cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mnspal : spalde test program                       cc
cc    evaluation of a spline function through its polynomial          cc
cc            representation in each knot interval.                   cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(21),y(21),t(20),c(20),d(6),cof(6)
      integer i,i1,i2,ier,j,jj,k,k1,l,l1,m,n,nk1
      real ai,aj,arg,fac,pol,tt,xx
c  set up the points where the splines will be evaluated.
      m = 21
      do 10 i=1,m
        ai = i-1
        x(i) = ai*0.5e-01
  10  continue
c  main loop for the different spline degrees.
      do 100 k=3,5,2
        k1 = k+1
c  n denotes the total number of knots.
        n = 2*k1+4
c  set up the knots of the spline
        j = n
c  the boundary knots
        do 20 i=1,k1
          t(i) = 0.
          t(j) = 0.1e+01
          j = j-1
  20    continue
c  the interior knots
        t(k1+1) = 0.1e+0
        t(k1+2) = 0.3e+0
        t(k1+3) = 0.4e+0
        t(k1+4) = 0.8e+0
c  generate the b-spline coefficients.
        nk1 = n-k1
        do 30 i=1,nk1
          ai = i
          c(i) = 0.1e-01*ai*(ai-0.5e01)
  30    continue
c  print the data for the spline.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        write(6,915)
        write(6,920) (c(i),i=1,nk1)
        l = k
        l1 = k1
c  main loop for the different points of evaluation.
        do 80 i=1,m
          arg = x(i)
c  search for knot interval t(l)<=x(i)<t(l+1).
  40      if(arg.lt.t(l1) .or. l.eq.nk1) go to 60
c  a new knot interval.
          l = l1
          l1 = l+1
          if(t(l).eq.t(l1)) go to 40
          write(6,925) t(l),t(l1)
c  calculate the spline derivatives at the midpoint tt of the interval
          tt = (t(l)+t(l1))*0.5e0
          call spalde(t,n,c,k1,tt,d,ier)
          write(6,930)
          write(6,935) (d(j),j=1,k1)
c  calculate the coefficients cof in the polynomial representation of
c  the spline in the current knot interval,i.e.
c    s(x) = cof(1)+cof(2)*(x-tt)+...+cof(k1)*(x-tt)**k
          fac = 0.1e01
          do 50 j=1,k1
            cof(j) = d(j)/fac
            aj = j
            fac = fac*aj
  50      continue
          write(6,940)
          write(6,935) (cof(j),j=1,k1)
          go to 40
c  evaluate the polynomial
  60      xx = arg-tt
          pol = cof(k1)
          jj = k1
          do 70 j=1,k
            jj = jj-1
            pol = pol*xx+cof(jj)
  70      continue
          y(i) = pol
  80    continue
        write(6,945)
        i2 = 0
        do 90 j=1,7
          i1 = i2+1
          i2 = i1+2
          write(6,950) (i,x(i),y(i),i=i1,i2)
  90    continue
 100  continue
      stop
c  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,15f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(16h0knot interval (,f4.1,1h,,f4.1,1h))
 930  format(1x,49hderivative values at the midpoint of the interval)
 935  format(2x,6e13.5)
 940  format(1x,45hcoefficients in the polynomial representation)
 945  format(1h0,3(7x,1hi,3x,4hx(i),4x,7hs(x(i))))
 950  format(1x,3(i8,f7.2,f11.5))
      end
