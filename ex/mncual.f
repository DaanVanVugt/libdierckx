cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mncual : cualde test program                       cc
cc         evaluation of a closed planar spline curve                 cc
cc                    x = sx(u) , y = sy(u)                           cc
cc            through its polynomial representation                   cc
cc                    in each knot interval.                          cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real t(20),c(40),u(20),sp(40),d(12),cof(2,6)
      integer i,idim,ier,ii,ip,i1,i2,j,jj,jn,j1,j2,j3,j4,k,kk,k1,
     * l,l1,m,n,nc,nd,nk,nk1
      real ai,aj,arg,fac,per,pol,tt,uu
c  we have a planar curve
      idim = 2
c  set up the dimension information
      nc = 40
      nd = 12
c  set up the points where the curve will be evaluated.
      m = 20
      do 10 i=1,m
        ai = i-1
        u(i) = ai*0.5e-01
  10  continue
c  main loop for the different spline degrees.
      do 120 k=3,5,2
c  the order of the spline.
        k1 = k+1
c  n denotes the total number of knots.
        n = 2*k1+4
c  set up the knots of the spline
        t(k1) = 0.
        t(k1+1) = 0.1e0
        t(k1+2) = 0.3e0
        t(k1+3) = 0.4e0
        t(k1+4) = 0.8e0
        t(k1+5) = 0.1e+01
c  fetch the b-spline coefficients for sx(u)
        c(1) = 0.1e+01
        c(2) = 0.3e+01
        c(3) = 0.4e+01
        c(4) = 0.5e+01
        c(5) = -0.1e+01
c  fetch the b-spline coefficients for sy(u)
        c(n+1) = 0.1e+01
        c(n+2) = 0.2e+01
        c(n+3) = -0.3e+01
        c(n+4) = 0.2e+01
        c(n+5) = 0.4e+01
c  incorporate the boundary conditions for periodic splines
        nk = n-k
        per = t(nk)-t(k1)
        do 20 j=1,k
c  the boundary knots
          i1 = nk+j
          i2 = nk-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
c  the boundary coefficients
          jn = j+n
          c(j+5) = c(j)
          c(jn+5) = c(jn)
  20    continue
c  print the data for the spline.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        write(6,915)
        nk1 = n-k1
        write(6,920) (c(i),i=1,nk1)
        write(6,925)
        i1 = n+1
        i2 = n+nk1
        write(6,920) (c(i),i=i1,i2)
        l = k
        l1 = k1
        kk = k1*idim
c  main loop for the different points of evaluation.
        ip = 0
        do 100 i=1,m
          arg = u(i)
c  search for knot interval t(l)<=u(i)<t(l+1).
  40      if(arg.lt.t(l1) .or. l.eq.nk1) go to 70
c  a new knot interval.
          l = l1
          l1 = l+1
          if(t(l).eq.t(l1)) go to 40
          write(6,930) t(l),t(l1)
c  calculate the spline derivatives at the midpoint tt of the interval
          tt = (t(l)+t(l1))*0.5e0
          call cualde(idim,t,n,c,nc,k1,tt,d,nd,ier)
          write(6,935)
          write(6,940) (d(j),j=1,kk)
c  calculate the coefficients cof in the polynomial representation of
c  the spline curve in the current knot interval,i.e.
c    sx(u) = cof(1,1)+cof(1,2)*(u-tt)+...+cof(1,k1)*(u-tt)**k
c    sy(u) = cof(2,1)+cof(2,2)*(u-tt)+...+cof(2,k1)*(u-tt)**k
          fac = 0.1e01
          jj = 0
          do 60 j=1,k1
            do 50 ii=1,idim
              jj = jj+1
              cof(ii,j) = d(jj)/fac
  50        continue
            aj = j
            fac = fac*aj
  60      continue
          write(6,945)
          write(6,950) (cof(1,j),j=1,k1)
          write(6,955)
          write(6,950) (cof(2,j),j=1,k1)
          go to 40
c  evaluate the polynomial curve
  70      uu = arg-tt
          do 90 ii=1,idim
            pol = cof(ii,k1)
            jj = k1
            do 80 j=1,k
              jj = jj-1
              pol = pol*uu+cof(ii,jj)
  80        continue
            ip = ip+1
            sp(ip) = pol
  90      continue
 100    continue
        write(6,960)
        i2 = 0
        j4 = 0
        do 110 j=1,10
          i1 = i2+1
          i2 = i1+1
          j1 = j4+1
          j2 = j1+1
          j3 = j2+1
          j4 = j3+1
          write(6,965) i1,u(i1),sp(j1),sp(j2),i2,u(i2),sp(j3),sp(j4)
 110    continue
 120  continue
      stop
c  format statements.
 900  format(31h0degree of the spline curve k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,12f6.1)
 915  format(1x,30hb-spline coefficients of sx(u))
 920  format(5x,14f5.0)
 925  format(1x,30hb-spline coefficients of sy(u))
 930  format(16h0knot interval (,f4.1,1h,,f4.1,1h))
 935  format(1x,49hcurve derivatives at the midpoint of the interval)
 940  format(1x,3(1x,2e12.4))
 945  format(1x,50hcoefficients in the polynomial represent. of sx(u))
 950  format(2x,6e13.5)
 955  format(1x,50hcoefficients in the polynomial represent. of sy(u))
 960  format(1x,2(5x,1hi,3x,4hu(i),4x,8hsx(u(i)),4x,8hsy(u(i))))
 965  format(1x,2(i6,f7.2,2f12.5))
      end
