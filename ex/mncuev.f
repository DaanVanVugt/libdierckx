cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mncuev : curev test program                        cc
cc             evaluation of a closed planar curve                    cc
cc                    x = sx(u) , y = sy(u)                           cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real u(20),t(20),c(40),sp(40)
      integer i,idim,i1,i2,ier,j,jn,j1,j2,j3,j4,k,k1,m,mx,n,nc,nk,nk1
      real ai,per
c  we have a planar curve
      idim = 2
c  set up the dimension information
      nc = 40
      mx = 40
c  set up the points where the curve will be evaluated.
      m = 20
      do 10 i=1,m
        ai = i-1
        u(i) = ai*0.5e-01
  10  continue
c  main loop for the different spline degrees.
      do 50 k=1,5
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
c  evaluate the spline curve
        call curev(idim,t,n,c,nc,k,u,m,sp,mx,ier)
c  print the results.
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
        write(6,930)
        i2 = 0
        j4 = 0
        do 40 j=1,10
          i1 = i2+1
          i2 = i1+1
          j1 = j4+1
          j2 = j1+1
          j3 = j2+1
          j4 = j3+1
          write(6,935) i1,u(i1),sp(j1),sp(j2),i2,u(i2),sp(j3),sp(j4)
  40    continue
  50  continue
      stop
c  format statements.
 900  format(31h0degree of the spline curve k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,12f6.1)
 915  format(1x,30hb-spline coefficients of sx(u))
 920  format(5x,14f5.0)
 925  format(1x,30hb-spline coefficients of sy(u))
 930  format(1x,2(5x,1hi,3x,4hu(i),4x,8hsx(u(i)),4x,8hsy(u(i))))
 935  format(1x,2(i6,f7.2,2f12.5))
      end
