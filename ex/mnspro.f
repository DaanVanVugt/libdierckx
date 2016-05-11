cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mnspro : sproot test program                       cc
cc      application : to find the intersection of a planar            cc
cc      cubic spline curve   x = sx(u)   y = sy(u)   with             cc
cc      a straight line   alfa*x + beta*y = gamma                     cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real t(13),c(26),zero(20),sp(40),cc(13)
      integer i,idim,ier,is,i1,i2,j,k,k1,l1,l2,m,mest,n,nc,nk1
      real alfa,beta,gamma,per
c  we have a planar curve
      idim = 2
c  we have a cubic spline curve.
      k = 3
      k1 = k+1
c  set up the dimension information
      nc = 26
      mest = 20
c  n denotes the total number of knots.
      n = 13
c  set up the knots of the spline curve
      t(4) = 0.
      t(5) = 0.2e0
      t(6) = 0.3e0
      t(7) = 0.5e0
      t(8) = 0.6e0
      t(9) = 0.7e0
      t(10) = 0.1e+01
c  fetch the b-spline coefficients for sx(u)
      c(1) = 0.1e+01
      c(2) = 0.3e+01
      c(3) = 0.4e+01
      c(4) = 0.5e+01
      c(5) = 0.3e+01
      c(6) = -0.1e+01
c  fetch the b-spline coefficients for sy(u)
      c(14) = 0.1e+01
      c(15) = 0.2e+01
      c(16) = -0.3e+01
      c(17) = 0.2e+01
      c(18) = 0.1e+01
      c(19) = 0.4e+01
c  we have a closed curve.
c  incorporate the boundary conditions for periodic splines
      per = t(10)-t(4)
      do 10 i=1,3
c  the boundary knots
        t(i) = t(i+6)-per
        t(i+10) = t(i+4)+per
c  the boundary coefficients
        c(i+6) = c(i)
        c(i+19) = c(i+13)
  10  continue
c  print the data of the spline curve.
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
c  loop for the different lines.
      do 200 is=1,5
        go to (110,120,130,140,150),is
c  fetch the parameters of the straight line.
 110    alfa = 0.
        beta = 0.1e+01
        gamma = 0.
        go to 160
 120    alfa = 0.1e+01
        beta = 0.
        go to 160
 130    beta = -0.1e+01
        go to 160
 140    alfa = 0.4e0
        beta = 0.3e0
        gamma = 0.12e+01
        go to 160
 150    beta = 0.4e0
        gamma = 0.
c  print the parameters of the straight line.
 160    write(6,930) alfa,beta,gamma
c  calculate the coefficients of s(u) = sx(u)*alfa + sy(u)*beta - gamma
        do 170 i=1,nk1
          j = i+n
          cc(i) = alfa*c(i)+beta*c(j)-gamma
 170    continue
c  find the zeros of s(u)
        call sproot(t,n,cc,zero,mest,m,ier)
        write(6,935) m
        if(m.eq.0) go to 200
c  find the intersection points
        call curev(idim,t,n,c,nc,k,zero,m,sp,nc,ier)
c  print the intersection points
        write(6,940)
        l2 = 0
        do 180 i=1,m
          l1 = l2+1
          l2 = l1+1
          write(6,945) i,zero(i),sp(l1),sp(l2)
 180    continue
 200  continue
      stop
c  format statements.
 900  format(31h0degree of the spline curve k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,7f6.1)
 915  format(1x,30hb-spline coefficients of sx(u))
 920  format(5x,14f5.0)
 925  format(1x,30hb-spline coefficients of sy(u))
 930  format(18h0intersection with,f6.1,5h *x +,f5.1,5h *y =,f5.1)
 935  format(1x,33hnumber of intersection points m =,i3)
 940  format(6x,1hi,7x,4hu(i),5x,8hsx(u(i)),4x,8hsy(u(i)))
 945  format(1x,i6,3f12.5)
      end
