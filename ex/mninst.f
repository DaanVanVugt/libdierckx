cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mninst : insert test program                       cc
cc      application : to find the sum of two periodic splines         cc
cc                 with different sets of knots                       cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real t1(30),c1(30),t2(30),c2(30),x(21),y(21),y1(21),y2(21)
      integer i,ier,iopt,ip,i1,i2,j,j1,j2,k,k1,m,nest,nk,n1,n1k1,n2,n2k1
      real ai,per
c  set up the points where the splines will be evaluated.
      m = 21
      do 10 i=1,m
        ai = i-1
        x(i) = ai*0.5e-01
  10  continue
c  set up the dimension information
      nest = 30
c  main loop for the different spline degrees.
      do 800 k=3,5,2
        k1 = k+1
        do 700 ip=1,2
c  if iopt = 1 the splines will be considered as periodic splines.
c  if iopt = 0 they will be considered as ordinary splines.
        iopt = 2-ip
        if(iopt.eq.0) write(6,900)
        if(iopt.ne.0) write(6,905)
        write(6,910) k
c  fetch the knots and b-spline coefficients of the first spline s1(x).
        n1 = 2*k1+5
        n1k1 = n1-k1
        t1(k1) = 0.
        t1(k1+1) = 0.2e0
        t1(k1+2) = 0.3e0
        t1(k1+3) = 0.4e0
        t1(k1+4) = 0.7e0
        t1(k1+5) = 0.9e0
        t1(k1+6) = 0.1e+01
        c1(1) = 0.1e+01
        c1(2) = 0.2e+01
        c1(3) =-0.1e+01
        c1(4) = 0.3e+01
        c1(5) = 0.3e+01
        c1(6) =-0.3e+01
c  fetch the knots and b-spline coefficients of the second spline s2(x).
        n2 = 2*k1+6
        n2k1 = n2-k1
        t2(k1) = 0.
        t2(k1+1) = 0.1e0
        t2(k1+2) = 0.2e0
        t2(k1+3) = 0.3e0
        t2(k1+4) = 0.4e0
        t2(k1+5) = 0.7e0
        t2(k1+6) = 0.8e0
        t2(k1+7) = 0.1e+01
        c2(1) = 0.2e+01
        c2(2) =-0.2e+01
        c2(3) = 0.1e+01
        c2(4) =-0.3e+01
        c2(5) = 0.4e+01
        c2(6) = 0.4e+01
        c2(7) = 0.4e+01
c  incorporate the boundary conditions for periodic splines.
        per = 0.1e+01
        nk = n1-k
        do 20 j=1,k
          i1 = nk+j
          i2 = nk-j
          j1 = k1+j
          j2 = k1-j
c  the boundary knots
          t1(i1) = t1(j1)+per
          t1(j2) = t1(i2)-per
          t2(i1+1) = t2(j1)+per
          t2(j2) = t2(i2+1)-per
c  the boundary coefficients
          c1(j+6) = c1(j)
          c2(j+7) = c2(j)
  20    continue
        if(iopt.ne.0) go to 100
c  if iopt=0 we insert k knots at the boundaries of the interval to
c  find the representation with coincident boundary knots
        do 40 j=1,k
          call insert(iopt,t1,n1,c1,k,0.1e+01,t1,n1,c1,nest,ier)
          n1 = n1-1
          call insert(iopt,t1,n1,c1,k,0.,t1,n1,c1,nest,ier)
          n1 = n1-1
          do 30 i=1,n1
            t1(i) = t1(i+1)
            c1(i) = c1(i+1)
  30      continue
  40    continue
        do 60 j=1,k
          call insert(iopt,t2,n2,c2,k,0.1e+01,t2,n2,c2,nest,ier)
          n2 = n2-1
          call insert(iopt,t2,n2,c2,k,0.,t2,n2,c2,nest,ier)
          n2 = n2-1
          do 50 i=1,n2
            t2(i) = t2(i+1)
            c2(i) = c2(i+1)
  50      continue
  60    continue
c  print knots and b-spline coefficients of the two splines.
 100    write(6,915)
        write(6,920) (t1(i),i=1,n1)
        write(6,925)
        write(6,930) (c1(i),i=1,n1k1)
        write(6,935)
        write(6,920) (t2(i),i=1,n2)
        write(6,940)
        write(6,930) (c2(i),i=1,n2k1)
c  evaluate the two splines
        call splev(t1,n1,c1,k,x,y1,m,ier)
        call splev(t2,n2,c2,k,x,y2,m,ier)
c  insert the knots of the second spline into those of the first one
        call insert(iopt,t1,n1,c1,k,0.1e0,t1,n1,c1,nest,ier)
        call insert(iopt,t1,n1,c1,k,0.8e0,t1,n1,c1,nest,ier)
c  insert the knots of the first spline into those of the second one
        call insert(iopt,t2,n2,c2,k,0.9e0,t2,n2,c2,nest,ier)
c  print the knots and coefficients of the splines in their new
c  representation
        n1k1 = n1-k1
        write(6,945)
        write(6,920) (t1(i),i=1,n1)
        write(6,925)
        write(6,930) (c1(i),i=1,n1k1)
        write(6,940)
        write(6,930) (c2(i),i=1,n1k1)
c  find the coefficients of the sum of the two splines.
        do 110 i=1,n1k1
          c1(i) = c1(i)+c2(i)
 110    continue
        write(6,950)
        write(6,930) (c1(i),i=1,n1k1)
c  evaluate this new spline and compare results
        call splev(t1,n1,c1,k,x,y,m,ier)
        write(6,955)
        do 200 i=1,m
          write(6,960) i,x(i),y1(i),y2(i),y(i)
 200    continue
 700  continue
 800  continue
      stop
c  format statements.
 900  format(41h0insertion algorithm for ordinary splines)
 905  format(41h0insertion algorithm for periodic splines)
 910  format(1x,25hdegree of the splines k =,i2)
 915  format(1x,30hposition of the knots of s1(x))
 920  format(5x,15f5.1)
 925  format(1x,30hb-spline coefficients of s1(x))
 930  format(5x,8f9.5)
 935  format(1x,30hposition of the knots of s2(x))
 940  format(1x,30hb-spline coefficients of s2(x))
 945  format(1x,37hposition of the knots after insertion)
 950  format(1x,33hb-spline coefficients of s1+s2(x))
 955  format(3h0 i,6x,1hx,7x,5hs1(x),7x,5hs2(x),6x,8hs1+s2(x))
 960  format(1x,i2,f8.2,3f12.5)
      end
