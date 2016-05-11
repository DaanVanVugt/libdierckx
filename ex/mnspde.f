cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mnspde : splder test program                       cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(7),y(42),t(20),c(20),wrk(20),d(6)
      integer i,ier,j,k,k1,l,m,n,nk1,nu
      real ai
c  set up the points where the spline derivatives will be evaluated.
      m = 7
      x(1) = 0.
      ai = 0.5e-01
      do 10 i=2,m
        x(i) = x(i-1)+ai
        ai = ai+0.5e-01
  10  continue
      x(m) = 0.1e+01
c  main loop for the different spline degrees.
      do 70 k=1,5
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
c  evaluate the spline derivatives.
        j = 1
        do 40 i=1,k1
c  nu denotes the order of the derivative
          nu = i-1
          call splder(t,n,c,k,nu,x,y(j),m,wrk,ier)
          j = j+m
  40    continue
c  print the results.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        write(6,915)
        write(6,920) (c(i),i=1,nk1)
        write(6,925) (i,i=1,5)
        write(6,930)
        do 60 i=1,m
          j = i
          do 50 l=1,k1
            d(l) = y(j)
            j = j+m
  50      continue
          write(6,935) i,x(i),(d(l),l=1,k1)
  60    continue
  70  continue
      stop
c  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,15f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(21x,5(i4,8x))
 930  format(1x,1hi,6h  x(i),3x,8hs(x(i)) ,5(4x,8hs (x(i))))
 935  format(1x,i1,f6.2,6e12.4)
      end
