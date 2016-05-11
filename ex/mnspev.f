cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mnspev : splev test program                        cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real x(21),y(21),t(20),c(20)
      integer i,i1,i2,ier,j,k,k1,m,n,nk1
      real ai
c  set up the points where the splines will be evaluated.
      m = 21
      do 10 i=1,m
        ai = i-1
        x(i) = ai*0.5e-01
  10  continue
c  main loop for the different spline degrees.
      do 50 k=1,5
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
c  evaluate the spline.
        call splev(t,n,c,k,x,y,m,ier)
c  print the results.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        write(6,915)
        write(6,920) (c(i),i=1,nk1)
        write(6,925) ier
        write(6,930)
        i2 = 0
        do 40 j=1,7
          i1 = i2+1
          i2 = i1+2
          write(6,935) (i,x(i),y(i),i=i1,i2)
  40    continue
  50  continue
      stop
c  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,15f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(1x,16herror flag ier =,i3)
 930  format(1x,3(7x,1hi,3x,4hx(i),4x,7hs(x(i))))
 935  format(1x,3(i8,f7.2,f11.5))
      end
