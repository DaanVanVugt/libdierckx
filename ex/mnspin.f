cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                 mnspin : splint test program                       cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real t(20),c(20),wrk(20)
      integer i,j,k,k1,n,nk1,ier
      real a,aint,ak,b,exint,splint
c  as an example we calculate some integrals of the form
c          / b
c         !     (1-x)**k  dx
c      a /
c
c  main loop for the different spline degrees.
      do 30 k=1,5
        k1 = k+1
        ak = k1
c  find the b-spline representation of the polynomial (1-x)**k.
        n = 2*k1
        j = n
        do 10 i=1,k1
          c(i) = 0.
          t(i) = 0.
          t(j) = 0.1e+01
          j = j-1
  10    continue
        c(1) = 0.1e+01
c  insert a number of knots
        call insert(0,t,n,c,k,0.8e0,t,n,c,20,ier)
        call insert(0,t,n,c,k,0.4e0,t,n,c,20,ier)
        call insert(0,t,n,c,k,0.3e0,t,n,c,20,ier)
        call insert(0,t,n,c,k,0.1e0,t,n,c,20,ier)
c  print the data for the spline.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        nk1 = n-k1
        write(6,915)
        write(6,920) (c(i),i=1,nk1)
c  loop for the different integration limits a and b.
        a = 0.
        b = 0.1e+01
        write(6,925)
        do 20 j=1,4
c  calculate the value of the spline integral
          aint = splint(t,n,c,k,a,b,wrk)
c  calculate the exact value of the integral
          exint = ((0.1e01-a)**k1-(0.1e01-b)**k1)/ak
          write(6,930) a,b,aint,exint
          a = a+0.1e0
          b = b-0.3e0
  20    continue
  30  continue
      stop
c  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,15f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(1h0,5x,1ha,6x,1hb,6x,6hsplint,7x,5hexint)
 930  format(1x,2f7.1,2f12.5)
      end
