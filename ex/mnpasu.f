cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc               mnpasu : parsur test program                         cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real u(21),v(11),f(693),tu(27),tv(17),c(900),wrk(2000),z(693),
     * wk(128)
      integer iwrk(80),iw(32),ipar(2)
      real ai,fp,s
      integer kwrk,lwrk,m,mu,mv,j0,j1,j2,j3,nc,nu,nuest,nv,nvest,
     * i,idim,ier,is,iopt,j,l
c  we generate the u-coordinates of the grid.
      mu = 21
      do 10 i=1,mu
        u(i) = i-1
  10  continue
c  we generate the v-coordinates of the grid.
      mv = 11
      do 20 i=1,mv
        v(i) = u(2*i-1)
  20  continue
c  the dimension of the surface
      idim = 3
c  we fetch and print the surface co-ordinates at the grid points
      write(6,900)
      write(6,905) (v(i),i=1,mv)
      write(6,910)
      m = mu*mv
      j0 = 0
      do 40 i=1,mu
        write(6,915) u(i)
        j1 = j0
        do 30 l=1,idim
          j2 = j1+1
          j3 = j1+mv
          read(5,920) (f(j),j=j2,j3)
          write(6,925) (f(j),j=j2,j3)
          j1 = j1+m
  30    continue
        j0 = j0+mv
  40  continue
c  we set up the dimension information
      nuest = 27
      nvest = 17
      lwrk = 2000
      kwrk = 80
c  main loop for the different spline approximations
      do 300 is=1,4
        go to (110,120,130,140),is
c  a smoothing surface with no periodicity conditions
 110    iopt = 0
        s = 0.07
        ipar(1) = 0
        ipar(2) = 0
        go to 200
c  a smoothing surface periodic in the v-variable
 120    ipar(2) = 1
        go to 200
c  a smoothing surface periodic in both variables
 130    ipar(1) = 1
        go to 200
c  finally we also calculate a least-squares spline surface
c  with specified knots.
 140    iopt = -1
        nu = 11
        nv = 11
        j = 5
        do 150 i=1,3
          ai = 5*i
          tu(j) = ai
          tv(j) = tu(j)
          j = j+1
 150    continue
c  determination of the spline surface.
 200    call parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,
     *   nvest,nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c  printing of the fitting results.
        if(iopt.ge.0) go to 210
        write(6,935) ipar(1),ipar(2)
        go to 220
 210    write(6,940) ipar(1),ipar(2)
        write(6,945) s
 220    write(6,950) fp,ier
        write(6,955) nu
        write(6,960)
        write(6,965) (tu(i),i=1,nu)
        write(6,970) nv
        write(6,960)
        write(6,965) (tv(i),i=1,nv)
        nc = (nu-4)*(nv-4)
        write(6,975)
        j1 = 0
        do 230 l=1,idim
          j0 = j1+1
          j1 = j1+nc
          write(6,980) (c(j),j=j0,j1)
 230    continue
c  evaluation of the spline surface.
        call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,z,693,
     *   wk,128,iw,32,ier)
        write(6,985)
        write(6,905) (v(i),i=1,mv,2)
        write(6,910)
        j0 = 0
        do 250 i=1,mu,4
          write(6,915) u(i)
          j1 = j0
          do 240 l=1,idim
            j2 = j1+1
            j3 = j1+mv
            write(6,925) (z(j),j=j2,j3,2)
            j1 = j1+m
 240      continue
          j0 = j0+mv*4
 250    continue
 300  continue
      stop
c  format statements.
 900  format(15h1the input data)
 905  format(1h0,2x,1hv,11(3x,f4.1))
 910  format(1h ,1x,1hu)
 915  format(1h ,f4.1)
 920  format(11f7.3)
 925  format(5x,11f7.3)
 935  format(37h0least-squares surface of periodicity,2i3)
 940  format(33h0smoothing surface of periodicity,2i3)
 945  format(20h smoothing factor s=,f8.2)
 950  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 955  format(1x,42htotal number of knots in the u-direction =,i3)
 960  format(1x,22hposition of the knots )
 965  format(5x,10f6.2)
 970  format(1x,42htotal number of knots in the v-direction =,i3)
 975  format(23h0b-spline coefficients )
 980  format(5x,8f9.4)
 985  format(1h0,37hspline values at selected grid points)
      end
