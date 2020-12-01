!AC -- read KORAL hdf5 files 11/20
!note: in koral convention, x3 is slowest changing, x2 is fastest
!make sure this works in 2d and in 3d!
module fluid_model_koralh5

  use, intrinsic :: iso_c_binding
  use HDF5

  use class_four_vector
  use interpolate, only: interp
  use kerr, only: kerr_metric, lnrf_frame, uks2ubl, bl2ks
  use phys_constants, only: pi,c2,k,mp,m
  use math, only: zbrent
  implicit none

  namelist /harm/  dfile, hfile, nt, indf

  logical :: shortfile
  character(len=500) :: dfile, hfile
  integer :: nx1, nx2, nx3,n, ndumps, nt, indf, has_electrons
  integer :: ndim=3
  integer :: minpolecell=2
  integer, dimension(:), allocatable :: dumps
  real :: tstep=10., toffset=0.
  real(8) :: dx1, dx2, dx3, gam, startx1, startx2, startx3
  real :: asim
  real :: r0, h0, my1, my2, pp
  real :: rmin, rmax, x1min, x1max, rbrk,fjet,fdisk,runi,rcj,rcd,rdj,rdd,al1,al2,rcyl,ncyl
  character(len=50) :: runmet, outmet
  real(8) :: scalefac
  real, dimension(:), allocatable :: t
  real, dimension(:), allocatable :: uniqx1, uniqx2, uniqx2_ext, uniqx3, uniqr, uniqph
  real, dimension(:), allocatable :: r_arr, th_arr, ph_arr
  real, dimension(:,:), allocatable :: r_th_arr
  real, dimension(:), allocatable :: rho_arr, p_arr, Be_arr
  real, dimension(:), allocatable :: te_arr, ti_arr
  real, dimension(:), allocatable :: u0_arr, vrl_arr, vpl_arr, vtl_arr
  real, dimension(:), allocatable :: b0_arr, br_arr, bth_arr, bph_arr

  interface init_koralh5_data
    module procedure init_koralh5_data
  end interface

  interface del_koralh5_data
    module procedure del_koralh5_data
  end interface

  interface initialize_koralh5_model
    module procedure initialize_koralh5_model
  end interface

  interface koralh5_vals
    module procedure koralh5_vals
  end interface

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Rowan 2019 temperature prescription 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine michael_e(rho,ttotcgs,bhl,tecgs) ! Rowan heating with guide field.
        real(kind=8), intent(in), dimension(:) :: rho,ttotcgs,bhl ! all in cgs, except bsq in HL!
        real(kind=8), intent(inout), dimension(size(rho)) :: tecgs
        real(kind=8), dimension(size(rho)) :: tiguess, teguess,thetai, thetae, betai,sigmai,sigmaw,wfaci, wface
        real(kind=8), dimension(size(rho)) :: adiabi, adiabe, delta,urat,tratout,betamax,bsq
        real(kind=8) :: bg, tratio

        bg = 0.3 ! fixed for now
        tratio = 1 ! electron to ion, fixed input for now

        tiguess = ttotcgs / (1. + tratio)
        teguess = tiguess*tratio
        thetai = tiguess*k/(mp*c2)
        thetae = teguess*k/(m*c2)
        adiabi = (10.+20*thetai)/(6+15.*thetai)
        adiabe = (10.+20*thetae)/(6+15.*thetae)
        wfaci = adiabi/(adiabi-1)
        wface = adiabe/(adiabe-1)

        bsq = bhl*bhl
        betai= 2*rho*k*tiguess/mp/(bsq)
        sigmai= bsq/(rho*c2)
        sigmaw = sigmai/(1.+ thetai*(wfaci + wface*tratio))
        betamax = 0.5/(sigmaw*(1+tratio))

        delta = 0.5*1.7*(tanh(0.33*bg)-0.4)*tanh(((1-betai/betamax)**1.5)/((0.42+tratio)*sigmaw**0.3)) + 0.5
        where(isnan(delta))
          delta = 0.5
        else where(delta<0)
          delta = 0.
        end where

        urat = 1./delta - 1.
        tratout =urat*(adiabi-1.)/(adiabe-1) !ion to electron temperature ratio
        !write(6,*) 'trat',minval(tratout),maxval(tratout)
        tecgs = ttotcgs/(1.+tratout)

    end subroutine michael_e

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Zhdankin 2019 temperature prescription )
    ! assume equal temperatures as input -- not strictly consistent!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine vladimir_e(ttotcgs,tecgs) 
        real(kind=8), intent(in), dimension(:) :: ttotcgs
        real(kind=8), intent(inout), dimension(size(ttotcgs)) :: tecgs
        real(kind=8), dimension(size(ttotcgs)) :: tiguess, thetai,thetae, gyrorat
        real(kind=8), dimension(size(ttotcgs)) :: urat, tratout,adiabi,adiabe
        real(kind=8) :: tratio

        tratio = 1 ! electron to ion, fixed input for now

        tiguess = ttotcgs / (1. + tratio)
        thetai = tiguess*k/(mp*c2)
        thetae = tratio*thetai*1836.15267
        adiabi = (10.+20*thetai)/(6+15.*thetai)
        adiabe = (10.+20*thetae)/(6+15.*thetae)

        gyrorat = 1836.15267*sqrt((meangamma(thetai)**2-1.)/(meangamma(thetae)**2-1.))
        urat = (gyrorat)**(2./3.) !ion-to-electron energy ratio
        tratout = urat*(adiabi-1.)/(adiabe-1.) !ion to electron temperature ratio

        !write(6,*) 'trat',minval(tratout),maxval(tratout)
        tecgs = ttotcgs/(1.+tratout)

    end subroutine vladimir_e

    function meangamma(theta) result(gammas)
        real(8), intent(in), dimension(:) :: theta
        real(8), dimension(size(theta)) :: gammas,xvals
        real(8), dimension(61) :: logthetaarr, gammarr

        logthetaarr = (/-3., -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2., &
        -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1., -0.9, &
        -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, &
        0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, &
        1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3./)

        gammarr = (/1.0015, 1.00189, 1.00238, 1.003, 1.00378, 1.00476, 1.006, 1.00756, &
        1.00954, 1.01203, 1.01519, 1.01918, 1.02424, 1.03066, 1.03883, &
        1.04925, 1.06257, 1.07966, 1.10165, 1.13008, 1.16699, 1.2151, 1.2781, &
        1.36087, 1.46995, 1.6139, 1.80386, 2.05412, 2.38271, 2.81217, &
        3.37044, 4.09199, 5.01943, 6.20561, 7.71639, 9.63422, 12.0626, &
        15.1319, 19.006, 23.8917, 30.0494, 37.8071, 47.5782, 59.8828, &
        75.3764, 94.8841, 119.445, 150.366, 189.295, 238.305, 300.005, &
        377.682, 475.471, 598.581, 753.568, 948.685, 1194.32, 1503.56, &
        1892.87, 2382.99, 3000./)

        where(theta.lt.1.d-3)
            gammas = 1.5*theta+1.
        else where(theta.gt.1.d3)
            gammas = 3*theta
        elsewhere
            xvals = log10(theta)
            gammas = interp1d_arr(logthetaarr,gammarr,xvals) 
        end where
    end function meangamma

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate gammagas from Te, Ti (in K)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_gammagas(te, ti, gammagas)
        real(8), intent(in), dimension(:) :: te,ti
        real(8), intent(out), dimension(size(te)) :: gammagas
        real(8), dimension(size(te)) :: the, thi, dxe, dxi, lne,lni, ge, gi, tratio

        the = te*1.68636949d-10  
        thi = ti*9.1842553d-10 ! AC might need to account for non-hydrogen plasma eventually

        dxe = 2.5*the
        dxi = 2.5*thi
        where(dxe.le.1d-8) ! AC hacky log1p function -- import better?
           lne = dxe*(1d0-0.5*dxe + dxe*dxe/3.- 0.25*dxe*dxe*dxe)
        elsewhere
           lne = log(1d0+dxe)
        end where

        where(dxi.le.1d-8) ! AC hacky log1p function -- import better?
           lni = dxi*(1d0-0.5*dxi + dxi*dxi/3. - 0.25*dxi*dxi*dxi)
        elsewhere
           lni = log(1d0+dxi)
        end where

        ge = 1 + the/(3*the - 0.6*lne)
        gi = 1 + thi/(3*thi - 0.6*lni)

        tratio = ti/te
        gammagas = 1. + (1. + tratio) / (1./(ge-1.) + tratio/(gi-1.)) 
    end subroutine calc_gammagas

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Coordinate transformations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! transform BL (r, theta) coordinates to modified Kerr-Schild mks3 (x1,x2)
    subroutine transformbl2mks3(r,th, x1, x2,r0, h,a,b,p)

        real, intent(in), dimension(:) :: r,th
        real, intent(out), dimension(size(th)) :: x1,x2
        real, intent(in) :: r0, h,a,b,p

        x1 = log(dble(r)-r0)
        x2 = 0.5*(1 + ((r**p)/(h*pi)) * ( atan(tan(0.5*h*pi)*(1-2*th/pi)) / ((b-a)*(2**p)+(a-0.5)*(r**p)) ))
    end subroutine transformbl2mks3

    ! transform modified Kerr-Schild mks3 (x1,x2) to Boyer-Lindquist (r,theta)
    subroutine transformmks32bl(x1, x2,r,th,r0, h,a,b,p)
        real, intent(in), dimension(:) :: x1, x2
        !real, intent(in), dimension(:) :: x2,r
        real, intent(out), dimension(size(x2)) :: r, th
        real, intent(in) :: r0, h,a,b,p
        
        r = exp(dble(x1)) + r0
        th = 0.5*pi*(1 + tan(h*pi*(-0.5+x2+(1-2*x2)*(a+(2**p)*(b-a)/(r**p))))/tan(0.5*h*pi))
    end subroutine transformmks32bl

    ! transform BL (r,theta) coordinates to modified Kerr-Schild mks2 (x1,x2)
    subroutine transformbl2mks2(r,th, x1, x2,r0, h)

        real, intent(in), dimension(:) :: r,th
        real, intent(out), dimension(size(th)) :: x1,x2
        real, intent(in) :: r0, h

        x1 = log(dble(r)-r0)
        x2 = 0.5 + atan(tan(0.5*pi*h)*(2*th/pi-1))/(h*pi)
    end subroutine transformbl2mks2

    ! transform modified Kerr-Schild mks2 (x1,x2) to Boyer-Lindquist (r,theta) 
    subroutine transformmks22bl(x1, x2,r,th,r0, h)
        real, intent(in), dimension(:) :: x1, x2
        !real, intent(in), dimension(:) :: x2,r
        real, intent(out), dimension(size(x2)) :: r, th
        real, intent(in) :: r0, h
        
        r = exp(dble(x1)) + r0
        th = 0.5*pi*(1 + tan(h*pi*(x2-0.5))/tan(0.5*h*pi))
    end subroutine transformmks22bl

    ! transform Boyer-Lindquist (r,theta) to jet coordinates (x1,x2) 
    ! this requires numerical inversion
    ! for simplicity use the global variables for parameters
    ! ignore SLOW cylindrification for now -- TODO add it back? 
    subroutine transformbl2jet(r,th,x1,x2)
        real, intent(in), dimension(:) :: r,th
        real, intent(out), dimension(size(th)) :: x1,x2
        !real, intent(in) :: r0,rmin, rmax, rbrk,fjet,fdisk,runi,rcj,rcd,rdj,rdd,al1,al2

        ! transform x1 --> r
        call bl2jet_x1(r, x1) 

        ! transform x2 --> theta
        call bl2jet_x2(r, th, x2) 

    end subroutine transformbl2jet

    ! transform jet coordinates (x1,x2) to Boyer-Lindquist (r,theta)
    ! for simplicity use the global variables for parameters
    ! ignore SLOW cylindrification for now -- TODO add it back? 
    subroutine transformjet2bl(x1, x2, r,th) !,r0, rmin, rmax, rbrk,fjet,fdisk,runi,rcj,rcd,rdj,rdd,al1,al2)
        real, intent(in), dimension(:) :: x1, x2
        real, intent(out), dimension(size(x2)) :: r, th
        !real, intent(in) :: r0,rmin, rmax, rbrk,fjet,fdisk,runi,rcj,rcd,rdj,rdd,al1,al2

        ! transform x1 --> r
        call jet2bl_r(x1, r) 

        ! transform x2 --> theta
        call jet2bl_th(r, x2, th) 
    end subroutine transformjet2bl

    ! jetcoords x1 --> r
    subroutine jet2bl_r(x1, r) 
        real, intent(in), dimension(:) :: x1
        real, intent(out), dimension(size(x1)) :: r
        real, dimension(size(x1)) :: x1scl 

        !AC TODO put in correct hyperexponential behavior for large-radius
        if(rbrk.lt.rmax) then
            write(6,*) 'current jetcoords code only supports grids with no hyperexponential (rbrk>=rmax)!'
            call exit(-1)
        end if
        x1scl = x1min + x1*(x1max-x1min) ! scale from (0,1) range to (x1min,x1max) range
        r = exp(dble(x1scl)) + r0 
    end subroutine jet2bl_r

    ! jetcoords r --> x1
    subroutine bl2jet_x1(r, x1) 
        real, intent(in), dimension(:) :: r
        real, intent(out), dimension(size(r)) :: x1
        real, dimension(size(x1)) :: x1scl 

        !AC TODO put in correct hyperexponential behavior for large-radius
        if(rbrk.lt.rmax) then
            write(6,*) 'current jetcoords code only supports grids with no hyperexponential (rbrk>=rmax)!'
            call exit(-1)
        end if
        x1scl = log(dble(r)-r0)
        x1 = (x1scl-x1min)/(x1max-x1min) ! scale from (x1min, x1max) to (0,1)

    end subroutine bl2jet_x1

    ! jetcoords x2 --> theta, given r
    subroutine jet2bl_th(r, x2, th) 
        real, intent(in), dimension(:) :: r, x2
        real, intent(out), dimension(size(x2)) :: th
        real, dimension(size(x2)) :: thdisk, thjet
        real(kind=8), dimension(size(x2)) :: wfrac
        real, dimension(size(x2)) :: r1, r2, xxjet, yyjet, xxdisk, yydisk, yyy, delta
        real, dimension(size(th)) :: mone, one
        mone(:)=-1d0; one(:)=1d0

        ! calculate theta for a jet-only grid
        r1 = minn(r, rdj*one, 0.5*rdj*one)/runi
        r2 = minn(r/(r1*runi), rcj*one/rdj, 0.5*rcj*one/rdj)
        xxjet = r1**al1
        yyjet = r2**al2
        yyy = yyjet*tan(0.5*pi*x2)
        thjet = 0.5*pi + atan2(yyy,xxjet)

        ! calculate theta for a disk-only grid
        r1 = minn(r, rdd*one, 0.5*rdd*one)/runi
        r2 = minn(r/(r1*runi), rcd*one/rdd, 0.5*rcd*one/rdd)
        xxdisk = r1**al1
        yydisk = r2**al2
        yyy = yydisk*tan(0.5*pi*x2)
        thdisk = 0.5*pi + atan2(yyy,xxdisk)

        ! calculate the fraction of the grid in the disk vs the jet
        delta = 2*(abs(x2)-fdisk)/(1-fjet-fdisk)-1
        wfrac = thetasmooth(dble(delta))

        ! combine the two values of theta
        th = wfrac*thjet + (1-wfrac)*thdisk

    end subroutine jet2bl_th 

    ! ks theta --> jetcoords x2, given r
    ! numerically inverts using zbrent
    subroutine bl2jet_x2(r, th, x2) 
        real, intent(in), dimension(:) :: r, th
        real, intent(out), dimension(size(th)) :: x2
        real, dimension(size(th)) :: thdisk, thjet, wfrac
        real, dimension(size(th)) :: r1, r2, xxjet, yyjet, xxdisk, yydisk, delta
        real, dimension(size(th)) :: mone, one
        real(kind=8), dimension(size(th)) :: dmone, done, x22
        real(kind=8), dimension(size(th), 5) :: searchargs
        mone(:)=-1d0; one(:)=1d0
        dmone(:)=-1d0; done(:)=1d0

        ! we can precompute several constant factors
        ! that depend only on r and globals
        r1 = minn(r, rdj*one, 0.5*rdj*one)/runi
        r2 = minn(r/(r1*runi), rcj*one/rdj, 0.5*rcj*one/rdj)
        xxjet = r1**al1
        yyjet = r2**al2

        r1 = minn(r, rdd*one, 0.5*rdd*one)/runi
        r2 = minn(r/(r1*runi), rcd*one/rdd, 0.5*rcd*one/rdd)
        xxdisk = r1**al1
        yydisk = r2**al2

        searchargs(:,1) = dble(th)
        searchargs(:,2) = dble(xxjet)
        searchargs(:,3) = dble(yyjet)
        searchargs(:,4) = dble(xxdisk)
        searchargs(:,5) = dble(yydisk)

        ! now find the solution numerically with zbrent
        x22 = zbrent(findx2diff,dmone,done,searchargs,1d-6)
        x2 = real(x22)
    end subroutine bl2jet_x2

    ! function for zbrent to find x2(theta)
    function findx2diff(x2, searchargs) result (diff)
        real(kind=8), intent(in), dimension(:) :: x2
        real(kind=8), intent(in), dimension(:,:) :: searchargs
        real(kind=8), dimension(size(x2)) :: th, theta, diff
        real(kind=8), dimension(size(x2)) :: xxjet, yyjet, xxdisk, yydisk, yyy, thjet, thdisk
        real(kind=8), dimension(size(x2)) :: wfrac, delta

        th = searchargs(:,1)
        xxjet = searchargs(:,2)
        yyjet = searchargs(:,3)
        xxdisk = searchargs(:,4)
        yydisk = searchargs(:,5)

        yyy = yyjet*tan(0.5*pi*x2)
        thjet = 0.5*pi + atan2(yyy,xxjet)
        yyy = yydisk*tan(0.5*pi*x2)
        thdisk = 0.5*pi + atan2(yyy,xxdisk)

        delta = 2*(abs(x2)-fdisk)/(1-fjet-fdisk)-1.
        wfrac = thetasmooth(delta)
        theta = wfrac*thjet + (1-wfrac)*thdisk

        diff = th - theta
    end function findx2diff

    ! smooth minimum
    function minn(a, b, df) result(res)
        real, intent(in), dimension(:) :: a, b, df
        real, dimension(size(a)) :: res
        real, dimension(size(a)) :: x,psi

        x = (b-a)/df
        where(x.lt.-1.)
           psi=0.
        elsewhere(x.gt.1.)
           psi=x
        elsewhere
           psi = ((-35*cos(0.5*pi*x)-(5./6.)*cos(1.5*pi*x) + 0.1*cos(2.5*pi*x))/(32.*pi)) + 0.5*(x+1.)
        end where
        res = b - psi*df

    end function minn

    ! smooth step
    function thetasmooth(x) result(res)
        real(kind=8), intent(in), dimension(:) :: x
        real(kind=8), dimension(size(x)) :: res
        where(x.lt.-1.)
           res=0.
        elsewhere(x.gt.1.)
           res=1.
        elsewhere
           res=0.5 + (70.*sin(0.5*pi*x) + 5*sin(1.5*pi*x) - sin(2.5*pi*x))/128.
        end where
    end function thetasmooth

    ! Locate a value in a sorted array
    function locate(xx,x) result(idx)
        implicit none
        real(8), dimension(:), intent(in) :: xx ! input arry
        real(8), intent(in) :: x ! point to locate
        integer :: n,jl,jm,ju,idx
        logical :: ascnd
        n=size(xx)
        ascnd = (xx(n) >= xx(1)) !is the array ascending or descending? 

        jl=0
        ju=n+1
        do
           if (ju-jl <= 1) exit
           jm=(ju+jl)/2
           if (ascnd.eqv.(x>=xx(jm))) then
              jl=jm
           else
              ju=jm
           end if
        end do

        if (x == xx(1)) then
           idx = 1
        else if (x == xx(n)) then
           idx = n-1
        else if(ascnd.and. (x < xx(1))) then ! over left edge
           idx = 0
        else if(.not.ascnd.and. (x > xx(1))) then ! over left edge
           idx = 0
        else if(ascnd.and. (x > xx(n))) then ! over right edge
           idx = n
        else if(.not.ascnd.and. (x < xx(n))) then ! over right edge
           idx = n
        else
           idx = jl
        end if

    end function locate

    ! 1d interpolation to single value
    function interp1d_single(x,y,xval) result(yval)
        implicit none
        integer :: n ! the size of the array
        real(8),dimension(:),intent(in) :: x,y ! the x and y arrays -- don't have to be equally spaced
        real(8),intent(in) :: xval ! the value at which to interpolate y
        real(8) :: yval ! the result
        integer :: ipos ! position of x value in x array
        real(8) :: frac
        
        n = size(x)
        ipos = locate(x,xval)

        if(ipos<n.and.ipos>0) then
           frac = (xval - x(ipos)) / (x(ipos+1) - x(ipos))
           yval = y(ipos) + frac * (y(ipos+1) - y(ipos))
        else if(ipos == n) then
           yval = y(n)
        else if(ipos == 0) then
           yval = y(1)
        else
           write(6,*) 'Error in interp1d_single!'
           stop
        end if

    end function interp1d_single

    function interp1d_arr(x,y,xvals) result(yvals) ! interpolate multiple xvals on same sampling
        implicit none
        integer :: nvals ! the size of the array
        real(8),dimension(:),intent(in) :: x,y ! the x and y arrays -- don't have to be equally spaced
        real(8),dimension(:),intent(in) :: xvals ! the value at which to interpolate y
        real(8),dimension(size(xvals)) :: yvals ! the result
        integer :: i 
        nvals = size(xvals)
        do i=1,nvals
            yvals(i) = interp1d_single(x,y,xvals(i))
        end do

    end function interp1d_arr
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Interpolates KORAL data to input coordinates
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine koralh5_vals(x0,a,rho,p,te,ti,b,u,bmag,Be,type,electrons)

        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        integer, intent(in) :: type
        integer, intent(out) :: electrons
        real(kind=8), dimension(size(x0)) :: done,nfac,pfac,bfac
        real, dimension(size(x0)) :: x3, x2,x1,zm,theta,fone, &
         vpl0,vrl0,vtl0,rd,td,pd,rttd,zr,dzero, &
         vr0,vth0,vph0,bth,dummy,tt,ttd,zt,zphi,zpp

        real, dimension(size(x0),2**(ndim+1)) :: ppi,Bei,rhoi,vrli,vtli, &
             vpli,bri,bthi,bphi,b0i,u0i,tei,tii

        real, dimension(size(x0)) :: dth,thu,thl,thone,minph, rdummy, x1trans
        integer, dimension(size(x0)*2*(2**ndim)) :: indx
        integer, dimension(size(x0)) :: lx1,lx2,lx3,ux1,ux2,ux3,x1l,x1u,x2l,x2u, x3l,x3u, &
         one,tindx,lx2l,lx2u,ux2l,ux2u,lx3l,lx3u,ux3l,ux3u,low,high
        real :: x2mintrust, x2maxtrust,thmintrust,thmaxtrust
        integer :: npts,i,j,maxtindx,ie
        real, dimension(size(x0)) :: b1tmp,b2tmp,b3tmp,b4tmp,u1tmp,ntmp
        real, dimension(size(x0)), intent(out) :: rho,p,bmag,Be,te,ti
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b
        integer :: jj
        real, dimension(nx2+2) :: th_arr_l, th_arr_u
        real(8) :: x2_a, x2_b, th_a, th_b

        ! for JETCOORDS, this determines if we use explicit transforms or interpolation
        ! inverse JETCOORDS transform is super expensive, and we do NOT use cylindrification
        logical, parameter :: jetcoords_interp=.True.

        one=1; dzero=0d0; done=1d0; fone=1.
        npts=size(x0)

        ! tell calling program if we have electrons or not
        electrons = has_electrons
   
        ! Transform from BL to Kerr-Schild to code x1,x2,x3 coordinates
        zt=x0%data(1)
        zr=x0%data(2)
        theta=x0%data(3)
        zpp=x0%data(4)

        ! Transform time from BL -->KS
        tt = bl2ks(dble(zr),dble(zt),dble(a),1) - bl2ks(dble(zr(1)),0d0,dble(a),1)
        tt = -tt + 1d-8

        ! Transform phi-coordinate from BL -->KS
        zphi=bl2ks(dble(zr), dble(zpp), dble(a))
        zphi=mod(zphi,(2.*pi))
        if(any(isnan(zphi))) then
           write(6,*) 'nan in zphi!'
        end if
        where(zphi.lt.0.)
            zphi = zphi + 2.*pi
        end where
        where(zphi.gt.pi) ! koral range is (-pi,pi) 
            zphi=zphi - 2.*pi
        end where

        ! Get code coordinate (x1,x2,x3) from KS coordinate (zr, theta, zphi)
        if(runmet.eq.'MKS2') then
            call transformbl2mks2(zr,theta, x1, x2,r0, h0)
        else if(runmet.eq.'MKS3') then 
            call transformbl2mks3(zr,theta, x1, x2,r0, h0,my1,my2,pp)
        else if(runmet.eq.'JETCOORDS') then
            if(jetcoords_interp) then
            ! for JETCOORDS, we will interpolate x2 instead of directly transforming)
                ! first, get x1 here
                call bl2jet_x1(zr, x1) ! TODO -- might want to interpolate here too if we use hyperexp
                ! we get x2 below after we find r indices
            else
                call transformbl2jet(zr,theta, x1, x2) 
            endif

        else
            write(6,*) 'cannot determine run metric!'
            call exit(-1)
        endif
        x3 = zphi ! we always leave the phi-coordinate unchanged

        ! Find the cell numbers that bracket the current points in each dimension
        lx1=floor((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        lx1=merge(lx1,one,lx1.ge.1) ! this enforces lx1>=1
        lx1=merge(lx1,nx1-1,lx1.le.nx1-1) ! this enforces lx1<=nx1-1
        ux1=lx1+1 

        ! need to get x2 by interpolation
        if((runmet.eq.'JETCOORDS').and.jetcoords_interp) then 
            th_arr_l(1) = 0; th_arr_u(1) = 0; ! extend range of theta to physical bounds
            th_arr_l(nx2+2) = pi; th_arr_u(nx2+2) = pi;

            do jj=1,npts
               ! get x2 by interpolating grid values at both lower and upper radius,
               ! then average
               th_arr_l(2:nx2+1) = r_th_arr(lx1(jj),:)
               th_arr_u(2:nx2+1) = r_th_arr(ux1(jj),:)
               x2_a = interp1d_single(dble(th_arr_l),dble(uniqx2_ext), dble(theta(jj)))
               x2_b = interp1d_single(dble(th_arr_u),dble(uniqx2_ext), dble(theta(jj)))
               x2(jj) = real(0.5*(x2_a + x2_b)) ! average the upper and lower x2 points (OK?)
            end do
        endif

        lx2=floor((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
        ux2=lx2+1
        lx2=merge(merge(lx2,one,lx2.ge.1),nx2-1,lx2.le.(nx2-1)) ! this enforces 1<=lx2<=nx2-1
        ux2=merge(ux2,one,ux2.ge.1)   ! AC TODO not sure why/if we want this sequence for lx2,ux2 -- poles?
        ux2=merge(ux2,nx2,ux2.le.nx2) ! AC took this from fluid_model_harm3d.f90

        lx3=floor((x3-uniqx3(1))/(uniqx3(nx3)-uniqx3(1))*(nx3-1))+1
        ux3=lx3+1

        if(nx3.eq.1) then ! 2D
            lx3 = 1
            ux3 = 1
            pd(:) = 1.

        else
            ! enforce periodic in phi
            minph=uniqph(lx3)
            where(ux3.gt.nx3) 
                !lx3 = nx3 ??
                !minph = uniqph(lx3)
                ux3=1
            end where
            where(lx3.lt.1) 
                lx3=nx3
                minph=uniqph(lx3)-2.*pi ! this should never be used since we push x3 range to (-pi,pi)
            end where
            !write(6,*) 'indices'
            ! phi-distance
            pd=((zphi-minph)/(uniqph(2)-uniqph(1))) ! grid is uniform in phi so we can use constant deltaph
        endif

        ! r-distance
        rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))

        ! theta-distance 
        ! need to transform x2 points dependent on r
        ! then deal with poles 
        x1trans = 0.5*(uniqx1(lx1) + uniqx1(ux1)) ! do coordinate transformations between lx1,ux1
        if(runmet.eq.'MKS2') then
            call transformmks22bl(x1trans, uniqx2(lx2), rdummy,thl,r0,h0) 
            call transformmks22bl(x1trans, uniqx2(ux2), rdummy,thu,r0,h0)
            call transformmks22bl(x1trans, uniqx2(one), rdummy,thone,r0,h0)
        else if(runmet.eq.'MKS3') then 
            call transformmks32bl(x1trans, uniqx2(lx2), rdummy,thl,r0,h0,my1,my2,pp) 
            call transformmks32bl(x1trans, uniqx2(ux2), rdummy,thu,r0,h0,my1,my2,pp)
            call transformmks32bl(x1trans, uniqx2(one), rdummy,thone,r0,h0,my1,my2,pp)
        else if(runmet.eq.'JETCOORDS') then
            if(jetcoords_interp) then
                ! get theta by interpolating grid values at both lower and upper radius,
                th_arr_l(1) = 0; th_arr_u(1) = 0; ! extend range of theta to physical bounds
                th_arr_l(nx2+2) = pi; th_arr_u(nx2+2) = pi;
                do jj=1,npts
                   th_arr_l(2:nx2+1) = r_th_arr(lx1(jj),:)
                   th_arr_u(2:nx2+1) = r_th_arr(ux1(jj),:)

                   th_a = interp1d_single(dble(uniqx2_ext), dble(th_arr_l),dble(uniqx2(lx2(jj))))
                   th_b = interp1d_single(dble(uniqx2_ext), dble(th_arr_u),dble(uniqx2(lx2(jj))))
                   thl(jj) = real(0.5*(th_a + th_b)) ! average the upper and lower theta points (OK?)

                   th_a = interp1d_single(dble(uniqx2_ext), dble(th_arr_l),dble(uniqx2(ux2(jj))))
                   th_b = interp1d_single(dble(uniqx2_ext), dble(th_arr_u),dble(uniqx2(ux2(jj))))
                   thu(jj) = real(0.5*(th_a + th_b)) ! average the upper and lower theta points (OK?)

                   th_a = interp1d_single(dble(uniqx2_ext), dble(th_arr_l),dble(uniqx2(one(jj))))
                   th_b = interp1d_single(dble(uniqx2_ext), dble(th_arr_u),dble(uniqx2(one(jj))))
                   thone(jj) = real(0.5*(th_a + th_b)) ! average the upper and lower theta points (OK?)
                end do
            else
                ! use direct transformation -- no cylindrification
                call transformjet2bl(x1trans, uniqx2(lx2), rdummy,thl) 
                call transformjet2bl(x1trans, uniqx2(ux2), rdummy,thu)
                call transformjet2bl(x1trans, uniqx2(one), rdummy,thone)
            endif
        else
            write(6,*) 'cannot determine run metric!'
            call exit(-1)
        endif
        
        where(ux2.ne.lx2)
           td=abs((theta-thl)/(thu-thl))
        elsewhere
           td=abs((theta-thl)/(2*thone))
        end where

        ! for nearest neighbor, use this instead of above
        !rd(:)=1.; td(:)=1.; pd(:)=1.

        ! When the geodesic is inside the innermost zone center located outside the horizon, 
        ! zero out the primitives
        where((uniqr(lx1).le.(1.+sqrt(1.-a**2.))).or.(lx1.eq.1))
           rd=1. ! use nearest neighbor rather than interpolate
           pfac=1e-6
           nfac=1e-6
           bfac=1e-6
        elsewhere
           pfac=1.
           nfac=1.
           bfac=1.
        end where

        ! Indices for nearest neighbors in 3D
        ! x2 is fastest changing index
        x2l=lx2-1
        x2u=ux2-1
        x1l=(lx1-1)*nx2 
        x1u=(ux1-1)*nx2
        x3l=(lx3-1)*nx2*nx1
        x3u=(ux3-1)*nx2*nx1

        ! find timestep indices:
        maxtindx=nt-2
        where(floor(tt/tstep).le.maxtindx)
           tindx=floor(tt/tstep)
           ttd=(tt-tindx*tstep)/tstep
        elsewhere
           tindx=maxtindx+1
           ttd=0.
        end where

        ! Create index array across times and 3D nearest neighbors (first half at first time, second half at second time):
        ! Need to make sure there aren't integer problems when these indices get very large
        indx=(/(/x1l+x2l+x3l+1+n*tindx/),(/x1l+x2u+x3l+1+n*tindx/),(/x1u+x2l+x3l+1+n*tindx/), &
               (/x1u+x2u+x3l+1+n*tindx/),(/x1l+x2l+x3u+1+n*tindx/), &
               (/x1l+x2u+x3u+1+n*tindx/),(/x1u+x2l+x3u+1+n*tindx/),(/x1u+x2u+x3u+1+n*tindx/), &
               (/x1l+x2l+x3l+1+n*(tindx+1)/),(/x1l+x2u+x3l+1+n*(tindx+1)/),(/x1u+x2l+x3l+1+n*(tindx+1)/), &
               (/x1u+x2u+x3l+1+n*(tindx+1)/),(/x1l+x2l+x3u+1+n*(tindx+1)/),(/x1l+x2u+x3u+1+n*(tindx+1)/), &
               (/x1u+x2l+x3u+1+n*(tindx+1)/),(/x1u+x2u+x3u+1+n*(tindx+1)/)/)

        ! keep indx from going out of bounds when loading one time slice
        where(indx.gt.size(rho_arr))
           indx=indx-n
        end where
     
        ! Get data values at array points
        rhoi=reshape(rho_arr(indx),(/npts,2**(ndim+1)/))
        ppi=reshape(p_arr(indx),(/npts,2**(ndim+1)/))
        Bei=reshape(Be_arr(indx),(/npts,2**(ndim+1)/))

        u0i=reshape(u0_arr(indx),(/npts,2**(ndim+1)/))
        vrli=reshape(vrl_arr(indx),(/npts,2**(ndim+1)/))
        vtli=reshape(vtl_arr(indx),(/npts,2**(ndim+1)/))
        vpli=reshape(vpl_arr(indx),(/npts,2**(ndim+1)/))

        b0i=reshape(b0_arr(indx),(/npts,2**(ndim+1)/))
        bri=reshape(br_arr(indx),(/npts,2**(ndim+1)/))
        bthi=reshape(bth_arr(indx),(/npts,2**(ndim+1)/))
        bphi=reshape(bph_arr(indx),(/npts,2**(ndim+1)/))

        ! interpolate data
        rttd=0.
        rho=interp(rhoi,rttd,pd,rd,td)
        p=interp(ppi,rttd,pd,rd,td)
        Be=interp(Bei,rttd,pd,rd,td)

        u1tmp=interp(u0i,rttd,pd,rd,td)
        vrl0=interp(vrli,rttd,pd,rd,td)
        vtl0=interp(vtli,rttd,pd,rd,td)
        vpl0=interp(vpli,rttd,pd,rd,td)

        b1tmp=interp(b0i,rttd,pd,rd,td)
        b2tmp=interp(bri,rttd,pd,rd,td)
        b3tmp=interp(bthi,rttd,pd,rd,td)
        b4tmp=interp(bphi,rttd,pd,rd,td)

        !mask data outside trusted region
        rho=merge(rho,dzero,x1.gt.uniqx1(1))*nfac
        p=merge(p,fone,x1.gt.uniqx1(1))*pfac 
        Be=merge(Be,dzero,x1.gt.uniqx1(1)) !AC -- does using zero make sense?

        u1tmp=merge(dble(u1tmp), done,x1.gt.uniqx1(1))
        vrl0=merge(vrl0,dzero,x1.gt.uniqx1(1))
        vtl0=merge(vtl0,dzero,x1.gt.uniqx1(1))
        vpl0=merge(vpl0,dzero,x1.gt.uniqx1(1))

        b1tmp=merge(b1tmp,fone,x1.gt.uniqx1(1))
        b2tmp=merge(b2tmp,fone,x1.gt.uniqx1(1))
        b3tmp=merge(b3tmp,fone,x1.gt.uniqx1(1))
        b4tmp=merge(b4tmp,fone,x1.gt.uniqx1(1))

        ! polar angle cuts
        ! AC -- assumes symmetric grid!!
        x2mintrust = uniqx2(minpolecell)
        x2maxtrust = uniqx2(size(uniqx2)-(minpolecell-1))
        thmintrust = 0.
        thmaxtrust = pi

        ! by default cut out 5 degrees away from poles for disk (also cut out high Be in fluid.f90)
        if (type.eq.1) then
           thmintrust = 0.09
           thmaxtrust = pi - 0.09

        ! cut out background entirely -- all theta > pi/2 
        else if(type.eq.2) then 
           thmaxtrust = 0.5*pi

        ! cut out foreground entirely -- all theta > pi/2 !AC -- assumes symmetric grid!!       
        else if(type.eq.3) then
           thmintrust = 0.5*pi
        end if

        ! theta based cuts for jet
        if((type.eq.1).or.(type.eq.2).or.(type.eq.3)) then
            rho=merge(rho,dzero,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
            p=merge(p,fone,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
            Be=merge(Be,dzero,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))

            u1tmp=merge(dble(u1tmp),done,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
            vrl0=merge(vrl0,dzero,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
            vtl0=merge(vtl0,dzero,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
            vpl0=merge(vpl0,dzero,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))

            b1tmp=merge(b1tmp,fone,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
            b2tmp=merge(b2tmp,fone,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
            b3tmp=merge(b3tmp,fone,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
            b4tmp=merge(b4tmp,fone,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
        end if

        ! x2 based cuts for removing polar axis cells
        rho=merge(rho,dzero,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        p=merge(p,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust)) 
        Be=merge(Be,dzero,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))

        u%data(1)=merge(dble(u1tmp),done,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        vrl0=merge(vrl0,dzero,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        vtl0=merge(vtl0,dzero,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        vpl0=merge(vpl0,dzero,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))

        b%data(1)=merge(b1tmp,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        b%data(2)=merge(b2tmp,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        b%data(3)=merge(b3tmp,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        b%data(4)=merge(b4tmp,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))

        ! te,ti only if electons are present
        if(has_electrons.eq.1) then
          tei=reshape(te_arr(indx),(/npts,2**(ndim+1)/))
          te=interp(tei,rttd,pd,rd,td)
          te=merge(te,fone,x1.gt.uniqx1(1))*pfac
          te=merge(te,fone,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
          te=merge(te,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust)) 

          tii=reshape(ti_arr(indx),(/npts,2**(ndim+1)/))
          ti=interp(tii,rttd,pd,rd,td)
          ti=merge(ti,fone,x1.gt.uniqx1(1))*pfac
          ti=merge(ti,fone,(theta.gt.thmintrust).and.(theta.lt.thmaxtrust))
          ti=merge(ti,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust)) 
        endif

        ! Compute magnitude of interpolated b-field and force b^2 > 0 (need to look into numerical issues here):
        call assign_metric(b,transpose(kerr_metric(zr, real(x0%data(3)),a)))
        bmag=b*b 
        bmag=merge(bmag,dzero,bmag.ge.0d0)
        bmag=sqrt(bmag)*bfac

        ! Get 4-velocity and Protect azimuthal velocities at poles
        call lnrf_frame(vrl0,vtl0,vpl0,zr,a,real(x0%data(3)),vr0,vth0,vph0,1)
        u%data(2)=u%data(1)*dble(vr0)
        u%data(3)=u%data(1)*dble(vth0)
        u%data(4)=u%data(1)*dble(vph0)
        call assign_metric(u,transpose(kerr_metric(zr,real(x0%data(3)) ,a)))
    end subroutine koralh5_vals

    ! Read input filenames from .in file
    subroutine read_koralh5_inputs(ifile)
        character(len=500), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=harm) !\harm\ dfile, hfile, nt, indf
        write(6,*) 'read: ', dfile, hfile, nt, indf
        close(unit=8)
    end subroutine read_koralh5_inputs

    ! Read header values from hdf5
    subroutine read_koralh5_data_header()

        integer(HID_T) :: file, dset, dspace, memtype, filetype
        integer :: hdferr ! error code

        integer(HSIZE_T), dimension(3) :: dims = (/1,1,1/) ! size of hdf5 read buffer, dummy
        integer(HSIZE_T), dimension(3) :: maxdims ! size of hdf5 read buffer, dummy
        integer(SIZE_T), parameter :: maxsize = 50 
        integer(SIZE_T) :: strsize
        real(8) :: rdataSd ! read buffer for real scalars
        integer :: rdataSi ! read buffer for integers
        character(len=maxsize), target :: rdataS ! read buffer for string
        type(C_PTR) :: f_ptr

        ! open file
        write(6,*) '---------------------------------'
        write(6,*) 'read koral HDF5 header from ', trim(hfile)
        write(6,*) '---------------------------------'
        call h5open_f(hdferr) 
        call h5fopen_f(hfile, H5F_ACC_RDONLY_F, file, hdferr)

        ! dimensions
        call h5dopen_f(file, '/header/n1', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_INTEGER, rdataSi, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        nx1 = rdataSi

        call h5dopen_f(file, '/header/n2', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_INTEGER, rdataSi, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        nx2 = rdataSi

        call h5dopen_f(file, '/header/n3', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_INTEGER, rdataSi, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        nx3 = rdataSi

        ! spin
        call h5dopen_f(file, '/header/bhspin', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        asim = rdataSd

        ! adiabatic index
        call h5dopen_f(file, '/header/gam', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        gam = rdataSd

        ! electrons or no
        call h5dopen_f(file, '/header/has_electrons', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_INTEGER, rdataSi, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        has_electrons = rdataSi

        ! grid parameters 
        call h5dopen_f(file, '/header/geom/startx1', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        startx1 = rdataSd

        call h5dopen_f(file, '/header/geom/startx2', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        startx2 = rdataSd

        call h5dopen_f(file, '/header/geom/startx3', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        startx3 = rdataSd

        call h5dopen_f(file, '/header/geom/dx1', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        dx1 = rdataSd

        call h5dopen_f(file, '/header/geom/dx2', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        dx2 = rdataSd

        call h5dopen_f(file, '/header/geom/dx3', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
        call h5dclose_f(dset, hdferr)
        dx3 = rdataSd

        ! run metric ! AC is there an easier way to get this string? 
        call h5dopen_f(file, '/header/metric_run', dset, hdferr)
        call H5Dget_type_f(dset, filetype, hdferr)
        call H5Tget_size_f(filetype, strsize, hdferr)
        call H5Dget_space_f(dset, dspace, hdferr)
        call H5Sget_simple_extent_dims_f(dspace, dims, maxdims, hdferr)
        call H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
        call H5Tset_size_f(memtype, maxsize, hdferr)
        f_ptr = C_LOC(rdataS)!C_LOC(rdataS(1:1))
        call H5Dread_f(dset, memtype, f_ptr, hdferr, dspace)
        runmet = rdataS
        write(6,*) 'metric_run ', runmet

        call H5Dclose_f(dset, hdferr)
        call H5Sclose_f(dspace, hdferr)
        call H5Tclose_f(memtype, hdferr)

        ! out metric ! AC is there an easier way to get this string? 
        call h5dopen_f(file, '/header/metric_out', dset, hdferr)
        call H5Dget_type_f(dset, filetype, hdferr)
        call H5Tget_size_f(filetype, strsize, hdferr)
        call H5Dget_space_f(dset, dspace, hdferr)
        call H5Sget_simple_extent_dims_f(dspace, dims, maxdims, hdferr)
        call H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
        call H5Tset_size_f(memtype, maxsize, hdferr)
        f_ptr = C_LOC(rdataS)!C_LOC(rdataS(1:1))
        call H5Dread_f(dset, memtype, f_ptr, hdferr, dspace)
        outmet = rdataS
        write(6,*) 'metric_out ', outmet

        call H5Dclose_f(dset, hdferr)
        call H5Sclose_f(dspace, hdferr)
        call H5Tclose_f(memtype, hdferr)

        if(.not.(outmet.eq.'KS' .or. outmet.eq.'BL')) then
           write(6,*) 'koral hdf5 currently only works with KS or BL output metric!' 
           call exit(-1)
        endif

        write(6,*) 'read hdf5 header dims: ', nx1, nx2, nx3
        write(6,*) 'read hdf5 header grid start: ', startx1, startx2, startx3
        write(6,*) 'read hdf5 header grid delta: ', dx1, dx2, dx3
        write(6,*) 'read hdf5 header spin: ', asim
        write(6,*) 'read hdf5 header gamma: ', gam
        write(6,*) 'read hdf5 header has_electrons: ', has_electrons

        ! coordinate parameters
        if(runmet.eq.'JETCOORDS') then

            call h5dopen_f(file, '/header/geom/jetcoords/mksr0', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            r0 = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/rmin', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            rmin = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/rmax', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            rmax = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/rbrk', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            rbrk = rdataSd

            if(rbrk.lt.rmax) then
                write(6,*) 'current jetcoords code only supports grids with no hyperexponential (rbrk>=rmax)!'
                call exit(-1)
            end if
      
            x1min = log(dble(rmin)-r0) !x1 coordinate minimum and maximum values
            x1max = log(dble(rmax)-r0)

            call h5dopen_f(file, '/header/geom/jetcoords/fjet', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            fjet = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/fdisk', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            fdisk = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/runi', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            runi = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/rcoll_jet', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            rcj = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/rcoll_disk', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            rcd = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/rdecoll_jet', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            rdj = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/rdecoll_disk', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            rdd = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/alpha1', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            al1 = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/alpha2', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            al2 = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/rcyl', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            rcyl = rdataSd

            call h5dopen_f(file, '/header/geom/jetcoords/ncyl', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            ncyl = rdataSd

            write(6,*) 'read jetcoords header vals 1: ', r0,rmin,rmax,rbrk
            write(6,*) 'read jetcoords header vals 2: ', fjet,fdisk,al1, al2
            write(6,*) 'read jetcoords header vals 3: ', runi,rcj,rcd,rdj,rdd
            write(6,*) 'read jetcoords header vals 4: ', rcyl,ncyl
        else if(runmet.eq.'MKS3') then
            call h5dopen_f(file, '/header/geom/mks3/mksr0', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            r0 = rdataSd

            call h5dopen_f(file, '/header/geom/mks3/mksh0', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            h0 = rdataSd

            call h5dopen_f(file, '/header/geom/mks3/mksmy1', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            my1 = rdataSd

            call h5dopen_f(file, '/header/geom/mks3/mksmy2', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            my2 = rdataSd

            call h5dopen_f(file, '/header/geom/mks3/mksp0', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            pp = rdataSd

            write(6,*) 'read mks3 header vals: ', r0,h0, my1, my2, pp
        else if(runmet.eq.'MKS2') then
            call h5dopen_f(file, '/header/geom/mks2/mksr0', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            r0 = rdataSd

            call h5dopen_f(file, '/header/geom/mks2/mksh0', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataSd, dims, hdferr)
            call h5dclose_f(dset, hdferr)
            h0 = rdataSd

            write(6,*) 'read mks2 header vals: ', r0,h0
        else
           write(6,*) 'koral hdf5 currently only works with JETCOORDS or MKS3 or MKS2 run metric!' 
           call exit(-1)
        endif

        ! close file
        call h5fclose_f(file , hdferr)
        call h5close_f(hdferr) 
    end subroutine read_koralh5_data_header

    ! Read data values from hdf5
    subroutine read_koralh5_data_file(data_file,tcur,rho,p,te,ti,u,b,Be)
        character(len=500), intent(in) :: data_file
        real(8), intent(out) :: tcur
        real(8), dimension(:), allocatable, intent(out) :: p,rho,Be,te,ti
        real(8), dimension(:), allocatable :: t00, rhoucon, bsq,gammagas
        type(four_vector), dimension(:), allocatable, intent(out) :: u,b
        type(four_vector), dimension(:), allocatable :: ucov,bcov
        integer :: i, n
        real(8) :: minrho

        integer(HID_T) :: file, dset, dspace 
        integer :: hdferr ! error code
        integer(HSIZE_T), dimension(3) :: dims ! size of hdf5 read buffer
        real(8), dimension(:,:,:), allocatable :: rdata ! read buffer for real arrays
        real(8) :: rdataS ! read buffer for real scalars

        ! for testing (x2, theta) interpolation
!        real, dimension(:), allocatable :: th_arr_l
!        integer :: testidx
!        real :: testth
!        real, dimension(1) :: tone, x2val

        ! allocate read buffer and data arrays
        dims = (/nx3,nx2,nx1/)  ! hdf5 array size -- NEED BACKWARDS ORDERING!!?? (z,y,x)
                                ! AC TODO -- can we make this more intiuitive? 
        allocate(rdata(nx3,nx2,nx1)) 
        n = nx1*nx2*nx3

        allocate(rho(n)); allocate(p(n)); allocate(Be(n))
        allocate(u(n)); allocate(b(n))
        allocate(ucov(n)); allocate(bcov(n))
        allocate(bsq(n)); allocate(rhoucon(n)); allocate(t00(n))
        if(has_electrons.eq.1) then
          allocate(te(n)); allocate(ti(n));allocate(gammagas(n))
        endif 
        ! open file

        call h5open_f(hdferr) 
        call h5fopen_f(data_file, H5F_ACC_RDONLY_F, file, hdferr)

        ! current time
        call h5dopen_f(file, '/t', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdataS, dims, hdferr)!, H5S_ALL_F, dspace, H5P_DEFAULT_F)
        tcur = rdataS

        write(6,*) '---------------------------------'
        write(6,*) 'read koral HDF5 data from ', trim(data_file), ' tcur: ', tcur
        write(6,*) '---------------------------------'

        ! create dataspace for reading arrays
        call h5screate_simple_f(3, dims, dspace, hdferr, dims) 

        ! read code grid -- we no longer need this since we can use left edge information!
!        ! x1 coordinate
!        call h5dopen_f(file, '/grid_run/x1', dset, hdferr, H5S_ALL_F)
!        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
!        call h5dclose_f(dset , hdferr)
!        do i=1, nx3  ! x2 (x3) is the fastest (slowest) changing index in koral convention 
!           x1_arr(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
!        end do

!        ! x2 coordinate
!        call h5dopen_f(file, '/grid_run/x2', dset, hdferr, H5S_ALL_F)
!        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
!        call h5dclose_f(dset , hdferr)
!        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
!           x2_arr(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
!        end do

!        ! x3 coordinate
!        call h5dopen_f(file, '/grid_run/x3', dset, hdferr, H5S_ALL_F)
!        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
!        call h5dclose_f(dset , hdferr)
!        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
!           x3_arr(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
!        end do
!   
!        uniqx2=x2_arr(1:nx2) ! fastest changing
!        uniqx1=x1_arr(1:nx2*(nx1-1)+1:nx2)
!        uniqx3=x3_arr(1:nx2*(nx1-1)*nx3+1:nx2*nx1)  ! slowest changing

        ! calculate code grid from left edge
        do i=1, nx1
          uniqx1(i) = startx1 + (float(i)-0.5)*dx1
        end do
        do i=1, nx2
          uniqx2(i) = startx2 + (float(i)-0.5)*dx2
        end do
        do i=1, nx3
          uniqx3(i) = startx3 + (float(i)-0.5)*dx3
        end do

        ! calculate unique r, ph points
        ! there are not uniq th points because th(x2) depends on r!
        if(runmet.eq.'JETCOORDS') then
           uniqr = exp(x1min + dble(uniqx1)*(x1max-x1min)) + r0 ! AC TODO add hyperexp case 
           uniqx2_ext(2:nx2+1) = uniqx2 ! this is extended to include end points (so is no longer equally spaced!) 
           uniqx2_ext(1) = -1.  ! these edge values hold specifically for JETCOORDS (only place this is used) 
           uniqx2_ext(nx2+2) = 1. 
        else
           uniqr = exp(dble(uniqx1)) + r0
        endif
        uniqph = uniqx3 ! the phi coordinate is always unchanged

        ! read output grid
        ! r coordinate
        call h5dopen_f(file, '/grid_out/r', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           r_arr(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
        end do

        ! th coordinate
        call h5dopen_f(file, '/grid_out/th', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           th_arr(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
        end do

        ! store grid of th for use in interpolation: 
        ! we ASSUME (r,th) is independent of phi (axisymmetric) 
        r_th_arr = transpose(rdata(1,:,:))

        ! test the interpolation
!        testidx = 100
!        testth = 3.
!        allocate(th_arr_l(nx2+2))
!        th_arr_l(1) = 0
!        th_arr_l(2:nx2+1) = r_th_arr(testidx,:)
!        th_arr_l(nx2+2) = pi

!        write(6,*) uniqx2_ext
!        write(6,*)
!        write(6,*) th_arr_l

!        tone = 1.
!        call bl2jet_x2(uniqr(testidx)*tone, testth*tone, x2val)
!        write(6,*) , x2val,interp1d_single(dble(th_arr_l), dble(uniqx2_ext), dble(testth))
!        deallocate(th_arr_l)
!        call exit(-1)


        ! phi coordinate
        call h5dopen_f(file, '/grid_out/ph', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           ph_arr(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
        end do

        ! code quantities
        ! rho
        call h5dopen_f(file, '/quants/rho', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           rho(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
        end do

        ! pgas
        call h5dopen_f(file, '/quants/pgas', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           p(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
        end do

        ! velocities are in KS!
        ! u0
        call h5dopen_f(file, '/quants/u0', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           u(1+(i-1)*nx2*nx1: i*nx2*nx1)%data(1) = pack(rdata(i,:,:),.true.)
        end do

        ! u1
        call h5dopen_f(file, '/quants/u1', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           u(1+(i-1)*nx2*nx1: i*nx2*nx1)%data(2) = pack(rdata(i,:,:),.true.)
        end do

        ! u2
        call h5dopen_f(file, '/quants/u2', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           u(1+(i-1)*nx2*nx1: i*nx2*nx1)%data(3) = pack(rdata(i,:,:),.true.)
        end do

        ! u3
        call h5dopen_f(file, '/quants/u3', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           u(1+(i-1)*nx2*nx1: i*nx2*nx1)%data(4) = pack(rdata(i,:,:),.true.)
        end do

        ! b0
        call h5dopen_f(file, '/quants/b0', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           b(1+(i-1)*nx2*nx1: i*nx2*nx1)%data(1) = pack(rdata(i,:,:),.true.)
        end do

        ! b1
        call h5dopen_f(file, '/quants/b1', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           b(1+(i-1)*nx2*nx1: i*nx2*nx1)%data(2) = pack(rdata(i,:,:),.true.)
        end do

        ! b2
        call h5dopen_f(file, '/quants/b2', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           b(1+(i-1)*nx2*nx1: i*nx2*nx1)%data(3) = pack(rdata(i,:,:),.true.)
        end do

        ! b3
        call h5dopen_f(file, '/quants/b3', dset, hdferr, H5S_ALL_F)
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
        call h5dclose_f(dset , hdferr)
        do i=1, nx3  ! x2 is the fastest changing index in koral convention 
           b(1+(i-1)*nx2*nx1: i*nx2*nx1)%data(4) = pack(rdata(i,:,:),.true.)
        end do

        ! electron data
        if(has_electrons.eq.1) then
            call h5dopen_f(file, '/quants/te', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
            call h5dclose_f(dset , hdferr)
            do i=1, nx3  ! x2 is the fastest changing index in koral convention 
               te(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
            end do

            call h5dopen_f(file, '/quants/ti', dset, hdferr, H5S_ALL_F)
            call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, dims, hdferr,H5S_ALL_F, dspace, H5P_DEFAULT_F)
            call h5dclose_f(dset , hdferr)
            do i=1, nx3  ! x2 is the fastest changing index in koral convention 
               ti(1+(i-1)*nx2*nx1: i*nx2*nx1) = pack(rdata(i,:,:),.true.)
            end do

            call calc_gammagas(te, ti, gammagas) !calculate variable gamma
        endif 

        ! close file
        call h5fclose_f(file , hdferr)
        call h5close_f(hdferr) 

        write(6,*) 'WARNING: ignoring minpolecell=',minpolecell, ' cells on polar axis!!'


        ! convert velocities, ks2bl, if necessary
        if(outmet.eq.'KS') then
           write(6,*) 'transforming velocities: ks2bl'
           u = uks2ubl(u,dble(r_arr),dble(asim))
           b = uks2ubl(b,dble(r_arr),dble(asim))
        endif

        call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,real(asim))))
        call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,real(asim))))

        ! calculate the bernoulli parameter
        ! AC TODO -- should this be in another function?
        bsq = b*b
        ucov = lower(u)
        bcov = lower(b)
        rhoucon = rho* u%data(1)

        if(has_electrons.eq.1) then  ! AC -- variable gamma case
          t00 = (rho + p*(gammagas/(gammagas-1.)) + bsq)*u%data(1)*ucov%data(1) + (p+0.5*bsq) - b%data(1)*bcov%data(1)
        else
          t00 = (rho + p*(gam/(gam-1.)) + bsq)*u%data(1)*ucov%data(1) + (p+0.5*bsq) - b%data(1)*bcov%data(1)
        endif
        Be = (-t00-rhoucon)/rhoucon

        write(6,*) 'read koral grid sizes', size(r_arr), size(th_arr), size(ph_arr)
!        write(6,*) 'read koral Be ', minval(Be), maxval(Be)
        write(6,*) 'read koral rho ', minval(rho), maxval(rho)
        write(6,*) 'read koral pgas ', minval(p), maxval(p)
!        write(6,*) 'read koral b0', minval(b%data(1)), maxval(b%data(1))
!        write(6,*) 'read koral b1', minval(b%data(2)), maxval(b%data(2))
!        write(6,*) 'read koral b2', minval(b%data(3)), maxval(b%data(3))
!        write(6,*) 'read koral b3', minval(b%data(4)), maxval(b%data(4))
!        write(6,*) 'read koral u0 ', minval(u%data(1)), maxval(u%data(1))
!        write(6,*) 'read koral u1 ', minval(u%data(2)), maxval(u%data(2))
!        write(6,*) 'read koral u2 ', minval(u%data(3)), maxval(u%data(3))
!        write(6,*) 'read koral u3 ', minval(u%data(4)), maxval(u%data(4))
        write(6,*) 'read koral bsq ', minval(bsq), maxval(bsq)
        write(6,*) 'read koral usq ', minval(u*u), maxval(u*u)

!        write(6,*) 'read koral transform coords r ', nx1, minval(r_arr), maxval(r_arr)
!        write(6,*) 'read koral transform coords theta  ', nx2, minval(th_arr), maxval(th_arr)
!        write(6,*) 'read koral transform coords phi', nx3, minval(ph_arr), maxval(ph_arr)

        if(any(rho.eq.0.)) then
           write(6,*) 'rho=0 regions found!'
           minrho = 1.d-6 * minval(merge(rho,1d100,rho.gt.0.))
              write(6,*) '  masking with ', minrho
           rho = merge(rho,minrho,rho.gt.0.)
        end if

        if(any(isnan(rho))) then
           write(6,*) 'nan in rho!'
        end if
        if(any(isnan(p))) then
           write(6,*) 'nan in p!'
        end if
        if(any(isnan(b%data(1)))) then
           write(6,*) 'nan in b0!'
        end if
        if(any(isnan(b%data(2)))) then
           write(6,*) 'nan in b1!'
        end if
        if(any(isnan(b%data(2)))) then
           write(6,*) 'nan in b2!'
        end if
        if(any(isnan(b%data(3)))) then
           write(6,*) 'nan in b3!'
        end if

        if(any(isnan(u%data(1)))) then
           write(6,*) 'nan in u0!'
        end if
        if(any(isnan(u%data(2)))) then
           write(6,*) 'nan in u1!'
        end if
        if(any(isnan(u%data(2)))) then
           write(6,*) 'nan in u2!'
        end if
        if(any(isnan(u%data(3)))) then
           write(6,*) 'nan in u3!'
        end if

        deallocate(ucov); deallocate(bcov)
        deallocate(rhoucon); deallocate(t00); deallocate(bsq);
        if(has_electrons.eq.1) then
          deallocate(gammagas)
        endif 

    end subroutine read_koralh5_data_file

    ! Initialize model parameters & load data from file
    subroutine initialize_koralh5_model(a,ifile,df,hf,ntt,indft,sfac)
        real(kind=8), intent(in) :: a
        real(8),intent(in),optional :: sfac
        character(len=500), intent(in), optional :: ifile
        character(len=500) :: default_ifile='koral.in' !AC
        character(len=500), intent(in), optional :: df,hf
        integer, intent(in), optional :: ntt,indft

        real,dimension(3) :: x1test,x2test,rtest,thtest,x1test2,x2test2

        !write(6,*) 'initialize koral HDF5 model ', trim(df), ' ', trim(hf), ntt, indft
        ! if all inputs are given, assign variables rather than reading from file
        if (present(df)) then
           dfile = df
           hfile = hf
           nt = ntt
           indf = indft
        else
           if (present(ifile)) then
              call read_koralh5_inputs(ifile)
           else
              call read_koralh5_inputs(default_ifile)
           endif
        endif

        ! get scale factor
        if (present(sfac)) then
            scalefac=sfac
        else
            scalefac=1.
        endif

        ! get simulation parameters from datafile header
        call read_koralh5_data_header() 

        if (abs(asim-a).gt.1e-4) then 
           write(6,*) 'ERROR -- Different simulation and grtrans spin values!'
           return
        endif

        ! initialize and load data arrays from file
        n=nx1*nx2*nx3
        call init_koralh5_data(n,n*nt)
        call load_koralh5_data(nt)

        ! test the jetcoords transform
!        x1test=(/0.014,0.11,0.49/)
!        x2test=(/0.01,-0.87,0.44/)
!        write(6,*) 'x1min,x1max',x1min,x1max


!        call transformjet2bl(x1test, x2test, rtest,thtest)
!        call transformbl2jet(rtest,thtest,x1test2,x2test2)
!        write(6,*) 'x1test in', x1test
!        write(6,*) 'r out', rtest
!        write(6,*) 'x1test back', x1test2
!        write(6,*)
!        write(6,*) 'x2test in', x2test
!        write(6,*) 'th out', thtest
!        write(6,*) 'x2test back', x2test2
!        call exit(-1)
    end subroutine initialize_koralh5_model

    ! Read data from file and then transfer to global arrays
    subroutine load_koralh5_data(nt)
        real(8), dimension(:), allocatable :: rho,p,Be,te,ti
        real, dimension(:), allocatable :: vrl, vtl, vpl
        type (four_vector), dimension(:), allocatable :: u, b
        integer, intent(in) :: nt
        character(len=500) :: append
        character(len=500) :: data_file
        integer :: k
        real(8) :: tcur, tstep_test

        allocate(rho(n)); allocate(p(n)); allocate(Be(n))
        allocate(vrl(n)); allocate(vtl(n)); allocate(vpl(n))
        allocate(u(n)); allocate(b(n))
        if(has_electrons.eq.1) then
          allocate(te(n)); allocate(ti(n));
        endif 

        if(has_electrons.eq.1) then
          allocate(te(n)); allocate(ti(n));
        endif 

        do k=1,nt
           !read data from  file
           write(append, fmt='(I4.4)') indf-(k-1)
           if(nt>1) then
            data_file = trim(dfile) // append
           else
            data_file = hfile
           endif
           call read_koralh5_data_file(data_file,tcur,rho,p,te,ti,u,b,Be)

           t(k)=tcur
            
           ! transform velocities to LNRF frame
           call lnrf_frame(real(u%data(2)/u%data(1)),real(u%data(3)/u%data(1)), & 
                           real(u%data(4)/u%data(1)), &
                           r_arr,asim,th_arr,vrl,vtl,vpl)

           write(6,*) 'koral scalefac: ', scalefac
           ! now assign to data arrays
           rho_arr((k-1)*n+1:k*n)=rho*scalefac
           p_arr((k-1)*n+1:k*n)=p*scalefac
           Be_arr((k-1)*n+1:k*n)=Be

           b0_arr((k-1)*n+1:k*n)=b%data(1)*sqrt(scalefac)
           br_arr((k-1)*n+1:k*n)=b%data(2)*sqrt(scalefac)
           bth_arr((k-1)*n+1:k*n)=b%data(3)*sqrt(scalefac)
           bph_arr((k-1)*n+1:k*n)=b%data(4)*sqrt(scalefac)

           u0_arr((k-1)*n+1:k*n)=u%data(1)
           vrl_arr((k-1)*n+1:k*n)=real(vrl)
           vpl_arr((k-1)*n+1:k*n)=real(vpl)
           vtl_arr((k-1)*n+1:k*n)=real(vtl)
           if(has_electrons.eq.1) then
             te_arr((k-1)*n+1:k*n)=te
             ti_arr((k-1)*n+1:k*n)=ti
           endif 

        end do

        ! Need to be able to do light curves even when nload (nt) = 1... Maybe read headers for first two files.
        if (nt.gt.1) then
           ! use default simulation time step unless told otherwise
           tstep=t(1)-t(2)
           tstep_test=t(nt-1)-t(nt)
           if (abs((tstep-tstep_test)/tstep).gt.0.1) then 
              write(6,*) 'WARNING -- Time step changes over dumps!'
           endif
        endif

        deallocate(rho); deallocate(p); deallocate(Be);
        deallocate(vrl); deallocate(vtl); deallocate(vpl)
        deallocate(u); deallocate(b)
        if(has_electrons.eq.1) then
          deallocate(te); deallocate(ti);
        endif 

    end subroutine  load_koralh5_data

    ! Advance to next timestep file
    subroutine advance_koralh5_timestep(dtobs)
        real(8), intent(in) :: dtobs
        integer :: nupdate
        nupdate = floor(dtobs / tstep)
        toffset = toffset+dtobs-tstep*nupdate
        nupdate = nupdate+floor(toffset/tstep)
        toffset = toffset-floor(toffset/tstep)*tstep
        indf = indf+nupdate
        write(6,*) 'advance koral timestep: ',dtobs,tstep,toffset,nupdate,indf 
        call update_koralh5_data(nupdate)
    end subroutine advance_koralh5_timestep

    ! Load new data for next timestep
    subroutine update_koralh5_data(nupdate)
        integer, intent(in) :: nupdate
        integer :: nshift
        if(nupdate.gt.nt) then ! just load all new data
           call load_koralh5_data(nt)

        else ! shift existing data back and then load more
           nshift=n*nupdate+1
           rho_arr(nshift:n*nt)=rho_arr(1:n*nt-nshift)
           p_arr(nshift:n*nt)=p_arr(1:n*nt-nshift)
           Be_arr(nshift:n*nt)=Be_arr(1:n*nt-nshift)
           br_arr(nshift:n*nt)=br_arr(1:n*nt-nshift)
           bth_arr(nshift:n*nt)=bth_arr(1:n*nt-nshift)
           bph_arr(nshift:n*nt)=bph_arr(1:n*nt-nshift)
           b0_arr(nshift:n*nt)=b0_arr(1:n*nt-nshift)
           vrl_arr(nshift:n*nt)=vrl_arr(1:n*nt-nshift)
           vtl_arr(nshift:n*nt)=vtl_arr(1:n*nt-nshift)
           vpl_arr(nshift:n*nt)=vpl_arr(1:n*nt-nshift)
           if(has_electrons.eq.1) then
             te_arr(nshift:n*nt)=te_arr(1:n*nt-nshift)
             ti_arr(nshift:n*nt)=ti_arr(1:n*nt-nshift)
           endif 

           call load_koralh5_data(nupdate)
        end if
    end subroutine update_koralh5_data

    ! Allocate arrays
    subroutine init_koralh5_data(nx, n)
        integer, intent(in) :: nx,n

        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(Be_arr(n)); allocate(u0_arr(n))
        allocate(vrl_arr(n)); allocate(vtl_arr(n)); allocate(vpl_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n))
        allocate(uniqx1(nx1)); allocate(uniqx2(nx2)); allocate(uniqx3(nx3))
        allocate(r_arr(nx)); allocate(th_arr(nx)); allocate(ph_arr(nx))
        allocate(t(nt))
        allocate(r_th_arr(nx1,nx2)); allocate(uniqx2_ext(nx2+2));
        if(has_electrons.eq.1) then
          allocate(te_arr(n)); allocate(ti_arr(n)); 
        endif 
    end subroutine init_koralh5_data

    ! Deallocate arrays
    subroutine del_koralh5_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(Be_arr)
        deallocate(u0_arr); deallocate(vrl_arr); deallocate(vtl_arr); deallocate(vpl_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr); deallocate(bph_arr)
        deallocate(uniqx1); deallocate(uniqx2); deallocate(uniqx3) 
        deallocate(r_arr); deallocate(th_arr); deallocate(ph_arr)
        deallocate(t)
        deallocate(r_th_arr); deallocate(uniqx2_ext)
        if(has_electrons.eq.1) then
          deallocate(te_arr); deallocate(ti_arr);
        endif 
    end subroutine del_koralh5_data

end module fluid_model_koralh5
