! Anantua 2019 semianalytic jet model
! Check conversion to cgs units!

module fluid_model_rrjet

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, calc_u0
      use phys_constants, only: pi, c2, m
      implicit none

      namelist /fluiddata/ betaeconst, betaecrit, ximax, bscl, pscl, pegasratio

      real :: betaeconst, betaecrit, ximax, bscl, pscl, pegasratio

      interface init_rrjet
        module procedure init_rrjet
      end interface
 
      interface rrjet_vals
        module procedure rrjet_vals
      end interface

      contains

        subroutine rrjet_vals(x0,a,rho,p,b,u,bmag)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real :: zcutoff
        real(kind=8), dimension(size(x0)) :: done
        real, dimension(size(x0)) :: x2,x1,zm,zr,theta,s,z,logz,xi,fone
        real, dimension(size(x0)) :: omega, phi, phiprime, current, Bs, Bz 
        real, dimension(size(x0)) :: Bphi,Vs,Vz,Vphi,V0z,Vr,Vth,Br,Bth
        real, dimension(size(x0)) :: u0,b0,vpl0,vrl0,vtl0,p0,rho0,bph0,bth0
        real, dimension(size(x0)) :: rd,td,dzero,vr0,vth0,vph0,dummy
        real, dimension(size(x0)) :: pgas,ww,betagas,betae
        real(8), dimension(size(x0),10) :: metric
        integer :: npts,nx1,nx2
        real, dimension(size(x0)), intent(out) :: rho,p,bmag
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b

        dzero=0d0; done=1d0; fone=1.; 
        
        ! Get spherical coordinates
        npts=size(x0)
        theta=x0%data(3)
        zm=cos(theta)
        
        ! Jet solution is symmetric about xy plane
        ! AC TODO ?? -- is this ok?
        x2=acos(abs(zm))
        theta=x2
        zr=x0%data(2)
        
        ! Get cylindrical coordinates
        s = zr*sin(theta)
        z = zr*cos(theta)
        xi = s*s / z
        where(z.eq.0)
          xi=ximax
        endwhere
                
        ! Get magnetic flux, fluid line angular velocity, and current 
        omega = 0.15*exp(-0.3*xi*xi)
        phi = tanh(0.3*xi) 
        phiprime = 0.3*(1.- phi*phi) 
        current = -2.*xi*omega*phiprime

        ! Get the B-field vector in cylindrical coordinates
        ! !AC - in dimensionless units, the flux through the horizon is 1
        ! !AC - we convert to physical units in  convert_fluidvars_rrjet in fluid.f90

        ! B-field scale = horizon flux / Rg^2
        !bscl = 1.02e4 ! old B-field scale

        Bs = bscl * s * phiprime / (2.*pi*z*z)
        Bz = bscl * phiprime / (pi*z)
        Bphi = bscl * current / (2*pi*s)
        Bphi = Bphi / s !AC -- put in coordinate basis
        
        ! Get the three-velocity in cylindrical coordinates
        logz = log10(z)
        V0z = -1.2 + 2.2*tanh(0.84*logz)
        Vz = V0z * exp(-0.001 * xi**4)
        Vs = s*Vz / (2.* z)
        Vphi = s*omega*(1.-Vz)
        Vphi = Vphi / s !AC -- put in coordinate basis
        
        ! Transform the velocity and B-field to spherical coordinates
        Br = s*Bs/zr + z*Bz/zr
        Bth = z*Bs/(zr*zr) - s*Bz/(zr*zr)
        Vr = s*Vs/zr + z*Vz/zr
        Vth = z*Vs/(zr*zr) - s*Vz/(zr*zr)

        ! AC swap theta component below equatorial plane !TODO -right? 
        where(zm.lt.0.)
            Bth = -Bth
            Vth = -Vth
        end where

        ! Find the timelike u0 and convert V -> u
        metric = kerr_metric(zr,real(x0%data(3)),a)
        u0 = calc_u0(metric,dble(Vr),dble(Vth),dble(Vphi))
        
        u%data(1)=u0
        u%data(2)=Vr*u0
        u%data(3)=Vth*u0
        u%data(4)=Vphi*u0
        call assign_metric(u,transpose(kerr_metric(zr,real(x0%data(3)),a)))

        ! Find b0 and convert B -> b
        b%data(1)=dzero
        b%data(2)=Br
        b%data(3)=Bth
        b%data(4)=Bphi
        call assign_metric(b,transpose(kerr_metric(zr,real(x0%data(3)),a)))

        b0 = b*u
        
        b%data(1)=b0
        b%data(2)=(Br + b0*u%data(2))/u0
        b%data(3)=(Bth + b0*u%data(3))/u0
        b%data(4)=(Bphi + b0*u%data(4))/u0

        ! Protect azimuthal velocities at poles
        ! Compute magnitude of interpolated b-field
        ! and force b^2 > 0 (need to look into numerical issues here):
        call assign_metric(b,transpose(kerr_metric(zr,real(x0%data(3)),a)))
        bmag=b*b;
        bmag=merge(bmag,dzero,bmag.ge.0d0)
        bmag=sqrt(bmag)

        ! Gas pressure fitting form from Richard : TODO what should pscl be? paramatrize better?
        ww = 2.5*sqrt(z)
        pgas = pscl*(z**(-4.5)) * ((0.7*s/ww)**3) * exp(-(0.45*s/ww)**4)
        betagas = 8*pi*pgas/(bmag*bmag + 1.e-20)
        
        if ((maxval(bmag*bmag).ge.0.).and.(minval(xi).lt.ximax)) then
            !write(6,*) 'max pgas', maxval(pgas)
            !write(6,*) 'max bsq', maxval(bmag*bmag + 1.e-20)
            !write(6,*) 'max pgas/bsq', (8*pi*maxval(pgas))/maxval(bmag*bmag + 1.e-20)
        endif

        ! Find the nonthermal electron pressure (betae * magnetic pressure)
        !p = betaeconst * (bmag*bmag) / (8*pi) !TODO -- 8pi or not? 

        if (pegasratio.gt.0) then
            p = pegasratio*pgas
        else if(betaeconst.eq.-1) then
            ! use power law for beta
            p = 1.e-6 * (0.08*z)**4
        else
            betae = betaeconst*fone
            if(betaecrit.gt.0) then
                where(betagas.gt.betaecrit)
                    betae = betaeconst * exp(-(betagas/betaecrit-1))
                end where
            endif
           
            p = betae * (bmag*bmag) / (8*pi)

        endif
        
        rho = p ! Convert pressure to number density in  convert_fluidvars_rrjet

        ! Zero everything where xi>ximax and z<2
        zcutoff=0.d0
        where((xi.gt.ximax).or.(abs(z).le.zcutoff))
           !u%data(1) = dzero;
           !u%data(2) = dzero;
           !u%data(3) = dzero;
           !u%data(4) = dzero;
           !b%data(1) = dzero;
           !b%data(2) = dzero;
           !b%data(3) = dzero;
           !b%data(4) = dzero;
           rho = dzero;
           bmag = dzero;
           p = dzero;
        end where
        
        ! !AC TODO ?? 
        ! !AC Cut off emission from below the eq. plane:
        !rho=merge(rho,rho*0.,zm.ge.0)
        !p = merge(p,p*0.,zm.ge.0)
        end subroutine rrjet_vals
    
        subroutine init_rrjet(betaeconst0, betaecrit0, ximax0, bscl0,pscl0,pegasratio0)
          real(8), intent(in) :: betaeconst0, betaecrit0, ximax0, bscl0,pscl0,pegasratio0
          betaeconst = betaeconst0
          betaecrit = betaecrit0
          ximax=ximax0
          bscl=bscl0
          pscl=pscl0
          pegasratio=pegasratio0
        write(6,*) 'read: ',betaeconst,betaecrit,ximax,bscl,pscl,pegasratio
        close(unit=8)
        end subroutine init_rrjet

      end module fluid_model_rrjet
