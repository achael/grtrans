! Anantua 2019 semianalytic jet model
! Check conversion to cgs units!

module fluid_model_rrjet

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, calc_u0
      use phys_constants, only: pi, c2, m
      implicit none

      namelist /fluiddata/ betaeconst, ximax

      real :: betaeconst, ximax

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
        real(kind=8), dimension(size(x0)) :: done
        real, dimension(size(x0)) :: x2,x1,zm,zr,theta,s,z,logz,xi,fone
        real, dimension(size(x0)) :: omega, phi, phiprime, current, Bs, Bz, Bphi,Vs,Vz,Vphi,V0z,Vr,Vth,Br,Bth
        real, dimension(size(x0)) ::  u0,b0,vpl0,vrl0,vtl0,p0,rho0,bph0,bth0,rd,td,dzero,vr0,vth0,vph0,dummy   
        real(8), dimension(size(x0),10) :: metric
        integer :: npts,nx1,nx2
        real, dimension(size(x0)), intent(out) :: rho,p,bmag
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b

!        write(6,*) 'rrjet: ',size(x0)
        dzero=0d0; done=1d0; fone=1.; 
        
        ! Get spherical coordinates
        npts=size(x0)
        theta=x0%data(3)
        zm=cos(theta)
        
        ! Jet solution is symmetric about xy plane
        ! ANDREW ?? -- is this ok?
        x2=acos(abs(zm)); theta=x2
        zr=x0%data(2)
        
        ! Get cylindrical coordinates
        s = zr*sin(theta)
        z = zr*cos(theta)
        xi = s*s / z
        where(z.eq.0)
          xi=1.e10
        endwhere
                
        ! Get magnetic flux, fluid line angular velocity, and current 
        omega = 0.15*exp(-0.3*xi*xi)
        phi = tanh(0.3*xi) 
        phiprime = 0.3*(1.- phi*phi) 
        current = -2.*xi*omega*phiprime

        ! Get the B-field vector in cylindrical coordinates
        ! !AC - in dimensionless units, the flux through the horizon is 1
        ! !AC - we convert to physical units in  convert_fluidvars_rrjet in fluid.f90
        Bs = 1.02e4 * s * phiprime / (2.*pi*z*z)
        Bz = 1.02e4 * phiprime / (pi*z)
        !Bz = 1.02e4 * phiprime / (2.*pi*z)
        Bphi = 1.02e4 * current / (2*pi*s)
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
        !b%data(2)=u0*(Br + b0*u%data(2))
        !b%data(3)=u0*(Bth + b0*u%data(3))
        !b%data(4)=u0*(Bphi + b0*u%data(4))
                        
        ! Protect azimuthal velocities at poles
        ! Compute magnitude of interpolated b-field
        ! and force b^2 > 0 (need to look into numerical issues here):
        call assign_metric(b,transpose(kerr_metric(zr,real(x0%data(3)),a)))
        bmag=b*b;
        bmag=merge(bmag,dzero,bmag.ge.0d0)
        bmag=sqrt(bmag)
!        write(6,*) 'maxr: ',maxval(uniqr),maxval(zr)
!        write(6,*) 'bmag: ',maxval(bmag),minval(bmag)
!        write(6,*) 'bmag: ',bmag

        ! Find the nonthermal electron pressure (betae * magnetic pressure)
        !p = 0.5 * betaeconst * (bmag*bmag) !/ (8*pi)

        if(betaeconst.eq.-1) then
            ! use power law for beta
            p = 1.e-6 * (0.08*z)**4
        else
            p = betaeconst * (bmag*bmag) !/ (8*pi)
        endif
        
        rho = p ! Convert pressure to number density in  convert_fluidvars_rrjet using emis params

        ! Zero everything where xi>ximax and z<2
        where((xi.gt.ximax).or.(abs(z).le.2d0))
           !u%data(1) = dzero;
           !u%data(2) = dzero;
           !u%data(3) = dzero;
           !u%data(4) = dzero;
           !b%data(1) = dzero;
           !b%data(2) = dzero;
           !b%data(3) = dzero;
           !b%data(4) = dzero;
           rho = dzero;
           bmag=dzero;
           p=dzero;
        endwhere
        
        ! !AC Cut off emission from below the eq. plane:
        rho=merge(rho,rho*0.,zm.ge.0)
        p = merge(p,p*0.,zm.ge.0)
        end subroutine rrjet_vals
    
        subroutine init_rrjet(betaeconst0, ximax0)
        !character(len=20), intent(in) :: ifile
        !open(unit=8,file=ifile,form='formatted',status='old')
          !read(8,nml=fluiddata)
          real(8), intent(in) :: betaeconst0,  ximax0
          betaeconst = betaeconst0
          ximax=ximax0
        write(6,*) 'read: ',betaeconst,ximax
        close(unit=8)
        end subroutine init_rrjet

      end module fluid_model_rrjet
