module fluid_model

  use class_four_vector
  use phys_constants, GC=>G
  use interpolate, only: interp, get_weight
  use kerr, only: kerr_metric, lnrf_frame, calc_rms, krolikc, calc_polvec, calc_u0, rms_vel
  use fluid_model_sphacc, only: sphacc_vals, init_sphacc, del_sphacc
  use fluid_model_sariaf, only: sariaf_vals, init_sariaf, del_sariaf
  use fluid_model_toy, only: toy_vals, init_toy, del_toy
  use fluid_model_powerlaw, only: powerlaw_vals, init_powerlaw, del_powerlaw
  use fluid_model_ffjet, only: initialize_ffjet_model, del_ffjet_data, &
       ffjet_vals
  use fluid_model_rrjet, only: init_rrjet, rrjet_vals      
  use fluid_model_phatdisk, only: phatdisk_vals, init_phatdisk, del_phatdisk, freq_tab
  use fluid_model_thindisk, only: thindisk_vals, init_thindisk
  use fluid_model_numdisk, only: initialize_numdisk_model, del_numdisk_data, &
       numdisk_vals
  use fluid_model_hotspot, only: hotspot_vals, init_hotspot, advance_hotspot_timestep
  use fluid_model_hotspot_schnittman, only: init_schnittman_hotspot, & 
       advance_schnittman_hotspot_timestep, schnittman_hotspot_vals
  use fluid_model_harm, only: initialize_harm_model, del_harm_data, harm_vals, &
       advance_harm_timestep
  use fluid_model_harm3d, only: initialize_harm3d_model, del_harm3d_data, harm3d_vals, &
      advance_harm3d_timestep
  use fluid_model_harmpi, only: initialize_harmpi_model, del_harmpi_data, harmpi_vals, &
      advance_harmpi_timestep
  use fluid_model_iharm, only: initialize_iharm_model, del_iharm_data, iharm_vals, &
      advance_iharm_timestep
  use fluid_model_koral, only: initialize_koral_model, del_koral_data, koral_vals, &
       advance_koral_timestep           
  use fluid_model_koral3d, only: initialize_koral3d_model, del_koral3d_data, koral3d_vals, &
       advance_koral3d_timestep
  use fluid_model_koralh5, only: initialize_koralh5_model, del_koralh5_data, koralh5_vals, &
       advance_koralh5_timestep
  use fluid_model_thickdisk, only: initialize_thickdisk_model, del_thickdisk_data, thickdisk_vals, &
       advance_thickdisk_timestep
  use fluid_model_mb09, only: initialize_mb09_model, del_mb09_data, mb09_vals, &
       advance_mb09_timestep
  use calc_gmin, only: calc_gmin_subroutine

  implicit none

  ! mnemonics for different nonthermal electron source models
  integer, parameter :: CONST=0,TAIL=1 

  ! mnemonics for different fluid models
  integer, parameter :: DUMMY=0,SPHACC=1,THINDISK=2,RIAF=3,HOTSPOT=4,PHATDISK=5,SCHNITTMAN=6, &
                        COSMOS=10,MB=11,HARM=12,FFJET=13,NUMDISK=14,THICKDISK=15,MB09=16, &
                        SARIAF=17,POWERLAW=18,HARM3D=19,HARMPI=20,TOY=21,KORAL=22,KORALNTH=23, &
                        SHELL=24,KORAL3D=25,KORAL3D_DISK=26,KORAL3D_TOPJET=27,KORAL3D_BOTJET=28, & 
                        IHARM=29, RRJET=30, KORALH5=31

  ! placeholder values for nothermal electrons 
  integer :: nrelbin=0
  real :: bingammamin=1., bingammamax=1.
 
  ! flag for if we have electrons or not -- !AC TODO handle this better
  integer :: has_electrons=0


  ! fluid data type
  type fluid
    integer :: model, nfreq, nrelbin
    real :: rin,bingammamin,bingammamax,sigcut
    real, dimension(:), allocatable :: rho,p,te,ti,bmag,rho2, &
         kela,kelb,kelc,keld,Be
    real, dimension(:,:), allocatable :: fnu
    real, dimension(:,:), allocatable :: nnth
    type (four_vector), dimension(:), allocatable :: u,b
  end type

  ! fluid args data type
  type fluid_args
     character(len=500) :: dfile,hfile,gfile,sim
     integer :: nt,indf,nfiles,jonfix,nw,nfreq_tab,nr,offset, &
          dindf,magcrit,bl06
     real(8) :: rspot,r0spot,n0spot,tscl,rscl,wmin,wmax,fmin, &
          fmax,rmax,sigt,fcol,mdot,mbh,nscl,nnthscl,nnthp,beta, &
          np,tp,rin,rout,thin,thout,phiin,phiout,scalefac,sigcut, &
          betaeconst,betaecrit,ximax,bscl,pscl,pegasratio
  end type

  ! source params data type
  type source_params
    real(kind=8) :: nfac,bfac,mbh,mdot,p1,p2,gmax,fpositron,gminval, &
         jetalphaval,muval,sigcut,ximax,betaeconst,betaecrit,bscl,pscl,pegasratio
    real(kind=8), dimension(:), allocatable :: gmin,jetalpha,mu
    integer :: type
  end type

  ! interfaces
  ! TODO: ultimately all source params stuff should probably go in own file / module
  interface initialize_source_params
     module procedure initialize_source_params
  end interface

  interface del_source_params
     module procedure del_source_params
  end interface
  
  interface assign_source_params
     module procedure assign_source_params
  end interface

  interface get_fluid_vars
!    module procedure get_fluid_vars_single
    module procedure get_fluid_vars_arr
  end interface

  interface convert_fluid_vars
    module procedure convert_fluid_vars_arr
  end interface convert_fluid_vars

  interface initialize_fluid_model
    module procedure initialize_fluid_model
  end interface

  interface assign_fluid_args
     module procedure assign_fluid_args
  end interface

  interface advance_fluid_timestep
     module procedure advance_fluid_timestep
  end interface

  interface load_fluid_model
    module procedure load_fluid_model
  end interface

  interface unload_fluid_model
    module procedure unload_fluid_model
  end interface unload_fluid_model

  contains

    ! Assign fluid model arguments to fluid_args data type
    subroutine assign_fluid_args(fargs,dfile,hfile,gfile,sim,nt,indf,nfiles,jonfix, &
                                 nw,nfreq_tab,nr,offset,dindf,magcrit, &
                                 rspot,r0spot, n0spot,tscl,rscl, &
                                 wmin,wmax,fmin,fmax,rmax,sigt,fcol, &
                                 mdot,mbh,nscl,nnthscl,nnthp,beta,bl06,np,tp, &
                                 rin,rout,thin,thout,phiin,phiout,scalefac,sigcut, &
                                 betaeconst,betaecrit,ximax,bscl,pscl,pegasratio)

          type (fluid_args), intent(inout) :: fargs
          character(len=500), intent(in) :: dfile,hfile,gfile,sim
          integer, intent(in) :: nt,indf,nfiles,jonfix,nw,nfreq_tab,nr,offset,dindf, &
               magcrit,bl06
          real(8), intent(in) :: rspot,r0spot,n0spot,tscl,rscl,wmin,wmax,fmin, &
               fmax,rmax,sigt,fcol,mdot,mbh,nscl,nnthscl,nnthp,beta,np,tp, &
               rin,rout,thin,thout,phiin,phiout,scalefac,sigcut,betaeconst,betaecrit,ximax,&
               bscl,pscl,pegasratio
          fargs%dfile = dfile; fargs%hfile = hfile; fargs%gfile=gfile
          fargs%sim = sim; fargs%nt = nt; fargs%indf = indf; fargs%nfiles = nfiles
          fargs%jonfix = jonfix; fargs%nw = nw; fargs%nfreq_tab = nfreq_tab
          fargs%nr = nr; fargs%offset = offset; fargs%dindf = dindf
          fargs%magcrit = magcrit; fargs%rspot = rspot; fargs%r0spot = r0spot
          fargs%n0spot = n0spot; fargs%tscl = tscl; fargs%rscl = rscl
          fargs%wmin = wmin; fargs%wmax = wmax; fargs%fmin = fmin
          fargs%fmax = fmax; fargs%rmax = rmax; fargs%sigt = sigt
          fargs%mbh = mbh; fargs%fcol = fcol; fargs%mdot = mdot
          fargs%nscl = nscl; fargs%nnthscl = nnthscl; fargs%nnthp = nnthp
          fargs%beta = beta; fargs%bl06 = bl06; fargs%np = np; fargs%tp=tp
          fargs%rin = rin; fargs%rout = rout; fargs%thin = thin
          fargs%thout = thout; fargs%phiin = phiin; fargs%phiout = phiout
          fargs%scalefac=scalefac; fargs%sigcut=sigcut;
          fargs%betaeconst=betaeconst; fargs%ximax=ximax
          fargs%betaecrit=betaecrit; fargs%bscl=bscl; fargs%pscl=pscl 
          fargs%pegasratio=pegasratio
    end subroutine assign_fluid_args

    ! Load fluid data from file and/or set parameters
    subroutine load_fluid_model(fname,a,fargs)
        real(kind=8), intent(in) :: a
        character(len=500), intent(in) :: fname
        character(len=500) :: ifile
        type (fluid_args) :: fargs

        if(fname=='COSMOS') then
            write(6,*) 'COSMOS fluid model removed!'
!           call init_cosmos(a,fargs)
        elseif(fname=='MB') then
            write(6,*) 'MB fluid model removed!'
!           call initialize_mb_model(a)
        elseif(fname=='RIAF') then
            write(6,*) 'RIAF fluid model removed! try SARIAF'
!           call init_riaf(a)
        elseif(fname=='SARIAF') then
            call init_sariaf(real(fargs%nscl),real(fargs%tscl),real(fargs%nnthscl), &
                             real(fargs%nnthp),real(fargs%beta),fargs%bl06) 
        elseif(fname=='POWERLAW') then
            call init_powerlaw(real(fargs%nscl),real(fargs%tscl),real(fargs%nnthscl), &
                               real(fargs%nnthp),real(fargs%beta),real(fargs%np),real(fargs%tp), &
                               real(fargs%rin),real(fargs%rout),real(fargs%thin),real(fargs%thout), &
                               real(fargs%phiin),real(fargs%phiout)) 
        elseif(fname=='TOY') then
            call init_toy(fargs%nscl,fargs%tscl,fargs%beta,fargs%bl06)

        elseif(fname=='THICKDISK') then
            call initialize_thickdisk_model(a,1,ifile,fargs%gfile,fargs%dfile,fargs%nt, &
                                            fargs%nfiles,fargs%indf,fargs%jonfix,fargs%offset,fargs%sim, &
                                            fargs%dindf,fargs%magcrit)
        elseif(fname=='MB09') then
            call initialize_mb09_model(a,1,ifile,fargs%gfile,fargs%dfile,fargs%nt, &
                                       fargs%nfiles,fargs%indf,fargs%jonfix,fargs%sim)
        elseif(fname=='HARM') then
            call initialize_harm_model(a,ifile,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf)
        elseif(fname=='HARM3D') then
            call initialize_harm3d_model(a,ifile,fargs%dfile,fargs%hfile,fargs%gfile,fargs%nt,fargs%indf)
        elseif(fname=='HARMPI') then
            call initialize_harmpi_model(a,ifile,fargs%dfile,fargs%hfile,fargs%gfile,fargs%nt,fargs%indf)
        elseif(fname=='IHARM') then
            call initialize_iharm_model(a,ifile,fargs%dfile,fargs%hfile,fargs%gfile,fargs%nt,fargs%indf)
        elseif(fname=='KORAL') then
            call initialize_koral_model(a,ifile,.false.,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf, & 
                                        fargs%scalefac,nrelbin,bingammamin,bingammamax)
        elseif(fname=='KORALNTH') then
            call initialize_koral_model(a,ifile,.true.,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf, &
                                        fargs%scalefac,nrelbin,bingammamin,bingammamax)            
        elseif(fname=='KORAL3D') then
            call initialize_koral3d_model(a,ifile,.false.,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf, & 
                                          fargs%scalefac,nrelbin,bingammamin,bingammamax)
        elseif(fname=='KORAL3D_DISK') then
            call initialize_koral3d_model(a,ifile,.false.,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf, & 
                                          fargs%scalefac,nrelbin,bingammamin,bingammamax)                      
        elseif(fname=='KORAL3D_TOPJET') then
            call initialize_koral3d_model(a,ifile,.false.,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf, & 
                                          fargs%scalefac,nrelbin,bingammamin,bingammamax)                      
        elseif(fname=='KORAL3D_BOTJET') then
            call initialize_koral3d_model(a,ifile,.false.,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf, & 
                                          fargs%scalefac,nrelbin,bingammamin,bingammamax)
        elseif(fname=='KORALH5') then
            call initialize_koralh5_model(a,ifile,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf,fargs%scalefac)                      
        elseif(fname=='SPHACC') then
            call init_sphacc()
        elseif(fname=='FFJET') then
            call initialize_ffjet_model(a,ifile,fargs%dfile)
        elseif(fname=='THINDISK') then
            call init_thindisk(real(a),ifile,real(fargs%mdot),real(fargs%mbh),real(fargs%rin),real(fargs%rout))
        elseif(fname=='PHATDISK') then
            call init_phatdisk(real(a),ifile,fargs%nw,real(fargs%wmin), &
                               real(fargs%wmax),fargs%nfreq_tab,real(fargs%fmin), &
                               real(fargs%fmax),fargs%nr, real(fargs%sigt), &
                               real(fargs%fcol))
        elseif(fname=='NUMDISK') then
            call initialize_numdisk_model(ifile,fargs%dfile,real(fargs%tscl), real(fargs%rscl))
        elseif(fname=='HOTSPOT') then
            call init_hotspot(ifile,real(fargs%rspot),real(fargs%r0spot), &
                              real(fargs%n0spot),fargs%bl06)
        elseif(fname=='SCHNITTMAN') then
            call init_schnittman_hotspot(ifile,real(fargs%rspot),real(fargs%r0spot), real(fargs%n0spot))
        elseif(fname=='RRJET') then
            call init_rrjet(fargs%betaeconst,fargs%betaecrit, fargs%ximax,fargs%bscl,fargs%pscl,fargs%pegasratio)
           
        endif
    end subroutine load_fluid_model

    ! Advance numeric simulation timestep
    subroutine advance_fluid_timestep(fname,dt)
        real(kind=8), intent(in) :: dt
        character(len=500), intent(in) :: fname

        if(fname=='COSMOS') then
!           call advance_cosmos_timestep(dt)
        elseif(fname=='MB') then
!           call advance_mb_timestep(dt)
        elseif(fname=='HARM') then
            call advance_harm_timestep(dt)
        elseif(fname=='HARM3D') then
            call advance_harm3d_timestep(dt)
        elseif(fname=='HARMPI') then
            call advance_harmpi_timestep(dt)
        elseif(fname=='IHARM') then
            call advance_iharm_timestep(dt)
        elseif(fname=='THICKDISK') then 
            call advance_thickdisk_timestep(dt)
        elseif(fname=='MB09') then 
            call advance_mb09_timestep(dt)
        elseif(fname=='KORAL') then
            call advance_koral_timestep(dt)
        elseif(fname=='KORALNTH') then
            call advance_koral_timestep(dt)
        elseif(fname=='KORAL3D') then
            call advance_koral3d_timestep(dt)
        elseif(fname=='KORAL3D_DISK') then
            call advance_koral3d_timestep(dt)
        elseif(fname=='KORAL3D_TOPJET') then
            call advance_koral3d_timestep(dt)
        elseif(fname=='KORAL3D_BOTJET') then
            call advance_koral3d_timestep(dt)
        elseif(fname=='KORALH5') then
            call advance_koralh5_timestep(dt)
        elseif(fname=='HOTSPOT') then
            call advance_hotspot_timestep(real(dt))
        elseif(fname=='SCHNITTMAN') then
            call advance_schnittman_hotspot_timestep(real(dt))
        endif
    end subroutine advance_fluid_timestep

    ! Specify the fluid model and allocate necessary arrays
    subroutine initialize_fluid_model(f,fname,a,nup)
        character(len=500), intent(in) :: fname
        type (fluid), intent(out) :: f
        integer, intent(in) :: nup
        real(kind=8), intent(in) :: a

        if(fname=='PHATDISK') then
           f%model=PHATDISK; f%nfreq=size(freq_tab)
           allocate(f%u(nup)); allocate(f%fnu(nup,f%nfreq))
           allocate(f%b(nup))
        else
           allocate(f%u(nup))
           allocate(f%rho(nup)) 
           allocate(f%b(nup))
           if(fname=='THINDISK') then
              f%model=THINDISK ! note: rho is used to store T for thindisk.
           elseif(fname=='NUMDISK') then
              f%model=NUMDISK
           else
              allocate(f%bmag(nup))
              allocate(f%p(nup))

              if(fname=='COSMOS') then
                 f%model=COSMOS
              elseif(fname=='MB') then
                 f%model=MB
              elseif(fname=='HARM') then
                 f%model=HARM
              elseif(fname=='HARM3D') then
                 f%model=HARM3D
              elseif(fname=='HARMPI') then
                 f%model=HARMPI
                 allocate(f%kela(nup))
                 allocate(f%kelb(nup))
                 allocate(f%kelc(nup))
                 allocate(f%keld(nup))
              elseif(fname=='IHARM') then
                 f%model=IHARM
              elseif(fname=='THICKDISK') then
                 f%model=THICKDISK
              elseif(fname=='MB09') then
                 f%model=MB09
              elseif(fname=='KORAL') then
                 f%model=KORAL
                 allocate(f%Be(nup))
                 allocate(f%te(nup))
                 allocate(f%ti(nup))
              elseif(fname=='KORALNTH') then
                 f%model=KORALNTH
                 allocate(f%Be(nup))
                 allocate(f%te(nup))
                 allocate(f%ti(nup))
                 allocate(f%nnth(nup, nrelbin))
                 f%nrelbin=nrelbin
                 f%bingammamin=bingammamin
                 f%bingammamax=bingammamax 
              elseif(fname=='KORAL3D') then
                 f%model=KORAL3D
                 allocate(f%Be(nup))
                 allocate(f%te(nup))
                 allocate(f%ti(nup))
              elseif(fname=='KORAL3D_DISK') then
                 f%model=KORAL3D_DISK
                 allocate(f%Be(nup))
                 allocate(f%te(nup))
                 allocate(f%ti(nup))
              elseif(fname=='KORAL3D_BOTJET') then
                 f%model=KORAL3D_BOTJET
                 allocate(f%Be(nup))
                 allocate(f%te(nup))
                 allocate(f%ti(nup))
              elseif(fname=='KORAL3D_TOPJET') then
                 f%model=KORAL3D_TOPJET
                 allocate(f%Be(nup))
                 allocate(f%te(nup))
                 allocate(f%ti(nup))
              elseif(fname=='KORALH5') then
                 f%model=KORALH5
                 allocate(f%Be(nup))
                 allocate(f%te(nup))
                 allocate(f%ti(nup))
              elseif(fname=='FFJET') then
                 f%model=FFJET
              elseif(fname=='RRJET') then
                 f%model=RRJET                 
              elseif(fname=='SPHACC') then
                 f%model=SPHACC
              elseif(fname=='RIAF') then
                 f%model=RIAF
              elseif(fname=='HOTSPOT') then
                 f%model=HOTSPOT
              elseif(fname=='SCHNITTMAN') then
                 f%model=SCHNITTMAN
              elseif(fname=='SARIAF') then
                 f%model=SARIAF
                 allocate(f%rho2(nup))
              elseif(fname=='TOY') then
                 f%model=TOY
              elseif(fname=='POWERLAW') then
                 f%model=POWERLAW
                 allocate(f%rho2(nup))
              else
                 write(6,*) 'WARNING: Unsupported fluid model -- using DUMMY'
                 f%model=DUMMY
              endif
           endif
        endif
    end subroutine initialize_fluid_model
 
    ! Deallocate model-specific arrays when we're done loading data
    subroutine unload_fluid_model(fname)
        character(len=500), intent(in) :: fname
        if(fname=='COSMOS') then
!           call del_cosmos()
        elseif(fname=='MB') then
!           call del_mb(a)
        elseif(fname=='HARM') then
            call del_harm_data()
        elseif(fname=='HARM3D') then
            call del_harm3d_data()
        elseif(fname=='HARMPI') then
            call del_harmpi_data()
        elseif(fname=='IHARM') then
            call del_iharm_data()
        elseif(fname=='THICKDISK') then
            call del_thickdisk_data()
        elseif(fname=='MB09') then
            call del_mb09_data()
        elseif(fname=='KORAL') then
            call del_koral_data()
        elseif(fname=='KORALNTH') then
            call del_koral_data() 
        elseif(fname=='KORAL3D') then
            call del_koral3d_data()
        elseif(fname=='KORAL3D_DISK') then
            call del_koral3d_data()
        elseif(fname=='KORAL3D_TOPJET') then
            call del_koral3d_data()
        elseif(fname=='KORAL3D_BOTJET') then
            call del_koral3d_data()
        elseif(fname=='KORALH5') then
            call del_koralh5_data()
        elseif(fname=='FFJET') then
            call del_ffjet_data()          
        elseif(fname=='PHATDISK') then
            call del_phatdisk()
        elseif(fname=='NUMDISK') then
            call del_numdisk_data()
        elseif(fname=='SPHACC') then
            call del_sphacc()
        elseif(fname=='SARIAF') then
            call del_sariaf()
        elseif(fname=='TOY') then
            call del_toy()
        elseif(fname=='POWERLAW') then
            call del_powerlaw()
        endif
    end subroutine unload_fluid_model

    ! Deallocate general arrays in fluid data when we're done loading data
    subroutine del_fluid_model(f)
        type (fluid), intent(inout) :: f
        if(f%model==PHATDISK) then
           deallocate(f%u); deallocate(f%fnu)
           deallocate(f%b)
        else
           deallocate(f%u); deallocate(f%rho); deallocate(f%b)
           if(f%model.ne.THINDISK.and.f%model.ne.NUMDISK) then
              deallocate(f%p)
              deallocate(f%bmag)
           endif
           if(f%model==SARIAF.or.f%model==POWERLAW) deallocate(f%rho2)
           if(f%model==KORALNTH) deallocate(f%nnth)
           if(f%model==KORAL.or.f%model==KORAL3D.or.f%model==KORAL3D_TOPJET) then
              deallocate(f%Be)
              deallocate(f%ti)
              deallocate(f%te)
           endif !AC TODO too dumb to figure out line continuation at 5 am
           if(f%model==KORAL3D_BOTJET.or.f%model==KORAL3D_DISK.or.f%model==KORALH5) then
              deallocate(f%Be)
              deallocate(f%ti)
              deallocate(f%te)
           endif
           if(f%model==HARMPI) then
              deallocate(f%kela); deallocate(f%kelb)
              deallocate(f%kelc); deallocate(f%keld)
           endif
        endif
        f%model=-1
    end subroutine del_fluid_model
   
    ! Get fluid variables at  an array of points x0
    subroutine get_fluid_vars_arr(x0,k0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0, k0
        type (fluid), intent(inout) :: f
        real(kind=8), intent(in) :: a

        SELECT CASE(f%model)
          CASE (SPHACC)
            call get_sphacc_fluidvars(x0,f)
          CASE (FFJET)
             call get_ffjet_fluidvars(x0,real(a),f)
          CASE (RRJET)
            call get_rrjet_fluidvars(x0,real(a),f)             
          CASE (THINDISK)
            call get_thindisk_fluidvars(x0,k0,real(a),f)
          CASE (PHATDISK)
            call get_phat_fluidvars(x0,k0,real(a),f)
          CASE (NUMDISK)
            call get_numdisk_fluidvars(x0,k0,real(a),f)
          CASE (HOTSPOT)
            call get_hotspot_fluidvars(x0,k0,real(a),f)
          CASE (SCHNITTMAN)
            call get_schnittman_hotspot_fluidvars(x0,real(a),f)
          CASE (HARM)
            call get_harm_fluidvars(x0,real(a),f)
          CASE (HARM3D)
            call get_harm3d_fluidvars(x0,real(a),f)
          CASE (HARMPI)
            call get_harmpi_fluidvars(x0,real(a),f)
          CASE (IHARM)
            call get_iharm_fluidvars(x0,real(a),f)
          CASE (THICKDISK)
            call get_thickdisk_fluidvars(x0,real(a),f)
          CASE (MB09)
            call get_mb09_fluidvars(x0,real(a),f)
          CASE (KORAL)
            call get_koral_fluidvars(x0,real(a),f)
          CASE (KORALNTH)
            call get_koral_fluidvars(x0,real(a),f)             
          CASE (KORAL3D)
            call get_koral3d_fluidvars(x0,real(a),f,0)
          CASE (KORAL3D_DISK)
            call get_koral3d_fluidvars(x0,real(a),f,1)
          CASE (KORAL3D_TOPJET)
            call get_koral3d_fluidvars(x0,real(a),f,2)
          CASE (KORAL3D_BOTJET)
            call get_koral3d_fluidvars(x0,real(a),f,3)
          CASE (KORALH5)
            call get_koralh5_fluidvars(x0,real(a),f,0)
          CASE (SARIAF)
            call get_sariaf_fluidvars(x0,real(a),f)
          CASE (TOY)
            call get_toy_fluidvars(x0,a,f)
          CASE (POWERLAW)
            call get_powerlaw_fluidvars(x0,real(a),f)
          CASE (DUMMY)
        END SELECT
    end subroutine get_fluid_vars_arr

    ! After sampling fluid variables at points x0, convert them to cgs units
    subroutine convert_fluid_vars_arr(f,ncgs,ncgsnth,nnthcgs,bcgs,tcgs,fnuvals,freqvals,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(inout) :: sp
        real(kind=8), intent(out), dimension(size(f%rho)) :: ncgs,ncgsnth,bcgs,tcgs
        real(kind=8), intent(out), dimension(size(f%rho),f%nrelbin) :: nnthcgs !nnthcgs is binned, ncgsnth is not
        real(kind=8), intent(out), dimension(:), allocatable :: freqvals
        real(kind=8), intent(out), dimension(:,:), allocatable :: fnuvals
        ncgs=0d0; ncgsnth=0d0; bcgs=0d0; tcgs=0d0

        SELECT CASE(f%model)
          CASE (SPHACC)
            call convert_fluidvars_sphacc(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (FFJET)
            call convert_fluidvars_ffjet(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (RRJET)
            call convert_fluidvars_rrjet(f,ncgs,ncgsnth,bcgs,tcgs,sp)             
          CASE (THINDISK)
            call convert_fluidvars_thindisk(f,tcgs,ncgs)
          CASE(PHATDISK)
            call convert_fluidvars_phatdisk(f,fnuvals,freqvals)
          CASE(NUMDISK)
            call convert_fluidvars_numdisk(f,tcgs,ncgs)
          CASE (HOTSPOT)
            call convert_fluidvars_hotspot(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (SCHNITTMAN)
            call convert_fluidvars_schnittman_hotspot(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (HARM)
            call convert_fluidvars_harm(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (HARM3D)
            call convert_fluidvars_harm3d(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (HARMPI)
            call convert_fluidvars_harmpi(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (IHARM)
            call convert_fluidvars_iharm(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (THICKDISK)
            call convert_fluidvars_thickdisk(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (MB09)
            call convert_fluidvars_mb09(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (KORAL)
            call convert_fluidvars_koral(f,ncgs,ncgsnth,bcgs,tcgs,nnthcgs,sp,0,0)
          CASE (KORALNTH)
            call convert_fluidvars_koral(f,ncgs,ncgsnth,bcgs,tcgs,nnthcgs,sp,0,0)             
          CASE (KORAL3D)
            call convert_fluidvars_koral(f,ncgs,ncgsnth,bcgs,tcgs,nnthcgs,sp,0,0)
          CASE (KORAL3D_DISK)
            call convert_fluidvars_koral(f,ncgs,ncgsnth,bcgs,tcgs,nnthcgs,sp,1,0)
          CASE (KORAL3D_TOPJET)
            call convert_fluidvars_koral(f,ncgs,ncgsnth,bcgs,tcgs,nnthcgs,sp,2,0)
          CASE (KORAL3D_BOTJET)
            call convert_fluidvars_koral(f,ncgs,ncgsnth,bcgs,tcgs,nnthcgs,sp,3,0)
          CASE (KORALH5) 
            call convert_fluidvars_koral(f,ncgs,ncgsnth,bcgs,tcgs,nnthcgs,sp,0,1)
          CASE (SARIAF)
            call convert_fluidvars_sariaf(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (TOY)
            call convert_fluidvars_toy(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (POWERLAW)
            call convert_fluidvars_powerlaw(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        END SELECT

        call assign_source_params(sp,ncgs,tcgs,ncgsnth)
    end subroutine convert_fluid_vars_arr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !get_fluidvars functions for specific models
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_thindisk_fluidvars(x0,k0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0, k0
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        real, dimension(size(x0)) :: T,omega
        real, dimension(size(x0),10) :: metric
        call thindisk_vals(real(x0%data(2)),real(x0%data(3)),a,T,omega)
        f%rho=T
        f%u%data(2)=0.; f%u%data(3)=0.
        f%b%data(1)=0.; f%b%data(2)=0.; f%b%data(3)=0.; f%b%data(4)=0.

        metric=kerr_metric(real(x0%data(2)),real(x0%data(3)),a)
!        f%b%data(2)=cos(x0%data(4))/metric(:,5)
!        f%b%data(4)=-sin(x0%data(4))/metric(:,10)

        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))
!       TEST SET OMEGA = 0 to look at pol
!       omega(:)=0d0
        f%u%data(1)=sqrt(-1./(metric(:,1)+2d0*metric(:,4)*omega+metric(:,10)* &
         omega*omega))
        f%u%data(4)=omega*f%u%data(1)
        ! Assign normal vector as magnetic field for comoving_ortho:
        f%b = calc_polvec(x0%data(2),cos(x0%data(3)),k0,dble(a),asin(1d0))
!       f%b%data(1)=-f%b%data(1)
!       f%b%data(2)=-f%b%data(2)
!       f%b%data(3)=-f%b%data(3)
!       f%b%data(4)=-f%b%data(4)
    end subroutine get_thindisk_fluidvars

    subroutine get_numdisk_fluidvars(x0,k0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0, k0
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        real, dimension(size(x0)) :: T,omega,phi
        real, dimension(size(x0),10) :: metric
        phi=x0%data(4); phi=phi+12.*acos(-1.)
        phi=mod(phi,(2.*acos(-1.)))
        call numdisk_vals(real(x0%data(2)),phi,a,T,omega)

        f%rho=T
        f%u%data(2)=0d0; f%u%data(3)=0d0
        f%b%data(1)=0d0; f%b%data(2)=0d0; f%b%data(3)=0d0; f%b%data(4)=0d0

        metric=kerr_metric(real(x0%data(2)),real(x0%data(3)),a)
        f%b%data(2)=cos(x0%data(4))/metric(:,5)
        f%b%data(4)=-sin(x0%data(4))/metric(:,10)

        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))
        f%u%data(1)=sqrt(-1d0/(metric(:,1)+2d0*metric(:,4)*omega+metric(:,10)*omega*omega))
        f%u%data(4)=omega*f%u%data(1)
        f%b = calc_polvec(x0%data(2),cos(x0%data(3)),k0,dble(a),0d0)
    end subroutine get_numdisk_fluidvars

    subroutine get_phat_fluidvars(x0,k0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0, k0
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        real, dimension(size(x0)) :: omega
        real, dimension(size(x0),10) :: metric
        real, dimension(size(x0),size(freq_tab)) :: fnu

        call phatdisk_vals(real(x0%data(2)),a,fnu,omega)
        f%fnu=fnu
        f%u%data(2)=0.; f%u%data(3)=0.
        f%b%data(1)=0.; f%b%data(2)=0.; f%b%data(3)=0.; f%b%data(4)=0.

        metric=kerr_metric(real(x0%data(2)),real(x0%data(3)),a)
        call assign_metric(f%u,transpose(metric))
!        call assign_metric(f%b,transpose(metric))

        f%u%data(1)=sqrt(-1d0/(metric(:,1)+2d0*metric(:,4)*omega+metric(:,10)*omega*omega))
        f%u%data(4)=omega*f%u%data(1)
        f%b = calc_polvec(x0%data(2),cos(x0%data(3)),k0,dble(a),0.d0)
    end subroutine get_phat_fluidvars

    subroutine get_ffjet_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f

        ! Computes properties of jet solution from Broderick & Loeb (2009)
        ! JAD 4/23/2010, fortran 3/30/2011
        call ffjet_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
    end subroutine get_ffjet_fluidvars

    subroutine get_rrjet_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f
        ! Computes properties of jet model from Anantua et al (2019)
        ! AC 4/2020
        call rrjet_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
    end subroutine get_rrjet_fluidvars
       
    subroutine get_harm_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f

        call harm_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
    end subroutine get_harm_fluidvars

    subroutine get_harm3d_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f

        call harm3d_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
    end subroutine get_harm3d_fluidvars

    subroutine get_harmpi_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f

        call harmpi_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag,f%kela, &
                         f%kelb,f%kelc,f%keld)
    end subroutine get_harmpi_fluidvars

    subroutine get_iharm_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f

        call iharm_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
    end subroutine get_iharm_fluidvars

    subroutine get_thickdisk_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f

        call thickdisk_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
    end subroutine get_thickdisk_fluidvars

    subroutine get_mb09_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f

        call mb09_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
    end subroutine get_mb09_fluidvars

    subroutine get_koral_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f !AC TODO merge with 3D

        call koral_vals(x0,a,f%rho,f%te,f%b,f%u,f%bmag,f%nnth) !AC replaced p-->te
    end subroutine get_koral_fluidvars

    subroutine get_koral3d_fluidvars(x0,a,f,type)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        integer, intent(in) :: type
        type (fluid), intent(inout) :: f

        call koral3d_vals(x0,a,f%rho,f%te, f%ti, f%b,f%u,f%bmag,f%Be,f%nnth,type) !AC replaced p-->te
    end subroutine get_koral3d_fluidvars

    subroutine get_koralh5_fluidvars(x0,a,f,type)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        integer, intent(in) :: type
        type (fluid), intent(inout) :: f
        integer :: electrons

        call koralh5_vals(x0,a,f%rho,f%p, f%te, f%ti, f%b,f%u,f%bmag,f%Be,type,has_electrons)
        has_electrons = electrons

    end subroutine get_koralh5_fluidvars


    subroutine get_sphacc_fluidvars(x0,f)
        type (four_vector), intent(in), dimension(:) :: x0
        type (fluid), intent(inout) :: f
        real, dimension(size(x0)) :: u,grr,g00,n,B,T,ur

        u=1./x0%data(2)
        call sphacc_vals(u,n,B,T,ur) ! Equipartition B field


        g00=-(1.-2.*u)
        grr=-1d0/g00
        f%u%data(2)=-ur
        f%u%data(3)=0d0
        f%u%data(4)=0d0

        f%u%data(1)=sqrt((-grr*f%u%data(2)*f%u%data(2)-1)/g00)


        ! use u dot b = 0, b dot b = B^2 to get components:
        f%b%data(1)=sqrt(f%u%data(2)**2*grr*B**2/ &
                         (f%u%data(2)**2*g00*grr+f%u%data(1)**2*g00*g00))
        f%b%data(2)=-sqrt(B**2/grr-f%b%data(1)**2*g00/grr)
        f%b%data(3)=0d0
        f%b%data(4)=0d0

        f%rho=n
        f%p=T
        f%bmag=B

    end subroutine get_sphacc_fluidvars

    subroutine get_hotspot_fluidvars(x0,k0,a,f)
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        type (four_vector), intent(in), dimension(:) :: x0,k0
        type (four_Vector), dimension(size(x0)) :: x,kin,xout
        real(kind=8), dimension(size(x0)) :: n

        x=x0; kin=k0
        ! transform geodesic phi, t
        x%data(4) = -1*acos(0.) - x%data(4)
        x%data(1) = -1.*x%data(1)

        call hotspot_vals(x,kin,dble(a),n,f%b,f%u,xout)
        f%rho = dble(n)
        f%bmag = sqrt(f%b*f%b)
    end subroutine get_hotspot_fluidvars

    subroutine get_schnittman_hotspot_fluidvars(x0,a,f)
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        type (four_vector), intent(in), dimension(:) :: x0
        type (four_Vector), dimension(size(x0)) :: x, xout
        real, dimension(size(x0),10) :: metric
        real, dimension(size(x0)) :: n,omega,gfac,omt,ut,lc,hc,om,ar,d,safe
        real :: rms

        call schnittman_hotspot_vals(x0,a,n)

        omega=1./(x0%data(2)**(3./2.)+a)
        f%rho = dble(n)
        f%u%data(2)=0d0; f%u%data(3)=0d0
        f%b%data(1)=0d0; f%b%data(2)=0d0; f%b%data(3)=0d0; f%b%data(4)=0d0
        metric=kerr_metric(real(x0%data(2)),real(x0%data(3)),a)
        f%u%data(1)=sqrt(-1d0/(metric(:,1)+2d0*metric(:,4)*omega+metric(:,10)*omega*omega))
        f%u%data(4)=omega*f%u%data(1)

        rms=calc_rms(a)
        d=x0%data(2)*x0%data(2)-2d0*x0%data(2)+a*a
        lc=(rms*rms-2.*a*sqrt(rms)+a*a)/(rms**1.5d0-2d0*sqrt(rms)+a)
        hc=(2.*x0%data(2)-a*lc)/d
        ar=(x0%data(2)*x0%data(2)+a*a)**2d0-a*a*d*sin(x0%data(3))**2d0
        om=2.*a*x0%data(2)/ar
        
        where(x0%data(2).gt.rms)
           omt=max(1d0/(x0%data(2)**(3d0/2d0)+a),om)
        elsewhere
           omt=max((lc+a*hc)/(x0%data(2)*x0%data(2)+2d0*x0%data(2)*(1d0+hc)),om)
        end where

        ut=sqrt(-1./(metric(:,1)+2.*metric(:,4)*omt+metric(:,10)*omt*omt))
        safe=metric(:,1)+2d0*metric(:,4)*omega+metric(:,10)*omega*omega
        f%u%data(1)=merge(f%u%data(1),dble(ut),safe.lt.0d0)
        f%u%data(4)=merge(f%u%data(4),omt*f%u%data(1),safe.lt.0d0)

        ! Toroidal magnetic field:
        gfac=1d0/sqrt((metric(:,10)*metric(:,1)-metric(:,4)*metric(:,4))* & 
           (metric(:,10)*f%u%data(4)*f%u%data(4)+f%u%data(1)* & 
           (2d0*metric(:,4)*f%u%data(4)+metric(:,1)*f%u%data(1))))
        f%b%data(1)=gfac*abs(metric(:,10)*f%u%data(4)+metric(:,4)*f%u%data(1))
        f%b%data(4)=-sign(1d0,metric(:,10)*f%u%data(4)+metric(:,4)*f%u%data(1)) &
           *(f%u%data(1)*metric(:,1)+metric(:,4)*f%u%data(4))*gfac

        call assign_metric(f%b,transpose(metric))
        call assign_metric(f%u,transpose(metric))


        f%bmag = sqrt(f%b*f%b)
    end subroutine get_schnittman_hotspot_fluidvars

    subroutine get_sariaf_fluidvars(x0,a,f)

        type (four_vector), intent(in), dimension(:) :: x0
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        real, dimension(size(x0)) :: u,ctheta
        real, dimension(size(x0)) :: riaf_vr,riaf_vth,riaf_omega, bmag,n,t,nnth
        real, dimension(size(x0)) :: ub,aleph,gfac,bb
        real, dimension(size(x0)) :: rr, rho2,psi4,stheta
        real :: rms

        real :: lambdae,game 
        real, dimension(size(x0)) :: delta,hhh 
        real(8), dimension(size(x0),10) :: metric
        real(8), dimension(size(x0)) :: hc,lc,d,ar,om,omtest,zero,vth,vphi,vr,u0
        integer :: bl06
        integer :: i
        type (four_vector), dimension(size(x0)) :: fu

        rr = x0%data(2)
        rms = calc_rms(a)
        zero = 0d0

        lambdae = (rms**2. - 2.*a*sqrt(rms) + a**2.)/(rms**(3./2.)-2.*sqrt(rms)+a)
        game = sqrt(1.-2./3./rms)
        delta = rr*rr - 2.*rr + a*a
        hhh = (2.*rr-a*lambdae)/delta

        u=1./rr
        ctheta =cos(x0%data(3))
        stheta = sin(x0%data(3))
        rho2 = rr**2. + a**2. * (ctheta)**2.
        psi4 = 2.*rr/rho2

        call sariaf_vals(a,ctheta,u,n,t,bmag,riaf_vr,riaf_vth,riaf_omega,nnth,bl06)

        ! construct four-velocity from three-velocity
        metric=kerr_metric(real(rr),real(x0%data(3)),a)
        lc = (rms**2d0-2d0*a*sqrt(rms)+a**2d0)/(rms**1.5d0-2d0*sqrt(rms)+a)
        d = rr**2d0-2d0*rr+a**2d0
        ar = (rr**2d0+a**2d0)**2d0-a**2d0*d*sin(x0%data(3))
        om = 2d0*a*rr/ar
        hc = (2d0*rr-a*lc)/d

        u0 = calc_u0(metric,dble(riaf_vr),dble(riaf_vth),dble(riaf_omega))
        fu = rms_vel(dble(a),x0%data(3),x0%data(2))

        where(rr.lt.rms)
           f%u%data(1) = fu%data(1)
           f%u%data(2) = fu%data(2)
           f%u%data(3) = fu%data(3)
           f%u%data(4) = fu%data(4)
        elsewhere
           f%u%data(1) = u0
           f%u%data(2) = riaf_vr*f%u%data(1)
           f%u%data(3) = riaf_vth*f%u%data(1)
           f%u%data(4) = riaf_omega*f%u%data(1)
        end where

        ! Radial magnetic field
        !aleph = -1.*(metric(:,4)*f%u%data(1)+metric(:,10)*f%u%data(4)) &
        !     /(metric(:,1)*f%u%data(1)+metric(:,4) * f%u%data(4))
        !bb = metric(:,1)*aleph*aleph + metric(:,10) + 2.*metric(:,4)*aleph !I hope this is never negative
        !f%b%data(4) = bmag/sqrt(bb) !what sign?
        !f%b%data(3) = 0d0
        !f%b%data(2) = 0d0
        !f%b%data(1) = aleph * f%b%data(4)s

        ! Toroidal magnetic field (!AC -- copied from fluid_model_hotspot.f90)
        gfac=1d0/sqrt((metric(:,10)*metric(:,1)-metric(:,4)*metric(:,4))* & 
           (metric(:,10)*f%u%data(4)*f%u%data(4)+f%u%data(1)* & 
           (2d0*metric(:,4)*f%u%data(4)+metric(:,1)*f%u%data(1))))
        f%b%data(1)=bmag*gfac*abs(metric(:,10)*f%u%data(4)+metric(:,4)*f%u%data(1))
        f%b%data(4)=-bmag*sign(1d0,metric(:,10)*f%u%data(4)+metric(:,4)*f%u%data(1)) &
                    *(f%u%data(1)*metric(:,1)+metric(:,4)*f%u%data(4))*gfac
        f%b%data(2)=0d0
        f%b%data(3)=0d0

        call assign_metric(f%b,transpose(metric))
        call assign_metric(f%u,transpose(metric))
        
        f%bmag = sqrt(f%b*f%b)
        f%rho = n
        f%p = t
        f%bmag = bmag
        f%rho2 = nnth

        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))

    end subroutine get_sariaf_fluidvars

    subroutine get_toy_fluidvars(x0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0
        type (fluid), intent(inout) :: f
        real(kind=8), intent(in) :: a
        type (four_vector), dimension(size(x0)) :: ulower
        real(kind=8), dimension(size(x0)) :: ctheta
        real(kind=8), dimension(size(x0)) :: bmag,n,l
        real(kind=8), dimension(size(x0)) :: ubar,aleph,bb
        real(kind=8), dimension(size(x0)) :: rr
        real(kind=8) :: rms
        real(kind=8), dimension(size(x0),10) :: metric,con_metric
        type (four_vector), dimension(size(x0)) :: fu

        rr = x0%data(2)
        ctheta =cos(x0%data(3))

        call toy_vals(a,ctheta,rr,n,l)
        
        ! construct four-velocity from three-velocity
        metric=kerr_metric(rr,x0%data(3),a)
        con_metric=kerr_metric(rr,x0%data(3),a,1)

        aleph = -1.*(metric(:,4)*f%u%data(1)+metric(:,10)*f%u%data(4)) &
             /(metric(:,1)*f%u%data(1)+metric(:,4) * f%u%data(4))
        ubar = sqrt(-1./(con_metric(:,1) + l*l*con_metric(:,10) - &
             2*l*con_metric(:,4)))
        ulower%data(1)=-ubar; ulower%data(2)=0d0; ulower%data(3)=0d0
        ulower%data(4)=l*ubar
        
        ! u^t = g^t\mu u_\mu = g^t\phi u_\phi + g^tt u_t
        f%u%data(2)=0d0; f%u%data(3)=0d0
        f%u%data(1)=con_metric(:,4)*ulower%data(4)+con_metric(:,1)*ulower%data(1)
        f%u%data(4)=con_metric(:,4)*ulower%data(1)+con_metric(:,10)*ulower%data(4)

        ! radial magnetic field
        bmag=1d0
        f%b%data(4) = bmag/sqrt(bb)
        f%b%data(3) = 0d0
        f%b%data(2) = 0d0
        f%b%data(1) = aleph * f%b%data(4)

        f%rho = n
        f%bmag = bmag

        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))

    end subroutine get_toy_fluidvars

    subroutine get_powerlaw_fluidvars(x0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        real, dimension(size(x0)) :: u,ctheta
        real, dimension(size(x0)) :: vr,vth,omega,bmag,n,t,nnth
        real, dimension(size(x0)) :: ub,aleph,bb
        real, dimension(size(x0)) :: rr, rho2,psi4,stheta
        real :: rms
        real :: lambdae,game 
        real, dimension(size(x0)) :: delta,hhh 
        real(8), dimension(size(x0),10) :: metric
        real(8), dimension(size(x0)) :: hc,lc,d,ar,om,omtest,zero,u0
        type (four_vector), dimension(size(x0)) :: fu

        rr = x0%data(2)
        rms = calc_rms(a)
        zero = 0d0

        lambdae = (rms**2. - 2.*a*sqrt(rms) + a**2.)/(rms**(3./2.)-2.*sqrt(rms)+a)
        game = sqrt(1.-2./3./rms)
        delta = rr*rr - 2.*rr + a*a
        hhh = (2.*rr-a*lambdae)/delta

        u=1./rr
        ctheta =cos(x0%data(3))
        stheta = sin(x0%data(3))
        rho2 = rr**2. + a**2. * (ctheta)**2.
        psi4 = 2.*rr/rho2

        call powerlaw_vals(a,ctheta,u,n,t,bmag,vr,vth,omega,nnth)

        ! construct four-velocity from three-velocity
        metric = kerr_metric(x0%data(2),x0%data(3),dble(a))
        u0 = calc_u0(metric,dble(vr),dble(vth),dble(omega))
!        fu = rms_vel(dble(a),x0%data(3),x0%data(2))
!        where(rr.lt.rms)
!           f%u%data(1) = fu%data(1)
!           f%u%data(2) = fu%data(2)
!           f%u%data(3) = fu%data(3)
!           f%u%data(4) = fu%data(4)
!        elsewhere
!           f%u%data(1) = u0
!           f%u%data(2) = vr*f%u%data(1)
!           f%u%data(3) = vth*f%u%data(1)
!           f%u%data(4) = omega*f%u%data(1)
!        end where

        f%u%data(1) = u0
        f%u%data(2) = vr*u0
        f%u%data(3) = vth*u0
        f%u%data(4) = omega*u0

        aleph = -1.*(metric(:,4)*f%u%data(1)+metric(:,10)*f%u%data(4)) &
             /(metric(:,1)*f%u%data(1)+metric(:,4) * f%u%data(4))

        bb = metric(:,1)*aleph*aleph + metric(:,10) + 2.*metric(:,4)*aleph ! hope this is never negative

        f%b%data(4) = bmag/sqrt(bb) !what sign?
        f%b%data(3) = 0d0
        f%b%data(2) = 0d0
        f%b%data(1) = aleph * f%b%data(4)
        f%rho = n
        f%p = t
        f%bmag = bmag
        f%rho2 = nnth

        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))

    end subroutine get_powerlaw_fluidvars

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !convert_fluidvars functions for specific models
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! scale code to cgs units given mbh and simulation and cgs mdot values
    subroutine scale_sim_units(mbh,mdotcgs,mdot,rho,p,bmag,ncgs,bcgs,tempcgs)
        real(kind=8), intent(in) :: mbh,mdot,mdotcgs
        real(kind=8) :: lcgs,tcgs
        real, dimension(:) :: rho,p,bmag
        real(kind=8), dimension(size(rho)) :: rhocgs,pcgs
        real(kind=8), dimension(size(rho)), intent(inout) :: ncgs,bcgs,tempcgs

        ! JAD 11/26/2012 adapted from IDL code
        ! Black hole mass sets time and length scales:
        lcgs=GC*mbh*msun/c**2; tcgs=lcgs/c
        rhocgs=mdotcgs/mdot/lcgs**3*tcgs*rho; ncgs=rhocgs/mp

        ! Use this to convert pressure:
        pcgs=p*rhocgs/rho*c**2d0

        ! Ideal gas temperature for single fluid (i.e., no separate e-/p):
        tempcgs=pcgs/ncgs/k

        ! finally, bfield conversion is just square root of this:
        bcgs=bmag*sqrt(rhocgs/rho)*c

        ! Convert HL units to Gaussian cgs:
        bcgs=bcgs*sqrt(4d0*pi)

    end subroutine scale_sim_units

    ! zero emissivity below a sigmacut threshold
    subroutine andrew_sigcut(bcgs,rhocgs,tempcgs,ncgs,sigcut)
        real(kind=8), dimension(:), intent(inout) :: bcgs,rhocgs,tempcgs,ncgs
        real(kind=8), intent(in) :: sigcut
        real(kind=8), dimension(size(bcgs)) :: sigmacgs
        sigmacgs = (bcgs*bcgs) / (rhocgs*8.988e20*4*pi)
        if(any(sigmacgs.ge.sigcut)) then

          ! JD changing values to see if we can avoid Koral NaNs
          where(sigmacgs.ge.sigcut)
             rhocgs = 0d0
             ncgs = 0d0
             tempcgs = 1d9
             bcgs = 1d-4
          end where
        end if
    end subroutine andrew_sigcut

    ! EHT theory notes formulae based on Moscibrodzka+2016
    ! assumes inputs are in cgs units
    subroutine charles_e(rho,p,u,b,beta_trans,rlow,rhigh,tcgs)
        real(kind=8), intent(in), dimension(:) :: rho,p,u,b
        real(kind=8), intent(in) :: rlow,rhigh,beta_trans
        real(kind=8), intent(inout), dimension(size(rho)) :: tcgs
        real(kind=8), dimension(size(rho)) :: beta,b2,trat

        ! here p = T_p + T_e and u = T_p + 2 T_e (assumes gammae=4/3, gammai=5/3)
        ! b is in HL
        beta=2d0*rho*k*p/mp/(b*b)
        b2=beta*beta

        where(b.gt.0d0)
         trat=rhigh*b2/(1d0+b2)+rlow/(1d0+b2)
        elsewhere
         trat=rhigh
        end where

        ! notes formula 9 with u = 2u/3nk = Tp+2Te from Koral
        tcgs=u/(2d0+trat)
    end subroutine charles_e

    ! as written in equation 1 of Moscibrodzka+2016 assuming cgs units
    subroutine monika_e_orig(rho,p,b,beta_trans,rlow,rhigh,trat)
        real(kind=4), intent(in), dimension(:) :: rho,p,b
        real(kind=8), intent(in) :: rlow,rhigh,beta_trans
        real(kind=8), intent(inout), dimension(size(rho)) :: trat
        real(kind=8), dimension(size(rho)) :: beta,b2

        ! here p is the total *temperature* Tp+2Te
        ! 2 is here to account for adiabatic indices
        ! b is in HL
        beta=2d0*rho*k*p/mp/(b*b)
        b2=beta*beta
        where(b.gt.0d0)
         trat=rhigh*b2/(1d0+b2)+rlow/(1d0+b2)
        elsewhere
         trat=rhigh
        end where

    end subroutine monika_e_orig

    ! as written in equation 1 of Moscibrodzka+2016 assuming code units
    subroutine monika_e(rho,p,b,beta_trans,rlow,rhigh,trat)
        real(kind=4), intent(in), dimension(:) :: rho,p,b
        real(kind=8), intent(in) :: rlow,rhigh,beta_trans
        real(kind=8), intent(inout), dimension(size(rho)) :: trat
        real(kind=8), dimension(size(rho)) :: beta,b2

        beta=p/(b*b)/0.5d0
        b2=beta*beta

        where(b.gt.0d0)
           trat=rhigh*b2/(1d0+b2)+rlow/(1d0+b2)
        elsewhere
           trat=rhigh
        end where

    end subroutine monika_e

    subroutine ressler_e(rho,kel,tcgs)
        real(kind=4), intent(in), dimension(:) :: kel,rho
        real(kind=8), intent(inout), dimension(size(rho)) :: tcgs
        real(kind=8), dimension(size(rho)) :: thetae
        real(kind=8) :: gamma
        gamma=4d0/3d0
        thetae = mp/m*kel*rho**(gamma-1.)
        tcgs = thetae*m*c2/k
    end subroutine ressler_e

    subroutine werner_e(rho,bmag,deltae)
        real(kind=4), intent(in), dimension(:) :: rho,bmag
        real(kind=8), intent(inout), dimension(size(rho)) :: deltae
        ! formula from figure 20 of Werner+2018 with sigma_i = b^2 / rho
        deltae=1d0/4d0+1d0/4d0*sqrt(bmag**2d0/rho/5d0/(2d0+bmag**2d0/rho/5d0))
    end subroutine werner_e

    ! non-thermal e- where jet energy density is high (e.g. Broderick & McKinney 2010, Dexter+2012)
    subroutine nonthermale_b2(alpha,gmin,p1,p2,bmagrho,bcgs,ncgsnth)
        real(kind=8), intent(in) :: alpha,gmin,p1,p2
        real(kind=8), intent(in), dimension(:) :: bmagrho,bcgs
        real(kind=8), intent(inout), dimension(:) :: ncgsnth
        where(bmagrho.gt.1d0)
           ncgsnth=alpha*bcgs**2d0/8d0/pi/gmin*(p1-2d0)/(p1-1d0)/8.2d-7
        elsewhere
           ncgsnth=0d0
        end where
    end subroutine nonthermale_b2

    subroutine convert_fluidvars_thickdisk(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: mdot,beta_trans
        mdot=0.0013; beta_trans=1d0
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs, &
             bcgs,tempcgs)
        call monika_e(f%rho,f%p,f%bmag,beta_trans, & 
                      1d0/sp%muval-1d0, sp%gminval*(1d0/sp%muval-1d0),trat)
        tempcgs = tempcgs/(1d0+trat)
        call nonthermale_b2(sp%jetalphaval,sp%gminval,sp%p1,sp%p2, &
             f%bmag**2d0/f%rho,bcgs,ncgsnth)
    end subroutine convert_fluidvars_thickdisk

    subroutine convert_fluidvars_mb09(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: mdot,beta_trans
        mdot=.0013; beta_trans=1d0

        ! code --> cgs units
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs,bcgs,tempcgs)

        ! thermal part
        call monika_e(f%rho,f%p,f%bmag,beta_trans, &
                      1d0/sp%muval-1d0, sp%gminval*(1d0/sp%muval-1d0),trat)
        tempcgs = tempcgs/(1d0+trat)

        ! non-thermal e- by hand
        ncgsnth=ncgs
    end subroutine convert_fluidvars_mb09

    subroutine convert_fluidvars_harm(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: mdot,beta_trans
        mdot=.003; beta_trans=1d0

        ! code --> cgs units
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs,bcgs,tempcgs)

        ! thermal part
        ! Moscibrodzka+2016 e- model with rlow = T_p / T_e from muval, rhigh = gmin*rlow
        ! reduces to T_p / T_e = const when gmin = 1 (should change input used for this)
        ! CHANGED TO USE MU INSTEAD OF TRAT to allow mu > 1/2 models to make sense 3/21/2017
        call monika_e(f%rho,f%p,f%bmag,beta_trans, &
                      sp%muval,sp%muval/sp%gminval,trat)
        tempcgs = tempcgs*trat
        ncgsnth=ncgs
    end subroutine convert_fluidvars_harm

    subroutine convert_fluidvars_harm3d(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat,rhocgs
        real(kind=8) :: mdot,beta_trans
        mdot=GC*sp%mbh*msun/c**3; beta_trans=1d0
        ! code --> cgs units
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs,bcgs,tempcgs)

        ! thermal part (rlow = 1/mu-1, rhigh=gmin*rlow)
        call monika_e(f%rho,f%p,f%bmag,beta_trans, &
                      1d0/sp%muval-1d0, sp%gminval*(1d0/sp%muval-1d0),trat)
        tempcgs = tempcgs/(1d0+trat)

        ! nonthermal part
        call nonthermale_b2(sp%jetalphaval,sp%gminval,sp%p1,sp%p2, f%bmag**2d0/f%rho,bcgs,ncgsnth)

        ! adding option to remove emission from highly magnetized regions \sigma > sigcut
        rhocgs=ncgs*mp
        call andrew_sigcut(bcgs,rhocgs,tempcgs,ncgs,dble(sp%sigcut))
    end subroutine convert_fluidvars_harm3d

    subroutine convert_fluidvars_iharm(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: mdot,beta_trans
        mdot=GC*sp%mbh*msun/c**3; beta_trans=1d0
        ! code --> cgs units
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs, bcgs,tempcgs)

        ! thermal part(rlow = 1/mu-1, rhigh=gmin*rlow)
        call monika_e(f%rho,f%p,f%bmag,beta_trans, &
                      1d0/sp%muval-1d0, sp%gminval*(1d0/sp%muval-1d0),trat)
        tempcgs = tempcgs/(1d0+trat)
    
        ! nonthermal part
        call nonthermale_b2(sp%jetalphaval,sp%gminval,sp%p1,sp%p2,f%bmag**2d0/f%rho,bcgs,ncgsnth)
    end subroutine convert_fluidvars_iharm
        
    subroutine convert_fluidvars_harmpi(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat,rhocgs
        real(kind=8) :: mdot,beta_trans
        mdot=GC*sp%mbh*msun/c**3; beta_trans=1d0

        ! code --> cgs units
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs, bcgs,tempcgs)
        rhocgs=ncgs*mp

        ! thermal part
        ! defaults to monika_e (rlow = 1/mu-1, rhigh=gmin*rlow)
        ! to use ressler_e call grtrans with gmin = -1
        if(sp%gminval.ge.1d0) then
           call monika_e(f%rho,f%p,f%bmag,beta_trans,&
                         1d0/sp%muval-1d0,sp%gminval*(1d0/sp%muval-1d0),trat)
           tempcgs = tempcgs/(1d0+trat)
        else if(sp%gminval.lt.0d0) then
           ! gmin=-1,2,3,4 to use kel4abcd
           ! not yet reading capability for ktot, kel4, kel4e, kel5, keldis
           if(sp%gminval.eq.-1d0) then
              call ressler_e(f%rho,f%kela,tempcgs)
           else if(sp%gminval.eq.-2d0) then
              call ressler_e(f%rho,f%kelb,tempcgs)
           else if(sp%gminval.eq.-3d0) then
              call ressler_e(f%rho,f%kelc,tempcgs)
           else
              call ressler_e(f%rho,f%keld,tempcgs)
           endif
        else
           ! this is the Werner+2018 heating model with an additional muval scaling
           ! if you want *just* Werner+ alone then use muval=1
           call werner_e(f%rho,f%bmag,trat)
           tempcgs = sp%muval*trat*tempcgs

        end if

        ! nonthermal part
        call nonthermale_b2(sp%jetalphaval,sp%gminval,sp%p1,sp%p2,f%bmag**2d0/f%rho,bcgs,ncgsnth)

        ! sigma cut
        call andrew_sigcut(bcgs,rhocgs,tempcgs,ncgs,dble(sp%sigcut))

    end subroutine convert_fluidvars_harmpi

    ! AC: KORAL quantities should ALREADY be in CGS 
    ! AC: be careful with HL vs Gaussian Bfield!
    subroutine convert_fluidvars_koral(f,ncgs,ncgsnth,bcgs,tempcgs,nnthcgs,sp,cuttype,vartype)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        integer, intent(in) :: cuttype, vartype
        real(kind=8) :: sigcut, mbh,lcgs,dcgs
        real(kind=8), dimension(size(f%rho)) :: sigmacgs
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho), f%nrelbin), intent(out) :: nnthcgs        
        real(kind=8), dimension(size(f%rho)) :: rhocgs, pcgs, tecgs, ticgs, bhl
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: beta_trans

        sigcut=sp%sigcut
        mbh=sp%mbh

        ! convert to cgs if necessary
        if(vartype.eq.1) then ! convert to cgs -- fix mdot

          lcgs=GC*mbh*msun/c**2
          dcgs=mbh*msun/(lcgs**3)

          rhocgs=dcgs*f%rho; 
          ncgs=rhocgs/mp
          pcgs=f%p*dcgs*(c**2d0)

          if(has_electrons.eq.1) then ! has_electrons is a global set in get_fluidvars_koralh5
            tecgs = f%te
            ticgs = f%ti
          else
            tecgs =pcgs/(2*ncgs)/k ! single electron-ion fluid, assumed equal temperature her
          endif

          bhl=f%bmag*sqrt(dcgs)*c
          bcgs=bhl*sqrt(4d0*pi) !HL to Gaussian

          nnthcgs=f%nnth ! binned non-thermal e- (TODO doesn't exist in h5 yet, so no units)
          ncgsnth=0. !AC TODO: add other koral nonthermal prescriptions

        else ! already in cgs

          rhocgs=f%rho
          ncgs=rhocgs/mp !TODO AC: handle non-hydrogen plasma
          tecgs=f%te ! AC: finally changed p-->te array
          ticgs = f%ti

          !convert HL B-field to gaussian
          bhl = f%bmag
          bcgs= bhl*sqrt(4*pi)

          nnthcgs=f%nnth ! binned non-thermal e-, already in cgs if present
          ncgsnth=0. !AC TODO: add other koral nonthermal prescriptions
        endif

        ! JD: allow scaling w/ Mdot
!       write(6,*) 'convert koral nfac: ', sp%nfac,sp%sigcut,sp%gminval
        bcgs=bcgs*sqrt(sp%nfac)
        ncgs=ncgs*sp%nfac; 
        rhocgs=rhocgs*sp%nfac;

        ! AC: by default the p variable stores electron temperature for KORAL
        ! JD: changed to allow postprocessing Monika model if gmin>1 (rhigh=gmin,rlow=1)
        if(sp%gminval.ge.1d0) then
           beta_trans = 1d0

           ! f%Be is confusingly the proton temperature here (AC: changed to electron temperature!)
           ! JD: CHANGING koral Rhigh method to align with EHT theory notes
!           call charles_e(f%rho,f%p+f%Be,2.*f%p+f%Be,f%bmag,beta_trans,1d0,sp%gminval,tempcgs)
           call charles_e(rhocgs,tecgs+ticgs,2.*tecgs+ticgs,bhl,beta_trans,1d0,sp%gminval,tempcgs)
        else 
!           tempcgs=f%p 
           tempcgs=f%te ! AC: finally changed p-->te array
        endif

        ! apply sigma cut
        call andrew_sigcut(bcgs,rhocgs,tempcgs,ncgs,sigcut)

        !ANDREW -- zero out density to zero out emissivity in disk or jet
        if(cuttype.eq.1) then  !zero out jet
            !AC --  old -- sigma cutoff defines jet
            !where(sigmacgs.ge.1)
            !   rhocgs = 0.
            !   tempcgs = 10.
            !   bcgs = 0.
            !end where
            !ANDREW --  new -- Be cutoff defines jet
            if(any(f%Be.ge.0.05)) then 
                where((f%Be.ge.0.05))
                   rhocgs = 0.
                   tempcgs = 10.
                   bcgs = 0.
                end where
            end if
       
        elseif((cuttype.eq.2).or.(cuttype.eq.3)) then !zero out disk
            !ANDREW --  old -- sigma cutoff defines jet
            !where(sigmacgs.le.1)
            !   rhocgs = 0.
            !   tempcgs = 10.
            !   bcgs = 0.
            !end where
            !ANDREW --  new -- Be cutoff defines jet
            if(any((f%Be.le.0.05).and.(sigmacgs.le.1))) then 
                where((f%Be.le.0.05).and.(sigmacgs.le.1))
                   rhocgs = 0.
                   tempcgs = 10.
                   bcgs = 0.
                end where
            end if
        endif
        
        
        ! warnings
        if(any(isnan(tempcgs))) then
            write(6,*) 'KORAL NAN temp: '
        endif
        if(any(isnan(ncgs))) then
            write(6,*) 'KORAL NAN n: '
        endif
        if(any(isnan(bcgs))) then
            write(6,*) 'KORAL NAN b: '
        endif
    end subroutine convert_fluidvars_koral

    subroutine convert_fluidvars_ffjet(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), &
          intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp

        ncgsnth=f%rho*sp%nfac; bcgs=f%bmag*sp%bfac; ncgs=0.; tcgs=0.
    end subroutine convert_fluidvars_ffjet

    subroutine convert_fluidvars_rrjet(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), &
             intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        real ::  prefac, prefacOLD, conv1, conv2
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)) :: pcgs, rhocgs
        
        ! pressure and magnetic field already in cgs
        bcgs = f%bmag
        pcgs = f%p
        
        ! nonthemal number density - get from total pressure
        prefac = 3.*(sp%p2-2.) / (sp%p2-1.)
        prefac = prefac * (sp%gminval**(1.-sp%p2) - sp%gmax**(1.-sp%p2))
        prefac = prefac / (sp%gminval**(2.-sp%p2) - sp%gmax**(2.-sp%p2)) 
        prefac=prefac/(m*c2)

        !from total pressure - ignoring gamma max
        !prefac = 3.*(sp%p2-2.) / ((sp%p2-1.)*sp%gminval) / (m*c2)

        ! nonthermal number density - from partial pressure, p=2 ONLY
        !prefac = 3.* (sp%gminval**(1.-sp%p2) - sp%gmax**(1.-sp%p2)) / ((sp%p2-1.) * (m*c2))
        
        ! calculate number density and mass density
        rhocgs = prefac * pcgs * mp
        ncgsnth= prefac * pcgs 
        
        ! the temperature and thermal energy density are both zero
        ncgs=ncgsnth;
        tcgs=0.

    end subroutine convert_fluidvars_rrjet

    subroutine convert_fluidvars_sphacc(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(inout) :: ncgs,ncgsnth,bcgs,tcgs

        ncgs=f%rho; bcgs=f%bmag; tcgs=f%p; ncgsnth=0.
    end subroutine convert_fluidvars_sphacc

    subroutine convert_fluidvars_hotspot(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), &
          intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp
    
        ncgs=f%rho; bcgs=f%bmag; ncgsnth=f%rho
    end subroutine convert_fluidvars_hotspot

    subroutine convert_fluidvars_schnittman_hotspot(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), &
          intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp
      
        ncgs=f%rho; bcgs=1d0
    end subroutine convert_fluidvars_schnittman_hotspot

    subroutine convert_fluidvars_thindisk(f,tcgs,ncgs)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), intent(out) :: tcgs,ncgs

        tcgs=f%rho; ncgs=1.
    end subroutine convert_fluidvars_thindisk

    subroutine convert_fluidvars_numdisk(f,tcgs,ncgs)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), intent(out) :: tcgs,ncgs

        tcgs=f%rho; ncgs=1.
    end subroutine convert_fluidvars_numdisk

    subroutine convert_fluidvars_phatdisk(f,fnu,nu)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(:,:), allocatable, intent(out) :: fnu
        real(kind=8), dimension(:), allocatable, intent(out) :: nu

        allocate(fnu(size(f%fnu,1),size(f%fnu,2)))
        allocate(nu(size(freq_tab)))
        fnu=f%fnu; nu=freq_tab
    end subroutine convert_fluidvars_phatdisk

    subroutine convert_fluidvars_sariaf(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8) :: riaf_n0, riaf_beta
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp

        ncgs = f%rho
        bcgs = f%bmag
        ! new var rho2 is used to hold a second density, in this non-thermal e-
        ncgsnth= f%rho2
        tcgs= f%p
    end subroutine convert_fluidvars_sariaf

    subroutine convert_fluidvars_toy(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp
    
        ncgs = f%rho
        bcgs = f%bmag
        tcgs=0d0
    end subroutine convert_fluidvars_toy

    subroutine convert_fluidvars_powerlaw(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8) :: n0, beta, beta_trans
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        type (source_params), intent(in) :: sp

        beta_trans = 1d0
        ncgs = f%rho
        bcgs = f%bmag
        ncgsnth= f%rho2
        call monika_e(f%rho,f%rho,f%bmag,beta_trans, &
                      1d0/sp%muval-1d0, sp%gminval*(1d0/sp%muval-1d0),trat)
        tcgs= f%p/(1d0+trat)
    end subroutine convert_fluidvars_powerlaw

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! source param routines (for nonthermal: constant or tail) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    subroutine assign_source_params_type(sp,type)
        character(len=500), intent(in) :: type
        type (source_params), intent(inout) :: sp
        if(type=='const') then
         sp%type=CONST
        elseif(type=='tail') then
         sp%type=TAIL
        else
         write(6,*) 'ERROR in assign_source_params_type: type not recognized', type
        endif
    end subroutine assign_source_params_type

    subroutine initialize_source_params(sp,nup)
        type (source_params), intent(inout) :: sp
        integer, intent(in) :: nup
        allocate(sp%gmin(nup))
        allocate(sp%jetalpha(nup))
        allocate(sp%mu(nup))
    end subroutine initialize_source_params

    subroutine del_source_params(sp)
        type (source_params), intent(inout) :: sp
        deallocate(sp%gmin)
        deallocate(sp%jetalpha)
        deallocate(sp%mu)
    end subroutine del_source_params
        
    subroutine assign_source_params(sp,ncgs,tcgs,ncgsnth)
        type (source_params), intent(inout) :: sp
        real(kind=8), dimension(:), intent(in) :: ncgs,tcgs
        real(kind=8), dimension(:), intent(inout) :: ncgsnth
        real(kind=8), dimension(size(ncgs)) :: x,one,gmin,gmax,zero,factor
        zero=0d0
        one=1d0
        gmax=sp%gmax

        !write(6,*) 'fluid assign source params: ',sp%type,CONST,TAIL
        select case(sp%type)
         case (CONST)
            sp%gmin=sp%gminval
            sp%jetalpha=sp%jetalphaval
            sp%mu=sp%muval
         case (TAIL)
            sp%jetalpha=sp%jetalphaval
            sp%mu=sp%muval
            ! Fluid -> calc_gmin -> mu*tcgs -> emis means that calc_gmin is given a tcgs
            ! pre-mu correction and needs to be corrected here.                 
            call calc_gmin_subroutine(sp%p2,sp%mu*k*tcgs/m/c/c,sp%jetalpha,gmin,x)
            sp%gmin=merge(gmin,gmax/2d0,gmin.le.gmax)
            factor=merge(one,(gmax/2d0/gmin)**(sp%p2 - 2.),gmin.le.gmax)

            !when gmin is corrected for being too large, multiply ncgsnth by a corrective factor.
            !The correction (1-p) is already applied, so the correction p-2 is needed.
            ncgsnth=factor * merge(x*ncgs*sp%gmin**(1.-sp%p2),zero,x.gt.0d0)

        end select
    end subroutine assign_source_params

end module fluid_model
