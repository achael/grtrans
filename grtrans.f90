subroutine grtrans_run(ifile,outfile)
   use omp_lib
   use grtrans_inputs
   use pgrtrans
   use grtrans, only: grtrans_driver
   use ray_trace, only: ray_set,initialize_raytrace_camera, &
        kwrite_raytrace_camera, del_raytrace_camera
   use fluid_model, only: load_fluid_model, unload_fluid_model, &
        advance_fluid_timestep, source_params, assign_source_params_type, &
        fluid_args, assign_fluid_args
   use emissivity, only: emis_params
   use geodesics, only: initialize_pixels, geokerr_args, &
     del_geokerr_args,initialize_geokerr_args, initialize_geo_tabs
   use chandra_tab24, only: load_chandra_tab24, del_chandra_tab24

   implicit none

   character(len=500), intent(in) :: outfile,ifile
   integer :: nextra=0, inum, gunit, i, ncams, j, m, l, nparams
   integer :: nthreads, threadnum, iii
   real (kind=8) :: wtime

   type (ray_set), dimension(:), allocatable :: c
   type (geokerr_args) :: gargs
   type (fluid_args) :: fargs
   type (source_params), dimension(:), allocatable :: sparams
   type (emis_params) :: eparams
   character(len=500), dimension(3) :: knames,kdescs

   knames(1)='nx'; knames(2)='ny'; kdescs(1)='# x pixels'
   kdescs(2)='# y pixels'
   knames(3)='nu'; kdescs(3)='Frequency (Hz)'
   gunit=12

   ! read input file (read_inputs.f90)
   call read_inputs(ifile)
   write(6,*) 'grtrans ifile: ',trim(ifile)

   ! grtrans_main is in pgrtrans.f90
   call grtrans_main(standard,mumin,mumax,nmu,phi0,spin,&
        uout,uin, rcut, nrotype, gridvals, nn,i1,i2,fname, dt, nt, nload, &
        nmdot, mdotmin, mdotmax,ename, mbh, nfreq, fmin, fmax, muval,&
        gmin, gmax,p1, p2, fpositron,jetalpha, stype,use_geokerr, nvals, iname,&
        cflag, extra, debug,outfile,fdfile,fhfile,fgfile,fsim,fnt,findf,fnfiles,fjonfix, &
        fnw,fnfreq_tab,fnr,foffset,fdindf,fmagcrit,frspot,fr0spot,fn0spot,ftscl,frscl, &
        fwmin,fwmax,ffmin,ffmax,frmax,fsigt,ffcol,fmdot,fnscl,fnnthscl,fnnthp,fbeta, &
        fbl06,fnp,ftp,frin,frout,fthin,fthout,fphiin,fphiout,fscalefac,sigcut, &
        betaeconst, betaecrit, ximax, bscl,pscl,pegasratio,&
        coefindx, epotherargs,nepotherargs)

   ! deallocate
   call del_pgrtrans_data()

   return

end subroutine grtrans_run
