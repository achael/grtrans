       program call_grtrans
       implicit none

       namelist /files/ ifile, ofile

       character(len=500) :: ifile, ofile
       character(len=500) :: readfile
       readfile='files.in'
       open(unit=8,file=readfile)
       read(8,nml=files)
       close(unit=8)
       call grtrans_run(ifile,ofile)
       end program
