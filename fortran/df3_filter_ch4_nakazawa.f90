!this one is same as main_wk.f90, but uses daily data from the model
! Prabir K Patra, May 02, 2021
  program df3_filter_ch4_nakazawa
    implicit double precision (a-h,o-z) !none
    external lsfbf4

    ! --- Runtime input/output paths, no recompilation for CYC / INCA
    character(len=1024) :: i_file, o_file, config_file
    character(len=16)   :: hcase
    integer             :: ios, recl_in, recl_out
    namelist /d3_filter_ch4_nakazawa/ hcase, i_file, o_file

    ![parameters]                      
    integer, parameter          :: s_year = 1999, f_year=2025, nyear=f_year-s_year+1
    integer, parameter          :: nmon   = 12
    integer, parameter          :: nstn   = 245
    integer, parameter          :: tlen   = 9862 ! from df3_combine_obs_actm.f90
    integer, parameter          :: ntrac  = 2, ntype=2, ndat2=nstn
    integer, parameter          :: KD     = 30000, KO=60000
    real*8, parameter           :: outvalue=-999.999

    integer                     :: n1,mon,yearin,monthin,dayin,nmis(nstn)       !for input to the fitting routine
    dimension                   :: mon(nmon), yearin(kd),monthin(kd),dayin(kd)  !these two line shouldn't be
    dimension                   :: valuein(kd), tim(KO),yout(KO),ytrn(KO)       !touched for proper running...
    real, dimension(ko)         :: curve,trend,resid,year
    real                        :: conc(nstn,ntrac,tlen), conct(tlen), time(tlen)
    real *4                     :: conc3h(nstn,ntrac,tlen), fits(3,nstn,ntrac,tlen), sigm
    character*98                :: filename,dirname
    character*9                 :: site(nstn)
    data                           mon / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
    

    ! --- Read runtime configuration passed by Python
    hcase = ''
    i_file = ''
    o_file = ''
    config_file = ''
    call get_command_argument(1, config_file)
    if (len_trim(config_file).eq.0) then
      print *, 'ERROR: usage: df3_filter_ch4_nakazawa <d3_filter_ch4_nakazawa.nml>'
      stop 1
    endif

    open(11, file=trim(config_file), status='old', action='read', iostat=ios)
    if (ios.ne.0) then
      print *, 'ERROR: cannot open configuration: ', trim(config_file)
      stop 2
    endif
    read(11, nml=d3_filter_ch4_nakazawa, iostat=ios)
    close(11)
    if (ios.ne.0) then
      print *, 'ERROR: bad configuration: ', trim(config_file)
      stop 3
    endif
    if (trim(hcase).ne.'CYC' .and. trim(hcase).ne.'INCA') then
      print *, 'ERROR: invalid hcase: ', trim(hcase)
      stop 4
    endif
    if (len_trim(i_file).eq.0 .or. len_trim(o_file).eq.0) then
      print *, 'ERROR: input/output filenames must be specified'
      stop 5
    endif
    print *, 'Case   : ', trim(hcase)
    print *, 'Input  : ', trim(i_file)
    print *, 'Output : ', trim(o_file)

    fits = -999.0
    inquire(iolength=recl_in) conc(:,:,1)
    open(50,file=i_file,access='direct',form='unformatted',status='old', recl=recl_in)
    do n = 1, tlen 
      read (50, rec=n) ((conc(i, j, n), i = 1, nstn), j = 1, ntrac)
    end do; 
    close (50); 
    print *, 'i/p file ', trim(i_file), n

    do jj = 1, ntrac
      nmis = 0
      do ii = 1, nstn
        print *, '  '
        print *, ' ---------------------------------------------- '      
        print *, ' Processing for jj, ii', jj, ii
        nval=0; valuein=outvalue; yearin=0; monthin=0; dayin=0; tim=0.; yout=0.; ytrn=0.; nreal=0
        do i = 1, nyear
          iyear = s_year + i - 1
          do j = 1, nmon
            monday = mon(j)
            if ( ((iyear/4).eq.(iyear/float(4))) .and. (j.eq.2) ) then 
                  monday = 29 !; print *, 'Leap year Feb_day', iyear, monday
            end if
            do l = 1, monday
              nval = nval+1; yearin(nval) = iyear; monthin(nval) = j; dayin(nval) = l
              timein = iyear + (j-1)/12.0 + ((l-1)/float(monday))/12.0
              time(nval) = timein
              if ( (conc(ii, jj, nval) .ge. 900.0) .and. (conc(ii, jj, nval) .le. 5000.0) ) then
                valuein(nval) = conc(ii, jj, nval)
                nreal = nreal+1
              else
                valuein(nval) = outvalue
                !   if (nreal.le.0) nmis(ns) = nmis(ns)+1
                !   nmis(ns) = nmis(ns)+1
              end if
    !         print '(4i5,4f12.4)', nval, yearin(nval), monthin(nval),dayin(nval),time(nval),timein, valuein(nval)
            end do
          end do
        end do
        ! pause 'is this the way you want to input the data, type GO to confirm'
        if ( nreal .le. 48 ) go to 97  ! no fitting


        print '(a,4i6)', ' ndata before  calling lsfbf4 = ', nval, kk, jj, ii
        call lsfbf4(nval,yearin,monthin,dayin,valuein,tim,yout,ytrn,ndat,outvalue)
        print *, ' ndata after calling lsfbf4 = ', nval, ndat
     
        !  print *, tim, ' fit time'; pause; print *, time, ' input time'; pause
        do i1 = 1, tlen+30
          do i2 = 1, tlen
          if ( abs(tim(i1)-time(i2)) .le. 0.0013 ) then
            fits (1, ii, jj, i2) = conc(ii,jj,i2); 
            fits (2, ii, jj, i2) = yout(i1) 
            fits (3, ii, jj, i2) = ytrn(i1) 
!           print '(2f12.4,9f9.2)', time(i2), tim(i1), fits (1, ii, jj, i2), fits (2, ii, jj, i2)
         end if; end do; end do
!        print '(3i5,9f9.2)', kk, jj, ii, fits (kk, jj, ii, 1) !; pause
     
      97 continue
      end do
    end do

    ! write synoptic variations here ! fits:ctl, fits_ndg5:WN
    inquire(iolength=recl_out) (((fits(k,i,j,1), i = 1, nstn), k = 1, 3), j = 1, ntrac)
    open(50,file=o_file,access='direct',form='unformatted',status='replace', recl=recl_out)
    do n = 1, tlen
      write (50, rec=n) (((fits(k, i, j, n), i = 1, nstn), k = 1, 3), j = 1, ntrac)
    end do
    close (50)
    print *, 'o/p file ', trim(o_file), 'records ', n

    stop
  end program df3_filter_ch4_nakazawa
