!  converting GCP 2025 prior conc at ACTM grids to obs locations
!     *.gt3 -> *.grd
!  Dmitry Belikov, 2025/09/17

  Program f2_actm_ch4_at_obs_sites
    implicit none

    ![runtime configuration]
    character(len=512)         :: obs_file, i_dir, o_dir, o_fil, config_file
    character(len=16)          :: hcase
    integer                    :: s_year, f_year, nargs
    namelist /actm_f2/ hcase, obs_file, i_dir, o_dir, o_fil, s_year, f_year

    ![parameter]
    integer,parameter          :: nsite = 245
    integer,parameter          :: ntrac = 1
    integer,parameter          :: imax = 128, jmax = 64, kmax = 67
    integer,parameter          :: nmon = 12
    integer,parameter          :: iunit = 99

    integer,save               :: nday(nmon) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    character(len=*),parameter :: tracn(ntrac) = ['ch4_inv1']

    ![vars]
    real*4, allocatable        :: conc(:, :, :) !(nsite,ntrac,ntime)
    character                  :: sitem(nsite)*60, fname*512, ayr*15, msg*512
    real*4, dimension(nsite)   :: slat, slon, slev, lat, lon, lev, shgt
    integer, dimension(nsite)  :: ilat, ilon, ilev, land
    real*4                     :: oval(imax, jmax, kmax)
    integer                    :: iyear, nt, ntime, ios, nrecl
    integer                    :: i,j,k,n,ii,jj,kk
    integer                    :: hstart(3), hend(3), hsize

    ![GTOOL3]
    integer       idate(6)    !! date (year, month, day, hour, min., sec.)
    integer       itime       !! forecast time
    ! character*(*) hunit       !! history unit
    character*16  head(64)    !! GTOOL3 header

    ![read runtime configuration]
    nargs = command_argument_count()
    if (nargs /= 1) then
      print *, 'Usage: f2_actm_ch4_at_obs_sites <actm_f2.nml>'
      stop 1
    end if

    call get_command_argument(1, config_file, status=ios)
    if (ios /= 0 .or. len_trim(config_file) == 0) then
      print *, 'Invalid configuration argument'
      stop 1
    end if

    hcase = ''
    obs_file = ''
    i_dir = ''
    o_dir = ''
    o_fil = ''
    s_year = -1
    f_year = -1

    open(unit=iunit, file=trim(config_file), form='formatted', status='old', action='read', &
         iostat=ios, iomsg=msg)
    if (ios /= 0) then
      print *, 'Cannot open configuration: ', trim(config_file)
      print *, trim(msg)
      stop 1
    end if
    read(iunit, nml=actm_f2, iostat=ios, iomsg=msg)
    close(iunit)
    if (ios /= 0) then
      print *, 'Cannot read actm_f2 namelist: ', trim(msg)
      stop 1
    end if

    if ((trim(hcase) /= 'CYC' .and. trim(hcase) /= 'INCA') .or. &
        len_trim(obs_file) == 0 .or. len_trim(i_dir) == 0 .or. &
        len_trim(o_dir) == 0 .or. len_trim(o_fil) == 0 .or. &
        s_year < 1900 .or. f_year < s_year) then
      print *, 'Invalid actm_f2 configuration'
      stop 1
    end if

    print *, 'ACTM case     : ', trim(hcase)
    print *, 'Site list     : ', trim(obs_file)
    print *, 'Input GT3 dir : ', trim(i_dir)
    print *, 'Output dir    : ', trim(o_dir)
    print *, 'Output prefix : ', trim(o_fil)
    print *, 'Years         : ', s_year, f_year

    ! read site list
    open(unit=iunit, file=obs_file, form='formatted', status='old', action='read', iostat=ios, iomsg=msg)
    if (ios /= 0) then
      print *, 'Cannot open site list: ', trim(obs_file)
      print *, trim(msg)
      stop 1
    end if

    print *, '  List of obs files'
    ! Skip 2 lines
    do k = 1, 2
      read(iunit,*,iostat=ios,iomsg=msg)
      if (ios /= 0) then
        print *, 'Cannot read site-list header: ', k
        print *, trim(msg)
        stop 1
      end if
    end do

    do n = 1, nsite
      read(iunit,*,iostat=ios,iomsg=msg) i,j,sitem(n),slat(n),slon(n),slev(n), &
                                      lat(n),lon(n),lev(n),ilon(n),ilat(n),ilev(n),land(n)
      if (ios /= 0) then
        print *, 'Cannot read site ', n, ' of ', nsite
        print *, trim(msg)
        stop 1
      end if

      if (ilon(n) < 1 .or. ilon(n) > imax) then
        print *, 'Invalid longitude index: ', trim(sitem(n)), ilon(n)
        stop 1
      end if
      if (ilat(n) < 1 .or. ilat(n) > jmax) then
        print *, 'Invalid latitude index: ', trim(sitem(n)), ilat(n)
        stop 1
      end if
      if (ilev(n) < 1 .or. ilev(n) > kmax) then
        print *, 'Invalid level index: ', trim(sitem(n)), ilev(n)
        stop 1
      end if

      print '(i5,A,A60,2f10.2)', i,'  ',sitem(n),slat(n),slon(n)
    end do
    close(iunit)
    print *, ''
    !print '(20a4)', sitem(:)(1:3)
    !print '(20i4)', ilev, nsite

    ! read model data
    print*,'=== Year loop ==========================================='
    print *, '  Run from ', s_year, ' to ', f_year, ' for ', f_year-s_year+1, ' years'
    print*,'========================================================='

    DO iyear = s_year, f_year
      write(ayr,'(i4.4)') iyear

      nt = 365
      if (mod(iyear,4) == 0 .and. (mod(iyear,100) /= 0 .or. mod(iyear,400) == 0)) nt = 366
      ntime = nt * 24 ! 1-hourly files

      print '(A,3i6)', '  ACTM processing for ', iyear, nt, ntime
      allocate(conc(1:nsite,1:ntrac,1:ntime),stat=ios,errmsg=msg)
      if (ios /= 0) then
        print *, 'Cannot allocate concentration array: ', trim(msg)
        stop 1
      end if
!
      do i = 1, ntrac
        oval = 0.0

        ! input
        fname = trim(i_dir) // trim(ayr) // '/' // trim(tracn(i)) // '.gt3'
        open(unit=iunit, file=fname, form='unformatted', status='old', action='read', &
             convert='big_endian', iostat=ios, iomsg=msg)
        if (ios /= 0) then
          print *, 'Cannot open input: ', trim(fname)
          print *, trim(msg)
          stop 1
        end if

        do j = 1, ntime
          !--- if grd file
          ! read (iunit,rec=j) (((oval (ii,jj,kk), ii=1,imax), jj=1,jmax), kk=1,kmax)

          !--- if gt3 file
          read(iunit,iostat=ios,iomsg=msg) head
          if (ios /= 0) then
            print *, 'Cannot read GT3 header: ', trim(fname), ' record ', j
            print *, trim(msg)
            stop 1
          end if

          !print *, ' header written, oval next'
          !print *, head(29),  head(32), head(35) !, imax, jmax, kmax
          !read (head(25), '(I16)') itime
          !read (head(27), 900) (idate(ii),ii=1,6)
          !900 format( I4.4,I2.2,i2.2,' ',I2.2,I2.2,I2.2 )

          if (trim(adjustl(head(38))) /= 'UR4') then
            print *, 'Expected GT3 UR4 data, found: ', head(38), ' record ', j
            stop 1
          end if

          do k = 1, 3
            read(head(30+3*(k-1)),*,iostat=ios) hstart(k)
            if (ios /= 0) then
              print *, 'Invalid GT3 axis start: ', k, ' record ', j
              stop 1
            end if
            read(head(31+3*(k-1)),*,iostat=ios) hend(k)
            if (ios /= 0) then
              print *, 'Invalid GT3 axis end: ', k, ' record ', j
              stop 1
            end if
          end do

          if (any(hstart /= 1) .or. any(hend /= (/imax,jmax,kmax/))) then
            print *, 'GT3 grid mismatch at record ', j
            print *, 'Axis starts: ', hstart, ' ends: ', hend
            stop 1
          end if

          read(head(64),*,iostat=ios) hsize
          if (ios /= 0) then
            print *, 'Invalid GT3 data size at record ', j
            stop 1
          end if
          if (hsize /= imax*jmax*kmax) then
            print *, 'GT3 data size mismatch: ', hsize, ' record ', j
            stop 1
          end if

          read(iunit,iostat=ios,iomsg=msg) (((oval(ii,jj,kk), ii=1,imax), jj=1,jmax), kk=1,kmax)
          if (ios /= 0) then
            print *, 'Cannot read GT3 field: ', trim(fname), ' record ', j
            print *, trim(msg)
            stop 1
          end if

          do n = 1, nsite !conc(nsite,ntrac,ntime)
            conc(n,i,j) = oval(ilon(n),ilat(n),ilev(n))
          end do
        end do

        ! Check for unexpected additional records
        read(iunit,iostat=ios,iomsg=msg) head
        if (ios == 0) then
          print *, 'Input contains more than ', ntime, ' hourly records: ', trim(fname)
          stop 1
        else if (ios > 0) then
          print *, 'Error checking end of input: ', trim(msg)
          stop 1
        end if
        close(iunit)

        print '(A,i3,A,i5,2A)', '    Processing input >>    Tracer: ', i, ' Ntime: ', ntime, ' File: ', trim(fname)
      end do
!
      ! output
      fname = trim(o_dir) // trim(o_fil) // ayr(1:4) // '.grd'
      print '(2A)', '      Writing output >>    File: ', trim(fname)
      inquire(iolength=nrecl) ((conc(n,i,1), n=1,nsite), i=1,ntrac)
      open(unit=iunit, file=fname, form='unformatted', access='direct', status='replace', &
           action='write', recl=nrecl, iostat=ios, iomsg=msg)
      if (ios /= 0) then
        print *, 'Cannot open output: ', trim(fname)
        print *, trim(msg)
        stop 1
      end if

      do j = 1, ntime
        write(iunit,rec=j,iostat=ios,iomsg=msg) ((conc(n,i,j), n=1,nsite), i=1,ntrac)
        if (ios /= 0) then
          print *, 'Cannot write output record ', j, ': ', trim(fname)
          print *, trim(msg)
          stop 1
        end if
      end do
      close(iunit,iostat=ios,iomsg=msg)
      if (ios /= 0) then
        print *, 'Cannot close output: ', trim(fname)
        print *, trim(msg)
        stop 1
      end if
!
      deallocate(conc)
    ENDDO
!
  End program f2_actm_ch4_at_obs_sites