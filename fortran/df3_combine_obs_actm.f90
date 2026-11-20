!  converting GCP 2025 model prior conc and obs 
!     *.gt3 -> *.grd
  Module df3_combine_obs_actm_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none 
    
    ![runtime configuration]
    character(len=512) :: obs_file, obs_dir, i_dir, o_dir, i_fil, o_yr, o_mm, o_dy
    character(len=512) :: config_file
    character(len=16)  :: hcase
    integer            :: s_year, f_year, nyear
    namelist /d3_combine_obs_actm/ hcase, obs_file, obs_dir, i_dir, o_dir, i_fil, &
                       o_yr, o_mm, o_dy, s_year, f_year

    ![parameter]
    integer, parameter         :: nsite  = 245
    integer, parameter         :: ntrac  = 2
    integer, parameter         :: imax   = 128, jmax = 64, kmax = 67
    integer, parameter         :: nmon   = 12
    integer, parameter         :: nhr    = 24
    integer,parameter          :: iun1   = 97
    integer,parameter          :: iun2   = 98
    integer,parameter          :: oun1   = 99
    real,parameter             :: undf   = -999.0
    integer,save               :: nday(nmon) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    character(len=10),parameter :: tracn(ntrac) = ['obs       ', 'ch4_inv1  ']
    integer,parameter          :: obs_hour_start = 0 ! 0 for hours 0:23; 1 for hours 1:24

    ![vars]
    real*4, allocatable          :: conc (:,:,:,:,:,:), conc_dy (:,:,:,:,:)
    real*4, allocatable        :: conc_mm(:,:,:,:), conc_mf(:,:,:,:), conc_yr(:,:,:)
    character                  :: siten(nsite)*98, fname*512, flg*3, subset*4, cyr*9, msg*512
    real*4, dimension(nsite)   :: slat
    real*4                     :: xx, tracer
    integer                    :: iy, jrec, ierr, n_st, n_en, n, mday, indx(nsite)
    integer                    :: ist, jyr, jmn, jdy, jhr, jtr, iyear

    
  contains      

    ! --- Small helper for validating case-specific filenames
    pure function lowercase(value) result(out)
      character(len=*), intent(in) :: value
      character(len=len(value)) :: out
      integer :: j, ich
      out = value
      do j = 1, len(value)
        ich = iachar(value(j:j))
        if (ich >= iachar('A') .and. ich <= iachar('Z')) out(j:j) = achar(ich+32)
      end do
    end function lowercase

    !*****************************************************************
    !*** get_obs
    Subroutine get_obs()
      integer :: site_id, group_id, ih, nvalid, ndup, ios
      integer :: nhit(nyear,nmon,31,nhr)
      character*1024 :: line
      logical :: exists

      if (obs_hour_start /= 0 .and. obs_hour_start /= 1) stop 'Invalid obs_hour_start'
      do ist = n_st, n_en
        read(iun1,*,iostat=ierr,iomsg=msg) site_id, group_id, siten(ist), slat(ist)
        call check_io(ierr,'Read site list')
        if (.not.ieee_is_finite(slat(ist))) stop 'Non-finite site latitude'
        fname = trim(obs_dir) // trim(siten(ist)) // '.txt'
        print '(i5,2x,a)', ist, trim(fname)
        inquire(file=fname,exist=exists)
        if (.not.exists) then
          print *, '  File not found; observations remain missing: ', trim(fname)
          cycle
        end if
        open(iun2,file=fname,status='old',action='read',iostat=ierr,iomsg=msg)
        call check_io(ierr,'Open observations')
        nhit = 0
        nvalid = 0
        ndup = 0
        jrec = 0
        do
          read(iun2,'(A)',iostat=ierr,iomsg=msg) line
          if (is_iostat_end(ierr)) exit
          call check_io(ierr,'Read observation line')
          jrec = jrec+1
          line = adjustl(line)
          if (len_trim(line) == 0) cycle
          if (line(1:1) == '#' .or. line(1:1) == '!') cycle
          read(line,*,iostat=ios) jyr, jmn, jdy, jhr, tracer
          if (ios /= 0) then
            print *, 'Invalid observation record: ', trim(fname), ' line ', jrec
            print *, trim(line)
            stop 1
          end if
          if (jyr < s_year .or. jyr > f_year) cycle
          if (jmn < 1 .or. jmn > nmon) stop 'Invalid observation month'
          if (jdy < 1 .or. jdy > month_days(jyr,jmn)) stop 'Invalid observation day'
          ! Skip observations with an invalid or missing hour
          ih = jhr-obs_hour_start+1
          if (ih < 1 .or. ih > nhr) cycle
          if (.not.ieee_is_finite(tracer)) cycle
          if (tracer <= 0.1) cycle
          iy = jyr-s_year+1
          nhit(iy,jmn,jdy,ih) = nhit(iy,jmn,jdy,ih)+1
          if (nhit(iy,jmn,jdy,ih) == 1) then
            conc(1,ist,iy,jmn,jdy,ih) = tracer
          else
            ndup = ndup+1
            xx = conc(1,ist,iy,jmn,jdy,ih)
            conc(1,ist,iy,jmn,jdy,ih) = xx+(tracer-xx)/real(nhit(iy,jmn,jdy,ih))
          end if
          nvalid = nvalid+1
        end do
        close(iun2,iostat=ierr,iomsg=msg)
        call check_io(ierr,'Close observations')
        print '(A,I8,A,I8)', '  Valid observations: ', nvalid, ' Additional records in same hour: ', ndup
      end do
      print *, '  Finished reading ', subset, ' CH4 data'

      ! [retained original comments and commented code]
          !if ((subset .eq. 'NOAA') .or. (subset .eq. 'WDCG') .or. (subset .eq. 'ICOS')) read (iun2, *) flg
              !if (subset .eq. 'ObsP') tracer = tracer*1e9
              ! print*, jyr, iy, jmn, jdy, jhr, tracer, subset
    End subroutine get_obs

    !*****************************************************************
    !*** get_actm
    Subroutine get_actm()
      real*4 :: record(nsite,ntrac-1)
      integer :: nrecl

      inquire(iolength=nrecl) record
      do iyear = s_year, f_year
        iy = iyear-s_year+1
        write(cyr,'(i4.4)') iyear
        n = 0
        fname = trim(i_dir) // trim(i_fil) // trim(cyr) // '.grd'
        open(iun2,file=fname,form='unformatted',access='direct',status='old', &
             action='read',recl=nrecl,iostat=ierr,iomsg=msg)
        call check_io(ierr,'Open ACTM')
        do jmn = 1, nmon
          mday = month_days(iyear,jmn)
          do jdy = 1, mday
            do jhr = 1, nhr
              n = n+1
              read(iun2,rec=n,iostat=ierr,iomsg=msg) record
              if (ierr /= 0) print *, 'ACTM record ', n, ' year/month/day/hour ', iyear,jmn,jdy,jhr
              call check_io(ierr,'Read ACTM')
              do jtr = 2, ntrac
                conc(jtr,:,iy,jmn,jdy,jhr) = record(:,jtr-1)
              end do
            end do
          end do
        end do
        close(iun2,iostat=ierr,iomsg=msg)
        call check_io(ierr,'Close ACTM')
        print *, iyear, ' Time ', n
      end do

      ! [retained original comments and commented code]
              !print*, iyear, iy, jmn, jdy, jhr, n
              !jtr = 2
              !ist = 2
              !print'(6i6, f12.2)', jtr,ist,iy,jmn,jdy,jhr, conc(jtr,ist,iy,jmn,jdy,jhr)
        !print *, iyear, ' Time ', n, ' SiteN ', nsite, ' Ntrac ', ntrac
    End subroutine get_actm

    !*****************************************************************
    !*** Preparing monthly/daily means
    Subroutine get_mean()
      real*8 :: sum_m(ntrac), sum_d(ntrac)
      integer :: nm, nd

      conc_dy = undf
      conc_mm = undf
      do ist = 1, nsite
      do jyr = 1, nyear
      do jmn = 1, nmon
        sum_m = 0.D0
        nm = 0
        do jdy = 1, month_days(s_year+jyr-1,jmn)
          sum_d = 0.D0
          nd = 0
          do jhr = 1, nhr
            ! Use the same valid hours for observations and every model tracer
            if (.not.all(ieee_is_finite(conc(:,ist,jyr,jmn,jdy,jhr)))) cycle
            if (.not.all(conc(:,ist,jyr,jmn,jdy,jhr) > 0.1)) cycle
            sum_d = sum_d+dble(conc(:,ist,jyr,jmn,jdy,jhr))
            sum_m = sum_m+dble(conc(:,ist,jyr,jmn,jdy,jhr))
            nd = nd+1
            nm = nm+1
          end do
          if (nd >= 1) conc_dy(:,ist,jyr,jmn,jdy) = real(sum_d/nd,4)
        end do
        if (nm >= 1) conc_mm(:,ist,jyr,jmn) = real(sum_m/nm,4)
      end do
      end do
      end do

      ! [retained original comments and commented code]
        ! time = s_year+jyr+(jmn-0.5)/12.0
              !comment for 'no real data selection'
                ! Monthly averages
                ! Daily averages
                !print*, sum_m, sum_d, nm, nd
          !print*, 'daily', conc_dy(jtr,ist,jyr,jmn,jdy)
        !print'(a, 4i6, f12.2)', 'monthly ', jtr,ist,jyr,jmn, conc_mm(jtr,ist,jyr,jmn)
        !write(51,'(3i5,15f12.4)')k+s_year,l,15,ch4mm(i,k,l),time
        !print*, jtr,ist,jyr,jmn, conc_mm(jtr,ist,jyr,jmn)
    End subroutine get_mean

    !*****************************************************************
    !*** write output
    Subroutine write_out()
      integer :: nrecl, doy, days_year
      real*4 :: record(nsite,ntrac)

      inquire(iolength=nrecl) record
      ! yr
      fname = trim(o_dir) // trim(o_yr)
      print *, '   Print yearly mean for 1-st latitude-sorted site: ', trim(siten(indx(1)))
      print *, '   year   ', '   obs    ', ' tracer1  '
      call open_output(nrecl)
      do jyr = 1, nyear
        ! writing arranged by latitude
        record = transpose(conc_yr(:,indx,jyr))
        write(iun2,rec=jyr,iostat=ierr,iomsg=msg) record
        call check_io(ierr,'Write annual means')
        print '(I6,20F10.2)', s_year+jyr-1, conc_yr(:,indx(1),jyr)
      end do
      close(iun2,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close annual means')
      print *, '  Nyear ', nyear

      ! mon
      print *, ''
      print *, '=== MONTHLY OBSERVATIONS: BEGIN =========================='
      print *, 'Each row is one month; first column = year + (month-0.5)/12.'
      print *, 'Following columns: CH4 observations for the first sites in the site-list order.'
      print *, 'Units are unchanged from the observation files; no conversion is applied.'
      print *, 'Means use common valid observation/model hours; monthly gaps are not filled.'
      print *, '-999.00 means no valid paired hours in that month.'
      do ist = 1, min(12,nsite)
        print '(A,I3,A,A)', '  site ', ist, ': ', trim(siten(ist))
      end do
      write(*,'(a10)',advance='no') 'year'
      do ist = 1, min(12,nsite)
        write(*,'(1x,a7,i0)',advance='no') 'site', ist
      end do
      write(*,*)
      fname = trim(o_dir) // trim(o_mm)
      call open_output(nrecl)
      n = 0
      do jyr = 1, nyear
      do jmn = 1, nmon
        n = n+1
        record = transpose(conc_mm(:,:,jyr,jmn))
        write(iun2,rec=n,iostat=ierr,iomsg=msg) record
        call check_io(ierr,'Write monthly means')
        print '(20F10.2)', s_year+jyr-1+(jmn-0.5)/12.0, conc_mm(1,1:min(12,nsite),jyr,jmn)
      enddo
      enddo
      close(iun2,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close monthly means')
      print *, '  Monthly Recs ', n, ' SiteN ', nsite
      print *, '=== MONTHLY OBSERVATIONS: END ============================'
      print *, 'Console table: observations only, first ', min(12,nsite), ' sites.'
      print *, 'Binary file: observations and model tracers, all ', nsite, ' sites.'
      print *, 'Monthly output file: ', trim(fname)
      print *, 'Next block: daily MODEL means for the first sites, using paired valid hours.'
      print *, ''

      ! day
      ! The shared validity mask is applied in get_mean; no extra daily-only threshold
      ! where (conc_dy < 100.0)
      !   conc_dy = undf
      ! end where
      fname = trim(o_dir) // trim(o_dy)
      call open_output(nrecl)
      n = 0
      do jyr = 1, nyear
        doy = 0
        days_year = 337+month_days(s_year+jyr-1,2)
        do jmn = 1, nmon
          mday = month_days(s_year+jyr-1,jmn)
          do jdy = 1, mday
            n = n+1
            doy = doy+1
            record = transpose(conc_dy(:,:,jyr,jmn,jdy))
            write(iun2,rec=n,iostat=ierr,iomsg=msg) record
            call check_io(ierr,'Write daily means')
            print '(20F10.2)', s_year+jyr-1+(doy-0.5)/days_year, conc_dy(2,1:min(12,nsite),jyr,jmn,jdy)
          enddo
        enddo
      enddo
      close(iun2,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close daily means')
      print *, '  Daily Recs ', n
    End subroutine write_out
    

    !*****************************************************************
    !*** filled
    Subroutine filled()
      real*8 :: sum_y(ntrac)
      real*4 :: sorted_lat(nsite)
      integer :: ny

      conc_mf = conc_mm
      conc_yr = undf
      ! Data fill, first for monthly gaps, then for annual gaps
      do ist = 1, nsite
      do jyr = 2, nyear-1
      do jmn = 1, nmon
        if (all(conc_mm(:,ist,jyr,jmn) > 0.1)) cycle
        if (.not.all(conc_mm(:,ist,jyr-1,jmn) > 0.1)) cycle
        if (.not.all(conc_mm(:,ist,jyr+1,jmn) > 0.1)) cycle
        conc_mf(:,ist,jyr,jmn) = (conc_mm(:,ist,jyr-1,jmn)+conc_mm(:,ist,jyr+1,jmn))/2.0
      enddo
      enddo
      enddo

      ! sort
      indx = 0
      sorted_lat = slat
      call hpsort_eps_epw(nsite,sorted_lat,indx,0.01)

      ! Preparing annual means and growth rates
      do ist = 1, nsite
      do jyr = 1, nyear
        sum_y = 0.D0
        ny = 0
        do jmn = 1, nmon
          if (.not.all(ieee_is_finite(conc_mf(:,ist,jyr,jmn)))) cycle
          if (.not.all(conc_mf(:,ist,jyr,jmn) > 0.1)) cycle
          sum_y = sum_y+dble(conc_mf(:,ist,jyr,jmn))
          ny = ny+1
        enddo
        ! at least 10months of data needed
        if (ny >= 10) conc_yr(:,ist,jyr) = real(sum_y/ny,4)
      enddo
      enddo

      ! [retained original comments and commented code]
        !ch4m3(j,i,(k-1)*12+l) = conc_mf(j,i,k,l)
        !ns = 0
              !ns = ns+1
          ! at least 10months of data needed
        !if (jtr .eq. ntrac) print *, i, trim(site(i)), n2, nmon*nyear
    End subroutine filled
    
    ! source : http://epw.org.uk/uploads/Ford2/sourcefile/sort.f90.html                                                                            
    ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
    ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
    !---------------------------------------------------------------------
    Subroutine hpsort_eps_epw (n, ra, ind, eps)
      !---------------------------------------------------------------------
      !
      ! use kinds, ONLY : DP
      implicit none  
      !-input/output variables
      integer, intent(in)  :: n  
      real(4), intent(in)  :: eps
      integer              :: ind (n)  
      real(4)              :: ra (n)
      !-local variables
      integer              :: i, ir, j, l, iind  
      real(4)              :: rra  

      ! initialize index array
      IF (ind (1) .eq.0) then  
         DO i = 1, n  
            ind (i) = i  
         ENDDO
      ENDIF
      ! nothing to order
      IF (n.lt.2) return  
      ! initialize indices for hiring and retirement-promotion phase
      l = n / 2 + 1  

      ir = n  

      sorting: do 
        ! still in hiring phase
        IF ( l .gt. 1 ) then  
          l    = l - 1  
          rra  = ra (l)  
          iind = ind (l)  
          ! in retirement-promotion phase.
        ELSE  
          ! clear a space at the end of the array
          rra  = ra (ir)  
          !
          iind = ind (ir)  
          ! retire the top of the heap into it
          ra (ir) = ra (1)  
          !
          ind (ir) = ind (1)  
          ! decrease the size of the corporation
          ir = ir - 1  
          ! done with the last promotion
          IF ( ir .eq. 1 ) then  
            ! the least competent worker at all !
            ra (1)  = rra  
            !
            ind (1) = iind  
            exit sorting  
          ENDIF
        ENDIF
        ! wheter in hiring or promotion phase, we
        i = l  
        ! set up to place rra in its proper level
        j = l + l  
        !
        DO while ( j .le. ir )  
          IF ( j .lt. ir ) then  
            ! compare to better underling
            IF ( hslt( ra (j),  ra (j + 1), eps ) ) then  
              j = j + 1  
              !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
              ! this means ra(j) == ra(j+1) within tolerance
              !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
            ENDIF
          ENDIF
          ! demote rra
          IF ( hslt( rra, ra (j), eps ) ) then  
              ra (i) = ra (j)  
              ind (i) = ind (j)  
              i = j  
              j = j + j  
           !else if ( .not. hslt ( ra(j) , rra ) ) then
              !this means rra == ra(j) within tolerance
              ! demote rra
             ! if (iind.lt.ind (j) ) then
             !    ra (i) = ra (j)
             !    ind (i) = ind (j)
             !    i = j
             !    j = j + j
             ! else
                 ! set j to terminate do-while loop
             !    j = ir + 1
             ! endif
              ! this is the right place for rra
          ELSE
              ! set j to terminate do-while loop
              j = ir + 1  
          ENDIF
        ENDDO
        ra (i) = rra  
        ind (i) = iind  
      END DO sorting    
    End subroutine hpsort_eps_epw

    ! internal function 
    ! compare two real number and return the result
    logical function hslt( a, b, eps)
      REAL(4)        :: a, b, eps
      IF( abs(a-b) <  eps ) then
        hslt = .false.
      ELSE
        hslt = ( a < b )
      end if
    end function hslt

    !*****************************************************************
    !*** month_days
    integer Function month_days(year, month) result(days)
      integer, intent(in) :: year, month

      if (month < 1 .or. month > nmon) stop 'Invalid month'
      days = nday(month)
      if (month == 2) then
        if (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) days = 29
      end if
    End function month_days

    !*****************************************************************
    !*** check_io
    Subroutine check_io(status, action)
      integer, intent(in) :: status
      character(len=*), intent(in) :: action

      if (status == 0) return
      print *, trim(action), ': ', trim(fname)
      print *, trim(msg)
      stop 1
    End subroutine check_io

    !*****************************************************************
    !*** open_output
    Subroutine open_output(nrecl)
      integer, intent(in) :: nrecl

      open(iun2,file=fname,form='unformatted',access='direct',status='replace', &
           action='write',recl=nrecl,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open output')
    End subroutine open_output

  End module df3_combine_obs_actm_mod


  Program df3_combine_obs_actm
    use df3_combine_obs_actm_mod
    implicit none

    ! --- Read case-specific runtime namelist
    if (command_argument_count() /= 1) then
      print *, 'Usage: df3_combine_obs_actm <d3_combine_obs_actm.nml>'
      stop 1
    end if
    call get_command_argument(1, config_file, status=ierr)
    if (ierr /= 0 .or. len_trim(config_file) == 0) stop 'Invalid configuration argument'

    hcase = ''
    obs_file = ''
    obs_dir = ''
    i_dir = ''
    o_dir = ''
    i_fil = ''
    o_yr = ''
    o_mm = ''
    o_dy = ''
    s_year = -1
    f_year = -1
    open(iun1,file=trim(config_file),status='old',action='read',iostat=ierr,iomsg=msg)
    if (ierr /= 0) then
      print '(A,A)', 'Cannot open D3 combine config: ', trim(config_file)
      print '(A)', trim(msg)
      stop 1
    end if
    read(iun1,nml=d3_combine_obs_actm,iostat=ierr,iomsg=msg)
    close(iun1)
    if (ierr /= 0) then
      print '(A,A)', 'Cannot read D3 combine namelist: ', trim(msg)
      stop 1
    end if
    if (trim(hcase) /= 'CYC' .and. trim(hcase) /= 'INCA') stop 'Invalid D3 combine case'
    if (s_year < 1900 .or. f_year < s_year) stop 'Invalid D3 combine years'
    if (len_trim(obs_file) == 0 .or. len_trim(obs_dir) == 0 .or. &
        len_trim(i_dir) == 0 .or. len_trim(o_dir) == 0 .or. &
        len_trim(i_fil) == 0 .or. len_trim(o_yr) == 0 .or. &
        len_trim(o_mm) == 0 .or. len_trim(o_dy) == 0) stop 'Incomplete D3 combine namelist'
    if (index(i_fil,trim(lowercase(hcase))) == 0) stop 'D3 combine model prefix does not match case'
    if (index(o_yr,trim(lowercase(hcase))) == 0 .or. &
        index(o_mm,trim(lowercase(hcase))) == 0 .or. &
        index(o_dy,trim(lowercase(hcase))) == 0) stop 'D3 combine output filename does not match case'
    if (i_dir(len_trim(i_dir):len_trim(i_dir)) /= '/') i_dir = trim(i_dir)//'/'
    if (o_dir(len_trim(o_dir):len_trim(o_dir)) /= '/') o_dir = trim(o_dir)//'/'
    if (obs_dir(len_trim(obs_dir):len_trim(obs_dir)) /= '/') obs_dir = trim(obs_dir)//'/'
    nyear = f_year-s_year+1

    print '(A,A)', '      ACTM case     : ', trim(hcase)
    print '(A,A)', '      Site list     : ', trim(obs_file)
    print '(A,A)', '      Obs dir       : ', trim(obs_dir)
    print '(A,A)', '      Input ACTM dir: ', trim(i_dir)
    print '(A,A)', '      Output dir    : ', trim(o_dir)
    print '(A,A)', '      Model prefix  : ', trim(i_fil)
    print '(A,2I6)', '      Years         : ', s_year, f_year

    fname = ''
    msg = ''
    if (ntrac < 2) stop 'Need observations and at least one model tracer'
    allocate(conc(1:ntrac,1:nsite,1:nyear,1:nmon,1:31,1:nhr),stat=ierr,errmsg=msg)
    call check_io(ierr,'Allocate hourly concentrations')
    allocate(conc_dy(1:ntrac,1:nsite,1:nyear,1:nmon,1:31),stat=ierr,errmsg=msg)
    call check_io(ierr,'Allocate daily concentrations')
    
    allocate(conc_mm(ntrac,nsite,nyear,nmon),stat=ierr,errmsg=msg)
    call check_io(ierr,'Allocate monthly means')
    allocate(conc_mf(ntrac,nsite,nyear,nmon),stat=ierr,errmsg=msg)
    call check_io(ierr,'Allocate filled monthly means')
    allocate(conc_yr(ntrac,nsite,nyear),stat=ierr,errmsg=msg)
    call check_io(ierr,'Allocate annual means')

    ! site list file
    print*,''
    print *, '=== Read obs_file ======================================='
    print*, '   ', obs_file
    print*,''
    
    fname = obs_file
    open(iun1,file=obs_file,status='old',action='read',iostat=ierr,iomsg=msg)
    call check_io(ierr,'Open site list')
    ! Skip 2 lines
    read(iun1,*,iostat=ierr,iomsg=msg)
    call check_io(ierr,'Read site-list header')
    read(iun1,*,iostat=ierr,iomsg=msg)
    call check_io(ierr,'Read site-list header')

    ! init conc
    conc = undf
    conc_dy = undf
    indx = 0
    slat = undf
    siten = ''
    conc_mm = undf
    conc_yr = undf
    conc_mf = undf

    ! get obs data from All 
    print*,''
    print *, '=== Run get_obs ======================================='
    n_st = 1
    ! n_en = 94
    n_en = nsite
    subset = 'obs'
    call get_obs()
    close(iun1,iostat=ierr,iomsg=msg)
    call check_io(ierr,'Close site list')

    ! get ACTM
    print *, ''
    print *, '=== Run get_actm ======================================'
    call get_actm()
    
    ! get mean
    print *, ''
    print *, '=== Run get_mean ======================================'
    call get_mean()
    
    ! get conc_mf
    print *, ''
    print *, '=== Run filled ========================================'
    call filled()

    ! write out
    print *, ''
    print *, '=== Run write_out ====================================='
    call write_out()
    deallocate(conc,conc_dy,conc_mm,conc_mf,conc_yr)
    
  End program df3_combine_obs_actm