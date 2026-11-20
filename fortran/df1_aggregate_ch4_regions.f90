  Module df1_aggregate_ch4_regions
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    ![input parameter]
    ! --- Runtime paths; fixed 1998-2026 period is intentional.
    character(len=1024)        :: flux_file, mask_file, uncr_file, config_file
    character(len=16)          :: hcase
    integer                    :: cfg_s_year, cfg_f_year, cfg_ios
    namelist /d1_aggregate_ch4_regions/ hcase, flux_file, mask_file, uncr_file, cfg_s_year, cfg_f_year

    ![parameter]
    integer,parameter          :: s_year = 1998, f_year = 2026, nyear = f_year-s_year+1
    integer,parameter          :: r8 = 8
    real(r8),parameter         :: pi = 3.1415926535897932384626433832795_r8
    real(r8),parameter         :: r0 = 6371.D3                      ! radius [m]
    integer,parameter          :: imax = 128, jmax = 64, kmax = 1
    integer,parameter          :: nemis = 1, nr54 = 54
    integer,parameter          :: nmon = 12
    integer,parameter          :: ounit = 99, munit = 11, funit = 98

    integer,save               :: nday(nmon) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    ![grid vars]
    real(r8)                   :: alon(imax)          ! longitude for GCM
    real(r8)                   :: dlon                ! diff. of longitude for GCM grid
    real(r8)                   :: alat(0:jmax+1)      ! latitude
    real(r8)                   :: blat(jmax+1)        ! boundary latitude
    real(r8)                   :: dlat(jmax)          ! latitude weight
    real(r8)                   :: gar_t42(jmax)       ! grid size jmax [m**2]
    real*4                     :: lnd54(imax,jmax), oval(imax,jmax,kmax)
    real(r8)                   :: garea(jmax), emis(imax,jmax,nmon,nemis)
    real(r8)                   :: ps54(nr54,nyear,nmon,nemis)
    real(r8)                   :: summ, monthly, cellmass, scale, total_all(nyear,nmon)
    real(r8)                   :: outside54(nyear,nmon)
    character*16               :: head(64), cyear
    character*512              :: fname, msg
    integer                    :: i, j, k, iy, imn, iyear, ifl, ierr, nrecl, ndays
    integer                    :: id54(imax,jmax)

    ! One total net flux field, including negative soil uptake
    character(len=*),parameter :: etype(nemis) = ['total']
    character(len=3),parameter :: regn54(nr54) = ['CA1','CA2','ALK','CA3','CA4','US1', &
                                                'US2','US3','US4','MEX','TA1','TA2', &
                                                'TA3','TA4','SA1','SA2','SA3','SA4', &
                                                'NA1','NA2','NA3','NA4','NA5','SF1', &
                                                'SF2','SF3','SF4','SB1','SB2','SB3', &
                                                'SB4','SB5','FSU','ARB','WIN','NEI', &
                                                'SIN','SWC','SEC','KAZ','MGL','NEC', &
                                                'JPK','SE1','SE2','AU1','AU2','AU3', &
                                                'AU4','EUS','UKR','SCN','RUS','SPN']
  contains

    !*****************************************************************
    !*** read runtime configuration, no changes to regional flux calculations
    Subroutine read_config()
      hcase = ''
      flux_file = ''
      mask_file = ''
      uncr_file = ''
      config_file = ''
      cfg_s_year = s_year
      cfg_f_year = f_year
      call get_command_argument(1, config_file)
      if (len_trim(config_file) == 0) error stop 'Usage: df1_aggregate_ch4_regions <d1_aggregate_ch4_regions.nml>'
      open(10,file=trim(config_file),status='old',action='read',iostat=cfg_ios)
      if (cfg_ios /= 0) error stop 'Cannot open D1 regions namelist'
      read(10,nml=d1_aggregate_ch4_regions,iostat=cfg_ios)
      close(10)
      if (cfg_ios /= 0) error stop 'Invalid D1 regions namelist'
      if (trim(hcase) /= 'CYC' .and. trim(hcase) /= 'INCA') error stop 'Invalid D1 regions case'
      if (cfg_s_year /= s_year .or. cfg_f_year /= f_year) error stop 'D1 regions period must remain 1998-2026'
      if (len_trim(flux_file) == 0 .or. len_trim(mask_file) == 0 .or. &
          len_trim(uncr_file) == 0) error stop 'D1 regions input/output path is empty'
      print '(2A)', '      Case          : ', trim(hcase)
      print '(2A)', '      Input flux GT3: ', trim(flux_file)
      print '(2A)', '      Regional mask : ', trim(mask_file)
      print '(2A)', '      Regional table: ', trim(uncr_file)
    End subroutine read_config

    !*****************************************************************
    !*** setgrid
    Subroutine setgrid_actm()
!     [internal]
      real(r8) :: ctheta(jmax)        !! Gauss Latitude (theta)
      real(r8) :: gw(jmax)            !! Gauss weight
      real(r8) :: edge
      integer :: i, j

      do i = 1, imax
        alon(i) = 360.D0*(i-1)/imax
      end do
      dlon = 360.D0/imax
      call gauss(ctheta,gw)
      dlat = gw
      alat(0) = 90.D0
      alat(jmax+1) = -90.D0
      do j = 1, jmax
        alat(j) = asin(ctheta(j))*180.D0/pi
      end do
      edge = 1.D0
      blat(1) = 90.D0
      do j = 1, jmax
        edge = edge-gw(j)
        blat(j+1) = asin(max(-1.D0,min(1.D0,edge)))*180.D0/pi
      end do
      blat(jmax+1) = -90.D0
      !print*, alon; print*, alat
      return
    End subroutine setgrid_actm

    !*****************************************************************
    !*** set gauss grid
    Subroutine gauss(ctheta, gw)
    !
    ! [output]
      real(r8),intent(out) :: ctheta(jmax), gw(jmax)
    !
    ! [work]
      real(r8) :: p0, p1, pn, dp
    !
    ! [internal parm]
      real(r8),parameter :: deltp = 1.D-14
      integer,parameter :: itrmax = 100
    !
    ! [internal work]
      integer :: n, j, iter
      real(r8) :: x0, x1

      do j = 1, (jmax+1)/2
        ! step 1
        x0 = cos(pi*(real(j,r8)-0.25D0)/(real(jmax,r8)+0.5D0))
        do iter = 1, itrmax
          ! step 2 - Pn
          p0 = 1.D0
          p1 = x0
          do n = 2, jmax
            pn = ((2.D0*n-1.D0)*x0*p1-(n-1.D0)*p0)/n
            p0 = p1
            p1 = pn
          enddo
          ! step 3 - d/d Pn
          dp = jmax*(x0*p1-p0)/(x0*x0-1.D0)
          ! step 4 - newton
          x1 = x0-p1/dp
          if (abs(x1-x0) < deltp) exit
          x0 = x1
        enddo
        if (iter > itrmax) stop 'Gaussian latitude iteration failed'
        x0 = x1
        p0 = 1.D0
        p1 = x0
        do n = 2, jmax
          pn = ((2.D0*n-1.D0)*x0*p1-(n-1.D0)*p0)/n
          p0 = p1
          p1 = pn
        enddo
        dp = jmax*(x0*p1-p0)/(x0*x0-1.D0)
        ctheta(j) = x0
        ctheta(jmax-j+1) = -x0
        ! step 5
        gw(j) = 2.D0/((1.D0-x0*x0)*dp*dp)
        gw(jmax-j+1) = gw(j)
      enddo
      if (abs(sum(gw)-2.D0) > 1.D-12) stop 'Invalid Gaussian weight sum'
      return
    End subroutine gauss

    !*****************************************************************
    !*** for reading to gt3 formatted files
    Subroutine gt3read()
      integer :: i, j, k, starts(3), ends(3), nvalues, idate(6), ios
      real*4 :: cell_missing

      read(funit,iostat=ierr,iomsg=msg) head
      call check_io(ierr,'Read GT3 header')
      if (trim(adjustl(head(38))) /= 'UR4') stop 'Expected GT3 UR4 field'
      if (trim(adjustl(head(29))) /= 'GLON128') stop 'Unexpected GT3 longitude axis'
      if (trim(adjustl(head(32))) /= 'GGLA64') stop 'Unexpected GT3 latitude axis'
      if (trim(adjustl(head(35))) /= 'SFC1') stop 'Expected GT3 surface field'
      if (trim(adjustl(head(16))) /= 'kg-CH4/m2/s') then
        print *, 'Unexpected flux units: ', head(16)
        stop 1
      end if
      do k = 1, 3
        read(head(30+3*(k-1)),*,iostat=ios) starts(k)
        if (ios /= 0) stop 'Invalid GT3 axis start'
        read(head(31+3*(k-1)),*,iostat=ios) ends(k)
        if (ios /= 0) stop 'Invalid GT3 axis end'
      end do
      if (any(starts /= 1) .or. any(ends /= (/imax,jmax,kmax/))) stop 'GT3 dimensions mismatch'
      read(head(64),*,iostat=ios) nvalues
      if (ios /= 0) stop 'Invalid GT3 field size'
      if (nvalues /= imax*jmax*kmax) stop 'GT3 field size mismatch'
      read(head(27),'(I4,I2,I2,1X,I2,I2,I2)',iostat=ios) idate
      if (ios /= 0) stop 'Invalid GT3 date'
      if (idate(1) /= iyear .or. idate(2) /= imn) then
        print *, 'Expected year/month ', iyear, imn, ' found ', idate
        stop 'GT3 time sequence mismatch'
      end if
      if (idate(3) < 1 .or. idate(3) > month_days(iyear,imn)) stop 'Invalid GT3 day'
      if (idate(4) < 0 .or. idate(4) > 23) stop 'Invalid GT3 hour'
      if (any(idate(5:6) < 0) .or. any(idate(5:6) > 59)) stop 'Invalid GT3 minute or second'
      read(funit,iostat=ierr,iomsg=msg) (((oval(i,j,k), i=1,imax), j=1,jmax), k=1,kmax)
      call check_io(ierr,'Read GT3 field')
      if (.not.all(ieee_is_finite(oval))) stop 'Non-finite input flux'
      read(head(39),*,iostat=ios) cell_missing
      if (ios == 0) then
        if (any(oval == cell_missing)) stop 'Missing values in input flux'
      end if
      return
    End subroutine gt3read

    !*****************************************************************
    !*** month_days
    integer Function month_days(year, month) result(days)
      integer,intent(in) :: year, month

      if (month < 1 .or. month > nmon) stop 'Invalid month'
      days = nday(month)
      if (month == 2) then
        if (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) days = 29
      end if
    End function month_days

    !*****************************************************************
    !*** check_io
    Subroutine check_io(status, action)
      integer,intent(in) :: status
      character(len=*),intent(in) :: action

      if (status == 0) return
      print *, trim(action), ': ', trim(fname)
      print *, 'Year/month: ', iyear, imn
      print *, trim(msg)
      stop 1
    End subroutine check_io

  End module df1_aggregate_ch4_regions
!
!*****************************************************************
!
  Program mainprog
    use df1_aggregate_ch4_regions
    implicit none

    print *, 'Start program'
    call read_config()
    iyear = 0
    imn = 0
    !iyear = 2000
    !write (cyear, '(i4.4)') iyear
    if (nemis /= 1) stop 'This program reads one total flux field; nemis must be 1'

    fname = mask_file
    inquire(iolength=nrecl) lnd54
    open(munit,file=mask_file,status='old',form='unformatted',access='direct', &
         action='read',recl=nrecl,iostat=ierr,iomsg=msg)
    call check_io(ierr,'Open regional mask')
    read(munit,rec=1,iostat=ierr,iomsg=msg) lnd54
    call check_io(ierr,'Read 54-region mask')
    close(munit,iostat=ierr,iomsg=msg)
    call check_io(ierr,'Close regional mask')
    if (.not.all(ieee_is_finite(lnd54))) stop 'Non-finite 54-region mask'

    ! Keep cells outside the named regions separate; never silently discard their budget
    id54 = 0
    do j = 1, jmax
    do i = 1, imax
      do k = 1, nr54
        if (abs(dble(lnd54(i,j))-k) <= 0.01D0) id54(i,j) = k
      end do
    end do
    end do
    print *, 'Mask 54: assigned/unassigned cells ', count(id54 > 0), count(id54 == 0)
    print *, 'Mask 54 min/max: ', minval(lnd54), maxval(lnd54)
    if (count(id54 > 0) == 0) stop 'No recognized regions in mask 54; check mask format/endian'

    call setgrid_actm()
    ! <set Facter from (1degx1deg/1yr) to (1mx1m/1s)>
    ! Gaussian weights sum to 2 in this version; use 2*pi, not 4*pi
    do j = 1, jmax
      gar_t42(j) = 2.D0*pi*r0**2*dlat(j)/imax
    end do
    garea = gar_t42
    !print *, garea
    print '(A,ES20.12)', 'Earth area [m2]: ', sum(garea)*imax

    ps54 = 0.D0
    total_all = 0.D0
    outside54 = 0.D0
    fname = flux_file
    open(funit,file=flux_file,form='unformatted',convert='big_endian',status='old', &
         action='read',iostat=ierr,iomsg=msg)
    call check_io(ierr,'Open flux file')
    print '(2A)', 'Input ', trim(flux_file)
    print *, 'Monthly net CH4 mass [Tg/month]; negative soil uptake is retained.'
    print *, 'Year Mon      Global       Reg54     Other54    Year-to-date'

    DO iy = 1, nyear
      iyear = s_year+iy-1
      !write (cyear, '(i4.4)') iyear
      emis = 0.D0
      summ = 0.D0
      ! read
      do imn = 1, nmon
        call gt3read !( 51, oval, im, jm, km, head )
        emis(:,:,imn,1) = dble(oval(:,:,1)) !; print *, head (27)
        ndays = month_days(iyear,imn)
        monthly = 0.D0
        scale = 0.D0
        ! convert
        do j = 1, jmax
        do i = 1, imax
          ! from Kg-CH4/m2/sec  --> Kg-CH4/grid/day
          emis(i,j,imn,1) = emis(i,j,imn,1)*garea(j)*86400.D0
          cellmass = emis(i,j,imn,1)*ndays/1.D9
          monthly = monthly+cellmass
          scale = scale+abs(cellmass)
          ! map to regions
          !54 inversion regions
          k = id54(i,j)
          if (k > 0) then
            ps54(k,iy,imn,1) = ps54(k,iy,imn,1)+cellmass
          else
            outside54(iy,imn) = outside54(iy,imn)+cellmass
          end if
          !12 aggregated regions
          ! Not used: only the first mask record (regions 1:54) is processed.
        end do
        end do
        total_all(iy,imn) = monthly
        summ = summ+monthly
        if (abs(sum(ps54(:,iy,imn,1))+outside54(iy,imn)-monthly) > 1.D-10*max(scale,1.D-30)) stop '54 budget mismatch'
        write(*,'(I5,I4,4F12.5)') iyear,imn,monthly,sum(ps54(:,iy,imn,1)),outside54(iy,imn),summ
      end do
      print '(A,I5,A,F14.6)', 'Year ', iyear, ' global net CH4 [Tg/year]: ', summ
    ENDDO

    ! Exact expected period: reject additional GT3 records
    read(funit,iostat=ierr,iomsg=msg) head
    if (ierr == 0) stop 'Additional GT3 records after f_year'
    if (.not.is_iostat_end(ierr)) call check_io(ierr,'Check end of GT3 input')
    close(funit,iostat=ierr,iomsg=msg)
    call check_io(ierr,'Close flux file')

    fname = uncr_file
    open(ounit,file=uncr_file,status='replace',action='write',iostat=ierr,iomsg=msg)
    call check_io(ierr,'Open monthly regional output')
    ! Preserve the existing table layout for downstream readers
    write(ounit,'(A5,15I9)',iostat=ierr,iomsg=msg) 'reg#', (imn, imn=1,nmon)
    call check_io(ierr,'Write regional header')
    do k = 1, nr54                                   ! for FIX case seasonality
      write(ounit,'(I5,15F9.2)',iostat=ierr,iomsg=msg) k, (sum(ps54(k,:,imn,1))/dble(nyear), imn=1,nmon)
      call check_io(ierr,'Write regional monthly climatology')
    end do
    close(ounit,iostat=ierr,iomsg=msg)
    call check_io(ierr,'Close monthly regional output')

    print *, 'Long-term mean annual net CH4 [Tg/year], ', s_year, ' to ', f_year
    print '(A,F14.6)', 'Global total:          ', sum(total_all)/dble(nyear)
    print '(A,F14.6)', 'All 54 regions:        ', sum(ps54(:,:,:,1))/dble(nyear)
    print '(A,F14.6)', 'Outside 54 regions:    ', sum(outside54)/dble(nyear)
    print *, 'Output: ', trim(uncr_file)
    print *, 'Values: monthly climatological net CH4 masses [Tg/month], not uncertainty estimates.'
    print *, 'Program finish OK'
  End program mainprog
