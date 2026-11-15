!  converting GCP 2025 emission to ACTM grids, Dmitry Belikov, 2025/09/17
!  old: ifort f0_fch4_gcp25_nc2gt3_t42.f90 -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
!       ifort f0_fch4_gcp25_nc2gt3_t42.f90 $(nf-config --fflags) $(nf-config --flibs)


  Module f0_fch4_gcp25_nc2gt3_t42_mod
    use netcdf
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none 

    ! [input/output files]
    ! --- Runtime paths and input file, prior GCP flux common to both OH cases.
    character(len=1024)           :: inp_dir, out_dir, input_nc, config_file
    character(len=32)             :: hcase, input_case
    integer                       :: cfg_s_year, cfg_f_year, cfg_ios
    namelist /actm_f0/ hcase, inp_dir, out_dir, input_nc, input_case, cfg_s_year, cfg_f_year
    integer,parameter             :: outun   = 6

    ! [parameter]
    integer,parameter             :: s_year = 1999, f_year = 2025, nyear = f_year-s_year+1
    integer,parameter             :: r8    = 8
    real(r8),parameter            :: pi    = 3.1415926535897932384626433832795_r8
    real(r8),parameter            :: r0    = 6371.D3
    real(r8),parameter            :: coeff = 86400.D0*1.D-9                  ! to convert to Tg/yr in print output
    integer,parameter             :: nmon  = 12
    integer,parameter             :: imax  = 128, jmax = 64, kmax = 1, ndims = 3

    ! [Input *.nc]
    integer,parameter             :: immax = 360 !3600          ! -180:179 (lower left corner of grid)
    integer,parameter             :: jmmax = 180 !1800          !  -90: 89 (lower left corner of grid)
    real(r8),parameter            :: dmdeg = 180.D0/jmmax
    real(r8)                      :: mlon (immax), mlat (jmmax) ! longitude/latitude for EDGAR

    integer,parameter             :: ounit = 98
    integer,save,dimension(nmon)  :: nday = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    integer,save,dimension(nmon)  :: mday = (/16, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16/)
    integer,save,dimension(nmon)  :: nhr  = (/12, 00, 12, 00, 12, 00, 12, 12, 00, 12, 00, 12/)

    ! [grid vars]
    real(r8)                      :: alon (imax)          ! longitude for GCM
    real(r8)                      :: dlon                 ! diff. of longitude for GCM grid
    real(r8)                      :: alat (0:jmax+1)      ! latitude
    real(r8)                      :: blat (jmax+1)        ! boundary latitude
    real(r8)                      :: dlat (jmax)          ! latitude weight
    real(r8)                      :: gar_t42(jmax)        ! grid size jmax [m**2]
    real(r8)                      :: gar_inp(jmmax)       ! grid size jmmax [m**2]
    integer                       :: iloc (immax)         !
    integer                       :: jloc (jmmax)         !
    real(r8)                      :: wlon(imax, immax), wlat(jmax, jmmax)

    ! [gt3 writing]
    integer                       :: jdate(6)  !! date (year, month, day, hour, min., sec.)
    integer                       :: itime
    character(len=*), parameter   :: hunit = 'kg-CH4/m2/s'
    real*4                        :: budg_mn(immax, jmmax, nmon), temp_mn(immax, jmmax, nmon)
    real*4                        :: budg_yr(immax, jmmax, nmon*nyear), temp_yr(immax, jmmax, nmon*nyear)
    real*4                        :: budg_vz(immax, jmmax, nmon*nyear), temp_vz(immax, jmmax, nmon*nyear)
    real*4                        :: temp_1(jmmax, immax, nmon*nyear),  temp_slice(immax, jmmax)
    real(r8)                      :: gtots (20, nyear)
    real*4                        :: oval (imax,jmax,1)
    character*25                  :: spcnm                ! GTOOL3 header
    character*255                 :: viz_dir, inp_file, out_file
    integer                       :: i, j, k, l, m!, iy

  contains

    !*****************************************************************
    !*** Runtime config (GT3 calculations are unchanged)
    Subroutine read_config()
      hcase = ''
      inp_dir = ''
      out_dir = ''
      input_nc = ''
      input_case = ''
      config_file = ''
      cfg_s_year = s_year
      cfg_f_year = f_year
      call get_command_argument(1,config_file)
      if (len_trim(config_file) == 0) error stop 'Usage: f0_fch4_gcp25_nc2gt3_t42 <actm_f0.nml>'
      open(10,file=trim(config_file),status='old',action='read',iostat=cfg_ios)
      if (cfg_ios /= 0) error stop 'Cannot open F0 namelist'
      read(10,nml=actm_f0,iostat=cfg_ios)
      close(10)
      if (cfg_ios /= 0) error stop 'Invalid F0 namelist'
      if (trim(hcase) /= 'CYC' .and. trim(hcase) /= 'INCA') error stop 'Invalid F0 case'
      if (cfg_s_year /= s_year .or. cfg_f_year /= f_year) error stop 'F0 period must remain 1999-2025'
      if (len_trim(inp_dir) == 0 .or. len_trim(out_dir) == 0 .or. &
          len_trim(input_nc) == 0 .or. len_trim(input_case) == 0) error stop 'F0 path is empty'
      if (trim(input_case) /= 'inv1') error stop 'F0 expected inversion prior inv1'
      print '(2A)', '      Case        : ', trim(hcase)
      print '(2A)', '      Input NC    : ', trim(inp_dir)//trim(input_nc)
      print '(2A)', '      Output dir  : ', trim(out_dir)
    End subroutine read_config

    !*****************************************************************
    !*** setgrid
    Subroutine setgrid_actm()
!     [internal]
      real(r8) :: ctheta(jmax), gw(jmax), edge
      integer  :: i, j

      do i = 1, imax
        alon(i) = 360.D0 * (i-1) / imax
      end do
      dlon = 360.D0 / imax
      call gauss(ctheta, gw)
      dlat = gw
      alat(0) = 90.D0
      alat(jmax+1) = -90.D0
      alat(1:jmax) = asin(ctheta) * 180.D0 / pi

      ! Gaussian cell boundaries from cumulative weights; sum(gw) = 2
      edge = 1.D0
      blat(1) = 90.D0
      do j = 1, jmax
        edge = edge - gw(j)
        blat(j+1) = asin(max(-1.D0, min(1.D0, edge))) * 180.D0 / pi
      end do
      blat(jmax+1) = -90.D0
      return

      ! [retained original comments and commented code]
      !print*, alon;      print*, alat  stop'alat'
    End subroutine setgrid_actm
  
    !*****************************************************************
    !*** set gauss grid
    Subroutine gauss(ctheta, gw)
    ! [output]
      real(r8), intent(out) :: ctheta(jmax), gw(jmax)
    ! [internal parm]
      real(r8), parameter   :: deltp = 1.D-14
      integer, parameter    :: itrmax = 100
    ! [internal work]
      integer               :: n, j, iter
      real(r8)              :: x0, x1, p0, p1, pn, dp

      do j = 1, (jmax+1)/2
        x0 = cos(pi * (real(j,r8)-0.25D0) / (real(jmax,r8)+0.5D0))
        do iter = 1, itrmax
          p0 = 1.D0
          p1 = x0
          do n = 2, jmax
            pn = ((2.D0*n-1.D0)*x0*p1 - (n-1.D0)*p0) / n
            p0 = p1
            p1 = pn
          end do
          dp = jmax * (x0*p1-p0) / (x0*x0-1.D0)
          x1 = x0 - p1/dp
          if (abs(x1-x0) < deltp) exit
          x0 = x1
        end do
        if (iter > itrmax) stop 'Gaussian latitude iteration failed'
        x0 = x1

        ! Recalculate derivative at the converged root
        p0 = 1.D0
        p1 = x0
        do n = 2, jmax
          pn = ((2.D0*n-1.D0)*x0*p1 - (n-1.D0)*p0) / n
          p0 = p1
          p1 = pn
        end do
        dp = jmax * (x0*p1-p0) / (x0*x0-1.D0)
        ctheta(j) = x0
        ctheta(jmax-j+1) = -x0
        gw(j) = 2.D0 / ((1.D0-x0*x0)*dp*dp)
        gw(jmax-j+1) = gw(j)
      end do
      if (abs(sum(gw)-2.D0) > 1.D-12) stop 'Invalid Gaussian weights'
      return

      ! [retained original comments and commented code]
    !
    !
    ! [work]
    !
    !
        ! step 1
          ! step 2 - Pn
          ! step 3 - d/d Pn
          ! step 4 - newton
        ! step 5
    End subroutine gauss

    !*****************************************************************
    !*** in/out settings
    Subroutine set_inout()
!     [internal]
      integer  :: i, j, ii, jj, shift
      real(r8) :: dx, dy, west, east, south, north, left, right, low, high, overlap
      real(r8) :: target, tol

      if (.not.all(ieee_is_finite(mlon))) stop 'Non-finite input longitude'
      if (.not.all(ieee_is_finite(mlat))) stop 'Non-finite input latitude'
      dx = 360.D0 / immax
      dy = 180.D0 / jmmax
      tol = 1.D-8

      ! Input must be a complete regular global cell-centred grid
      do i = 1, immax
        ii = mod(i,immax) + 1
        if (abs(modulo(mlon(ii)-mlon(i),360.D0)-dx) > tol) stop 'Invalid longitude spacing'
      end do
      do j = 2, jmmax
        if (abs(abs(mlat(j)-mlat(j-1))-dy) > tol) stop 'Invalid latitude spacing'
        if ((mlat(j)-mlat(j-1))*(mlat(2)-mlat(1)) <= 0.D0) stop 'Non-monotonic latitude'
      end do
      if (abs(minval(mlat)+90.D0-dy/2.D0) > tol) stop 'Invalid south latitude edge'
      if (abs(maxval(mlat)-90.D0+dy/2.D0) > tol) stop 'Invalid north latitude edge'
      mlon = modulo(mlon,360.D0)

      do j = 1, jmax
        gar_t42(j) = 2.D0 * pi * r0**2 * dlat(j) / imax
      end do
      do j = 1, jmmax
        south = (mlat(j)-dy/2.D0)*pi/180.D0
        north = (mlat(j)+dy/2.D0)*pi/180.D0
        gar_inp(j) = r0**2 * dx*pi/180.D0 * (sin(north)-sin(south))
      end do
      print '(A, f12.4)', 'Area T42   = ', sum(gar_t42)*imax*1.D-12
      print '(A, f12.4)', 'Area input = ', sum(gar_inp)*immax*1.D-12

      ! Longitude overlap [radians], including the cyclic boundary
      wlon = 0.D0
      do i = 1, immax
        west = mlon(i)-dx/2.D0
        east = mlon(i)+dx/2.D0
        do ii = 1, imax
          do shift = -1, 1
            left = alon(ii)-dlon/2.D0+360.D0*shift
            right = alon(ii)+dlon/2.D0+360.D0*shift
            overlap = max(0.D0, min(east,right)-max(west,left))
            wlon(ii,i) = wlon(ii,i) + overlap*pi/180.D0
          end do
        end do
        iloc(i) = maxloc(wlon(:,i),dim=1)
      end do

      ! Latitude overlap in sin(latitude); supports either latitude order
      wlat = 0.D0
      do j = 1, jmmax
        south = mlat(j)-dy/2.D0
        north = mlat(j)+dy/2.D0
        do jj = 1, jmax
          low = max(south,blat(jj+1))
          high = min(north,blat(jj))
          if (high > low) wlat(jj,j) = sin(high*pi/180.D0)-sin(low*pi/180.D0)
        end do
        jloc(j) = maxloc(wlat(:,j),dim=1)
      end do

      ! Check that overlap weights cover every input and output cell
      do i = 1, immax
        if (abs(sum(wlon(:,i))-dx*pi/180.D0) > tol) stop 'Input longitude coverage failed'
      end do
      do ii = 1, imax
        if (abs(sum(wlon(ii,:))-dlon*pi/180.D0) > tol) stop 'Output longitude coverage failed'
      end do
      do j = 1, jmmax
        target = gar_inp(j)/(r0**2*dx*pi/180.D0)
        if (abs(sum(wlat(:,j))-target) > tol) stop 'Input latitude coverage failed'
      end do
      do jj = 1, jmax
        if (abs(sum(wlat(jj,:))-dlat(jj)) > tol) stop 'Output latitude coverage failed'
      end do

      ! [retained original comments and commented code]
      !print*, 'mlon swap', mlon
      !print *, 'iLoc ', iloc; print *, 'jLoc ', jloc
    End subroutine set_inout

    !*****************************************************************
    !*** read_nc with 12 month records
    Subroutine read_nc_12mon()
      print '(2A)', 'Read input   ::', trim(inp_file)
      call read_flux_nc(trim(inp_file), 'fch4', 'longitude', 'latitude', nmon, temp_mn)

      ! [retained original comments and commented code]
      !print *, ndims_in, nvars_in, ngatts_in, unlimdimid_in
      !return fch4
    End subroutine read_nc_12mon

    !*****************************************************************
    !*** read_nc with yearly records
    Subroutine read_nc_yr()
      print '(2A)', 'Read input   ::', trim(inp_file)
      call read_flux_nc(trim(inp_file), 'flux_ch4_prior', 'lon', 'lat', nmon*nyear, temp_yr)

      ! [retained original comments and commented code]
      ! call set_inout () ! mlon, mlat -> iloc, jloc
      !! Rearrange data from temp_1 to temp_yr using transpose
      !do k = 1, nmon*nyear
      !   temp_slice = transpose(temp_1(:,:,k))
      !   temp_yr(:,:,k) = temp_slice
      !end do     
      !
      !! Convert units "Kg-CH4/m2/mon" -> "Kg-CH4/m2/s"
      !do l = 1, nyear
      !do m = 1, nmon
      !  kk = m + (l-1)*nmon
      !  temp_yr(:,:,kk) = temp_yr(:,:,kk) / nday(m) / 86400.0
      !end do
      !end do  
    End subroutine read_nc_yr
    
    !*****************************************************************
    !*** check
    Subroutine check(status)
      integer, intent ( in) :: status
      
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 1
      end if
    End subroutine check 

    !*****************************************************************
    !*** write to gt3
    Subroutine write_t42()
      real(r8)                               :: vol_inp, vol_t42, rate_inp, rate_t42, scale
      real(r8)                               :: work(imax,jmmax), gdz(imax,jmax)
      real*4                                 :: rdz(imax,jmax)
      integer                                :: i, j, l, m, kk, days, ios
      character*255                          :: ayear, fname, msg

      fname = trim(out_file)//'t42.gt3'
      print '(2A)', ' Write output ::', trim(fname)
      open(unit=ounit, file=fname, form='unformatted', status='replace', action='write', &
           convert='big_endian', iostat=ios, iomsg=msg)
      if (ios /= 0) then
        print *, trim(msg)
        stop 'Cannot open GT3 output'
      end if

      ! loop yearly
      do l = 0, nyear-1
        write(ayear,'(i4.4)') s_year+l
        vol_inp = 0.D0
        vol_t42 = 0.D0
        do m = 1, nmon
          kk = m + l*nmon
          days = month_days(s_year+l,m)

          ! Conservative remapping by cell overlap; preserve negative soil uptake
          work = matmul(wlon,real(budg_yr(:,:,kk),r8))
          gdz = r0**2 * matmul(work,transpose(wlat))
          do j = 1, jmax
            gdz(:,j) = gdz(:,j) / gar_t42(j)
          end do
          if (.not.all(ieee_is_finite(gdz))) stop 'Non-finite remapped flux'
          rdz = real(gdz,kind=kind(rdz))

          rate_inp = 0.D0
          rate_t42 = 0.D0
          scale = 0.D0
          do j = 1, jmmax
            rate_inp = rate_inp + sum(real(budg_yr(:,j,kk),r8))*gar_inp(j)
            scale = scale + sum(abs(real(budg_yr(:,j,kk),r8)))*gar_inp(j)
          end do
          do j = 1, jmax
            rate_t42 = rate_t42 + sum(real(rdz(:,j),r8))*gar_t42(j)
          end do
          if (abs(rate_t42-rate_inp) > 1.D-6*max(scale,1.D-30)) stop 'Flux conservation check failed'
          vol_inp = vol_inp + rate_inp*days*coeff
          vol_t42 = vol_t42 + rate_t42*days*coeff

          itime = kk-1
          jdate = (/s_year+l, m, 1+days/2, mod(days,2)*12, 0, 0/)
          oval(:,:,1) = rdz
          call gt3write()
        end do
        print '(A, A5, 2f12.4)', '   Year Grid Inp/T42 :: ', ayear, vol_inp, vol_t42
      end do
      close(unit=ounit,iostat=ios,iomsg=msg)
      if (ios /= 0) then
        print *, trim(msg)
        stop 'Cannot close GT3 output'
      end if

      ! [retained original comments and commented code]
      ! total budget
      ! re-map at t42 grid
      ! negative numbers from interpolation
        ! print*, ''
        !print *, 1 + l*nmon, (l+1)*nmon
    End subroutine write_t42
    
    !*****************************************************************
    !*** gt3write
    Subroutine gt3write()
!     [internal]
      integer                              :: i, j, ios
      character*255                        :: msg
      character*16                         :: head(64)    ! GTOOL3 header

      ! make GTOOL3 header
      print '(A, 2i5, 5i3, 2x, A)', '     gt3write info :: ', itime, jdate, hunit
      !print '(8i3)', itime, jdate

      head(:) = '' ! reset all headers
      write (head(1), '(I16)') 9010
      head(2) = 'GCP-Global'           !FOS only
      head(3) = 'QFLUX07'
      write (head(12), '(I16)') 0
      write (head(13), '(I16)') 1
      !head(14) = 'tracer  ' !trim(spcnm)
      head(14) = trim(spcnm)
      head(16) = hunit
      write (head(25), '(I16)') itime
      head(26) = 'MONTH'
      write (head(27), 900) (jdate(i),i=1,6)
      write (head(28), '(I16)') 0
      if (imax.lt.100) then    ! X-axis
         write (head(29), '(A,I2)') 'GLON',imax
      else
         write (head(29), '(A,I3)') 'GLON',imax
      end if
      write (head(30), '(I16)') 1
      write (head(31), '(I16)') imax
      if (jmax.lt.100) then    ! Y-axis
         write (head(32), '(A,I2)') 'GGLA',jmax
      else
         write (head(32), '(A,I3)') 'GGLA',jmax
      end if
      write (head(33), '(I16)') 1
      write (head(34), '(I16)') jmax
      if (kmax.eq.32) then    ! for l32(troposphere)
         write (head(35), '(A,I2.2,A)') 'CSIG',kmax,'t'
      else if (kmax.eq.1) then  !! for surface level
         write (head(35), '(A)') 'SFC1'
      else
         write (head(35), '(A,I2.2)') 'CSIG',kmax
      end if
      write (head(36), '(I16)') 1
      write (head(37), '(I16)') kmax
      head(38) = 'UR4'
      write (head(39), '(E16.7)') -0.9990000e+03
      write (head(40), '(E16.7)') -0.9990000e+03
      write (head(41), '(E16.7)') -0.9990000e+03
      write (head(42), '(E16.7)') -0.9990000e+03
      write (head(43), '(E16.7)') -0.9990000e+03
      write (head(44), '(I16)'  ) 1
      write (head(46), '(I16)'  ) 0
      write (head(47), '(E16.7)') 0.
      write (head(48), '(I16)'  ) 0
      write (head(63), '(A)'    ) 'Prabir_Patra'
      write (head(64), '(I16)'  ) imax*jmax*kmax

      write(ounit,iostat=ios,iomsg=msg) head !; print *, ' header written, oval next'
      if (ios /= 0) then
        print *, trim(msg)
        stop 'Cannot write GT3 header'
      end if
      write(ounit,iostat=ios,iomsg=msg) ((oval(i,j,1), i = 1, imax), j = 1, jmax)
      if (ios /= 0) then
        print *, trim(msg)
        stop 'Cannot write GT3 field'
      end if

  900 format(I4.4,I2.2,i2.2,' ',I2.2,I2.2,I2.2)

      return
    End subroutine gt3write
    
    !*****************************************************************
    !*** get_gtotal
    Subroutine get_gtotal(idx, cyl)
      real(r8)                               :: sum1
      integer, intent(in)                    :: idx, cyl
      integer                                :: j, l, m, kk, days

      if (idx < 1 .or. idx > size(gtots,1)) stop 'Invalid budget index'
      if (cyl /= 0 .and. cyl /= 1) stop 'Invalid climatology flag'
      ! emis in kg-ch4/m2/s
      do l = 1, nyear
        sum1 = 0.D0
        do m = 1, nmon
          kk = m + (l-1)*nmon
          days = month_days(s_year+l-1,m)
          do j = 1, jmmax
            if (cyl == 1) then
              sum1 = sum1 + sum(real(temp_mn(:,j,m),r8))*gar_inp(j)*days
            else
              sum1 = sum1 + sum(real(budg_yr(:,j,kk),r8))*gar_inp(j)*days
            end if
          end do
        end do
        gtots(idx,l) = sum1*coeff
      end do

      ! [retained original comments and commented code]
            ! cyl or iav fluxes
            !if (cyl .eq. 1) then
            !  kk = m
            !  sum1 = sum1 + temp_mn(i,j,kk)*cellarea * nday(m)
            !  !sum1 = sum1 + temp_mn(i,j,kk)*gar_inp(j) * 86400.0 * mday(m)
            !  !print*, idx, l, m, i, j, temp_yr(i,j,kk), sum1
            !else
            !  kk = m + (l-1)*nmon
            !  sum1 = sum1 + temp_yr(i,j,kk)*cellarea * nday(m)
            !  !sum1 = sum1 + temp_mn(i,j,kk)*gar_inp(j) * 86400.0 * mday(m)
            !endif
        !if (idx.eq.11) 
        !print*, gtots (idx, l), idx, l, sum1*coeff
    End subroutine get_gtotal
    
    
    !*****************************************************************
    !*** visit_flx Cao/WH
    Subroutine read_visit_flx()
      real*4                                 :: temis(immax,jmmax)
      integer                                :: irec, kk, l, m, nbytes, ios, days
      integer, parameter                     :: visit_s_year = 1950
      character*255                          :: msg

      print '(2A)', 'Read input   ::', trim(inp_file)
      ! Assumes the source starts in January 1950; set visit_s_year to the actual start
      if (s_year < visit_s_year) stop 'VISIT start year precedes source'
      irec = (s_year-visit_s_year)*nmon
      inquire(iolength=nbytes) temis
      open(50,file=inp_file,form='unformatted',access='direct',status='old',recl=nbytes, &
           action='read',iostat=ios,iomsg=msg)
      if (ios /= 0) then
        print *, trim(msg)
        stop 'Cannot open VISIT input'
      end if
      temp_vz = 0.0
      do l = 1, nyear
        do m = 1, nmon
          irec = irec + 1
          read(50,rec=irec,iostat=ios,iomsg=msg) temis
          if (ios /= 0) then
            print *, trim(msg)
            stop 'Cannot read VISIT record'
          end if
          if (.not.all(ieee_is_finite(temis))) stop 'Non-finite VISIT input'
          kk = m + (l-1)*nmon
          days = month_days(s_year+l-1,m)
          ! Original conversion assumes mg-CH4/m2/month, not kg-CH4/month
          temp_vz(:,:,kk) = temis*1.0e-6/86400.0/days
        end do
      end do
      close(50,iostat=ios,iomsg=msg)
      if (ios /= 0) then
        print *, trim(msg)
        stop 'Cannot close VISIT input'
      end if
    End subroutine read_visit_flx
    
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
    !*** days_before_year
    integer Function days_before_year(year) result(days)
      integer, intent(in) :: year
      integer            :: y

      if (year < 1 .or. year > 9999) stop 'Unsupported calendar year'
      y = year-1
      days = 365*y + y/4 - y/100 + y/400
    End function days_before_year

    !*****************************************************************
    !*** check_time
    Subroutine check_time(ncid, nrec)
      integer, intent(in) :: ncid, nrec
      integer            :: varid, status, p, yr, mon, day, hour, minute, ios, j, kk, y, m
      real(r8)           :: times(nrec), factor, origin, begin_month, end_month, stamp, second
      character*255      :: units, calendar_name, origin_text

      call check(nf90_inq_varid(ncid,'time',varid))
      call check(nf90_get_var(ncid,varid,times))
      if (.not.all(ieee_is_finite(times))) stop 'Non-finite time coordinate'
      if (any(times(2:nrec) <= times(1:nrec-1))) stop 'Time must increase strictly'
      call check(nf90_get_att(ncid,varid,'units',units))
      calendar_name = 'standard'
      status = nf90_get_att(ncid,varid,'calendar',calendar_name)
      if (status == nf90_enotatt) then
        calendar_name = 'standard'
      else
        call check(status)
      end if
      select case(trim(calendar_name))
      case('standard','gregorian','proleptic_gregorian')
      case default
        stop 'Only Gregorian monthly input is supported'
      end select

      p = index(units,' since ')
      if (p == 0) stop 'Missing time reference date'
      select case(trim(units(:p-1)))
      case('days','day')
        factor = 86400.D0
      case('hours','hour')
        factor = 3600.D0
      case('minutes','minute')
        factor = 60.D0
      case('seconds','second')
        factor = 1.D0
      case default
        stop 'Unsupported time units'
      end select
      origin_text = adjustl(units(p+7:))
      if (len_trim(origin_text) < 10) stop 'Invalid time reference date'
      read(origin_text(1:4),*,iostat=ios) yr
      if (ios /= 0) stop 'Invalid reference year'
      read(origin_text(6:7),*,iostat=ios) mon
      if (ios /= 0) stop 'Invalid reference month'
      read(origin_text(9:10),*,iostat=ios) day
      if (ios /= 0) stop 'Invalid reference day'
      if (yr < 1583) stop 'Reference date must use the Gregorian calendar'
      if (day < 1 .or. day > month_days(yr,mon)) stop 'Invalid reference day'
      hour = 0
      minute = 0
      second = 0.D0
      if (len_trim(origin_text) > 10) then
        if (len_trim(origin_text) < 19) stop 'Expected reference time HH:MM:SS'
        read(origin_text(12:13),*,iostat=ios) hour
        if (ios /= 0) stop 'Invalid reference hour'
        read(origin_text(15:16),*,iostat=ios) minute
        if (ios /= 0) stop 'Invalid reference minute'
        read(origin_text(18:19),*,iostat=ios) second
        if (ios /= 0) stop 'Invalid reference second'
        if (len_trim(origin_text) > 19) then
          if (trim(origin_text(20:)) /= 'Z' .and. trim(origin_text(20:)) /= '+00:00') then
            stop 'Unsupported reference timezone or fractional second'
          end if
        end if
      end if
      if (hour < 0 .or. hour > 23 .or. minute < 0 .or. minute > 59) stop 'Invalid reference time'
      if (second < 0.D0 .or. second >= 60.D0) stop 'Invalid reference second'
      origin = real(days_before_year(yr)+day-1,r8)*86400.D0 + hour*3600.D0 + minute*60.D0 + second
      do j = 1, mon-1
        origin = origin + month_days(yr,j)*86400.D0
      end do

      ! One record must fall within each expected calendar month
      begin_month = real(days_before_year(s_year),r8)*86400.D0
      do kk = 1, nrec
        y = s_year+(kk-1)/nmon
        m = mod(kk-1,nmon)+1
        end_month = begin_month + month_days(y,m)*86400.D0
        stamp = origin + times(kk)*factor
        if (stamp < begin_month .or. stamp >= end_month) then
          print *, 'Unexpected input month at record ', kk, ' expected ', y, m
          stop 'NetCDF period does not match s_year/f_year'
        end if
        begin_month = end_month
      end do
    End subroutine check_time

    !*****************************************************************
    !*** read_flux_nc
    Subroutine read_flux_nc(path, name, lon_name, lat_name, nrec, values)
      character(len=*), intent(in) :: path, name, lon_name, lat_name
      integer, intent(in)          :: nrec
      real*4, intent(out)          :: values(immax,jmmax,nrec)
      integer                     :: ncid, varid, lonid, latid, timeid, dims(nf90_max_var_dims), rank
      integer                     :: n, kk, status, coordid, ndims_in, nvars_in, ngatts_in, unlimdimid_in
      real(r8)                    :: field(immax,jmmax,1), fill_value, missing, packing, offset
      logical                     :: has_fill, has_missing
      character*255               :: units

      call check(nf90_open(path,NF90_NOWRITE,ncid))
      call check(nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,unlimdimid_in))
      print *, ndims_in, nvars_in, ngatts_in, unlimdimid_in
      call check(nf90_inq_dimid(ncid,lon_name,lonid))
      call check(nf90_inquire_dimension(ncid,lonid,len=n))
      if (n /= immax) stop 'Input longitude size mismatch'
      call check(nf90_inq_dimid(ncid,lat_name,latid))
      call check(nf90_inquire_dimension(ncid,latid,len=n))
      if (n /= jmmax) stop 'Input latitude size mismatch'
      call check(nf90_inq_dimid(ncid,'time',timeid))
      call check(nf90_inquire_dimension(ncid,timeid,len=n))
      if (n /= nrec) stop 'Input time size mismatch'

      call check(nf90_inq_varid(ncid,lon_name,coordid))
      call check(nf90_inquire_variable(ncid,coordid,ndims=rank,dimids=dims))
      if (rank /= 1 .or. dims(1) /= lonid) stop 'Invalid longitude coordinate dimensions'
      call check(nf90_get_var(ncid,coordid,mlon))
      call check(nf90_inq_varid(ncid,lat_name,coordid))
      call check(nf90_inquire_variable(ncid,coordid,ndims=rank,dimids=dims))
      if (rank /= 1 .or. dims(1) /= latid) stop 'Invalid latitude coordinate dimensions'
      call check(nf90_get_var(ncid,coordid,mlat))
      call check(nf90_inq_varid(ncid,'time',coordid))
      call check(nf90_inquire_variable(ncid,coordid,ndims=rank,dimids=dims))
      if (rank /= 1 .or. dims(1) /= timeid) stop 'Invalid time coordinate dimensions'
      if (nrec == nmon*nyear) call check_time(ncid,nrec)

      call check(nf90_inq_varid(ncid,name,varid))
      call check(nf90_inquire_variable(ncid,varid,ndims=rank,dimids=dims))
      if (rank /= 3) stop 'Flux must have three dimensions'
      ! Fortran order lon,lat,time corresponds to Python order time,lat,lon
      if (any(dims(1:3) /= (/lonid,latid,timeid/))) stop 'Expected NetCDF flux dimensions time,lat,lon'
      units = ''
      status = nf90_get_att(ncid,varid,'units',units)
      if (status /= nf90_noerr .and. status /= nf90_enotatt) call check(status)
      if (status == nf90_noerr) then
        select case(trim(units))
        case('kg m-2 s-1','kg m^-2 s^-1','kg m**-2 s**-1','kg/m2/s','kg-CH4/m2/s','kg CH4 m-2 s-1')
        case default
          print *, 'Unrecognized flux units: ', trim(units)
          stop 'Input must be kg-CH4/m2/s'
        end select
      else
        print *, 'Flux units absent: assuming kg-CH4/m2/s from the Python conversion'
      end if

      status = nf90_get_att(ncid,varid,'_FillValue',fill_value)
      has_fill = status == nf90_noerr
      if (status /= nf90_noerr .and. status /= nf90_enotatt) call check(status)
      status = nf90_get_att(ncid,varid,'missing_value',missing)
      has_missing = status == nf90_noerr
      if (status /= nf90_noerr .and. status /= nf90_enotatt) call check(status)
      ! Set defaults after reading absent attributes, not before nf90_get_att
      status = nf90_get_att(ncid,varid,'scale_factor',packing)
      if (status == nf90_enotatt) then
        packing = 1.D0
      else
        call check(status)
      end if
      status = nf90_get_att(ncid,varid,'add_offset',offset)
      if (status == nf90_enotatt) then
        offset = 0.D0
      else
        call check(status)
      end if
      print '(A,2ES16.7)', ' scale_factor/add_offset :: ', packing, offset
      if (.not.ieee_is_finite(packing)) stop 'Non-finite scale_factor attribute'
      if (.not.ieee_is_finite(offset)) stop 'Non-finite add_offset attribute'
      do kk = 1, nrec
        call check(nf90_get_var(ncid,varid,field,start=(/1,1,kk/),count=(/immax,jmmax,1/)))
        if (.not.all(ieee_is_finite(field))) stop 'NaN or infinity in input flux'
        if (has_fill) then
          if (any(field == fill_value)) stop 'Fill values in input flux; fix sector masks before summing'
        end if
        if (has_missing) then
          if (any(field == missing)) stop 'Missing values in input flux'
        end if
        if (any(abs(field) >= real(nf90_fill_float,r8))) stop 'Default fill value or overflow in input'
        field = field*packing + offset
        if (.not.all(ieee_is_finite(field))) then
          print *, 'Unpacking failed at record ', kk, ' scale/offset = ', packing, offset
          stop 'Invalid unpacked flux'
        end if
        if (any(abs(field) > real(huge(values),r8))) stop 'Flux overflows real*4'
        values(:,:,kk) = real(field(:,:,1),kind=kind(values))
      end do
      call check(nf90_close(ncid))
    End subroutine read_flux_nc

  End module f0_fch4_gcp25_nc2gt3_t42_mod
!
!*****************************************************************
!
  Program mainprog
    use f0_fch4_gcp25_nc2gt3_t42_mod
    implicit none
   
    character                               :: tracers(20)*40, flx_viz(2)*40, run*7
    integer                                 :: jf, ji, maxiav, maxcyl, maxviz, itr, ntrac

    call read_config()
    call setgrid_actm ()

    ntrac = 1
    if (ntrac > size(gtots,1)) stop 'Too many tracers'
    gtots = 0.D0

    ! scenarious
    do ji=1, ntrac
      run = trim(input_case) 
      budg_yr = 0.0

      spcnm = trim(run)
      out_file = trim(out_dir)//'fch4_gcp25_'//trim(run)//'_'
      print*, spcnm, trim(out_file)
      print*
      
      ! fluxes
      inp_file = trim(inp_dir)//trim(input_nc)
      call read_nc_yr()

      budg_yr = budg_yr + temp_yr
      call set_inout()
      call get_gtotal(ji, 0)
      tracers(ji) = trim(input_case)

      call write_t42()
    enddo

    ! check g-total
    if (outun .ne. 6) &
      open (outun , file=trim(out_dir)//'budget_ch4_gcp25_.txt', status = 'unknown')
    write (outun, '(a16,a20,i3)') 'Flux for gcp25 ', ' Tracer number: ', ntrac
    write (outun, '(a5,5x,29a8)') 'Year', (tracers(itr), itr=1,ntrac)!, 'Total'
    do l = 1, nyear
      write (outun, '(i5,5x,29f8.2)')  s_year +l-1, (gtots(itr,l), itr=1,ntrac)
    end do
    if (outun .ne. 6) &
      close (outun)


  End program mainprog
