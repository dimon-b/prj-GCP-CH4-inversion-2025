! ifort s3_ch4_grd2gt3.f90 -o s3_ch4_grd2gt3.exe

  Module pf4_ch4_grd2gt3_mod
    ! use netcdf
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none 

    ! [files, dirs]
    character(len=*),parameter    :: work_dir = '../results2025/INCA/'
    character(len=*),parameter   :: inout_dir = work_dir//'flux2d/'
    character(len=*),parameter   :: icase = 's060'
    character(len=*),parameter   :: flux_suffix = ''

    ! [parameter]
    integer,parameter            :: s_year = 1999, f_year = 2025, nyear=f_year-s_year+1
    integer,parameter            :: r8    = 8
    real(r8),parameter           :: pi    = 3.1415926535897932384626433832795_r8
    real(r8),parameter           :: r0    = 6371.D3
    integer,parameter            :: nmon=12, imax=128, jmax=64, kmax=1

    ! [Input *.nc]
    integer,parameter            :: immax = 360 !3600          ! -180:179 (lower left corner of grid)
    integer,parameter            :: jmmax = 180 !1800          !  -90: 89 (lower left corner of grid)
    real(r8),parameter           :: dmdeg = 180.D0/jmmax
    real(r8)                     :: mlon (immax), mlat (jmmax) ! longitude/latitude for EDGAR

    integer,parameter            :: ounit = 98
    integer,save,dimension(nmon) :: nday = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    integer,save,dimension(nmon) :: mday = (/16, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16/)
    integer,save,dimension(nmon) :: nhr  = (/12, 00, 12, 00, 12, 00, 12, 12, 00, 12, 00, 12/)

    ! [grid vars]
    real(r8)                     :: alon (imax)          ! longitude for GCM
    real(r8)                     :: dlon                 ! diff. of longitude for GCM grid
    real(r8)                     :: alat (0:jmax+1)      ! latitude
    real(r8)                     :: blat (jmax+1)        ! boundary latitude
    real(r8)                     :: dlat (jmax)          ! latitude weight
    real(r8)                     :: gar_t42(jmax)        ! grid size jmax [m**2]
    real(r8)                     :: gar_inp(jmmax)       ! grid size jmmax [m**2]
    integer                      :: iloc (immax)         !
    integer                      :: jloc (jmmax)         !

    ! [flux vars]
    integer,parameter            :: ninv=1, nunp=3, nunx=3
    character(len=25),parameter  :: invcases(ninv) = ['inv1']
    character(len=25),parameter  :: flxcases(ninv) = ['fch4_gcp2025_prior_inv1']
    character(len=15),parameter  :: unxcases(nunx) = ['ctl', 'ux2', 'ux4']
    character(len=15),parameter  :: unpcases(nunp) = ['p30', 'p50', 'p99']

    ! [gt3 writing]
    integer                      :: jdate(6)  !! date (year, month, day, hour, min., sec.)
    integer                      :: itime
    character(len=*), parameter  :: hunit = 'kg-CH4/m2/s'
    real*4                       :: flux(immax, jmmax, nmon*nyear), soil(immax, jmmax, nmon*nyear)
    real(r8)                     :: gtots (nyear)
    real*4                       :: oval (imax,jmax,1)
    character*16                 :: spcnm                ! GTOOL3 header
    character*255                :: inp_dir, out_dir, inp_file, out_file
    character                    :: invc*5, unp*5, unx*5
    integer                      :: i, j, k, l, m, iu, ic, iy, jp, jx

    real(r8)                     :: wlon(imax,immax), wlat(jmax,jmmax)
    real(r8)                     :: sedge(jmax+1)
    integer                      :: ierr
    character*512                :: io_file, msg

  contains      

    !*****************************************************************
    !*** setgrid
    Subroutine setgrid_actm()
!     [internal]
      real(r8) :: ctheta(jmax)        !! Gauss Latitude (theta)
      real(r8) :: gw(jmax)            !! Gauss weight
      integer :: i,j
      do i = 1, imax
        alon(i) = 360.D0*(i-1)/imax
      end do
      dlon = 360.D0/imax
      call gauss(ctheta,gw)
      dlat = gw
      alat(0) = 90.D0
      alat(jmax+1) = -90.D0
      sedge(1) = 1.D0
      do j = 1, jmax
        alat(j) = asin(ctheta(j))*180.D0/pi
        sedge(j+1) = sedge(j)-gw(j)
      end do
      sedge(jmax+1) = -1.D0
      blat = asin(max(-1.D0,min(1.D0,sedge)))*180.D0/pi
      return
      ! print*, alon
      ! print*, alat
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
        x0 = cos(pi*(dble(j)-0.25D0)/(dble(jmax)+0.5D0))
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
    !*** in/out settings
    Subroutine set_inout()
!     [internal]
      integer :: i,j,ii,jj,shift
      real(r8) :: dx,sl,sr,tl,tr,overlap,lo,hi
      dx = 360.D0/immax
      do i = 1, immax
        mlon(i) = -180.D0+(dble(i)-0.5D0)*dx
      end do
      do j = 1, jmmax
        mlat(j) = -90.D0+(dble(j)-0.5D0)*dmdeg
        lo = sin((-90.D0+dble(j-1)*dmdeg)*pi/180.D0)
        hi = sin((-90.D0+dble(j)*dmdeg)*pi/180.D0)
        gar_inp(j) = r0**2*(dx*pi/180.D0)*(hi-lo)
      end do
      gar_t42 = 2.D0*pi*r0**2*dlat/imax
      print '(A,ES20.12)', 'Area T42   [m2] = ',sum(gar_t42)*imax
      print '(A,ES20.12)', 'Area input [m2] = ',sum(gar_inp)*immax
      ! Longitude overlap includes periodic copies around the dateline.
      wlon = 0.D0
      do ii = 1, imax
        tl = alon(ii)-0.5D0*dlon
        tr = alon(ii)+0.5D0*dlon
        do i = 1, immax
          do shift = -1, 1
            sl = mlon(i)-0.5D0*dx+360.D0*shift
            sr = mlon(i)+0.5D0*dx+360.D0*shift
            overlap = max(0.D0,min(sr,tr)-max(sl,tl))
            wlon(ii,i) = wlon(ii,i)+overlap/dlon
          end do
        end do
      end do
      ! Latitude overlap is measured in sin(latitude), proportional to area.
      wlat = 0.D0
      do jj = 1, jmax
        do j = 1, jmmax
          lo = sin((-90.D0+dble(j-1)*dmdeg)*pi/180.D0)
          hi = sin((-90.D0+dble(j)*dmdeg)*pi/180.D0)
          overlap = max(0.D0,min(hi,sedge(jj))-max(lo,sedge(jj+1)))
          wlat(jj,j) = overlap/dlat(jj)
        end do
      end do
      if (maxval(abs(sum(wlon,dim=2)-1.D0)) > 1.D-12) stop 'Longitude coverage mismatch'
      if (maxval(abs(sum(wlat,dim=2)-1.D0)) > 1.D-12) stop 'Latitude coverage mismatch'
      ! where (mlon > 0.0)
          ! mlon = mlon
      ! elsewhere
          ! mlon = mlon + 360.0
      ! end where
      !print*, 'mlon swap', mlon
      !print *, 'iLoc ', iloc; print *, 'jLoc ', jloc
    End subroutine set_inout

    !*****************************************************************
    !*** write to gt3
    Subroutine write_t42()
      ! real*8                                 :: budg (immax, jmmax, nmon*nyear)
      real(r8) :: vol_inp,vol_t42,month_inp,month_t42,abs_inp,secs,relerr
      real(r8) :: tmp(imax,jmmax),gdz(imax,jmax)
      real*4 :: rdz(imax,jmax)
      integer :: kk,days
      character*255 :: ayear,fname
      ! re-map at t42 grid
      fname = trim(out_file)
      io_file = fname
      print '(2A)', ' Write output ::',trim(fname)
      open(ounit,file=fname,form='unformatted',convert='big_endian',status='replace',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open GT3 output')
      ! loop yearly
      do iy = 0, nyear-1
        write(ayear,'(I4.4)') s_year+iy
        vol_inp = 0.D0
        vol_t42 = 0.D0
        do k = 1, nmon
          kk = k+iy*nmon
          days = month_days(s_year+iy,k)
          secs = dble(days)*86400.D0
          tmp = matmul(wlon,dble(flux(:,:,kk)))
          gdz = matmul(tmp,transpose(wlat))
          month_inp = 0.D0
          abs_inp = 0.D0
          do j = 1, jmmax
            month_inp = month_inp+sum(dble(flux(:,j,kk)))*gar_inp(j)
            abs_inp = abs_inp+sum(abs(dble(flux(:,j,kk))))*gar_inp(j)
          end do
          month_t42 = 0.D0
          do j = 1, jmax
            month_t42 = month_t42+sum(gdz(:,j))*gar_t42(j)
          end do
          if (abs(month_inp-month_t42) > 1.D-11*max(abs_inp,1.D-30)) stop 'Conservative remap budget mismatch'
          rdz = real(gdz)
          if (.not.all(ieee_is_finite(rdz))) stop 'Non-finite output flux'
          month_t42 = 0.D0
          do j = 1, jmax
            month_t42 = month_t42+sum(dble(rdz(:,j)))*gar_t42(j)
          end do
          if (abs(month_inp-month_t42) > 2.D-7*max(abs_inp,1.D-30)) stop 'UR4 rounding budget mismatch'
          vol_inp = vol_inp+month_inp*secs*1.D-9
          vol_t42 = vol_t42+month_t42*secs*1.D-9
          itime = kk-1
          jdate = (/s_year+iy,k,1+days/2,mod(days,2)*12,0,0/)
          oval(:,:,1) = rdz
          call gt3write()
        end do
        gtots(iy+1) = vol_inp
        relerr = 0.D0
        if (abs(vol_inp) > 1.D-30) relerr = (vol_t42-vol_inp)/abs(vol_inp)*100.D0
        print '(A,A4,2F14.6,ES14.5)', '   Year Input/T42 [Tg CH4/year], difference [%]: ', &
          ayear(1:4),vol_inp,vol_t42,relerr
      end do
      close(ounit,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close GT3 output')
          !print*, iy, k, kk, mday(k), nhr(k), rdz (imax/2, 10, kk)
          !stop
    End subroutine write_t42
    
    !*****************************************************************
    !*** gt3write
    Subroutine gt3write()
!     [internal]
      character*16                         :: head(64)    ! GTOOL3 header

      ! make GTOOL3 header
      !print '(A, 2i5, 5i3, 2x, A)', '     gt3write info :: ', itime, jdate, hunit
      !print '(8i3)', itime, jdate

      head(:) = '' ! reset all headers

      write (head(1), '(I16)') 9010
      head(2) = 'GCP2025 post'           !FOS only
      head(3) = 'QFLUX07'      
      write (head(12), '(I16)') 0
      write (head(13), '(I16)') 1
      head(14) = trim(spcnm)
      head(16) = hunit
      write (head(25), '(I16)') itime
      head(26) = 'MONTH'
      write (head(27), 900) (jdate(i),i=1,6)
      write (head(28), '(I16)') 1
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

      write (ounit,iostat=ierr,iomsg=msg) head !; print *, ' header written, oval next'
      call check_io(ierr,'Write GT3 header')
      write (ounit,iostat=ierr,iomsg=msg) ((oval(i,j,1), i = 1, imax), j = 1, jmax)
      call check_io(ierr,'Write GT3 field')

  900 format(I4.4,I2.2,i2.2,' ',I2.2,I2.2,I2.2)

      return
    End subroutine gt3write
    
    !*****************************************************************
    !*** get_gtotal
    Subroutine get_gtotal()
      real(r8) :: sum1
      integer :: kk
      !temis in kg-ch4/m2/s
      do l = 1, nyear
        sum1 = 0.D0
        do m = 1, nmon
          kk = m+(l-1)*nmon
          do j = 1, jmmax
            sum1 = sum1+sum(dble(flux(:,j,kk)))*gar_inp(j)*86400.D0*month_days(s_year+l-1,m)
          end do
        end do
        gtots(l) = sum1*1.D-9 ! Tg CH4/year
      end do
        !print*, gtots (idx, l), idx, l, sum1
    End subroutine get_gtotal
    
    
    !*****************************************************************
    !*** read grd
    Subroutine read_grd()
      real*4                                 :: temis(immax, jmmax)
      integer                                :: irec, kk, nrecl
    
      print '(2A)', ' Read input   ::', trim(inp_file)
      io_file = inp_file
      irec = 0
      flux = 0.0
      inquire(iolength=nrecl) temis
      open (50,file=inp_file,form='unformatted',access='direct',status='old', recl=nrecl,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open GRD input')
      do l = 1, nyear
        do m = 1, nmon
          temis=0.0
          irec = irec + 1 
          read (50,rec=irec,iostat=ierr,iomsg=msg) ((temis(i,j), i = 1, immax), j = 1, jmmax)
          if (ierr /= 0) print *, 'Record: ',irec
          call check_io(ierr,'Read GRD record')
          if (.not.all(ieee_is_finite(temis))) stop 'Non-finite GRD flux'
          if (any(temis == -999.0)) stop 'Missing GRD flux (-999)'
          do i = 1, immax
          do j = 1, jmmax
            kk = m + (l-1)*nmon
            flux(i,j,kk) = flux(i,j,kk) + temis(i,j)!*1.0e-6/86400.0/nday(m)  !emis in Kg-CH4/mon
            ! print*, flux(i,j,kk)
          end do
          end do   ! annual flux/12.0 + monthly flux
        end do
      end do
      close(50,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close GRD input')
    End subroutine read_grd

    !*****************************************************************
    !*** read grd soil
    Subroutine read_grd_soil()
      real*4                                 :: temis(immax, jmmax)
      integer                                :: irec, kk, nrecl
      character*255                          :: inp_file_soil

      inp_file_soil = work_dir//'priors/' // 'fch4_gcp2025_soil0_inv1.grd'
      print '(2A)', ' Read input   ::', trim(inp_file_soil)
      io_file = inp_file_soil
      irec = 0
      soil = 0.0
      inquire(iolength=nrecl) temis
      open (50,file=inp_file_soil,form='unformatted',access='direct',status='old', recl=nrecl,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open GRD input')
      do l = 1, nyear
        do m = 1, nmon
          temis=0.0
          irec = irec + 1 
          read (50,rec=irec,iostat=ierr,iomsg=msg) ((temis(i,j), i = 1, immax), j = 1, jmmax)
          if (ierr /= 0) print *, 'Record: ',irec
          call check_io(ierr,'Read GRD record')
          if (.not.all(ieee_is_finite(temis))) stop 'Non-finite GRD flux'
          if (any(temis == -999.0)) stop 'Missing GRD flux (-999)'
          do i = 1, immax
          do j = 1, jmmax
            kk = m + (l-1)*nmon
            soil(i,j,kk) = soil(i,j,kk) + temis(i,j)
            ! print*, soil(i,j,kk)
          end do
          end do
        end do
      end do
      close(50,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close GRD input')
    End subroutine read_grd_soil    
    !*****************************************************************
    integer Function month_days(year,month) result(days)
      integer,intent(in) :: year,month
      days = nday(month)
      if (month == 2) then
        if (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) days = 29
      end if
    End function month_days

    Subroutine check_io(status,action)
      integer,intent(in) :: status
      character(len=*),intent(in) :: action
      if (status == 0) return
      print '(3a)', trim(action), ': ', trim(io_file)
      print '(a)', trim(msg)
      stop 1
    End subroutine check_io
  End module pf4_ch4_grd2gt3_mod
!
!*****************************************************************
!
  Program mainprog
    use pf4_ch4_grd2gt3_mod
    implicit none
    

    integer                                 :: outun

    ! dirs
    inp_dir = inout_dir
    out_dir = inout_dir

    call setgrid_actm ()
    call set_inout ()
    !call read_grd_soil ()
    
    outun = 50 !6
    
    ! scenarious
    DO jp = 1, nunp
      unp = unpcases(jp)
      DO jx = 1, nunx
        unx = unxcases(jx)
          DO ic = 1, ninv, 1
            invc = invcases(ic)
            
            inp_file = trim(inp_dir) // trim(unp) // '/'//icase//'_'// trim(unx) //'_' // trim(invc) // flux_suffix // '.grd'
            out_file = trim(out_dir) // trim(unp) // '/'//icase//'_'// trim(unx) //'_' // trim(invc) // flux_suffix // '.gt3'
            ! print*, 'inp_file ', trim(inp_file)
            ! print*, 'out_file ', trim(out_file)
            spcnm = 'fch4_'//trim(invc)
            call read_grd
            !flux = flux + soil
            call write_t42 ()
            ! call get_gtotal()
            
            ! check g-total
            ! if (outun .ne. 6) &
              ! open (outun , file='./budget'//'budget_ch4_gcp21_'//trim(invc)//'.txt', status = 'unknown')
            ! write (outun, '(a6,a9,a20,i3)') ' Flux: ', invc, ' Componets number: ', ntrac
            ! write (outun, '(a5,5x,29a8)') 'Year', (invc(itr)(6:12), itr=1,ntrac), 'Total'
            ! do iy = 1, nyear
              ! write (outun, '(29f8.2)')  1999.5+iy, (gtots(itr,iy), itr=1,ntrac), sum(gtots(:,iy))
            ! end do
            ! if (outun .ne. 6) &
              ! close (50)
        ENDDO
      ENDDO
    ENDDO

    print *, 'Program finish OK'
  End program mainprog