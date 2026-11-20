! ifort s2_ch4_grd.f90 

  Module pf3_ch4_grd_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none 

    ![input parameter]
    character(len=*),parameter   :: work_dir = '../results2025/INCA/'
    character(len=*),parameter   :: icase = 's060'
    character(len=*),parameter   :: prd_dir  = work_dir//icase//'_'
    character(len=*),parameter   :: prd_file = 'predicted.sources'
    character(len=*),parameter   :: flx_dir  = '../fluxes/gcp2025_flux_prior/'
    character(len=*),parameter   :: lsc_dir  = work_dir//'losscorr/'
    character(len=*),parameter   :: lsc_file = 'gcp_LC_'
    character(len=*),parameter   :: msk_file = '../inv_dir/mask/mask80.bin'
    character(len=*),parameter   :: lnb_file = '../inv_dir/mask/ch4_regn_pulse.grd'

    ![output parameter]
    character(len=*),parameter   :: out_dir = work_dir//'flux2d/'
    character(len=*),parameter   :: cgt_file = 'check_GT_'
    logical,parameter            :: apply_lc = .true.
    character(len=*),parameter   :: flux_suffix = ''
    
    ! [parameter]
    integer,parameter            :: s_year = 1999, f_year = 2025, nyear=f_year-s_year+1
    integer, parameter           :: nyrsk  = 0
    integer,parameter            :: r8 = 8
    real(r8),parameter           :: pi     = 3.1415926535897932384626433832795_r8
    real(r8),parameter           :: r0     = 6371.D3
    integer, parameter           :: nmon   = 12
    integer,parameter            :: nlon   = 360, nlat=180, nreg=54
    integer,save                 :: nday(nmon) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    ! [vars]
    real*4                       :: binreg(nlat,nlon)!, lnd54(nlat,nlon),lnd12(nlat,nlon)
    real(r8)                     :: tdi (nyear,nmon,nreg), time (nyear,nmon)
    real*4                       :: ocntr(nmon,nlat,nlon)
    real(r8)                     :: lc_fact(nyear), gtotal(nyear)
    real(r8)                     :: angle(nlat), cellarea(nlat)
    real(r8)                     :: totl(nreg), totc(nreg), totlc(nreg)!, totcc(nreg)
    real(r8)                     :: ofset, sum1, sum2, pstf
    character                    :: srcnm (nreg)*20, fname*210, dname*210, invc*5, unp*5, unx*5, flxk*25
    integer                      :: ireg, irec, jrec, ngrd(nreg), ngrdc(nreg)
    integer                      :: i, j, k, m, iy, im, ierr, jp, jx, ic
    integer                      :: p

    ! [flux vars]
    integer,parameter            :: ninv=1, nunp=3, nunx=3
    character(len=25),parameter  :: invcases(ninv) = ['inv1']
    character(len=25),parameter  :: flxcases(ninv) = ['fch4_gcp2025_prior_inv1']
    character(len=15),parameter  :: unxcases(nunx) = ['ctl', 'ux2', 'ux4']
    character(len=15),parameter  :: unpcases(nunp) = ['p30', 'p50', 'p99']

    real*4,allocatable,dimension (:,:,:,:) :: flxc, flux


    real(r8),parameter           :: basis_days = 365.D0
    integer                      :: regid(nlat,nlon), nrecl, cmdstat, exitstat
    character*512                :: io_file, msg, budget_file
    real(r8)                     :: expected_tdi(nreg), norm_seconds

  contains      

    !*****************************************************************
    !*** read_mask
    Subroutine read_mask()
      io_file = msk_file
      inquire(iolength=nrecl) binreg
      open(11,file=msk_file,status='old',form='unformatted',access='direct',recl=nrecl,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open mask')
      read(11,rec=1,iostat=ierr,iomsg=msg) ((binreg(j,i),i=1,nlon),j=1,nlat)
      call check_io(ierr,'Read mask')
      close(11,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close mask')
      if (.not.all(ieee_is_finite(binreg))) stop 'Non-finite mask'
      regid = 0
      do j = 1, nlat
      do i = 1, nlon
        do k = 1, nreg
          if (abs(dble(binreg(j,i))-k) <= 0.01D0) regid(j,i) = k
        end do
      end do
      end do
      print *, 'Mask assigned/unassigned cells: ',count(regid > 0),count(regid == 0)
      if (count(regid > 0) == 0) stop 'No recognized regions in mask'
    End subroutine read_mask
    
    
    !*****************************************************************
    !*** read land basis
    Subroutine read_landb()
      io_file = lnb_file
      inquire(iolength=nrecl) ((ocntr(1,j,i),i=1,nlon),j=1,nlat)
      open(50,file=lnb_file,status='old',access='direct',form='unformatted',recl=nrecl,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open land basis')
      ocntr = 0.0
      read(50,rec=1,iostat=ierr,iomsg=msg) ((ocntr(1,j,i),i=1,nlon),j=1,nlat)
      call check_io(ierr,'Read land basis')
      close(50,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close land basis')
      if (.not.all(ieee_is_finite(ocntr(1,:,:)))) stop 'Non-finite land basis'
      if (any(ocntr(1,:,:) < 0.0)) stop 'Negative land basis; check input'
    End subroutine read_landb

    !*****************************************************************
    !*** read_pred_srcs
    Subroutine read_pred_srcs()
      character*512 :: line
      real(r8) :: xtime,expected_time
      io_file = fname
      open(12,file=fname,status='old',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open predicted sources')
      tdi = 0.D0
      do m = 1, nreg
        read(12,'(A20)',iostat=ierr,iomsg=msg) srcnm(m)
        call check_io(ierr,'Read source name')
        if (len_trim(srcnm(m)) == 0) stop 'Empty predicted source name'
        if (m == 1) print *, 'Reg name: ',srcnm(m)
        read(12,'(A)',iostat=ierr,iomsg=msg) line
        call check_io(ierr,'Read source header')
        do i = 1, nyear
        do j = 1, nmon
          read(12,*,iostat=ierr,iomsg=msg) xtime,tdi(i,j,m)
          if (ierr /= 0) print *, 'Region/year/month: ',m,s_year+i-1,j
          call check_io(ierr,'Read predicted correction')
          if (.not.ieee_is_finite(xtime) .or. .not.ieee_is_finite(tdi(i,j,m))) stop 'Non-finite predicted source'
          expected_time = s_year+i-1+(dble(j)-0.5D0)/nmon
          if (abs(xtime-expected_time) > 0.05D0) stop 'Predicted source time mismatch; expected decimal year'
          if (m == 1) then
            time(i,j) = xtime
          else
            if (abs(time(i,j)-xtime) > 1.D-5) stop 'Different time axes across regions'
          end if
        end do
        end do
        ! repeat last years
        ! Disabled: retain the supplied corrections for every year, including 2024 and 2025.
      end do
      close(12,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close predicted sources')
      print *, 'Read ',nreg,' regions x ',nyear*nmon,' monthly corrections; no years replaced.'
          ! debug
          !if (m.eq.1) print '(2x, 4i6, 2f10.4)', m, i, i+s_year-1, j-1, time (i, j-1), tdi(i,j-1,m)
          !if (i > nyear-5) print '(2x, 4i6, 2f10.4)', m, i, i+s_year-1, j-1, time (i, j-1), tdi(i,j-1,m)
    End subroutine read_pred_srcs

    !*****************************************************************
    !*** read_flx
    Subroutine read_flx()
      io_file = fname
      inquire(iolength=nrecl) ((flux(1,1,j,i),i=1,nlon),j=1,nlat)
      open(50,file=fname, status='old', access='direct', form='unformatted', recl=nrecl, iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open prior flux')

      irec = 0
      do iy = 1, nyear
      do im = 1, nmon
        irec = irec + 1
        read (50, rec=irec,iostat=ierr,iomsg=msg) ((flux(iy, im, j, i), i = 1, nlon), j = 1, nlat)
        call check_io(ierr,'Read prior flux')
        if (.not.all(ieee_is_finite(flux(iy,im,:,:)))) stop 'Non-finite prior flux'
        if (any(flux(iy,im,:,:) == -999.0)) stop 'Missing prior flux (-999)'
      end do
      end do
      close(50,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close prior flux')
      print *,'Got init CH4 fluxes from : ', trim(fname)
      print *,'        # of recs : ', irec
    End subroutine read_flx

    !*****************************************************************
    !*** get_offset
    Subroutine get_offset()
      real(r8) :: lo,hi
      ofset = 0.D0
      do j = 1, nlat
        angle(j) = (-90.D0+(dble(j)-0.5D0)*180.D0/nlat)*pi/180.D0
        lo = (-90.D0+dble(j-1)*180.D0/nlat)*pi/180.D0
        hi = (-90.D0+dble(j)*180.D0/nlat)*pi/180.D0
        cellarea(j) = r0**2*(2.D0*pi/nlon)*(sin(hi)-sin(lo))
        ofset = ofset+sum(dble(ocntr(1,j,:)))*cellarea(j)
      end do
      norm_seconds = basis_days*86400.D0
      print *, 'Earth area [m2]: ',sum(cellarea)*nlon
      print *, 'Basic total offset [Tg/365-day year]: ',ofset*1.D-9*norm_seconds
    End subroutine get_offset    

    !*****************************************************************
    !*** convert from grid-point data to area weighted values
    Subroutine convert_flx()
      real(r8) :: value
      io_file = fname
      open(51,file=fname,status='replace',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open regional budget')
      flxc = 0.0
      totc = 0.D0
      ngrd = 0
      do j = 1, nlat
      do i = 1, nlon
        ireg = regid(j,i)
        if (ireg == 0) cycle
        ngrd(ireg) = ngrd(ireg)+1
        totc(ireg) = totc(ireg)+dble(ocntr(1,j,i))*cellarea(j)*norm_seconds*1.D-9
      end do
      end do
      do iy = 1, nyear
      do im = 1, nmon
        !applying corrections here
        do ireg = 1, nreg
          if (totc(ireg) <= 0.D0 .and. tdi(iy,im,ireg) /= 0.D0) then
            print *, 'Zero regional basis with nonzero correction: ',ireg,s_year+iy-1,im,tdi(iy,im,ireg)
            stop 1
          end if
        end do
        ngrdc = 0
        totlc = 0.D0
        write(51,*,iostat=ierr,iomsg=msg) s_year+iy-1,im
        call check_io(ierr,'Write budget date')
        do j = 1, nlat
        do i = 1, nlon
          ireg = regid(j,i)
          if (ireg == 0) cycle
          ngrdc(ireg) = ngrdc(ireg)+1
          if (totc(ireg) <= 0.D0) cycle
          value = dble(ocntr(1,j,i))*tdi(iy,im,ireg)/totc(ireg)
          flxc(iy,im,j,i) = real(value)
          if (.not.ieee_is_finite(flxc(iy,im,j,i))) stop 'Non-finite gridded correction'
          totlc(ireg) = totlc(ireg)+dble(flxc(iy,im,j,i))*cellarea(j)*norm_seconds*1.D-9
        end do
        end do
        do i = 1, nreg
          if (abs(totlc(i)-tdi(iy,im,i)) > 2.D-7*max(abs(tdi(iy,im,i)),1.D-30)) then
            print *, 'Region/target/actual: ',i,tdi(iy,im,i),totlc(i)
            stop 'Regional correction budget mismatch'
          end if
          write(51,'(2I9,5E16.7)',iostat=ierr,iomsg=msg) i,ngrdc(i),totc(i),totlc(i),tdi(iy,im,i)
          call check_io(ierr,'Write regional budget')
        end do
      enddo
      enddo
      close(51,iostat=ierr,iomsg=msg)!; pause
      call check_io(ierr,'Close regional budget')
              !totc(ireg) = totc(ireg) + 1.0 * cellarea(j) * 365.0           !old style
        !print *, ' totc ', totc
        !totcc=0.0
              ! debug
              !print'(i5, 10e12.2)', ireg, tdi(iy,im,ireg), totc(ireg), flxc(iy,im,j,i), flux(iy,im,j,i)
    End subroutine convert_flx

    !*****************************************************************
    !*** write_flx
    Subroutine write_flx()
      real*4 :: field(nlon,nlat)
      io_file = fname
      inquire(iolength=nrecl) field
      open(50,file=fname,access='direct',form='unformatted',status='replace',recl=nrecl,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open posterior flux')
      irec = 0
      do iy = 1, nyear
      do im = 1, nmon
        irec = irec+1
        ! Loss correction retains its original global scope, including unassigned cells.
        do j = 1, nlat
        do i = 1, nlon
          field(i,j) = real((dble(flux(iy,im,j,i))+dble(flxc(iy,im,j,i)))*lc_fact(iy))
        end do
        end do
        if (.not.all(ieee_is_finite(field))) stop 'Non-finite posterior flux'
        write(50,rec=irec,iostat=ierr,iomsg=msg) field
        call check_io(ierr,'Write posterior flux')
      end do
      enddo
      close(50,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close posterior flux')
        ! write (50, rec=irec) (( flxc(iy,im,j,i), i = 1, nlon), j =1, nlat)   !if priors not added
      !print *, '    Write TDI results for ', s_year+iy, irec 
    End subroutine write_flx

    !*****************************************************************
    !*** check_gtotal
    Subroutine check_gtotal()
      print '(7x, 10a)', '  year ', '   flxc   ',  '   prior   ', '   post   ', 'post*lc_f', '  LC_dif  ', ' invc  '
      do iy = 1, nyear
        sum1 = 0.0
        sum2 = 0.0
        do im = 1, nmon
          do j = 1, nlat
            do i = 1, nlon
              sum1 = sum1 + dble(flxc(iy,im,j,i)) * cellarea(j) * dble(month_days(s_year+iy-1,im)) * 86400.D0 * 1.D-9
              sum2 = sum2 + dble(flux(iy,im,j,i)) * cellarea(j) * dble(month_days(s_year+iy-1,im)) * 86400.D0 * 1.D-9
            end do
          end do
        end do
        pstf = sum1+sum2
        print '(7x, i6, 5f10.4, 2x,  a)', s_year+iy-1, sum1, sum2, pstf, pstf*lc_fact(iy), pstf*(lc_fact(iy) - 1.), trim(invc)
        io_file = budget_file
        write(71, '(i6, 5f10.4, 2x,  a)',iostat=ierr,iomsg=msg) s_year+iy-1, sum1, sum2, pstf, pstf*lc_fact(iy), &
          pstf*(lc_fact(iy) - 1.), trim(invc)
        call check_io(ierr,'Write global budget')
        gtotal(iy) = sum1+sum2
      end do
    End subroutine check_gtotal

    !*****************************************************************
    !*** read_lc_fact
    Subroutine read_lc_fact()
      integer :: ll,jj
      real(r8) :: xx
      character*512 :: fline
      logical :: seen(nyear)
      io_file = fname
      seen = .false.
      lc_fact = 0.D0
      open(61,file=fname,status='old',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open loss correction')
      read(61,'(A)',iostat=ierr,iomsg=msg) fline
      call check_io(ierr,'Read loss correction header')
      do
        read(61,'(A)',iostat=ierr,iomsg=msg) fline
        if (is_iostat_end(ierr)) exit
        call check_io(ierr,'Read loss correction row')
        if (len_trim(fline) == 0) cycle
        read(fline,*,iostat=ierr,iomsg=msg) jj,xx
        call check_io(ierr,'Parse loss correction row')
        if (jj < s_year .or. jj > f_year) cycle
        ll = jj-s_year+1
        if (seen(ll)) stop 'Duplicate loss correction year'
        if (.not.ieee_is_finite(xx) .or. xx <= 0.D0) stop 'Invalid loss correction factor'
        seen(ll) = .true.
        lc_fact(ll) = xx
      end do
      close(61,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close loss correction')
      do ll = 1, nyear
        if (.not.seen(ll)) then
          print *, 'Missing loss correction year: ',s_year+ll-1
          stop 1
        end if
      end do
        !print*, jj, xx, ierr, ll
          !print*,  'the total lines ', ll, s_year
            ! debug
            !print*, lc_fact(ll), ll
    End subroutine read_lc_fact

    
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
      print *, trim(action), ': ',trim(io_file)
      print *, trim(msg)
      stop 1
    End subroutine check_io

    Subroutine make_dir(path)
      character(len=*),intent(in) :: path
      if (index(path,"'") /= 0) stop 'Single quote in directory path'
      call execute_command_line("mkdir -p '"//trim(path)//"'",exitstat=exitstat,cmdstat=cmdstat,cmdmsg=msg)
      if (cmdstat /= 0) then
        print *, trim(msg)
        stop 'Cannot execute mkdir'
      end if
      if (exitstat /= 0) stop 'Cannot create output directory'
    End subroutine make_dir
  End module pf3_ch4_grd_mod
!
!*****************************************************************
!
  Program mainprog
    use pf3_ch4_grd_mod
    implicit none
    
    print*, 'Start program'
    call make_dir(out_dir)
    print *, 'TDI normalization: Tg CH4 per 365-day year; LC multiplier applies globally.'
    
    allocate (flux(1:nyear, 1:nmon, 1: nlat, 1: nlon) ); flux = 0.0 ! pkp / PKP
    allocate (flxc(1:nyear, 1:nmon, 1: nlat, 1: nlon) ); flxc = 0.0 ! 0.0 for invfiles; 9.9 for GrADS

    DO jp = 1, nunp
      unp = unpcases(jp)
      DO jx = 1, nunx
        unx = unxcases(jx)

        fname = out_dir // cgt_file // trim(unp) // '_'// trim(unx) //'.txt'
        budget_file = fname
        io_file = fname
        open(71,file=fname,status='replace',iostat=ierr,iomsg=msg)
        call check_io(ierr,'Open global budget')
        write(71, '(10a)',iostat=ierr,iomsg=msg) '  year ', '   flxc   ',  '   prior   ', '   post   ', &
          'post*lc_f', '  LC_dif  ', ' invc  '
        call check_io(ierr,'Write global budget header')

        DO ic = 1, ninv, 1
          invc = invcases(ic)
          flxk = flxcases(ic)
          
          print*,''
          print*,'=== BEGIN ================================================'
          print*,''
          print*, 'Run for : ', trim(invc), ' ', trim(unp), ' ', trim(unx)

          print*,''
          print*,'=== Read mask ============================================'
          print*, '      Read mask '
          call read_mask()

          print*,''
          print*,'=== Read predicted ======================================='
          fname = prd_dir // trim(unp) // '/'//icase//'_'// trim(unx) //'_' // trim(invc) //'/'//prd_file
          print'(A, A90)', '       ', fname
          call read_pred_srcs()
    
          print*,''
          print*,'=== Read landb ==========================================='
          call read_landb()
    
          print*,''
          print*,'=== Read prior ==========================================='
          fname = flx_dir // trim(flxk) // '.grd'
          print'(A, A90)', '       ', fname
          jrec = 0
          call read_flx
          
          print*,''
          print*,'=== Get offset ==========================================='
          call get_offset()
          
          print*,''
          print*,'=== Read loss correction ================================='
          fname = lsc_dir // trim(lsc_file) // trim(invc) // '.txt'
          print'(A, A90)', '       ', fname
          if (apply_lc) then
            call read_lc_fact()
          else
            lc_fact = 1.D0
            print *, 'Identity LC: posterior flux without loss correction'
          end if
    
          print*,''
          print*,'=== Convert flux ========================================='
          fname = out_dir//icase//'_'//trim(unp)//'_'//trim(unx)//'_'//trim(invc)//flux_suffix//'_fluxlnreg.dat'
          call convert_flx()
          
          print*,''
          print*,'=== Write output ========================================='
          fname = out_dir // trim(unp) // '/'//icase//'_'// trim(unx) //'_' // trim(invc) // flux_suffix // '.grd'
          p = index(fname, '/', back=.true.)
          if (p > 0) then
            dname = fname(:p-1)
          else
            dname = '.'
          end if
          call make_dir(dname)
          print'(A, A90)', '       ', fname
          call write_flx()
          
          print*,'===  Check G-total ======================================='
          call check_gtotal()
    
        ENDDO
        io_file = budget_file
        close(71,iostat=ierr,iomsg=msg)
        call check_io(ierr,'Close global budget')
      ENDDO
    ENDDO
    
    print*, 'End program'
  
  End program mainprog


