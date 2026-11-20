
  Module pf1_ch4_burloss_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none 

    ! [input files]
    character(len=*),parameter    :: grz_file = '/S/home01/G5070/y0715/MIROC4.0/data/a_data/t42/grz'
    ! CYC
    character(len=*),parameter    :: ch4_dir  = '../trout_MIROC/CH4_t42l67_INCA_H_apr_261104/'
    character(len=*),parameter    :: work_dir = '../results2025/INCA/'

    logical,parameter             :: kg_label_is_ppb = .true.
    logical,parameter             :: hourly_input = .true.

    ! [output files]
    ! dataset for loss correction
    character(len=*),parameter    :: burd_file = work_dir // 'losscorr/frt_burd_'


    ! [params]
    integer,parameter             :: s_year = 1999, f_year = 2025, nyear=f_year-s_year+1
    integer,parameter             :: istyr = s_year-1
    integer,parameter             :: nspcs = 3, ntrs=1
    integer,parameter             :: r8    = 8
    real(r8),parameter            :: pi    = 3.1415926535897932384626433832795_r8
    real(r8),parameter            :: avn   = 6.02214179D23  ! no. of atms per mole; Avogadro's Number
    real(r8),parameter            :: r0    = 6371.D3        ! radius [m]
    integer,parameter             :: imax  = 128, jmax= 64, kmax= 67
    integer,parameter             :: nmon  = 12
    integer,parameter             :: ounit = 30

    ! [grid vars]                
    real(r8)                      :: gar_t42(jmax)        ! grid size jmax [m**2]
    real(r8)                      :: dlat (jmax)          ! latitude weight
    real                          :: oval(imax,jmax,kmax)
    real(r8)                      :: garea(jmax), grz(imax,jmax)
    real(r8)                      :: bloss_y(3,nyear,nspcs), bloss_m(3,nyear*nmon,nspcs), bloss_gcp(4,nyear*nmon,nspcs)
    real(r8)                      :: ch4_ref(nyear*nmon)
    character*16                  :: head(64), ayear*19, tracer*6
    character*512                 :: fname, msg, output_name
    real(r8)                      :: summ, sum1, sum2, sum_sh, sum_nh, sum_tr, sum_st, tloss, burden
    integer                       :: i, j, k, iy, is, it, ierr, trec, tr

  ![data variables]
    integer,save                  :: mday(nmon)  = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    character(len=4),parameter    :: spcs(nspcs) = ['loh', 'lcl', 'lod']
    character(len=4),parameter    :: trs(ntrs)   = ['inv1']
    real(r8),save                 :: fvmr(nspcs) = (/1.D-9, 1.D-9, 1.D-9/)
    real(r8),parameter            :: molar_air = 28.97D0 ! dry air [g/mol]
    real(r8),save                 :: wtm(nspcs)  = (/16.D0, 16.D0, 16.D0/)
    integer,save                  :: lla_ref(3)  = (/120, 64, 2/) ! SPO coordinates in ACTM grid

    real(r8),dimension (imax,jmax,kmax,nmon)        :: ch4cc,ch4ls,ch4lt,zgt,airm !ch4tl,

    ! Original integration stops at kmax-1: no upper interface is supplied for kmax.
    integer,parameter             :: kint = kmax-1, ktrop = 30
    ! GGLA64 indexing assumed north to south; change only if the input axis differs.
    logical,parameter             :: north_first = .true.
    integer                       :: sh_first, sh_last, nh_first, nh_last, nrecl, ndays_year
    real(r8)                      :: dz(imax,jmax,kmax,nmon), total_weight, check_total
    real(r8),parameter            :: missing = -999.D0
    character*16                  :: ndm_units, ch4_units, loss_units, z_units, grz_units

  contains      

    !*****************************************************************
    !*** setgrid
    Subroutine setgrid_actm()

!     [internal]
      real(r8) :: ctheta (jmax)        !! Gauss Latitude (theta)
      real(r8) :: gw     (jmax)        !! Gauss weight
      real(r8) :: dlon                 ! diff. of longitude for GCM grid
      real(r8) :: alon (imax)          ! longitude for GCM
      real(r8) :: alat (0:jmax+1)      ! latitude
      real(r8) :: blat (jmax+1)        ! boundary latitude

      do i=1,imax
        alon(i) = 360.D0 * (i-1)/imax
      end do
      dlon = alon(2)-alon(1)

      call gauss (ctheta, gw)
      dlat = gw

      do j=1, jmax
        alat(j) = 180.D0 / pi *asin(ctheta(j))
      end do

      blat(1) = 90.D0
      blat(jmax+1) =-90.D0
      do j = 2, jmax
        blat(j) = (alat(j) + alat(j-1)) / 2.D0
      end do
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
    !*** for reading to gt3 formatted files
    Subroutine gt3read(klev,unit_in)
      integer,intent(in) :: klev        !! history levels
      integer :: ia, ja, ka, starts(3), ends(3), axis, status, idate(6)
      integer,optional,intent(in) :: unit_in
      integer :: iu
      real :: fill
      iu = ounit
      if (present(unit_in)) iu = unit_in
      read(iu,iostat=ierr,iomsg=msg) head
      call check_io(ierr,'Read GT3 header')
      if (trim(adjustl(head(38))) /= 'UR4') stop 'Expected GT3 UR4 data'
      do axis = 1, 3
        read(head(30+3*(axis-1)),*,iostat=status) starts(axis)
        if (status /= 0) stop 'Invalid GT3 axis start'
        read(head(31+3*(axis-1)),*,iostat=status) ends(axis)
        if (status /= 0) stop 'Invalid GT3 axis end'
      end do
      if (any(starts /= 1) .or. any(ends /= (/imax,jmax,klev/))) stop 'Unexpected GT3 dimensions'
      if (klev > 1 .and. .not.present(unit_in)) then
        read(head(27),'(I4,I2,I2,1X,I2,I2,I2)',iostat=status) idate
        if (status /= 0) stop 'Invalid GT3 date'
        if (idate(1) /= istyr+iy .or. idate(2) /= it) then
          print '(a,2(1x,i0),a,6(1x,i0))', 'Expected year/month ',istyr+iy,it,' found ',idate
          stop 'GT3 time mismatch'
        end if
      end if
      oval = 0.0
      read(iu,iostat=ierr,iomsg=msg) (((oval(ia,ja,ka),ia=1,imax),ja=1,jmax),ka=1,klev)
      call check_io(ierr,'Read GT3 field')
      !if (.not.all(ieee_is_finite(oval(:,:,1:klev)))) stop 'Non-finite GT3 data'
      if (.not.all(ieee_is_finite(oval(:,:,1:min(klev,kint))))) then
        print '(2A)', 'File: ',trim(fname)
        print '(2A)', 'GT3 date: ',trim(head(27))
        print '(A,I0)', 'Non-finite cells: ',count(.not.ieee_is_finite(oval(:,:,1:min(klev,kint))))
        do ka = 1, min(klev,kint)
        do ja = 1, jmax
        do ia = 1, imax
          if (.not.ieee_is_finite(oval(ia,ja,ka))) then
            print '(A,3I6)', 'First invalid i,j,k: ',ia,ja,ka
            stop 'Non-finite GT3 data'
          end if
        end do
        end do
        end do
      end if
      
      read(head(39),*,iostat=status) fill
      if (status == 0) then
        if (any(oval(:,:,1:min(klev,kint)) == fill)) stop 'Missing GT3 data'
      end if
      return
    End subroutine gt3read
    
    !*****************************************************************
    !*** horizontal grid area in m^2
    Subroutine get_garea()
      do j = 1, jmax
        garea(j) = 2.d0 * pi * r0**2 * 1.d0 / dble(imax) * dlat(j)
      end do
      gar_t42 = garea
      print '(a,f0.3,a)', '  Earth area ', sum(garea)*imax, ' m2'
      print '(a,f0.3,a)', '  Earth area ', sum(garea)*imax*1.D-6, ' km2'
      print '(*(a))', ''
    End subroutine get_garea
    
    !*****************************************************************
    !*** read ndm
    Subroutine read_ndm()
      open(ounit,file=fname,form='unformatted',convert='big_endian',status='old',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open input')
      do it = 1, nmon
        oval = 0.0
        call gt3read(kmax)
        if (any(oval <= 0.0)) stop 'Non-positive air number density'
        airm (:, :, :, it) = dble(oval (:, :, :))
        if (it == 1) then
          ndm_units = head(16)
          call check_units(ndm_units,'density')
          print '(*(a))', 'NDM units (expected molecules/cm3): ', trim(ndm_units)
        else
          if (head(16) /= ndm_units) stop 'NDM units changed within year'
        end if
      end do
      close(ounit,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close input')
    End subroutine read_ndm
    
    !*****************************************************************
    !*** read loss
    Subroutine read_loss()
      open(ounit,file=fname,form='unformatted',convert='big_endian',status='old',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open input')
      do it = 1, nmon
        oval = 0.0
        call gt3read(kmax)
        if (any(oval < 0.0)) stop 'Negative loss: check sign convention or missing values'
        ch4lt (:, :, :, it) = dble(oval (:, :, :))
        if (it == 1) then
          loss_units = head(16)
          call check_units(loss_units,'loss')
          print '(*(a))', trim(spcs(is)), ' units (expected molecules/cm3/s): ', trim(loss_units)
          if (trim(adjustl(loss_units)) == 's-1' .or. trim(adjustl(loss_units)) == '1/s') then
            stop 'First-order loss rate requires multiplication by CH4 number density; check input units'
          end if
        else
          if (head(16) /= loss_units) stop 'Loss units changed within year'
        end if
      end do
      close(ounit,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close input')
    End subroutine read_loss

    !*****************************************************************
    !*** read ch4 gt3
    Subroutine read_ch4_gt3()
      real(r8) :: mix_scale
      open(ounit,file=fname,form='unformatted',convert='big_endian',status='old',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open input')
      do it = 1, nmon
        oval = 0.0
        call gt3read(kmax)
        call check_io(ierr,'Read CH4 field')
        if (.not.all(ieee_is_finite(oval)) .or. any(oval < 0.0)) stop 'Invalid CH4 concentration'

        if (it == 1) then
          ch4_units = head(16)
          call check_units(ch4_units,'mixing',mix_scale)
          print '(*(a))', 'CH4 input units: ', trim(ch4_units)
        else
          if (head(16) /= ch4_units) stop 'CH4 units changed within year'
        end if
        ch4cc (:, :, :, it) = dble(oval (:, :, :)) * airm (:, :, :, it) * mix_scale
        trec = (iy -1)*12 + it
        ch4_ref(trec) = dble(oval(lla_ref(1),lla_ref(2),lla_ref(3)))*(mix_scale*1.D9)
      end do
      close(ounit,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close input')
    End subroutine read_ch4_gt3

    !*****************************************************************
    !*** read ch4 grd
    Subroutine read_ch4_grd()
      inquire(iolength=nrecl) oval
      open (ounit, file=fname, form='unformatted', access='direct', &
                               status='old', recl=nrecl, iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open input')
      do it = 1, nmon
        oval = 0.0
        read (ounit, rec=it,iostat=ierr,iomsg=msg) ((( oval(i,j,k), i = 1, imax), j = 1, jmax), k = 1, kmax)
        call check_io(ierr,'Read CH4 field')
        if (.not.all(ieee_is_finite(oval)) .or. any(oval < 0.0)) stop 'Invalid CH4 concentration'
        ch4cc (:, :, :, it) = dble(oval (:, :, :)) * airm (:, :, :, it) * fvmr(1)
        
        trec = (iy -1)*12 + it
        ch4_ref(trec) = oval (lla_ref(1), lla_ref(2), lla_ref(3))
      end do
      close(ounit,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close input')
    End subroutine read_ch4_grd

    !*****************************************************************
    !*** read Z
    Subroutine read_z()
      open(ounit,file=fname,form='unformatted',convert='big_endian',status='old',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open input')
      do it = 1, nmon
        oval = 0.0
        call gt3read(kmax)
        zgt (:, :, :, it) = dble(oval (:, :, :))
        if (it == 1) then
          z_units = head(16)
          call check_units(z_units,'height')
          print '(*(a))', 'Z units (expected m): ', trim(z_units)
        else
          if (head(16) /= z_units) stop 'Z units changed within year'
        end if
      end do
      close(ounit,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close input')
     End subroutine read_z

    !*****************************************************************
    !*** read GRZ [surface height, m]
    Subroutine read_grz()
      open(ounit,file=fname,form='unformatted',convert='big_endian',status='old',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open input')
      call gt3read(1)
      grz_units = head(16)
      call check_units(grz_units,'height')
      grz (:, :) = dble(oval (:, :, 1))
      print '(*(a))', 'GRZ units (expected m): ', trim(grz_units)
      close(ounit,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close input')    
    End subroutine read_grz

    !*****************************************************************
    integer Function month_days(year,month) result(days)
      integer,intent(in) :: year,month
      days = mday(month)
      if (month == 2) then
        if (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) days = 29
      end if
    End function month_days

    Subroutine get_dz()
      integer :: im, lev
      dz = 0.D0
      do im = 1, nmon
        if (any(zgt(:,:,1,im) < grz)) stop 'Lowest level below surface; check height units'
        do lev = 2, kmax
          if (any(zgt(:,:,lev,im) <= zgt(:,:,lev-1,im))) stop 'Non-increasing model height'
        end do
        dz(:,:,1,im) = 0.5D0*(zgt(:,:,1,im)+zgt(:,:,2,im))-grz
        do lev = 2, kint
          dz(:,:,lev,im) = 0.5D0*(zgt(:,:,lev+1,im)-zgt(:,:,lev-1,im))
        end do
        if (any(dz(:,:,1:kint,im) <= 0.D0)) stop 'Non-positive layer thickness'
      end do
    End subroutine get_dz

    Subroutine check_units(unit_text,quantity,mix_scale)
      character(len=*),intent(in) :: unit_text,quantity
      real(r8),optional,intent(out) :: mix_scale ! input CH4 to mol/mol
      character*64 :: unit_norm
      integer :: ix, ic, nc
      logical :: accepted
      unit_norm = ''
      nc = 0
      do ix = 1, len_trim(unit_text)
        if (unit_text(ix:ix) == ' ' .or. unit_text(ix:ix) == '^' .or. unit_text(ix:ix) == '*') cycle
        nc = nc+1
        ic = iachar(unit_text(ix:ix))
        if (ic >= iachar('A') .and. ic <= iachar('Z')) ic = ic+32
        unit_norm(nc:nc) = achar(ic)
      end do
      ! ACTM GT3 writes cubic centimetres as cm|3".
      ix = index(unit_norm,'cm|3"')
      if (ix > 0) unit_norm = unit_norm(:ix-1)//'cm3'//unit_norm(ix+5:)
      accepted = .false.
      select case(quantity)
      case('height')
        accepted = unit_norm == 'm' .or. unit_norm == 'meter' .or. unit_norm == 'metre'
      case('mixing')
        select case(trim(unit_norm))
        case('ppb','ppbv','nmol/mol')
          accepted = .true.
          if (present(mix_scale)) mix_scale = 1.D-9
        case('kg/kg')
          accepted = .true.
          if (present(mix_scale)) mix_scale = merge(1.D-9,molar_air/wtm(1),kg_label_is_ppb)
        end select
      case('density')
        select case(trim(unit_norm))
        case('cm-3','1/cm3','molec/cm3','molecules/cm3','molecule/cm3','moleccm-3','moleculescm-3','moleculecm-3')
          accepted = .true.
        end select
      case('loss')
        select case(trim(unit_norm))
        case('cm-3s-1','1/cm3/s','molec/cm3/s','molecules/cm3/s','molecule/cm3/s')
          accepted = .true.
        case('moleccm-3s-1','moleculescm-3s-1','moleculecm-3s-1','molec/cm3s','molecules/cm3s')
          accepted = .true.
        end select
      end select
      if (.not.accepted) then
        print '(*(a))', 'File: ',trim(fname)
        print '(*(a))', 'Unconfirmed units for ',quantity,': ',trim(unit_text)
        stop 'Check units before applying the conversion factors'
      end if
    End subroutine check_units

    Subroutine check_io(status,action)
      integer,intent(in) :: status
      character(len=*),intent(in) :: action
      if (status == 0) return
      print '(*(a))', trim(action), ': ', trim(fname)
      print '(*(a))', trim(msg)
      stop 1
    End subroutine check_io

    Subroutine process_hourly_year()
      integer,parameter :: nf = 3+nspcs
      integer :: units(nf), f, day, hour, nhours, year_hours, stamp(6), ref_stamp(6)
      character*512 :: files(nf), prefix
      character*16 :: labels(nf), extra(64)
      character(len=7),parameter :: quantity(3) = ['height ','density','mixing ']
      real(r8) :: volume(imax,jmax,kint), density(imax,jmax,kint), heights(imax,jmax,kmax)
      real(r8) :: hemis(2), annual_b(2), annual_l(2,nspcs), vertical(2,nspcs), parts(2,2)
      real(r8) :: rate, row, factor, ref_sum, mix_scale, raw_ref
      integer :: lev, lat, h, v, species

      prefix = trim(ch4_dir)//trim(ayear)
      files(1) = trim(prefix)//'z.gt3'
      files(2) = trim(prefix)//'ndm.gt3'
      files(3) = trim(prefix)//'ch4_'//trim(tracer)//'.gt3'
      do species = 1, nspcs
        files(3+species) = trim(prefix)//trim(spcs(species))//'_'//trim(tracer)//'.gt3'
      end do
      do f = 1, nf
        fname = files(f)
        open(newunit=units(f),file=fname,form='unformatted',convert='big_endian',status='old',action='read',iostat=ierr,iomsg=msg)
        call check_io(ierr,'Open hourly GT3')
      end do
      annual_b = 0.D0
      annual_l = 0.D0
      vertical = 0.D0
      year_hours = 0
      do it = 1, nmon
        trec = (iy-1)*nmon+it
        nhours = 24*month_days(istyr+iy,it)
        hemis = 0.D0
        ref_sum = 0.D0
        bloss_m(2,trec,:) = 0.D0
        do day = 1, month_days(istyr+iy,it)
        do hour = 0, 23
          ! The same sequential field number in every annual file represents the same model hour.
          do f = 1, nf
            fname = files(f)
            call gt3read(kmax,units(f))
            read(head(27),'(i4,2i2,1x,3i2)',iostat=ierr,iomsg=msg) stamp
            call check_io(ierr,'Read hourly GT3 date')
            if (f == 1) then
              ref_stamp = stamp
              if (any(stamp(1:3) /= [istyr+iy,it,day])) then
                print '(a,3(1x,i0),a,i0,a,a)', 'Expected date',istyr+iy,it,day,' record ',hour+1,'; found ',trim(head(27))
                stop 1
              end if
            else if (any(stamp /= ref_stamp)) then
              print '(a,6(1x,i0))', 'z.gt3 time:',ref_stamp
              print '(2a,6(1x,i0))', trim(fname),' time:',stamp
              stop 1
            end if
            if (it == 1 .and. day == 1 .and. hour == 0) then
              labels(f) = head(16)
              if (f == 3) then
                call check_units(labels(f),'mixing',mix_scale)
              else if (f <= 2) then
                call check_units(labels(f),trim(quantity(f)))
              else
                call check_units(labels(f),'loss')
              end if
              print '(1x,3a)', trim(fname),'  units: ',trim(adjustl(labels(f)))
            else if (labels(f) /= head(16)) then
              print '(3a,4(1x,i0))', 'Units changed in ',trim(fname),' at model hour ',istyr+iy,it,day,hour
              stop 1
            end if

            select case(f)
            case(1)
              heights = dble(oval)
              if (any(heights(:,:,1) < grz)) error stop 'Lowest level below surface'
              do lev = 2, kmax
                if (any(heights(:,:,lev) <= heights(:,:,lev-1))) error stop 'Non-increasing model height'
              end do
              volume(:,:,1) = 0.5D0*(heights(:,:,1)+heights(:,:,2))-grz
              do lev = 2, kint
                volume(:,:,lev) = 0.5D0*(heights(:,:,lev+1)-heights(:,:,lev-1))
              end do
              do lev = 1, kint
              do lat = 1, jmax
                volume(:,lat,lev) = volume(:,lat,lev)*garea(lat)
              end do
              end do
            case(2)
              if (any(oval <= 0.0)) error stop 'Non-positive air number density'
              density = dble(oval(:,:,1:kint))
            case(3)
              if (any(oval < 0.0)) error stop 'Negative CH4 concentration'
              factor = mix_scale*1.D6/avn*wtm(1)*1.D-12
              if (day == 1 .and. hour == 0 .and. it == 1) then
                raw_ref = dble(oval(lla_ref(1),lla_ref(2),lla_ref(3)))
                print '(a,3es16.6)', '    First CH4: raw min/max/SPO ',minval(oval),maxval(oval),raw_ref
                print '(a,2es16.6)', '    CH4 scale to mol/mol; SPO [ppb] ',mix_scale,raw_ref*mix_scale*1.D9
                print '(a,2es16.6)', '    First NDM min/max [molecules/cm3] ',minval(density),maxval(density)
                print '(a,es16.6)', '    Integrated dry air mass [kg] ',sum(density*volume)*(1.D6/avn)*(molar_air*1.D-3)
              end if
              ref_sum = ref_sum+dble(oval(lla_ref(1),lla_ref(2),lla_ref(3)))*(mix_scale*1.D9)
              do lev = 1, kint
              do lat = 1, jmax
                h = merge(1,2,lat >= sh_first .and. lat <= sh_last)
                row = sum(dble(oval(:,lat,lev))*density(:,lat,lev)*volume(:,lat,lev))
                hemis(h) = hemis(h)+row*factor
              end do
              end do
              if (it == 1 .and. day == 1 .and. hour == 0) print '(a,es16.6)', '    First CH4 burden [Tg] ',sum(hemis)
            case default
              !if (any(oval < 0.0)) error stop 'Negative chemical loss'
              if (any(oval(:,:,1:kint) < 0.0)) then
                print '(2a)', 'File: ',trim(fname)
                print '(2a)', 'GT3 date: ',trim(head(27))
                print '(a,2es16.7)', 'Loss min/max, integrated levels: ', &
                  minval(oval(:,:,1:kint)),maxval(oval(:,:,1:kint))
                print '(a,3i6)', 'Minimum at i,j,k: ',minloc(oval(:,:,1:kint))
                print '(a,i0)', 'Negative cells: ',count(oval(:,:,1:kint) < 0.0)
                !error stop 'Negative chemical loss'
                if (day > 30) error stop 'Negative chemical loss'
              end if
              
              species = f-3
              factor = 3600.D0*1.D6/avn*wtm(species)*1.D-12
              parts = 0.D0
              do lev = 1, kint
                v = merge(1,2,lev <= ktrop)
                do lat = 1, jmax
                  h = merge(1,2,lat >= sh_first .and. lat <= sh_last)
                  rate = sum(dble(oval(:,lat,lev))*volume(:,lat,lev))
                  parts(h,v) = parts(h,v)+rate*factor
                end do
              end do
              annual_l(:,species) = annual_l(:,species)+sum(parts,dim=2)
              vertical(:,species) = vertical(:,species)+sum(parts,dim=1)
              bloss_m(2,trec,species) = bloss_m(2,trec,species)+sum(parts)
            end select
          end do
        end do
        end do
        annual_b = annual_b+hemis
        year_hours = year_hours+nhours
        bloss_m(1,trec,:) = sum(hemis)/nhours
        ch4_ref(trec) = ref_sum/nhours
        print '(a,i4,1x,i2.2,a,i4)', '    Completed ',istyr+iy,it,'  hours: ',nhours
      end do
      do f = 1, nf
        fname = files(f)
        read(units(f),iostat=ierr,iomsg=msg) extra
        if (.not.is_iostat_end(ierr)) then
          call check_io(ierr,'Check end of hourly file')
          print '(*(a))', 'Extra record after December 31 in ',trim(fname)
          stop 1
        end if
        close(units(f),iostat=ierr,iomsg=msg)
        call check_io(ierr,'Close hourly GT3')
      end do
      annual_b = annual_b/year_hours
      burden = sum(annual_b)
      print '(a,3es16.6)', '  NH/SH/Total Annual burden ',annual_b(2),annual_b(1),burden
      do species = 1, nspcs
        tloss = sum(annual_l(:,species))
        check_total = sum(vertical(:,species))
        if (abs(check_total-tloss) > 1.D-10*max(tloss,1.D-30)) error stop 'Vertical loss budget mismatch'
        print '(2a,3f12.4)', '    NH/SH Annual loss ',spcs(species),annual_l(2,species),annual_l(1,species),tloss
        print '(2a,3f12.4)', '    Level 1:30 / 31:66 / total loss ',spcs(species),vertical(:,species),check_total
        bloss_y(1,iy,species) = burden
        bloss_y(2,iy,species) = tloss
        if (tloss > 0.D0) bloss_y(3,iy,species) = burden/tloss
        do it = 1, nmon
          trec = (iy-1)*nmon+it
          if (bloss_m(2,trec,species) <= 0.D0) cycle
          total_weight = dble(year_hours)/(24*month_days(istyr+iy,it))
          bloss_m(3,trec,species) = bloss_m(1,trec,species)/(bloss_m(2,trec,species)*total_weight)
        end do
      end do
    End subroutine process_hourly_year

  End module pf1_ch4_burloss_mod

  Program mainprog
    use pf1_ch4_burloss_mod
    implicit none
    
    print '(*(a))',''
    print '(*(a))', '=== Start program ======================================='
    print '(*(a))',''
    print '(*(a))', '  Input model   dir  ', ch4_dir
    print '(*(a))', '  Output working dir  ', work_dir
    print '(*(a))',''
          
    iy = 0
    it = 0
    sh_first = 1
    sh_last = jmax/2
    nh_first = jmax/2+1
    nh_last = jmax
    if (north_first) then
      nh_first = 1
      nh_last = jmax/2
      sh_first = jmax/2+1
      sh_last = jmax
    end if
    print '(a,i0,a)', 'Integrated levels: 1 to ',kint,'; top level excluded (upper interface unavailable).'
    print '(*(a))', 'Level 30 split is fixed; it is not a diagnosed tropopause.'
    print '(*(a))', 'Units: CH4 ppb or kg/kg dry air; NDM molec/cm3; loss molec/cm3/s; height m.'
    print '(*(a))', 'Annual burden is day-weighted; output ch4_ref remains the December reference.'
    if (hourly_input) print '(a)', 'Hourly GT3: 24 matched records per calendar day; 3600 s per field.'
    call setgrid_actm
    call get_garea

    fname = trim(grz_file)
    call read_grz

    DO tr = 1, ntrs
      tracer = trs(tr)
      print '(*(a))', '  Tracer = ', tracer
      
      bloss_y = missing
      bloss_m = missing
      bloss_gcp = missing
      ch4_ref = missing
      output_name = burd_file//trim(tracer)//'.txt'
      fname = output_name
      open(60,file=output_name,status='replace',iostat=ierr,iomsg=msg)
      call check_io(ierr,'Open burden output')
      write (60, '(20a10)',iostat=ierr,iomsg=msg) 'year    ', 'ch4_ref ', 'bch4    ', (spcs(is), is = 1, nspcs)
      call check_io(ierr,'Write burden header')
      !open (63, file = blss_file//trim(tracer)//'.txt', status = 'unknown')
      !write (63, '(20a10)') 'year ', 'month ', 'ch4_ref ', 'bch4  ', (spcs(is), is = 1, nspcs)
      !open (62, file = ltim_file, status = 'unknown')
      !write (62, '(20a10)') 'year ', (spcs(is), is = 1, nspcs)
      
      !
      ! ------ loop years
      !
      DO iy = 1, nyear
        print '(*(a))', ''
        print '(a, i6)', '    Year = ', istyr+iy
        ! get dirs
        write (ayear, '(i4.4,a1)') istyr+iy, '/'

        if (hourly_input) then
          call process_hourly_year()
        else
        ch4cc = 0.0
        zgt = 0.0
        sum_sh = 0.D0
        sum_nh = 0.D0
        is = 1
        !
        ! ------ reading files
        !
        ! ndm
        fname = trim(ch4_dir) // trim(ayear) // 'ndm.gt3'
        call read_ndm()

        ! ch4
        fname = trim(ch4_dir) // trim(ayear) // 'ch4_'//trim(tracer)//'.gt3'
        call read_ch4_gt3()

        ! z
        fname = trim(ch4_dir) // trim(ayear) // 'z.gt3'
        call read_z
        call get_dz
        ndays_year = 0

        !
        ! ------ estimating atmospheric burden of CH4
        !
        do it = 1, nmon
            ! === For Southern Hemisphere; surface layer first ===
            summ = 0.D0
            do j = sh_first, sh_last; 
            do i = 1, imax 
              if (ch4cc(i,j,1,it) .gt. 0.0) then
                summ = summ + ch4cc(i,j,1,it)*(1.D6/avn)*garea(j) *   &
                       dz(i,j,1,it) ! layer 1
              end if

              !for all other layers
              do k = 2, kmax-1 
                if (ch4cc(i,j,k,it) .gt. 0.0) then   
                  summ = summ + ch4cc(i,j,k,it)*(1.D6/avn)*garea(j) * &   ! 1e6 for cc==>>m3
                         dz(i,j,k,it)
                end if 
              end do
            end do
            end do
            sum1 = summ*1.0d-12*wtm(is)   ! in Tg-CH4

            ! === For Northern Hemisphere; surface layer first ===
            summ = 0.D0
            do j = nh_first, nh_last
            do i = 1, imax 
              if (ch4cc(i,j,1,it) .gt. 0.0) then
                 summ = summ + ch4cc(i,j,1,it)*(1.D6/avn)*garea(j) *   &
                        dz(i,j,1,it) ! layer 1
              end if

              !for all other layers
              do k = 2, kmax-1 
                if (ch4cc(i,j,k,it) .gt. 0.0) then   
                   summ = summ + ch4cc(i,j,k,it)*(1.D6/avn)*garea(j) *   &
                     dz(i,j,k,it)
                end if 
              end do
            end do
            end do
            
            sum2 = summ*1.0d-12*wtm(is)   ! in Tg-CH4

            sum_sh = sum_sh + sum1*dble(month_days(istyr+iy,it))
            sum_nh = sum_nh + sum2*dble(month_days(istyr+iy,it))
            ndays_year = ndays_year+month_days(istyr+iy,it)
            trec = (iy -1)*12 + it
            bloss_m (1, trec, :) = sum1+sum2
            !print'(a30,i5,3f12.4)','Monthly burden ', it, sum_nh, sum_sh, sum_sh+sum_nh
        end do
        burden = (sum_sh+sum_nh)/dble(ndays_year)
        print '(a,3es16.6)','  NH/SH/Total Annual burden ',sum_nh/dble(ndays_year), sum_sh/dble(ndays_year), burden

        !
        ! ------ loop chemical speaces
        !
        DO is = 1, nspcs
          ch4ls = 0.0
          ch4lt = 0.0

          ! spcs
          fname = trim(ch4_dir) // trim(ayear) //  trim(spcs(is))//'_'//trim(tracer)//'.gt3'
          call read_loss()

          ! ------ estimating atmospheric loss of CH4 in northern/southern hemispheres
          sum_sh = 0.D0
          sum_nh = 0.D0
          ch4ls = ch4lt * 86400.D0  ! molecules/cm3/s ==>> molecules/cm3/day
          ! print *, 'CH4ls ', ch4ls (1,32,3,:); print *, 'CH4lt ', ch4lt (1,32,3,:); print *, 'CH4cc ', ch4cc (1,32,3,:)

          do it = 1, nmon
            ! For Southern Hemisphere; surface leyer first 
            summ = 0.D0
            do j = sh_first, sh_last
            do i = 1, imax
              if (ch4ls(i,j,1,it) .gt. 0.0) then
                summ = summ + ch4ls(i,j,1,it)*(1.D6/avn)*garea(j) *   &
                      dz(i,j,1,it) ! layer 1
              end if

              !for all other layers
              do k = 2, kmax-1 
                if (ch4ls(i,j,k,it) .gt. 0.0) then
                  summ = summ + ch4ls(i,j,k,it)*(1.D6/avn)*garea(j) *  &   ! 1e6 for cc==>>m3
                     dz(i,j,k,it)
                end if
              end do
            end do
            end do
            sum1 = summ*1.0d-12*wtm(is)*dble(month_days(istyr+iy,it))   ! in Tg-CH4/month

            ! For Northern Hemisphere; surface leyer first 
            summ = 0.D0
            do j = nh_first, nh_last
            do i = 1, imax
              if (ch4ls(i,j,1,it) .gt. 0.0) then
                summ = summ + ch4ls(i,j,1,it)*(1.D6/avn)*garea(j) *   &
                   dz(i,j,1,it) ! layer 1
              end if 
              
              !for all other layers
              do k = 2, kmax-1
                if (ch4ls(i,j,k,it) .gt. 0.0) then
                   summ = summ + ch4ls(i,j,k,it)*(1.D6/avn)*garea(j) *  &
                     dz(i,j,k,it)
                end if
              end do
            end do
            end do
            sum2 = summ*1.0d-12*wtm(is)*dble(month_days(istyr+iy,it))   ! in Tg-CH4/month

            sum_sh = sum_sh + sum1
            sum_nh = sum_nh + sum2 
            trec = (iy -1)*12 + it
            bloss_m (2, trec, is) = sum1+sum2
          end do
          print '(2a,3f12.4)','    NH/SH Annual loss   ', spcs(is), sum_nh, sum_sh, sum_sh+sum_nh   ! for months
          tloss = sum_sh+sum_nh

          ! print '(a,f12.4, a)', '      Lifetime, burden/loss ', burden/tloss !, ' year' !; print *, ' '
          bloss_y (1, iy, is) = burden
          bloss_y (2, iy, is) = tloss
          if (tloss > 0.D0) bloss_y (3, iy, is) = burden/tloss
          do it = 1, nmon
            trec = (iy-1)*nmon+it
            if (bloss_m(2,trec,is) > 0.D0) then
              total_weight = dble(ndays_year)/dble(month_days(istyr+iy,it))
              bloss_m(3,trec,is) = bloss_m(1,trec,is)/(bloss_m(2,trec,is)*total_weight)
            end if
          end do

          ! ------ estimating atmospheric loss of CH4 for troposphere/stratosphere
          sum_tr = 0.D0
          sum_st = 0.D0
          do it = 1, nmon 

            ! For troposphere
            summ = 0.D0
            do j = 1, jmax; do i = 1, imax; do k = 1, ktrop
            if (ch4ls(i,j,k,it) .gt. 0.0) then
              summ = summ + ch4ls(i,j,k,it)*(1.D6/avn)*garea(j) *  &   ! 1e6 for cc==>>m3
                dz(i,j,k,it)
            end if
            end do; end do; end do
            sum1 = summ*1.0d-12*wtm(is)*dble(month_days(istyr+iy,it))   !  in Tg-CH4/month

            ! For stratosphere
            summ = 0.D0
            do j = 1, jmax; do i = 1, imax; do k = ktrop+1, kint
            if (ch4ls(i,j,k,it) .gt. 0.0) then
              summ = summ + ch4ls(i,j,k,it)*(1.D6/avn)*garea(j) *  &
                dz(i,j,k,it)
            end if
            end do; end do; end do
            sum2 = summ*1.0d-12*wtm(is)*dble(month_days(istyr+iy,it))   !  in Tg-CH4/month

            sum_tr = sum_tr + sum1
            sum_st = sum_st + sum2
            trec = (iy -1)*12 + it

          end do
          check_total = sum_tr+sum_st
          if (abs(check_total-tloss) > 1.D-10*max(abs(tloss),1.D-30)) stop 'Vertical loss budget mismatch'
          print '(2A,3F12.4)', '    Level 1:30 / 31:66 / total loss ',spcs(is),sum_tr,sum_st,check_total
          !print '(2a,3f12.4)', '      Tr/St Annual loss ', spcs(is), sum_tr, sum_st, sum_tr+sum_st
          !print *, ''

        END DO   ! spes loop
        end if
        fname = output_name
        if (max(abs(ch4_ref(iy*12)),abs(bloss_y(1,iy,1)),maxval(abs(bloss_y(2,iy,:)))) >= 9999999.D0) then
          write(60,'(I5,5(1X,ES16.6E3))',iostat=ierr,iomsg=msg) istyr+iy,ch4_ref(iy*12),bloss_y(1,iy,1),bloss_y(2,iy,:)
        else
          write(60,'(I5,20F12.4)',iostat=ierr,iomsg=msg) istyr+iy,ch4_ref(iy*12),bloss_y(1,iy,1),bloss_y(2,iy,:)
        end if
        call check_io(ierr,'Write annual burden/loss')
        !write (62, '(i5, 20f12.4)') istyr+iy, (bloss_y (3, iy, k), k = 1, nspcs)
        !do it = 1, nmon 
        !  trec = (iy -1)*12 + it
        !  write (63, '(2i5, 20f12.4)') istyr+iy, it, ch4_ref(trec), &
        !                               bloss_m (1, trec, 1), (bloss_m (2, trec, k), k = 1, nspcs)
        !enddo
      END DO   ! year loop
      close(60,iostat=ierr,iomsg=msg)
      call check_io(ierr,'Close burden output')
      !close(63)
    END DO   ! tracer loop
    

    print '(*(a))', 'Program finish OK'
  End
