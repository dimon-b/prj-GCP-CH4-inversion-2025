! =============================================================
! === INCA OH NetCDF -> GTOOL3 (T42L67)
! Original converter: i1_INCA_OH_nc2gt3_t42.f90
! Revised for GCP 2025: uses command-line paths, reads one month at a time.
! The Python launcher validates calendar, dimensions, units and GT3 output.
!
! Build: ifort -O2 -traceback df0_inca_oh_nc2gt3_t42.f90 -o df0_inca_oh_nc2gt3_t42 $(nf-config --fflags) $(nf-config --flibs)
! Run:   ./df0_inca_oh_nc2gt3_t42 INPUT.nc OUTPUT.gt3 START_YEAR
! =============================================================

module df0_inca_oh_nc2gt3_t42_mod
    use netcdf
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: imax=128, jmax=64, kmax=67
    integer, parameter :: nmon=12, hlen=64
    integer, parameter :: mday(12) = (/16,15,16,16,16,16,16,16,16,16,16,16/)
    integer, parameter :: nhr(12)  = (/12, 0,12, 0,12, 0,12,12, 0,12, 0,12/)
    real(real64), parameter :: lonstep=360.0_real64/real(imax,real64)

contains

    ! =============================================================
    ! === NetCDF call check
    subroutine check(status, context)
        integer, intent(in) :: status
        character(*), intent(in) :: context
        if (status /= nf90_noerr) then
            print '(3A)', 'NETCDF ERROR [',trim(context),']: '//trim(nf90_strerror(status))
            stop 1
        endif
    end subroutine check

    ! =============================================================
    ! === One GTOOL3 record (original ACTM header convention)
    subroutine write_gt3_record(unit, vals, year, month, itime)
        integer, intent(in) :: unit, year, month, itime
        real(real32), intent(in) :: vals(imax,jmax,kmax)
        character(16) :: head(hlen)
        integer :: jdate(6), i, ios
        character(512) :: msg

        jdate = (/year,month,mday(month),nhr(month),0,0/)
        head = ''
        write(head(1),'(I16)') 9010
        head(2) = 'GCP-INCA-OH'
        head(3) = 'NDOH'
        write(head(12),'(I16)') 0
        write(head(13),'(I16)') 1
        head(14) = 'GCP-INCA-OH'
        head(16) = 'molec/cm3'
        write(head(25),'(I16)') itime
        head(26) = 'MONTH'
        write(head(27),'(I4.4,I2.2,I2.2,1X,I2.2,I2.2,I2.2)') (jdate(i),i=1,6)
        write(head(28),'(I16)') 0
        head(29) = 'GLON128'
        write(head(30),'(I16)') 1
        write(head(31),'(I16)') imax
        head(32) = 'GGLA64'
        write(head(33),'(I16)') 1
        write(head(34),'(I16)') jmax
        head(35) = 'CSIG67'
        write(head(36),'(I16)') 1
        write(head(37),'(I16)') kmax
        head(38) = 'UR4'
        do i=39,43
            write(head(i),'(E16.7)') -999.0_real32
        enddo
        write(head(44),'(I16)') 1
        write(head(46),'(I16)') 0
        write(head(47),'(E16.7)') 0.0_real32
        write(head(48),'(I16)') 0
        head(63) = 'Dmitry_Belikov'
        write(head(64),'(I16)') imax*jmax*kmax

        write(unit,iostat=ios,iomsg=msg) head
        if (ios /= 0) then
            print *, 'GT3 header write error: ',trim(msg)
            stop 1
        endif
        write(unit,iostat=ios,iomsg=msg) vals
        if (ios /= 0) then
            print *, 'GT3 data write error: ',trim(msg)
            stop 1
        endif
    end subroutine write_gt3_record

    ! =============================================================
    ! === Convert complete NetCDF file, one month at a time
    subroutine convert(input_file,output_file,start_year)
        character(*), intent(in) :: input_file,output_file
        integer, intent(in) :: start_year
        integer :: ncid, varid, nd, nmonths, i, j, t, year, month
        integer :: xdid,ydid,zdid,tdid,dimids(nf90_max_var_dims), ndims, ios, outun
        real(real32) :: lon(imax),lat(jmax),sigma(kmax)
        real(real32), allocatable :: frame(:,:,:,:), vals(:,:,:)
        logical :: reverse_lat
        character(512) :: msg
        character(64) :: units

        print '(2A)', 'Input NetCDF : ',trim(input_file)
        print '(2A)', 'Output GT3   : ',trim(output_file)
        print '(A,I0)', 'Start year   : ',start_year

        call check(nf90_open(trim(input_file),nf90_nowrite,ncid),'open')
        call check(nf90_inq_dimid(ncid,'x',xdid),'x dimension')
        call check(nf90_inq_dimid(ncid,'y',ydid),'y dimension')
        call check(nf90_inq_dimid(ncid,'level',zdid),'level dimension')
        call check(nf90_inq_dimid(ncid,'time',tdid),'time dimension')
        call check(nf90_inquire_dimension(ncid,xdid,len=nd),'x length')
        if (nd /= imax) stop 'Expected 128 longitudes'
        call check(nf90_inquire_dimension(ncid,ydid,len=nd),'y length')
        if (nd /= jmax) stop 'Expected 64 latitudes'
        call check(nf90_inquire_dimension(ncid,zdid,len=nd),'level length')
        if (nd /= kmax) stop 'Expected 67 levels'
        call check(nf90_inquire_dimension(ncid,tdid,len=nmonths),'time length')
        if (nmonths <= 0 .or. mod(nmonths,nmon) /= 0) stop 'Time must contain full years'

        call check(nf90_inq_varid(ncid,'x',varid),'longitude variable')
        call check(nf90_get_var(ncid,varid,lon),'longitude values')
        call check(nf90_inq_varid(ncid,'y',varid),'latitude variable')
        call check(nf90_get_var(ncid,varid,lat),'latitude values')
        call check(nf90_inq_varid(ncid,'level',varid),'sigma variable')
        call check(nf90_get_var(ncid,varid,sigma),'sigma values')
        do i=1,imax
            if (abs(real(lon(i),real64)-lonstep*(i-1)) > 0.001_real64) then
                stop 'Longitude grid mismatch: expected 0 .. 357.1875'
            endif
        enddo
        if (.not.all(ieee_is_finite(lat))) stop 'Invalid latitude coordinate'
        if (.not.all(ieee_is_finite(sigma))) stop 'Invalid sigma coordinate'
        if (any(sigma <= 0.0_real32) .or. any(sigma > 1.0_real32)) stop 'Invalid sigma levels'
        if (any(abs(lat) > 90.0_real32)) stop 'Invalid latitudes'
        reverse_lat = lat(1) < lat(jmax)
        if (.not.reverse_lat .and. .not.(lat(1) > lat(jmax))) stop 'Invalid latitude ordering'
        do j=2,jmax
            if (reverse_lat .and. lat(j) <= lat(j-1)) stop 'Non-monotonic latitude grid'
            if (.not.reverse_lat .and. lat(j) >= lat(j-1)) stop 'Non-monotonic latitude grid'
        enddo

        call check(nf90_inq_varid(ncid,'OH',varid),'OH variable')
        dimids=0
        call check(nf90_inquire_variable(ncid,varid,ndims=ndims,dimids=dimids),'OH dimensions')
        if (ndims /= 4) stop 'Expected OH(time, level, y, x)'
        if (any(dimids(1:4) /= (/xdid,ydid,zdid,tdid/))) stop 'Unexpected OH dimension order'
        units=''
        call check(nf90_get_att(ncid,varid,'units',units),'OH units')
        select case(trim(adjustl(units)))
        case ('molecules cm-3','molec/cm3','molecules/cm3')
            continue
        case default
            print *, 'Unexpected OH units: ',trim(units)
            stop 1
        end select

        allocate(frame(imax,jmax,kmax,1), vals(imax,jmax,kmax))

        open(newunit=outun,file=trim(output_file),form='unformatted', &
             status='new',action='write',convert='big_endian',iostat=ios,iomsg=msg)
        if (ios /= 0) then
            print *, 'Cannot create GT3: ',trim(msg)
            stop 1
        endif

        do t=1,nmonths
            year=start_year+(t-1)/12
            month=mod(t-1,12)+1
            call check(nf90_get_var(ncid,varid,frame, &
                 start=(/1,1,1,t/),count=(/imax,jmax,kmax,1/)), 'read monthly OH')
            if (reverse_lat) then
                vals=frame(:,jmax:1:-1,:,1)
            else
                vals=frame(:,:,:,1)
            endif
            if (.not.all(ieee_is_finite(vals))) then
                print *, 'NaN/Inf found: ',year,month
                stop 1
            endif
            if (any(vals < 0.0_real32)) then
                print *, 'Negative OH found: ',year,month
                stop 1
            endif
            call write_gt3_record(outun,vals,year,month,t-1)
            print '(A,I4.4,A,I2.2,A,ES11.3)', 'GT3: ',year,'-',month, ' mean OH=',sum(real(vals,real64))/size(vals)
        enddo
        close(outun)
        call check(nf90_close(ncid),'close')
        print '(A,I0)', 'GT3 records written: ',nmonths
    end subroutine convert

end module df0_inca_oh_nc2gt3_t42_mod

program main
    use df0_inca_oh_nc2gt3_t42_mod
    implicit none
    character(1024) :: input_file, output_file, arg
    integer :: yr, ios
    if (command_argument_count() /= 3) then
        print *, 'Usage: df0_inca_oh_nc2gt3_t42 input.nc output.gt3 start_year'
        stop 1
    endif
    call get_command_argument(1,input_file)
    call get_command_argument(2,output_file)
    call get_command_argument(3,arg)
    read(arg,*,iostat=ios) yr
    if (ios/=0 .or. yr<1900 .or. yr>2200) stop 'Invalid start year'
    call convert(trim(input_file),trim(output_file),yr)
end program main
