! ============================================================
! === Monthly GlobalView data: count and compact valid rows
!     Both passes use gv_load and exactly the same validity rule.
!     Existing station selection is by last matching 45-char name.
!     The order of active stations and GF indices is unchanged.
module gv_valid_months
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: gv_load
contains

    logical function gv_valid(obs, sigma)
        real, intent(in) :: obs, sigma
        gv_valid = .false.
        if (.not.ieee_is_finite(obs)) return
        if (.not.ieee_is_finite(sigma)) return
        if (obs == -999.0 .or. sigma == -999.0) return
        if (sigma <= 0.0) return
        gv_valid = .true.
    end function gv_valid

    subroutine gv_load(datf,nobs,neqns,monperyr,datnm,dat,datind, &
                       datyr,datmnth,cdat,ndat,datav,firstdat,lastdat, &
                       nread,count_only)
        implicit none
        include 'netcdf.inc'
        character(len=*), intent(in) :: datf
        integer, intent(in) :: nobs,neqns,monperyr,ndat,firstdat,lastdat
        character(len=80), intent(in) :: datnm(:)
        integer, intent(in) :: datav(:)
        real, intent(inout) :: dat(:),cdat(:)
        integer, intent(inout) :: datind(:),datyr(:),datmnth(:)
        integer, intent(out) :: nread
        logical, intent(in) :: count_only
        integer :: ncid,statid,timeid,typeid,nameid,namelenid
        integer :: obsid,sitevarid,tvarid,dimids(3),rc
        integer :: nsites,nmonths,ntypes,namelen,nstart
        integer :: i,j,k,l,m,isite,nmatches,ntime,nexpect
        real :: expected_time
        real, allocatable :: time(:), values(:,:,:)
        character(len=45), allocatable :: names(:)

        rc = nf_open(trim(datf),0,ncid)
        call gv_check(rc,'open '//trim(datf))
        rc = nf_inq_dimid(ncid,'station_count',statid)
        call gv_check(rc,'dim station_count')
        rc = nf_inq_dimlen(ncid,statid,nsites)
        call gv_check(rc,'len station_count')
        rc = nf_inq_dimid(ncid,'tlen',timeid)
        call gv_check(rc,'dim tlen')
        rc = nf_inq_dimlen(ncid,timeid,nmonths)
        call gv_check(rc,'len tlen')
        rc = nf_inq_dimid(ncid,'data_type',typeid)
        call gv_check(rc,'dim data_type')
        rc = nf_inq_dimlen(ncid,typeid,ntypes)
        call gv_check(rc,'len data_type')
        rc = nf_inq_dimid(ncid,'namelen',namelenid)
        call gv_check(rc,'dim namelen')
        rc = nf_inq_dimlen(ncid,namelenid,namelen)
        call gv_check(rc,'len namelen')
        if (nsites < nobs .or. ntypes /= 2 .or. namelen /= 45) then
            write(*,*) 'GV dimensions:',nsites,nmonths,ntypes,namelen
            error stop 'Invalid F7 NetCDF dimensions'
        endif

        rc = nf_inq_varid(ncid,'site_name',nameid)
        call gv_check(rc,'site_name')
        rc = nf_inq_varid(ncid,'tlen',tvarid)
        call gv_check(rc,'tlen variable')
        rc = nf_inq_varid(ncid,'tracer_data',obsid)
        call gv_check(rc,'tracer_data')
        rc = nf_inq_vardimid(ncid,obsid,dimids)
        call gv_check(rc,'tracer_data dimensions')
        if (dimids(1) /= typeid .or. dimids(2) /= statid .or. &
            dimids(3) /= timeid) error stop 'F7 tracer_data dimension order'

        allocate(time(nmonths),names(nsites),values(2,nsites,nmonths))
        rc = nf_get_var_real(ncid,tvarid,time)
        call gv_check(rc,'read time')
        rc = nf_get_var_text(ncid,nameid,names)
        call gv_check(rc,'read station names')
        rc = nf_get_var_real(ncid,obsid,values)
        call gv_check(rc,'read tracer_data')
        rc = nf_close(ncid)
        call gv_check(rc,'close F7 file')

        if (monperyr /= 12) error stop 'GV supports monthly observations only'
        nexpect = (lastdat-firstdat+1)*monperyr
        if (nexpect <= 0) error stop 'Invalid inversion period'
        nstart = 0
        do i=1,nmonths
            if (.not.ieee_is_finite(time(i))) cycle
            if (int(time(i)) /= firstdat) cycle
            if (abs(time(i)-(real(firstdat)+0.5/12.0)) < 0.015) then
                nstart = i
                exit
            endif
        enddo
        if (nstart == 0 .or. nstart+nexpect-1 > nmonths) then
            write(*,*) 'GV firstdat,lastdat,time length:',firstdat,lastdat,nmonths
            error stop 'F7 NetCDF does not cover inversion months'
        endif

        ! --- Fail on shifted or non-monthly time coordinates.
        do m=1,nexpect
            ntime = nstart+m-1
            expected_time = real(firstdat)+(real(m)-0.5)/12.0
            if (.not.ieee_is_finite(time(ntime))) error stop 'Nonfinite F7 time'
            if (abs(time(ntime)-expected_time) > 0.015) then
                write(*,*) 'Incorrect F7 time at month ',m, &
                           time(ntime),expected_time
                error stop 'F7 time index shifted'
            endif
        enddo

        nread=0
        do l=1,nobs
            if (datav(l) /= 9999) error stop 'Only monthly data allowed'
            isite=0
            nmatches=0
            do i=1,nsites
                if (datnm(l)(1:45) == names(i)) then
                    isite=i
                    nmatches=nmatches+1
                endif
            enddo
            if (isite == 0) then
                write(*,*) 'GV absent station: ',trim(datnm(l))
                error stop 'Requested station missing from F7 NetCDF'
            endif
            if (nmatches > 1 .and. count_only) then
                print '(A,A,A,I0,A,I0)', 'GV duplicate name: ', &
                      trim(datnm(l)), ' matches=', nmatches, &
                      ' last row=',isite
            endif
            do m=1,nexpect
                ntime=nstart+m-1
                if (.not.gv_valid(values(1,isite,ntime), &
                                  values(2,isite,ntime))) cycle
                nread=nread+1
                if (count_only) cycle
                if (nread > neqns) error stop 'GV count exceeds allocated rows'
                datind(nread)=l
                datyr(nread)=firstdat+(m-1)/12
                datmnth(nread)=mod(m-1,12)+1
                dat(nread)=values(1,isite,ntime)
                cdat(nread)=values(2,isite,ntime)**2
            enddo
        enddo
        if (count_only) then
            print '(A,I0,A,I0)', 'GV valid station-months=',nread, &
                                ' excluded=',nobs*nexpect-nread
        else
            print '(A,I0)', 'GV read valid station-months=',nread
        endif
    end subroutine gv_load

    subroutine gv_check(status,where)
        implicit none
        include 'netcdf.inc'
        integer, intent(in) :: status
        character(len=*), intent(in) :: where
        if (status == nf_noerr) return
        write(*,*) 'NetCDF error at ',trim(where),': ',nf_strerror(status)
        error stop 'NetCDF F7 read failed'
    end subroutine gv_check
end module gv_valid_months

! ============================================================
! === Count exactly the observations that read_globalview accepts
subroutine count_globalview(datf,datflist,nobs,monperyr, &
                            firstdat,lastdat,nvalid)
    use gv_valid_months, only: gv_load
    implicit none
    character(len=80), intent(in) :: datf,datflist
    integer, intent(in) :: nobs,monperyr,firstdat,lastdat
    integer, intent(out) :: nvalid
    integer :: nall,i,flag,station_id,ns,ios
    real :: lat,lon
    character(len=80) :: sname
    character(len=80), allocatable :: names(:)
    integer, allocatable :: av(:)
    real :: dummy_data(1),dummy_cdat(1)
    integer :: dummy_ind(1),dummy_year(1),dummy_month(1)

    allocate(names(nobs),av(nobs))
    open(unit=77,file=trim(datflist),status='old',action='read', &
         iostat=ios)
    if (ios /= 0) error stop 'Cannot open inversion datflist'
    read(77,*,iostat=ios)
    if (ios /= 0) error stop 'Bad datflist header'
    read(77,*,iostat=ios) nall
    if (ios /= 0 .or. nall < nobs) error stop 'Bad datflist count'
    ns=0
    do i=1,nall
        read(77,*,iostat=ios) station_id,flag,sname,lat,lon
        if (ios /= 0) error stop 'Malformed datflist row'
        if (flag == 0) cycle
        if (flag /= 9999) error stop 'Nonmonthly active station'
        ns=ns+1
        if (ns > nobs) error stop 'Active station count changed'
        names(ns)=sname
        av(ns)=flag
    enddo
    close(77)
    if (ns /= nobs) error stop 'Active station count mismatch'
    call gv_load(datf,nobs,1,monperyr,names,dummy_data,dummy_ind, &
                 dummy_year,dummy_month,dummy_cdat,nobs,av, &
                 firstdat,lastdat,nvalid,.true.)
end subroutine count_globalview

! ============================================================
! === Compact only valid station-month data into contiguous arrays
subroutine read_globalview(datf,nobs,neqns,monperyr,datnm,dat, &
                           datind,datyr,datmnth,cdat,ndat,datav, &
                           firstdat,lastdat,nread)
    use gv_valid_months, only: gv_load
    implicit none
    character(len=80), intent(in) :: datf
    integer, intent(in) :: nobs,neqns,monperyr,ndat,firstdat,lastdat
    character(len=80), intent(in) :: datnm(ndat)
    real, intent(inout) :: dat(neqns),cdat(neqns)
    integer, intent(inout) :: datind(neqns),datyr(neqns),datmnth(neqns)
    integer, intent(in) :: datav(ndat)
    integer, intent(out) :: nread
    call gv_load(datf,nobs,neqns,monperyr,datnm,dat,datind, &
                 datyr,datmnth,cdat,ndat,datav,firstdat,lastdat, &
                 nread,.false.)
end subroutine read_globalview
