! ifort df1_write_prior_uncertainty.f90 $(nf-config --fflags) $(nf-config --flibs)
! rm ../inv_dir/sourceprior/sourcepriorch4_p*
!
  program write_sourcetdi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    include 'netcdf.inc'

    ![input/output parameter]
    ! --- Runtime paths; prefix must not include an uncertainty percentage.
    character(len=1024)                 :: un_file, nc_file, config_file
    character(len=16)                   :: hcase
    integer                             :: cfg_s_year, cfg_f_year, cfg_ios
    namelist /d1_write_prior_uncertainty/ hcase, un_file, nc_file, cfg_s_year, cfg_f_year

    ![parameter]
    integer,parameter                  :: s_year = 1998, f_year = 2026, nyear=f_year-s_year+1
    integer,parameter                  :: ndattype=2, lreg=54, mtot=12, npresub=4, oreg=1
    real, dimension(:,:,:), allocatable :: s_presub
    real, dimension(:,:,:), allocatable :: s_land!,s_ocean
    real, dimension(:),     allocatable :: s_offset,time
    integer                            :: ncid,ok,npresubid,ndatid,lregid,monthid,nameid,lenid
    integer                            :: npresubp1id,spresubid,slandid,vdims1(1),vdims2(2),vdims3(3) !soceanid
    integer                            :: start(3),count(3),timeid,soffsetid,yearid !oregid
    character*10                       :: presubname(npresub+1), dum, un_percent_str
    character*512                      :: fname,msg,line
    real                               :: srcl(lreg),uncl(lreg,mtot),row(mtot)
    integer                            :: un_percent,i,j,isy,iy,im,ju,ios
    logical                            :: seen(lreg)
    integer,save,dimension(3)           :: un_percents = (/30, 50, 99/)
    data                                  presubname/' anth',' biob',' othr',' wetl', ' ch4off'/
    data                                  srcl / lreg * 0.0 /

    ! --- Read runtime configuration
    hcase = ''
    un_file = ''
    nc_file = ''
    config_file = ''
    cfg_s_year = s_year
    cfg_f_year = f_year
    call get_command_argument(1,config_file)
    if (len_trim(config_file) == 0) error stop 'Usage: df1_write_prior_uncertainty <d1_write_prior_uncertainty.nml>'
    open(10,file=trim(config_file),status='old',action='read',iostat=cfg_ios)
    if (cfg_ios /= 0) error stop 'Cannot open D1 uncertainty namelist'
    read(10,nml=d1_write_prior_uncertainty,iostat=cfg_ios)
    close(10)
    if (cfg_ios /= 0) error stop 'Invalid D1 uncertainty namelist'
    if (trim(hcase) /= 'CYC' .and. trim(hcase) /= 'INCA') error stop 'Invalid D1 uncertainty case'
    if (cfg_s_year /= s_year .or. cfg_f_year /= f_year) error stop 'D1 uncertainty period must remain 1998-2026'
    if (len_trim(un_file) == 0 .or. len_trim(nc_file) == 0) error stop 'D1 uncertainty input/output path is empty'
    print '(2A)', '      Case          : ', trim(hcase)
    print '(2A)', '      Regional table: ', trim(un_file)
    print '(2A)', '      NetCDF prefix : ', trim(nc_file)

    allocate(s_presub(npresub,nyear,ndattype))
    s_presub=0.0
    allocate(s_land(lreg,mtot*nyear,ndattype))
    s_land=0.0
    allocate(s_offset(ndattype))
    s_offset=0.0
    allocate(time(mtot*nyear))
    time=0.0

!   need to read in or otherwise set values
!   need to fix fossil to increase with year
    isy = s_year-1980
    print *, 'Years parameters ', isy, isy+nyear-1, nyear
    print *, 'Calendar years ', s_year, f_year
    s_presub(:,:,1) = 0.0
    s_offset(1) = 0.0
    s_presub(1,:,2) = 0.01
    s_presub(2,:,2) = 0.01
    s_presub(3,:,2) = 0.01
    s_presub(4,:,2) = 0.01
    s_offset(2) = 0.0
!
    ! f4 input: monthly climatological net CH4 masses [Tg/month].
    uncl = 0.0
    seen = .false.
    fname = un_file
    open(50,file=un_file,status='old',form='formatted',action='read',iostat=ios,iomsg=msg)
    call check_io(ios,'Open regional input')
    read(50,'(A)',iostat=ios,iomsg=msg) line
    call check_io(ios,'Read regional header')
    do i = 1, lreg
      read(50,'(A)',iostat=ios,iomsg=msg) line
      call check_io(ios,'Read regional row')
      read(line,*,iostat=ios,iomsg=msg) j,row
      call check_io(ios,'Parse regional row')
      if (j < 1 .or. j > lreg) stop 'Invalid region number'
      if (seen(j)) stop 'Duplicate region number'
      if (.not.all(ieee_is_finite(row))) stop 'Non-finite regional flux'
      seen(j) = .true.
      uncl(j,:) = row
      print '(i5,15f9.2)', j, row
    end do
    do
      read(50,'(A)',iostat=ios,iomsg=msg) line
      if (is_iostat_end(ios)) exit
      call check_io(ios,'Check end of regional input')
      if (len_trim(line) /= 0) stop 'Unexpected extra regional input row'
    end do
    close(50,iostat=ios,iomsg=msg)
    call check_io(ios,'Close regional input')
    if (.not.all(seen)) stop 'Missing region number'
    ! Preserve the original annual-equivalent convention: Tg/month * 12.
    ! Uncertainty is non-negative even where the net regional flux is negative.
    uncl = abs(uncl)*12.0

!      do i = 1, lreg; j = 0; do iy = 1, nyear; do im = 1, mtot
!         j = j + 1
!!        s_land (i,j,1) = srcl(i); s_land (i,j,2) = uncl(i,im) * 1.0 !Inv07_ case
!         s_land (i,j,1) = srcl(i); s_land (i,j,2) = uncl(i,im) * 0.7 !Inv08_ case
!      end do; end do; end do

    iy = 0
    do i = 1, nyear
    do j = 1, mtot
      iy = iy+1
      time(iy) = s_year+float(i-1)+(float(j)-.5)/float(mtot)
    enddo
    enddo
    print *, 'Time points: ', size(time), ' first/last: ', time(1),time(size(time))

    do ju = 1, size(un_percents)
      un_percent = un_percents(ju)
      if (un_percent < 0 .or. un_percent > 100) stop 'Invalid uncertainty percentage'
      print *, 'Run for ', ju, un_percent
      do i = 1, lreg
        j = 0
        do iy = 1, nyear
        do im = 1, mtot
          j = j+1
          s_land(i,j,1) = srcl(i)
          s_land(i,j,2) = uncl(i,im)*real(un_percent)/100.0
        end do
        end do
      end do
      print *, 'Land prior correction min/max: ', minval(s_land(:,:,1)),maxval(s_land(:,:,1))
      print *, 'Land uncertainty min/max [Tg/year equivalent]: ', minval(s_land(:,:,2)),maxval(s_land(:,:,2))
      write(un_percent_str,'(I0)') un_percent
      fname = trim(nc_file)//trim(un_percent_str)//'.nc'
      ! Preserve the original protection against overwriting existing files.
      call check_nc(nf_create(trim(fname),nf_noclobber,ncid),'Create file')
!
      call check_nc(nf_def_dim(ncid,'presub_count',npresub,npresubid),'Define presub_count')
      call check_nc(nf_def_dim(ncid,'presubp1',npresub+1,npresubp1id),'Define presubp1')
      call check_nc(nf_def_dim(ncid,'datatype_count',ndattype,ndatid),'Define datatype_count')
      call check_nc(nf_def_dim(ncid,'land_region',lreg,lregid),'Define land_region')
  !   ok = nf_def_dim(ncid,'ocean_region',oreg,oregid)
  !   if (ok.ne.nf_noerr) write(6,*) 'def_dim 4: ',nf_strerror(ok)
      call check_nc(nf_def_dim(ncid,'month',mtot*nyear,monthid),'Define month')
      call check_nc(nf_def_dim(ncid,'year',nyear,yearid),'Define year')
      call check_nc(nf_def_dim(ncid,'namelen',10,lenid),'Define namelen')
!
      vdims2(1)=lenid
      vdims2(2)=npresubp1id
      call check_nc(nf_def_var(ncid,'presub_name',nf_char,2,vdims2,nameid),'Define presub_name')
      vdims1(1)=monthid
      call check_nc(nf_def_var(ncid,'time',nf_float,1,vdims1,timeid),'Define time')
      vdims1(1)=ndatid
      call check_nc(nf_def_var(ncid,'sources_offset',nf_float,1,vdims1,soffsetid),'Define sources_offset')
      vdims3(1)=npresubid
      vdims3(2)=yearid
      vdims3(3)=ndatid
      call check_nc(nf_def_var(ncid,'sources_presub',nf_float,3,vdims3,spresubid),'Define sources_presub')
!
      vdims3(1)=lregid
      vdims3(2)=monthid
      vdims3(3)=ndatid
      call check_nc(nf_def_var(ncid,'sources_land',nf_float,3,vdims3,slandid),'Define sources_land')
  !   vdims3(1)=oregid
  !   ok = nf_def_var(ncid,'sources_ocean',nf_float,3,vdims3,soceanid)
  !   if (ok.ne.nf_noerr) write(6,*) 'def_var 3:',nf_strerror(ok)
!
      call put_text(timeid,'long_name','decimal year at nominal mid-month')
      call put_text(slandid,'units','Tg CH4/year')
      call put_text(slandid,'datatype_order','1: prior correction; 2: prior uncertainty')
      call put_text(slandid,'uncertainty_basis','abs(monthly climatological net mass) * 12 * percentage / 100')
      call check_nc(nf_enddef(ncid),'End definitions')
!
      call check_nc(nf_put_var_text(ncid,nameid,presubname),'Write presub_name')
      call check_nc(nf_put_var_real(ncid,timeid,time),'Write time')
      call check_nc(nf_put_var_real(ncid,soffsetid,s_offset),'Write sources_offset')
      start=1
      count(1)=npresub
      count(2)=nyear
      count(3)=ndattype
      call check_nc(nf_put_vara_real(ncid,spresubid,start,count,s_presub),'Write sources_presub')
      count(1)=lreg
      count(2)=mtot*nyear
      count(3)=ndattype
      call check_nc(nf_put_vara_real(ncid,slandid,start,count,s_land),'Write sources_land')
  !   count(1)=oreg
  !   ok = nf_put_vara_real(ncid,soceanid,start,count,s_ocean)
  !   if (ok.ne.nf_noerr) write(6,*) 'put_vara 3: ',nf_strerror(ok)
  !
      call check_nc(nf_close(ncid),'Close file')
      print *, 'Written: ', trim(fname)
    enddo
    deallocate(s_presub,s_land,s_offset,time)
    print *, 'Program finish OK'
  contains
    Subroutine check_nc(status,action)
      integer,intent(in) :: status
      character(len=*),intent(in) :: action
      if (status == nf_noerr) return
      print *, trim(action), ': ', trim(fname)
      print *, nf_strerror(status)
      stop 1
    End subroutine check_nc

    Subroutine check_io(status,action)
      integer,intent(in) :: status
      character(len=*),intent(in) :: action
      if (status == 0) return
      print *, trim(action), ': ', trim(fname)
      print *, trim(msg)
      stop 1
    End subroutine check_io

    Subroutine put_text(varid,name,value)
      integer,intent(in) :: varid
      character(len=*),intent(in) :: name,value
      call check_nc(nf_put_att_text(ncid,varid,name,len_trim(value),trim(value)),'Write '//name)
    End subroutine put_text
  end program write_sourcetdi
