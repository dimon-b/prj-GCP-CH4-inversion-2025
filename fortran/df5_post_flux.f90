! ============================================================
! === D5: posterior NetCDF -> GT3 conservative conversion
!     Python streams monthly posterior 360x180 float32 fields
!     to Fortran stdin. No intermediate prior/post GRD is written.
!     Output GT3: 128x64 Gaussian, kg CH4 m-2 s-1.
! ============================================================
module d5_post_mod
  use iso_fortran_env, only: real32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: sp=real32, dp=real64
  integer, parameter :: nlon=360, nlat=180, imax=128, jmax=64
  real(dp), parameter :: pi=3.1415926535897932384626433832795_dp
  real(dp), parameter :: radius=6371000.0_dp, tgscale=1.0e-9_dp
  character(len=1024) :: out_gt3='', budget_file=''
  integer :: initial_year=1999, nmonths=324
  namelist /post_config/ out_gt3, budget_file, initial_year, nmonths
  real(dp) :: area_in(nlat), area_t42(jmax), latweight(jmax), sinedge(jmax+1)
  real(dp) :: wlon(imax,nlon), wlat(jmax,nlat), gaussian_lon(imax), gaussian_lat(jmax)
  real(dp) :: temp(imax,nlat), mapped(imax,jmax)
  real(sp) :: post(nlon,nlat), map32(imax,jmax)
contains
  subroutine die(message)
    character(len=*), intent(in) :: message
    print '(2a)', 'D5 ERROR: ',trim(message)
    error stop 1
  end subroutine

  subroutine io_check(status,what)
    integer, intent(in) :: status
    character(len=*), intent(in) :: what
    if(status/=0) call die(what)
  end subroutine

  integer function days_in_month(year,month) result(nd)
    integer,intent(in) :: year,month
    integer,parameter :: days(12)=[31,28,31,30,31,30,31,31,30,31,30,31]
    nd=days(month)
    if(month==2.and.mod(year,4)==0.and.(mod(year,100)/=0.or.mod(year,400)==0)) nd=29
  end function

  subroutine read_settings()
    character(len=1024) :: path
    integer :: ios,u
    call get_command_argument(1,path,status=ios)
    if(ios/=0.or.len_trim(path)==0) call die('Usage: df5_post_flux <config.nml>')
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    call io_check(ios,'Open config')
    read(u,nml=post_config,iostat=ios)
    call io_check(ios,'Read post_config')
    close(u)
    if(len_trim(out_gt3)==0.or.len_trim(budget_file)==0) call die('Missing output path')
    if(nmonths<1.or.nmonths>348.or.initial_year<1900) call die('Invalid period')
  end subroutine

  subroutine gauss_grid()
    integer::j,k,it,i,sh
    real(dp)::x0,x1,p0,p1,pn,deriv,lo,hi,overlap,dlon,dx,sl,sr,tl,tr
    real(dp)::weight(jmax)
    do j=1,(jmax+1)/2
      x0=cos(pi*(real(j,dp)-0.25_dp)/(real(jmax,dp)+0.5_dp))
      do it=1,100
        p0=1.0_dp
        p1=x0
        do k=2,jmax
          pn=((2.0_dp*k-1.0_dp)*x0*p1-(k-1.0_dp)*p0)/k
          p0=p1
          p1=pn
        end do
        deriv=real(jmax,dp)*(x0*p1-p0)/(x0*x0-1.0_dp)
        x1=x0-p1/deriv
        if(abs(x1-x0)<1.0e-14_dp) exit
        x0=x1
      end do
      if(it>100) call die('Gaussian grid iteration failed')
      x0=x1
      p0=1.0_dp
      p1=x0
      do k=2,jmax
        pn=((2.0_dp*k-1.0_dp)*x0*p1-(k-1.0_dp)*p0)/k
        p0=p1
        p1=pn
      end do
      deriv=real(jmax,dp)*(x0*p1-p0)/(x0*x0-1.0_dp)
      gaussian_lat(j)=x0
      gaussian_lat(jmax-j+1)=-x0
      weight(j)=2.0_dp/((1.0_dp-x0*x0)*deriv**2)
      weight(jmax-j+1)=weight(j)
    end do
    if(abs(sum(weight)-2.0_dp)>1.0e-12_dp) call die('Gaussian weights invalid')
    sinedge(1)=1.0_dp
    do j=1,jmax
      latweight(j)=weight(j)
      sinedge(j+1)=sinedge(j)-weight(j)
      area_t42(j)=2.0_dp*pi*radius**2*weight(j)/imax
    end do
    sinedge(jmax+1)=-1.0_dp
    wlat=0.0_dp
    do j=1,jmax
      do i=1,nlat
        lo=sin((-90.0_dp+real(i-1,dp))*pi/180.0_dp)
        hi=sin((-90.0_dp+real(i,dp))*pi/180.0_dp)
        overlap=max(0.0_dp,min(hi,sinedge(j))-max(lo,sinedge(j+1)))
        wlat(j,i)=overlap/weight(j)
      end do
    end do
    dlon=360.0_dp/imax
    dx=360.0_dp/nlon
    wlon=0.0_dp
    do j=1,imax
      gaussian_lon(j)=360.0_dp*real(j-1,dp)/imax
      tl=gaussian_lon(j)-0.5_dp*dlon
      tr=gaussian_lon(j)+0.5_dp*dlon
      do i=1,nlon
        do sh=-1,1
          sl=-180.0_dp+(real(i,dp)-1.0_dp)*dx+360.0_dp*sh
          sr=sl+dx
          wlon(j,i)=wlon(j,i)+max(0.0_dp,min(sr,tr)-max(sl,tl))/dlon
        end do
      end do
    end do
    if(maxval(abs(sum(wlat,dim=2)-1.0_dp))>1.0e-12_dp) call die('Latitude remap not conservative')
    if(maxval(abs(sum(wlon,dim=2)-1.0_dp))>1.0e-12_dp) call die('Longitude remap not conservative')
  end subroutine

  subroutine make_gt3_header(unit,year,month)
    integer,intent(in)::unit,year,month
    character(len=16)::head(64)
    integer::i,day,hr
    head=''
    write(head(1),'(i16)') 9010
    head(2)='GCP2025 post'
    head(3)='QFLUX07'
    write(head(12),'(i16)') 0
    write(head(13),'(i16)') 1
    head(14)='fch4_inv1'
    head(16)='kg-CH4/m2/s'
    write(head(25),'(i16)') (year-initial_year)*12+month-1
    head(26)='MONTH'
    day=1+days_in_month(year,month)/2
    hr=mod(days_in_month(year,month),2)*12
    write(head(27),'(i4.4,i2.2,i2.2,a1,i2.2,i2.2,i2.2)') year,month,day,' ',hr,0,0
    write(head(28),'(i16)') 1
    write(head(29),'(a,i3)') 'GLON',imax
    write(head(30),'(i16)') 1
    write(head(31),'(i16)') imax
    write(head(32),'(a,i2)') 'GGLA',jmax
    write(head(33),'(i16)') 1
    write(head(34),'(i16)') jmax
    head(35)='SFC1'
    write(head(36),'(i16)') 1
    write(head(37),'(i16)') 1
    head(38)='UR4'
    do i=39,43
      write(head(i),'(e16.7)') -999.0_dp
    end do
    write(head(44),'(i16)') 1
    write(head(46),'(i16)') 0
    write(head(47),'(e16.7)') 0.0_dp
    write(head(48),'(i16)') 0
    head(63)='GCP2025'
    write(head(64),'(i16)') imax*jmax
    write(unit) head
  end subroutine

  subroutine process_flux()
    integer :: ios,u_in,u_gt3,u_csv,m,j,year,month
    real(dp) :: seconds,source_mass,gt3_mass,abs_total,rel_error
    real(dp) :: longitude_area
    longitude_area=2.0_dp*pi*radius**2/real(nlon,dp)
    do j=1,nlat
      area_in(j)=longitude_area*(sin((-90.0_dp+real(j,dp))*pi/180.0_dp)- &
                     sin((-90.0_dp+real(j-1,dp))*pi/180.0_dp))
    end do
    open(newunit=u_in,file='/dev/stdin',form='unformatted',access='stream', &
         status='old',action='read',iostat=ios)
    call io_check(ios,'Open stdin float32 stream')
    open(newunit=u_gt3,file=trim(out_gt3),form='unformatted',convert='big_endian', &
         status='new',action='write',iostat=ios)
    call io_check(ios,'Open GT3 output')
    open(newunit=u_csv,file=trim(budget_file),status='new',action='write',iostat=ios)
    call io_check(ios,'Open budget CSV')
    write(u_csv,'(a)') 'year,month,input_Tg_month,gt3_Tg_month,gt3_rel_error'
    call gauss_grid()
    do m=1,nmonths
      year=initial_year+(m-1)/12
      month=mod(m-1,12)+1
      read(u_in,iostat=ios) post
      call io_check(ios,'Reading streamed monthly flux (short stream)')
      if(.not.all(ieee_is_finite(post))) call die('Nonfinite input posterior field')
      seconds=real(days_in_month(year,month),dp)*86400.0_dp
      source_mass=0.0_dp
      abs_total=0.0_dp
      do j=1,nlat
        source_mass=source_mass+sum(real(post(:,j),dp))*area_in(j)*seconds*tgscale
        abs_total=abs_total+sum(abs(real(post(:,j),dp)))*area_in(j)*seconds*tgscale
      end do
      temp=matmul(wlon,real(post,dp))
      mapped=matmul(temp,transpose(wlat))
      map32=real(mapped,sp)
      if(.not.all(ieee_is_finite(map32))) call die('Nonfinite remapped field')
      gt3_mass=0.0_dp
      do j=1,jmax
        gt3_mass=gt3_mass+sum(real(map32(:,j),dp))*area_t42(j)*seconds*tgscale
      end do
      rel_error=(gt3_mass-source_mass)/max(abs_total,1.0e-6_dp)
      if(abs(rel_error)>3.e-7_dp) call die('Conservative GT3 mass budget failure')
      call make_gt3_header(u_gt3,year,month)
      write(u_gt3,iostat=ios) map32
      call io_check(ios,'Write GT3 field')
      write(u_csv,'(i0,",",i0,3(",",es24.15))') year,month,source_mass,gt3_mass,rel_error
      if(month==12.or.m==nmonths) print '(a,i0,a,i0)', 'D5 GT3 year ',year,' month ',month
    end do
    close(u_in)
    close(u_gt3)
    close(u_csv)
  end subroutine
end module

program df5_post_flux
  use d5_post_mod
  implicit none
  call read_settings()
  call process_flux()
  print '(a)', 'D5 NETCDF -> GT3 (NO LOSS CORRECTION): PASS'
end program
