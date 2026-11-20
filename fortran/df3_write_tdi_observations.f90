! ifort df3_write_tdi_observations.f90 $(nf-config --fflags) $(nf-config --flibs)
!
  Module df3_write_tdi_observations_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    include 'netcdf.inc'
    
    ![i/o files]
    ! --- Runtime paths; the output summary must be case-specific to avoid collisions.
    character(len=1024)                     :: obs_file, o_obsrv, i_infit, o_prior, o_pr_nc, config_file
    character(len=16)                       :: hcase
    integer                                 :: cfg_s_year, cfg_f_year, cfg_ios
    namelist /d3_write_tdi_observations/ hcase, obs_file, o_obsrv, i_infit, o_prior, o_pr_nc, cfg_s_year, cfg_f_year

    ![parameters]
    integer, parameter                      :: s_year = 1999, f_year = 2025, nyear = f_year-s_year+1
    integer, parameter                      :: nmon   = 12
    integer, parameter                      :: nsite  = 245
    integer, parameter                      :: ntime  = 9862 ! from df3_combine_obs_actm.f90
    integer, parameter                      :: n_mlo  = 34 
    integer, parameter                      :: ntrac  = 2, ntype=2, ndat2=nsite
    integer, parameter                      :: tlen   = nmon*nyear

    character(len=8),parameter              :: tracer(ntrac) = [character(len=8) :: 'ch4obs','ch4_inv1']
    character(len=3),parameter              :: unxcases(3)   = ['ctl','ux2','ux4']
    integer,save                            :: scales(3)     = (/1, 2, 4/)
    integer,save                            :: nday(nmon)    = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    ![vars]
    real, dimension(:),allocatable          :: lat, lon, alt
    real, dimension(:,:,:),allocatable      :: newdat1, newdat
    integer, dimension(:), allocatable      :: isiteno, itime
    character*45, dimension(:), allocatable :: sitenames

    real                                    :: tvect(tlen), off1
    real*8                                  :: sum2(nmon), summ
    real                                    :: stdev(ndat2,nmon), fits(3,nsite,ntrac,ntime), tr_mm(3,ndat2,ntrac,tlen)
    integer                                 :: vdims3(3), start3(3), count3(3), vdims2(2), typeid, nameid, lenid
    integer                                 :: sid, tid, tracdatid, ncid, statcntid, ierr, tlenid
    integer                                 :: nmis(3,ndat2,ntrac), nd2(nmon), ndat, ndat1, n1
    integer                                 :: iyear,junc,modl,n,l,k,i,j,it,nd,nm,iy,im,nda,mday,md,nfl
    character(len=512)                      :: fileout
    character(len=65)                       :: siten(nsite)

    character*512                           :: filename, msg, prior_file, obsrv_file
    integer                                 :: ios, nrecl, expected_days, nvalid, vdims1(1)
    real,parameter                          :: missing = -999.0
    real                                    :: summary(3)
    logical                                 :: valid_month(tlen)



  contains

    !*****************************************************************
    !*** Read runtime configuration (the inversion data calculations are unchanged)
    Subroutine read_config()
      hcase = ''
      obs_file = ''
      o_obsrv = ''
      i_infit = ''
      o_prior = ''
      o_pr_nc = ''
      config_file = ''
      cfg_s_year = s_year
      cfg_f_year = f_year
      call get_command_argument(1,config_file)
      if (len_trim(config_file) == 0) error stop 'Usage: df3_write_tdi_observations <d3_write_tdi_observations.nml>'
      open(10,file=trim(config_file),status='old',action='read',iostat=cfg_ios)
      if (cfg_ios /= 0) error stop 'Cannot open D3 observations namelist'
      read(10,nml=d3_write_tdi_observations,iostat=cfg_ios)
      close(10)
      if (cfg_ios /= 0) error stop 'Invalid D3 observations namelist'
      if (trim(hcase) /= 'CYC' .and. trim(hcase) /= 'INCA') error stop 'Invalid D3 observations case'
      if (cfg_s_year /= s_year .or. cfg_f_year /= f_year) error stop 'D3 observations period must remain 1999-2025'
      if (len_trim(obs_file) == 0 .or. len_trim(o_obsrv) == 0 .or. &
          len_trim(i_infit) == 0 .or. len_trim(o_prior) == 0 .or. &
          len_trim(o_pr_nc) == 0) error stop 'D3 observations input/output path is empty'
      print '(2A)', '      Case             : ', trim(hcase)
      print '(2A)', '      Site list        : ', trim(obs_file)
      print '(2A)', '      Filtered daily   : ', trim(i_infit)
      print '(2A)', '      Observation stats: ', trim(o_obsrv)
      print '(2A)', '      Prior difference : ', trim(o_prior)
      print '(2A)', '      NetCDF prefix    : ', trim(o_pr_nc)
    End subroutine read_config

    !*****************************************************************
    !*** write_nc
    Subroutine write_nc

      fileout = trim(o_pr_nc)//unxcases(junc)//tracer(modl)(5:8)//'.nc'
      print*, '  Write nc file ', fileout
      
      filename = fileout
      ierr = nf_create(trim(fileout),nf_clobber,ncid)
      call check_nc(ierr,'nf_create')
!
      ierr = nf_def_dim(ncid,'station_count',ndat1,statcntid)
      call check_nc(ierr,'nf_def_dim')
      ierr = nf_def_dim(ncid,'tlen',tlen,tlenid)
      call check_nc(ierr,'nf_def_dim')
      ierr = nf_def_dim(ncid,'data_type',ntype,typeid)
      call check_nc(ierr,'nf_def_dim')
      ierr = nf_def_dim(ncid,'namelen',45,lenid)
      call check_nc(ierr,'nf_def_dim')

      print *, '    nc variables assigned '

      vdims2(1)=lenid; vdims2(2)=statcntid
      ierr = nf_def_var(ncid,'site_name',nf_char,2,vdims2,nameid)
      call check_nc(ierr,'nf_def_var')
      vdims1(1)=statcntid
      ierr = nf_def_var(ncid,'station_count',NF_SHORT,1,vdims1,sid)
      call check_nc(ierr,'nf_def_var')
      vdims1(1)=tlenid
      ierr = nf_def_var(ncid,'tlen',NF_FLOAT,1,vdims1,tid)
      call check_nc(ierr,'nf_def_var')
      vdims3(2) = statcntid; vdims3(3) = tlenid; vdims3(1) = typeid
      ierr = nf_def_var(ncid,'tracer_data',NF_FLOAT,3,vdims3,tracdatid)
      call check_nc(ierr,'nf_def_var')
!
      ierr = NF_PUT_ATT_TEXT(ncid,sid,'long_name',14,'station number')
      call check_nc(ierr,'NF_PUT_ATT_TEXT')
      ierr = nf_put_att_text(ncid,tracdatid,'long_name',26,'monthly OBS-ACTM and error')
      call check_nc(ierr,'nf_put_att_text')
      ierr = nf_put_att_text(ncid,tracdatid,'units',3,'ppb')
      call check_nc(ierr,'nf_put_att_text')
!
      ierr = nf_put_att_real(ncid,tracdatid,'_FillValue',nf_float,1,missing)
      call check_nc(ierr,'Set missing value')
      ierr = nf_put_att_text(ncid,tracdatid,'data_type_order',27,'1: OBS-ACTM; 2: uncertainty')
      call check_nc(ierr,'Set data types')
      ierr = nf_enddef(ncid)
      call check_nc(ierr,'nf_enddef')
!
      ierr = nf_put_var_text(ncid,nameid,sitenames)
      call check_nc(ierr,'nf_put_var_text')

      ierr = nf_put_var_int(ncid,sid,isiteno)
      call check_nc(ierr,'nf_put_var_int')
      ierr = nf_put_var_real(ncid,tid,tvect)
      call check_nc(ierr,'nf_put_var_real')
!
      start3(1) = 1; start3(2) = 1; start3(3) = 1
      count3(1) = 2; count3(2) = ndat1; count3(3) = tlen
      print *, '    nc counts ', count3
      allocate (newdat1(2,ndat1,tlen))
      newdat1=0.0
      ! it is safer to reinitialize the array before writing
      newdat1 (:,:,:) = newdat (1:2,1:ndat1,:)
      ierr = nf_put_vara_real(ncid,tracdatid,start3,count3,newdat1)
      call check_nc(ierr,'nf_put_vara_real')

      deallocate (newdat1)
!
      print *, '    nc data written'

      ierr = nf_close(ncid)
      call check_nc(ierr,'nf_close')
    End subroutine write_nc

    !*****************************************************************
    !*** read_sites
    Subroutine read_sites
      character*512 :: dum
      integer :: n, i, j
      filename = obs_file
      open(51,file=obs_file,status='old',action='read',iostat=ios,iomsg=msg)
      call check_io(ios,'Open site list')
      read(51,'(A)',iostat=ios,iomsg=msg) dum
      call check_io(ios,'Read site header 1')
      read(51,'(A)',iostat=ios,iomsg=msg) dum
      call check_io(ios,'Read site header 2')
      do n = 1, ndat2
        read(51,*,iostat=ios,iomsg=msg) i,j,siten(n),lat(n),lon(n),alt(n)
        call check_io(ios,'Read site row')
        if (i <= 0 .or. i > 32767) stop 'Site number outside NetCDF short range'
        if (any(isiteno(1:n-1) == i)) stop 'Duplicate site number'
        if (.not.all(ieee_is_finite((/lat(n),lon(n),alt(n)/)))) stop 'Non-finite site coordinates'
        if (abs(lat(n)) > 90.0) stop 'Invalid site latitude'
        isiteno(n) = i
      end do
      close(51,iostat=ios,iomsg=msg)
      call check_io(ios,'Close site list')
        !print*, i,j,siten(n)
      !print *, '  Read sitelist len = ', ndat2 
      !pause
    End subroutine read_sites

    !*****************************************************************
    !*** write_prior
    Subroutine write_prior
      filename = prior_file
      open(52,file=prior_file,status='replace',iostat=ios,iomsg=msg) !CH4 prior file in ascii
      call check_io(ios,'Open prior differences')
      do n = 1, ndat1
        write(52,'(i5,2x,a45,2x,9f9.2)',iostat=ios,iomsg=msg) n,sitenames(n),lat(n),lon(n),alt(n)
        call check_io(ios,'Write prior station')
        valid_month = ieee_is_finite(newdat(1,n,:)) .and. ieee_is_finite(newdat(2,n,:))
        valid_month = valid_month .and. newdat(1,n,:) /= missing
        valid_month = valid_month .and. newdat(2,n,:) > 0.0 .and. newdat(2,n,:) <= 20.0
        nvalid = count(valid_month)
        write(52,*,iostat=ios,iomsg=msg) nvalid
        call check_io(ios,'Write prior record count')
        do nm = 1, tlen
          if (.not.valid_month(nm)) cycle
          write(52,'(i5,f12.5,4f12.4)',iostat=ios,iomsg=msg) nm,tvect(nm),(newdat(j,n,nm),j=1,3)
          call check_io(ios,'Write prior monthly row')
        end do
      end do
      close(52,iostat=ios,iomsg=msg)
      call check_io(ios,'Close prior differences')
      print *, '  Prior data written: ', trim(prior_file)
    End subroutine write_prior

    !*****************************************************************
    !*** read_fits (synoptic variations - daily data)) 
    Subroutine read_fits
      expected_days = 0
      do iy = s_year, f_year
      do im = 1, nmon
        expected_days = expected_days+month_days(iy,im)
      end do
      end do
      if (expected_days /= ntime) stop 'ntime does not match calendar period'
      filename = i_infit
      inquire(iolength=nrecl) (((fits(k,n,j,1),n=1,nsite),k=1,3),j=1,ntrac)
      open(50,file=i_infit,access='direct',form='unformatted',status='old',recl=nrecl,iostat=ios,iomsg=msg)
      call check_io(ios,'Open daily fits')
      do l = 1, ntime       !  Daily values    ! j = 1: for Observations - is of our interests here
        read(50,rec=l,iostat=ios,iomsg=msg) (((fits(k,n,j,l),n=1,nsite),k=1,3),j=1,ntrac)
        if (ios /= 0) print *, 'Daily record: ', l
        call check_io(ios,'Read daily fits')
      end do
      close(50,iostat=ios,iomsg=msg)
      call check_io(ios,'Close daily fits')
      print *, 'Read daily fits: ', ntime, ' records, ', nsite, ' sites'
      !print *, '  Read fitted file ', i_infit, '  # of records ', l-1
      !print*, fits(1, 5, 1, :)
    End subroutine read_fits

    !*****************************************************************
    !*** monthly-means and monthly-stdevs for each sites
    Subroutine get_mean2std
      ! monthly-means for each sites
      nmis = 0
      tr_mm = missing
      do k = 1, 3
      do n = 1, ndat2
        nd = 0
        nm = 0
        do iy = 1, nyear
          iyear = s_year+iy-1
          do im = 1, nmon
            mday = month_days(iyear,im)
            nm = nm+1
            do j = 1, ntrac
              summ = 0.D0
              nda = 0
              do md = 1, mday
                !for valid obs
                ! Use the same valid days for observations and model.
                if (.not.all(valid_conc(fits(k,n,:,nd+md)))) cycle
                nda = nda+1
                summ = summ+dble(fits(k,n,j,nd+md))
              end do
              if (nda >= 1) then
                tr_mm(k,n,j,nm) = real(summ/dble(nda))
              else
                nmis(k,n,j) = nmis(k,n,j)+1
              end if
            end do
            nd = nd+mday
          end do
        end do
        if (nd /= ntime) stop 'Monthly averaging calendar mismatch'
      end do
      end do
      print *, '    monthly paired means done, nm= ', nm
      ! monthly-stdevs for each sites
      ! This is RMS(raw observation - fitted observation), pooled by calendar month.
      stdev = missing
      do n = 1, ndat2
        nd = 0
        sum2 = 0.D0
        nd2 = 0
        do iy = 0, nyear-1
          iyear = s_year+iy
          do im = 1, nmon
            mday = month_days(iyear,im)
            do md = 1, mday
              nd = nd+1
              if (.not.valid_conc(fits(1,n,1,nd))) cycle
              if (.not.valid_conc(fits(2,n,1,nd))) cycle
              nd2(im) = nd2(im)+1
              
              if (n == 42 .and. im == 12) then
                if (abs(fits(1,n,1,nd)-fits(2,n,1,nd)) > 20.0) then
                  print '(A,3I5,3F12.3)', 'SPO obs raw/fit/dif: ', iyear, im, md, &
                    fits(1,n,1,nd), fits(2,n,1,nd), fits(1,n,1,nd)-fits(2,n,1,nd)
                end if
              end if
              
              sum2(im) = sum2(im)+(dble(fits(1,n,1,nd))-dble(fits(2,n,1,nd)))**2
            end do
          end do
        end do
        if (nd /= ntime) stop 'Residual RMS calendar mismatch'
        do im = 1, nmon
          if (nd2(im) >= 1) stdev(n,im) = real(sqrt(sum2(im)/dble(nd2(im))))
        end do
      end do
      print *, '    monthly residual RMS done - measurements only, in ppb'
            !print '(7i5,f10.2)', k, n, j, iy, im, id, nd, fits(k, n, j, nd)
            !for valid obs
      !n = 1
      !print*, siten(n), tr_mm (1,n,1,:)
      !pause
        !print '(a9,2x,9i5)', siten(n)(1:8),k, n, j, nd, nm, nmis(k,n,j)
        !print *, siten(n)(1:25), stdev(n,:), nd2
      !pause
    End subroutine get_mean2std
    !*****************************************************************
    integer Function month_days(year,month) result(days)
      integer,intent(in) :: year,month
      days = nday(month)
      if (month == 2) then
        if (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) days = 29
      end if
    End function month_days

    elemental logical Function valid_conc(value)
      real,intent(in) :: value
      valid_conc = ieee_is_finite(value) .and. value >= 1201.0 .and. value <= 5001.0
    End function valid_conc

    Subroutine check_nc(status,action)
      integer,intent(in) :: status
      character(len=*),intent(in) :: action
      if (status == nf_noerr) return
      print *, trim(action), ': ', trim(filename)
      print *, nf_strerror(status)
      stop 1
    End subroutine check_nc

    Subroutine check_io(status,action)
      integer,intent(in) :: status
      character(len=*),intent(in) :: action
      if (status == 0) return
      print *, trim(action), ': ', trim(filename)
      print *, trim(msg)
      stop 1
    End subroutine check_io
  End module df3_write_tdi_observations_mod
!
!*****************************************************************
!
  Program df3_write_tdi_observations
    use df3_write_tdi_observations_mod
    implicit none
    
    print*, 'Start program'
    call read_config()
    allocate(lat(ndat2)); allocate(lon(ndat2)); allocate(alt(ndat2))
    lat=0.0; lon=0.0; alt=0.0
    allocate(sitenames(ndat2))
    allocate(newdat(3,ndat2,tlen))
    newdat=missing
    allocate(isiteno(ndat2)); allocate(itime(ndat2)); isiteno=0; itime=0

    call read_sites()
    call read_fits()
    
    print '(2x, a, i5, a, i5, a, i5)', 'nyear = ', nyear, '  f_year - s_year = ', f_year - s_year, '  ntime = ', ntime

    call get_mean2std
    

    DO junc = 1,3 
    DO modl = 2, ntrac
      print *, ''
      print *, 'Start processing for ', junc, modl

      newdat = missing
      itime = 0
      ! Keep legacy ASCII filenames for ctl; other cases have distinct suffixes.
      prior_file = trim(o_prior)
      obsrv_file = trim(o_obsrv)
      if (junc /= 1) then
        prior_file = trim(o_prior)//'.'//unxcases(junc)
        obsrv_file = trim(o_obsrv)//'.'//unxcases(junc)
      end if
      filename = obsrv_file
      open(55,file=obsrv_file,status='replace',iostat=ios,iomsg=msg) !to be used in inversion
      call check_io(ios,'Open observation summary')
      write(55,'(A,I4,A,I4,A5,A10)',iostat=ios,iomsg=msg) ' CH4 data ',s_year,'-',f_year,unxcases(junc),tracer(modl)
      call check_io(ios,'Write summary title')
      write (55, '(2x,a,2x,a10,2x,a45,2x,9a10)',iostat=ios,iomsg=msg) ' # ', 'itime', 'siteName', &
                        'Latitude', 'Longitude', 'Altitude', &
                        'RealDataFrac', 'newdat1', 'newdat2', 'newdat3'
      call check_io(ios,'Write summary columns')

      !for the Station in the filelist
      n1 = 0
      do n = 1, ndat2
        n1 = n1 + 1 
        nm = 0
        do iy = 1, nyear
        iyear = s_year + iy - 1
        
        do im = 1, nmon 
          nm = nm + 1
          !Observed - Model
          if (valid_conc(tr_mm(2,n,1,nm)) .and. valid_conc(tr_mm(2,n,modl,nm))) then
            if (ieee_is_finite(stdev(n,im)) .and. stdev(n,im) >= 0.0) then
              newdat(1,n1,nm) = tr_mm(2,n,1,nm)-tr_mm(2,n,modl,nm)
              newdat(2,n1,nm) = stdev(n,im)*scales(junc)+2.00+3.00
              if (.not.valid_conc(tr_mm(1,n,1,nm))) then
                newdat(2,n1,nm) = newdat(2,n1,nm)*2.0 !when there is no observed values - avoid extrapolation
              end if
            end if
          end if
          !1:Observed; 2:fitted in tr_mm(1,), we use fitted
          if (valid_conc(tr_mm(1,n,1,nm)) .and. valid_conc(tr_mm(1,n,modl,nm))) then
            newdat(3,n1,nm) = tr_mm(1,n,1,nm)-tr_mm(1,n,modl,nm)
          end if

          ! to check one site
          !if ((n.eq.5) .and. ((iy .le. 3) .or. (iy .ge. nyear - 5))) print'(a5,2i6,10f12.2)', &
          if ((n.eq.42) .and. (iyear.le.2003) .and. (im.eq.12)) print'(a5,2i6,10f12.2)', &
            siten(n)(5:7), iyear, im, &
            tr_mm (1,n,1,nm), tr_mm (1,n,modl,nm), tr_mm (2,n,1,nm), tr_mm (2,n,modl,nm), &
            tr_mm (3,n,1,nm), tr_mm (3,n,modl,nm), &
            newdat(1,n1,nm), newdat(2,n1,nm), newdat(3,n1,nm)

        end do
        end do
        sitenames(n1) = siten(n)
        lat(n1) = lat(n)
        lon(n1) = lon(n)
        alt(n1) = alt(n)
        if (count(newdat(2,n1,:) > 0.0) >= tlen-10) itime(n1) = 9999  ! 'tlen':fitted data without gaps at the end
        do k = 1, 3
          valid_month = ieee_is_finite(newdat(k,n1,:)) .and. newdat(k,n1,:) /= missing
          nvalid = count(valid_month)
          summary(k) = missing
          if (nvalid > 0) summary(k) = real(sum(dble(newdat(k,n1,:)),mask=valid_month)/dble(nvalid))
        end do
        write(55,'(2I9,2X,A45,2X,9F10.2)',iostat=ios,iomsg=msg) n1,itime(n1),sitenames(n1), &
          lat(n1),lon(n1),alt(n1),real(tlen-nmis(1,n,1))/real(tlen),summary
        call check_io(ios,'Write observation summary row')
        !if ( nmis(1,n,1) .ge.  9 ) itime(n1) = 9999  ! 'tlen':when there are obseved value for all months
        !write(55,'(2i9,2x,a45,2x,9f10.2)') n1,itime(n1),sitenames(n1),lat(n1),lon(n1),alt(n1),(tlen-nmis(1,n,1))/float(tlen), &
        !                   sum(newdat(1,n1,:))/float(tlen), sum(newdat(2,n1,:))/float(tlen), sum(newdat(3,n1,:))/float(tlen)
      end do
      close(55,iostat=ios,iomsg=msg)
      call check_io(ios,'Close observation summary')
      ! pause

      ! time vector
      do nm = 1, tlen
        tvect(nm)=s_year+(float(nm)-0.5)/12.
      enddo
      ndat1 = n1
      print *, '  ndat1, ndat2 ', ndat1, ndat2

      ! where (newdat .le. -200.) newdat = -999.0
      call write_prior()
      call write_nc()
      newdat=missing
      itime=0
      print *, 'All done for junc, modl ', junc, modl
      
      !do n1 = 1, ndat1
      !  write(*,'(i5, 9f10.2)') n1, sum(newdat(1,n1,:))/float(tlen), sum(newdat(2,n1,:))/float(tlen), sum(newdat(3,n1,:))/float(tlen)
      !enddo
    END DO  ! For each modl
    END DO  ! For each junc
!
    print *, 'Program finish OK'
  End program df3_write_tdi_observations
