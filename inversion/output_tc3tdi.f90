      subroutine output_tc3tdi( outputf, nunks, neqns, ngroups, nsources, ntotdat, &
          ndat, nlandreg, npresub, firstsrc,lastsrc,source, sourcetype, &
          sourceind, sourceyr, sourcemnth, sourcenm, csrc, dat, datind, datyr,     &
          datmnth, datav, datnm, cdat, grouplst, groupnm, grouplen, firstdat,      &
          lastdat, prefname, header)
! this routine handles the output for the tdi
! dumping values, variances and covariances for source and data
! see tdi.tex for general documentation
!     parameters:

! PKP has made some changes in the writing lately, search for PKP to see those
      implicit none
      character*80 outputf ! master output file name
      integer nunks,    &  ! total number of unknowns
           neqns,       &  ! number of equationsinten
           ngroups,     &  ! total number of groups
           nsources,    &  ! number of different source types (files)
           ntotdat,     &  ! total number of data (neqns - constraints)
           ndat            ! number of different types of data (files)
      real source(nunks)   ! source values
      integer sourceind(nunks)     ! index into source description arrays
      integer sourcetype(nsources) ! type of source
      integer sourceyr(nunks)      ! year of pulse
      integer sourcemnth(nunks)    ! month of pulse
      character*80 sourcenm(nsources) ! names
      real, dimension(1) :: csrc   ! source covar, now packed
      real dat(neqns)              ! data values
      integer datind(neqns)        ! index into data descriptions
      integer datyr(neqns)         ! year of datum
      integer datmnth(neqns)       ! month of datum
      integer datav(ndat)          ! averaging period of datum
      character*80 datnm(ndat)     ! data name
      real cdat(neqns)             ! covariance
      integer grouplst(nunks, ngroups) ! lists of sources in each group
      character*80 groupnm(ngroups) ! names of groups
      integer grouplen(ngroups)    ! sourcer per group
      integer firstdat             ! year of first datum
      integer lastdat              ! year of last datum
      integer firstsrc,lastsrc
      integer nlandreg,npresub
      character*(*) prefname       ! prefix name for output
!     character*80 prefname
      character*80 header          ! for output files
      character*80 curname
      logical done(nunks)
      integer kount,kount2
      include 'tdi.h'
      include 'common.h'

! name list
      character*80 blockf                     ! master filename for covariance blocks
      character*80 basis_file                 ! name of file for flux basis funcs
      integer max_disallowed                  ! maximum number of fluxes to exclude
      parameter(max_disallowed = 5)
      character*80 disallowed(max_disallowed) ! list of excluded flux names
      logical dosources, dodata, dogroups, docovar, donetcdf
      namelist /outputnml/ &
           basis_file, & ! file of flux basis functions
           blockf,     & ! name of master file for covariance blocks
           dosources,  & ! whether to dump sources
           dodata,     & ! flag to whether to dump data
           dogroups,   & ! flag for whether to dump groups
           docovar,    & ! flag for whether to dump covar
           donetcdf,   & ! do we want netcdf output
           disallowed    ! list of excluded fluxes
      
! by dfault, dump everything
      data dosources/.true./;  data dodata   /.true./
      data dogroups /.true./; data docovar  /.true./
      data donetcdf /.true./

      integer curind     ! current index for names etc
      character*80 line  ! read from block master file
      character*80 subnm ! subsidiary file name for covar blocks
      character*80 code  ! description of covar block
      integer i,j,k,n,i1,j1,ind ! index variables
      integer m,l,maxm,lyear
      integer nb         ! index variable for covar blocks
      real grouptot      ! summed source over group
      real groupvar      ! variance of grouptotal
      integer blkunit    ! unit number for block files
      parameter(blkunit=40)
      integer nblocks    ! number of covar blocks
      integer fnblnk     ! first nonblank function
      integer lnblnk     ! last nonblank function
      integer blksz      ! number of groups in block
      integer getblklst  ! returns names of groups in block
      character*80 blklst(maxperline)    ! list of group names in block
      integer blknum(maxperline)         ! numbers of group in block
      real blkcov(maxperline,maxperline) ! block covariance
      real getpack       ! function for getting value from packed array
      integer findname   ! finds name in list of names
! getting words from line
      character*80 words(maxwords), temp_name
      integer split      ! function to split line into words

! variables for netcdf
      integer nlats, nlons                     ! number of lats and longs
      real, dimension(:), allocatable :: lats, lons, times ! dimensions 
      real, dimension(:,:), allocatable :: map, dcsrc ! really dummy
      real, dimension(:,:,:), allocatable :: fluxes, basis_funcs
      character*80, dimension(:), allocatable :: basis_names
      include 'netcdf.inc'
      integer start(3),count(3),timeid,ncid,sfid,ok
      integer longitudedim, latitudedim, timedim
      integer dims(3),iret,longitudeid,latitudeid
      integer getmap     ! function for getting cif file
      integer nfluxes    ! total number of flux basis functions
      integer num_months ! number of data months
      integer num_stats  ! number of stations
      integer num_obs    ! number of station observations
      integer, dimension(:), allocatable :: datid ! data ids for netcdf

! read namelist
      open(outunit,file=outputf,status='old',err=99999)
      rewind outunit
      read(outunit,nml=outputnml,err=99998)
      close(outunit)

!     temp_name  = prefname(1:index(prefname, ' ')-1)//'.txt'
!     open(outunit,file=temp_name,   err=99997)
!     write(outunit,'(a80)') header
! first do sources,almost graphz format
      if(dosources) then
!       temp_name  = prefname(1:lnblnk(prefname))//'.sources'
        temp_name = prefname//'.sources'
        open(outunit,file=temp_name,   err=99997)
!       write(outunit,*) 'sources'
!       assume known ordering of sources in arrays
        kount=1
        kount2 = 0
        do l=1,(nlandreg+npresub)
          if (l.le.(nlandreg)) then
            maxm=monperyr
          else
            maxm=1
          endif
          if (l.eq.(nlandreg+npresub)) then
            lyear=1
          else
            lyear=lastsrc-firstsrc+1
          endif
          write(outunit,*,err=99994) sourcenm(kount)(1:60)
          curname = sourcenm(kount)
          write(outunit,*,err=99993) lyear*maxm,',',2
          do n=1,lyear
            do m=1,maxm
              k = kount2 + (m-1)*lyear + n
              if (sourcenm(sourceind(k)).ne.curname) stop 'check output source code'
              if (maxm.ne.1) then
                  write(outunit,999,err=99992) float(sourceyr(k))+ &
                  (float(sourcemnth(k))-0.5)/monperyr, source(k), sqrt(getpack(csrc,k,k))
 999      format(f9.4,2(1pe14.5))
              else
                  write(outunit,999,err=99992)float(sourceyr(k)), source(k), sqrt(getpack(csrc,k,k))
              endif
            enddo
          enddo
          kount  = kount + maxm
          kount2 = kount2 + maxm*lyear
        enddo

        close(outunit)
      endif

! now data
      if(dodata) then
!       temp_name  = prefname(1:lnblnk(prefname))//'.data'
        temp_name = prefname//'.data'
        open(outunit,file=temp_name,   err=99997)
!       write(outunit,*) 'data'
        curind = -1 ! impossible initial value
        do  n =1, ntotdat
          if(datind(n) .ne. curind) then ! onto new datum
            curind = datind(n)
            write(outunit,*,err=99991) datnm(curind)(1:60)
! now how many data of this type are there
            if(curind .eq. ndat) then ! last data type
              k = ntotdat - n + 1
            else
              do k = n,ntotdat
                if(datind(k) .ne. curind) goto 2000 ! break
              end do
 2000         k = k - n
            endif
            write(outunit,*,err=99993) k 
          endif
! construct decimal year and print out time, value and std-dev
          if (datav(datind(n)).eq.9999) then
            write(outunit,998,err=99990)float(datyr(n)) + &
            (float(datmnth(n))-0.5)/12.,dat(n), sqrt(cdat(n))
          else
            write(outunit,998,err=99990) float(datyr(n)) +    &
           (float(datmnth(n))-0.5)*float(datav(datind(n)))/   &
             float(pointsperyear), dat(n), sqrt(cdat(n))
 998      format(f10.5,2(1pe13.5))
          endif
        end do
        close(outunit)
      endif
      
      if(dogroups) then
!       temp_name  = prefname(1:lnblnk(prefname))//'.groups'
        temp_name = prefname//'.groups'
        open(outunit,file=temp_name,   err=99997)
        write(outunit,*) 'group av, std dev, sum, sum std dev, name'
        write(outunit,*) ngroups,',',4
        do n = 1,ngroups
          grouptot = 0.
          groupvar = 0.
          do i=1,grouplen(n) ! number of sources in group
            grouptot = grouptot + source(grouplst(i,n))
            do j = 1,grouplen(n)
              groupvar = groupvar +getpack(csrc,grouplst(i,n), grouplst(j,n))
            end do
          end do
          write(outunit,997,err=99989) n,grouptot/(grouplen(n)+1e-6),  & ! 1e-6 deals with empty groups 0/0
          sqrt(groupvar)/(grouplen(n)+1e-6),grouptot,sqrt(groupvar),groupnm(n)
 997      format(i5,4(1pe12.3),3x,a20)
        end do
        close(outunit)
      endif
      
! and now the really yucky bit where
! covariances of groups are written out
      if(docovar) then
        open(blkunit,file=blockf,status='old',err=99988)
        rewind blkunit
        read(blkunit,*,err=99987)
        read(blkunit,*,err=99986) nblocks
        
        do nb =1,nblocks

          read(blkunit,'(a80)',err=99985) line !  read whole thing now split
!         subnm first field, code (description) second field
          if(split(line, words, maxwords).lt.2) then
            write(stderr,*) 'output:format error line ',nb,' in ',blockf
            stop
          endif
          subnm = words(1)
          code = words(2)
          blksz = getblklst(subnm,blklst,maxperline)
          if(blksz .lt. 0) then
            write(stderr,*) 'error getting block list at data line ', nb,' in ',blockf
            stop
          endif
! now turn group names into numbers
          do j = 1,blksz
            blknum(j) = findname(groupnm, ngroups, blklst(j))
            if(blknum(j) .eq. 0) then
              write(stderr,*) 'error, ',blklst(j),' not found in ', 'group list, skipping whole block'
              write(outunit,*) 'error, ',blklst(j),' not found in ','group list, skipping whole block'
              goto 3000 ! really continue
            endif
            
          end do
            
! now finally in a position to do some output
          write(outunit,*,err=99983) code     ! explanatory header for block
          write(outunit,996,err=99983) (blklst(i),i=1,blksz) ! matrix header
 996      format(10x,7a10)
! now compute covariance matrix
! which is still the sum of the relevant submatrix of csrc
          do i=1,blksz
            do j=1,blksz
              blkcov(i,j) = 0.
              do i1 = 1,grouplen(blknum(i))
                do j1 = 1,grouplen(blknum(j))
                  blkcov(i,j) = blkcov(i,j) + getpack(csrc,grouplst(i1, blknum(i)), &
                            grouplst(j1,blknum(j)))
                end do
              end do
            end do
          end do
          do j=1, blksz
            write(outunit,995,err=99982) blklst(j), (blkcov(i,j),i=1,blksz)
 995        format(a10,7f10.3)
          end do
 3000     continue              ! target of break kludge
        end do
      endif
      close(outunit)
      

      if(donetcdf) then 
  ! this section takes a series of flux maps
  !corresponding to c fluxes from each region and for each year and month c
  !multiplies them by the relevant scaling factor and sums them c to
  !produce the total flux for each month c fluxes on the dissallowed list
  !are not included c the routine associates fluxes by name so that a
 !source called c ??ocean-A will match a flux map called ocean-A c this is
!to allow for the m-ocean-A, n-ocean-A etc names c Author: Peter Rayner c
!date: june 1999
        
! first find the basis functions, count them and allocate space
        open(outunit, file=basis_file, form='unformatted',status='old', err=99981)
        rewind outunit
        read(outunit) nlats
        allocate(lats(nlats))
        read(outunit) lats
        read(outunit) nlons
        allocate(lons(nlons))
        allocate(map(nlons, nlats))
        rewind outunit
                                ! now count maps
        nfluxes = 0
        do while(.true.)
          i = getmap(outunit, nlats, nlons, nlats, lats, nlons, lons,temp_name, map)
          if(i .eq. 1) then
            write(stderr,*) 'output: error getting flux map ',nfluxes+1, 'from ',basis_file
            stop
          else if(i .eq. -1)then ! EOF
            exit
          else
            nfluxes = nfluxes + 1
          endif
        end do
                                ! counted them, now allocate
        allocate(basis_funcs(nlons, nlats, nfluxes))
        allocate(basis_names(nfluxes))
                                ! and now read them in
        rewind outunit
        do i = 1, nfluxes
          j =  getmap(outunit, nlats, nlons, nlats, lats, nlons, lons, basis_names(i), basis_funcs(:,:,i))
        end do
                                ! ok, got all the information we need, create and fill flux series
        num_months = (lastdat - firstdat +1)*pointsperyear
        allocate(fluxes(nlons, nlats, num_months))
        fluxes = 0.
        do n = 1,nunks
          ind = sourceind(n)
          if((sourcetype(ind) <= 3) .and. ((sourceyr(n) .eq. everyyear).or.  &
               (sourceyr(n) .eq. linearsource) .or. (sourceyr(n) .ge. firstdat)) &
               .and. (findname(disallowed, max_disallowed, sourcenm(ind)).eq.0)) then ! carbon flux and not on the blacklist
            temp_name = sourcenm(ind)(3:)
            j = findname(basis_names, nfluxes, temp_name)
            if(j .eq. 0) then
              write(stderr,*) 'output:could not find source ',n,'named ',sourcenm(ind),' in basis function list'
            else
!              print*,n,j,sourcenm(ind),source(n)
              if(sourceyr(n) .eq. everyyear) then
                do i = firstdat, lastdat
                  k = (i - firstdat)* pointsperyear + sourcemnth(n)
                  fluxes(:,:,k) = fluxes(:,:,k) + source(n)*basis_funcs(:,:, j)
                end do
              else
                
                k = (sourceyr(n) - firstdat)* pointsperyear + sourcemnth(n)
                fluxes(:,:,k) = fluxes(:,:,k) + source(n)*basis_funcs(:,:, j)
              end if
              
            endif
            
          endif
        end do
                                ! now write netcdf output
        allocate(times(num_months))
        temp_name = prefname(1:index(prefname, ' ')-1)//'.cdf'
        ok = nf_create(temp_name,NF_NOCLOBBER,ncid)
        if (ok.ne.nf_noerr) write(6,*) 'create: ',nf_strerror(ok)

! define dimensions
        ok = nf_def_dim(ncid,'longitude',nlons,longitudedim)
        if (ok.ne.nf_noerr) write(6,*) 'defdim 1: ',nf_strerror(ok)
        ok = nf_def_dim(ncid,'latitude',nlats,latitudedim)
        if (ok.ne.nf_noerr) write(6,*) 'defdim 2: ',nf_strerror(ok)
        ok = nf_def_dim(ncid,'time',num_months,timedim)
        if (ok.ne.nf_noerr) write(6,*) 'defdim 3: ',nf_strerror(ok)
! define variables and attributes
        ok = nf_def_var(ncid,'longitude',NF_FLOAT,1,longitudedim,longitudeid)
        if (ok.ne.nf_noerr) write(6,*) 'defvar 1: ',nf_strerror(ok)
        ok = nf_put_att_text(ncid,longitudeid,'units',7,'degrees')
        if (ok.ne.nf_noerr) write(6,*) 'defatt 1: ',nf_strerror(ok)
        ok = nf_def_var(ncid,'latitude',NF_FLOAT,1,latitudedim,latitudeid)
        if (ok.ne.nf_noerr) write(6,*) 'defvar 2: ',nf_strerror(ok)
        ok = nf_put_att_text(ncid,latitudeid,'units',7,'degrees')
        if (ok.ne.nf_noerr) write(6,*) 'defatt 2: ',nf_strerror(ok)
        ok = nf_def_var(ncid,'time',NF_FLOAT,1,timedim,timeid)
        if (ok.ne.nf_noerr) write(6,*) 'defvar 3: ',nf_strerror(ok)
        ok = nf_put_att_text(ncid,timeid,'units',4,'year')
        if (ok.ne.nf_noerr) write(6,*) 'defatt 3: ',nf_strerror(ok)
        dims(1) = longitudedim
        dims(2) = latitudedim
        dims(3) = timedim
        ok = nf_def_var(ncid,'fluxes',NF_FLOAT,3,dims,sfid)
        if (ok.ne.nf_noerr) write(6,*) 'defvar 4: ',nf_strerror(ok)
        ok = nf_put_att_text(ncid,sfid,'units',8,'g/m^2/yr')
        if (ok.ne.nf_noerr) write(6,*) 'defatt 4: ',nf_strerror(ok)
        ! now count number of stations, unfortunately not
! ndat since that has constraints in it too
        do i =1,neqns
          if((datyr(i)<firstdat).or.(datyr(i)>lastdat)) exit
        end do
        num_obs = i -1
        num_stats = datind(num_obs)
! now define a variable for each station and constraint
! with time as its dimension
        allocate(datid(num_stats))
        dims(1) = timedim
        do i = 1, num_stats
          ok = nf_def_var(ncid,datnm(i),NF_FLOAT,1,timedim,datid(i))
         if (ok.ne.nf_noerr) write(6,*) 'defvar 5: ',nf_strerror(ok)
          ok = nf_put_att_text(ncid,datid(i),'units',4,'ppmv')
         if (ok.ne.nf_noerr) write(6,*) 'defatt 5: ',nf_strerror(ok)
        end do
        
!       call ncendf(ncid, iret) ! leaving define mode 
        ok = nf_enddef(ncid)
        if (ok.ne.nf_noerr) write(6,*) 'enddef: ',nf_strerror(ok)
! write longitude
        ok = nf_put_var_real(ncid,longitudeid,lons)
        if (ok.ne.nf_noerr) write(6,*) 'putvar 1: ',nf_strerror(ok)
! write latitude
        ok = nf_put_var_real(ncid,latitudeid,lats)
        if (ok.ne.nf_noerr) write(6,*) 'putvar 2: ',nf_strerror(ok)
        do k=1,num_months
          times(k) = float(firstdat) +(k -0.5)/pointsperyear
        enddo
        ok = nf_put_var_real(ncid,timeid,times)
        if (ok.ne.nf_noerr) write(6,*) 'putvar 3: ',nf_strerror(ok)
! now write fluxes
        start(1) = 1; start(2) = 1; start(3) = 1
        count(1) = nlons; count(2) = nlats; count(3) = num_months
        ok = nf_put_vara_real(ncid,sfid,start,count,fluxes)
        if (ok.ne.nf_noerr) write(6,*) 'putvar 4: ',nf_strerror(ok)
! now try writing data
        do i =1, num_obs
          ind = (datyr(i) -firstdat)*pointsperyear +datmnth(i)
           ok = nf_put_var1_real(ncid,datid(datind(i)),ind,dat(i))
        if (ok.ne.nf_noerr) write(6,*) 'putvar 5: ',nf_strerror(ok)
        end do
        
!       call ncclos (ncid, iret)
        ok = nf_close(ncid)
        if (ok.ne.nf_noerr) write(6,*) 'close: ',nf_strerror(ok)
        deallocate(lats,lons,map,basis_funcs,basis_names,fluxes,times,datid)
      endif
            
      return
99999 write(stderr,*) 'could not open name list file ',outputf,' in output'; stop
99998 write(stderr,*) 'error reading namelist from ',outputf,' in output'; stop
99997 write(stderr,*) 'could not open file ',prefname(1:index(prefname, ' ')-1)//'.text',' in output'; stop
99995 write(stderr,*) 'error writing header in output'; stop
99994 write(stderr,*) 'error writing source name ',curind,' in output'; stop
99993 write(stderr,*) 'error writing number of elements in output'; stop
99992 write(stderr,*) 'error writing source ',n,' in output'; stop
99991 write(stderr,*) 'error writing  data name ',curind,' in output'; stop
99990 write(stderr,*) 'error writing datum ',n,' in output'; stop
99989 write(stderr,*) 'error writing group total ',n,' in output'; stop
99988 write(stderr,*) 'error opening block file ',blockf; stop
99987 write(stderr,*) 'error reading header from ',blockf; stop
99986 write(stderr,*) 'error reading nblocks from ',blockf; stop
99985 write(stderr,*) 'error reading data line ',nb,' in ',blockf; stop
99984 write(stderr,*) 'error writing header for line ',n,' in ',blockf; stop
99983 write(stderr,*) 'eror writing matrix header for line ',n,' in ',blockf; stop
99982 write(stderr,*) 'error writing matrix for data line ',n,' in ',blockf; stop
99981 write(stderr,*) 'output: could not open basis file',basis_file; stop
      
      end
