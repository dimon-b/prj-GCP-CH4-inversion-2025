! This program performs a time-dependent inversion of trace gas
! records in terms of a series of pulsed sources with predetermined
! spatial structure.

! The input data is the Green's function response for each source, a 
! series of prior estimates for the time histories of each source (along
! with covariances) and a time history of various types of data.
! Addtional information is a set of constraints on groups of sources.

! The Green's functions are precalculated by the insertion of pulsed
! sources into an atmospheric transport model.

! The computational core of the code is the routine bayesinv.f which is
! a general-purpose Bayesian inverse solver.

! Author: Peter Rayner

! Date commenced: March 1996


! The implicit call on the next line is commented out since it crashes
! the compiler on SGI-R8000 machines as of 96/05/17. The code is very
! likely wrong but compiles without error on a R4400.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    program tdi

      implicit none

      include 'tdi.h'    ! parameters
      include 'common.h' ! common blocks 

! the name list specifies control files and some flags

      character*80 :: title,sourcef,datf, outputf,groupf,matdump,greenf,svdf, &
           greendump,srccorf, datcorf,recipef, constf,datflist

! dummy values for  name list

      data title   /'untitled'/;    data sourcef   /'dummysource'/
      data datf    /'dummydata'/;   data groupf    /'dummygroup'/
      data constf  /'dummyconst'/;  data outputf   /'dummyoutput'/
      data greenf  /'dummygreen'/ ; data srccorf   /'dummysrccor'/
      data datcorf /'dummydatcor'/; data recipef   /'dummyrecipe'/
      data matdump /'dummydump'/;   data greendump /'greendump'/
!     dumpmats=.false.;concheck=.false.;doinfluence=.false.;dodiff=.false. !pkp
!     doref=.false.;dosvd=.true.;useoldgreens=.false.;commentary=.false. !moved down

      integer numsvdvecs

      data numsvdvecs /-1/ ! that will mean dump everything
      data svdf       /'svdout'/ 

      real srcscale, datscale ! see namelist 

      data srcscale /1./ ! don't change it
      data datscale /1./ ! don't change it 

      integer firstsrc, lastsrc, firstdat, lastdat

! defaults for period of inversion, make it select everything

      data firstsrc/-10000/; data lastsrc /10000/ 
      data firstdat/-10000/; data lastdat /10000/
      integer nlandreg, npresub
      namelist     /control/  &
           title,       & ! title for output files
           sourcef,     & ! name of source file
           datf,        & ! name of data control file
           outputf,     & ! name of output control file
           groupf,      & ! name  of group file
           constf,      & ! name of constraint file
           greenf,      & ! name of green's function master file
           recipef,     & ! name of src-data coefficient recipe file
           srccorf,     & ! name of source correlation file
!          datcorf,       ! name of data correlation file
           dumpmats,    & ! flag to say whether to dump matrices
!          do_extend,     ! use old solution to extend data
           matdump,     & ! name of output file
           dumpgreens,  & ! flat to say wehterh to dump green's funcs
!          greendump,     ! name of green's func dump file
           concheck,    & ! flag to check condition number in inversion
           doinfluence, & ! flag to  do jacobian influence check
           dodiff,      & ! flag to do differential rather than integral form
           doref,       & ! flag to do trend + difference calc
           dosvd,       & ! flag to do an svd of frechet
           numsvdvecs,  & ! number of singular vectors
!          svdf,          ! file to dump svd
!          useoldgreens,  ! grab old GF's from previous greendump
           dumpgroups,  & ! flag to dump group lists
!          commentary,    ! be chatty about progress
           srcscale,    & ! factor to scale source variance for tuning chi-sq
           datscale,    & ! factor to scale data variance for tuning chi-sq
           firstsrc,    & ! year of first source
           lastsrc,     & ! year of last source
           firstdat,    & ! year of first datum
           lastdat,     & ! year of last datum
!     rml 29/9/00
           tc3steady    & ! flag for tc3 style input
          ,pseudodata   & ! flag for pseudodata (2190 timestep, 228 sites)
          ,nlandreg,npresub  & ! no. of land,ocean,presub   
          ,datflist     &
!     rml 5/3/01
          ,tc3tdi

! dimensions of everything
      integer nsources   ! number of source components
      integer nunks      ! n-unknowns nsources*time-history-length
      integer nobs       ! number of sets of observations
      integer ndat       ! number of sets of data and constraints
      integer ntotdat    ! total data elements
      integer nconsts    ! number of constraint equations
      integer ngroups    ! number of groups
      integer neqns      ! total number of equations for dimensioning
      integer ngreens    ! number of distinct GF's
      integer lgreen     ! length of GF record
! rml 29/9/00


      integer nobsmax   ! number of sites in input file

!  flag
      logical shipdata

! some functions
      integer getnum       ! function for counting objects
      integer getlgreen    ! function regurning length of GF record
      integer defgrps      ! both counts and assigns default groups
      integer noffsets     ! number of global offsets
      integer getoffsets   ! function to count offsets

! some local variables, surprisingly few
      integer status       ! returned from some subroutines
      integer i,j          ! needed index variable
      integer nvecs        ! number of SV's
      integer kount        ! counter for number of sites
      integer nread        ! count actually read from F7 NetCDF
      character(len=80) :: arg1
      integer index, itime ! temporary variables
      integer nyear        ! no. of years in tdi, set to 1 for csi
      real                                      :: templat,templon
      character                                 :: tempname*80

! md covariance, we need to allocate it here
      real, dimension(:), allocatable :: var  ! model covariance *and* LSV's
      integer dummy                           ! useless argument for some routines
      real, dimension(:), allocatable :: tmpmat1, tmpmat2,scrdg ! one row of a matrix
      real, dimension(:,:), pointer   :: dvar ! data variance

      include 'pointer.h'

      nsources = 0; nunks = 0; nobs = 0; ndat = 0; ntotdat= 0
      nconsts = 0; ngroups = 0; neqns = 0; nobsmax = 0
!     dumpmats=.false.;concheck=.false.;doinfluence=.false. !;dodiff=.false. !pkp
!     doref=.false.;dosvd=.true.;useoldgreens=.false.;commentary=.false. !data satement

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (command_argument_count().gt.0) then
          call get_command_argument(1,arg1)
          if (trim(arg1).eq.'--valid-months-version') then
             write(*,'(a)') 'TDI_VALID_MONTHS_V1'
             stop
          endif
      endif
      Open(mstunit,file='control')
      Rewind mstunit
      Read(mstunit,nml=control)
      Close(mstunit)

      Write(stdout,nml=control)

! first get numbers of data, sources and constraints for allocating
! arrays number of types of sources directly, it's so simple
      ntotdat = 0

      If (tc3steady) then
! sources, currently 268
         nsources = nlandreg*monperyr + npresub
         nunks    = nsources; nyear = 1
! data
         Open(inunit,file=datflist,status='old',err=99988)
         Rewind inunit; Read(inunit,*); Read(inunit,*,err=99998) nobsmax
         kount = 0
         Do i = 1, nobsmax            !in our list file, itime = 9999 or 0 only
            Read(inunit,*) index,itime  !if itime=0 then, it doesnt add to nobs
            If (itime == 9999) then; kount = kount + 1  ! kount becomes nobs..
               ntotdat = ntotdat + monperyr
            Elseif (itime /= 0) then;kount = kount + 1  !we dont have such cases - pkp
               ntotdat = ntotdat + pointsperyear/itime
            End if
         End do; Close(inunit)
         nobs = kount; Write(6,*) 'nobs nyear = ',nobs, nyear
!        ntotdat = nobs*pointsperyear  ! true while itime=1
      Elseif (tc3tdi) then
         nsources = (nlandreg)*monperyr + npresub
! offset isn't multiplied by the number of years
         nunks = (nsources-1)*(lastsrc-firstsrc+1) + 1; nyear=(lastsrc-firstsrc+1)
         write(6,*) 'nunks= ',nunks, nyear
! data
         open(inunit, file=datflist,status='old',err=99988)
         rewind inunit; read(inunit,*); read(inunit,*,err=99998) nobsmax
         kount = 0
         do i = 1, nobsmax
            read(inunit,*) index,itime
            if (itime == 9999) then; kount = kount + 1
               ! Count valid months after the active names are known.
            elseif (itime.ne.0) then                    !we dont have such cases - pkp
               stop 'Only monthly data allowed with tc3tdi'
            endif !; print *, ' i, ntotdat ', i, kount, nobsmax, ntotdat
         enddo; close(inunit)
          nobs = kount
          call count_globalview(datf,datflist,nobs,monperyr, &
               firstdat,lastdat,ntotdat)
          if (ntotdat.le.0) stop 'TDI: no valid observations'
          write(6,*) 'nobs, valid ntotdat, excluded = ', &
               nobs,ntotdat,nobs*monperyr*(lastdat-firstdat+1)-ntotdat
      Else
         open(inunit, file=sourcef, status='old',err=99999)
         rewind inunit; read(inunit,*) 
         read(inunit,*,err=99998) nsources; close(inunit)
! now data, for both tc3steady and tc3tdi
         open(inunit,file=datf,status='old',err=99996)
         rewind inunit; read(inunit,*) ! skip header
         read(inunit,*,err=99995) nobs; close(inunit)
      Endif

! now number of constraint equations

      open(inunit,file=constf,status='old',err=99994)
      rewind inunit; read(inunit,*) ! skip header
      read(inunit,*,err=99993) nconsts; close(inunit)
      if (nconsts < 0) then
         write(stderr,*) 'error: nconsts= ',nconsts,' in file ',constf
         stop
      endif

      if (.not.(tc3steady.or.tc3tdi)) then
! now get other numbers from recursive reading via master files
! yucky enough to use routine
         nunks = getnum(sourcef, firstsrc, lastsrc)
         if (nunks .le. 0) then
            Write(stderr,*) 'error, nsrouces=',nunks
            Write(stderr,*) ' from file ',sourcef
            stop
         endif
         ntotdat = getnum(datf, firstdat, lastdat)
         if (ntotdat .lt. 0) then
            write(stderr,*) 'error: ntotdat=',ntotdat,' in file ',datf
            stop
         endif
      endif

! now  global offsets, one for each type of data, find how many?      
!      noffsets = getoffsets(datf)

! must postpone numbers of groups until we can compute
! default groups from source names etc
! but at least read explicit groups at this stage

      open(inunit,file=groupf,status='old',err=99992)
      rewind inunit; read(inunit,*)    ! skip header
      read(inunit,*,err=99991) ngroups ! will increase this later
      close(inunit)
      if (ngroups < 0) then
         write(stderr,*) 'error: ngroups= ',ngroups,' in file ',groupf
         stop
      endif

! constraints are really extra data make space for them
      neqns = ntotdat + nconsts    ! total number of equations
      ndat = nobs+nconsts
      print*,'ndat,ntotdat,nconsts,neqns', ndat,ntotdat,nconsts,neqns

! now about Green's functions
      lgreen = getlgreen(tc3tdi) ! almost certainly dummy function
      Write(*,*) 'lgreen',lgreen
      if (tc3steady.or.tc3tdi) then
         ngreens = nsources
      else
         open(inunit, file=greenf, status='old', err=99990)
         rewind inunit; read(inunit,*) ! skip header record
         read(inunit,*,err=99989) ngreens ! total number of distinct green's funcs
         close(inunit)  ! that's all we need at the moment
      endif

! now do space allocation, totally SGI specific
! can't do groups yet since don't know how many until have sources
! commented out in f90 version with changing allocation rules

      allocate(src(nunks));         src = 0.0
      allocate(srcnm(nsources))
      allocate(srcind(nunks));      srcind = 0
      allocate(greenind(nsources)); greenind = 0
      allocate(srctype(nsources));  srctype = 0
      allocate(csrc(nunks));        csrc = 0.0
      allocate(srcyr(nunks));       srcyr = 0
      allocate(srcmnth(nunks));     srcmnth = 0
      allocate(dat(neqns));         dat = 0.0
      allocate(datcopy(neqns));     datcopy = 0.0
      allocate(cdat(neqns));        cdat = 0.0
      allocate(siteindex(nobs));    siteindex=0
      allocate(datnm(ndat))
      do i=1,ndat
        write(datnm(i),'(80x)')
      enddo
      allocate(datlat(ndat));       datlat = 0.0
      allocate(datlon(ndat));       datlon = 0.0
      allocate(dattype(ndat));      dattype = 0
      allocate(datav(ndat));        datav = 0
      allocate(datyr(neqns));       datyr = 0
      allocate(datmnth(neqns));     datmnth = 0
      allocate(datind(neqns));      datind = 0
      allocate(frech(neqns,nunks)); frech = 0.0

      if(nconsts.gt. 0) then
        allocate(constnm(nconsts))
      endif

! note the extra space for global offsets in source arrays
! now get source information
      if (tc3steady) then
        call read_sourcep  (sourcef,monperyr,nlandreg,npresub, &
             nsources,nunks,srcnm,greenind,srctype,src,srcind,srcyr,srcmnth,csrc)
      elseif (tc3tdi) then
        call read_sourcetdi(sourcef,monperyr,nlandreg,npresub, &
             nsources,nunks,srcnm,greenind,srctype,src,srcind,srcyr,srcmnth,csrc,firstsrc,lastsrc)
      else
        call getsrc(sourcef,nsources,nunks,ngreens,srcnm,greenind,srctype,src, &
             srcind, srcyr,srcmnth, csrc, firstsrc, lastsrc)
      endif


      If (tc3steady.or.tc3tdi) then
         Open(inunit,file=datflist,status='old',err=99988)
         Rewind inunit; Read(inunit,*); Read(inunit,*,err=99998) nobsmax
         kount = 0; Write(*,*) 'Sites used'
         Do i = 1, nobsmax
            Read(inunit, *) index,itime,tempname,templat,templon  ! free format pkp / 20050113
 8001       Format(2i6,3x,a25,f7.2,f9.2)
            If (itime /= 0) then
               kount = kount + 1; siteindex(kount) = index
               If (.not.pseudodata.and.itime < 9999) then
                  Write(*,*) 'inversion from observations must use '
                  Write(*,*) 'monthly data: itime=9999'
                  stop
               Endif
               datav(kount) = itime; datnm(kount)(1:45) = tempname
               Write(6,*) siteindex(kount),datav(kount),tempname  !pkp
               datlat(kount) = templat; datlon(kount) = templon
               dattype(kount) = 1
            Endif
         Enddo; Close(inunit)
         If (tc3steady) then
            if (pseudodata) then
               stop 'change pointsperyear in tdi.h'
               call read_pseudodata(datf,nobsmax,nobs,neqns,pointsperyear, &
                    siteindex,dat,datind,datyr,datmnth,cdat,ndat,datav)
            else
               call read_monthlyobs(datf,nobs,neqns,monperyr,datnm,dat,    &
                    datind,datyr,datmnth,cdat,ndat,datav)
            Endif
         Else
! tdi case
            if (pseudodata) then
               stop 'pseudodata not set up for tdi'
            Else
               call read_globalview(datf,nobs,neqns,monperyr,datnm,dat,   &
                     datind,datyr,datmnth,cdat,ndat,datav, &
                     firstdat,lastdat,nread)
                if (nread.ne.ntotdat) then
                   write(stderr,*) 'TDI count/read mismatch:', &
                        ntotdat,nread
                   stop 2
                endif
            Endif
         Endif
      Else
        call getdat(datf,nobs,neqns,datnm,datlat,datlon,dattype,dat,      &
                    datind,datyr,datmnth,cdat,firstdat,lastdat)
      Endif


! now find how many groups
! we already have explicit ones, now count default groups

      ngroups = ngroups + defgrps(ngroups, nunks, nsources, srcnm, srcyr,   &
        srcmnth, srctype, srcind,dummy, dummy, dummy, firstdat, lastdat, .true.)

! the last .true. means here that we're only counting
! this is imperative since the grouplst hasn't been allocated yet
      allocate(groupnm(ngroups))
      allocate(grouplst(nunks, ngroups)); grouplst = 0
      allocate(grouplen(ngroups));        grouplen = 0

      call getgrps(groupf, ngroups, nunks, nsources, srcnm, srcyr, srcmnth, &
           srctype, srcind, groupnm, grouplen, grouplst, firstdat, lastdat)
      write(stdout,*) 'nunks=',nunks,' neqns=',neqns, ' ngroups=',ngroups,  &
           'max size=', (nunks+neqns+ngroups)*nunks/1e6

! build GF's
! we want to allocate them here so when  we throw them away
! they don't just leave a hole

!-------------------------start of the problem here--------------------------------
      allocate(greens(nyear, lgreen, nobs, ngreens)); greens = 0.0; print*,'Greens allocated '

      if (tc3steady .or. tc3tdi) then !; nyear = 1
        call tc3greens(greenf,greens,nyear,lgreen,nobs,ngreens,nlandreg,  &
          npresub,siteindex,ndat,datav,greendump,tc3tdi)
      else
        call bldGreens(greenf,greens,nyear,lgreen,nobs,ngreens,datnm,greendump) !call for the monthly mean version
      endif
      
!     can use standard bldfrech as long as pointsperyear=12 in tdi.h
      call bldfrech(recipef, nunks, neqns, nsources, ngreens, lgreen, &
           nobs, ntotdat, frech, srcyr, srcmnth, srcind, srctype,     &
           greenind, datyr, datmnth, datind, dattype, greens, firstsrc, nyear)

!-------------------------end of the problem here--------------------------------

! can now throw away GFs
      deallocate(greens)
      call getconsts( constf, neqns, nunks, ntotdat, nobs,           &
           frech, dat, cdat, datyr, datmnth, datind, datnm,          &
           ngroups, groupnm, grouplen, grouplst, firstdat, lastdat)

! now tune variances if necessary
      call tunevar(nunks, csrc, srcscale, neqns, cdat, datscale)

! now correlations for both source and data
      call srccor(srccorf, nunks, ngroups, grouplst, grouplen, groupnm, srcyr, csrc)
      call datcor( datcorf, neqns, ntotdat, cdat)

      if(doinfluence) call influence(nunks, neqns, nsources, ntotdat, ndat, &
           srcind, srcyr, srcmnth, srcnm, dat, datind, datyr, datmnth, datnm, cdat, frech)

! because of the lack of garbage collection in f90, we want to use
! the same space for the source covariance and the left singular 
! vectors in the SVD which follows.
! we-d like to allocate and deallocate as we need but
! it turns out we're frequently keeping two copies that way,
! somewhat unpredictably.  So, we'll allocate it here,
! make sure it's big enough for both jobs and use pointers
!  within svdcalc.f
! first find out how big it needs to be
! the covariance is a packed symmetric array order nunks
! the LSV's have length neqns and it depends on how many SV's there are
      nvecs = min(neqns, nunks) ! number of s-vecs
! now make sure it's big enough for either job
      allocate(var(max((nunks*(nunks+1))/2,neqns*nvecs)))
      var(:) = 0.0
! now populate with initial covariance
      do i=1,nunks; call putpack(var, i, i, csrc(i)); end do

      if(do_extend) then ! use old sources to build extended data rcrd
        open(sbunit,file=matdump,form='unformatted')
        rewind sbunit
        call getmats(sbunit, nunks, nsources, src, srcnm, srcyr, srcmnth, srcind, var)
        close( sbunit)
! now calculate data based on Jacobian and these sources
        call dgemv('n', neqns, nunks, 1.0, frech, neqns, src, 1, 0.0, dat, 1)
! now do data covariance, we do it one element at a time
! to avoid storing the whole data covariance matrix
! we only want the diagonal anyway
! the expression is frech * csrc * frech^t
! we calculate each row of frech* csrc individually
! then multiply it by a column of frech^t
! first allocate space for single row of frech * csrc
        allocate(tmpmat1(nunks)); allocate(tmpmat2(nunks))
        do i = 1,neqns
          call sspmv('u', nunks, 1., var, frech(i,1), neqns, 0.0, tmpmat1, 1)
          cdat(i) = dot_product(frech(i,:), tmpmat1)
        end do
!  now output it
        call output( outputf, nunks, neqns, ngroups, nsources, neqns, ndat,      &
             src, srctype, srcind, srcyr, srcmnth, srcnm, var, dat, datind, datyr, &
             datmnth, datnm, cdat, grouplst, groupnm, grouplen, firstdat, lastdat, &
             'extended ', title)                  ! space at end of name mandatory
        stop ' just data extension' 
      endif
      
      if (tc3steady) then
        call output_tc3( outputf, nunks, neqns, ngroups, nsources, neqns, ndat,  &
           src, srctype, srcind, srcyr, srcmnth, srcnm, var, dat, datind, datyr, &
           datmnth, datav, datnm, cdat, grouplst, groupnm, grouplen, firstdat,   &
           lastdat, 'prior', title) ! space at end of name mandatory
      elseif (tc3tdi) then
        call output_tc3tdi( outputf, nunks, neqns, ngroups, nsources, neqns, ndat, &
           nlandreg,npresub,firstsrc,lastsrc, src, srctype, srcind, srcyr, &
           srcmnth, srcnm, var, dat, datind, datyr, datmnth, datav, datnm, cdat,   &
           grouplst, groupnm, grouplen, firstdat, lastdat, 'prior', title) ! space at end of name mandatory
      else
      call output( outputf, nunks, neqns, ngroups, nsources, neqns, ndat, src,   &
           srctype, srcind, srcyr, srcmnth, srcnm, var, dat, datind, datyr,      &
           datmnth, datnm, cdat, grouplst, groupnm, grouplen, firstdat, lastdat, &
           'prior ', title) ! space at end of name mandatory
!      stop 'just want prior data fit'
      endif

      if(dosvd) then
        call svdcalc(svdf, nunks, neqns, neqns, frech, neqns, cdat, csrc,      &
             numsvdvecs, nsources, ntotdat, ndat, src, srcind, srcyr, srcmnth, & 
             srcnm, dat, datind, datyr, datmnth, datnm, var)
      else

! create and populate data covar as well
        allocate(dvar(neqns,neqns)); dvar = 0.0
        do i=1,neqns; dvar(i,i) =  cdat(i); end do
        
        call bayesinv( nunks, neqns, neqns, frech, src, dat, nunks, var, neqns, dvar, status)
! now copy data variance back for output
        do i=1,neqns; cdat(i) = dvar(i,i); end do
      endif
      
      if(dumpmats) then
        open(sbunit,file=matdump,form='unformatted'); rewind sbunit
        call putmats(sbunit, nunks, nsources, src, srcnm, srcyr, srcmnth, srcind, var)
        close( sbunit)
      endif
      
      if (tc3steady) then
        call output_tc3( outputf, nunks, neqns, ngroups, nsources, neqns, ndat, &
           src, srctype, srcind, srcyr, srcmnth, srcnm, var, dat, datind, datyr, &
           datmnth, datav, datnm, cdat, grouplst, groupnm, grouplen, firstdat, &
           lastdat, 'predicted', title) ! space at end of name mandatory
      elseif (tc3tdi) then
        call output_tc3tdi( outputf, nunks, neqns, ngroups, nsources, neqns, ndat, &
           nlandreg,npresub,firstsrc,lastsrc, src, srctype, srcind, srcyr, &
           srcmnth, srcnm, var, dat, datind, datyr, datmnth, datav, datnm, cdat,   &
           grouplst, groupnm, grouplen, firstdat, lastdat, 'predicted', title) ! space at end of name mandatory
      else
      call output(outputf, nunks, neqns, ngroups, nsources, neqns, ndat, src,  &
           srctype, srcind, srcyr, srcmnth, srcnm, var, dat, datind, datyr,    &
           datmnth, datnm, cdat, grouplst, groupnm, grouplen, firstdat, lastdat, &
           'predicted ', title) ! space at end of name mandatory
      endif

      stop 
99999 write(stderr,*) 'could not open file',sourcef; stop 
99998 write(stderr,*) 'error: nsources in file',sourcef; stop 
99997 write(stderr,*) 'cannot open file ',matdump; stop  
99996 write(stderr,*) 'cannot open file ',datf; stop  
99995 write(stderr,*) 'error: ndat from file ',datf; stop  
99994 write(stderr,*) 'error opening ',constf; stop  
99993 write(stderr,*) 'error: nconsts in file ',constf; stop  
99992 write(stderr,*) 'error opening ',groupf; stop  
99991 write(stderr,*) 'error reading ngroups in ',groupf; stop  
99990 write(stderr,*) 'error opening GF file ',greenf; stop  
99989 write(stderr,*) 'could not read ngreens ',greenf; stop  
99988 write(stderr,*) 'cannot open file ',datflist
      stop; end
