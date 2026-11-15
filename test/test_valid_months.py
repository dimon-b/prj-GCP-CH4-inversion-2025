# -*- coding: utf-8 -*-
"""Compile and execute the compacted F7 reader with a synthetic NetCDF-API stub.

No real inversion or scientific data files are touched.
"""
import shutil
import subprocess
import tempfile
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / 'inversion' / 'read_globalview.f90'

STUB_INC = '''
integer, parameter :: nf_noerr = 0
integer :: nf_open, nf_inq_dimid, nf_inq_dimlen, nf_inq_varid
integer :: nf_inq_vardimid, nf_get_var_text, nf_get_var_real, nf_close
character(len=80) :: nf_strerror
'''

STUB = r'''
integer function nf_open(filename,mode,ncid)
  implicit none
  character(len=*) filename
  integer mode,ncid
  ncid=1
  nf_open=0
end function
integer function nf_inq_dimid(ncid,name,dimid)
  implicit none
  integer ncid,dimid
  character(len=*) name
  select case (trim(name))
  case('station_count'); dimid=1
  case('tlen'); dimid=2
  case('data_type'); dimid=3
  case('namelen'); dimid=4
  case default; error stop 'Unknown mock dimension'
  end select
  nf_inq_dimid=0
end function
integer function nf_inq_dimlen(ncid,dimid,n)
  implicit none
  integer ncid,dimid,n
  select case (dimid)
  case(1); n=3
  case(2); n=12
  case(3); n=2
  case(4); n=45
  case default; error stop 'Unknown mock dimension length'
  end select
  nf_inq_dimlen=0
end function
integer function nf_inq_varid(ncid,name,varid)
  implicit none
  integer ncid,varid
  character(len=*) name
  select case(trim(name))
  case('tracer_data'); varid=10
  case('site_name'); varid=11
  case('tlen'); varid=12
  case default; error stop 'Unknown mock variable'
  end select
  nf_inq_varid=0
end function
integer function nf_inq_vardimid(ncid,varid,dims)
  implicit none
  integer ncid,varid,dims(*)
  dims(1)=3
  dims(2)=1
  dims(3)=2
  nf_inq_vardimid=0
end function
integer function nf_get_var_text(ncid,varid,names)
  implicit none
  integer ncid,varid
  character(len=45) names(*)
  names(1)='ALT'
  names(2)='ALT' ! Last-name-match should select this row
  names(3)='MLO'
  nf_get_var_text=0
end function
integer function nf_get_var_real(ncid,varid,arr)
  implicit none
  integer ncid,varid,i,m,idx
  real arr(*)
  if (varid==12) then
      do m=1,12
          arr(m)=1999.0+(real(m)-0.5)/12.0
      enddo
  else
      do m=1,12
          do i=1,3
              idx = 2*((m-1)*3+(i-1))+1
              arr(idx) = real(i*100 + m)
              arr(idx+1) = 2.0
              if (i==3) arr(idx+1) = 3.0
          enddo
      enddo
      arr(2*((2-1)*3+(2-1))+1) = -999.0 ! ALT Feb invalid
      arr(2*((3-1)*3+(2-1))+2) = 0.0    ! ALT Mar invalid
      arr(2*((12-1)*3+(3-1))+1) = -999.0 ! MLO Dec invalid
  endif
  nf_get_var_real=0
end function
integer function nf_close(ncid)
  implicit none
  integer ncid
  nf_close=0
end function
character(len=80) function nf_strerror(status)
  implicit none
  integer status
  nf_strerror='synthetic error'
end function
'''

DRIVER = r'''
program test_reader
    implicit none
    character(len=80) :: names(2), mockf, listf
    integer :: av(2),nvalid,nread
    integer, allocatable :: ind(:),yrs(:),mons(:)
    real, allocatable :: dat(:),cdat(:)

    open(11,file='datflist',status='replace')
    write(11,'(A)') 'Station list'
    write(11,'(I0)') 3
    write(11,'(A)') '1 9999 ALT 82 -62'
    write(11,'(A)') '2 0 OTHER 0 0'
    write(11,'(A)') '3 9999 MLO 19 -155'
    close(11)

    mockf='mock.nc'; listf='datflist'
    call count_globalview(mockf,listf,2,12,1999,1999,nvalid)
    if (nvalid /= 21) error stop 'Incorrect count: expected 21'
    allocate(ind(nvalid+1),yrs(nvalid+1),mons(nvalid+1))
    allocate(dat(nvalid+1),cdat(nvalid+1))
    ind=-1; yrs=-1; mons=-1; dat=-1; cdat=-1
    names(1)='ALT'; names(2)='MLO'
    av=9999
    call read_globalview(mockf,2,nvalid+1,12,names,dat,ind, &
                         yrs,mons,cdat,2,av,1999,1999,nread)
    if (nread /= nvalid) error stop 'Count differs from read'
    if (ind(1) /= 1 .or. mons(1) /= 1) error stop 'Bad first row'
    if (abs(dat(1)-201.0)>1.e-6) error stop 'Not last duplicate station'
    if (mons(2) /= 4 .or. dat(2) /= 204.0) error stop 'Invalid months not skipped'
    if (ind(11) /= 2 .or. mons(11) /= 1) error stop 'Bad second site index'
    if (ind(21) /= 2 .or. mons(21) /= 11) error stop 'Bad final month'
    if (cdat(1) /= 4.0 .or. cdat(11) /= 9.0) error stop 'Covariance not squared'
    if (ind(22) /= -1) error stop 'Unexpected extra row'
    print '(A)', 'F7 compacted-reader runtime test: PASS (21/24 months)'
end program
'''


def run_tests():
    compiler = shutil.which('gfortran')
    if not compiler:
        print('F7 reader runtime test: SKIP (gfortran unavailable)')
        return
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        (directory / 'netcdf.inc').write_text(STUB_INC)
        (directory / 'mock_netcdf.f90').write_text(STUB)
        (directory / 'driver.f90').write_text(DRIVER)
        cmd = [compiler, '-I', str(directory), '-fcheck=all', str(SRC),
               str(directory / 'mock_netcdf.f90'), str(directory / 'driver.f90'),
               '-o', str(directory / 'test_reader')]
        subprocess.run(cmd, check=True, cwd=directory)
        result = subprocess.run([str(directory / 'test_reader')], capture_output=True,
                                text=True, check=False, cwd=directory)
        if result.returncode: raise AssertionError('STDOUT:'+result.stdout+'\nSTDERR:'+result.stderr)
        print(result.stdout.strip())


if __name__ == '__main__':
    run_tests()
