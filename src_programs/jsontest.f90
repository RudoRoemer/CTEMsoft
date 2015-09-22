!*****************************************************************************************
program jsontest

!! this is a test to convert a small nml file into a json file
    
use local
use io
use NameListHandlers
use JSONsupport

IMPLICIT NONE

integer(kind=irg)       :: n_errors
character(fnlen)        :: nmlfile, jsonname
type(KosselNameListType):: knl
type(MCCLNameListType)  :: mcnl


! read a MCOpenCL input file
nmlfile = 'EMMCOpenCL.nml'
call GetMCCLNameList(nmlfile, mcnl)

! and write it to json format
nmlfile = 'EMMCOpenCL.nml'
jsonname = 'EMMCOpenCL.json'
call JSONwriteMCCLNameList(mcnl, jsonname, n_errors)


stop


! first we read a Kossel namelist file
!nmlfile = 'Kosseltest.nml'
!call GetKosselNameList(nmlfile, knl)

! and then we write it to json format
!jsonname = 'Kossel.json'
!n_errors = 0
!call JSONwriteKosselNameList(knl, jsonname, n_errors)

!if (n_errors /= 0) stop 1


! next, we need to do the inverse (json->nml) and compare the two nml files
!n_errors = 0
!!call JSONreadKosselNameList(knl, jsonname, n_errors)

!write (*,*) knl%stdout 
!write (*,*) knl%numthick
!write (*,*) knl%npix
!write (*,*) knl%maxHOLZ
!write (*,*) knl%nthreads
!write (*,*) knl%k
!write (*,*) knl%fn
!write (*,*) knl%voltage
!write (*,*) knl%dmin
!write (*,*) knl%convergence
!write (*,*) knl%startthick
!write (*,*) knl%thickinc
!write (*,*) knl%minten
!write (*,*) knl%xtalname
!write (*,*) knl%outname
!
!write (*,*) 'number of errors = ',n_errors

!if (n_errors /= 0) stop 1

end program jsontest
!*****************************************************************************************

