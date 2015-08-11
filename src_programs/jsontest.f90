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


! first we read a Kossel namelist file
nmlfile = 'Kosseltest.nml'
call GetKosselNameList(nmlfile, knl)

! and then we write it to json format
jsonname = 'Kossel.json'
n_errors = 0
call JSONwriteKosselNameList(knl, jsonname, n_errors)

if (n_errors /= 0) stop 1

call Message('nml -> json conversion completed')

end program jsontest
!*****************************************************************************************

