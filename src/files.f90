! ###################################################################
! Copyright (c) 2014, Marc De Graef/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list 
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this 
!        list of conditions and the following disclaimer in the documentation and/or 
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
!        of its contributors may be used to endorse or promote products derived from 
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

!--------------------------------------------------------------------------
! EMsoft:files.f90
!--------------------------------------------------------------------------
!
! MODULE: files
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief everything that has to do with file-based input-output
! 
!> @version
!
!> @date    1/ 5/99 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   01/10/14 MDG 4.0 update after new cell type
!> @date   03/29/15 MDG 5.0 reformatted xtal files in HDF5 format; removed obsolete routines
!> @date   05/05/15 MDG 5.1 removed all getenv() calls; replaced with global path strings
!--------------------------------------------------------------------------

module files

use local
use typedefs

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: DumpXtalInfo
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Write a brief summary of the crystal structure on the screen
! 
!> @date    1/ 5/99 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   01/10/14 MDG 4.0 update after new cell type
!> @date   06/06/14 MDG 4.1 added cell pointer as argument, corrected Message routine
!--------------------------------------------------------------------------
subroutine DumpXtalInfo(cell, stdout)    

use constants
use io
use symmetry

IMPLICIT NONE

type(unitcell), pointer                 :: cell
integer(kind=irg),INTENT(IN),OPTIONAL   :: stdout

integer(kind=irg)                       :: i, j, oi_int(3), std
real(kind=dbl)                          :: oi_real(5)

 std = 6
 if (PRESENT(stdout)) std = stdout

 call Message('', frm = "(A/)", stdout = std)
 call Message('Crystal Structure Information', frm = "('-->',A,'<--')", stdout = std)
 oi_real(1) = cell%a
 call WriteValue('  a [nm]             : ', oi_real, 1, "(F9.5)", stdout = std)
 oi_real(1) = cell%b
 call WriteValue('  b [nm]             : ', oi_real, 1, "(F9.5)", stdout = std)
 oi_real(1) = cell%c
 call WriteValue('  c [nm]             : ', oi_real, 1, "(F9.5)", stdout = std)
 oi_real(1) = cell%alpha
 call WriteValue('  alpha [deg]        : ', oi_real, 1, "(F9.5)", stdout = std)
 oi_real(1) = cell%beta
 call WriteValue('  beta  [deg]        : ', oi_real, 1, "(F9.5)", stdout = std)
 oi_real(1) = cell%gamma
 call WriteValue('  gamma [deg]        : ', oi_real, 1, "(F9.5)", stdout = std)
 oi_real(1) = cell%vol
 call WriteValue('  Volume [nm^3]      : ', oi_real, 1, "(F12.8)", stdout = std)
 oi_int(1) = cell%SYM_SGnum
 call WriteValue('  Space group #      : ', oi_int, 1, "(1x,I3)", stdout = std)
 call WriteValue('  Space group symbol : ', trim(SYM_SGname(cell%SYM_SGnum)) , stdout = std)
 call WriteValue('  Generator String   : ',  trim(SYM_GL(cell%SYM_SGnum)) , stdout = std)
 if ((cell%SYM_SGset.eq.2).AND.(cell%xtal_system.ne.5)) then 
  call Message('   Using second origin setting', frm = "(A)", stdout = std)
 endif
 if ((cell%SYM_SGset.eq.2).AND.(cell%xtal_system.eq.5)) then 
  call Message('   Using rhombohedral parameters', frm = "(A)", stdout = std)
 endif
  if (cell%SG%SYM_centrosym) then 
    call Message('   Structure is centrosymmetric', frm = "(A)", stdout = std)
 else 
   call Message('   Structure is non-centrosymmetric', frm = "(A)", stdout = std)
 end if
! generate atom positions and dump output  
 call Message('', frm = "(A/)", stdout = std)
 call CalcPositions(cell,'v')
 oi_int(1) = cell%ATOM_ntype
 call WriteValue('  Number of asymmetric atom positions ', oi_int, 1, stdout = std)
 do i=1,cell%ATOM_ntype
  oi_int(1:3) = (/i, cell%ATOM_type(i), cell%numat(i)/)
  call WriteValue('  General position / atomic number / multiplicity :', oi_int, 3,"(1x,I3,'/',I2,'/',I3,$)", stdout = std)
  call Message(' ('//ATOM_sym(cell%ATOM_type(i))//')', frm = "(A)", stdout = std)
  call Message('   Equivalent positions  (x y z  occ  DWF) ', frm = "(A)", stdout = std)
  do j=1,cell%numat(i)
    oi_real(1:5) = (/cell%apos(i, j,1:3),dble(cell%ATOM_pos(i,4:5))/)
    call WriteValue('         > ', oi_real, 5,"(2x,4(F9.5,','),F9.5)", stdout = std)
  end do
end do
call Message('', frm = "(A/)", stdout = std)

end subroutine DumpXtalInfo

!--------------------------------------------------------------------------
!
! SUBROUTINE: CrystalData
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief load or generate crystal data
! 
!> @param cell unit cell pointer
!> @param verbose (OPTIONAL)
!> @param fname file name (optional)
!
!> @date    1/ 5/99 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   01/10/14 MDG 4.0 update after new cell type
!> @date   06/06/14 MDG 4.1 added cell pointer and loadingfile arguments
!> @date   03/30/15 MDG 5.0 changed file format to HDF; always assume that the file exists
!--------------------------------------------------------------------------
subroutine CrystalData(cell,verbose)

use io
use crystal
use symmetry

IMPLICIT NONE

type(unitcell), pointer                 :: cell
logical,INTENT(IN),OPTIONAL             :: verbose

integer(kind=irg)                       :: i, ipg, isave

call ReadDataHDF(cell)

! strucdef = .TRUE.
 cell%hexset = .FALSE.
 if (cell%xtal_system.eq.4) cell%hexset = .TRUE.
 if ((cell%xtal_system.eq.5).AND.(cell%SYM_SGset.ne.2)) cell%hexset = .TRUE.

! compute the metric matrices
 call CalcMatrices(cell)

! First generate the point symmetry matrices, then the actual space group.
! Get the symmorphic space group corresponding to the point group
! of the actual space group
 ipg=0
 do i=1,32
  if (SGPG(i).le.cell%SYM_SGnum) ipg=i
 end do

! if the actual group is also the symmorphic group, then both 
! steps can be done simultaneously, otherwise two calls to 
! GenerateSymmetry are needed.
 if (SGPG(ipg).eq.cell%SYM_SGnum) then
  call GenerateSymmetry(cell,.TRUE.)
 else
  isave = cell%SYM_SGnum
  cell%SYM_SGnum = SGPG(ipg)
  call GenerateSymmetry(cell,.TRUE.)
  cell%SYM_SGnum = isave
  call GenerateSymmetry(cell,.FALSE.)
 end if

! and print the information on the screen
if (present(verbose)) then
 if (verbose) then
   call DumpXtalInfo(cell)
 end if
end if 

end subroutine CrystalData

!--------------------------------------------------------------------------
!
! SUBROUTINE: SaveDataHDF
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief save crystal structure data to an HDF file
! 
!> @param cell unit cell pointer
!
!> @date    1/ 5/99 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   01/10/14 MDG 4.0 update after new cell type
!> @date   06/06/14 MDG 4.1 added cell pointer argument
!> @date   03/29/15 MDG 5.0 branched from old version; HDF support 
!> @date   11/07/15 MDG 5.1 correction to writing of SEM_SGset variable
!--------------------------------------------------------------------------
subroutine SaveDataHDF(cell)

use io
use crystal
use HDF5
use HDFsupport
 
IMPLICIT NONE

type(unitcell), pointer, INTENT(IN)     :: cell

type(HDFobjectStackType),pointer        :: HDF_head

character(11)                           :: dstr
character(15)                           :: tstr
character(fnlen)                        :: progname = 'EMmkxtal.f90', groupname, dataset, fname
integer(kind=irg)                       :: hdferr
real(kind=dbl)                          :: cellparams(6)
integer(kind=irg),allocatable           :: atomtypes(:)
real(kind=sgl),allocatable              :: atompos(:,:)

nullify(HDF_head)

! Initialize FORTRAN interface.
!
CALL h5open_f(hdferr)

call timestamp(datestring=dstr, timestring=tstr)

fname = trim(xtalpathname)//trim(cell%fname)
hdferr =  HDF_createFile(fname, HDF_head)

groupname = 'CrystalData'
hdferr = HDF_createGroup(groupname, HDF_head)
dataset = 'ProgramName'
hdferr = HDF_writeDatasetStringArray(dataset, progname, 1, HDF_head)

dataset = 'CreationDate'
hdferr = HDF_writeDatasetStringArray(dataset, dstr, 1, HDF_head)

dataset = 'CreationTime'
hdferr = HDF_writeDatasetStringArray(dataset, tstr, 1, HDF_head)

dataset = 'Creator'
hdferr = HDF_writeDatasetStringArray(dataset, username, 1, HDF_head)

dataset = 'CrystalSystem'
hdferr = HDF_writeDatasetInteger(dataset, cell%xtal_system, HDF_head)

dataset = 'LatticeParameters'
cellparams = (/ cell%a, cell%b, cell%c, cell%alpha, cell%beta, cell%gamma /)
hdferr = HDF_writeDatasetDoubleArray1D(dataset, cellparams, 6, HDF_head)

dataset = 'SpaceGroupNumber'
hdferr = HDF_writeDatasetInteger(dataset, cell%SYM_SGnum, HDF_head)

! make sure we do not write a '0' for the SGset variable; it must be either 1 or 2
if (cell%SYM_SGset.eq.0) cell%SYM_SGset = 1
dataset = 'SpaceGroupSetting'
hdferr = HDF_writeDatasetInteger(dataset, cell%SYM_SGset, HDF_head)

dataset = 'Natomtypes'
hdferr = HDF_writeDatasetInteger(dataset, cell%ATOM_ntype, HDF_head)

allocate(atomtypes(cell%ATOM_ntype))
atomtypes(1:cell%ATOM_ntype) = cell%ATOM_type(1:cell%ATOM_ntype)
dataset = 'Atomtypes'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, atomtypes, cell%ATOM_ntype, HDF_head)
deallocate(atomtypes)

allocate(atompos(cell%ATOM_ntype,5))
atompos(1:cell%ATOM_ntype,1:5) = cell%ATOM_pos(1:cell%ATOM_ntype,1:5)
dataset = 'AtomData'
hdferr = HDF_writeDatasetFloatArray2D(dataset, atompos, cell%ATOM_ntype, 5, HDF_head)
deallocate(atompos)

call HDF_pop(HDF_head,.TRUE.)

call h5close_f(hdferr)

end subroutine SaveDataHDF

!--------------------------------------------------------------------------
!
! SUBROUTINE: ReadDataHDF
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read crystal structure data from an HDF file
! 
!> @param cell unit cell pointer
!
!> @date    1/ 5/99 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   01/10/14 MDG 4.0 update after new cell type
!> @date   06/06/14 MDG 4.1 added cell pointer argument
!> @date   03/29/15 MDG 5.0 branched from old version; HDF support 
!> @date   11/07/15 MDG 5.1 corrected reading of SYM_SGset for older xtal files
!--------------------------------------------------------------------------
subroutine ReadDataHDF(cell)

use io
use crystal
use HDF5
use HDFsupport
 
IMPLICIT NONE

type(unitcell), pointer, INTENT(IN)     :: cell

type(HDFobjectStackType),pointer        :: HDF_head

character(fnlen)                        :: dataset, groupname, fname
integer(HSIZE_T)                        :: dims(1), dims2(2)
integer(kind=irg)                       :: hdferr
real(kind=dbl)                          :: cellparams(6)
integer(kind=irg),allocatable           :: atomtypes(:)
real(kind=sgl),allocatable              :: atompos(:,:)

nullify(HDF_head)

call h5open_f(hdferr)

fname = trim(xtalpathname)//trim(cell%fname)
hdferr =  HDF_openFile(fname, HDF_head)

groupname = 'CrystalData'
hdferr = HDF_openGroup(groupname, HDF_head)

dataset = 'CrystalSystem'
cell%xtal_system = HDF_readDatasetInteger(dataset, HDF_head)

dataset = 'LatticeParameters'
cellparams = HDF_readDatasetDoubleArray1D(dataset, dims, HDF_head)
cell%a = cellparams(1)
cell%b = cellparams(2)
cell%c = cellparams(3)
cell%alpha = cellparams(4)
cell%beta = cellparams(5)
cell%gamma = cellparams(6)

dataset = 'SpaceGroupNumber'
cell%SYM_SGnum = HDF_readDatasetInteger(dataset, HDF_head)

dataset = 'SpaceGroupSetting'
cell%SYM_SGset = HDF_readDatasetInteger(dataset, HDF_head)
! this parameter must be either 1 or 2, but is initialized to 0;
! some older .xtal files may still have 0 in them, so we correct this here
if (cell%SYM_SGset.eq.0) cell%SYM_SGset = 1

dataset = 'Natomtypes'
cell%ATOM_ntype = HDF_readDatasetInteger(dataset, HDF_head)

dataset = 'Atomtypes'
atomtypes = HDF_readDatasetIntegerArray1D(dataset, dims, HDF_head)
cell%ATOM_type(1:cell%ATOM_ntype) = atomtypes(1:cell%ATOM_ntype) 
deallocate(atomtypes)

dataset = 'AtomData'
atompos = HDF_readDatasetFloatArray2D(dataset, dims2, HDF_head)
cell%ATOM_pos(1:cell%ATOM_ntype,1:5) = atompos(1:cell%ATOM_ntype,1:5) 
deallocate(atompos)

call HDF_pop(HDF_head,.TRUE.)

call h5close_f(hdferr)

! for trigonal space groups we need to set SYM_trigonal to .TRUE.
if ((cell%SYM_SGnum.ge.143).and.(cell%SYM_SGnum.le.167)) then
  cell%SG%SYM_trigonal = .TRUE.
else
  cell%SG%SYM_trigonal = .FALSE.
end if 

! we have not yet implemented the rhombohedral setting of the trigonal 
! space groups, so this needs to remain at .FALSE. always.
cell%SG%SYM_second = .FALSE.
!if (cell%SYM_SGset.ne.0) cell%SG%SYM_second=.TRUE.

end subroutine ReadDataHDF


!--------------------------------------------------------------------------
!
! SUBROUTINE: ReadDataFile ---> obsolete and removed 03/30/15
!
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE: SafeOpenFile  ---> obsolete and removed 03/30/15
!
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE: SafeCloseFile ---> obsolete and removed 03/30/15
!
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE: CopyTemplateFiles
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  copy template files into local folder
!
!> @details In the resources folder, there is a text file called templatecodes.txt
!> in which each template file is given a unique ID number.  The present routine
!> receives the requested numbers and then looks into the file to figure out 
!> which ones need to be copied.
!
!> @param nt number of files to copy
!> @param templatelist integer array identifying the template files to be copied
! 
!> @date   06/26/13 MDG 1.0 first attempt
!> @date   06/08/14 MDG 2.0 added stdout argument
!> @date   10/06/14 MDG 2.1 corrected off-by-one error in templatelist
!> @date   05/05/15 MDG 2.2 removed getenv() call; replaced by global path string
!--------------------------------------------------------------------------
subroutine CopyTemplateFiles(nt,templatelist,stdout)

use io

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: nt
integer(kind=irg),INTENT(IN)            :: templatelist(*)
integer(kind=irg),INTENT(IN),OPTIONAL   :: stdout

character(fnlen)                :: templates(maxnumtemplates)
character(fnlen)                :: input_name, output_name
integer(kind=sgl)               :: ios, i, j, std
character(255)                  :: line

std = 6
if (PRESENT(stdout)) std = stdout

! first open and read the resources/templatecodes.txt file

open(UNIT=dataunit,FILE=trim(templatecodefilename), &
        STATUS='old', FORM='formatted',ACCESS='sequential')

templates = ''
do
 read(dataunit,'(I3.3,A)',iostat=ios) j, line
! write (*,*) j,' ->',trim(line),'<-'
 if (ios.ne.0) then 
  exit
 end if
 templates(j+1) = trim(line)
end do
CLOSE(UNIT=dataunit, STATUS='keep')

do i=1,nt
 input_name = trim(templatepathname)//trim(templates(templatelist(i)+1))
 output_name = templates(templatelist(i)+1)
 OPEN(UNIT=dataunit,FILE=trim(input_name), STATUS='old', FORM='formatted',ACCESS='sequential')
 OPEN(UNIT=dataunit2,FILE=trim(output_name), STATUS='unknown', FORM='formatted',ACCESS='sequential')
 do
        read(dataunit,'(A)',iostat=ios) line
        if (ios.ne.0) then 
          exit
        end if
        write(dataunit2,'(A)') trim(line)
  end do
 CLOSE(UNIT=dataunit, STATUS='keep')
 CLOSE(UNIT=dataunit2, STATUS='keep')
 call Message('  -> created template file '//trim(templates(templatelist(i)+1)), frm = "(A)", stdout = std)
end do
 
end subroutine CopyTemplateFiles


!--------------------------------------------------------------------------
!
! SUBROUTINE: Interpret_Program_Arguments
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  interpret the command line arguments
!
!> @details In the resources folder, there is a text file called templatecodes.txt
!> in which each template file is given a unique ID number.  The present routine
!> receives the requested numbers and then looks into the file to figure out 
!> which ones need to be copied.
!
!> @param nmldefault default nml file name on input, actual file name on output
!> @param numt number of files to (potentially) copy
!> @param templatelist integer array identifying the template files to be copied
! 
!> @date   06/26/13 MDG 1.0 first attempt (this might belong elsewhere...)
!> @date   09/04/13 MDG 1.1 minor modifications to arguments
!> @date   06/08/14 MDG 2.0 added stdout argument
!--------------------------------------------------------------------------
subroutine Interpret_Program_Arguments(nmldefault,numt,templatelist,progname,stdout)

use io

IMPLICIT NONE

character(fnlen),INTENT(INOUT)          :: nmldefault
integer(kind=irg),INTENT(IN)            :: numt
integer(kind=irg),INTENT(IN)            :: templatelist(*)
character(fnlen),INTENT(IN)             :: progname
integer(kind=irg),INTENT(IN),OPTIONAL   :: stdout

integer(kind=irg)                       :: numarg       !< number of command line arguments
integer(kind=irg)                       :: iargc        !< external function for command line
character(fnlen)                        :: arg          !< to be read from the command line
character(fnlen)                        :: nmlfile      !< nml file name
integer(kind=irg)                       :: i, std
logical                                 :: haltprogram

std = 6
if (PRESENT(stdout)) std = stdout

numarg = iargc()
nmlfile = ''
nmlfile = trim(nmldefault)

haltprogram = .FALSE.
if (numarg.ge.1) haltprogram = .TRUE.

if (numarg.gt.0) then ! there is at least one argument
  do i=1,numarg
    call getarg(i,arg)
!    mess = 'Found the following argument: '//trim(arg); call Message("(/A/)")
! does the argument start with a '-' character?    
    if (arg(1:1).eq.'-') then
        if (trim(arg).eq.'-h') then
         call Message(' Program should be called as follows: ', frm = "(/A)", stdout = std)
         call Message('        '//trim(progname)//' -h -t [nmlfile]', frm = "(A)", stdout = std)
         call Message(' where nmlfile is an optional file name for the namelist file;', frm = "(A/)", stdout = std)
         call Message(' If absent, the default name '''//trim(nmldefault)//''' will be used.', frm = "(A)", stdout = std)
         call Message(' To create templates of all possible input files, type '//trim(progname)//' -t', frm = "(A)", stdout = std)
         call Message(' To produce this message, type '//trim(progname)//' -h', frm = "(A)", stdout = std)
         call Message(' All program arguments can be combined in the same order;  ', frm = "(A)", stdout = std)
         call Message(' the argument without - will be interpreted as the input file name.', frm = "(A/)", stdout = std)
        end if
        if (trim(arg).eq.'-t') then
! with this option the program creates template namelist files in the current folder so that the 
! user can edit them (file extension will be .template; should be changed by user to .nml)
                call Message('Creating program template files:', frm = "(/A)", stdout = std)
                call CopyTemplateFiles(numt,templatelist)
        end if
    else
! no, the first character is not '-', so this argument must be the filename
! if it is present, but any of the other arguments were present as well, then
! we stop the program. 
        nmlfile = arg
        if (numarg.eq.1) haltprogram = .FALSE.
    end if
  end do
end if

if (haltprogram) then
  call Message('To execute program, remove all flags except for nml input file name', frm = "(/A/)", stdout = std)
  stop
end if

nmldefault = nmlfile

end subroutine Interpret_Program_Arguments




end module files
