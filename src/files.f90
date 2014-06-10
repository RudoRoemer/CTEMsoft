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
! CTEMsoft:files.f90
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
!> @param loadingfile logical
!> @param fname file name (optional)
!
!> @date    1/ 5/99 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   01/10/14 MDG 4.0 update after new cell type
!> @date   06/06/14 MDG 4.1 added cell pointer and loadingfile arguments
!--------------------------------------------------------------------------
subroutine CrystalData(cell, loadingfile, fname, stdout)

use io
use crystal
use symmetry

IMPLICIT NONE

type(unitcell), pointer                 :: cell
logical, INTENT(INOUT)                  :: loadingfile
character(fnlen),OPTIONAL,INTENT(IN)    :: fname        !< optional file name
integer(kind=irg),INTENT(IN),OPTIONAL   :: stdout

integer(kind=irg)                       :: io_int(1), std
logical                                 :: fr = .TRUE.

 std = 6
 if (PRESENT(stdout)) std = stdout

 loadingfile = .FALSE.
 if (PRESENT(fname)) then 
  cell%fname = fname
  loadingfile = .TRUE.
  call ReadDataFile(cell,fr)
 else
  call ReadValue(' Load file (0) or new data (1) ? ', io_int, 1)
  call Message('', frm = "(A/)", stdout = std)
  if (io_int(1).ne.0) then
   cell%SYM_SGset=0
   call GetLatParm(cell)
   call CalcMatrices(cell)
   call GetSpaceGroup(cell)
   call GenerateSymmetry(cell,.TRUE.)
   call GetAsymPos(cell)
   call SaveData(cell)
  else
   loadingfile = .TRUE.
   call ReadDataFile(cell)
  end if
 end if

end subroutine CrystalData

!--------------------------------------------------------------------------
!
! SUBROUTINE: SaveData
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief save crystal structure data to file
! 
!> @todo This should be replaced with a text-based editable file format (.txt. or .xml)
!
!> @param cell unit cell pointer
!
!> @date    1/ 5/99 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   01/10/14 MDG 4.0 update after new cell type
!> @date   06/06/14 MDG 4.1 added cell pointer argument
!--------------------------------------------------------------------------
subroutine SaveData(cell)

use io
use crystal
 
IMPLICIT NONE

type(unitcell), pointer :: cell

! call SafeOpenFile('xt','unformatted',cell%fname)
 open (dataunit,file=trim(cell%fname),status='unknown',form='unformatted')
! save lattice parameters, crystal system, space group number and contents
! of the asymmetric unit.
 write (dataunit) cell%xtal_system, cell%a,cell%b,cell%c,cell%alpha,cell%beta,cell%gamma
 write (dataunit) cell%ATOM_type, cell%ATOM_pos, cell%ATOM_ntype, cell%SYM_SGnum, cell%SYM_SGset
 call SafeCloseFile('xt','keep', cell%fname)

end subroutine

!--------------------------------------------------------------------------
!
! SUBROUTINE: ReadDataFile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  load crystal data in memory
!
!> @details Also converts the older file format to the current one.
!
!> @param cell unit cell pointer
!> @param fr optional logical, passed on to SafeOpenFile routine
! 
!> @todo This should be replaced with a text-based editable file format (.txt. or .xml)
!
!> @date    1/ 5/99 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   01/10/14 MDG 4.0 update after new cell type
!> @date   06/06/14 MDG 4.1 added cell pointer argument
!--------------------------------------------------------------------------
subroutine ReadDataFile(cell,fr)

use io
use crystal
use symmetry

IMPLICIT NONE

type(unitcell), pointer                 :: cell
logical,optional,INTENT(IN)             :: fr                   !< logical 

integer(kind=irg)                       :: i,ipg,isave,iost
logical                                 :: fread = .TRUE.

! the following variables are used for the older (smaller) xtal file format
integer                                 :: ATOM_type(50),ATOM_ntype,SYM_SGnum,SYM_SGset
real                                    :: ATOM_pos(50,5)


 if (present(fr)) fread=.FALSE.
  
 call SafeOpenFile('xt','unformatted',cell%fname,fread)
 read (dataunit) cell%xtal_system, cell%a,cell%b,cell%c,cell%alpha,cell%beta,cell%gamma
 read (dataunit,iostat=iost) cell%ATOM_type, cell%ATOM_pos, cell%ATOM_ntype, cell%SYM_SGnum, cell%SYM_SGset

 if (iost.eq.0) then 
   call SafeCloseFile('xt','keep',cell%fname,.TRUE.)   ! ok, this is a modern file with 100 max atom positions in asymmetric unit
 else ! this a legacy input file with only 50 maximum atom positions in the asymmetric unit cell... conversion is needed
! close the file first, then re-open it and read with alternative variables
   call Message('Found a data file of the old format; converting to new format', frm = "(A)")
   call SafeCloseFile('xt','keep',cell%fname,.TRUE.)     
   call SafeOpenFile('xt','unformatted',cell%fname,.FALSE.)
   read (dataunit) cell%xtal_system, cell%a,cell%b,cell%c,cell%alpha,cell%beta,cell%gamma
   read (dataunit,iostat=iost) ATOM_type, ATOM_pos, ATOM_ntype, SYM_SGnum, SYM_SGset
   call SafeCloseFile('xt','keep',cell%fname)     
! and copy them into the new arrays
   cell%ATOM_type(1:50) = ATOM_type(1:50)
   cell%ATOM_pos(1:50,1:5) = ATOM_pos(1:50,1:5)
   cell%ATOM_ntype = ATOM_ntype
   cell%SYM_SGnum = SYM_SGnum
   cell%SYM_SGset = SYM_SGset
   call SaveData(cell)
   call Message('Conversion complete', frm = "(A)")
 end if
 
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
 call DumpXtalInfo(cell)

end subroutine

!--------------------------------------------------------------------------
!
! SUBROUTINE: SafeOpenFile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  open a data file
!
!> @details based on example 10-1 in Stephen Chapman''s Fortran 90/95 book.
!
!> @param ftyp file type identifier
!> @param frm string to indicate file format (formatted, unformatted)
!> @param fname file name
!> @param fread optional logical
!> @param stdout optional output unit identifier
! 
!> @todo This should be replaced with a text-based editable file format (.txt. or .xml)
!
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   06/06/14 MDG 4.0 added stdout argument
!--------------------------------------------------------------------------
subroutine SafeOpenFile(ftyp,frm,fname,fread, loadingfile, stdout)

use error
use io

IMPLICIT NONE

character(2),INTENT(IN)                 :: ftyp         !< selects the logical unit
character(fnlen),INTENT(INOUT)          :: fname        !< file name string
character(*),INTENT(IN)                 :: frm          !< formatting option
logical,OPTIONAL,INTENT(IN)             :: fread        !< read if .TRUE., write if .FALSE.
logical,INTENT(INOUT),OPTIONAL          :: loadingfile
integer(kind=irg),INTENT(IN),OPTIONAL	  :: stdout

character(1)                            :: ans
logical                                 :: lexist,lopen 
integer(kind=irg)                       :: iunit, std

 std = 6
 if (PRESENT(stdout)) std = stdout

 select case (ftyp)
  case('xt'); iunit = dataunit
  case('d1'); iunit = dataunit
  case('d2'); iunit = dataunit2
  case('d3'); iunit = dataunit3
  case('ps'); iunit = psunit
 end select

 lopen = .FALSE.
 openfile: do while(.not.lopen)  ! repeat this attempt until we have an open file

! get filename
  if (PRESENT(fread)) then
   if (fread.eqv..TRUE.) then ! get the input filename 
    call ReadValue(' Enter input file name : ', fname,"(A)", stdout = std)
   end if
  else
   call ReadValue(' Enter output file name : ', fname, "(A)", stdout = std)
  end if

! does file already exist ?
  inquire(file=trim(fname),exist=lexist)
  if ((.not.lexist).and.(loadingfile)) call FatalError('SafeOpenFile ',' input file does not exist ')
  
  exists: if (.not.lexist) then
! ok, file does not already exist, so open it
    open(unit=iunit,file=trim(fname),status='new',action='write',form=frm)
    lopen = .TRUE.
  else
! file exists.  Should we replace it or read it ?
   if (PRESENT(fread)) then  ! file exists and we op to read it
     open(unit=iunit,file=trim(fname),status='old',action='read',form = frm)
     lopen = .TRUE.
   else  ! file exists and we need to write, so ask what to do
    call ReadValue(' Output file exists.  Overwrite it ? (y/n) ', ans,"(A1)", stdout = std)
    replace: if ((ans == 'Y').or.(ans == 'y')) then
! open file for writing 
     open(unit=iunit,file=trim(fname),status='replace',action='write',form = frm)
     lopen = .TRUE.
    end if replace
   end if
  end if exists

 end do openfile


 if (lopen) then ! ok, file is open, so we print a message
  if (PRESENT(fread)) then
   call Message('  File '//trim(fname)//' open for input ', frm = "(A)", stdout = std)
  else
   call Message('  File '//trim(fname)//' open for output ', frm = "(A)", stdout = std)
  end if
 end if

end subroutine SafeOpenFile

!--------------------------------------------------------------------------
!
! SUBROUTINE: SafeCloseFile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  close a data file
!
!> @param ftyp file type identifier
!> @param stt string to indicate file status (discard or keep)
!> @param fname file name
!> @param quiet optional logical to prevent output
!> @param stdout output unit identifier
! 
!> @todo This should be replaced with a text-based editable file format (.txt. or .xml)
!
!> @date   01/05/99 MDG 1.0 original
!> @date   05/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!> @date   06/08/14 MDG 4.0 added stdout argument
!--------------------------------------------------------------------------
subroutine SafeCloseFile(ftyp,stt,fname,quiet,stdout)

use io
use error

IMPLICIT NONE

character(2),INTENT(IN)                 :: ftyp         !< file type identifier
character(*),INTENT(IN)                 :: stt          !< status string
character(fnlen),INTENT(IN)             :: fname        !< file name
logical,OPTIONAL,INTENT(IN)             :: quiet        !< logical verbose or not
integer(kind=irg),OPTIONAL,INTENT(IN)   :: stdout       !< output unit identifier 

integer(kind=irg)                       :: iunit,ier,std

std = 6
if (PRESENT(stdout)) std = stdout

 select case (ftyp)
  case('xt'); iunit = dataunit
  case('d1'); iunit = dataunit
  case('d2'); iunit = dataunit2
  case('d3'); iunit = dataunit3
  case('ps'); iunit = psunit
 end select

 close(unit=iunit,status=stt,iostat=ier)

 if (ier.ne.0) call FatalError('SafeCloseFile',' data file NOT saved', stdout = std)

 if (.not.present(quiet)) then
  call Message(' Data has been stored in file '//trim(fname), frm = "(A)", stdout = std)
 end if
 
end subroutine SafeCloseFile

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

call getenv("CTEMsoft2013templates",templatepathname)
call getenv("CTEMsoft2013resources",resourcepathname)
!write (*,*) 'template folder location: ',templatepathname
!write (*,*) 'resource folder location: ',resourcepathname

! first open and read the resources/templatecodes.txt file
open(UNIT=dataunit,FILE=trim(resourcepathname)//'/'//templatecodefilename, &
        STATUS='old', FORM='formatted',ACCESS='sequential')

templates = ''
do
 read(dataunit,'(I3.3,A)',iostat=ios) j, line
! write (*,*) j,' ->',trim(line),'<-'
 if (ios.ne.0) then 
  exit
 end if
 templates(j) = trim(line)
end do
CLOSE(UNIT=dataunit, STATUS='keep')

do i=1,nt
 input_name = trim(templatepathname)//'/'//trim(templates(templatelist(i)))
! write (*,*) i,'->',trim(input_name),'<-'
 output_name = templates(templatelist(i))
 OPEN(UNIT=dataunit,FILE=trim(input_name), STATUS='old', FORM='formatted',ACCESS='sequential')
 OPEN(UNIT=dataunit2,FILE=trim(output_name), STATUS='replace', FORM='formatted',ACCESS='sequential')
 do
        read(dataunit,'(A)',iostat=ios) line
        if (ios.ne.0) then 
          exit
        end if
        write(dataunit2,'(A)') trim(line)
  end do
 CLOSE(UNIT=dataunit, STATUS='keep')
 CLOSE(UNIT=dataunit2, STATUS='keep')
 call Message('  -> created template file '//trim(templates(templatelist(i))), frm = "(A)", stdout = std)
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
