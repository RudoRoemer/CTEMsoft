! ###################################################################
! Copyright (c) 2015, Marc De Graef/Carnegie Mellon University
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
! CTEMsoft2013:CTEMpedKIN.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMpedKIN 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Kinematical precession electron diffraction dictionary creation
!
!> @date 03/01/15 MDG 1.0 original
!--------------------------------------------------------------------------
program CTEMpedKIN

use local
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(PEDKINNameListType)          :: pednl

nmldeffile = 'CTEMpedKIN.nml'
progname = 'CTEMpedKIN.f90'
progdesc = 'Kinematical Precession Electron Diffraction Dictionary Generation'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 15 /), progname)

! deal with the namelist stuff
call GetPEDKINNameList(nmldeffile,pednl)

! print some information
call CTEMsoft(progname, progdesc)

! generate a set of master EBSD patterns
 call PEDKIN_dictionary(pednl,progname)

end program CTEMPEDKIN

!--------------------------------------------------------------------------
!
! SUBROUTINE:PEDKIN_dictionary
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a kinematical PED dictionary
!
!> @param nmlfile namelist file name
!
!> @date 03/02/15 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine PEDKIN_dictionary(pednl,progname)

use local
use typedefs
use dictmod
use crystal
use initializers
use gvectors
use io
use symmetry
use quaternions
use NameListTypedefs
use constants
use rotations
use so3
use math

type(PEDKINNameListType),INTENT(IN)     :: pednl
character(fnlen),INTENT(IN)             :: progname

integer(kind=irg)               :: FZcnt, pgnum
type(FZpointd),pointer          :: FZlist, FZtmp
real(kind=sgl)                  :: la, dval, dmin, glen, gmax, io_real(3), om(3,3), k(3), sgmax, FN(3)
integer(kind=irg)               :: gp(3), imh, imk, iml, nref, gg(3), ix, iy, iz, io_int(3)
logical                         :: verbose

real(kind=sgl),allocatable      :: pedpattern(:,:)

type(unitcell),pointer          :: cell
type(DynType),save              :: Dyn
type(gnode),save                :: rlp
type(reflisttype),pointer       :: reflist, nexts, rltmpa


sgmax = 0.10

!=============================================
!=============================================
! crystallography section
nullify(cell)
allocate(cell)

verbose = .TRUE.

call Initialize_Cell(cell,Dyn,rlp,pednl%xtalname, pednl%dmin, pednl%voltage, verbose)

! determine the point group number
j=0
do i=1,32
 if (SGPG(i).le.cell % SYM_SGnum) j=i
end do
pgnum = j

!=============================================
!=============================================
! generation of all potential reflections inside a reciprocal space sphere
! computed from the camera length and the detector size ...

! first set the maximum |g| value that can possibly give rise to a diffracted beam on the detector (diagonal)
gmax = (float(pednl%npix)*pednl%pixelsize/1000.0/sqrt(2.0)) / sngl(cell%mLambda) / pednl%camlen
io_real(1) = gmax
call WriteValue(' Length of longest g-vector : ', io_real, 1, "(F8.4)")

! this code is taken from the Initialize_ReflectionList routine, but we do not
! need everything from that routine 
! get the size of the lookup table
  gp = shape(cell%LUT)
  imh = (gp(1)-1)/4
  imk = (gp(2)-1)/4
  iml = (gp(3)-1)/4

write (*,*) 'shape of LUT = ',shape(cell%LUT)
  
  nullify(reflist)
  nullify(rltmpa)
  nref = 0
 
! transmitted beam has excitation error zero
  gg = (/ 0,0,0 /)
  call AddReflection(rltmpa, reflist, cell, nref, gg)   ! this guarantees that 000 is always the first reflection

! now compute |U_g|^2 for all allowed reflections; 
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
        if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
         gg = (/ ix, iy, iz /)
         glen = CalcLength(cell, float(gg), 'r' )

! find all reflections, ignoring double diffraction spots
         if ((IsGAllowed(cell,gg)).and.(glen.le.gmax)) then ! allowed by the lattice centering, if any
            call AddReflection(rltmpa, reflist, cell, nref, gg )
! we'll use the sangle field of the rltail structure to store |Ug|^2
            rltmpa%sangle = cdabs(cell%LUT(ix, iy, iz))**2
         end if ! IsGAllowed
        end if
       end do izl
      end do iyl
    end do ixl
    
  io_int(1) = nref
  call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")

!=============================================
!=============================================
! rotation sampling section
! we need to get a sampling of orientation space, starting from the 
! cubochoric representation
nullify(FZlist)
FZcnt = 0
call sampleRFZ(pednl%ncubochoric, pgnum, FZcnt, FZlist)


!=============================================
!=============================================
! create the output array
allocate(pedpattern(pednl%npix,pednl%npix))

!=============================================
!=============================================
! open the output file


!=============================================
!=============================================
! and loop over all orientations
FZtmp => FZlist                        ! point to the top of the list
do i = 1, FZcnt                        ! loop over all entries
  om = ro2om(FZtmp%rod)                ! convert to passive rotation matrix
! multiplication with (0,0,1) produces the normalized beam direction in a
! cartesian reference frame; so now we can compute the excitation errors 
! for every reflection and keep only the ones that are small
  k = (/ 0.0, 0.0, 1.0 /)
  k = matmul(om,k)
  FN = k
  k = k/sngl(mLambda)
! first we go through the entire reflection list and compute the excitation errors
  rltmpa => reflist%next
  nexts => rltmpa
  do i=1,nref
    gg = rltmpa%hkl
    rltmpa%sg = Calcsg(cell,float(gg),k,FN)
! should we consider this point any further ? If so, add it to the strong reflection linked list
    if (abs(rltmpa%sg).le.sgmax) then 
      nexts%nexts => rltmpa
    end if
    rltmpa => rltmpa%next
  end do






  FZtmp => FZtmp%next                  ! point to the next entry
end do




end subroutine PEDKIN_dictionary



