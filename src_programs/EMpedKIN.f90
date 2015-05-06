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
! EMsoft:EMpedKIN.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMpedKIN 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Kinematical precession electron diffraction dictionary creation
!
!> @date 03/01/15 MDG 1.0 original
!> @date 03/03/15 MDG 1.1 first tests with realistic parameters; reasonable results
!--------------------------------------------------------------------------
program EMpedKIN

use local
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(PEDKINNameListType)          :: pednl

nmldeffile = 'EMpedKIN.nml'
progname = 'EMpedKIN.f90'
progdesc = 'Kinematical Precession Electron Diffraction Dictionary Generation'

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 15 /), progname)

! deal with the namelist stuff
call GetPEDKINNameList(nmldeffile,pednl)

! generate a set of master EBSD patterns
 call PEDKIN_dictionary(pednl,progname)

end program EMPEDKIN

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
!> @date 03/06/15 MDG 1.1 testing complete; added EulerAngles.txt output
!--------------------------------------------------------------------------
subroutine PEDKIN_dictionary(pednl,progname)

use local
use typedefs
use dictmod
use crystal
use initializers
use gvectors
use io
use diffraction
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
real(kind=sgl)                  :: la, dval, dmin, glen, gmax, io_real(3), om(3,3), k(3), sgmax, FN(3), xgmin, Ig, Igmax, & 
                                   maxint, w, ku(3), kp(3), rnmpp, dx, dy, eu(3)
integer(kind=irg)               :: gp(3), imh, imk, iml, nref, gg(3), ix, iy, iz, io_int(3), ww, nsize, tdp, sx, sy
logical                         :: verbose

real(kind=sgl),allocatable      :: pedpattern(:,:), image(:,:), xx(:,:), yy(:,:), line(:), dot(:,:)

type(unitcell),pointer          :: cell
type(DynType),save              :: Dyn
type(gnode),save                :: rlp
type(reflisttype),pointer       :: reflist, nexts, rltmpa


sgmax = 0.50

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
write (*,*) 'point group number = ', pgnum

!=============================================
!=============================================
! rotation sampling section
! we need to get a sampling of orientation space, starting from the 
! cubochoric representation
nullify(FZlist)
FZcnt = 0
write (*,*) 'pgnum = ', pgnum
write (*,*) 'N = ',pednl%ncubochoric

call sampleRFZ(pednl%ncubochoric, pgnum, FZcnt, FZlist)

io_int(1) = FZcnt
call WriteValue(' Number of incident beam directions       : ', io_int, 1, "(I8)")

!=============================================
!=============================================
! generation of all potential reflections inside a reciprocal space sphere
! computed from the camera length and the detector size ...

! first set the maximum |g| value that can possibly give rise to a diffracted beam on the detector (diagonal)
  gmax = sqrt(2.0) * float(pednl%npix) * pednl%rnmpp
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
  rltmpa%xg = 0.0
  xgmin = 100000.0
  Igmax = 0.0

! now compute |U_g|^2 for all allowed reflections; 
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
        if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
         gg = (/ ix, iy, iz /)
         glen = CalcLength(cell, float(gg), 'r' )

! find all reflections, ignoring double diffraction spots
         if ((IsGAllowed(cell,gg)).and.(glen.le.gmax).and.(glen.gt.0.0)) then ! allowed by the lattice centering, if any
            call AddReflection(rltmpa, reflist, cell, nref, gg )
! we'll use the sangle field of the rltail structure to store |Ug|^2; we will also need the extinction distance
            rltmpa%sangle = cdabs(cell%LUT(ix, iy, iz))**2
            if (rltmpa%sangle.gt.Igmax) Igmax = rltmpa%sangle
            rltmpa%xg = 1.0/(cdabs(cell%LUT(ix,iy,iz))*cell%mLambda)
            if (rltmpa%xg.lt.xgmin) xgmin = rltmpa%xg
         end if ! IsGAllowed
        end if
       end do izl
      end do iyl
    end do ixl
    
io_int(1) = nref
call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")


!=============================================
!=============================================
! create the coordinate arrays for the Gaussian peaks
rnmpp = 1.0/pednl%rnmpp
ww = 4
tdp = 2*ww+1
allocate(xx(-ww:ww,-ww:ww), yy(-ww:ww,-ww:ww), line(-ww:ww), dot(-ww:ww,-ww:ww))
line = (/ (float(i),i=-ww,ww) /) * rnmpp
xx = spread(line,dim=1,ncopies=2*ww+1)
yy = transpose(xx)


!=============================================
!=============================================
! create the output array
nsize = pednl%npix/2 + ww 
allocate(pedpattern(-nsize:nsize,-nsize:nsize))
allocate(image(pednl%npix,pednl%npix))
maxint = Igmax 
write (*,*) ' Maximum intensity = ',maxint

!=============================================
!=============================================
! open the output files, one for the patterns, another one for the Euler angle triplets
open(unit=dataunit2,file=trim(pednl%eulername),status='unknown',form='formatted')
write (dataunit2,"(I6)") FZcnt
FZtmp => FZlist                        ! point to the top of the list
do i=1,FZcnt
  eu = ro2eu(FZtmp%rod)
  write (dataunit2,"(3F10.5)") eu(1), eu(2), eu(3)
  FZtmp => FZtmp%next                  ! point to the next entry
end do
close(unit=dataunit2,status='keep')

open(unit=dataunit,file=trim(pednl%outname),status='unknown',form='unformatted')
write (dataunit) pednl%npix, FZcnt

!=============================================
!=============================================
! and loop over all orientations
FZtmp => FZlist                        ! point to the top of the list
orientationloop: do i = 1, FZcnt       ! loop over all incident beam directions

! set the output pattern to zero
  pedpattern = 0.0
  image = 0.0

! convert the rodrigues vector to a passive rotation matrix.
  om = ro2om(FZtmp%rod)                

! multiplication with (0,0,1) produces the normalized beam direction in a
! cartesian reference frame; so now we can compute the excitation errors 
! for every reflection and keep only the ones that are sufficiently small
  k = (/ 0.0, 0.0, 1.0 /)
  ku = matmul(om,k)
  FN = ku
  k = ku/sngl(cell%mLambda)

! first we go through the entire reflection list and compute the excitation errors
! those points that satisfy the cutoff are linked via the nexts pointers
  rltmpa => reflist%next
  nexts => rltmpa
  do j=1,nref
    gg = rltmpa%hkl
    rltmpa%sg = Calcsg(cell,float(gg),k,FN)
! should we consider this point any further ? If so, add it to the strong reflection linked list
    if (abs(rltmpa%sg).le.sgmax) then 
      nexts%nexts => rltmpa
      nexts => rltmpa
    end if
    rltmpa => rltmpa%next
  end do

! then, for each point in the nexts list, we compute the components of k' = k+g+s
! and place them in the proper reference frame; we skip the incident beam since it is 
! meaningless in the kinematical approximation
  nexts => reflist%next%nexts
  do 
! determine the vector k'
    kp = k + float(nexts%hkl) + nexts%sg*ku
    kp = matmul(transpose(om),kp)

! get the intensity for each point
    w = sngl(cPi)*nexts%sg*pednl%thickness
    if (abs(w).lt.1.0e-6) then
      Ig = nexts%sangle  ! * (sngl(cPi)*pednl%thickness/nexts%xg)**2
    else 
      Ig = nexts%sangle * (sin(w)/w)**2 ! * (sngl(cPi)*pednl%thickness/nexts%xg)**2
    end if

! determine the spot coordinates on the detector
    x = rnmpp * kp(1)
    y = rnmpp * kp(2)

! and plot that spot as a small Gaussian in the pedpattern array, assuming it falls on the detector.
    if ((abs(x).le.nsize-ww).and.(abs(y).le.nsize-ww)) then
      sx = nint(x)
      sy = nint(y)
      dx = x-sx
      dy = y-sy
      dot = (Ig/Igmax)**0.2 * exp(-((xx-dx)**2+(yy-dy)**2)*0.003)
      pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) = pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) + dot(-ww:ww,-ww:ww)
    end if

! and repeat this until the end of the list
    if (.not. associated(nexts%nexts)) EXIT
    nexts => nexts%nexts
  end do

! save the pedpattern to file
  image(1:pednl%npix,1:pednl%npix) = pedpattern(-nsize+ww:nsize-ww,-nsize+ww:nsize-ww)
  write (dataunit) image

! and write the corresponding euler angles to a file as well


! reset the nexts linked list and start over
  nexts => reflist%next
  rltmpa => nexts%nexts
  do 
    nullify(nexts%nexts)
    if (.not. associated(rltmpa%nexts)) EXIT
    nexts => rltmpa
    rltmpa => rltmpa%nexts
  end do

  FZtmp => FZtmp%next                  ! point to the next entry
end do orientationloop

close(unit=dataunit,status='keep')



end subroutine PEDKIN_dictionary



