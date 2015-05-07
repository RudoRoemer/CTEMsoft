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
! EMsoft:initializers.f90
!--------------------------------------------------------------------------
!
! MODULE: initializers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief several basic initialization routines
!
!> @date 01/10/14 MDG 1.0 new version
!--------------------------------------------------------------------------

module initializers

public :: Initialize_ReflectionList
interface Initialize_ReflectionList
        module procedure Initialize_ReflectionList
        module procedure Initialize_ReflectionList_EwaldSweep
end interface


contains

!--------------------------------------------------------------------------
!
! subroutine: Initialize_Cell
!
!> @author Marc De Graef
!
!> @brief perform all steps to initialize a unit cell type variable
!
!> @param cell unit cell pointer
!> @param xtalname file name for crystal structure
!> @param dmin smallest d-spacing to consider
!> @param voltage accelerating voltage (needed to compute relativistic scattering factors)
!
!> @date 01/10/14 MDG 1.0 original
!> @date 06/10/14 MDG 2.0 rewrite without global variables
!--------------------------------------------------------------------------
subroutine Initialize_Cell(cell,Dyn,rlp,xtalname, dmin, voltage, verbose)

use local
use typedefs
use crystal
use symmetry
use files
use io
use error
use gvectors
use diffraction

IMPLICIT NONE

type(unitcell),pointer                     :: cell
type(DynType),INTENT(INOUT)                :: Dyn
type(gnode),INTENT(INOUT)                  :: rlp
character(fnlen),INTENT(IN)                :: xtalname
real(kind=sgl),INTENT(IN)                  :: dmin
real(kind=sgl),INTENT(IN)                  :: voltage
logical,INTENT(IN),OPTIONAL                :: verbose

integer(kind=irg)                          :: istat, io_int(3), skip
integer(kind=irg)                          :: imh, imk, iml, gg(3), ix, iy, iz
real(kind=sgl)                             :: dhkl, io_real(3), ddt
logical                                    :: loadingfile


! clear the cell variable (set everything to zero)
 call ResetCell(cell)

! load the crystal structure file, which also computes all the important 
! matrices as well as all the symmetry arrays
 cell%SG%SYM_reduce=.TRUE.
 cell%fname = xtalname
 call CrystalData(cell,verbose)
 cell%voltage = dble(voltage)

 skip = 3        ! always use Weickenmeier&Kohl scattering coefficients, including absorptive form factors
 call CalcWaveLength(cell,rlp,skip,verbose)

! generate all atom positions
! if the cell is dostorted, then this is not exactly correct, but it should be close for small distortions
 call CalcPositions(cell,'v')

! compute the range of reflections for the lookup table and allocate the table
! The master list is easily created by brute force
 imh = 1
 do 
   dhkl = 1.0/CalcLength(cell,  (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imh = imh + 1
 end do
 imk = 1
 do 
   dhkl = 1.0/CalcLength(cell, (/0.0_sgl,float(imk),0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imk = imk + 1
 end do
 iml = 1
 do 
   dhkl = 1.0/CalcLength(cell, (/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.dmin) EXIT
   iml = iml + 1
 end do

 if (present(verbose)) then
  if (verbose) then
    io_int = (/ imh, imk, iml /)
    call WriteValue(' Range of reflections along a*, b* and c* = ',io_int,3)
  end if
 end if
  
! the LUT array stores all the Fourier coefficients, so that we only need to compute them once... i.e., here and now
 allocate(cell%LUT(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
 if (istat.ne.0) call FatalError('InitializeCell:',' unable to allocate cell%LUT array')
 cell%LUT = dcmplx(0.D0,0.D0)
 
! allocate an array that keeps track of potential double diffraction reflections
 allocate(cell%dbdiff(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
 if (istat.ne.0) call FatalError('InitializeCell:',' unable to allocate cell%dbdiff array')
 cell%dbdiff = .FALSE.
 ddt = 1.0e-10


! next, we compute the overall lookup table cell%LUT; we do not, at this point, create a 
! list of linked reflections; in the old code, this was done at the same time, but it appears
! it is better to decouple these two computations. In this new approach, we'll compute a much
! shorter linked list based on the incident wave vector direction.

! first, we deal with the transmitted beam
 gg = (/ 0,0,0 /)
 call CalcUcg(cell,rlp,gg)   
 Dyn%Upz = rlp%Vpmod         ! U'0 normal absorption parameter 
 if (present(verbose)) then
  if (verbose) then
   io_real(1) = rlp%xgp
   call WriteValue(' Normal absorption length [nm] = ', io_real, 1)
  end if
 end if
 
! and add this reflection to the look-up table
 cell%LUT(0,0,0) = rlp%Ucg

 if (present(verbose)) then
  if (verbose) then
   call Message('Generating Fourier coefficient lookup table ... ', frm = "(/A,$)")
  end if
 end if
 
! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections,
! since the dynamical matrix uses g-h as its index !!!  
ixl: do ix=-2*imh,2*imh
iyl:  do iy=-2*imk,2*imk
izl:   do iz=-2*iml,2*iml
        gg = (/ ix, iy, iz /)
        if (IsGAllowed(cell,gg)) then  ! is this reflection allowed by lattice centering ?
! add the reflection to the look up table
           call CalcUcg(cell,rlp,gg )
           cell%LUT(ix, iy, iz) = rlp%Ucg
! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
           if (cabs(rlp%Ucg).le.ddt) then 
             cell%dbdiff(ix,iy,iz) = .TRUE.
           end if
        end if ! IsGAllowed
       end do izl
      end do iyl
    end do ixl

  if (present(verbose)) then
   if (verbose) then
    call Message('Done', frm = "(A/)")
   end if
  end if
  

! that's it
end subroutine Initialize_Cell



!--------------------------------------------------------------------------
!
! SUBROUTINE: Initialize_ReflectionList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize the potential reflection list for a given wave vector
!
!> @param cell unit cell pointer
!> @param BetheParameter Bethe potential structure
!> @param FN  foil normal
!> @param k zone axis direction cosines in direct Bravais lattice
!> @param dmin smallest lattice d-spacing to consider
!> @param listroot pointer to top of list (could be cell%reflist)
!> @param nref number of reflections in main list (used to be DynNbeams)
!> @param verbose (optional) used for debugging purposes mostly
!
!> @date 01/10/14 MDG 1.0 original, based on old Compute_ReflectionList
!> @date 01/13/14 MDG 1.1 update for new cell type definition and new Bethe potential criterion
!> @date 06/15/14 MDG 2.0 update for removal of all globals
!> @date 06/16/14 MDG 2.1 added recursive
!> @date 06/23/14 MDG 2.2 replaced Dyn structure by FN
!--------------------------------------------------------------------------
recursive subroutine Initialize_ReflectionList(cell, listroot, BetheParameter, FN, k, dmin, nref, verbose)

use local
use typedefs
use io
use crystal
use constants
use gvectors
use diffraction
use symmetry

IMPLICIT NONE

type(unitcell),pointer                          :: cell
type(reflisttype),pointer                       :: listroot
type(BetheParameterType),INTENT(INOUT)          :: BetheParameter
real(kind=sgl),INTENT(IN)                       :: FN(3)
real(kind=sgl),INTENT(IN)                       :: k(3)
real(kind=sgl),INTENT(IN)                       :: dmin
integer(kind=irg),INTENT(INOUT)                 :: nref
logical,INTENT(IN),OPTIONAL                     :: verbose

integer(kind=irg)                               :: imh, imk, iml, gg(3), ix, iy, iz, i, minholz, RHOLZ, im, istat, N, &
                                                   ig, numr, ir, irsel
real(kind=sgl)                                  :: dhkl, io_real(9), H, g3(3), g3n(3), FNg(3), ddt, s, kr(3), exer, &
                                                   rBethe_i, rBethe_d, sgp, r_g, la, dval
integer(kind=irg)                               :: io_int(3), gshort(3), gp(3)
type(reflisttype),pointer                       :: rltail

! set the truncation parameters
  rBethe_i = BetheParameter%c3          ! if larger than this value, we ignore the reflection completely
  rBethe_d = BetheParameter%sgdbdiff    ! excitation error cutoff for double diffraction reflections
  la = 1.0/sngl(cell%mLambda)
  
! get the size of the lookup table
  gp = shape(cell%LUT)
  imh = (gp(1)-1)/4
  imk = (gp(2)-1)/4
  iml = (gp(3)-1)/4
  
  nullify(listroot)
  nullify(rltail)
 
! transmitted beam has excitation error zero
  gg = (/ 0,0,0 /)
  call AddReflection(rltail, listroot, cell, nref, gg )   ! this guarantees that 000 is always the first reflection
  rltail%sg = 0.0


! now compute |sg|/|U_g|/lambda for the other allowed reflections; if this parameter is less than
! the threshhold, rBethe_i, then add the reflection to the list of potential reflections
! note that this uses the older form of the Bethe Potential truncation parameters for now
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
        if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
         gg = (/ ix, iy, iz /)
         dval = 1.0/CalcLength(cell, float(gg), 'r' )

         if ((IsGAllowed(cell,gg)).and.(dval.gt.dmin)) then ! allowed by the lattice centering, if any
          sgp = Calcsg(cell,float(gg),k,FN)
          if (cell%dbdiff(ix, iy, iz)) then ! potential double diffraction reflection
            if (abs(sgp).le.rBethe_d) then 
              call AddReflection(rltail, listroot, cell, nref, gg )
              rltail%sg = sgp
              rltail%dbdiff = .TRUE.
            end if
          else
            r_g = la * abs(sgp)/cdabs(cell%LUT(ix, iy, iz))
            if (r_g.le.rBethe_i) then 
              call AddReflection(rltail, listroot, cell, nref, gg )
              rltail%sg = sgp
              rltail%dbdiff = .FALSE.
            end if
          end if
         end if ! IsGAllowed
        end if
       end do izl
      end do iyl
    end do ixl
    
  if (present(verbose)) then 
    if (verbose) then 
      io_int(1) = nref
      call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
    end if
  end if

end subroutine Initialize_ReflectionList


!--------------------------------------------------------------------------
!
! SUBROUTINE: Initialize_ReflectionList_EwaldSweep
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize the potential reflection list for a given precession electron diffraction geometry
!
!> @param cell unit cell pointer
!> @param listroot main pointer for reflection linked list (typically reflist in calling routine)
!> @param FN  foil normal
!> @param k zone axis direction cosines in direct Bravais lattice
!> @param nref number of reflections in main list (used to be DynNbeams)
!> @param pedangle precession cone angle in mrad
!> @param goffset offset parameter for reflections to be included outside the Ewald sphere precession volume
!> @param verbose (optional) used for debugging purposes mostly
!
!> @date 01/10/14 MDG 1.0 original, based on old Compute_ReflectionList
!> @date 01/13/14 MDG 1.1 update for new cell type definition and new Bethe potential criterion
!> @date 06/15/14 MDG 2.0 update for removal of all globals
!> @date 06/16/14 MDG 2.1 added recursive
!> @date 06/23/14 MDG 2.2 replaced Dyn structure by FN
!> @date 11/30/14 MDG 3.0 forked from original to accommodate PED case
!--------------------------------------------------------------------------
recursive subroutine Initialize_ReflectionList_EwaldSweep(cell, listroot, FN, k, nref, pedangle, goffset, verbose)

use local
use typedefs
use io
use crystal
use constants
use gvectors
use diffraction
use symmetry

IMPLICIT NONE

type(unitcell),pointer                          :: cell
type(reflisttype),pointer                       :: listroot
real(kind=sgl),INTENT(IN)                       :: FN(3)
real(kind=sgl),INTENT(IN)                       :: k(3)
integer(kind=irg),INTENT(INOUT)                 :: nref
real(kind=sgl),INTENT(IN)                       :: pedangle
real(kind=sgl),INTENT(IN)                       :: goffset
logical,INTENT(IN),OPTIONAL                     :: verbose

integer(kind=irg)                               :: imh, imk, iml, gg(3), ix, iy, iz, io_int(3), gp(3)
real(kind=sgl)                                  :: FNg(3), c, s, kr(3), sgp, la, kstar(3), gperp(3), gpara(3), bup, blo, y, z, &
                                                   gdk, glen, gplen
type(reflisttype),pointer                       :: rltail

! init a couple of parameters
  la = 1.0/sngl(cell%mLambda)
  c = la * cos(pedangle/1000.0)
  s = 2.0 * la * sin(pedangle/1000.0)

! reciprocal space wave vector
  call TransSpace(cell, k, kstar, 'd', 'r')
  call NormVec(cell, kstar, 'r')
  kstar = la * kstar
  
! get the size of the lookup table
  gp = shape(cell%LUT)
  imh = (gp(1)-1)/4
  imk = (gp(2)-1)/4
  iml = (gp(3)-1)/4
  
  nullify(listroot)
  nullify(rltail)
 
! transmitted beam has excitation error zero, and set xg to zero; xg will store the accumulated intensity for each reflection
  gg = (/ 0,0,0 /)
  call AddReflection(rltail, listroot, cell, nref, gg )   ! this guarantees that 000 is always the first reflection
  rltail%sg = 0.0
  rltail%xg = 0.0

! scan through the reciprocal lattice volume corresponding to the dmin value
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
        if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
         gg = (/ ix, iy, iz /)
         if (IsGAllowed(cell,gg)) then ! allowed by the lattice centering, if any
! first we need to determine the parallel and perpendicular components of this g vector with respect to the beam direction in reciprocal space
          gdk = CalcDot(cell,float(gg),kstar,'r')       ! projection of gg onto k*
          gpara = gdk * kstar
          gperp = float(gg) - gpara
! then get the length of the perpendicular and parallel components, including sign of parallel component
          glen = CalcLength(cell, gperp, 'r')
          gplen = CalcLength(cell, gpara, 'r')
        ! sign of length depends on dot product sign of gg onto k*
          if (gdk.le.0.0) gplen = -gplen
! compute the upper and lower bounds for this value of glen
          y = glen*s
          z = c*c-glen*glen
          bup = goffset+c-sqrt(z-y)
          blo = -goffset+c-sqrt(z+y)
! and check whether or not this point should be taken into account
          if ((blo.le.gplen).and.(gplen.le.bup)) then 
            sgp = Calcsg(cell,float(gg),k,FN)
! note that we are not applying any Bethe parameter conditions here since those will be applied for each beam orientation separately
            if (cell%dbdiff(ix, iy, iz)) then ! potential double diffraction reflection
                call AddReflection(rltail, listroot, cell, nref, gg )
                rltail%sg = sgp
                rltail%xg = 0.0
                rltail%dbdiff = .TRUE.
            else
                call AddReflection(rltail, listroot, cell, nref, gg )
                rltail%sg = sgp
                rltail%xg = 0.0
                rltail%dbdiff = .FALSE.
            end if
          end if ! reflection inside precession-swept Ewald sphere volume
         end if ! IsGAllowed
        end if ! not the origin
       end do izl
      end do iyl
    end do ixl
    
  if (present(verbose)) then 
    if (verbose) then 
      io_int(1) = nref
      call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
    end if
  end if

end subroutine Initialize_Reflectionlist_EwaldSweep


end module initializers
