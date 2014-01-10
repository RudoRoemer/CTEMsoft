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
! CTEMsoft:initializers.f90
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

contains

!--------------------------------------------------------------------------
!
! subroutine: InitializeCell
!
!> @author Marc De Graef
!
!> @brief perform all steps to initialize a unit cell type variable
!
!> @param xtalname file name for crystal structure
!> @param dmin smallest d-spacing to consider
!> @param voltage accelerating voltage (needed to compute relativistic scattering factors)
!
!> @date 01/10/14 MDG original
!--------------------------------------------------------------------------
subroutine InitializeCell(xtalname, dmin, voltage)

use local
use crystalvars
use crystal
use files
use io
use gvectors
use dynamical

IMPLICIT NONE

character(fnlen),INTENT(IN)                :: xtalname
real(kind=sgl),INTENT(IN)                  :: dmin
real(kind=sgl),INTENT(IN)                  :: voltage

integer(kind=irg)                          :: istat, io_int(3), skip
integer(kind=irg)	                    :: imh, imk, iml, gg(3), ix, iy, iz
real(kind=sgl)                             :: dhkl, io_real(3), ddt


! make sure the cell variable exists
 if (.not.associated(cell)) then
  allocate(cell,stat=istat)
  if (istat.ne.0) call FatalError('InitializeCell:',' unable to allocate cell pointer')
 end if

! clear the cell variable (set everthing to zero)
 call ResetCell

! load the crystal structure file, which also computes all the important 
! matrices as well as all the symmetry arrays
 cell%SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)

 skip = 3        ! always use Weickenmeier&Kohl scattering coefficients, including absorptive form factors
 call CalcWaveLength(dble(voltage),skip)

! compute the range of reflections for the lookup table and allocate the table
! The master list is easily created by brute force
 imh = 1
 do 
   imh = imh + 1
   dhkl = 1.0/CalcLength(  (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
 end do
 imk = 1
 do 
   imk = imk + 1
   dhkl = 1.0/CalcLength( (/0.0_sgl,float(imk),0.0_sgl/), 'r')
  if (dhkl.lt.dmin) EXIT
 end do
 iml = 1
 do 
   iml = iml + 1
   dhkl = 1.0/CalcLength( (/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.dmin) EXIT
 end do
 io_int = (/ imh, imk, iml /)
 call WriteValue(' Range of reflections along a*, b* and c* = ',io_int,3)

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
! shorted linked list based on the incident wave vector direction.

! first, we deal with the transmitted beam
 gg = (/ 0,0,0 /)
 call CalcUcg(gg)   
 DynUpz = rlp%Vpmod         ! U'0 normal absorption parameter 
 io_real(1) = rlp%xgp
 call WriteValue(' Normal absorption length [nm] = ', io_real, 1)

! and add this reflection to the look-up table
 cell%LUT(0,0,0) = rlp%Ucg

 mess = 'Generating Fourier coefficient lookup table ... '
 call Message("(/A,$)")

! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections,
! since the dynamical matrix uses g-h as its index !!!  
ixl: do ix=-2*imh,2*imh
iyl:  do iy=-2*imk,2*imk
izl:   do iz=-2*iml,2*iml
        gg = (/ ix, iy, iz /)
        if (IsGAllowed(gg)) then  ! is this reflection allowed by lattice centering ?
! add the reflection to the look up table
 	   call CalcUcg( gg )
           cell%LUT(ix, iy, iz) = rlp%Ucg
! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
           if (cabs(rlp%Ucg).le.ddt) then 
             cell%dbdiff(ix,iy,iz) = .TRUE.
           end if
        end if ! IsGAllowed
       end do izl
      end do iyl
    end do ixl
 mess = 'Done'
 call Message("(A/")

! generate all atom positions
 call CalcPositions('v')

! that's it
end subroutine InitializeCell




end module initializers
