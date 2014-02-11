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
! subroutine: Initialize_Cell
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
subroutine Initialize_Cell(xtalname, dmin, voltage, verbose)

use local
use crystalvars
use crystal
use symmetry
use files
use io
use error
use gvectors
use diffraction
use dynamical

IMPLICIT NONE

character(fnlen),INTENT(IN)                :: xtalname
real(kind=sgl),INTENT(IN)                  :: dmin
real(kind=sgl),INTENT(IN)                  :: voltage
logical,INTENT(IN),OPTIONAL	             :: verbose
integer(kind=irg)                          :: istat, io_int(3), skip
integer(kind=irg)	                    :: imh, imk, iml, gg(3), ix, iy, iz
real(kind=sgl)                             :: dhkl, io_real(3), ddt


! make sure the cell variable exists
! if (.not.associated(cell)) then
!  allocate(cell,stat=istat)
!  if (istat.ne.0) call FatalError('InitializeCell:',' unable to allocate cell pointer')
! end if

! clear the cell variable (set everything to zero)
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
   dhkl = 1.0/CalcLength(  (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imh = imh + 1
 end do
 imk = 1
 do 
   dhkl = 1.0/CalcLength( (/0.0_sgl,float(imk),0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imk = imk + 1
 end do
 iml = 1
 do 
   dhkl = 1.0/CalcLength( (/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.dmin) EXIT
   iml = iml + 1
 end do

 if (present(verbose)) then
    io_int = (/ imh, imk, iml /)
   call WriteValue(' Range of reflections along a*, b* and c* = ',io_int,3)
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
! shorted linked list based on the incident wave vector direction.

! first, we deal with the transmitted beam
 gg = (/ 0,0,0 /)
 call CalcUcg(gg)   
 DynUpz = rlp%Vpmod         ! U'0 normal absorption parameter 
 if (present(verbose)) then
   io_real(1) = rlp%xgp
   call WriteValue(' Normal absorption length [nm] = ', io_real, 1)
 end if
 
! and add this reflection to the look-up table
 cell%LUT(0,0,0) = rlp%Ucg

 if (present(verbose)) then
  mess = 'Generating Fourier coefficient lookup table ... '
  call Message("(/A,$)")
 end if
 
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

  if (present(verbose)) then
   mess = 'Done'
   call Message("(A/)")
  end if
  
! generate all atom positions
 call CalcPositions('v')

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
!> @param k zone axis direction cosines in direct Bravais lattice
!> @param dmin smallest lattice d-spacing to consider
!> @param verbose (optional) used for debugging purposes mostly
!
!> @date 01/10/14 MDG 1.0 original, based on old Compute_ReflectionList
!> @date 01/13/14 MDG 1.1 update for new cell type definition and new Bethe potential criterion
!--------------------------------------------------------------------------
subroutine Initialize_ReflectionList(k,dmin,verbose)

use local
use io
use crystal
use constants
use dynamical
use gvectors
use diffraction
use symmetry

IMPLICIT NONE

real(kind=sgl),INTENT(IN)			:: k(3)
real(kind=sgl),INTENT(IN)                     :: dmin
logical,INTENT(IN),OPTIONAL	               :: verbose

integer(kind=irg)				:: imh, imk, iml, gg(3), ix, iy, iz, i, minholz, RHOLZ, im, istat, N, &
						ig, numr, ir, irsel
real(kind=sgl)					:: dhkl, io_real(9), H, g3(3), g3n(3), FNg(3), ddt, s, kr(3), exer, &
                                                rBethe_i, rBethe_d, sgp, r_g, la, dval
integer(kind=irg)				:: io_int(3), gshort(3), gp(3)


! set the truncation parameters 
  rBethe_i = BetheParameter%c3          ! if larger than this value, we ignore the reflection completely
  rBethe_d = BetheParameter%sgdbdiff    ! excitation error cutoff for double diffraction reflections
  la = 1.0/sngl(mLambda)
  
! get the size of the lookup table
  gp = shape(cell%LUT)
  imh = (gp(1)-1)/4
  imk = (gp(2)-1)/4
  iml = (gp(3)-1)/4
  
! get rid of the existing linked list, if any 
  if (associated(cell%reflist)) then
    call Delete_gvectorlist()
    nullify(cell%reflist)
  end if 
  
! transmitted beam has excitation error zero
  gg = (/ 0,0,0 /)
  call AddReflection( gg )   ! this guarantees that 000 is always the first reflection

! now compute |sg|/|U_g|/lambda for the other allowed reflections; if this parameter is less than
! the threshhold, rBethe_i, then add the reflection to the list of potential reflections
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
        if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
         gg = (/ ix, iy, iz /)
         dval = 1.0/CalcLength( float(gg), 'r' )

         if ((IsGAllowed(gg)).and.(dval.gt.dmin)) then ! allowed by the lattice centering, if any
          sgp = Calcsg(float(gg),k,DynFN)
          if (cell%dbdiff(ix, iy, iz)) then ! potential double diffraction reflection
            if (abs(sgp).le.rBethe_d) then 
              call AddReflection( gg )
              rltail%sg = sgp
              rltail%dbdiff = .TRUE.
            end if
          else
            r_g = la * abs(sgp)/cabs(cell%LUT(ix, iy, iz))
            if (r_g.le.rBethe_i) then 
              call AddReflection( gg )
              rltail%sg = sgp
              rltail%dbdiff = .FALSE.
            end if
          end if
         end if ! IsGAllowed
        end if
       end do izl
      end do iyl
    end do ixl
    
  cell%DynNbeamsLinked = cell%DynNbeams

  if (present(verbose)) then 
    if (verbose) then 
      io_int(1) = cell%DynNbeams
      call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
    end if
  end if
  
end subroutine Initialize_ReflectionList




end module initializers
