! ###################################################################
! Copyright (c) 2013, Marc De Graef/Carnegie Mellon University
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
! CTEMsoft2013:gvectors.f90
!--------------------------------------------------------------------------
!
! MODULE: gvectors
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief variables and types needed to determine lists of reciprocal lattice vectors.  
!
!> @details This was part of the dynamical module, but was moved into a separate
!> module.  The new module includes a routine to delete the linked list, and also routines 
!> to allocate linked lists, which are used by almost all dynamical scattering codes.
!>
!> Due to some complicated module interdependencies the CalcDynMat routine is in
!> this module rather than in diffraction, where it would logically belong.  We may need
!> to figure out how to change that. 
! 
!> @date   04/29/13 MDG 1.0 original
!--------------------------------------------------------------------------
module gvectors

use local

IMPLICIT NONE

! linked list of reflections
type reflisttype  
  integer(kind=irg)          	:: num, &  		! sequential number
                              	hkl(3),&		! Miller indices
				HOLZN,& 		! belongs to this HOLZ layer
				famnum,&		! family number
				nab(2)  		! decomposition with respect to ga and gb
  character(1)               	:: ForB    		! 'F' for full, 'B' for Bethe
  logical                    	:: dbdiff  		! double diffraction reflection ?
  complex(kind=dbl)		:: Ucg,&   		! potential coefficient
                               amp     		! amplitude of beam
  real(kind=dbl)             	:: sg,xg,Ucgmod,sangle   ! excitation error, extinction distance, and modulus of potential coefficient; scattering angle (mrad)
  type(reflisttype),pointer 	:: next    		! connection to next entry
end type reflisttype

type(reflisttype),pointer 	:: reflist, & 	! linked list of reflections
                                rltail, &  		! end of linked list
                                rltmpa,rltmpb 	! temporary pointers

complex(kind=dbl),allocatable 	:: LUT(:,:,:)
logical,allocatable		:: dbdiff(:,:,:)

! define the cutoff parameters for the Bethe potential approach (and set to zero initially)
type BetheParameterType
	real(kind=sgl)			:: weakcutoff = 0.0_sgl
	real(kind=sgl)			:: cutoff = 0.0_sgl
	real(kind=sgl)			:: sgcutoff = 0.0_sgl
	integer(kind=irg)		:: nns
	integer(kind=irg)		:: nnw
	integer(kind=irg)		:: minweak
	integer(kind=irg)		:: minstrong
	integer(kind=irg)		:: maxweak
	integer(kind=irg)		:: maxstrong
	integer(kind=irg)		:: totweak
	integer(kind=irg)		:: totstrong
	integer(kind=irg),allocatable 	:: weaklist(:) 
	integer(kind=irg),allocatable 	:: stronglist(:)
	integer(kind=irg),allocatable 	:: weakhkl(:,:)
	integer(kind=irg),allocatable 	:: stronghkl(:,:)
	real(kind=sgl),allocatable	:: weaksg(:)
	real(kind=sgl),allocatable	:: strongsg(:)
	integer(kind=sgl),allocatable	:: reflistindex(:)		! used to map strong reflections onto the original reflist
end type BetheParameterType

type(BetheParameterType)	:: BetheParameter

! finally, create allocatable arrays with the same fields as the rlp list
! these may be needed for multithreaded applications, since linked lists
! cause OpenMP runtime problems.  This may not be a real problem but caused
! by me misunderstanding some aspects of OpenMP ... it's happened before ...
integer(kind=irg),allocatable	:: Reflist_hkl(:,:)
complex(kind=dbl),allocatable	:: Reflist_Ucg(:)
real(kind=sgl),allocatable	:: Reflist_sg(:), Reflist_xg(:)

! same for the BetheParameter arrays; somehow OpenMP doesn't like this type of variable...
integer(kind=irg),allocatable	:: BPweaklist(:), BPstronglist(:), BPweakhkl(:,:), BPstronghkl(:,:), BPreflistindex(:)
real(kind=sgl),allocatable	:: BPweaksg(:), BPstrongsg(:)
integer(kind=irg)		:: BPnns, BPnnw


contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: MakeRefList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief allocate and initialize the linked reflection list
!
!> @date   10/20/98 MDG 1.0 original
!> @date    5/22/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date  03/26/13  MDG 3.0 updated IO
!--------------------------------------------------------------------------
subroutine MakeRefList

use local
use error
use dynamical

IMPLICIT NONE

integer(kind=irg)  :: istat

! create it if it does not already exist
if (.not.associated(reflist)) then
  DynNbeams = 0
  allocate(reflist,stat=istat)
  if (istat.ne.0) call FatalError('MakeRefList:',' unable to allocate pointer')
  rltail => reflist                 	! tail points to new value
  nullify(rltail%next)              	! nullify next in new value
end if

end subroutine MakeRefList


!--------------------------------------------------------------------------
!
! SUBROUTINE: AddReflection
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief add a reflection to the linked reflection list
!
!> @param hkl Miller indices
!
!> @date   10/20/98 MDG 1.0 original
!> @date    5/22/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date  03/26/13  MDG  3.0 updated IO
!--------------------------------------------------------------------------
subroutine AddReflection(hkl)

use local
use error
use dynamical
use diffraction

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: hkl(3)		!< Miller indices of reflection to be added to list
integer(kind=irg)  		:: istat

! create reflist if it does not already exist
 if (.not.associated(reflist)) call MakeRefList

! create a new entry
 allocate(rltail%next,stat=istat)  		! allocate new value
 if (istat.ne.0) call FatalError('AddReflection',' unable to add new reflection')

 rltail => rltail%next             		! tail points to new value
 nullify(rltail%next)              		! nullify next in new value

 DynNbeams = DynNbeams + 1         		! update reflection counter
 rltail%num = DynNbeams            		! store reflection number
 rltail%hkl = hkl                  		! store Miller indices
 call CalcUcg(hkl)                 		! compute potential Fourier coefficient
 rltail%Ucg = rlp%Ucg              		! and store it in the list
 rltail%xg = rlp%xg                		! also store the extinction distance
 rltail%famnum = 0				! init this value for Prune_ReflectionList
 rltail%Ucgmod = cabs(rlp%Ucg)   		! added on 2/29/2012 for Bethe potential computations
 rltail%sangle = 1000.0*dble(CalcDiffAngle(hkl(1),hkl(2),hkl(3)))    ! added 4/18/2012 for EIC project HAADF/BF tomography simulations

end subroutine AddReflection

!--------------------------------------------------------------------------
!
! SUBROUTINE: RankReflections
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief rank reflections
!
!> @param k 
!> @param ga 
!> @param gb
!> @param fcnt 
!> @param srza 
!> @param iorder 
!
!> @date   10/20/98 MDG 1.0 original
!> @date    5/22/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date  03/26/13  MDG  3.0 updated IO
!--------------------------------------------------------------------------
subroutine RankReflections(k,ga,gb,fcnt,srza,iorder)

use local
use multibeams
use diffraction
use symmetryvars
use symmetry
use io
use crystal
use error
use dynamical

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: k(3)	!< incident beam direction
integer(kind=irg),INTENT(IN)	:: ga(3)	!< first pattern vector
integer(kind=irg),INTENT(IN)	:: gb(3)	!< second pattern vector
integer(kind=irg),INTENT(OUT)	:: fcnt	!< 
character(2),INTENT(IN)	:: srza      !< = SR for systematic row and ZA for zone axis orientation
integer(kind=irg),INTENT(OUT)	:: iorder	!< 


character(*),parameter 	:: var1 = ''
integer(kind=irg)   		:: i,j,ii,jj, &    ! various counters
                      		gn(3), & ! reciprocal lattice vectors
                      		il(48), gg(3),rfcnt,inm,imm,num,jcnt,iproj,jproj,ier, io_int(2)
real(kind=sgl)      		:: kk(3), &
                      		M(2,2),X(2),D ! variables used to decompose reciprocal lattice vectors
logical             		:: a,sgg 
logical,allocatable 		:: z(:,:)


 kk = float(k)
 iorder = SG % SYM_NUMpt

 if (srza.eq.'ZA') then 
  mess = ' The program will list all independent families of reflections of the zone'; call Message("(/A)")
  mess = ' that can be written as linear combinations of the two independent reflections'; call Message("(/A)")
  call ReadValue(' Enter the maximum multiple of ga and gb that should be included [integers] : ', io_int, 2)
  imm = io_int(1)
  inm = io_int(2)
 else
  mess = ' The program will list all independent families of reflections of the systematic row.'; call Message("(/A)")
  call ReadValue(' Enter the maximum multiple of ga that should be included [I] : ', io_int, 1)
  inm = io_int(1)
 end if

 fcnt = 0

! systematic row or zone axis ?
 select case (srza)

 case('ZA');

! use a logical array z to exclude family members once one
! member of the family has been evaluated
  if (.not.allocated(z))  allocate(z(-inm:inm,-imm:imm),stat=ier)
  if (ier.ne.0) call FatalError('RankReflections',' unable to allocate memory for array z')
  z(-inm:inm,-imm:imm) = .FALSE.

  do i=-inm,inm
   do j=-imm,imm
    if (.not.z(i,j)) then
     fcnt = fcnt + 1
     if (fcnt.gt.numr) call FatalError('RankReflections ',' too many families ')
     gn = i*ga+j*gb
     call CalcFamily(gn,num,'r')
     call GetOrder(kk,il,num,jcnt)
     numfam(fcnt) = jcnt
    do jj = 1,jcnt
      do ii = 1,3
       gg(ii) = itmp(il(jj),ii)
      end do
      if (jj.eq.1) then
       glen(fcnt) = CalcLength(float(gg),'r')
      end if

! determine the components of all family members with
! respect to ga and gb      
      M(1,1) = CalcDot(float(gb),float(gb),'c')
      M(1,2) = -CalcDot(float(ga),float(gb),'c')
      M(2,1) = M(1,2)
      M(2,2) = CalcDot(float(ga),float(ga),'c')
      D = M(1,1)*M(2,2) - M(1,2)*M(2,1)
      X(1) = CalcDot(float(gg),float(ga),'c')
      X(2) = CalcDot(float(gg),float(gb),'c')
      X = matmul(M,X)/D
      iproj = int(X(1))
      jproj = int(X(2))
      if (((iproj.ge.-inm).and.(iproj.le.inm)).and.((jproj.ge.-imm).and.(jproj.le.imm))) then
       z(iproj,jproj) = .TRUE.
      endif
      do ii = 1,3
       family(fcnt,jj,ii) = itmp(il(jj),ii)
      end do
     end do
    end if
   end do
  end do
  deallocate(z)

 case('SR');   ! systematic row case
! In this case we should mostly check for centrosymmetry, since 
! the non-centrosymmetric case may give +g and -g reflections with
! different intensities.
! The most straightforward way to check this is to determine 
! whether or not -g belongs to the family {+g}.  If it does, then
! the systematic row will be symmetric in +-g, else it is not.
  call CalcFamily(ga,num,'r')
  sgg = .FALSE.
  do ii=1,num
   gg(1:3) = itmp(ii,1:3)
   if (sum(gg+ga).eq.0) then 
     sgg = .TRUE.
   end if
  end do
  fcnt = fcnt + 1
  numfam(fcnt) = 1
  glen(fcnt) = 0.0
  family(fcnt,1,1:3) = (/0,0,0/)

  if (sgg) then
! -g does belong to the same family; enumerate all families
   iorder = 2
   do i=1,inm
    fcnt = fcnt + 1
    numfam(fcnt) = 2
    glen(fcnt) = CalcLength(float(i*ga),'r')
    family(fcnt,1,1:3) = i*ga(1:3)
    family(fcnt,2,1:3) = -i*ga(1:3)
   end do
  else
! -g does not belong to the same family
   iorder = 1
   do i=-inm,inm
    fcnt = fcnt + 1
    numfam(fcnt) = 1
    glen(fcnt) = CalcLength(float(i*ga),'r')
    family(fcnt,1,1:3) = i*ga(1:3)
   end do
  end if

 end select   ! systematic row or zone axis

! allocate the variables defined in module multibeams
 if (.not.allocated(gm)) allocate(gm(fcnt),stat=ier)
 if (ier.ne.0) call FatalError(var1,'gm')
 if (.not.allocated(idx)) allocate(idx(fcnt),stat=ier)
 if (ier.ne.0) call FatalError(var1,'idx')
 if (.not.allocated(al)) allocate(al(fcnt),stat=ier)
 if (ier.ne.0) call FatalError(var1,'al')
 if (.not.allocated(V)) allocate(V(fcnt,4),stat=ier)
 if (ier.ne.0) call FatalError(var1,'V')
 gm=0.0
 gm(1:fcnt) = glen(1:fcnt)

! next, rank the families by increasing |g|
 call SPSORT(gm,fcnt,idx,1,ier)

! and print them out in the correct order, removing the ones that are not
! allowed by lattice centering !  Reflections with zero structure factor,
! but not forbidden by lattice centering may give rise to double diffracted
! beams, so they can not be excluded at this point.
 mess = ' List of independent reflections, ranked by increasing |g|'; call Message("(/A/)")
 rfcnt = 0
 do i=1,fcnt
  gn(1:3) = family(idx(i),1,1:3)
  a = IsGAllowed(gn)
  if (a) then
   call CalcUcg(gn)
   rfcnt=rfcnt+numfam(idx(i))
   al(idx(i)) = .TRUE.
   V(idx(i),1) = rlp%Vmod
   V(idx(i),2) = rlp%Vphase
   V(idx(i),3) = rlp%Vpmod
   V(idx(i),4) = rlp%Vpphase
    write (stdout,"('(',I3,I3,I3,'); |g| = ',F7.3,'; # = ',I2,'; nref = ',I3,'; |Vg| = ',F8.5,' Volt')") &
        (gn(j),j=1,3),gm(idx(i)),numfam(idx(i)),rfcnt,rlp%Vmod
  else
   al(idx(i)) = .FALSE.
  end if
 end do
 
end subroutine RankReflections

!--------------------------------------------------------------------------
!
! SUBROUTINE: SelectReflections
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determine a subset of beams to be included in multi-beam
!> simulation
!
!> @param fcnt
!> @param rfcnt
!> @param ccnt
!
!> @date    5/24/01 MDG 1.0 f90
!> @date   11/27/01 MDG 1.1 added kind support
!> @date  03/26/13  MDG 2.0 updated IO
!--------------------------------------------------------------------------
subroutine SelectReflections(fcnt,rfcnt,ccnt)

use local
use multibeams
use crystalvars
use symmetryvars
use symmetry
use io
use dynamical

IMPLICIT NONE

integer(kind=irg)    	:: rfcnt, gn(3), fcnt,i,j,igv, ccnt, io_int(1)
real(kind=sgl)       	:: gmax, io_real(1)

intent(IN) 		:: fcnt
intent(OUT)		:: rfcnt,ccnt

! ask for maximum value in terms of |g| or in terms of Vmod
 mess = ' The reciprocal lattice vectors contributing to the computation'; call Message("(/A)")
 mess = ' must now be selected.  You can use a maximum |g| value to '; call Message("(A)")
 mess = ' truncate reciprocal space, or you can use a |V| threshold.'; call Message("(A/)")
! is this a non-symmorphic space group ? If so, add a little note
 if (minval(abs(SGsym - cell % SYM_SGnum)).ne.0) then
  mess = ' The list above may show reflections with zero structure factor because'; call Message("(A)")
  mess = ' of non-symmorphic space group symmetry elements.  Those reflections'; call Message("(A)")
  mess = ' can be incorporated in the simulation only with the maximum |g| criterion.'; call Message("(A/)")
 end if
 call ReadValue('Use maximum |g| as criterion (1) or threshold |V| (2) :', io_int, 1)
 igv = io_int(1)

 if (igv.eq.1) then
  mess = ' You have selected the maximum |g| criterion.'; call Message("(/A)")
  call ReadValue('Enter the truncation value for |g| (strictly smaller) : ', io_real, 1)
  gmax = io_real(1)
  do i=1,fcnt
   if ((gm(idx(i)).gt.gmax).and.(al(idx(i)))) al(idx(i)) = .FALSE.
  end do 
 else
  mess = ' You have selected the threshold |V| criterion.'; call Message("(/A)")
  call ReadValue('Enter the minimum value for |V| (strictly larger) : ', io_real, 1)
  gmax = io_real(1)
  do i=1,fcnt
   if ((V(idx(i),1).lt.gmax).and.(al(idx(i)))) al(idx(i)) = .FALSE.
  end do 
 endif

 rfcnt = 0
 ccnt = 0
 mess = ' The following reflections will be included :'; call Message("(A/)")
 do i=1,fcnt
  if (al(idx(i))) then 
! place all reflections in the linked reflection list
   do j=1,numfam(idx(i))
    call AddReflection(family(idx(i),j,1:3))
    rltail%famnum = ccnt+1
   end do
! and print list of included families
   gn(1:3) = family(idx(i),1,1:3)
   ccnt = ccnt + 1
   rfcnt=rfcnt+numfam(idx(i))
    write (stdout,"('(',I3,I3,I3,'); |g| = ',F7.3,'; # = ',I2,'; nref = ',I3,'; |Vg| = ',F8.5,' Volt')") &
        (gn(j),j=1,3),gm(idx(i)),numfam(idx(i)),rfcnt,V(idx(i),1)
  end if
 end do

 io_int(1) = rfcnt
 call WriteValue('Total number of reflections contributing to computation : ', io_int, 1) 

end subroutine SelectReflections


!--------------------------------------------------------------------------
!
! SUBROUTINE: Printrlp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief output the contents of the rlp structure
!
!> @param first logical switch to provide long form output (optional)
!
!> @date   10/20/98 MDG 1.0 original
!> @date    5/22/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date  03/26/13  MDG  3.0 updated IO
!--------------------------------------------------------------------------
subroutine Printrlp(first)

use local
use io
use constants
use dynamical

IMPLICIT NONE

logical,optional,intent(INOUT) 	:: first		!< switch for long/short output
integer(kind=irg)			:: oi_int(3)
real(kind=sgl)				:: oi_real(7)
complex(kind=sgl)			:: oi_cmplx(1)

if (present(first)) then
 if (first) then
  mess = '     Scattering factors : '; call Message("(/A,$)")
  
  if (rlp%method.eq.'WK') then 
   if (rlp%absorption.eqv..TRUE.) then 
    mess = ' Weickenmeier-Kohl (with absorption)'; call Message("(A/)")
   else
    mess = ' Weickenmeier-Kohl'; call Message("(A/)")
   end if
  else
    mess = ' Doyle-Turner/Smith-Burge'; call Message("(A/)")
  end if

  if (rlp%absorption.eqv..TRUE.) then
    mess = '   h  k  l    |g|    Ucg_r  Ucg_i   |Ug|    phase   |Ugp|   phase   xi_g   xi_gp    ratio  Re-1/q_g-Im'
    call Message("(A)") 
  else
    mess = '   h  k  l    |g|    Ucg_r  |Ug|    phase    xi_g   1/q_g'
    call Message("(A)") 
  end if
  first = .FALSE.
 end if
end if

if (rlp%absorption.eqv..TRUE.) then
 oi_int(1:3) = rlp%hkl(1:3)
 call WriteValue('',oi_int, 3, "(1x,3I3,1x,$)")
 oi_real(1) = rlp%g
 call WriteValue('',oi_real, 1, "(F9.4,$)")
 oi_cmplx(1) = rlp%Ucg
 call WriteValue('',oi_cmplx, 1, "(2F7.3,1x,$)")
 oi_real(1:7)  = (/ rlp%Umod,rlp%Vphase*180.0/sngl(cPi),rlp%Upmod,rlp%Vpphase*180.0/sngl(cPi),rlp%xg,rlp%xgp,rlp%ar /)
 call WriteValue('',oi_real, 7, "(4F8.3,3F8.1,$)")
 oi_cmplx(1) = rlp%qg
 call WriteValue('',oi_cmplx, 1, "(2F8.3)")
else
 oi_int(1:3) = rlp%hkl(1:3)
 call WriteValue('',oi_int, 3, "(1x,3I3,1x,$)")
 oi_real(1) = rlp%g
 call WriteValue('',oi_real, 1, "(F9.4,$)")
 oi_real(1) = real(rlp%Ucg)
 call WriteValue('',oi_real, 1, "(F7.3,1x,$)")
 oi_real(1:3)  = (/ rlp%Umod,rlp%Vphase*180.0/sngl(cPi),rlp%xg /)
 call WriteValue('',oi_real, 3, "(2F8.3,F8.1,$)")
 oi_cmplx(1) = rlp%qg
 call WriteValue('',oi_cmplx, 1, "(2F8.3)")
end if

end subroutine Printrlp

!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcDynMat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the dynamical matrix
!
!> @details compute the dynamical matrix for either Bloch wave
!> or Darwin-Howie-Whelan approach. Bethe potentials
!> are not yet implemented in this routine, but should
!> be at some point in the future...
!
!> @param calcmode computational mode (D-H-W, BLOCH, DIAGH, DIAGB, BETHE)
!
!> @note This is the old routine; should be replaced everywhere by a call to 
!> Compute_DynMat
!
!> @todo Add support for Bethe potentials; probably needs a complete rewrite to include all cases
!
!> @date   10/20/98 MDG 1.0 original
!> @date    5/22/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date  03/26/13  MDG  3.0 updated IO
!--------------------------------------------------------------------------
subroutine CalcDynMat(calcmode)

use local
use constants
use error
use crystal
use dynamical
use diffraction

IMPLICIT NONE

character(5),INTENT(IN)       	:: calcmode  !< computational mode  (D-H-W, BLOCH, DIAGH, DIAGB, or BETHE)
integer(kind=irg) 		:: istat,ir,ic
real(kind=sgl)     		:: glen,exer
complex(kind=dbl)  		:: czero,pre

! has reflist been allocated ?
if (.not.associated(reflist)) call FatalError('CalcDynMat',' reflection list has not been allocated')

! initialize some parameters
czero = cmplx(0.0,0.0,dbl)
pre = cmplx(0.0,cPi,dbl)

! allocate DynMat if it hasn't already been allocated
if (.not.allocated(DynMat)) then
  allocate(DynMat(DynNbeams,DynNbeams),stat=istat)
  DynMat = czero
! get the absorption coefficient
  call CalcUcg((/0,0,0/))
  DynUpz = rlp%Vpmod
end if

! are we supposed to fill the off-diagonal part ?
 if ((calcmode.eq.'D-H-W').or.(calcmode.eq.'BLOCH')) then
  rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
  do ir=1,DynNbeams
   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
   do ic=1,DynNbeams
    if (ic.ne.ir) then  ! exclude the diagonal
! compute Fourier coefficient of electrostatic lattice potential 
     call CalcUcg(rltmpa%hkl - rltmpb%hkl)
     if (calcmode.eq.'D-H-W') then
      DynMat(ir,ic) = pre*rlp%qg
     else
      DynMat(ir,ic) = rlp%Ucg
     end if
    end if
    rltmpb => rltmpb%next  ! move to next column-entry
   end do
   rltmpa => rltmpa%next   ! move to next row-entry
  end do
 end if

! or the diagonal part ?
 if ((calcmode.eq.'DIAGH').or.(calcmode.eq.'DIAGB')) then
  rltmpa => reflist%next   ! point to the front of the list
! ir is the row index
  do ir=1,DynNbeams
   glen = CalcLength(float(rltmpa%hkl),'r')
   if (glen.eq.0.0) then
    DynMat(ir,ir) = cmplx(0.0,DynUpz,dbl)
   else  ! compute the excitation error
    exer = Calcsg(float(rltmpa%hkl),DynWV,DynFN)
    rltmpa%sg = exer
    if (calcmode.eq.'DIAGH') then  !
     DynMat(ir,ir) = cmplx(0.0,2.D0*cPi*exer,dbl)
    else
     DynMat(ir,ir) = cmplx(2.D0*exer/mLambda,DynUpz,dbl)
    end if
   endif
   rltmpa => rltmpa%next   ! move to next row-entry
  end do
 end if

! or are we making use of Bethe potentials ?
! This part is yet to be implemented

end subroutine CalcDynMat

!--------------------------------------------------------------------------
!
! SUBROUTINE: Prune_ReflectionList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief select from the reflection list those g-vectors that will be counted in an LACBED computation
!
!> @details This routine basicaly repeats a section form the Compute_DynMat routine
!> without actually computing the dynamical matrix; it simple keeps track of all the
!> beams that are at one point or another regarded as strong beams and write their 
!> information to the dataunit file (which must be opened in the calling program).
!> We'll use the famnum field in the rlp linked list to flag the strong reflections.
!
!> @param numk number of wave vectors to consider
!> @param nstrong total number of unique strong beams
!
!> @date   09/20/13 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine Prune_ReflectionList(numk,nstrong)

use local
use dynamical
use kvectors

IMPLICIT NONE

integer(kind=irg),INTENT(IN)		:: numk
integer(kind=irg),INTENT(OUT)		:: nstrong

integer(kind=irg)			:: ik, ig, istrong
real(kind=sgl)     			:: sgp, lUg, cut1, cut2
integer(kind=irg),allocatable		:: strongreflections(:,:)


! reset the value of DynNbeams in case it was modified in a previous call 
DynNbeams = DynNbeamsLinked
  	
nstrong = 0

! loop over all the incident beam directions
  ktmp => khead
! loop over all beam orientations, selecting them from the linked list
  do ik = 1,numk
  
! reset the reflection linked list
    rltmpa => reflist%next

! pick the first reflection since that is the transmitted beam (only on the first time)
    if (ik.eq.1) then 
      rltmpa%famnum = 1    
      nstrong = nstrong + 1
    end if

! loop over all reflections in the linked list    
    rltmpa => rltmpa%next
    reflectionloop: do ig=2,DynNbeamsLinked

! We compare the precomputed |sg| with two multiples of lambda |Ug|
!
!  |sg|>cutoff lambda |Ug|   ->  don't count reflection
!  cutoff lambda |Ug| > |sg| > weakcutoff lambda |Ug|  -> weak reflection
!  weakcutoff lambda |Ug| > |sg|  -> strong reflection
!
        sgp = abs(rltmpa%sg) 
        lUg = cabs(rltmpa%Ucg) * mLambda
        cut1 = BetheParameter%cutoff * lUg
        cut2 = BetheParameter%weakcutoff * lUg

! we have to deal separately with double diffraction reflections, since
! they have a zero potential coefficient !        
        if ( dbdiff(rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3)) ) then  ! it is a double diffraction reflection
          if ((sgp.le.BetheParameter%sgcutoff).and.(rltmpa%famnum.ne.1)) then         
	    nstrong = nstrong + 1
    	    rltmpa%famnum = 1    
          end if
        else   ! it is not a double diffraction reflection
          if (sgp.le.cut1) then  ! count this beam
! is this a weak or a strong reflection (in terms of Bethe potentials)? 
             	if ((sgp.le.cut2).and.(rltmpa%famnum.ne.1)) then ! it's a strong beam
	    		nstrong = nstrong + 1
		     	rltmpa%famnum = 1    
             	end if
          end if
        end if
! go to the next beam in the list
      rltmpa => rltmpa%next
    end do reflectionloop
! go to the next incident beam direction
    if (ik.ne.numk) ktmp => ktmp%next
  end do ! ik loop

! ok, now that we have the list, we'll go through it again to set sequential numbers instead of 1's
! at the same time, we'll extract the data that we need to write to the data file
  rltmpa => reflist%next
  rltmpa => rltmpa%next
  reflectionloop2: do ig=2,DynNbeamsLinked
    if (rltmpa%famnum.eq.1) then
	istrong = istrong + 1
      	rltmpa%famnum = istrong
    endif
! go to the next beam in the list
    rltmpa => rltmpa%next
  end do reflectionloop2

end subroutine Prune_ReflectionList


!--------------------------------------------------------------------------
!
! SUBROUTINE: Compute_DynMat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the entire dynamical matrix, including HOLZ and Bethe potentials
!
!> @details This is a complicated routine because it forms the core of the dynamical
!> scattering simulations.  The routine must be capable of setting up the dynamical 
!> matrix for systematic row and zone axis, with or without HOLZ reflections, and 
!> with Bethe potentials for the Bloch wave case.  The routine must also be able to
!> decide which reflections are weak, and which are strong (again for the Bloch wave
!> case, but potentially also for other cases? Further research needed...).
!
!> @param calcmode string that describes the particular matrix mode
!> @param kk incident wave vector
!> @param IgnoreFoilNormal switch for foil normal inclusion in sg computation
!> @param IncludeSecondOrder (optional) switch to include second order correction to Bethe potentials
!
!> @date   05/06/13 MDG 1.0 original
!> @date   08/30/13 MDG 1.1 correction of effective excitation error
!> @date   09/20/13 MDG 1.2 added second order Bethe potential correction switch
!--------------------------------------------------------------------------
recursive subroutine Compute_DynMat(calcmode,kk,IgnoreFoilNormal,IncludeSecondOrder)

use local
use dynamical
use error
use constants
use crystal
use diffraction
use io

IMPLICIT NONE

character(*),INTENT(IN)		:: calcmode		!< computation mode
real(kind=dbl),INTENT(IN)		:: kk(3)		!< incident wave vector
logical,INTENT(IN)			:: IgnoreFoilNormal	!< how to deal with the foil normal
logical,INTENT(IN),OPTIONAL		:: IncludeSecondOrder	!< second order Bethe potential correction switch

complex(kind=dbl)  			:: czero,pre, weaksum, ughp, uhph
integer(kind=irg) 		 	:: istat,ir,ic,nn, iweak, istrong, iw, ig, ll(3), gh(3), nnn
real(kind=sgl)     			:: glen,exer,gg(3), kpg(3), gplen, sgp, lUg, cut1, cut2
real(kind=dbl)				:: lsfour, weaksgsum 
logical					:: AddSecondOrder
 
AddSecondOrder = .FALSE.
if (present(IncludeSecondOrder)) AddSecondOrder = .TRUE.

! has the list of reflections been allocated ?
if (.not.associated(reflist)) call FatalError('Compute_DynMat',' reflection list has not been allocated')

! if the dynamical matrix has already been allocated, deallocate it first
! this is partially so that no program will allocate DynMat itself; it must be done
! via this routine only.
if (allocated(DynMat)) deallocate(DynMat)

! initialize some parameters
czero = cmplx(0.0,0.0,dbl)	! complex zero
pre = cmplx(0.0,cPi,dbl)		! i times pi

if (calcmode.ne.'BLOCHBETHE') then

! allocate DynMat
	  allocate(DynMat(DynNbeams,DynNbeams),stat=istat)
	  DynMat = czero
! get the absorption coefficient
	  call CalcUcg( (/0,0,0/) )
	  DynUpz = rlp%Vpmod

! are we supposed to fill the off-diagonal part ?
	 if ((calcmode.eq.'D-H-W').or.(calcmode.eq.'BLOCH')) then
	  rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
	  do ir=1,DynNbeams
	   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
	   do ic=1,DynNbeams
	    if (ic.ne.ir) then  ! exclude the diagonal
! compute Fourier coefficient of electrostatic lattice potential 
	     gh = rltmpa%hkl - rltmpb%hkl
	     if (calcmode.eq.'D-H-W') then
	      call CalcUcg(gh)
	      DynMat(ir,ic) = pre*rlp%qg
	     else
	      DynMat(ir,ic) = LUT(gh(1),gh(2),gh(3))
	     end if
	    end if
	    rltmpb => rltmpb%next  ! move to next column-entry
	   end do
	   rltmpa => rltmpa%next   ! move to next row-entry
	  end do
	 end if
	 
! or the diagonal part ?
	 if ((calcmode.eq.'DIAGH').or.(calcmode.eq.'DIAGB')) then
	  rltmpa => reflist%next   ! point to the front of the list
! ir is the row index
	  do ir=1,DynNbeams
	   glen = CalcLength(float(rltmpa%hkl),'r')
	   if (glen.eq.0.0) then
	    DynMat(ir,ir) = cmplx(0.0,DynUpz,dbl)
	   else  ! compute the excitation error
	    exer = Calcsg(float(rltmpa%hkl),sngl(kk),DynFN)
	    rltmpa%sg = exer
	    if (calcmode.eq.'DIAGH') then  !
	     DynMat(ir,ir) = cmplx(0.0,2.D0*cPi*exer,dbl)
	    else
	     DynMat(ir,ir) = cmplx(2.D0*exer/mLambda,DynUpz,dbl)
	    end if
	   endif
	   rltmpa => rltmpa%next   ! move to next row-entry
	  end do
	 end if

else  ! this is the Bloch wave + Bethe potentials initialization (originally implemented in the EBSD programs)

! we don't know yet how many strong reflections there are so we'll need to determine this first
! this number depends on some externally supplied parameters, which we will get from a namelist
! file (which should be read only once by the Set_Bethe_Parameters routine), or from default values
! if there is no namelist file in the folder.
	if (BetheParameter%cutoff.eq.0.0) call Set_Bethe_Parameters


! reset the value of DynNbeams in case it was modified in a previous call 
  	DynNbeams = DynNbeamsLinked
  	
! precompute lambda^2/4
	lsfour = mLambda**2*0.25D0
  
! first, for the input beam direction, determine the excitation errors of 
! all the reflections in the master list, and count the ones that are
! needed for the dynamical matrix (weak as well as strong)
        if (.not.allocated(BetheParameter%weaklist)) allocate(BetheParameter%weaklist(DynNbeams))
        if (.not.allocated(BetheParameter%stronglist)) allocate(BetheParameter%stronglist(DynNbeams))
        if (.not.allocated(BetheParameter%reflistindex)) allocate(BetheParameter%reflistindex(DynNbeams))

	BetheParameter%weaklist = 0
	BetheParameter%stronglist = 0
	BetheParameter%reflistindex = 0


    	rltmpa => reflist%next

! deal with the transmitted beam first
    nn = 1		! nn counts all the scattered beams that satisfy the cutoff condition
    nnn = 1		! nnn counts only the strong beams
    BetheParameter%stronglist(nn) = 1   ! make sure that the transmitted beam is always a strong beam ...
    BetheParameter%weaklist(nn) = 0
    BetheParameter%reflistindex(nn) = 1

    rltmpa%sg = 0.D0    

! loop over all reflections in the linked list    
    rltmpa => rltmpa%next
    reflectionloop: do ig=2,DynNbeamsLinked
      gg = float(rltmpa%hkl)        		! this is the reciprocal lattice vector 

! deal with the foil normal; if IgnoreFoilNormal is .TRUE., then assume it is parallel to the beam direction
     if (IgnoreFoilNormal) then 
! we're taking the foil normal to be parallel to the incident beam direction at each point of
! the standard stereographic triangle, so cos(alpha) = 1 always in eqn. 5.11 of CTEM
        kpg = kk+gg                		! k0 + g (vectors)
        gplen = CalcLength(kpg,'r')  	! |k0+g|
        rltmpa%sg = (1.0/mLambda**2 - gplen**2)*0.5/gplen
     else
	rltmpa%sg = Calcsg(gg,sngl(kk),DynFN)
     end if

! use the reflection num entry to indicate whether or not this
! reflection should be used for the dynamical matrix
! We compare |sg| with two multiples of lambda |Ug|
!
!  |sg|>cutoff lambda |Ug|   ->  don't count reflection
!  cutoff lambda |Ug| > |sg| > weakcutoff lambda |Ug|  -> weak reflection
!  weakcutoff lambda |Ug| > |sg|  -> strong reflection
!
        sgp = abs(rltmpa%sg) 
        lUg = cabs(rltmpa%Ucg) * mLambda
        cut1 = BetheParameter%cutoff * lUg
        cut2 = BetheParameter%weakcutoff * lUg

! we have to deal separately with double diffraction reflections, since
! they have a zero potential coefficient !        
        if ( dbdiff(rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3)) ) then  ! it is a double diffraction reflection
          if (sgp.le.BetheParameter%sgcutoff) then         
            	nn = nn+1
            	nnn = nnn+1
		BetheParameter%stronglist(ig) = 1
            	BetheParameter%reflistindex(ig) = nnn
          end if
        else   ! it is not a double diffraction reflection
          if (sgp.le.cut1) then  ! count this beam
             	nn = nn+1
! is this a weak or a strong reflection (in terms of Bethe potentials)? 
             	if (sgp.le.cut2) then ! it's a strong beam
              		nnn = nnn+1
	      		BetheParameter%stronglist(ig) = 1
              		BetheParameter%reflistindex(ig) = nnn
             	else ! it's a weak beam
              		BetheParameter%weaklist(ig) = 1
             	end if
          end if
        end if
! go to the next beam in the list
      rltmpa => rltmpa%next
    end do reflectionloop

! if we don't have any beams in this list (unlikely, but possible if the cutoff and 
! weakcutoff parameters have unreasonable values) then we abort the run
! and we report some numbers to the user 
	 if (nn.eq.0) then
	   mess = ' no beams found for the following parameters:'; call Message("(A)")
	   write (stdout,*) ' wave vector = ',kk,'  -> number of beams = ',nn
	   mess =  '   -> check cutoff and weakcutoff parameters for reasonableness'; call Message("(A)")
	   call FatalError('Compute_DynMat','No beams in list')
	end if

! next, we define nns to be the number of strong beams, and nnw the number of weak beams.
	 BetheParameter%nns = sum(BetheParameter%stronglist)
	 BetheParameter%nnw = sum(BetheParameter%weaklist)

! We may want to keep track of the total and average numbers of strong and weak beams  
	 BetheParameter%totweak = BetheParameter%totweak + BetheParameter%nnw
	 BetheParameter%totstrong = BetheParameter%totstrong + BetheParameter%nns
	 if (BetheParameter%nnw.lt.BetheParameter%minweak) BetheParameter%minweak=BetheParameter%nnw
	 if (BetheParameter%nnw.gt.BetheParameter%maxweak) BetheParameter%maxweak=BetheParameter%nnw
	 if (BetheParameter%nns.lt.BetheParameter%minstrong) BetheParameter%minstrong=BetheParameter%nns
	 if (BetheParameter%nns.gt.BetheParameter%maxstrong) BetheParameter%maxstrong=BetheParameter%nns

! allocate arrays for weak and strong beam information
	if (allocated(BetheParameter%weakhkl)) deallocate(BetheParameter%weakhkl)
	if (allocated(BetheParameter%weaksg)) deallocate(BetheParameter%weaksg)
	if (allocated(BetheParameter%stronghkl)) deallocate(BetheParameter%stronghkl)
	if (allocated(BetheParameter%strongsg)) deallocate(BetheParameter%strongsg)
	allocate(BetheParameter%weakhkl(3,BetheParameter%nnw),BetheParameter%weaksg(BetheParameter%nnw))
	allocate(BetheParameter%stronghkl(3,BetheParameter%nns),BetheParameter%strongsg(BetheParameter%nns))

! here's where we extract the relevant information from the linked list (much faster
! than traversing the list each time...)
	rltmpa => reflist%next    ! reset the a list
	iweak = 0
	istrong = 0
	do ir=1,DynNbeamsLinked
	     if (BetheParameter%weaklist(ir).eq.1) then
	        iweak = iweak+1
	        BetheParameter%weakhkl(1:3,iweak) = rltmpa%hkl(1:3)
	        BetheParameter%weaksg(iweak) = rltmpa%sg
	     end if
	     if (BetheParameter%stronglist(ir).eq.1) then
	        istrong = istrong+1
	        BetheParameter%stronghkl(1:3,istrong) = rltmpa%hkl(1:3)
	        BetheParameter%strongsg(istrong) = rltmpa%sg
	     end if
	   rltmpa => rltmpa%next
	end do

! now we are ready to create the dynamical matrix
	DynNbeams = BetheParameter%nns

! allocate DynMat if it hasn't already been allocated and set to complex zero
	  if (allocated(DynMat)) deallocate(DynMat)
	  allocate(DynMat(DynNbeams,DynNbeams),stat=istat)
	  DynMat = czero

! get the absorption coefficient
	  call CalcUcg( (/0,0,0/) )
	  DynUpz = rlp%Vpmod

! ir is the row index
       do ir=1,BetheParameter%nns
! ic is the column index
          do ic=1,BetheParameter%nns
! compute the Bethe Fourier coefficient of the electrostatic lattice potential 
              if (ic.ne.ir) then  ! not a diagonal entry
                 ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
                 DynMat(ir,ic) = LUT(ll(1),ll(2),ll(3)) 
! and subtract from this the total contribution of the weak beams
                 weaksum = czero
                 do iw=1,BetheParameter%nnw
                      ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%weakhkl(1:3,iw)
                      ughp = LUT(ll(1),ll(2),ll(3)) 
                      ll = BetheParameter%weakhkl(1:3,iw) - BetheParameter%stronghkl(1:3,ic)
                      uhph = LUT(ll(1),ll(2),ll(3)) 
                      weaksum = weaksum +  ughp * uhph *cmplx(1.D0/BetheParameter%weaksg(iw),0.0,dbl)
                 end do
! and correct the dynamical matrix element to become a Bethe potential coefficient
                 DynMat(ir,ic) = DynMat(ir,ic) - cmplx(0.5D0*mLambda,0.0D0,dbl)*weaksum
! do we need to add the second order corrections ?
		  if (AddSecondOrder) then 
		    weaksum = czero
		  end if
              else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! determine the total contribution of the weak beams
                 weaksgsum = 0.D0
		  do iw=1,BetheParameter%nnw
                      ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%weakhkl(1:3,iw)
                      ughp = LUT(ll(1),ll(2),ll(3)) 
                      weaksgsum = weaksgsum +  cabs(ughp)**2/BetheParameter%weaksg(iw)
                 end do
                 weaksgsum = weaksgsum * mLambda/2.D0
                 DynMat(ir,ir) = cmplx(2.D0*BetheParameter%strongsg(ir)/mLambda-weaksgsum,DynUpz,dbl)
! do we need to add the second order corrections ?
		  if (AddSecondOrder) then 
		    weaksum = czero
		  end if
               end if
          end do
        end do
! that should do it for the initialization of the dynamical matrix

end if	 ! Bethe potential initialization

end subroutine Compute_DynMat

!--------------------------------------------------------------------------
!
! SUBROUTINE: Set_Bethe_Parameters
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Read the Bethe potential parameters from a file, if it exists; otherwise take defaults
!
!> @details The parameters set in this routine determine the difference between strong and
!> weak reflections.  The routine checks for the presence of the BetheParameters.nml file
!> in the current folder.  If present, it will read the parameters, otherwise it will use 
!> defaults which have been determined to be reasonable based on dynamical EBSD runs. 
!
!> @date   05/08/13 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine Set_Bethe_Parameters()

use local
use io

IMPLICIT NONE

character(fnlen),parameter 	:: Bethefilename = 'BetheParameters.nml'
logical				:: fexist
real(kind=sgl)			:: weakcutoff, cutoff, sgcutoff

namelist /BetheList/ weakcutoff, cutoff, sgcutoff

! check for the presence of the namelist file in the current folder
inquire(file=trim(Bethefilename),exist=fexist)

! set all default values (must be done here, since nml file may not contain all of them)
weakcutoff = 50.0  	! dimensionless cutoff parameter, smaller = strong, larger = weak
cutoff = 100.0		! overall cutoff parameter
sgcutoff = 5.0		! sg cutoff for double diffraction reflections

if (fexist.eq..TRUE.) then ! check for the file in the local folder
! read the parameters from the file
 OPEN(UNIT=dataunit,FILE=trim(Bethefilename),DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=BetheList)
 CLOSE(UNIT=dataunit)
 mess = 'Read Bethe parameters from BetheParameters.nml'; call Message("(A)")
 write (stdout,nml=BetheList)
end if

BetheParameter % weakcutoff = weakcutoff
BetheParameter % cutoff = cutoff
BetheParameter % sgcutoff = sgcutoff

end subroutine Set_Bethe_Parameters


!--------------------------------------------------------------------------
!
! SUBROUTINE: Compute_ReflectionList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the entire reflection list for general conditions (including HOLZ)
!
!> @details also computes the LUT (LookUpTable) that stores all the scattering
!> potential coefficients that are needed to fill the dynamical matrix (only for the 
!> general case).
!
!> @param dmin minimum d-spacing to allow in the list
!> @param k zone axis indices
!> @param ga first reciprocal vector of zone
!> @param gb second reciprocal vector 
!> @param method  approach to follow (ALL or ZA)
!> @param ConvertList logical, determines whether or not conversion of list is needed
!> @param maxholz maximum range of HOLZ reflections to include
!
!> @date 04/29/13 MDG 1.0 original
!> @date 09/20/13 MDG 1.1 corrected handling of LUT
!--------------------------------------------------------------------------
subroutine Compute_ReflectionList(dmin,k,ga,gb,method,ConvertList,maxholz)

use local
use io
use crystal
use dynamical
use diffraction
use symmetry

IMPLICIT NONE

real(kind=sgl),INTENT(IN)				:: dmin
integer(kind=irg),INTENT(IN)				:: k(3)
integer(kind=irg),INTENT(IN)				:: ga(3)
integer(kind=irg),INTENT(IN)				:: gb(3)
character(*),INTENT(IN)				:: method
logical,INTENT(IN)					:: ConvertList
integer(kind=irg),INTENT(IN),OPTIONAL			:: maxholz

integer(kind=irg)					:: imh, imk, iml, gg(3), ix, iy, iz, i, minholz, RHOLZ, im, istat, N, ig
real(kind=sgl)						:: dhkl, io_real(9), H, g3(3), g3n(3), FNg(3), ddt
integer(kind=irg)					:: io_int(3), gshort(3), gp(3)

! determine the master list of reflections for the general case (i.e., a box in reciprocal space)
if (method.eq.'ALL') then 
! set threshold for double diffraction detection
  ddt = 1.0e-6

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

! the LUT array stores all the Fourier coefficients, so that we only need to compute them once...
  allocate(LUT(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
  LUT = dcmplx(0.D0,0.D0)

! allocate an array that keeps track of potential double diffraction reflections
  allocate(dbdiff(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
  dbdiff = .FALSE.
  
  DynNbeams = 1
  gg = (/ 0,0,0 /)
  call AddReflection( gg )   ! this guarantees that 000 is always the first reflection
  call CalcUcg(gg)   
  DynUpz = rlp%Vpmod
  io_real(1) = rlp%xgp
  call WriteValue(' Normal absorption length = ', io_real, 1)

! and add this reflection to the look-up table
  LUT(gg(1),gg(2),gg(3)) = rlp%Ucg

! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections,
! since the dynamical matrix uses g-h as its index !!!  However, the linked list of reflections
! should only contain the g, h refelctions separately, not their differences. (i.e., smaller box)
    do ix=-2*imh,2*imh
      do iy=-2*imk,2*imk
       do iz=-2*iml,2*iml
        gg = (/ ix, iy, iz /)
        if ((IsGAllowed(gg)).AND.(sum(abs(gg)).ne.0)) then
! if this g is inside the original box, then add it to the linked list
          if ((abs(ix).le.imh).and.(abs(iy).le.imk).and.(abs(iz).le.iml)) then 
            call AddReflection(gg)
          else
            call CalcUcg(gg)
          end if
! add the reflection to the look up table
          LUT(ix, iy, iz) = rlp%Ucg
! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
          if (cabs(rlp%Ucg).le.ddt) then 
            dbdiff(ix,iy,iz) = .TRUE.
          end if
        end if
       end do
      end do
    end do
  io_int(1) = DynNbeams
  call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
end if   ! method = ALL


! Determine the masterlist of reflections for the zone axis case with HOLZ;
if (method.eq.'ZA') then
! set threshold for double diffraction detection
  ddt = 1.0e-6

! distance between consecutive HOLZ layers in nm-1
  H = 1.0/CalcLength(float(k),'d')

! determine g3 basis vector (vector normal to ZOLZ, with length H)
  call CalcCross(float(ga),float(gb),g3,'r','r',1)
  call NormVec(g3,'r')
  g3n = g3
  g3 = H * g3
  FNg = (/ CalcDot(DynFN,float(ga),'r'), CalcDot(DynFN,float(gb),'r'), CalcDot(DynFN,g3,'r') /)
    
  io_real = (/ float(ga(1:3)), float(gb(1:3)), g3(1:3) /)
  call WriteValue('basis vectors for this computation: ', io_real, 9, "(/'ga = ',3f10.5,/'gb = ',3f10.5,/'g3 = ',3f10.5,/)")
  io_real(1) = H;
  call WriteValue('reciprocal interplanar spacing H = ', io_real, 1, "(F10.4,' nm^-1'/)")

  call ShortestGFOLZ(k,ga,gb,gshort,gp)
  io_int = gshort
  call WriteValue(' shortest vector to FOLZ = ', io_int, 3, "('(',3I3,')',/)")

! The master list is most easily created by brute force; we'll compute the 
! radius of the FOLZ ring, scale it by the length of ga or gb, turn that into an integer
! and make the range twice as large...  All reflections from FOLZ and ZOLZ within
! this range will become part of the list
  minholz = 0
  i = maxholz
  DynNbeams = 1
  RHOLZ = sqrt(2.0*H*float(i)/mLambda - (float(i)*H)**2)
  im = max( nint(RHOLZ/CalcLength(float(ga),'r')),nint(RHOLZ/CalcLength(float(gb),'r')) )

if (im.eq.0) im = 10

! allocate the look up table with potential coefficients
  allocate(LUT(-2*im:2*im,-2*im:2*im,-2*im:2*im),stat=istat)
  LUT = dcmplx(0.D0,0.D0)

! allocate an array that keeps track of potential double diffraction reflections
  allocate(dbdiff(-2*im:2*im,-2*im:2*im,-2*im:2*im),stat=istat)
  dbdiff = .FALSE.

  call AddReflection( (/0,0,0/) )   ! this guarantees that 000 is always the first reflection
  LUT(0,0,0) = rlp%Ucg
  do N=minholz,maxholz
    do ix=-2*im,2*im
      do iy=-2*im,2*im
        gg = ix*ga + iy*gb + N*gshort
        if ((IsGAllowed(gg)).AND.(CalcLength(float(gg),'r').ne.0.0)) then
	  if  ((abs(gg(1)).le.2*im).and.(abs(gg(2)).le.2*im).and.(abs(gg(3)).le.2*im) ) then
             call AddReflection(gg) 
          else
            call CalcUcg(gg)
          end if
! add the reflection to the look up table
          LUT(gg(1), gg(2), gg(3)) = rlp%Ucg
! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
          if (cabs(rlp%Ucg).le.ddt) then 
            dbdiff(ix,iy,iz) = .TRUE.
          end if
        end if
      end do
    end do
  end do
  io_int(1) = DynNbeams
  call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I5,/)")
end if

! copy the DynNbeams variable
DynNbeamsLinked = DynNbeams
 
! should we convert the linked list to a regular set of arrays ?
if (ConvertList) then
  allocate(Reflist_hkl(3,DynNbeams), Reflist_Ucg(DynNbeams), Reflist_sg(DynNbeams), Reflist_xg(DynNbeams), stat = istat)
  rltmpa => reflist%next
  do ig=1,DynNbeams
    Reflist_hkl(1:3,ig) = rltmpa%hkl(1:3)
    Reflist_Ucg(ig) = rltmpa%Ucg
    Reflist_sg(ig) = rltmpa%sg
    Reflist_xg(ig) = rltmpa%xg
    rltmpa => rltmpa%next
  end do
end if 
 
 
 
 
end subroutine Compute_ReflectionList




!--------------------------------------------------------------------------
!
! SUBROUTINE: ShortestGFOLZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brieffind the G vector (displacement FOLZ w.r.t. ZOLZ, Chapter 3)
!
!> @param k wavevector
!> @param ga first reflection
!> @param gb second reflection
!> @param  gshort shortest g-vector inclined to ga-gb plane
!> @param gp
!
!> @date  01/29/02 MDG 1.0 original
!> @date   04/29/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine ShortestGFOLZ(k,ga,gb,gshort,gp)

use local
use io
use crystal
use crystalvars
use error

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: k(3)
integer(kind=irg),INTENT(IN)	:: ga(3)
integer(kind=irg),INTENT(IN)	:: gb(3)
integer(kind=irg),INTENT(OUT)	:: gshort(3)
integer(kind=irg),INTENT(OUT)	:: gp(3)


real(kind=sgl)               	:: gmin,gam11,gam12,gam22,g1(3),g2(3),g3(3),glen
integer(kind=irg),parameter  	:: inm = 8
integer(kind=irg)            	:: ih,ik,il,NN, io_int(1)

! look for the shortest reflection satisfying hu+kv+lw = 1
! This could be replaced by code from Jackson's paper (1987),
! but it does essentially the same thing.
 gmin = 100.0
 NN=1
 g1 = float(ga)
 g2 = float(gb)
 do while((gmin.eq.100.0).and.(NN.lt.4))
  do ih=-inm,inm
   do ik=-inm,inm
    do il=-inm,inm
! does this reflection lie in the plane NN ?
     if ((ih*k(1)+ik*k(2)+il*k(3)).eq.NN) then
      glen = CalcLength(float((/ih,ik,il/)),'r')
      if (glen.lt.gmin) then
       gmin = glen
       gshort =  (/ ih,ik,il /) 
      end if
     end if
    end do
   end do
  end do
  if (gmin.eq.100.0) then 
    io_int(1) = NN
    call WriteValue(' Could not find any reflections with hu+kv+lw = ', io_int, 1, "(I2)")
    NN = NN+1
  end if
 end do
 if (gmin.eq.100.0) then ! for some reason there is no reflection with N<=3 ...
  call FatalError('ShortestGFOLZ','HOLZ: could not find any reflections with hu+kv+lw<=maxholz ...')
 end if
 g3 = float(gshort)
! projected components of G
 gam11 = CalcDot(g1,g1,'r')
 gam12 = CalcDot(g1,g2,'r')
 gam22 = CalcDot(g2,g2,'r')
 gmin = 1.0/(gam11*gam22-gam12**2)
 gp(1) = (CalcDot(g3,g1,'r')*gam22-CalcDot(g3,g2,'r')*gam12)*gmin
 gp(2) = (CalcDot(g3,g2,'r')*gam11-CalcDot(g3,g1,'r')*gam12)*gmin

end subroutine ShortestGFOLZ




!--------------------------------------------------------------------------
!
! SUBROUTINE: Delete_gvectorlist
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief delete the entire linked list
!
!> @date   04/29/13 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine Delete_gvectorlist()

IMPLICIT NONE

! deallocate the entire linked list before returning, to prevent memory leaks
rltail => reflist
rltmpa => rltail % next
do 
  deallocate(rltail)
  if (.not. associated(rltmpa)) EXIT
  rltail => rltmpa
  rltmpa => rltail % next
end do

end subroutine Delete_gvectorlist


end module gvectors
