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
! CTEMsoft:gvectors.f90 
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
!> @date 04/29/13 MDG 1.0 original 
!> #date 01/10/14 MDG 2.0 new version, now use-d in crystalvars
!> @date 06/09/14 MDG 3.0 removed all global variables and replaced them by arguments
!--------------------------------------------------------------------------
module gvectors

use local
use typedefs

IMPLICIT NONE

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: MakeRefList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief allocate and initialize the linked reflection list
!
!> @param cell unit cell pointer
!> @param rltail reflist pointer
!
!> @date  10/20/98 MDG 1.0 original
!> @date   5/22/01 MDG 2.0 f90
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  03/26/13 MDG 3.0 updated IO
!> @date  01/10/14 MDG 4.0 account for new version of cell type
!> @date  06/09/14 MDG 4.1 added cell and rltail as arguments
!--------------------------------------------------------------------------
subroutine MakeRefList(cell, rltail)

use error

IMPLICIT NONE

type(unitcell),pointer	        :: cell
type(reflisttype),pointer	:: rltail

integer(kind=irg)  :: istat

! create it if it does not already exist
if (.not.associated(cell%reflist)) then
  cell%DynNbeams = 0
  allocate(cell%reflist,stat=istat)
  if (istat.ne.0) call FatalError('MakeRefList:',' unable to allocate pointer')
  rltail => cell%reflist              ! tail points to new value
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
!> @param rltail reflisttype variable
!> @param cell unit cell pointer
!> @param hkl Miller indices
!
!> @date  10/20/98 MDG 1.0 original
!> @date   5/22/01 MDG 2.0 f90
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  03/26/13 MDG 3.0 updated IO
!> @date  01/10/14 MDG 4.0 account for new version of cell type
!> @date  06/09/14 MDG 4.1 added rltail and cell as arguments
!--------------------------------------------------------------------------
subroutine AddReflection(rltail,cell,hkl)

use error
use diffraction

IMPLICIT NONE

type(reflisttype),pointer	:: rltail
type(unitcell),pointer	        :: cell
integer(kind=irg),INTENT(IN)	:: hkl(3)		!< Miller indices of reflection to be added to list

integer(kind=irg)  		:: istat

! create reflist if it does not already exist
 if (.not.associated(cell%reflist)) call MakeRefList(cell,rltail)

! create a new entry
 allocate(rltail%next,stat=istat)  		! allocate new value
 if (istat.ne.0) call FatalError('AddReflection',' unable to add new reflection')

 rltail => rltail%next             		! tail points to new value
 nullify(rltail%next)              		! nullify next in new value

 cell%DynNbeams = cell%DynNbeams + 1         	! update reflection counter
 rltail%num = cell%DynNbeams            	! store reflection number
 rltail%hkl = hkl                  		! store Miller indices
 rltail%Ucg = cell%LUT( hkl(1), hkl(2), hkl(3) ) ! and store it in the list
 rltail%famnum = 0				! init this value for Prune_ReflectionList
! rltail%Ucgmod = cabs(rlp%Ucg)   		! added on 2/29/2012 for Bethe potential computations
! rltail%sangle = 1000.0*dble(CalcDiffAngle(hkl(1),hkl(2),hkl(3)))    ! added 4/18/2012 for EIC project HAADF/BF tomography simulations
! rltail%thetag = rlp%Vphase                   ! added 12/14/2013 for CTEMECCI program
 nullify(rltail%nextw)
 nullify(rltail%nexts)
 
end subroutine AddReflection

!--------------------------------------------------------------------------
!
! SUBROUTINE: Printrlp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief output the contents of the rlp structure
!
!> @param rlp gnode structure
!> @param first logical switch to provide long form output (optional)
!> @param stdout optional output unit identifier
!
!> @date  10/20/98 MDG 1.0 original
!> @date   5/22/01 MDG 2.0 f90
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  03/26/13 MDG 3.0 updated IO
!> @date  06/09/14 MDG 4.0 added rlp and stdout arguments
!--------------------------------------------------------------------------
subroutine Printrlp(rlp,first,stdout)

use io
use constants

IMPLICIT NONE

type(gnode),INTENT(IN)	               :: rlp
logical,optional,intent(INOUT) 	:: first		!< switch for long/short output
integer(kind=irg),OPTIONAL,INTENT(IN) :: stdout

integer(kind=irg)			:: oi_int(3), std
real(kind=sgl)				:: oi_real(7)
complex(kind=sgl)			:: oi_cmplx(1)

std = 6
if (PRESENT(stdout)) std = stdout

if (present(first)) then
 if (first) then
  call Message('     Scattering factors : ', frm = "(/A,$)", stdout = std)
    if (rlp%method.eq.'WK') then 
   if (rlp%absorption.eqv..TRUE.) then 
    call Message(' Weickenmeier-Kohl (with absorption)', frm = "(A/)", stdout = std)
   else
    call Message(' Weickenmeier-Kohl', frm = "(A/)", stdout = std)
   end if
  else
    call Message(' Doyle-Turner/Smith-Burge', frm = "(A/)", stdout = std)
  end if

  if (rlp%absorption.eqv..TRUE.) then
    call Message('   h  k  l    |g|    Ucg_r  Ucg_i   |Ug|    phase   |Ugp|   phase   xi_g   xi_gp    ratio  Re-1/q_g-Im', &
        frm = "(A)", stdout = std)
  else
    call Message('   h  k  l    |g|    Ucg_r  |Ug|    phase    xi_g   1/q_g', frm = "(A)", stdout = std)
  end if
  first = .FALSE.
 end if
end if

if (rlp%absorption.eqv..TRUE.) then
 oi_int(1:3) = rlp%hkl(1:3)
 call WriteValue('',oi_int, 3, "(1x,3I3,1x,$)", stdout = std)
 oi_real(1) = rlp%g
 call WriteValue('',oi_real, 1, "(F9.4,$)", stdout = std)
 oi_cmplx(1) = rlp%Ucg
 call WriteValue('',oi_cmplx, 1, "(2F7.3,1x,$)", stdout = std)
 oi_real(1:7)  = (/ rlp%Umod,rlp%Vphase*180.0/sngl(cPi),rlp%Upmod,rlp%Vpphase*180.0/sngl(cPi),rlp%xg,rlp%xgp,rlp%ar /)
 call WriteValue('',oi_real, 7, "(4F8.3,3F8.1,$)", stdout = std)
 oi_cmplx(1) = rlp%qg
 call WriteValue('',oi_cmplx, 1, "(2F8.3)", stdout = std)
else
 oi_int(1:3) = rlp%hkl(1:3)
 call WriteValue('',oi_int, 3, "(1x,3I3,1x,$)", stdout = std)
 oi_real(1) = rlp%g
 call WriteValue('',oi_real, 1, "(F9.4,$)", stdout = std)
 oi_real(1) = real(rlp%Ucg)
 call WriteValue('',oi_real, 1, "(F7.3,1x,$)", stdout = std)
 oi_real(1:3)  = (/ rlp%Umod,rlp%Vphase*180.0/sngl(cPi),rlp%xg /)
 call WriteValue('',oi_real, 3, "(2F8.3,F8.1,$)", stdout = std)
 oi_cmplx(1) = rlp%qg
 call WriteValue('',oi_cmplx, 1, "(2F8.3)", stdout = std)
end if

end subroutine Printrlp

!--------------------------------------------------------------------------
!
! SUBROUTINE: Apply_BethePotentials
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief tag weak and strong reflections in cell%reflist
!
!> @param cell unit cell pointer
!
!> @details This routine steps through the cell%reflist linked list and 
!> determines for each reflection whether it is strong or weak or should be
!> ignored.  Strong and weak reflections are then linked in a new list via
!> the nexts and nextw pointers, along with the nns and nnw counters.
!> This routine produces no output parameters, and makes use of the
!> BetheParameter variables.
!
!> @date  01/14/14 MDG 1.0 original version
!> @date  06/09/14 MDG 2.0 added cell and BetheParameter arguments
!--------------------------------------------------------------------------
subroutine Apply_BethePotentials(cell, BetheParameter)

use diffraction

IMPLICIT NONE

type(unitcell),pointer                         :: cell
type(BetheParameterType),INTENT(INOUT)         :: BetheParameter

integer(kind=irg),allocatable   :: glist(:,:)
real(kind=dbl),allocatable      :: rh(:)
type(reflisttype),pointer       :: rl, lastw, lasts
integer(kind=irg)               :: icnt, istat, gmh(3), ir, ih
real(kind=dbl)                  :: sgp, la, m


la = 1.D0/mLambda

nullify(cell%firstw)

! first we extract the list of g-vectors from reflist, so that we can compute 
! all the g-h difference vectors
allocate(glist(3,cell%DynNbeams),rh(cell%DynNbeams),stat=istat)
rl => cell%reflist%next
icnt = 0
do
  if (.not.associated(rl)) EXIT
  icnt = icnt+1
  glist(1:3,icnt) = rl%hkl(1:3)
  rl => rl%next
end do

! initialize the strong and weak reflection counters
cell%nns = 1
cell%nnw = 0

! the first reflection is always strong
rl => cell%reflist%next
rl%strong = .TRUE.
rl%weak = .FALSE.
lasts => rl
nullify(lasts%nextw)

! next we need to iterate through all reflections in glist and 
! determine which category the reflection belongs to: strong, weak, ignore
irloop: do ir = 2,icnt
  rl => rl%next
  rh = 0.D0
  sgp = la * abs(rl%sg)
  do ih = 1,icnt
   gmh(1:3) = glist(1:3,ir) - glist(1:3,ih)
   if (cell%dbdiff(gmh(1),gmh(2),gmh(3))) then  ! it is a double diffraction reflection with |U|=0
! to be written
     rh(ih) = 0.D0 
   else 
     rh(ih) = sgp/cdabs( cell%LUT(gmh(1), gmh(2), gmh(3)) )
   end if
  end do

! which category does reflection ir belong to ?
  m = minval(rh)
! write (*,*) glist(1:3,ir),'; minval(rh) = ',m

! m > c2 => ignore this reflection
  if (m.gt.BetheParameter%c2) then
!    write (*,*) '    -> ignore'
    rl%weak = .FALSE.
    rl%strong = .FALSE.
    CYCLE irloop
  end if
! c1 < m < c2 => weak reflection
  if ((BetheParameter%c1.lt.m).and.(m.le.BetheParameter%c2)) then
    if (cell%nnw.eq.0) then
      cell%firstw => rl
      lastw => rl
    else
      lastw%nextw => rl
      lastw => rl
      nullify(lastw%nexts)
    end if
    rl%weak = .TRUE.
    rl%strong = .FALSE.
    cell%nnw = cell%nnw + 1
!    write (*,*) '    -> weak ',cell%nnw
    CYCLE irloop
  end if

! m < c1 => strong
  if (m.le.BetheParameter%c1) then
    lasts%nexts => rl
    nullify(lasts%nextw)
    lasts => rl
    rl%weak = .FALSE.
    rl%strong = .TRUE.
    cell%nns = cell%nns + 1
!    write (*,*) '    -> strong ',cell%nns
  end if  
end do irloop

deallocate(glist, rh)

end subroutine Apply_BethePotentials


!--------------------------------------------------------------------------
!
! SUBROUTINE: Prune_ReflectionList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief select from the reflection list those g-vectors that will be counted in an LACBED computation
!
!> @details This routine basicaly repeats a section from the Compute_DynMat routine
!> without actually computing the dynamical matrix; it simply keeps track of all the
!> beams that are at one point or another regarded as strong or weak beams.
!> We'll use the famnum field in the rlp linked list to flag the strong reflections.
!> Linked list entries that are not used are instantly removed from the linked list.
!
!> @param cell unit cell pointer
!> @param khead start of reflisttype linked list
!> @param reflist reflection linked list
!> @param Dyn dynamical scattering structure
!> @param BetheParameter Bethe parameter structure 
!> @param numk number of wave vectors to consider
!> @param nbeams total number of unique beams
!
!> @date  09/20/13 MDG 1.0 original
!> @date  10/05/13 MDG 1.1 removal of unused reflections from linked list
!> @date  10/05/13 MDG 1.2 changed the order of nested loops to speed things up a bit
!> @date  10/07/13 MDG 1.3 added section to reset the famhkl entries after pruning
!> @date  01/10/14 MDG 4.0 account for new version of cell type
!> @date  06/09/14 MDG 4.1 added cell, reflist, BetheParameter and khead arguments
!--------------------------------------------------------------------------
subroutine Prune_ReflectionList(cell,khead,reflist,Dyn,BetheParameter,numk,nbeams)

use io
use crystal
use kvectors
use diffraction

IMPLICIT NONE

type(unitcell),pointer	                :: cell
type(kvectorlist),pointer              :: khead
type(reflisttype),pointer	        :: reflist
type(DynType),INTENT(INOUT)            :: Dyn
type(BetheParameterType),INTENT(INOUT):: BetheParameter
integer(kind=irg),INTENT(IN)		:: numk
integer(kind=irg),INTENT(OUT)		:: nbeams

integer(kind=irg)			:: ik, ig, istrong, curfam(3), newfam(3)
real(kind=sgl)     			:: sgp, lUg, cut1, cut2
!integer(kind=irg),allocatable		:: strongreflections(:,:)
type(kvectorlist),pointer	        :: ktmp
type(reflisttype),pointer	        :: rltmpa, rltmpb

! reset the value of DynNbeams in case it was modified in a previous call 
cell%DynNbeams = cell%DynNbeamsLinked

nbeams = 0

! reset the reflection linked list
  rltmpa => cell%reflist%next

! pick the first reflection since that is the transmitted beam (only on the first time)
  rltmpa%famnum = 1    
  nbeams = nbeams + 1

! loop over all reflections in the linked list    
!!!! this will all need to be changed with the new Bethe potential criteria ...  
  rltmpa => rltmpa%next
  reflectionloop: do ig=2,cell%DynNbeamsLinked
    lUg = cdabs(rltmpa%Ucg) * mLambda
    cut1 = BetheParameter%cutoff * lUg
    cut2 = BetheParameter%weakcutoff * lUg

! loop over all the incident beam directions
    ktmp => khead
! loop over all beam orientations, selecting them from the linked list
    kvectorloop: do ik = 1,numk
! We compare |sg| with two multiples of lambda |Ug|
!
!  |sg|>cutoff lambda |Ug|   ->  don't count reflection
!  cutoff lambda |Ug| > |sg| > weakcutoff lambda |Ug|  -> weak reflection
!  weakcutoff lambda |Ug| > |sg|  -> strong reflection
!
! 	sgp = abs(CalcsgHOLZ(float(rltmpa%hkl),sngl(ktmp%kt),sngl(mLambda)))
! 	write (*,*) rltmpa%hkl,CalcsgHOLZ(float(rltmpa%hkl),sngl(ktmp%kt), &
!			sngl(mLambda)),Calcsg(float(rltmpa%hkl),sngl(ktmp%k),DynFN)
        sgp = abs(Calcsg(cell,float(rltmpa%hkl),sngl(ktmp%k),Dyn%FN)) 
! we have to deal separately with double diffraction reflections, since
! they have a zero potential coefficient !        
        if ( cell%dbdiff(rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3)) ) then  ! it is a double diffraction reflection
          if (sgp.le.BetheParameter%sgcutoff) then         
	    nbeams = nbeams + 1
    	    rltmpa%famnum = 1    
	    EXIT kvectorloop	! this beam did contribute, so we no longer need to consider it
          end if
        else   ! it is not a double diffraction reflection
          if (sgp.le.cut1) then  ! count this beam, whether it is weak or strong
	    nbeams = nbeams + 1
	    rltmpa%famnum = 1    
	    EXIT kvectorloop	! this beam did contribute, so we no longer need to consider it
          end if
        end if
 
! go to the next incident beam direction
       if (ik.ne.numk) ktmp => ktmp%next
     end do kvectorloop  ! ik loop

! go to the next beam in the list
   rltmpa => rltmpa%next
  end do reflectionloop

  call Message(' Renumbering reflections', frm = "(A)")

! change the following with the new next2 pointer in the reflist type !!!
  
! ok, now that we have the list, we'll go through it again to set sequential numbers instead of 1's
! at the same time, we'll deallocate those entries that are no longer needed.
  rltmpa => reflist%next
  rltmpb => rltmpa
  rltmpa => rltmpa%next	! we keep the first entry, always.
  istrong = 1
  reflectionloop2: do ig=2,cell%DynNbeamsLinked
    if (rltmpa%famnum.eq.1) then
	istrong = istrong + 1
      	rltmpa%famnum = istrong
	rltmpa => rltmpa%next
	rltmpb => rltmpb%next
    else   ! remove this entry from the linked list
	rltmpb%next => rltmpa%next
	deallocate(rltmpa)
	rltmpa => rltmpb%next
    endif
! go to the next beam in the list
  end do reflectionloop2

! reset the number of beams to the newly obtained number
  cell%DynNbeamsLinked = nbeams
  cell%DynNbeams = nbeams

! go through the entire list once again to correct the famhkl
! entries, which may be incorrect now; famhkl is supposed to be one of the 
! reflections on the current list, but that might not be the case since
! famhkl was first initialized when there were additional reflections on
! the list... so we set famhkl to be the same as the first hkl in each family.
  rltmpa => reflist%next%next  ! no need to check the first one
reflectionloop3:  do while (associated(rltmpa))
    curfam = rltmpa%famhkl
    if (sum(abs(curfam-rltmpa%hkl)).ne.0) then
      newfam = rltmpa%hkl
      do while (sum(abs(rltmpa%famhkl-curfam)).eq.0)
        rltmpa%famhkl = newfam
        rltmpa => rltmpa%next
	if ( .not.associated(rltmpa) ) EXIT reflectionloop3
      end do
    else
      do while (sum(abs(rltmpa%famhkl-curfam)).eq.0)
        rltmpa => rltmpa%next
	if ( .not.associated(rltmpa) ) EXIT reflectionloop3
      end do
    end if
! go to the next beam in the list
  end do reflectionloop3

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
!> @param cell unit cell pointer
!> @param reflist reflection list pointer
!> @param Dyn dynamical scattering structure
!> @param BetheParameter Bethe parameter structure 
!> @param calcmode string that describes the particular matrix mode
!> @param kk incident wave vector
!> @param kt tangential component of incident wave vector (encodes the Laue Center)
!> @param IgnoreFoilNormal switch for foil normal inclusion in sg computation
!> @param IncludeSecondOrder (optional) switch to include second order correction to Bethe potentials
!
!> @date   05/06/13 MDG 1.0 original
!> @date   08/30/13 MDG 1.1 correction of effective excitation error
!> @date   09/20/13 MDG 1.2 added second order Bethe potential correction switch
!> @date   06/09/14 MDG 2.0 added cell, reflist and BetheParameter as arguments
!--------------------------------------------------------------------------
recursive subroutine Compute_DynMat(cell,reflist,Dyn,BetheParameter,calcmode,kk,kt,IgnoreFoilNormal,IncludeSecondOrder)

use error
use constants
use crystal
use diffraction
use io

IMPLICIT NONE

type(unitcell),pointer	                :: cell
type(reflisttype),pointer	        :: reflist
type(DynType),INTENT(INOUT)            :: Dyn
type(BetheParameterType),INTENT(INOUT) :: BetheParameter
character(*),INTENT(IN)		:: calcmode		!< computation mode
real(kind=dbl),INTENT(IN)		:: kk(3),kt(3)		!< incident wave vector and tangential component
logical,INTENT(IN)			:: IgnoreFoilNormal	!< how to deal with the foil normal
logical,INTENT(IN),OPTIONAL		:: IncludeSecondOrder	!< second order Bethe potential correction switch

complex(kind=dbl)  			:: czero,pre, weaksum, ughp, uhph
integer(kind=irg) 		 	:: istat,ir,ic,nn, iweak, istrong, iw, ig, ll(3), gh(3), nnn, nweak, io_int(1)
real(kind=sgl)     			:: glen,exer,gg(3), kpg(3), gplen, sgp, lUg, cut1, cut2, io_real(3)
real(kind=dbl)				:: lsfour, weaksgsum 
logical					:: AddSecondOrder
type(gnode)                           :: rlp
type(reflisttype),pointer	        :: rltmpa, rltmpb
 
AddSecondOrder = .FALSE.
if (present(IncludeSecondOrder)) AddSecondOrder = .TRUE.

! has the list of reflections been allocated ?
if (.not.associated(reflist)) call FatalError('Compute_DynMat',' reflection list has not been allocated')

! if the dynamical matrix has already been allocated, deallocate it first
! this is partially so that no program will allocate DynMat itself; it must be done
! via this routine only.
if (allocated(Dyn%DynMat)) deallocate(Dyn%DynMat)

! initialize some parameters
czero = cmplx(0.0,0.0,dbl)	! complex zero
pre = cmplx(0.0,cPi,dbl)		! i times pi

if (calcmode.ne.'BLOCHBETHE') then

! allocate DynMat
	  allocate(Dyn%DynMat(cell%DynNbeams,cell%DynNbeams),stat=istat)
	  Dyn%DynMat = czero
! get the absorption coefficient
	  call CalcUcg(cell, rlp, (/0,0,0/) )
	  Dyn%Upz = rlp%Vpmod

! are we supposed to fill the off-diagonal part ?
	 if ((calcmode.eq.'D-H-W').or.(calcmode.eq.'BLOCH')) then
	  rltmpa => cell%reflist%next    ! point to the front of the list
! ir is the row index
	  do ir=1,cell%DynNbeams
	   rltmpb => cell%reflist%next   ! point to the front of the list
! ic is the column index
	   do ic=1,cell%DynNbeams
	    if (ic.ne.ir) then  ! exclude the diagonal
! compute Fourier coefficient of electrostatic lattice potential 
	     gh = rltmpa%hkl - rltmpb%hkl
	     if (calcmode.eq.'D-H-W') then
	      call CalcUcg(cell, rlp,gh)
	      Dyn%DynMat(ir,ic) = pre*rlp%qg
	     else
	      Dyn%DynMat(ir,ic) = cell%LUT(gh(1),gh(2),gh(3))
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
	  do ir=1,cell%DynNbeams
	   glen = CalcLength(cell,float(rltmpa%hkl),'r')
	   if (glen.eq.0.0) then
	    Dyn%DynMat(ir,ir) = cmplx(0.0,Dyn%Upz,dbl)
	   else  ! compute the excitation error
	    exer = Calcsg(cell,float(rltmpa%hkl),sngl(kk),Dyn%FN)

	    rltmpa%sg = exer
	    if (calcmode.eq.'DIAGH') then  !
	     Dyn%DynMat(ir,ir) = cmplx(0.0,2.D0*cPi*exer,dbl)
	    else
	     Dyn%DynMat(ir,ir) = cmplx(2.D0*exer/mLambda,Dyn%Upz,dbl)
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
	if (BetheParameter%cutoff.eq.0.0) call Set_Bethe_Parameters(BetheParameter)


! reset the value of DynNbeams in case it was modified in a previous call 
  	cell%DynNbeams = cell%DynNbeamsLinked
  	
! precompute lambda^2/4
	lsfour = mLambda**2*0.25D0
  
! first, for the input beam direction, determine the excitation errors of 
! all the reflections in the master list, and count the ones that are
! needed for the dynamical matrix (weak as well as strong)
        if (.not.allocated(BetheParameter%weaklist)) allocate(BetheParameter%weaklist(cell%DynNbeams))
        if (.not.allocated(BetheParameter%stronglist)) allocate(BetheParameter%stronglist(cell%DynNbeams))
        if (.not.allocated(BetheParameter%reflistindex)) allocate(BetheParameter%reflistindex(cell%DynNbeams))
        if (.not.allocated(BetheParameter%weakreflistindex)) allocate(BetheParameter%weakreflistindex(cell%DynNbeams))

	BetheParameter%weaklist = 0
	BetheParameter%stronglist = 0
	BetheParameter%reflistindex = 0
	BetheParameter%weakreflistindex = 0

    	rltmpa => cell%reflist%next

! deal with the transmitted beam first
    nn = 1		! nn counts all the scattered beams that satisfy the cutoff condition
    nnn = 1		! nnn counts only the strong beams
    nweak = 0		! counts only the weak beams
    BetheParameter%stronglist(nn) = 1   ! make sure that the transmitted beam is always a strong beam ...
    BetheParameter%weaklist(nn) = 0
    BetheParameter%reflistindex(nn) = 1

    rltmpa%sg = 0.D0    
! write (*,*) 'DynNbeamsLinked = ',DynNbeamsLinked

! loop over all reflections in the linked list    
    rltmpa => rltmpa%next
    reflectionloop: do ig=2,cell%DynNbeamsLinked
      gg = float(rltmpa%hkl)        		! this is the reciprocal lattice vector 

! deal with the foil normal; if IgnoreFoilNormal is .TRUE., then assume it is parallel to the beam direction
     if (IgnoreFoilNormal) then 
! we're taking the foil normal to be parallel to the incident beam direction at each point of
! the standard stereographic triangle, so cos(alpha) = 1 always in eqn. 5.11 of CTEM
        kpg = kk+gg                		! k0 + g (vectors)
        gplen = CalcLength(cell,kpg,'r')  	! |k0+g|
        rltmpa%sg = (1.0/mLambda**2 - gplen**2)*0.5/gplen
     else
	rltmpa%sg = Calcsg(cell,gg,sngl(kk),Dyn%FN)
! here we need to determine the components of the Laue Center w.r.t. the g1 and g2 vectors
! and then pass those on to the routine; 
!	rltmpa%sg = CalcsgHOLZ(gg,sngl(kt),sngl(mLambda))
! write (*,*) gg, Calcsg(gg,sngl(kk),DynFN), CalcsgHOLZ(gg,sngl(kt),sngl(mLambda))
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
        lUg = cdabs(rltmpa%Ucg) * mLambda
        cut1 = BetheParameter%cutoff * lUg
        cut2 = BetheParameter%weakcutoff * lUg

! we have to deal separately with double diffraction reflections, since
! they have a zero potential coefficient !        
        if ( cell%dbdiff(rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3)) ) then  ! it is a double diffraction reflection
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
              		nweak = nweak+1
              		BetheParameter%weaklist(ig) = 1
              		BetheParameter%weakreflistindex(ig) = nweak
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
	   call Message(' no beams found for the following parameters:', frm = "(A)")
	   io_real(1:3) = kk(1:3)
           call WriteValue(' wave vector = ', io_real,3)
           io_int(1) = nn
           call WriteValue('  -> number of beams = ', io_int, 1)
	   call Message( '   -> check cutoff and weakcutoff parameters for reasonableness', frm = "(A)")
	   call FatalError('Compute_DynMat','No beams in list')
	end if

! next, we define nns to be the number of strong beams, and nnw the number of weak beams.
	 BetheParameter%nns = sum(BetheParameter%stronglist)
	 BetheParameter%nnw = sum(BetheParameter%weaklist)
	 
! add nns to the weakreflistindex to offset it; this is used for plotting reflections on CBED patterns
	do ig=2,cell%DynNbeamsLinked
	  if (BetheParameter%weakreflistindex(ig).ne.0) then
	    BetheParameter%weakreflistindex(ig) = BetheParameter%weakreflistindex(ig) + BetheParameter%nns
	  end if
	end do
	
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
	if (allocated(BetheParameter%strongID)) deallocate(BetheParameter%strongID)
	allocate(BetheParameter%weakhkl(3,BetheParameter%nnw),BetheParameter%weaksg(BetheParameter%nnw))
	allocate(BetheParameter%stronghkl(3,BetheParameter%nns),BetheParameter%strongsg(BetheParameter%nns))
        allocate(BetheParameter%strongID(BetheParameter%nns))

! here's where we extract the relevant information from the linked list (much faster
! than traversing the list each time...)
	rltmpa => cell%reflist%next    ! reset the a list
	iweak = 0
	istrong = 0
	do ir=1,cell%DynNbeamsLinked
	     if (BetheParameter%weaklist(ir).eq.1) then
	        iweak = iweak+1
	        BetheParameter%weakhkl(1:3,iweak) = rltmpa%hkl(1:3)
	        BetheParameter%weaksg(iweak) = rltmpa%sg
	     end if
	     if (BetheParameter%stronglist(ir).eq.1) then
	        istrong = istrong+1
	        BetheParameter%stronghkl(1:3,istrong) = rltmpa%hkl(1:3)
	        BetheParameter%strongsg(istrong) = rltmpa%sg
! make an inverse index list
		BetheParameter%strongID(istrong) = ir		
	     end if
	   rltmpa => rltmpa%next
	end do

! now we are ready to create the dynamical matrix
	cell%DynNbeams = BetheParameter%nns

! allocate DynMat if it hasn't already been allocated and set to complex zero
	  if (allocated(Dyn%DynMat)) deallocate(Dyn%DynMat)
	  allocate(Dyn%DynMat(cell%DynNbeams,cell%DynNbeams),stat=istat)
	  Dyn%DynMat = czero

! get the absorption coefficient
	  call CalcUcg(cell, rlp, (/0,0,0/) )
	  Dyn%Upz = rlp%Vpmod

! ir is the row index
       do ir=1,BetheParameter%nns
! ic is the column index
          do ic=1,BetheParameter%nns
! compute the Bethe Fourier coefficient of the electrostatic lattice potential 
              if (ic.ne.ir) then  ! not a diagonal entry
                 ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
                 Dyn%DynMat(ir,ic) = cell%LUT(ll(1),ll(2),ll(3)) 
        ! and subtract from this the total contribution of the weak beams
         weaksum = czero
         do iw=1,BetheParameter%nnw
              ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%weakhkl(1:3,iw)
              ughp = cell%LUT(ll(1),ll(2),ll(3)) 
              ll = BetheParameter%weakhkl(1:3,iw) - BetheParameter%stronghkl(1:3,ic)
              uhph = cell%LUT(ll(1),ll(2),ll(3)) 
              weaksum = weaksum +  ughp * uhph *cmplx(1.D0/BetheParameter%weaksg(iw),0.0,dbl)
         end do
        ! and correct the dynamical matrix element to become a Bethe potential coefficient
         Dyn%DynMat(ir,ic) = Dyn%DynMat(ir,ic) - cmplx(0.5D0*mLambda,0.0D0,dbl)*weaksum
! do we need to add the second order corrections ?
		  if (AddSecondOrder) then 
		    weaksum = czero
		  end if
              else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! determine the total contribution of the weak beams
                 weaksgsum = 0.D0
		  do iw=1,BetheParameter%nnw
                      ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%weakhkl(1:3,iw)
                      ughp = cell%LUT(ll(1),ll(2),ll(3)) 
                      weaksgsum = weaksgsum +  cdabs(ughp)**2/BetheParameter%weaksg(iw)
                 end do
                 weaksgsum = weaksgsum * mLambda/2.D0
                 Dyn%DynMat(ir,ir) = cmplx(2.D0*BetheParameter%strongsg(ir)/mLambda-weaksgsum,Dyn%Upz,dbl)
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
!> @param silent optional, if present, then no output
!> @date   05/08/13 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine Set_Bethe_Parameters(BetheParameter,silent)

use io

IMPLICIT NONE

type(BetheParameterType),INTENT(INOUT)        :: BetheParameter
logical,INTENT(IN),OPTIONAL	:: silent

character(fnlen),parameter 	:: Bethefilename = 'BetheParameters.nml'
logical				:: fexist
real(kind=sgl)			:: weakcutoff, cutoff, sgcutoff

namelist /BetheList/ weakcutoff, cutoff, sgcutoff

! check for the presence of the namelist file in the current folder
inquire(file=trim(Bethefilename),exist=fexist)

! set all default values (must be done here, since nml file may not contain all of them)
weakcutoff = 80.0  	! dimensionless cutoff parameter, smaller = strong, larger = weak
cutoff = 160.0		! overall cutoff parameter
sgcutoff = 0.05		! sg cutoff for double diffraction reflections

if (fexist) then ! check for the file in the local folder
! read the parameters from the file
 OPEN(UNIT=dataunit,FILE=trim(Bethefilename),DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=BetheList)
 CLOSE(UNIT=dataunit)
 if (.not.present(silent)) then
   call Message('Read Bethe parameters from BetheParameters.nml', frm = "(A)")
   write (6,nml=BetheList)
 end if
end if

BetheParameter % weakcutoff = weakcutoff
BetheParameter % cutoff = cutoff
BetheParameter % sgcutoff = sgcutoff

end subroutine Set_Bethe_Parameters



!--------------------------------------------------------------------------
!
! SUBROUTINE: ShortestGFOLZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brieffind the G vector (displacement FOLZ w.r.t. ZOLZ, Chapter 3)
!
!> @param cell unit cell pointer
!> @param k wavevector
!> @param ga first reflection
!> @param gb second reflection
!> @param  gshort shortest g-vector inclined to ga-gb plane
!> @param gp
!
!> @date  01/29/02 MDG 1.0 original
!> @date  04/29/13 MDG 2.0 rewrite
!> @date  06/09/14 MDG 2.1 added cell argument 
!--------------------------------------------------------------------------
subroutine ShortestGFOLZ(cell,k,ga,gb,gshort,gp)

use io
use crystal
use error

IMPLICIT NONE

type(unitcell),pointer	        :: cell
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
      glen = CalcLength(cell,float((/ih,ik,il/)),'r')
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
 gam11 = CalcDot(cell,g1,g1,'r')
 gam12 = CalcDot(cell,g1,g2,'r')
 gam22 = CalcDot(cell,g2,g2,'r')
 gmin = 1.0/(gam11*gam22-gam12**2)
 gp(1) = (CalcDot(cell,g3,g1,'r')*gam22-CalcDot(cell,g3,g2,'r')*gam12)*gmin
 gp(2) = (CalcDot(cell,g3,g2,'r')*gam11-CalcDot(cell,g3,g1,'r')*gam12)*gmin

end subroutine ShortestGFOLZ




!--------------------------------------------------------------------------
!
! SUBROUTINE: Delete_gvectorlist
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief delete the entire linked list
!
!> @param cell unit cell pointer
!
!> @date   04/29/13 MDG 1.0 original
!> @date   06/09/14 MDG 1.1 added cell argument
!--------------------------------------------------------------------------
subroutine Delete_gvectorlist(cell)

IMPLICIT NONE

type(unitcell),pointer	:: cell

type(reflisttype),pointer :: rltail, rltmpa

! deallocate the entire linked list before returning, to prevent memory leaks
rltail => cell%reflist
rltmpa => rltail % next
do 
  deallocate(rltail)
  if (.not. associated(rltmpa)) EXIT
  rltail => rltmpa
  rltmpa => rltail % next
end do

end subroutine Delete_gvectorlist


end module gvectors
