! ###################################################################
! Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
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
! CTEMsoft2013:kvectors.f90
!--------------------------------------------------------------------------
!
! MODULE: kvectors
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief variables and types needed to determine lists of wave vectors
!
!> @details This was part of the dynamical module, but was moved into a separate
!> module once I realized that I had several versions of the Calckvectors subroutine.
!> From now on, there is only one single version that can deal with all possible cases.
! 
!> @date   04/29/13 MDG 1.0 original
!--------------------------------------------------------------------------
module kvectors

use local

IMPLICIT NONE

! derived type definitions 

! linked list of wave vectors (used by all diffraction programs)
type kvectorlist  
  integer(kind=irg) 		:: i,j         		! image coordinates
  real(kind=dbl)    		:: kt(3)       	! tangential component of wavevector
  real(kind=dbl)    		:: kn          	! normal component
  real(kind=dbl)    		:: k(3)        	! full wave vector
  type(kvectorlist),pointer	:: next     		! connection to next wave vector
end type kvectorlist

contains

!--------------------------------------------------------------------------
!
! FUNCTION: Kdelta
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Kronecker delta function, returns 1 or 0
!
!> @param i first entry 
!> @param j second entry 
!
!> @date   04/29/13 MDG 1.0 original
!--------------------------------------------------------------------------
function Kdelta(i,j) result(res)

use local

IMPLICIT NONE

integer(kind=irg),INTENT(IN)  	:: i,j
integer(kind=irg)		:: res

 if (i.eq.j) then 
   res = 1
 else
   res = 0
 end if

end function Kdelta

!--------------------------------------------------------------------------
!
! SUBROUTINE: Calckvectors
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create a linked list of wave vectors
!
!> @details This is a new version that combines several older routines.  The most important 
!> aspects of this routine are a) linked list can use regular mapping or modified Lambert mapping;
!> b) list makes use of crystal symmetry (although that feature can be turned off); c) routine
!> has been cleaned up, and there is now a Delete_kvectorlist function as well.
!
!> @param khead head of linked list
!> @param cell unit cell pointer
!> @param k central wave vector
!> @param ga reciprocal lattice vector normal to k
!> @param ktmax maximum length of tangential component
!> @param npx number of vectors along x 
!> @param npy number of vectors along y 
!> @param numk total number of wave vectors in list
!> @param isym Laue group number
!> @param ijmax used for conical illumination
!> @param mapmode sets the mapping mode to be used ('Conical', 'Standard', 'StandardConical', 'RoscaLambert')
!
!> @todo The Standard and RoscaLambert mapmodes have different considerations of the 
!> Laue groups; this needs to be verified and, if necessary, simplified to a single set of 
!> conditions.  This might also allow Addkvector and Add_knode to become a single routine.
!
!> @date   04/29/13 MDG 1.0 original
!> @date   06/09/14 MDG 2.0 added khead and cell as arguments
!--------------------------------------------------------------------------
subroutine Calckvectors(khead,cell,k,ga,ktmax,npx,npy,numk,isym,ijmax,mapmode,usehex)

use io
use error
use constants
use diffraction
use crystal
use crystalvars

IMPLICIT NONE

type(kvectorlist),pointer              :: khead
type(unitcell),pointer	                :: cell
real(kind=dbl),INTENT(IN)		:: k(3)		!< initial wave vector
real(kind=dbl),INTENT(IN)		:: ga(3)	!< "horizontal" reciprocal lattice vector
real(kind=dbl),INTENT(IN)		:: ktmax	!< maximum length of tangential wave vector
integer(kind=irg),INTENT(IN)		:: npx		!< number of kvectors along x
integer(kind=irg),INTENT(IN)		:: npy		!< number of kvectors along y
integer(kind=irg),INTENT(OUT)		:: numk		!< total number of kvectors in linked list
integer(kind=irg),INTENT(IN)		:: isym		!< Laue symmetry group number 
integer(kind=irg),INTENT(INOUT)	:: ijmax	!< max parameter used for Conical and StandardConical modes
character(*),INTENT(IN)		:: mapmode 	!< controls the type of mapping used ('Standard' or 'RoscaLambert')
!real(kind=sgl),INTENT(IN)		:: klaue(2)	!< Laue center coordinates
logical,INTENT(IN),OPTIONAL		:: usehex	!< hexagonal mode for RoscaLambert mapmode

integer(kind=irg)       		:: istat,i,j,istart,iend,jstart,jend, imin, imax, jmin, jmax
real(kind=dbl)				:: glen, gan(3), gperp(3), kstar(3), delta
logical					:: hexgrid = .FALSE.
character(3)				:: grid
type(kvectorlist),pointer             :: ktail, ktmp

! first, if khead already exists, delete it
 if (associated(khead)) then    	 	! deallocate the entire linked list
    call Delete_kvectorlist(khead)
 end if   
 
! do we know this mapmode ?
if ( .not.( (mapmode.eq.'Conical').or.(mapmode.eq.'Standard').or.(mapmode.eq.'StandardConical').or. &
   (mapmode.eq.'RoscaLambert') ) ) then
  call FatalError('Calckvectors','mapmode unknown')
end if

if (mapmode.eq.'Conical') then ! used for CBED without symmetry application, including CTEMZAdefect
! compute geometrical factors 
 glen = CalcLength(cell,ga,'r')              		! length of ga
 gan = ga/glen                                 	! normalized ga
 delta = 2.0*ktmax*glen/(2.0*float(npx)+1.0)   	! grid step size in nm-1 
 call TransSpace(cell,k,kstar,'d','r')       		! transform incident direction to reciprocal space
 call CalcCross(cell,ga,kstar,gperp,'r','r',0)      	! compute g_perp = ga x k
 call NormVec(cell,gperp,'r')                       	! normalize g_perp
 call NormVec(cell,kstar,'r')                       	! normalize reciprocal beam vector

! allocate the head and tail of the linked list
 allocate(khead,stat=istat)   				! allocate new value
 if (istat.ne.0) call FatalError('Calckvectors','unable to allocate khead pointer')
 ktail => khead                      			! tail points to new value
 nullify(ktail%next)                			! nullify next in new value
 numk = 1                          			! keep track of number of k-vectors so far
 ktail%i = 0                             		! i-index of beam
 ktail%j = 0                             		! j-index of beam
 ktail%kt = (/0.0,0.0,0.0/)				! no tangential component for central beam direction
 ktail%k = kstar/mLambda				! divide by wavelength
 ktail%kn = CalcDot(cell,ktail%k,kstar,'r')			! normal component

! set the loop limits
 imin = -npx; imax = npx; jmin = -npy; jmax = npy; 

! and loop over the entire range (without symmetry considerations
 do i=imin,imax
  do j=jmin,jmax
   if (.not.((i.eq.0).and.(j.eq.0))) then  		! the point (0,0) has already been taken care of
    if ((i**2+j**2).le.ijmax) then 			! only directions inside the incident cone
     call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/)) 	! add k-vector to linked list
    end if
   end if
  end do
 end do

end if  ! mapmode = Conical


! standard or standard-conical kvector list, as used by CBED and other programs
if ( (mapmode.eq.'Standard').or.(mapmode.eq.'StandardConical') ) then

! for standard mode, we want to make sure that ijmax, which need not be defined by 
! the calling program for this mode, is set to a large value
  if (mapmode.eq.'Standard') then 
    ijmax = (5*npx)**2
  end if

! compute geometrical factors 
 glen = CalcLength(cell,ga,'r')              		! length of ga
 gan = ga/glen                                 	! normalized ga
 delta = 2.0*ktmax*glen/(2.0*float(npx)+1.0)   	! grid step size in nm-1 
 call TransSpace(cell,k,kstar,'d','r')       		! transform incident direction to reciprocal space
 call CalcCross(cell,ga,kstar,gperp,'r','r',0)      	! compute g_perp = ga x k
 call NormVec(cell,gperp,'r')                       	! normalize g_perp
 call NormVec(cell,kstar,'r')                       	! normalize reciprocal beam vector

! allocate the head and tail of the linked list
 allocate(khead,stat=istat)   				! allocate new value
 if (istat.ne.0) call FatalError('Calckvectors','unable to allocate khead pointer')
 ktail => khead                      			! tail points to new value
 nullify(ktail%next)                			! nullify next in new value
 numk = 1                          			! keep track of number of k-vectors so far
 ktail%i = 0                             		! i-index of beam
 ktail%j = 0                             		! j-index of beam
 ktail%kt = (/0.0,0.0,0.0/)				! no tangential component for central beam direction
 ktail%k = kstar/mLambda				! divide by wavelength
 ktail%kn = CalcDot(cell,ktail%k,kstar,'r')			! normal component

! implement symmetry Table 7.3 from CTEM book
  select case(isym)  ! negative values -> systematic row; positive -> zone axis
   case(-1)  ! centrosymmetric systematic row
     imin = 0; imax = npx; grid = 'srw'
   case(-2)  ! non-centrosymmetric systematic row
     imin = -npx; imax = npx; grid = 'srw'
   case(1)  ! 2D Group 1
     imin = -npx; imax = npx; jmin = -npy; jmax = npy; grid = 'sqa'
   case(2)  ! 2D Group 2
     imin = -npx; imax = npx; jmin = 0; jmax = npy; grid = 'sqb'
   case(3)  ! 2D Group m
     imin = -npx; imax = npx; jmin = 0; jmax = npy; grid = 'sqa'
   case(4)  ! 2D Group 2mm
     imin = 0; imax = npx; jmin = 0; jmax = npy; grid = 'sqa'
   case(5)  ! 2D Group 4
     imin = 1; imax = npx; jmin = 0; jmax = npy; grid = 'sqa'
   case(6)  ! 2D Group 4mm
     imin = 0; imax = npx; jmin = 0; jmax = npy; grid = 'sqc'
   case(7)  ! 2D Group 3   (cubic version)
     grid = 'hxa'; hexgrid=.TRUE.
   case(8)  ! 2D Group 31m  (cubic version)
     grid = 'hxb'; hexgrid=.TRUE.
   case(9)  ! 2D Group 6
     grid = 'hxe'; hexgrid=.TRUE.
   case(10)  ! 2D Group 6mm
     grid = 'hxf'; hexgrid=.TRUE.
   case(11)  ! 2D Group 3  (hexagonal setting)
     grid = 'hxc'; hexgrid=.TRUE.
   case(12)  ! 2D Group 31m (hexagonal setting)
     grid = 'hxd'; hexgrid=.TRUE.
   case(13)  ! 2D Group 3m1 (cubic setting)
     grid = 'hxg'; hexgrid=.TRUE.
   case(14)  ! 2D Group 3m1 (hexagonal setting)
     grid = 'hxh'; hexgrid=.TRUE.
   case default   ! we should never get here
     call FatalError('Calckvectors','unknown isym value')
  end select

! now do the real work for standard sets of wave vectors
  select case(grid)

   case('srw')		! systematic row incident beam orientations
     do i=imin,imax
      if (i.ne.0) then  				! the point (0,0) has already been taken care of
       call Add_knode(ktail,cell,i,0,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/)) 
      end if
     end do

   case('sqa')		! from here on, all orientations are zone axis cases for all Laue groups
     do i=imin,imax
      do j=jmin,jmax
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/)) 
        end if
       end if
      end do
     end do

   case('sqb')
     do i=imin,imax
      jloop_sqb:  do j=jmin,jmax
       if ((j.eq.0).and.(i.lt.0)) cycle jloop_sqb   	! skip the points  (i<0,0)
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/)) 
        end if
       end if
      end do jloop_sqb   
     end do

   case('sqc')
     do j=0,jmax
      do i=j,imax
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/)) 
        end if
       end if
      end do
     end do

   case('hxa')
     do j=0,npy
      do i=1-Kdelta(j,0),npx
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/),hexgrid) 
        end if
       end if
      end do
    end do

   case('hxb')
     do j=0,npy
      do i=j,npx
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/),hexgrid) 
        end if
       end if
      end do
    end do

   case('hxc')
     do j=0,npy
      do i=1-Kdelta(j,0)-j,npx-j
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/),hexgrid) 
        end if
       end if
      end do
     end do

   case('hxd')
     do j=0,npy
      do i=0,npx-j
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/),hexgrid) 
        end if
       end if
      end do
     end do

   case('hxe')
     do j=0,npy-1
      do i=1-Kdelta(j,0),npx-j
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/),hexgrid) 
       end if
      end if
      end do
     end do

   case('hxf')
     do j=0,npy/2
      do i=j,npx-j
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/),hexgrid) 
        end if
       end if
      end do
     end do

   case('hxg')
     do j=0,npy
      do i=j/2,min(2*j,npy)
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/),hexgrid) 
        end if
       end if
      end do
     end do

   case('hxh')
     do j=0,npy
      do i=-j/2,min(j,npy-1)
       if (.not.((i.eq.0).and.(j.eq.0))) then  	! the point (0,0) has already been taken care of
        if (i**2+j**2.le.ijmax) then
         call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/),hexgrid) 
        end if
       end if
      end do
     end do

  case default  ! we should never get here
    call FatalError('Calckvectors:','unknown grid type value')
    
  end select  ! grid value

end if ! mapmode.eq.'Standard' or 'StandardConical'

! the next type of grid is the one used for the modified Lambert maps in the dynamical EBSD 
! programs; this requires some special care, since these mappings are a little trickier than 
! those of the standard mapmode.  While it is possible to use a plain Lambert projection as
! well, here we only allow for the RoscaLambert mode.

if (mapmode.eq.'RoscaLambert') then 
   if (usehex) then  		! hexagonal grid step size
      delta =  2.D0*dsqrt(cPi)/3.0D0**0.75D0/dble(npx)
      hexgrid = .TRUE.
   else				! square grid step size
      delta = dsqrt(cPi*0.5D0)/dble(npx)
      hexgrid = .FALSE.
   end if

! allocate the head of the linked list
   allocate(khead,stat=istat)         		! allocate new value
   if (istat.ne.0) call FatalError('Calckvectors',' unable to allocate khead pointer')
   ktail => khead                      	! tail points to new value
   nullify(ktail%next)                		! nullify next in new value
   numk = 1                          		! keep track of number of k-vectors so far
   ktail%i = 0                             	! i-index of beam
   ktail%j = 0                             	! j-index of beam
   kstar = (/ 0.0, 0.0, 1.0 /)			! we always use c* as the center of the RoscaLambert projection
   call NormVec(cell,kstar,'c')                	! normalize incident direction
   kstar = kstar/mLambda           		! divide by wavelength
! and transform to reciprocal crystal space using the structure matrix
   ktail%k = matmul(transpose(cell%dsm),kstar)
   ktail%kn = 1.0/mLambda

! deal with each Laue group symmetry separately
 select case (isym)
 
   case (1)  ! triclinic symmetry
   	istart = -npx
	iend = npx
   	jstart = -npy
	jend = npy
	  do j=jstart,jend
	    do i=istart,iend   ! 
		call AddkVector(ktail,cell,numk,delta,i,j)
	    end do
	  end do

  case (2)   !  monoclinic symmetry
  	istart = -npx
	iend = npx
   	jstart = 0
	jend = npy
	  do j=jstart,jend
	   do i=istart,iend   ! 
		call AddkVector(ktail,cell,numk,delta,i,j)
	   end do
	  end do

  case (3,4,10)  ! orthorhombic mmm, tetragonal 4/m, cubic m-3
  	istart = 0
	iend = npx
   	jstart = 0
	jend = npy
	  do j=jstart,jend
	   do i=istart,iend   ! 
		call AddkVector(ktail,cell,numk,delta,i,j)
	   end do
	  end do

  case (5,11)  ! tetragonal 4/mmm, cubic m-3m
    	istart = 0
	iend = npx
   	jstart = 0
	jend = npy
	  do i=istart,iend
	   do j=jstart,i   ! 
		call AddkVector(ktail,cell,numk,delta,i,j)
	   end do
	  end do

  case (6,7)   ! npy is now npyhex !
   	istart = 0
	iend = npy
   	jstart = 0
	jend = npy
	  do j=jstart,jend
	    do i=istart,iend   ! 
		call AddkVector(ktail,cell,numk,delta,i,j,hexgrid)
	    end do
	  end do

  case (8)   ! npy is now npyhex !
   	istart = 0
	iend = npy
   	jstart = 0
	jend = npy
	  do j=jstart,jend
	    do i=j,iend   ! 
		call AddkVector(ktail,cell,numk,delta,i,j,hexgrid)
	    end do
	  end do

  case (9)   ! npy is now npyhex !
   	istart = 0
	iend = npx   ! was npy
   	jstart = 0
	jend = npx   ! was npy
	  do j=jstart,jend
	    do i=2*j,iend   ! 
		call AddkVector(ktail,cell,numk,delta,i,j,hexgrid)
	    end do
	  end do

  case (12)   ! npy is now npyhex !  This is the second setting of Laue group -3m.
   	istart = 0
	iend = npy
   	jstart = 0
	jend = npy/2
	  do j=jstart,jend
	    do i=2*j,iend   ! 
		call AddkVector(ktail,cell,numk,delta,i,j,hexgrid)
		call AddkVector(ktail,cell,numk,delta,i-j,-j,hexgrid)
	    end do
	  end do

 end select

end if

end subroutine Calckvectors


!--------------------------------------------------------------------------
!
! SUBROUTINE: CalckvectorsSymmetry
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create a linked list of wave vectors, using the whole pattern symmetry
!
!> @details This is a new version to test whether or not we can use the whole pattern
!> symmetry to determine the relevant list of incident wave vectors; this should be a 
!> general routine, so that we do not need to consider each symmetry case separately.
!> This will require a floating point version of the Apply2DPGSymmetry routine in symmetry.f90.
!
!> @param khead of head of linked list
!> @param cell unit cell pointer
!> @param TDPG 2D point group structure
!> @param k central wave vector
!> @param ga reciprocal lattice vector normal to k
!> @param ktmax maximum length of tangential component
!> @param npx number of vectors along x 
!> @param npy number of vectors along y 
!> @param numk total number of wave vectors in list
!> @param isym Laue group number
!> @param ijmax used for conical illumination
!> @param klaue fractional Laue center coordinates
!> @param debug when present, an output file with the kselected array is produced
!
!> @date   10/03/13 MDG 1.0 original
!> @date   06/09/14 MDG 2.0 added khead and cell arguments
!--------------------------------------------------------------------------
recursive subroutine CalckvectorsSymmetry(khead,cell,TDPG,k,ga,ktmax,npx,npy,numk,isym,ijmax,klaue,debug)

use io
use error
use constants
use diffraction
use crystal
use crystalvars
use Lambert

IMPLICIT NONE

type(kvectorlist),pointer              :: khead
type(unitcell),pointer	                :: cell
type(symdata2D),INTENT(INOUT)          :: TDPG
real(kind=dbl),INTENT(IN)		:: k(3)		!< initial wave vector
real(kind=dbl),INTENT(IN)		:: ga(3)	!< "horizontal" reciprocal lattice vector
real(kind=dbl),INTENT(IN)		:: ktmax	!< maximum length of tangential wave vector
integer(kind=irg),INTENT(IN)		:: npx		!< number of kvectors along x
integer(kind=irg),INTENT(IN)		:: npy		!< number of kvectors along y
integer(kind=irg),INTENT(OUT)		:: numk		!< total number of kvectors in linked list
integer(kind=irg),INTENT(IN)		:: isym		!< Laue symmetry group number 
integer(kind=irg),INTENT(INOUT)	:: ijmax	!< max parameter used for Conical and StandardConical modes
real(kind=sgl),INTENT(IN)		:: klaue(2)	!< fractional Laue center coordinates
logical,INTENT(IN),OPTIONAL		:: debug

integer(kind=irg),allocatable		:: kselected(:,:)	!< keeps track of which k-vectors have already been considered

integer(kind=irg)       		:: istat,i,j, iequiv(2,12), nequiv, jj, nx, ny
real(kind=dbl)				:: glen, gan(3), gperp(3), kstar(3), delta, Lauexy(2)
logical					:: hexgrid = .FALSE.
real(kind=sgl)				:: kt(3),kr(3)
real(kind=sgl)				:: ktlen
type(kvectorlist),pointer             :: ktail

nx = 2*npx
ny = 2*npy
allocate(kselected(-nx:nx,-ny:ny))

! initialize the kselected array to 0
kselected = 0

! compute geometrical factors 
 glen = CalcLength(cell,ga,'r')              		! length of ga
 Lauexy = glen * klaue					! scaled Laue center coordinates
 gan = ga/glen                                 	! normalized ga
 delta = 2.0*ktmax*glen/(2.0*float(npx)+1.0)   	! grid step size in nm-1 
 call TransSpace(cell,k,kstar,'d','r')       		! transform incident direction to reciprocal space
 call CalcCross(cell,ga,kstar,gperp,'r','r',0)      	! compute g_perp = ga x k
 call NormVec(cell,gperp,'r')                       	! normalize g_perp
 call NormVec(cell,kstar,'r')                       	! normalize reciprocal beam vector

! allocate the head and tail of the linked list
 allocate(khead,stat=istat)   				! allocate new value
 if (istat.ne.0) call FatalError('Calckvectors','unable to allocate khead pointer')
 ktail => khead                      			! tail points to new value
 nullify(ktail%next)                			! nullify next in new value
 numk = 1                          			! keep track of number of k-vectors so far
 ktail%i = 0                             		! i-index of beam
 ktail%j = 0                             		! j-index of beam
 
! use the Laue center coordinates to define the tangential component of the incident wave vector
 kt = - Lauexy(1)*gan - Lauexy(2)*gperp  		! tangential component of k
 ktail%kt = kt                    			! store tangential component of k
 ktlen = CalcLength(cell,kt,'r')**2      			! squared length of tangential component

 kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar 	! complete wave vector
 ktail%k = kr                     			! store in pointer list
 ktail%kn = CalcDot(cell,ktail%k,kstar,'r')    		! normal component of k
  
 kselected(0,0) = 2

 if (maxval(abs(klaue)).eq.0.0) then 			! zone axis orientation, so we should use symmetry
! we scan over the entire range of potential beam directions, defined by npx and npy along with
! the conical truncation parameter ijmax; for each point we check whether or not it has been considered
! before; it it has, we move on, if it hasn't, then we add this point to the linked list in the usual way.
! we do this by computing the equivalent (i,j) using the Whole Pattern symmetry.
     do i=-nx,nx
      do j=-ny,ny
  	if (kselected(i,j).eq.0) then
          if ((i*i+j*j).le.ijmax) then
! first of all, add the present point to the linked list
	   call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,(/ 0.0,0.0/)) 
! then compute the equivalent points and flag all of them in kselected
	   call Apply2DPGSymmetry(TDPG,i,j,isym,iequiv,nequiv)
	   kselected(iequiv(1,1),iequiv(2,1)) = 2
	   if (nequiv.gt.1) then 
	    do jj=2,nequiv
	     kselected(iequiv(1,jj),iequiv(2,jj)) = 1
	    end do
	   end if
        end if
       end if
      end do
     end do
  else							! not a zone axis, so no symmmetry
     do i=-nx,nx
      do j=-ny,ny
  	if (kselected(i,j).eq.0) then
          if ((i*i+j*j).le.ijmax) then
! first of all, add the present point to the linked list
	    call Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,sngl(Lauexy))
	    kselected(i,j) = 2
        end if
       end if
      end do
     end do
  end if

! for debugging purposes, we can write the kselected array to a file.
if (present(debug)) then
  open(unit=20,file='kselected.data',status='unknown',form='unformatted')
  write (20) 2*nx+1,2*ny+1
  write (20) kselected
  close(unit=20,status='keep')
end if

! and clean up the kselected array
deallocate(kselected)

end subroutine CalckvectorsSymmetry


!--------------------------------------------------------------------------
!
! SUBROUTINE: Add_knode
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief add one entry to the linked wave vector list (standard mode)
!
!> @param ktail current entry in linked list
!> @param cell unit cell pointer 
!> @param i x image coordinate 
!> @param j y image coordinate 
!> @param numk total number of wave vectors in list
!> @param delta scale parameter
!> @param gan normalized g-vector
!> @param gperp normalized perpendicular g-vector
!> @param kstar reciprocal components of wave vector
!> @param klaue fractional Laue center coordinates
!> @param hexgrid (optional) indicates hexagonal sampling if present
!
!> @todo implement Laue center coordinates for hexagonal grid
!
!> @date   04/29/13 MDG 1.0 original
!> @date   06/09/14 MDG 2.0 added ktail and cell as arguments
!--------------------------------------------------------------------------
subroutine Add_knode(ktail,cell,i,j,numk,delta,gan,gperp,kstar,klaue,hexgrid)

use local
use error
use crystal
use constants
use diffraction

IMPLICIT NONE

type(kvectorlist),pointer              :: ktail
type(unitcell),pointer	                :: cell
integer(kind=irg),INTENT(IN)		:: i
integer(kind=irg),INTENT(IN)		:: j
integer(kind=irg),INTENT(INOUT)	:: numk
real(kind=dbl),INTENT(IN)		:: delta
real(kind=dbl),INTENT(IN)		:: gan(3)
real(kind=dbl),INTENT(IN)		:: gperp(3)
real(kind=dbl),INTENT(IN)		:: kstar(3)
real(kind=sgl),INTENT(IN)		:: klaue(2)
logical,INTENT(IN),OPTIONAL		:: hexgrid

real(kind=sgl)				:: kt(3),kr(3)
real(kind=sgl)				:: ktlen
integer(kind=irg)			:: istat


allocate(ktail%next,stat=istat)  		! allocate new value
if (istat.ne.0) call FatalError('Add_knode:',' unable to allocate pointer')
ktail => ktail%next               		! tail points to new value
nullify(ktail%next)              		! nullify next in new value
numk = numk + 1                 		! keep track of number of k-vectors so far
ktail%i = i                      		! i-index of beam
ktail%j = j                      		! j-index of beam

! is it a square or hexagonal grid ?
if (present(hexgrid)) then
  kt = -(float(i)-float(j)*0.5)*delta*gan - float(j)*delta*gperp*0.5*sqrt(3.0)  ! tangential component of k
else
  kt = -(klaue(1)+float(i)*delta)*gan - (klaue(2)+float(j)*delta)*gperp  ! tangential component of k
end if

ktail%kt = kt                    		! store tangential component of k
ktlen = CalcLength(cell,kt,'r')**2      		! squared length of tangential component

kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar 	! complete wave vector
ktail%k = kr                     		! store in pointer list
ktail%kn = CalcDot(cell,ktail%k,kstar,'r')    	! normal component of k

end subroutine Add_knode


!--------------------------------------------------------------------------
!
! FUNCTION: GetSextant
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determines in which sextant a point (x,y) is located (used for RoscaLambert mapmode)
!
!> @param x x coordinate 
!> @param y y coordinate 
!
!> @date   11/21/12 MDG 1.0 original 
!> @date   04/29/13 MDG 1.1 modified for kvectors module
!--------------------------------------------------------------------------
function GetSextant(x,y) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN):: x, y 

real(kind=dbl),parameter	:: srt = 1.732050808   ! sqrt(3.D0)
integer(kind=irg)		:: res
real(kind=dbl)			:: xx

xx = dabs(x*srt)    	! |x| sqrt(3)

if (y.ge.0) then
  if (y.ge.xx) then
	res = 0
  else
	if (x.gt.0.D0) then
	  res = 1
	else
	  res = 5
	end if
  end if
else
  if (dabs(y).ge.xx) then
	res = 3
  else
	if (x.gt.0.D0) then
	  res = 2
	else
	  res = 4
	end if
  end if
end if

end function GetSextant


!--------------------------------------------------------------------------
!
! SUBROUTINE: Addkvector
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief add a k-vector for square or hexagonal grid sampling mode (used for RoscaLambert mapmode)
!
!> @todo replace the coordinate transformation with one defined in the Lambert module
!
!> @param ktail current entry in linked list
!> @param cell unit cell pointer
!> @param numk total number of wave vectors in list
!> @param delta scale parameter
!> @param i x image coordinate 
!> @param j y image coordinate 
!> @param usehex (optional) indicates hexagonal sampling if present
!
!> @date   11/21/12 MDG 1.0 original 
!> @date   04/29/13 MDG 1.1 modified for kvectors module
!> @date   06/09/14 MDG 2.0 added ktail as argument
!--------------------------------------------------------------------------
subroutine AddkVector(ktail,cell,numk,delta,i,j,usehex)

use local
use io
use constants
use error
use diffraction
use crystal
use crystalvars
use dynamical

IMPLICIT NONE

type(kvectorlist),pointer              :: ktail
type(unitcell),pointer	                :: cell
integer(kind=irg),INTENT(INOUT)	:: numk
real(kind=dbl),INTENT(IN)		:: delta
integer(kind=irg),INTENT(IN)		:: i
integer(kind=irg),INTENT(IN)		:: j
logical,INTENT(IN),OPTIONAL		:: usehex


integer(kind=irg)         		:: istat, ks
real(kind=dbl)                 	:: kstar(3), x, y, rr, q, iPi, XX, YY, xp, yp
logical 				:: goahead
real(kind=dbl),parameter		:: srt = 0.86602540D0   	! sqrt(3.D0)/2.D0
real(kind=dbl),parameter		:: isrt = 0.577350269D0   	! 1.D0/sqrt(3.D0)
real(kind=dbl),parameter		:: rtt = 1.7320508076D0   	!  sqrt(3)
real(kind=dbl),parameter		:: prea = 0.525037568D0   	!  3^(1/4)/sqrt(2pi)
real(kind=dbl),parameter		:: preb = 1.050075136D0   	!  3^(1/4)sqrt(2/pi)
real(kind=dbl),parameter		:: prec = 0.90689968D0   	!  pi/2sqrt(3)
real(kind=dbl),parameter		:: pred = 2.09439510D0   	!  2pi/3

! [06/09/14] all this needs to be replaced with calls to Lambert module !!!

! initalize some parameters
iPi = 1.D0/cPi  ! inverse of pi
goahead = .FALSE.

! hexagonal sampling or not?
if (present(usehex)) then
  x = (i - j*0.5)*delta
  y = j*delta*srt
else
  x = i*delta
  y = j*delta
end if

! r^2
rr = x*x+y*y
	
! this is the more correct mapping from a uniform square grid to a sphere via an equal area map
! to the 2D circle first; it is computationally slightly longer due to the trigonometric function calls
! but it is mathematically more correct.  We do distinguish here between the square and hexagon
! projections; Note that the hexagonal case must be mirrored x <-> y with respect to the 
! analytical derivation in the MSMSE paper due to a rotation of the hexagonal cell.

if ( .not.((i.eq.0).and.(j.eq.0)) ) then  ! skip (0,0) 
      if (usehex) then  ! we're projecting from a hexagonal array
! decide which sextant the point (i,j) is located in.
          ks = GetSextant(x,y)
	  select case (ks)
	  case (0,3)
	  	XX = preb*y*dcos(x*prec/y)
	  	YY = preb*y*dsin(x*prec/y)
	  case (1,4)
	  	xp = y+rtt*x
	  	yp = y*pred/xp
	  	XX = prea*xp*dsin(yp)
	  	YY = prea*xp*dcos(yp)
	  case (2,5)
	  	xp = y-rtt*x
	  	yp = y*pred/xp
	  	XX = prea*xp*dsin(yp)
	  	YY = -prea*xp*dcos(yp)	  
	  end select
	  q = XX**2+YY**2
	  kstar = (/ 0.5D0*XX*dsqrt(4.D0-q), 0.5D0*YY*dsqrt(4.D0-q),1.D0-0.5D0*q /)
          goahead = .TRUE.
      else   ! we're projecting from a square array
! decide which equation to use  [ (8) or (9) from Rosca's paper, with r=1 ]
	     if (dabs(x).le.dabs(y)) then
                 q = 2.D0*y*iPi*dsqrt(cPi-y*y)
  	          kstar = (/ q*dsin(x*cPi*0.25D0/y), q*dcos(x*cPi*0.25D0/y), 1.D0-2.D0*y*y*iPi /)  
	     else
                  q = 2.D0*x*iPi*dsqrt(cPi-x*x)
  	          kstar = (/ q*dcos(y*cPi*0.25D0/x), q*dsin(y*cPi*0.25D0/x), 1.D0-2.D0*x*x*iPi /)  
    	     end if
             goahead = .TRUE.
      end if
end if 
 
if (goahead) then 
     allocate(ktail%next,stat=istat)  			! allocate new value
     if (istat.ne.0) call FatalError('Addkvector:',' unable to allocate ktail pointer')
     ktail => ktail%next          			! tail points to new value
     nullify(ktail%next)          			! nullify next in new value
     numk = numk + 1       				! keep track of number of k-vectors so far
     if (usehex) then 					! transform the hex coordinates to square-array coordinates
       ktail%i = i - j/2+mod(j,2)/2           	! i-index of beam
       ktail%j = j                      		! j-index of beam
     else  						! leave the square coordinates unchanged
       ktail%i = i                      		! i-index of beam
       ktail%j = j                      		! j-index of beam
     end if 
     call NormVec(cell,kstar,'c')                		! normalize incident direction in cartesian space
     kstar = kstar/mLambda                 		! divide by wavelength
! and transform to reciprocal crystal space using the direct structure matrix
     ktail%k = matmul(transpose(cell%dsm),kstar)
     ktail%kn = 1.0/mLambda
end if

end subroutine AddkVector

!--------------------------------------------------------------------------
!
! SUBROUTINE: Delete_kvectorlist
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief delete the entire linked list
!
!> @date   04/29/13 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine Delete_kvectorlist(khead)

IMPLICIT NONE

type(kvectorlist),pointer        :: khead
type(kvectorlist),pointer        :: ktmp, ktail

! deallocate the entire linked list before returning, to prevent memory leaks
ktail => khead
ktmp => ktail % next
do 
  deallocate(ktail)
  if (.not. associated(ktmp)) EXIT
  ktail => ktmp
  ktmp => ktail % next
end do

end subroutine Delete_kvectorlist

end module kvectors
