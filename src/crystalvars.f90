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
! CTEMsoft:crystalvars.f90
!--------------------------------------------------------------------------
!
! MODULE: crystalvars
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Definition of all variables and types for crystallographic computations
!
!> @details  Defines the unitcell type and the orientation type, as well as 
!> the main cell variable used by all crystallographic computations
! 
!> @date     1/5/99 MDG 1.0 original
!> @date    7/16/99 MDG 1.1 added error handling and TransCoor
!> @date    4/ 5/00 MDG 1.2 modified TransCoor to include new mInvert
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   03/19/13 MDG 3.0 update to new version
!> @date   10/17/13 MDG 3.1 added HOLZentries type
!> @date    1/10/14 MDG 4.0 new version, many new entries in unitcell type
!> @date    6/ 5/14 MDG 4.1 removed variable declaration for cell
!--------------------------------------------------------------------------
module crystalvars

use local
use symmetryvars

!> [note added on 1/10/14]
!> first we define the reflisttype (all the information needed for a given reciprocal lattice point or rlp).
!> This used to be in the gvectors module, but the type definitions make more sense here.
!> In this linked list, we want to keep everything that might be needed to perform rlp-related 
!> simulations, except for the Fourier coefficient of the lattice potential, wich is kept in 
!> a lookup table.  Anything that can easily be derived from the LUT does not need to be stored.
!> [end note]
!
! linked list of reflections
type reflisttype  
  integer(kind=irg)          	:: num, &  		! sequential number
                              	   hkl(3),&		! Miller indices
                              	   famhkl(3),&		! family representative Miller indices
				   HOLZN,& 		! belongs to this HOLZ layer
				   strongnum,& 		! sequential number for strong beams
				   weaknum,& 		! sequential number for weak beams
				   famnum		! family number
! removed 1/10/14		   nab(2)  		! decomposition with respect to ga and gb
  logical                    	:: dbdiff  		! double diffraction reflection ?
  real(kind=dbl)             	:: sg, &               ! excitation error
                                  xg, &              ! extinction distance
! removed 1/10/14                 Ucgmod, &          ! modulus of Fourier coefficient
                                  sangle, &          ! scattering angle (mrad)
                                  thetag             ! phase angle, needed for ECCI simulations
  logical	                :: strong, weak       ! is this a strong beam or not; both .FALSE. means 'do not consider'
  complex(kind=dbl)            :: Ucg	               ! Fourier coefficient, copied from cell%LUT
  type(reflisttype),pointer 	:: next    		! connection to next entry in master linked list
  type(reflisttype),pointer 	:: nexts    		! connection to next strong entry in linked list
  type(reflisttype),pointer 	:: nextw    		! connection to next weak entry in linked list
end type reflisttype


! a structure that contains all the relevant HOLZ geometry information;
! it can be filled by calling the GetHOLZGeometry subroutine in crystal.f90
type HOLZentries
  real(kind=sgl)	:: g1(3),g2(3),g3(3),gx(3),gy(3),LC1,LC2,H,FNr(3),FNg(3),gp(2),gtoc(2,2),gshort(3)
  integer(kind=irg)	:: uvw(3),FN(3)
end type

! this used to be a global variable declaration, but in Release 3.0 we no longer have any 
! global variables...
! type (HOLZentries)	:: HOLZdata



!> [note added on 1/10/14]
!> To make this package functional for multi-phase materials, we must have
!> a way to load the structural information for multiple crystal structures
!> simultaneously.  And we must also be able to perform scattering computations
!> for any given phase. This means that it is probably best if we define a 
!> pointer type to a large collection of unit cell related parameters, so that
!> we can easily switch from one to the other.  As a consequence, we need to 
!> pass the cell pointer along with the other arguments for every single routine
!> that performs simulations that require cell information.  So, first we define a 
!> pointer to a unitcell, and then we can fill in all the information needed,
!> including symmetry operators, atom coordinates, the potential coefficient lookup
!> table, pointers to the linked g-vector list, etc. In order to reduce the number
!> of changes to be made to the source code, we will keep the current variable names
!> as much as possible.
!> [end note]  
!
!> The following are the components of the unitcell type:
!
!> lattice parameters
!> a		= a parameter [nm]
!> b		= b parameter [nm]
!> c		= c parameter [nm]
!> alpha	= alpha angle [deg]
!> beta	        = beta angle [deg]
!> gamma	= gamma angle [deg]
!> vol		= unit cell volume [nm^3]
!
!> metric information
!> dmt	        = direct space metric tensor
!> rmt		= reciprocal space metric tensor
!> dsm	        = direct space structure matrix
!> rsm	        = reciprocal space structure matrix
!> [removed on 1/10/14] krdel	= Kronecker delta (unit matrix)
!
!> asymmetric unit contents
!> ATOM_ntype	= actual number of occupied positions in asymmetric unit
!> ATOM_type	= atomic number for each atom in asymmetric unit
!> ATOM_pos	= fractional coordinates (x,y,z), occupation, Debye-Waller factor for each atom in asymmetric unit
!
!> the structure file name
!> fname	= crystal structure file name
!
!> use hexagonal or regular indices (comes from old local.f90 module)
!> hexset	= logical to determine whether to use 3(FALSE) or 4(TRUE) index notation 
!
!> atom coordinate array
!> apos        = allocatable array for atom coordinates
!
!> storage space for the potential Fourier coefficient lookup table
!> LUT	        = lookup table (allocatable)
!
!> double diffraction logical array
!> dbdiff	= indicates whether a reflection could be due to double diffraction
!
!> is this space group non-symmorphic or not ?
!> nonsymmorphic = logical .TRUE. or .FALSE.
!
!> space group symmetry
!> SG          = space group symmetry structure defined in symmetryvars.f90
!
!> linked g-vector list (used to be in gvectors.f90)
!> reflist     = starting point of linked list of g-vectors
!
!> firstw	= pointer to first weak beam entry in list
!
!> number of beams in linked list (used to be in dynamical.f90)
!> DynNbeams   = current number being considered
!> DynNbeamsLinked = total number
!> nns	= number of strong beams
!> nnw = number of weak beams
type unitcell
  real(kind=dbl)	                :: a,b,c,alpha,beta,gamma
  real(kind=dbl)	                :: dmt(3,3),rmt(3,3),dsm(3,3),rsm(3,3),vol
  integer(kind=irg)	                :: ATOM_type(maxpasym),ATOM_ntype,SYM_SGnum,xtal_system,SYM_SGset
  real(kind=sgl)	                :: ATOM_pos(maxpasym,5)
  character(fnlen)	                :: fname
  logical				 :: hexset
  real(kind=dbl),allocatable           :: apos(:,:,:)
  complex(kind=dbl),allocatable        :: LUT(:,:,:)
  logical,allocatable                  :: dbdiff(:,:,:)
  logical                              :: nonsymmorphic
  type(symdata)                        :: SG
  type(reflisttype),pointer            :: reflist
  type(reflisttype),pointer 	        :: firstw    		! connection to first weak entry in linked list
  integer(kind=irg)                    :: DynNbeams, DynNbeamsLinked, nns, nnw
end type

!> this type is used to define an orientation relation, i.e., two parallel
!> directions and two parallel planes
type orientation
  real(kind=sgl)  	:: tA(3), tB(3), gA(3), gB(3)
end type


!> cell is a pointer to the generic unit cell variable used in all programs.  
! This pointer is allocated by the InitializeCell routine in the initializers.f90 module
! 
! in Release 2.0, the cell variable was a global variable.  In Release 3 and beyond,
! we aim to have no global variables at all and, instead, pass all variables as 
! function and subroutine arguments.
! type(unitcell) :: cell



end module
