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
! CTEMsoft2013:crystalvars.f90
!--------------------------------------------------------------------------
!
! MODULE: crystalvars
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Definition of all variables and types for crystallographic computations
!
!> @details  Defines the unitcell type and the orientation type, as well as 
!> the main cell variable used by all crystallographic computations
! 
!> @date   1/5/99   MDG 1.0 original
!> @date    7/16/99 MDG 1.1 added error handling and TransCoor
!> @date    4/ 5/00 MDG 1.2 modified TransCoor to include new mInvert
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   03/19/13 MDG 3.0 update to new version
!> @date   10/17/13 MDG 3.1 added HOLZentries type
!--------------------------------------------------------------------------
module crystalvars

use local

!
!> lattice parameters
!> a		= a parameter [nm]
!> b		= b parameter [nm]
!> c		= c parameter [nm]
!> alpha	= alpha angle [deg]
!> beta	= beta  angle [deg]
!> gamma	= gamma angle [deg]
!> vol		= unit cell volume [nm^3]
!
!> metric information
!> dmt	= direct space metric tensor
!> rmt		= reciprocal space metric tensor
!> dsm	= direct space structure matrix
!> rsm	= reciprocal space structure matrix
!> krdel	= Kronecker delta (unit matrix)
!
!
!> asymmetric unit contents
!> ATOM_ntype	= actual number of occupied positions in asymmetric unit
!> ATOM_type		= atomic number for each atom in asymmetric unit
!> ATOM_pos		= fractional coordinates (x,y,z), occupation, Debye-Waller factor for each atom in asymmetric unit
!> fname			= crystal structure file name
type unitcell
  real	(kind=dbl)	:: a,b,c,alpha,beta,gamma
  real	(kind=dbl)	:: dmt(3,3),rmt(3,3),dsm(3,3),rsm(3,3),krdel(3,3),vol
  integer(kind=irg)	:: ATOM_type(maxpasym),ATOM_ntype,SYM_SGnum,xtal_system,SYM_SGset
  real(kind=sgl)	:: ATOM_pos(maxpasym,5)
  character(fnlen)	:: fname
  logical			:: SYM_reduce,SYM_second,SYM_trigonal
end type

!> this type is used to define an orientation relation, i.e., two parallel
!> directions and two parallel planes
type orientation
  real(kind=sgl)  	:: tA(3), tB(3), gA(3), gB(3)
end type

!> cell is the generic unit cell variable used in all programs.  
type (unitcell) 	:: cell

! next we add a structure that contains all the relevant HOLZ geometry information;
! it can be filled by calling the GetHOLZGeometry subroutine in crystal.f90
type HOLZentries
  real(kind=sgl)	:: g1(3),g2(3),g3(3),gx(3),gy(3),LC1,LC2,H,FNr(3),FNg(3),gp(2),gtoc(2,2),gshort(3)
  integer(kind=irg)	:: uvw(3),FN(3)
end type

type (HOLZentries)	:: HOLZdata


end module
