
!--------------------------------------------------------------------------
! CTEMsoft2013:dynamical.f90
!--------------------------------------------------------------------------
!
! MODULE: dynamical
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief variables and types needed by diffraction module
! 
!> @todo some of these variables are obsolete, so this needs to be cleaned up
!
!> @date   10/13/98 MDG 1.0 original
!> @date    5/22/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date    3/14/02 MDG 2.2 added CalcDynMat routine
!> @date    4/29/13 MDG 3.0 rewrite, moved stuff into gvectors and kvectors modules
!> @date   11/05/13 MDG 3.1 added DynMats for multiphase calculations
!> @date   01/10/14 MDG 4.0 adjusted for new cell type definition
!--------------------------------------------------------------------------
module dynamical

use local

IMPLICIT NONE

! allocatable arrays
complex(kind=dbl),allocatable :: W(:), &         	! eigenvalue vector for Bloch wave method
                                 CG(:,:), &      	! eigenvector matrix
                                 alpha(:), &     	! excitation amplitude vector
                                 DHWMz(:,:),&		! Darwin-Howie-Whelan matrix
                                 DynMat(:,:), &  	! dynamical matrix
                                 DynMat0(:,:), &  	! dynamical matrix (for programs that need two or more of them)
                                 DynMat1(:,:), &  	! dynamical matrix (for programs that need two or more of them)
                                 DynMat2(:,:), &  	! dynamical matrix (for programs that need two or more of them)
                                 DynMat3(:,:), &  	! dynamical matrix (for programs that need two or more of them)
                                 phiz(:),Az(:,:) 	! used for Taylor expansion of scattering matrix

! The parameters in gnode are computed by CalcUcg 
type gnode
  character(2)         	:: method   ! computation method (WK = Weickenmeier-Kohl, DT = Doyle-Turner/Smith-Burge, XR for XRD)
  logical              	:: absorption ! is absorption included or not ?
  integer(kind=irg) 		:: hkl(3)   ! Miller indices
  real(kind=sgl)       	:: xg, &    ! extinction distance [nm]
                         	   xgp, &   ! absorption length [nm]
                          	   ar, &    ! aborption ratio
                          	   g, &     ! length of reciprocal lattice vectors [nm^-1]
                          	   Vmod,Vpmod, & ! modulus of Vg and Vgprime [V]
                          	   Umod,Upmod, & ! modulus of Ug and Ugprime [nm^-2]
                          	   Vphase,Vpphase ! phase factors of Vg and Vgprime [rad]
  complex(kind=sgl)		:: Ucg, &   ! scaled potential Fourier coefficient [nm^-2]
                          	   Vg, &    ! potential Fourier coefficient [V]
                          	   qg       ! interaction parameter for Darwin-Howie-Whelan equations [nm^-1]
end type gnode

type(gnode)            	:: rlp      ! reciprocal lattice point

! other vectors needed for dynamical computations
real(kind=sgl)   		:: DynWV(3), &       ! wave vector expressed in reciprocal frame
                    		DynFN(3), &       ! Foil normal in reciprocal frame
		    		DynUpz            ! U'_0 normal absorption parameter

! moved to cell type in crystalvars.f90
! integer(kind=irg)		:: DynNbeams, DynNbeamsLinked      ! number of beams

end module dynamical


