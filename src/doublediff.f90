
!--------------------------------------------------------------------------
! CTEMsoft2013:doublediff.f90
!--------------------------------------------------------------------------
!
! MODULE: doublediff
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Anything related to double diffraction
! 
!> @details The following logical array is used to tag the potential double 
!> diffraction spots in non-symmorphic space groups.
!
!> @todo This may not need to be a separate module
!
!> @date   10/13/98 MDG 1.0 original
!> @date    5/22/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date    3/14/02 MDG 2.2 added CalcDynMat routine
!--------------------------------------------------------------------------
module doublediff

use local

logical,allocatable  	:: dbdiff(:)
logical              	:: nonsymmorphic

end module doublediff
