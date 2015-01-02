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
! CTEMsoft:dictmod.f90
!--------------------------------------------------------------------------
!
! MODULE: dictmod
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Dictionary indexing routines
!
!> @details This module contains all the routines that deal with the dictionary
!> indexing approach, both in terms of computing all the dot products (which uses
!> OpenCL kernels) and in terms of the subsequent indexing on the symmetrized 
!> quaternion unit sphere using the modified von Mises-Fisher distribution.  All
!> the details for this approach can be found in two papers:
!>
!> "A dictionary based approach for EBSD indexing", Yu Hui Chen, Se Un Park, Dennis Wei, 
!> Greg Newstadt, Michael Jackson, Jeff Simmons, Alfred Hero, and Marc De Graef, 
!> Microscopy & Microanalysis, under review (2015).
!>
!> "Parameter estimation in spherical symmetry groups", Yu Hui Chen, Dennis Wei,
!> Gregory Newstadt, Marc De Graef, Jeff Simmons, and Al Hero, IEEE Signal Processing 
!> Letters, in print (2015)
!>
!> Here is an example program showing how the routines can be called:
!>
!> program t
!> 
!> use local
!> use typedefs
!> use dictmod
!> use quaternions
!> use constants
!> 
!> integer(kind=irg)               :: nums, seed
!> real(kind=dbl),allocatable      :: samples(:,:)
!> type(dicttype),pointer          :: dict
!> real(kind=dbl)                  :: q0, q1, q2, q3, muhat(4), kappahat, qu(4)
!> 
!> ! this is a test of the dictionary indexing portion that deals with the 
!> ! modified von Mises-Fisher distribution; the results must be the same 
!> ! as those produced by the original Matlab code...
!> 
!> seed = 432514
!> nums = 1000
!> allocate(samples(4,nums))
!> 
!> allocate(dict)
!> dict%Num_of_init = 3
!> dict%Num_of_iterations = 30
!> dict%pgnum = 32
!> 
!> open(UNIT=dataunit,file='../VMF/Matlab/samples.txt',status='old',form='formatted')
!> do i=1,nums
!>   read(dataunit,"(4E16.7)") q1,q2,q3,q0
!>   samples(1:4,i) = (/q0, q1, q2, q3 /)
!> end do
!> close(UNIT=dataunit,status='keep')
!> 
!> 
!> call InitDictionaryIndexing(dict)
!> 
!> do i=1,10
!>   call EMforVMF(samples, dict, nums, seed, muhat, kappahat)
!> 
!>   write (*,*) '  '
!>   write (*,*) 'mu    = ',muhat
!>   write (*,*) 'kappa = ',kappahat
!>   write (*,*) 'equivalent angular precision : ',180.D0*dacos(1.D0-1.D0/kappahat)/cPi
!> end do
!> 
!> end program
!> 

! 
!> @date 12/31/14 MDG 1.0 original (based on UMich Matlab code and IDL intermediate version)
!> @date 01/02/15 MDG 1.1 debug of code; produces same result as Matlab code
!--------------------------------------------------------------------------

module dictmod

IMPLICIT NONE

public  :: InitDictionaryIndexing, EMforVMF
private :: logCp, VMF_Estep, VMF_Mstep, VMFDensity,  VMF_getQandL

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: InitDictionaryIndexing
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief initialize the dictionary indexing parameters (symmetry operators and precomputed Ap lookup table)
!
!> @details For all details, see following paper:
!
!> @param dict dictionary parameter pointer (must be declared in calling routine)
!
!> @date 12/31/14 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine InitDictionaryIndexing(dict)
! recursive subroutine InitDictionaryIndexing(dict)

use local
use typedefs
use error
use math

IMPLICIT NONE

type(dicttype),pointer,INTENT(INOUT)    :: dict

integer(kind=irg)                       :: i
real(kind=dbl)                          :: y1, y2


! here we need to analyze the rotational symmetry group, and copy the appropriate 
! quaternion symmetry operators into the dict%Pm array

! first get the number of the rotational point group that corresponds to the crystal point group
dict%prot = PGrot(dict%pgnum)
! possible values for dict%prot are: (/1,3,6,9,12,16,18,21,24,28,30/)
! corresponding to the point groups 1, 2, 222, 4, 422, 3, 32, 6, 622, 23, and 432, respectively

!------------
! IMPORTANT NOTE: the original von Mises-Fischer (VMF) approach requires that q and -q are considered to 
! be separate quaternions, so the original Matlab code included the negatives of all quaternion symmetry operators 
! as well, leading to a cardinality of twice the rotational point group order.  It appears that we do not have to
! do so if we replace the exponential in the VMF by a hyperbolic cosine function, which would account directly
! for the q, -q duplicity...
!------------

! identity operator is part of all point groups
dict%Pm = 0.D0                  ! initialize all entries to zero
dict%Pm(1:4,1) = SYM_Qsymop(1:4,1)

! select statement for each individual rotational point group (see typedefs.f90 for SYM_Qsymop definitions)
select case (dict%prot) 
        case(1)         ! 1 (no additional symmetry elements)
                dict%Nqsym = 1*2
                dict%Pm(1:4,2) = -dict%Pm(1:4,1)

        case(3)         ! 2  (we'll assume that the two-fold axis lies along the e_y-axis)
                dict%Nqsym = 2*2
                dict%Pm(1:4,2) = SYM_Qsymop(1:4,3)
                do i=3,4
                  dict%Pm(1:4,i) = -dict%Pm(1:4,i-2)
                end do

        case(6)         ! 222
                dict%Nqsym = 4*2
                do i=2,4
                  dict%Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
                do i=5,8
                  dict%Pm(1:4,i) = -dict%Pm(1:4,i-4)
                end do

        case(9)         ! 4
                dict%Nqsym = 4*2
                dict%Pm(1:4,2) = SYM_Qsymop(1:4,4)
                dict%Pm(1:4,3) = SYM_Qsymop(1:4,7)
                dict%Pm(1:4,4) = SYM_Qsymop(1:4,10)
                do i=5,8
                  dict%Pm(1:4,i) = -dict%Pm(1:4,i-4)
                end do

        case(12)        ! 422
                dict%Nqsym = 8*2
                dict%Pm(1:4,2) = SYM_Qsymop(1:4,4)
                dict%Pm(1:4,3) = SYM_Qsymop(1:4,7)
                dict%Pm(1:4,4) = SYM_Qsymop(1:4,10)
                dict%Pm(1:4,5) = SYM_Qsymop(1:4,2)
                dict%Pm(1:4,6) = SYM_Qsymop(1:4,3)
                dict%Pm(1:4,7) = SYM_Qsymop(1:4,11)
                dict%Pm(1:4,8) = SYM_Qsymop(1:4,12)
                do i=9,16
                  dict%Pm(1:4,i) = -dict%Pm(1:4,i-8)
                end do

        case(16)        ! 3
                dict%Nqsym = 2
                call FatalError('InitDictionaryIndexing','this symmetry has not yet been implemented (pg 3)')

        case(18)        ! 32 (needs special handling)
                dict%Nqsym = 2
                call FatalError('InitDictionaryIndexing','this symmetry has not yet been implemented (pg 32)')

        case(21)        ! 6
                dict%Nqsym = 2
                call FatalError('InitDictionaryIndexing','this symmetry has not yet been implemented (pg 6)')

        case(24)        ! 622
                dict%Nqsym = 2
                call FatalError('InitDictionaryIndexing','this symmetry has not yet been implemented (pg 622)')

        case(28)        ! 23
                dict%Nqsym = 12*2
                do i=2,4
                  dict%Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
                do i=17,24
                  dict%Pm(1:4,4+(i-16)) = SYM_Qsymop(1:4,i)
                end do
                do i=13,24
                  dict%Pm(1:4,i) = -dict%Pm(1:4,i-12)
                end do

        case(30)        ! 432
                dict%Nqsym = 24*2
                do i=2,24
                  dict%Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
                do i=25,48
                  dict%Pm(1:4,i) = -dict%Pm(1:4,i-24)
                end do

        case default    ! this should never happen ...
                call FatalError('InitDictionaryIndexing','unknown rotational point group number')
end select


! the next part of the initial Matlab code computes a lookup table for the parameter Ap(u) (Appendix in paper)
dict%Apnum = 1400000
allocate(dict%xAp(dict%Apnum), dict%yAp(dict%Apnum))

! define the xAp array
dict%xAp = (/ (0.001D0+dble(i-1)*0.001D0,i=1,dict%Apnum)  /)
do i=1,dict%Apnum/2
  dict%yAp(i) = BesselIn(dict%xAp(i), 2) / BesselI1(dict%xAp(i))
end do

! do the rest via linear interpolation
y1 = dict%yAp(699000)   
y2 = dict%yAp(dict%Apnum/2) 
y1 = y2-y1
do i=1,dict%Apnum/2
  dict%yAp(i+dict%Apnum/2) = y2 + y1*(dict%xAp(i+dict%Apnum/2)-700.D0)
end do

end subroutine InitDictionaryIndexing

!--------------------------------------------------------------------------
!
! SUBROUTINE: EMforVMF
!
!> @author Yu-Hui Chen, U. Michigan / Marc De Graef, Carnegie Mellon University
!
!> @brief Expectation maximization approach to maximum likelihood problem for mu and kappa
!
!> @details For all details, see following paper:
!
!> @param X list of input quaternions
!> @param dict dictionary parameter pointer (must be declared in calling routine)
!> @param nums number of input quaternions
!> @param seed for normal random number generator 
!> @param muhat output mean orientation
!> @param kappahat output concentration parameter
!
!> @date 01/01/15 MDG 1.0 original, based on Chen's Matlab version + simplifications
!--------------------------------------------------------------------------
recursive subroutine EMforVMF(X, dict, nums, seed, muhat, kappahat)

use local
use constants
use typedefs
use math! , only:r8vec_normal_01, r4_uniform_01          ! array of normal random numbers
use quaternions
use rotations, only:qu2ro               ! we only need to move to Rodrigues-Frank space
use so3, only:IsinsideFZ                ! we only need to do a test ...

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: X(4,nums)
type(dicttype),pointer,INTENT(INOUT)    :: dict
integer(kind=irg),INTENT(IN)            :: nums
integer(kind=irg),INTENT(INOUT)         :: seed
real(kind=dbl),INTENT(OUT)              :: muhat(4)
real(kind=dbl),INTENT(OUT)              :: kappahat

integer(kind=irg)                       :: i, j, N, Pmdims, init, dd
!integer(kind=irg),parameter             :: Num_of_init=3                        ! preset based on cubic case testing
!integer(kind=irg),parameter             :: Num_of_iterations=30                 ! preset based on cubic case testing
integer(kind=irg)                       :: FZtype, FZorder
real(kind=dbl),allocatable              :: Mu_All(:,:), Kappa_All(:), R_All(:,:,:), L_All(:), &
                                           R(:,:), Q(:), L(:)
real(kind=dbl)                          :: Mu(4), PmMu(4), MuKa(5), Qi, Li, rod(4), qu(4), Kappa

! In this routine, we perform the EM algorithm to obtain an estimate for the 
! mean direction and concentration parameter of the modified von Mises-Fisher (mVMF)
! distribution that models the statistics of the orientation point cloud.
! We made one simplification with respect to the paper: instead of using twice
! the number of symmetry operators to account for the q/-q symmetry of the quaternion
! unit sphere, we combine both of them in one step by replacing the exponential in the 
! mVMF expression by a hyperbolic cosine.


! array sizes (we use shorthand notations)
N = nums
Pmdims = dict%Nqsym

! initialize some auxiliary arrays
allocate(Mu_All(dict%Num_of_init,4), Kappa_All(dict%Num_of_init), &
         R_All(N,Pmdims,dict%Num_of_init),L_All(dict%num_of_init))
Mu_All = 0.D0
Kappa_All = 0.D0
R_All = 0.D0
L_All = 0.D0


! main loop (EM typically uses a few starting parameter sets to make sure we don't get stuck in a local maximum)
do init=1,dict%Num_of_init 
! create a vector to hold the results
  allocate(R(N,Pmdims))
  R = 0.D0

! generate a normal random vector and normalize it as a starting guess for Mu (i.e., a unit quaternion)
  call R8VEC_normal_01(4,seed,Mu)
  Mu = Mu/cabs(Mu)
! Mu = (/ -0.144499, -0.560162, 0.803906, -0.138108 /)     ! this is a test vector for Mu

! the CTEMsoft package only considers quaternions with positive first component, 
! so we may need to change all the signs
  if (Mu(1).lt.0.D0) Mu = -Mu

! starting value for Kappa
  Kappa = 30.D0

! define the number of iterations and the Q and L function arrays
  allocate (Q(dict%Num_of_iterations), L(dict%Num_of_iterations))
  Q = 0.D0
  L = 0.D0


! and here we go with the EM iteration...
! we use quaternion multiplication throughout instead of the matrix version in the Matlab version
! quaternion multiplication has been verified against the 4x4 matrix multiplication of the atlab code on 01/02/15 
  iloop: do i=1,dict%Num_of_iterations 
! E-step
    R = VMF_Estep(X,dict,Pmdims,N,Mu,Kappa)

! M-step
    MuKa = VMF_Mstep(X,dict,Pmdims,N,R)

! calculate the Q and Likelihood function values
    call VMF_getQandL(X,dict,Pmdims,nums,MuKa,R,Qi,Li)
    L(i) = Li
    Q(i) = Qi

! update the containers
    Mu_All(init,1:4) = MuKa(1:4)
    Kappa_All(init) = MuKa(5)
    R_All(1:N,1:Pmdims,init) = R(1:N,1:Pmdims)
    L_All(init) = L(i)
    Mu = MuKa(1:4)
    Kappa = MuKa(5)

! and terminate if necessary
    if (i.ge.2) then 
      if (abs(Q(i)-Q(i-1)).lt.0.05) then 
        EXIT iloop
      end if
    end if
  end do iloop
  deallocate(R,Q,L)
end do

dd = maxloc(L_All,1)
Mu = Mu_all(dd,1:4)
kappahat = Kappa_All(dd)

write (*,*) ' All solutions : '
do i=1,dict%Num_of_init
  write (*,*) Mu_All(i,1:4), '---> ',L_All(i)
end do

! the CTEMsoft package only considers quaternions with positive first component, 
! so we may need to change all the signs
if (Mu(1).lt.0.D0) Mu = -Mu

! the final step is to make sure that the resulting Mu lies in the same
! fundamental zone that the dictionary elements are located in; since we start
! the EM iterations from a random quaternion, there is no guarantee that the 
! result lies in the same fundamental zone. Therefore, we cycle through all the 
! equivalent quaternions, and stop as soon as we find one in the Rodrigues 
! fundamental zone, which requires routines from the rotations and so3 modules. 
FZtype = FZtarray(dict%pgnum)
FZorder = FZoarray(dict%pgnum)

FZloop: do i=1,Pmdims
  qu = quat_mult(Mu,dict%Pm(1:4,i))
  if (qu(1).lt.0.D0) qu = -qu
  rod = qu2ro(qu)
  if (IsinsideFZ(rod,FZtype,FZorder)) EXIT FZloop
end do FZloop
muhat = qu

deallocate(Mu_All, Kappa_All, R_All, L_All)

end subroutine EMforVMF


!--------------------------------------------------------------------------
!
! FUNCTION: VMF_Estep
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief computes the E step of the EM process, verified against Matlab code on 01/02/15
!
!> @param X list of input quaternions
!> @param dict pointer to dictionary type
!> @param Pmdims number of quaternion symmetry operators
!> @param nums number of samples
!> @param Mu current guess for mean quaternion
!> @param Kappa input parameter
!
!> @date 01/01/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function VMF_Estep(X,dict,Pmdims,nums,Mu,Kappa) result(R)

use local
use typedefs
use quaternions

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: X(4,nums)
type(dicttype),pointer,INTENT(IN)       :: dict
integer(kind=irg),INTENT(IN)            :: Pmdims
integer(kind=irg),INTENT(IN)            :: nums
real(kind=dbl),INTENT(IN)               :: Mu(4)
real(kind=dbl),INTENT(IN)               :: Kappa
real(kind=dbl)                          :: R(nums,Pmdims)

integer(kind=irg)                       :: j 
real(kind=dbl)                          :: Rdenom(nums), PmMu(4)


! this implements equation (15) of the appendix of the paper
do j=1,Pmdims
  PmMu = quat_mult(Mu,dict%Pm(1:4,j))
  R(1:nums,j) = VMFDensity(X, nums, PmMu, Kappa)
end do

! and normalize
Rdenom = 1.D0/sum(R,2)
do j=1,Pmdims 
  R(1:nums,j) = R(1:nums,j)*Rdenom(1:nums)
end do

end function VMF_Estep


!--------------------------------------------------------------------------
!
! FUNCTION: VMF_Mstep
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief computes the M step of the EM process
!
!> @param X list of input quaternions
!> @param dict pointer to dictionary type
!> @param Pmdims number of quaternion symmetry operators
!> @param nums number of samples
!> @param R weight factors form the E step
!
!> @date 01/01/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function VMF_Mstep(X,dict,Pmdims,nums,R) result(MuKa)

use local
use typedefs
use quaternions

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: X(4,nums)
type(dicttype),pointer,INTENT(IN)       :: dict
integer(kind=irg),INTENT(IN)            :: Pmdims
integer(kind=irg),INTENT(IN)            :: nums
real(kind=dbl),INTENT(IN)               :: R(nums,Pmdims)
real(kind=dbl)                          :: MuKa(5)

real(kind=dbl)                          :: tmpGamma(4), nGamma, diff(dict%Apnum), qu(4)
integer(kind=irg)                       :: minp, i, j


! this is simplified from the Matlab routine and uses straight summations and
! quaternion multiplication instead of arrays
tmpGamma = 0.D0
do j=1,Pmdims 
  do i=1,nums
    qu = quat_mult(X(1:4,i),conjg(dict%Pm(1:4,j)))
    tmpGamma = tmpGamma +  R(i,j) * qu
  end do
end do
nGamma = cabs(tmpGamma)
MuKa(1:4) = tmpGamma/nGamma


! find kappa corresponding to this value of gamma (equation 17 in appendix of paper)
diff = dabs( nGamma/dble(nums) - dict%yAp ) 
minp = minloc( diff, 1 )
if (minp.eq.1) minp = 2 
MuKa(5) = dict%xAp(minp)

end function VMF_Mstep


!--------------------------------------------------------------------------
!
! SUBROUTINE: VMF_getQandL
!
!> @author Yu-Hui Chen, U. Michigan / Marc De Graef, Carnegie Mellon University
!
!> @brief Computes the Q array and the log-likelihood array
!
!> @details For all details, see following paper:
!
!> @param X list of input quaternions
!> @param dict dictionary parameter pointer (must be declared in calling routine)
!> @param Pmdims number of quaternion symmetry operators to consider
!> @param number of input quaternions
!> @param MuKa  vector with Mu and Kappa
!> @param R output from the E step
!> @param Q output Q 
!> @param L output L 
!
!> @date 01/01/15 MDG 1.0 original, based on Chen's Matlab version + simplifications
!--------------------------------------------------------------------------
recursive subroutine VMF_getQandL(X,dict,Pmdims,nums,MuKa,R,Q,L)

use local
use typedefs
use quaternions

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: X(4,nums)
type(dicttype),pointer,INTENT(INOUT)    :: dict
integer(kind=irg),INTENT(IN)            :: Pmdims
integer(kind=irg),INTENT(IN)            :: nums
real(kind=dbl),INTENT(IN)               :: MuKa(5)
real(kind=dbl),INTENT(IN)               :: R(nums,Pmdims)
real(kind=dbl),INTENT(OUT)              :: Q
real(kind=dbl),INTENT(OUT)              :: L

real(kind=dbl)                          :: Phi(nums,Pmdims), PmMu(4), qu(4)
integer(kind=irg)                       :: j

! compute the auxiliary Phi array
Phi = 0.D0
qu = MuKa(1:4)
do j=1,Pmdims
  PmMu = quat_mult(dict%Pm(1:4,j), qu)
  Phi(1:nums,j) = VMFDensity(X, nums, PmMu, MuKa(5))
end do
Phi = Phi/dble(dict%Nqsym)

! and convert the array into the Q and L parameters.
L = sum(dlog(sum(Phi,2)))
Q = sum(R*dlog(Phi))

end subroutine VMF_getQandL


!--------------------------------------------------------------------------
!
! FUNCTION: VMFDensity
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief computes the VMF density function, using a cosh() expression rather than exp()
!
!> @details function used by the VMFDensity function
!> original in Matlab by Yu-Hui Chen, U. Michigan
!> converted to IDL by MDG, 12/18/14, simplified arguments
!> converted to f90 by MDG, 12/31/14, further simplifications
!> output validated against Matlab output on 12/31/14
!
!> @param X input quaternion samples
!> @param nums number of samples
!> @param kappa input parameter
!
!> @date 01/01/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function VMFDensity(X,nums,mu,kappa) result(vmf)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: X(4,nums)
integer(kind=irg),INTENT(IN)            :: nums
real(kind=dbl),INTENT(IN)               :: mu(4)
real(kind=dbl),INTENT(IN)               :: kappa
real(kind=dbl)                          :: vmf(nums)

integer(kind=irg)                       :: j
real(kind=dbl)                          :: C

C = logCp(kappa)

do j=1,nums
! vmf(j) = 2.D0*C*dcosh(kappa*dot_product(mu,X(1:4,j)))
  vmf(j) = dexp(C+kappa*dot_product(mu,X(1:4,j)))
end do

end function VMFDensity

!--------------------------------------------------------------------------
!
! FUNCTION: logCp
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief computes the logarithm of Cp
!
!> @details function used by the VMFDensity function
!> original in Matlab by Yu-Hui Chen, U. Michigan
!> converted to IDL by MDG, 12/18/14, simplified arguments
!> converted to f90 by MDG, 12/31/14, further simplifications
!> output validated against Matlab output on 12/31/14
!
!> @param kappa input parameter
!
!> @date 12/31/14 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function logCp(kappa) result(lCp)

use local
use constants
use math 

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: kappa
real(kind=dbl)                  :: lCp

! pre-computed constants
real(kind=dbl),parameter        :: C=0.025330295911D0           ! C = 1.D0/(2.D0*cPi)**2
! y2 = dlog( C * ( 700.D0 / BesselI1(700.D0) ) ) 
! y1 = dlog( C * ( 699.D0 / BesselI1(699.D0) ) ) 
real(kind=dbl),parameter        :: y2my1=-0.997856378284D0, y2=-692.929658999631D0       ! y2my1 = y2-y1

! for arguments larger than kappa=700, we use a straight line interpolation, otherwise we compute the full function 
if (kappa.gt.700.D0) then 
  lCp = y2 + y2my1*(kappa-700.D0)
else 
  lCp = dlog( C * ( kappa / BesselI1(kappa) ) ) 
end if

end function logCp




end module dictmod
