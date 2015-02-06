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
!> call DI_Init(dict)
!> 
!> do i=1,10
!>   call DI_EMforVMF(samples, dict, nums, seed, muhat, kappahat)
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
!> @date 01/04/15 MDG 1.2 trial implementation of model using hyperbolic functions instead of exponential
!> @date 01/06/15 MDG 1.3 changed public routine names with DI_ in front
!> @date 01/07/15 MDG 1.4 added VMF sampling routines 
!> @date 01/09/15 MDG 1.5 replaced several computations by numerically more stable versions
!> @date 02/05/15 MDG 1.6 added sampling for axial Watson distribution
!--------------------------------------------------------------------------

module dictmod

IMPLICIT NONE

public  :: DI_Init, DI_EMforVMF, VMFDensity, DI_Similarity_Classifier, logCp, DI_SampleVMF
private :: VMF_Estep, VMF_Mstep, VMF_getQandL, CardIntersection, randVMFMeanDir, randUniformSphere, &
           getVMFDensityLBM,  VMFMeanDirDensity

contains



!--------------------------------------------------------------------------
!
! FUNCTION: DI_SampleVMF
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief Sample the von Mises-Fisher distribution on the quaternion unit sphere
!
!> @param N number of samples to return
!> @param seed random number generator seed value
!> @param mu mean direction (unit quaternion)
!> @param kappa concentration
!
!> @date 01/07/15 MDG 1.0 original, based on Yu-Hui's Matlab code
!--------------------------------------------------------------------------
recursive function DI_SampleVMF(N, seed, mu, kappa) result(sVMF)

use local
use error

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: N
integer(kind=irg),INTENT(INOUT)         :: seed
real(kind=dbl),INTENT(IN)               :: mu(4)
real(kind=dbl),INTENT(IN)               :: kappa
real(kind=dbl)                          :: sVMF(4,N)

real(kind=dbl)                          :: nq, tmpmu(4), RandSphere(3,N), t(N), RS(4,N), lmu(4)
integer(kind=irg)                       :: i
! parameters for the singular value decomposition
integer(kind=irg)                       :: nr, LDA, LDU, LDVT, lwork, info
real(kind=dbl)                          :: mA(4,4), ss(4), u(4,4), vt, work(20)

! make sure the input quaternion is normalized
nq = dsqrt(sum(mu*mu))
if (nq.eq.0.D0) call FatalError('DI_SampleVMF','Input quaternion has zero length')
lmu = mu/nq

! initialize a bunch of parameters
tmpmu = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
sVMF = 0.D0
RS = 0.D0

! get the t-parameter and the distribution of random directions on the 2-sphere
t = randVMFMeanDir(N, kappa,seed)

RandSphere = randUniformSphere(N,seed)
RS(2:4,1:N) = RandSphere(1:3,1:N)

! merge these two parameters into the desired random variables
sVMF = transpose( spread(t,DIM=2,NCOPIES=4) * spread(tmpmu,DIM=1,NCOPIES=N) + &
                  spread(dsqrt(1.D0-t*t),DIM=2,NCOPIES=4) * transpose(RS) )

! Rotate the distribution along the desired mean direction mu
! In Matlab, one uses the null() operator which returns the null space of the argument
! This is then inserted into a 4x4 rotation matrix and multiplied with the quaternions
! from the random sample.  The null space of the input quaternion can be computed with
! singular value decomposition, which is done with the dgesvd Lapack routine. The matrix
! returned as u is the desired rotation matrix, except that the numbers in the first 
! column must have their signs reversed.

mA = 0.D0
mA(1:4,1) = lmu(1:4)
nr = 4
LDA = 4
LDVT = 1
LDU = 4
lwork = 20
call DGESVD('A','N',nr,nr,mA,LDA,ss, u, LDU, vt, LDVT, work, lwork, info)
u(1:4,1) = -u(1:4,1)

! next, apply this 4x4 rotation matrix to all of the generated quaternions to
! rotate them along the mean direction mu
do i=1,N
        sVMF(1:4,i) = matmul(u,sVMF(1:4,i))
end do

end function DI_SampleVMF


!--------------------------------------------------------------------------
!
! FUNCTION: randUniformSphere
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief 
!
!> @param N number of samples to return
!> @param seed random number generator seed value
!
!> @date 01/07/15 MDG 1.0 original, based on Yu-Hui's Matlab code, output transposed
!--------------------------------------------------------------------------
recursive function randUniformSphere(N,seed) result(ranSphere)

use local
use math

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: N
integer(kind=irg),INTENT(INOUT)         :: seed
real(kind=dbl)                          :: ranSphere(3,N)

real(kind=dbl)                          :: nq, NR(N*3), randNorm(3,N)
integer(kind=irg)                       :: i

ranSphere = 0.D0
call R8VEC_normal_01(N*3,seed,NR)
randNorm = reshape( NR, (/ 3, N /) )

! and normalize the three-vectors
do i=1,N
  nq = dsqrt(sum(randNorm(1:3,i)**2))
  RanSphere(1:3,i) = randNorm(1:3,i)/nq
end do

end function randUniformSphere


!--------------------------------------------------------------------------
!
! FUNCTION: randVMFMeanDir
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief 
!
!> @param N number of samples to return
!> @param k concentration
!> @param seed random number generator seed value
!
!> @date 01/07/15 MDG 1.0 original, based on Yu-Hui's Matlab code
!--------------------------------------------------------------------------
recursive function randVMFMeanDir(N, k, seed) result(t)

use local
use math

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: N
real(kind=dbl),INTENT(IN)               :: k
integer(kind=irg),INTENT(INOUT)         :: seed
real(kind=dbl)                          :: t(N)

real(kind=dbl)                          :: LBM(2), h, x, C
integer(kind=irg)                       :: i

! Find the left bound and maximum
C = 0.D0
LBM = getVMFDensityLBM(k,C)

! apply the rejection sampling algorithm
t = 0.D0
do i=1,N
  do 
    x = r8_uniform_01(seed)*(1.D0-LBM(1))+LBM(1)
    h = VMFMeanDirDensity(x, k, C)
    if (r8_uniform_01(seed)*LBM(2).le.h) EXIT 
  end do
  t(i) = x
end do 

end function randVMFMeanDir

!--------------------------------------------------------------------------
!
! FUNCTION: getVMFDensityLBM
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief determines the left bound and maximum for rejection sampling
!
!> @param k concentration
!> @param C constant prefactor of distribution function
!
!> @date 01/07/15 MDG 1.0 original, based on Yu-Hui's Matlab code
!--------------------------------------------------------------------------
recursive function getVMFDensityLBM(k,C) result(LBM)

use local
use error

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: k
real(kind=dbl),INTENT(INOUT)            :: C
real(kind=dbl)                          :: LBM(2)

real(kind=dbl),parameter                :: min_thresh=0.00001D0
real(kind=dbl)                          :: s, x
integer(kind=irg)                       :: f

! first we look for the left bound
f = 1
do 
   x = -1.D0+dble(f)*0.00001D0
   if (x.eq.1.D0) call FatalError('getVMFDensityLBM','reached +1 in leftbound determination')
   s = VMFMeanDirDensity(x,k,C)
   if (s.ge.min_thresh) EXIT
   f = f+1
end do
!
LBM(1) =  -1.D0 +dble(f)*0.00001D0

! for the simplified version of the density function, we have an analytical
! expression for where the maximum of the function occurs [convert the BesselI(3/2,x)
! to hyperbolic functions, then to exponential, and ignore the negative exponential
! which will be very small for reasonably sized k and t...]
x = (-1.D0+dsqrt(1.D0+4.D0*k*k))/(2.D0*k)

LBM(2) = VMFMeanDirDensity(x,k,C)


end function getVMFDensityLBM

!--------------------------------------------------------------------------
!
! FUNCTION: VMFMeanDirDensity
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief function to be sampled for VMF random sampling; we're using a close approximation
!
!> @param x argument value
!> @param k concentration
!> @param C constant prefactor
!
!> @date 01/07/15 MDG 1.0 original, based on Yu-Hui's Matlab code
!--------------------------------------------------------------------------
recursive function VMFMeanDirDensity(x, k, C) result(y)

use local
use constants
use math
use error

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: x
real(kind=dbl),INTENT(IN)       :: k
real(kind=dbl),INTENT(INOUT)    :: C
real(kind=dbl)                  :: y


if (dabs(x).gt.1.D0) call FatalError('VMFMeanDirDensity','argument must be in [-1,1]')

! explicit expression for p=4 (Gamma[3/2]Gamma[1/2] = pi/2)
! and the BesselI(3/2) function reduces to hyperbolic functions
! diverges for k->0, and becomes really small for large k
if (C.eq.0.D0) then
  C = 2.D0*k**(2.5D0)/dsqrt(2.D0*cPi)/(k-1.D0)
end if

! this is a close approximation, really good for larger values of k
! and numerically more stable than the original, which has problems for k>600 or so
y = C * dexp(k*(x-1.D0))*dsqrt(1.D0-x*x)

end function VMFMeanDirDensity

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! sampling routines for the axial Watson distribution
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! FUNCTION: DI_SampleWatson
!
!> @author Marc De Graef, Carnegie Mellon University 
!
!> @brief Sample the axial Watson distribution on the quaternion unit sphere
!
!> @param N number of samples to return
!> @param seed random number generator seed value
!> @param mu mean direction (unit quaternion)
!> @param kappa concentration
!
!> @date 02/05/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function DI_SampleWatson(N, seed, mu, kappa) result(sWatson)

use local
use error

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: N
integer(kind=irg),INTENT(INOUT)         :: seed
real(kind=dbl),INTENT(IN)               :: mu(4)
real(kind=dbl),INTENT(IN)               :: kappa
real(kind=dbl)                          :: sWatson(4,N)

real(kind=dbl)                          :: nq, tmpmu(4), RandSphere(3,N), t(N), RS(4,N), lmu(4)
integer(kind=irg)                       :: i

! parameters for the singular value decomposition
integer(kind=irg)                       :: nr, LDA, LDU, LDVT, lwork, info
real(kind=dbl)                          :: mA(4,4), ss(4), u(4,4), vt, work(20)

! make sure the input quaternion is normalized
nq = dsqrt(sum(mu*mu))
if (nq.eq.0.D0) call FatalError('DI_SampleWatson','Input quaternion has zero length')
lmu = mu/nq

! initialize a bunch of parameters
tmpmu = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
sWatson = 0.D0
RS = 0.D0

! get the t-parameter and the distribution of random directions on the 2-sphere
t = randWatsonMeanDir(N, kappa,seed)

RandSphere = randUniformSphere(N,seed)
RS(2:4,1:N) = RandSphere(1:3,1:N)

! merge these two parameters into the desired random variables
sWatson = transpose( spread(t,DIM=2,NCOPIES=4) * spread(tmpmu,DIM=1,NCOPIES=N) + &
                  spread(dsqrt(1.D0-t*t),DIM=2,NCOPIES=4) * transpose(RS) )

! Rotate the distribution along the desired mean direction mu
! In Matlab, one uses the null() operator which returns the null space of the argument
! This is then inserted into a 4x4 rotation matrix and multiplied with the quaternions
! from the random sample.  The null space of the input quaternion can be computed with
! singular value decomposition, which is done with the dgesvd Lapack routine. The matrix
! returned as u is the desired rotation matrix, except that the numbers in the first 
! column must have their signs reversed.

mA = 0.D0
mA(1:4,1) = lmu(1:4)
nr = 4
LDA = 4
LDVT = 1
LDU = 4
lwork = 20
call DGESVD('A','N',nr,nr,mA,LDA,ss, u, LDU, vt, LDVT, work, lwork, info)
u(1:4,1) = -u(1:4,1)

! next, apply this 4x4 rotation matrix to all of the generated quaternions to
! rotate them along the mean direction mu
do i=1,N
        sWatson(1:4,i) = matmul(u,sWatson(1:4,i))
end do

end function DI_SampleWatson

!--------------------------------------------------------------------------
!
! FUNCTION: randWatsonMeanDir
!
!> @author Marc De Graef, Carnegie Mellon University 
!
!> @brief 
!
!> @param N number of samples to return
!> @param k concentration
!> @param seed random number generator seed value
!
!> @date 02/05/15 MDG 1.0 original 
!--------------------------------------------------------------------------
recursive function randWatsonMeanDir(N, k, seed) result(t)

use local
use math

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: N
real(kind=dbl),INTENT(IN)               :: k
integer(kind=irg),INTENT(INOUT)         :: seed
real(kind=dbl)                          :: t(N)

real(kind=dbl)                          :: LBM(2), h, x, C
integer(kind=irg)                       :: i

! Find the left bound and maximum
C = 0.D0
LBM = getWatsonDensityLBM(k,C)

! apply the rejection sampling algorithm
t = 0.D0
do i=1,N
  do 
    x = r8_uniform_01(seed)*(1.D0-LBM(1))+LBM(1)
    h = WatsonMeanDirDensity(x, k, C)
    if (r8_uniform_01(seed)*LBM(2).le.h) EXIT 
  end do
  t(i) = x
end do 

end function randWatsonMeanDir

!--------------------------------------------------------------------------
!
! FUNCTION: getWatsonDensityLBM
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief determines the left bound and maximum for rejection sampling
!
!> @param k concentration
!> @param C constant prefactor of distribution function
!
!> @date 01/07/15 MDG 1.0 original, based on Yu-Hui's Matlab code
!--------------------------------------------------------------------------
recursive function getWatsonDensityLBM(k,C) result(LBM)

use local
use error

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: k
real(kind=dbl),INTENT(INOUT)            :: C
real(kind=dbl)                          :: LBM(2)

real(kind=dbl),parameter                :: min_thresh=0.00001D0
real(kind=dbl)                          :: s, x
integer(kind=irg)                       :: f

! first we look for the left bound
f = 1
do 
   x = dble(f)*0.00001D0
   if (x.eq.1.D0) call FatalError('getWatsonDensityLBM','reached +1 in leftbound determination')
   s = WatsonMeanDirDensity(x,k,C)
   if (s.ge.min_thresh) EXIT
   f = f+1
end do
!
LBM(1) = dble(f)*0.00001D0

! for the simplified version of the density function, we have an analytical
! expression for where the maximum of the function occurs 
x = dsqrt((2.D0*k-1.D0)/(2.D0*k))

LBM(2) = WatsonMeanDirDensity(x,k,C)

write (*,*) ' LBM = ',LBM

end function getWatsonDensityLBM


!--------------------------------------------------------------------------
!
! FUNCTION: WatsonMeanDirDensity
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief function to be sampled for Watson random sampling; we're using a close approximation
!
!> @param x argument value
!> @param k concentration
!> @param C constant prefactor
!
!> @date 02/05/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function WatsonMeanDirDensity(x, k, C) result(y)

use local
use constants
use math
use error

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: x
real(kind=dbl),INTENT(IN)       :: k
real(kind=dbl),INTENT(INOUT)    :: C
real(kind=dbl)                  :: y

real(kind=dbl),parameter        :: CC = 144.43253338822560946D0         ! 256/sqrt(pi)


if (dabs(x).gt.1.D0) call FatalError('WatsonMeanDirDensity','argument must be in [-1,1]')

! approximate expression for p=4 
if (C.eq.0.D0) then
  C = CC*k**4.5D0/(525.D0+4.D0*k*(45.D0+8.D0*k*(3.D0+4.D0*k)))
end if

! this is a close approximation, really good for larger values of k
! and numerically more stable than the original, which has problems for k>600 or so
y = C * dexp(k*(x*x-1.D0))*dsqrt(1.D0-x*x)

end function WatsonMeanDirDensity



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! below we have a series of functions/subroutines used to index patterns based
! on the set of 40 or so closest matches in the dictionary; these routines perform
! an averaging on the quaternion unit sphere and return the mean direction and
! the concentration parameter.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
! SUBROUTINE: DI_Init
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief initialize the dictionary indexing parameters (symmetry operators and precomputed Ap lookup table)
!
!> @details For all details, see following paper:
!
!> @param dict dictionary parameter pointer (must be declared in calling routine)
!> @param fmode  'VMF' or 'WAT' for von Mises-Fisher and Watson distributions, respectively.
!> @param full optional, if present duplicate all operators for +q and -q
!
!> @date 12/31/14 MDG 1.0 original
!> @date 01/06/15 MDG 1.1 added optional argument full
!--------------------------------------------------------------------------
recursive subroutine DI_Init(dict,fmode,full)

use local
use typedefs
use error
use math

IMPLICIT NONE

type(dicttype),pointer,INTENT(INOUT)    :: dict
character(3),INTENT(IN)                 :: fmode
logical,INTENT(IN),OPTIONAL             :: full

integer(kind=irg)                       :: i, M
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
! for the q, -q duplicity... This is taken care of by the optional argument full
!------------

M=1
dict%full = .FALSE.

if (present(full)) then
  M=2
  dict%full = .TRUE.
end if

! identity operator is part of all point groups
dict%Pm = 0.D0                  ! initialize all entries to zero
dict%Pm(1:4,1) = SYM_Qsymop(1:4,1)

! select statement for each individual rotational point group (see typedefs.f90 for SYM_Qsymop definitions)
select case (dict%prot) 
        case(1)         ! 1 (no additional symmetry elements)
                dict%Nqsym = 1*M
                dict%Pm(1:4,2) = -dict%Pm(1:4,1)

        case(3)         ! 2  (we'll assume that the two-fold axis lies along the e_y-axis)
                dict%Nqsym = 2*M
                dict%Pm(1:4,2) = SYM_Qsymop(1:4,3)
                if (M.eq.2) then 
                        do i=3,4
                                dict%Pm(1:4,i) = -dict%Pm(1:4,i-2)
                        end do
                end if

        case(6)         ! 222
                dict%Nqsym = 4*M
                do i=2,4
                  dict%Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
                if (M.eq.2) then 
                        do i=5,8
                                dict%Pm(1:4,i) = -dict%Pm(1:4,i-4)
                        end do
                end if

        case(9)         ! 4
                dict%Nqsym = 4*M
                dict%Pm(1:4,2) = SYM_Qsymop(1:4,4)
                dict%Pm(1:4,3) = SYM_Qsymop(1:4,7)
                dict%Pm(1:4,4) = SYM_Qsymop(1:4,10)
                if (M.eq.2) then 
                        do i=5,8
                                dict%Pm(1:4,i) = -dict%Pm(1:4,i-4)
                        end do
                end if

        case(12)        ! 422
                dict%Nqsym = 8*M
                dict%Pm(1:4,2) = SYM_Qsymop(1:4,4)
                dict%Pm(1:4,3) = SYM_Qsymop(1:4,7)
                dict%Pm(1:4,4) = SYM_Qsymop(1:4,10)
                dict%Pm(1:4,5) = SYM_Qsymop(1:4,2)
                dict%Pm(1:4,6) = SYM_Qsymop(1:4,3)
                dict%Pm(1:4,7) = SYM_Qsymop(1:4,11)
                dict%Pm(1:4,8) = SYM_Qsymop(1:4,12)
                if (M.eq.2) then 
                        do i=9,16
                                dict%Pm(1:4,i) = -dict%Pm(1:4,i-8)
                        end do
                end if

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
                dict%Nqsym = 12*M
                do i=2,4
                  dict%Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
                do i=17,24
                  dict%Pm(1:4,4+(i-16)) = SYM_Qsymop(1:4,i)
                end do
                if (M.eq.2) then 
                        do i=13,24
                                dict%Pm(1:4,i) = -dict%Pm(1:4,i-12)
                        end do
                end if

        case(30)        ! 432
                dict%Nqsym = 24*M
                do i=2,24
                  dict%Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
                if (M.eq.2) then 
                        do i=25,48
                                dict%Pm(1:4,i) = -dict%Pm(1:4,i-24)
                        end do
                end if

        case default    ! this should never happen ...
                call FatalError('InitDictionaryIndexing','unknown rotational point group number')
end select

! von Mises-Fisher mode:
! the next part of the initial Matlab code computes a lookup table for the parameter Ap(u) (Appendix in paper)
! this lookup table is only used when the ratio of the BesselI functions is between 0 and 0.95; for the 
! region between 0.95 and 1, we use an analytical approximation (see VMF_Mstep routine).
!
! Watson mode:
! we've used a similar approach to create a lookup table for values of kappa that are smaller than 35, in
! which case we use the standard ration of Kummer functions:  Kummer[3/2,3,k]/Kummer[1/2,2,k]/k.  For
! larger kappa values, we have an expansion using the large argument behavior of the modified Bessel functions.
! 
dict%Apnum = 35000
allocate(dict%xAp(dict%Apnum), dict%yAp(dict%Apnum))

  ! define the xAp array
dict%xAp = (/ (0.001D0+dble(i-1)*0.001D0,i=1,dict%Apnum)  /)

if (fmode.eq.'VMF') then ! von Mises-Fisher distribution
  do i=1,dict%Apnum
    dict%yAp(i) = BesselIn(dict%xAp(i), 2) / BesselI1(dict%xAp(i))
  end do
else ! Watson distribution
  do i=1,dict%Apnum
    y1 = BesselI1(dict%xAp(i)*0.5D0)
    y2 = BesselI0(dict%xAp(i)*0.5D0)
    dict%yAp(i) = y1 / (y2-y1) / dict%xAp(i)
  end do
end if

end subroutine DI_Init

!--------------------------------------------------------------------------
!
! SUBROUTINE: DI_EMforVMF
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
!> @date 01/06/15 MDG 1.1 added optional argument full
!--------------------------------------------------------------------------
recursive subroutine DI_EMforVMF(X, dict, nums, seed, muhat, kappahat)

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
integer(kind=irg)                       :: FZtype, FZorder
real(kind=dbl),allocatable              :: Mu_All(:,:), Kappa_All(:), R_All(:,:,:), L_All(:), &
                                           R(:,:), Q(:), L(:)
real(kind=dbl)                          :: Mu(4), PmMu(4), MuKa(5), Qi, Li, rod(4), qu(4), Kappa

! In this routine, we perform the EM algorithm to obtain an estimate for the 
! mean direction and concentration parameter of the modified von Mises-Fisher (mVMF)
! distribution that models the statistics of the orientation point cloud.

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
! quaternion multiplication has been verified against the 4x4 matrix multiplication of the Matlab code on 01/02/15 
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
      if (abs(Q(i)-Q(i-1)).lt.0.01) then 
        EXIT iloop
      end if
    end if
  end do iloop
  deallocate(R,Q,L)
end do

dd = maxloc(L_All,1)
Mu = Mu_all(dd,1:4)
kappahat = Kappa_All(dd)

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

end subroutine DI_EMforVMF


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
!> @date 01/06/15 MDG 1.1 added optional argument full
!> @date 01/09/15 MDG 1.1 removed optional argument full (incorporated in dicttype)
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

integer(kind=irg)                       :: j, i 
real(kind=dbl)                          :: Rdenom(nums), PmMu(4), coshsum(nums), arg(nums), qu(4), C


if (dict%full) then
  C = logCp(kappa)
! this implements equation (15) of the appendix of the paper
  do j=1,Pmdims
    PmMu = quat_mult(Mu,dict%Pm(1:4,j))
    R(1:nums,j) = VMFDensity(X, nums, PmMu, Kappa, C)
  end do
! and determine the normalization factors
  Rdenom = 1.D0/sum(R,2)
else
  coshsum = 0.D0
! this implements equation (15) of the appendix of the paper
  do j=1,Pmdims
    PmMu =  quat_mult(Mu,dict%Pm(1:4,j)) 
    do i=1,nums
      arg(i) = dot_product(PmMu,X(1:4,i))
    end do
    coshsum = coshsum + dcosh(Kappa*arg)
    R(1:nums,j) = dsinh(Kappa*arg)
  end do
! and determine the normalization factors
  Rdenom = 1.D0/coshsum
end if 

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
!> @date 01/06/15 MDG 1.1 added optional full parameter
!> @date 01/09/15 MDG 1.2 removed optional argument full (incorporated in dicttype)
!> @date 01/09/15 MDG 1.3 introduced accurate numerical approximation for kappa determination
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

real(kind=dbl)                          :: tmpGamma(4), nGamma, diff(dict%Apnum), qu(4), y
integer(kind=irg)                       :: minp, i, j


if (dict%full) then
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
else
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
end if

! find kappa corresponding to this value of gamma (equation 17 in appendix of paper)
! we split this into two regionds: 0<=y<0.95, for which we use the look-up table 
! approach, and 0.95<=y<=1, for which we have derived an analytical approximation
! that is pretty accurate in the relevant region of kappa>30.
y = nGamma/dble(nums)
if (y.ge.0.95D0) then
  MuKa(5) = (15.D0-3.D0*y+dsqrt(15.D0+90.D0*y+39.D0*y*y))/(16.D0*(1.0D0-y))
else
  diff = dabs( y - dict%yAp ) 
  minp = minloc( diff, 1 )
  if (minp.eq.1) minp = 2 
  MuKa(5) = dict%xAp(minp)
end if

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
!> @date 01.06/15 MDG 1.1 added optional full parameter
!> @date 01/09/15 MDG 1.1 removed optional argument full (incorporated in dicttype)
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

real(kind=dbl)                          :: Phi(nums,Pmdims), PmMu(4), qu(4), C
integer(kind=irg)                       :: j
real(kind=dbl),parameter                :: eps = 0.001D0

if (dict%full) then 
  C = logCp(MuKa(5))
! compute the auxiliary Phi array
  Phi = 0.D0
  qu = MuKa(1:4)
  do j=1,Pmdims
    PmMu = quat_mult(dict%Pm(1:4,j), qu)
    Phi(1:nums,j) = VMFDensity(X, nums, PmMu, MuKa(5), C)
  end do
  Phi = Phi/dble(dict%Nqsym)

! and convert the array into the Q and L parameters.
  L = sum(dlog(sum(Phi,2)))
  Q = sum(R*dlog(Phi+eps))
else
  C = dexp(logCp(MuKa(5)))
! compute the auxiliary Phi array
  Phi = 0.D0
  qu = MuKa(1:4)
  do j=1,Pmdims
    PmMu = quat_mult(dict%Pm(1:4,j), qu)
    Phi(1:nums,j) = VMFDensity(X, nums, PmMu, MuKa(5),C)
  end do
  Phi = Phi/dble(dict%Nqsym)

! and convert the array into the Q and L parameters.
  L = sum(dlog(sum(Phi,2)))
  Q = sum(R*dlog(Phi))
end if

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
!> @param mu mean direction
!> @param kappa concentration
!> @param C logCp(kappa) or exp(logCp(kappa) (precomputed in calling routine)
!
!> @date 01/01/15 MDG 1.0 original
!> @date 01/06/15 MDG 1.1 added C and full parameters
!--------------------------------------------------------------------------
recursive function VMFDensity(X,nums,mu,kappa,C) result(vmf)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: X(4,nums)
integer(kind=irg),INTENT(IN)            :: nums
real(kind=dbl),INTENT(IN)               :: mu(4)
real(kind=dbl),INTENT(IN)               :: kappa
real(kind=dbl),INTENT(IN)               :: C
real(kind=dbl)                          :: vmf(nums)

integer(kind=irg)                       :: j

do j=1,nums
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
!> @date 01/09/14 MDG 1.1 introduced more accurate numerical approximation for Cp
!> @date 01/09/14 MDG 1.2 moved some of the constants in front of the logarithm
!--------------------------------------------------------------------------
recursive function logCp(kappa) result(lCp)

use local
use constants
use math 

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: kappa
real(kind=dbl)                  :: lCp

! pre-computed constants
real(kind=dbl),parameter        :: C=-3.675754132818690967D0    ! C = ln(1.D0/(2.D0*cPi)**2)
real(kind=dbl),parameter        :: C2=4.1746562059854348688D0   ! C2 = ln(512/sqrt(2)/pi^(3/2))

! for arguments larger than kappa=30, we use a simple numerical approximation
if (kappa.gt.30.D0) then 
  lCp = kappa**4.5D0/(-105D0+8.D0*kappa*(-15.D0+16.D0*kappa*(-3.D0+8.D0*kappa)))
  lCp = C2 - kappa + dlog(lCp)
else 
  lCp = C + dlog( kappa / BesselI1(kappa) )
end if

end function logCp

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! the following routines have to do with the neighborhood similarity analysis
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
! SUBROUTINE: DI_similarity_classifier
!
!> @author Saransh Singh, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief classify the point as grain interior or anomalous point
!
!> @details takes the kNN neighbor information as input and returns the 
!  whether the point lies in the interior of the grain or lies on the 
!  grain boundary. Details in pg 11 of the dictionary indexing paper
!
!> @param array input array
!> @param k number of top matches for each pixel
!> @param npx number of pixels in the x direction
!> @param npy number of pixels in the y direction
!
!> @date 01/05/15 SS 1.0 original
!> @date 01/06/15 MDG 1.1 simplified summation loop and renamed routine
!--------------------------------------------------------------------------
subroutine DI_Similarity_Classifier(array,k,npx,npy,returnarr)

use local

IMPLICIT NONE

integer(kind=sgl),INTENT(IN)            :: array(npx,npy,k)
integer(kind=irg),INTENT(IN)            :: k
integer(kind=irg),INTENT(IN)            :: npx
integer(kind=irg),INTENT(IN)            :: npy
real(kind=sgl),INTENT(OUT)              :: returnarr(npx,npy)

integer(kind=irg)                       :: ii,jj,ki,kj,similarity_measure_sum,res
real(kind=sgl)                          :: similarity_measure


similarity_measure_sum = 0
similarity_measure = 0.0
returnarr = 0.0

do ii = 2,npx-1
    do jj = 2,npy-1
      do ki = -1, 1
        do kj = -1, 1
          if ((abs(ki)+abs(kj)).ne.0) then 
            call CardIntersection(array(ii+ki,jj+kj,1:k),array(ii,jj,1:k),k,res)
            similarity_measure_sum = similarity_measure_sum + res

            similarity_measure = float(similarity_measure_sum)/float(8*k)
            returnarr(ii,jj) = similarity_measure

            similarity_measure = 0.0
            similarity_measure_sum = 0
          end if
        end do
      end do
    end do
end do

end subroutine DI_Similarity_Classifier

!--------------------------------------------------------------------------
!
! SUBROUTINE: CardIntersection
!
!> @author Saransh Singh, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief calculate the cardinality of the intersection of two sets
!
!> @param set1
!> @param set2
!> @param k number of elements in each set
!
!> @date 01/05/15 SS 1.0 original
!> @date MDG 1.1 changed types to integer 
!--------------------------------------------------------------------------
recursive subroutine CardIntersection(set1,set2,k,res)

use local

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: set1(k)
integer(kind=irg),INTENT(IN)            :: set2(k)
integer(kind=irg),INTENT(IN)            :: k
integer(kind=irg),INTENT(OUT)           :: res

integer(kind=irg)                       :: ii,jj
jj = 1
res = 0

do ii = 1,k
    do jj = 1,k
        if (set1(ii) .eq. set2(jj)) then
            res = res + 1
            EXIT
        end if
    end do
end do

end subroutine CardIntersection


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! the following routines have to do with the axial Watson distribution
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
! FUNCTION: WatsonDensity
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief computes the axial Watson density function
!
!> @param X input quaternion samples
!> @param nums number of samples
!> @param mu mean direction
!> @param kappa concentration
!> @param CW logCpW(kappa) (precomputed in calling routine)
!
!> @date 01/24/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function WatsonDensity(X,nums,mu,kappa,CW) result(watson)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: X(4,nums)
integer(kind=irg),INTENT(IN)            :: nums
real(kind=dbl),INTENT(IN)               :: mu(4)
real(kind=dbl),INTENT(IN)               :: kappa
real(kind=dbl),INTENT(IN)               :: CW
real(kind=dbl)                          :: watson(nums)

integer(kind=irg)                       :: j

do j=1,nums
 watson(j) = dexp(CW+kappa*dot_product(mu,X(1:4,j))**2)
end do

end function WatsonDensity

!--------------------------------------------------------------------------
!
! FUNCTION: logCpW
!
!> @author Marc De Graef, Carnegie Mellon University 
!
!> @brief computes the logarithm of Cp for the Watson distribution function
!
!> @details For large value of kappa (>20), we use an approximation based on
!> the large argument expansions for the modified Bessel functions
!
!> @param kappa input parameter
!
!> @date 01/25/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function logCpW(kappa) result(lCpW)

use local
use constants
use math 

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: kappa
real(kind=dbl)                  :: lCpW

! pre-computed constants
real(kind=dbl),parameter        :: C2=5.4243952068443172530D0   ! C2 = ln(128*sqrt(pi))

! for arguments larger than kappa=20, we use a simple numerical approximation
if (kappa.gt.20.D0) then 
  lCpW = kappa**4.5D0/(525.D0 + 4.D0*kappa*(45.D0 + 8.D0*kappa*(3.D0 + 4.D0*kappa)))
  lCpW = C2 - kappa + dlog(lCpW)
else 
  lCpW = -kappa*0.5D0 - dlog( BesselI0(kappa*0.5D0) - BesselI1(kappa*0.5D0) )
end if

end function logCpW


!--------------------------------------------------------------------------
!
! SUBROUTINE: DI_EMforWatson
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Expectation maximization approach to maximum likelihood problem for mu and kappa, axial Watson distribution
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
!> @date 01/25/15 MDG 1.0 original, based on Chen's original VMF implementation
!--------------------------------------------------------------------------
recursive subroutine DI_EMforWatson(X, dict, nums, seed, muhat, kappahat)

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
integer(kind=irg)                       :: FZtype, FZorder
real(kind=dbl),allocatable              :: Mu_All(:,:), Kappa_All(:), R_All(:,:,:), L_All(:), &
                                           R(:,:), Q(:), L(:)
real(kind=dbl)                          :: Mu(4), PmMu(4), MuKa(5), Qi, Li, rod(4), qu(4), Kappa

! In this routine, we perform the EM algorithm to obtain an estimate for the 
! mean direction and concentration parameter of the axial Watson distribution
! that models the statistics of the orientation point cloud for the case where
! -q and q represent identical rotations (i.e., as is the case for quaternions)

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
! quaternion multiplication has been verified against the 4x4 matrix multiplication of the Matlab code on 01/02/15 
  iloop: do i=1,dict%Num_of_iterations 
! E-step
    R = Watson_Estep(X,dict,Pmdims,N,Mu,Kappa)
! M-step
    MuKa = Watson_Mstep(X,dict,Pmdims,N,R)
! calculate the Q and Likelihood function values
    call Watson_getQandL(X,dict,Pmdims,nums,MuKa,R,Qi,Li)
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
      if (abs(Q(i)-Q(i-1)).lt.0.01) then 
        EXIT iloop
      end if
    end if
  end do iloop
  deallocate(R,Q,L)
end do

dd = maxloc(L_All,1)
Mu = Mu_all(dd,1:4)
kappahat = Kappa_All(dd)

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

end subroutine DI_EMforWatson

!--------------------------------------------------------------------------
!
! FUNCTION: Watson_Estep
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief computes the E step of the EM process, 
!
!> @param X list of input quaternions
!> @param dict pointer to dictionary type
!> @param Pmdims number of quaternion symmetry operators
!> @param nums number of samples
!> @param Mu current guess for mean quaternion
!> @param Kappa input parameter
!
!> @date 01/01/15 MDG 1.0 original
!> @date 01/06/15 MDG 1.1 added optional argument full
!> @date 01/09/15 MDG 1.1 removed optional argument full (incorporated in dicttype)
!--------------------------------------------------------------------------
recursive function Watson_Estep(X,dict,Pmdims,nums,Mu,Kappa) result(R)

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

integer(kind=irg)                       :: j, i 
real(kind=dbl)                          :: Rdenom(nums), PmMu(4), coshsum(nums), arg(nums), qu(4), C


C = logCpW(kappa)
! 
do j=1,Pmdims
    PmMu = quat_mult(Mu,dict%Pm(1:4,j))
    R(1:nums,j) = WatsonDensity(X, nums, PmMu, Kappa, C)
end do
! and determine the normalization factors
Rdenom = 1.D0/sum(R,2)

do j=1,Pmdims 
  R(1:nums,j) = R(1:nums,j)*Rdenom(1:nums)
end do

end function Watson_Estep


!--------------------------------------------------------------------------
!
! FUNCTION: Watson_Mstep
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief computes the M step of the EM process for the Watson distribution
!
!> @param X list of input quaternions
!> @param dict pointer to dictionary type
!> @param Pmdims number of quaternion symmetry operators
!> @param nums number of samples
!> @param R weight factors form the E step
!
!> @date 01/01/15 MDG 1.0 original
!> @date 01/06/15 MDG 1.1 added optional full parameter
!> @date 01/09/15 MDG 1.2 removed optional argument full (incorporated in dicttype)
!> @date 01/09/15 MDG 1.3 introduced accurate numerical approximation for kappa determination
!--------------------------------------------------------------------------
recursive function Watson_Mstep(X,dict,Pmdims,nums,R) result(MuKa)

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

real(kind=dbl)                          :: diff(dict%Apnum), qu(4), y, Tscatt(4,4), tmp(4,4), D4k
integer(kind=irg)                       :: minp, i, j

! variables needed for the dsyev Lapack eigenvalue routine
CHARACTER                               :: JOBZ, UPLO
INTEGER                                 :: INFO, LDA, LWORK, NN
DOUBLE PRECISION                        :: A( 4 , 4 ), W( 4 ), WORK( 20 )

! here, we compute the modified scattering matrix Tscatt and compute its largest eigenvalue
! under the assumption that kappa will always be positive for the types of problems that we
! need to consider; we need to use the outer product, implemented using spread calls
Tscatt = 0.D0
do j=1,Pmdims 
  do i=1,nums
      qu = quat_mult(X(1:4,i),conjg(dict%Pm(1:4,j)))
      tmp = spread(qu,dim=2,ncopies=4)*spread(qu,dim=1,ncopies=4)
      Tscatt = Tscatt + R(i,j) * tmp
  end do
end do
Tscatt = Tscatt/dble(nums)

JOBZ = 'V'
UPLO = 'U'
NN = 4
LDA = 4
LWORK = 20

A = Tscatt

call DSYEV( JOBZ, UPLO, NN, A, LDA, W, WORK, LWORK, INFO )

qu(1:4) = A(1:4,4)
MuKa(1:4) = qu(1:4)

! this next step comes from Mardia&Jupp; may need to be modified using Yu Hui's expression...
D4k = dot_product(qu,matmul(Tscatt,qu))

! and here we invert D4k to get an estimate of kappa
if (D4k.ge.0.94D0) then
  MuKa(5) = (5.D0*D4k-11.D0-dsqrt(39.D0-12.D0*D4k+9.D0*D4k**2))/(8.D0*(D4k-1.D0))
else
  diff = dabs( D4k - dict%yAp ) 
  minp = minloc( diff, 1 )
  if (minp.eq.1) minp = 2 
  MuKa(5) = dict%xAp(minp)
end if

end function Watson_Mstep

!--------------------------------------------------------------------------
!
! SUBROUTINE: Watson_getQandL
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
!> @date 01.06/15 MDG 1.1 added optional full parameter
!> @date 01/09/15 MDG 1.1 removed optional argument full (incorporated in dicttype)
!> @date 02/04/15 MDG 2.0 new routine, derived from equivalent VMF routine
!--------------------------------------------------------------------------
recursive subroutine Watson_getQandL(X,dict,Pmdims,nums,MuKa,R,Q,L)

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

real(kind=dbl)                          :: Phi(nums,Pmdims), PmMu(4), qu(4), C
integer(kind=irg)                       :: j
real(kind=dbl),parameter                :: eps = 0.00001D0

  C = logCpW(MuKa(5))
! compute the auxiliary Phi array
  Phi = 0.D0
  qu = MuKa(1:4)
  do j=1,Pmdims
    PmMu = quat_mult(dict%Pm(1:4,j), qu)
    Phi(1:nums,j) = WatsonDensity(X, nums, PmMu, MuKa(5), C)
  end do
  Phi = Phi/dble(dict%Nqsym)

! and convert the array into the Q and L parameters.
  L = sum(dlog(sum(Phi,2)))
  Q = sum(R*dlog(Phi+eps))

end subroutine Watson_getQandL


end module dictmod
