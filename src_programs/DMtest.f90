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
! EMsoft:DMtest.f90
!--------------------------------------------------------------------------
!
! PROGRAM: DMtest
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief test the dynamical matrix for multiple space group settings
!
!> @details This program takes the silicon structure in both space group
!> settings and compares the dynamical matrix entries as well as some other
!> EBSD-related quantities in an attempt to identify a problem with phases.
! 
!> @date   09/06/15 MDG 1.0 original
!--------------------------------------------------------------------------
program DMtest

use local
use initializers
use diffraction
use gvectors
use kvectors
use symmetry
use MBmodule
use io

IMPLICIT NONE

type(unitcell),pointer          :: cell1, cell2
type(DynType)                   :: Dyn1, Dyn2
type(kvectorlist),pointer       :: khead, ktmp
type(gnode)                     :: rlp1, rlp2
type(symdata2D)                 :: TDPG
type(BetheParameterType)        :: BetheParameters
type(reflisttype),pointer       :: reflist1, reflist2, firstw1, firstw2, rltmp
complex(kind=dbl),allocatable   :: DynMat1(:,:), DynMat2(:,:)

real(kind=sgl)                  :: dmin, voltage, ga(3), k(3), ktmax, klaue(2), FN(3), thick(1), lambdaE(1)
real(kind=dbl)                  :: depthstep, th
integer(kind=irg)               :: npx, npy, isym, numk, ijmax, nref1, nref2, i, j, nns1, nns2, nnw1, nnw2, numset, nt, izz, gzero
logical                         :: verbose
character(fnlen)                :: progname, progdesc, xtalname1, xtalname2
complex(kind=dbl),allocatable   :: Sgh1(:,:,:), Sgh2(:,:,:)
complex(kind=dbl),allocatable   :: Lgh1(:,:), Lgh2(:,:)
integer(kind=irg),allocatable   :: nat1(:), nat2(:)
complex(kind=dbl)               :: sum1, sum2
real(kind=sgl),allocatable      :: inten1(:,:), inten2(:,:)

progname = 'DMtest.f90'
progdesc = 'test program for dynamical matrices'

! print some information
call EMsoft(progname, progdesc)

! cell1 is the Si structure in the first space group setting, cell2 in the second setting
nullify(cell1, cell2)
nullify(khead)
nullify(ktmp)

allocate(cell1, cell2)
xtalname1 = 'Si1.xtal'
xtalname2 = 'Si2.xtal'
voltage = 20000.0
dmin = 0.06

! load the structure files and initialize associated arrays
verbose = .TRUE.
call Initialize_Cell(cell1,Dyn1,rlp1, xtalname1, dmin, voltage, verbose)
call Initialize_Cell(cell2,Dyn2,rlp2, xtalname2, dmin, voltage, verbose)

! pick an incident beam direction 
isym = 9
k = (/ 0.0, 1.0, 1.0 /)
klaue = (/ 0.0, 0.0 /)
ga = (/ 1.0, 0.0, 0.0 /)
ktmax = 0.0
npx = 0
npy = 0 
ijmax = 0
call CalckvectorsSymmetry(khead,cell1,TDPG,dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue,.FALSE.)
write (*,*) 'number of beams selected : ',numk, khead%k

! determine g-vector lists 
write (*,*) 'Bethe = ', BetheParameters%c1, BetheParameters%c2, BetheParameters%c3

FN = sngl(khead%k)
dmin = 0.09
call Initialize_ReflectionList(cell1, reflist1, BetheParameters, FN, sngl(khead%k), dmin, nref1, verbose)
call Initialize_ReflectionList(cell2, reflist2, BetheParameters, FN, sngl(khead%k), dmin, nref2, verbose)
write (*,*) 'number of g-vectors in list ',nref1, nref2
write (*,*) 'Si setting 1'
rltmp => reflist1%next
do i=1, nref1
  call CalcUcg(cell1, rlp1, rltmp%hkl)
  write (*,*) rltmp%hkl, rltmp%sg, rltmp%Ucg, rlp1%qg
  rltmp => rltmp%next
end do

write (*,*) 'Si setting 2'
rltmp => reflist2%next
do i=1, nref2
  call CalcUcg(cell2, rlp2, rltmp%hkl)
  write (*,*) rltmp%hkl, rltmp%sg, rltmp%Ucg, rlp2%qg
  rltmp => rltmp%next
end do

! ok, so both structures produce the same list of reflections, but the Fourier coefficients are different
! next, we compute the Bloch wave dynamical matrices and compare them
! determine strong and weak reflections
call Apply_BethePotentials(cell1, reflist1, firstw1, BetheParameters, nref1, nns1, nnw1)
write (*,*) 'strong/weak 1: ',nns1,nnw1

call Apply_BethePotentials(cell2, reflist2, firstw2, BetheParameters, nref2, nns2, nnw2)
write (*,*) 'strong/weak 2: ',nns2,nnw2
allocate(DynMat1(nns1,nns1))

! generate the dynamical matrix
call GetDynMat(cell1, reflist1, firstw1, rlp1, DynMat1, nns1, nnw1)

write (*,*) 'Dynamical matrix 1'
do i=1,4
  write (*,*) (DynMat1(i,j),j=1,4)
end do

allocate(DynMat2(nns2,nns2))
call GetDynMat(cell2, reflist2, firstw2, rlp2, DynMat2, nns2, nnw2)

write (*,*) 'Dynamical matrix 2'
do i=1,4
  write (*,*) (DynMat2(i,j),j=1,4)
end do


! get the eigenvalues etc for the Bloch wave approach and compute intensities for some thickness
nt = 1
thick(1) = 50.0
allocate(inten1(nt,nns1))
allocate(inten2(nt,nns2))

call CalcPEDint(DynMat1,cell1,sngl(khead%kn),nns1,nt,thick,inten1)
call CalcPEDint(DynMat2,cell2,sngl(khead%kn),nns2,nt,thick,inten2)

rltmp => reflist1%next
do i=1,nns1
  write (*,*) rltmp%hkl, inten1(1,i), inten2(1,i), abs(inten1(1,i)-inten2(1,i))
  rltmp=>rltmp%nexts
end do


! compute the Sgh matrix for EBSD
numset = 1
allocate(Sgh1(nns1,nns1,numset), Sgh2(nns2,nns2,numset),nat1(numset), nat2(numset))

call CalcSgh(cell1,reflist1,nns1,numset,Sgh1,nat1)
call CalcSgh(cell2,reflist2,nns2,numset,Sgh2,nat2)

write (*,*) 'Sgh matrix 1'
do i=1,6
  write (*,*) (Sgh1(i,j,1),j=1,4)
end do

write (*,*) 'Sgh matrix 2'
do i=1,6
  write (*,*) (Sgh2(i,j,1),j=1,4)
end do

do i=1,nns1
  do j=1,nns1
    if (abs(Sgh1(i,j,1)-Sgh2(i,j,1)).gt.0.1) write (*,*) i,j,Sgh1(i,j,1),' <-> ',Sgh2(i,j,1)
  end do
end do

! compute the Lgh matrices
allocate(Lgh1(nns1,nns1), Lgh2(nns2,nns2))
gzero = 1
th = 50.D0
depthstep = 50.D0
lambdaE(1) = 1.0
izz = 1

call CalcLgh(DynMat1,Lgh1,th,khead%kn,nns1,gzero,depthstep,lambdaE,izz)
call CalcLgh(DynMat2,Lgh2,th,khead%kn,nns2,gzero,depthstep,lambdaE,izz)

write (*,*) 'Lgh matrix 1'
do i=1,6
  write (*,*) (Lgh1(i,j),j=1,4)
end do

write (*,*) 'Lgh matrix 2'
do i=1,6
  write (*,*) (Lgh2(i,j),j=1,4)
end do

sum1 = cmplx(0.D0,0.D0)
sum2 = cmplx(0.D0,0.D0)

do i=1,nns1
  do j=1,nns1
    if (abs(Lgh1(i,j)-Lgh2(i,j)).gt.0.0001) write (*,*) i,j,Lgh1(i,j),' <-> ',Lgh2(i,j)
    sum1 = sum1 + Sgh1(i,j,1)*Lgh1(i,j)
    sum2 = sum2 + Sgh2(i,j,1)*Lgh2(i,j)
  end do
end do

write (*,*) 'Total sum = ', sum1, sum2

end program DMtest

