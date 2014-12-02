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
! CTEMsoft2013:MBmodule.f90
!--------------------------------------------------------------------------
!
! MODULE: MBmodule
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routine for multi-beam diffraction calculations, used by more than one program
!
!> @date 7/10/13   MDG 1.0 original
!> @date 7/12/13   MDG 1.1 added forward cube to ball to quaternion mappings
!> @date 8/01/13   MDG 1.2 added standard Lambert projection
!> @date 8/12/13   MDG 1.3 added inverse Lambert projections for Ball to Cube
!> @date 9/20/13   MDG 1.4 added ApplyLaueSymmetry
!--------------------------------------------------------------------------
module MBmodule

use local
use typedefs

IMPLICIT NONE

contains


!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcBWint
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the scattered intensities for a range of thicknesses
!
!> @param Dyn dynamical scattering structure
!> @param cell unit cell pointer
!> @param ktmp wave vector structure
!> @param BetheParameter Bethe potential parameters
!> @param nn number of strong beams
!> @param nw number of weak beams
!> @param nt number of thickness values
!> @param thick thickness array
!> @param inten output intensity list, including weak beams
!
!> @date  10/13/98 MDG 1.0 original
!> @date   7/04/01 MDG 2.0 f90
!> @date  04/29/13 MDG 3.0 inclusion of Bethe weak beams
!> @date  06/10/14 MDG 4.0 added Dyn, cell, ktmp, and BetheParameter arguments
!--------------------------------------------------------------------------
subroutine CalcBWint(Dyn,cell,ktmp,BetheParameter,nn,nw,nt,thick,inten)

use io
use diffraction
use kvectors
use gvectors
use constants

IMPLICIT NONE

type(DynType),INTENT(INOUT)     :: Dyn
type(unitcell),pointer          :: cell
type(kvectorlist),pointer       :: ktmp
type(BetheParameterType),INTENT(IN) :: BetheParameter
integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
integer(kind=irg),INTENT(IN)    :: nw                   !< number of weak beams
integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
real(kind=sgl),INTENT(IN)       :: thick(nt)            !< thickness array
real(kind=sgl),INTENT(INOUT)    :: inten(nt,nn+nw)      !< output intensities (both strong and weak)

integer(kind=irg)               :: i,j,IPIV(nn), ll(3), jp
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), &
                                   lalpha(nn), delta(nn,nn), weak(nw,nn), Ucross(nw,nn), tmp(nw,nn), c
real(kind=sgl)                  :: th

! compute the eigenvalues and eigenvectors
 Minp = Dyn%DynMat
 IPIV = 0
 call BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
 lW = cPi*Wloc/cmplx(ktmp%kn,0.0)
 do i=1,nn
  lalpha(i) = CGinv(i,1)
 end do

! in preparation for the intensity computation, we need the prefactor array for the
! weak beam amplitude computation...
! we will also need a potential coefficient array for the cross coefficients, which we 
! compute in the same loop
 do j=1,nn   ! strong beam loop
   do jp=1,nw  ! weak beam loop
! prefactor value
     c = cmplx(2.D0*BetheParameter%weaksg(jp)/cell%mLambda) - 2.D0*ktmp%kn*Wloc(j)
     weak(jp,j) = cmplx(-1.D0,0.D0)/c
! cross potential coefficient
     ll(1:3) = BetheParameter%weakhkl(1:3,jp) - BetheParameter%stronghkl(1:3,j)
     Ucross(jp,j) = cell%LUT( ll(1),ll(2),ll(3) )
   end do
 end do

! compute the strong beam intensities, stored in the first nn slots of inten 
! we can also compute the weak beams, since they make use of the same diag(1:nn) expression
! as the strong beams, plus a few other factors (excitation error, wave length, Fourier coefficients)
 do i=1,nt
  th = thick(i)
  diag(1:nn)=exp(-th*imag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
! the delta array is common to the strong and weak beam intensity computation, so we compute it first
  do j=1,nn
   delta(j,1:nn) = lCG(j,1:nn)*diag(1:nn)
  end do
! strong beams
  do j=1,nn
   inten(i,j) = cdabs(sum(delta(j,1:nn)))**2
  end do 
! weak beams
  tmp = matmul(Ucross,delta)
  do jp=1,nw
   inten(i,nn+jp) = cdabs( sum(weak(jp,1:nn)*tmp(jp,1:nn)) )**2
  end do  
 end do
   
end subroutine CalcBWint


!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcPEDint
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the scattered intensities for the Precession Electron Diffraction mode
!
!> @param DynMat dynamical matrix
!> @param cell unit cell pointer
!> @param kn  normal component
!> @param nn number of strong beams
!> @param nt number of thickness values
!> @param thick thickness array
!> @param inten output intensity list
!
!> @date  10/13/98 MDG 1.0 original
!> @date   7/04/01 MDG 2.0 f90
!> @date  04/29/13 MDG 3.0 inclusion of Bethe weak beams
!> @date  06/10/14 MDG 4.0 added Dyn, cell, ktmp, and BetheParameter arguments
!> @date  11/28/14 MDG 4.1 forked from CalcBWint
!--------------------------------------------------------------------------
subroutine CalcPEDint(DynMat,cell,kn,nn,nt,thick,inten)

use io
use diffraction
use kvectors
use gvectors
use constants

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)    :: DynMat(nn,nn)
type(unitcell),pointer          :: cell
real(kind=sgl),INTENT(IN)       :: kn
integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
real(kind=sgl),INTENT(IN)       :: thick(nt)            !< thickness array
real(kind=sgl),INTENT(INOUT)    :: inten(nt,nn)         !< output intensities (both strong and weak)

integer(kind=irg)               :: i,j,IPIV(nn), ll(3), jp
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), &
                                   lalpha(nn), delta(nn,nn), s, c
real(kind=sgl)                  :: th

! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
 lW = cPi*Wloc/cmplx(kn,0.0)
 lalpha(1:nn) = CGinv(1:nn,1)

! make sure the alpha excitation coefficients are normalized 
! s = sum(cdabs(lalpha(1:nn))**2)
! if (s.ne.1.D0) then
!  s = dcmplx(1.D0/dsqrt(s),0.D0)
!  lalpha = lalpha*s
! endif 

! compute the strong beam intensities, stored in the first nn slots of inten 
 do i=1,nt
  th = thick(i)
  diag(1:nn)=exp(-th*imag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
! the delta array is common to the strong and weak beam intensity computation, so we compute it first
  do j=1,nn
   delta(j,1:nn) = lCG(j,1:nn)*diag(1:nn)
  end do
! strong beams
  do j=1,nn
   inten(i,j) = cdabs(sum(delta(j,1:nn)))**2
  end do 
 end do
   
end subroutine CalcPEDint



!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcKint
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the Kossel intensities for a range of thicknesses
!
!> @param DynMat dynamical matrix
!> @param ktmp wave vector structure
!> @param nn number of strong beams
!> @param nt number of thickness values
!> @param thick thickness array
!> @param inten output intensity list, including weak beams
!
!> @date  10/13/98 MDG 1.0 original
!> @date   7/04/01 MDG 2.0 f90
!> @date  04/29/13 MDG 3.0 inclusion of Bethe weak beams
!> @date  01/08/14 MDG 3.1 forked from CalcBWint, specialized for Kossel computation
!> @date  06/10/14 MDG 4.0 added Dyn, cell, and ktmp arguments
!> @date  06/15/14 MDG 4.1 removed global W, CG and alpha initializations
!> @date  06/16/14 MDG 4.2 made routine recursive for OPenMP
!--------------------------------------------------------------------------
recursive subroutine CalcKint(DynMat,kn,nn,nt,thick,Iz)

use local
use io
use diffraction
use kvectors
use gvectors
use constants

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)    :: DynMat(nn,nn)
real(kind=sgl),INTENT(IN)       :: kn
integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
real(kind=sgl),INTENT(IN)       :: thick(nt)            !< thickness array
real(kind=sgl),INTENT(INOUT)    :: Iz(nt)               !< output intensities

integer(kind=irg)               :: j, IPIV(nn), k
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn), Wloc(nn), lCG(nn,nn), lW(nn), lalpha(nn)
real(kind=dbl)                  :: s, q, t


! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)


! the alpha coefficients are in the first column of the inverse matrix
 lW = cPi*Wloc/cmplx(kn,0.0)
 lalpha(1:nn) = CGinv(1:nn,1)

! make sure the alpha excitation coefficients are normalized 
 s = sum(cdabs(lalpha(1:nn))**2)
 if (s.ne.1.D0) then
  s = dcmplx(1.D0/dsqrt(s),0.D0)
  lalpha = lalpha*s
 endif 
 
! compute the thickness array 
 Iz = 0.D0
 do j=1,nn
    q = -4.D0*cPi*aimag(lW(j))
    s = cdabs(lalpha(j))**2
    do k=1,nt
      t = q*thick(k)
      if (abs(t).lt.30.D0) Iz(k) = Iz(k) +  s * exp(t)
    end do
 end do   
   
end subroutine CalcKint

!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcSgh
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute structure factor-like array for EBSD, ECCI and ECP simulations
!
!> @param cell unit cell pointer
!> @param nn dimension of array
!> @param Sgh output array
!> @param nat normalization array
!
!> @date 03/05/14  MDG 1.0 original (used to be in-line in ECP and ECCI programs)
!> @date 03/11/14  MDG 1.1 converted to diagonal Sgh array only
!> @date 06/19/14  MDG 2.0 no globals, taken out of CTEMECCI.f90
!--------------------------------------------------------------------------
recursive subroutine CalcSgh(cell,reflist,nn,numset,Sgh,nat)

use local
use typedefs
use crystal
use gvectors
use constants
use symmetry

IMPLICIT NONE

type(unitcell),pointer                  :: cell
type(reflisttype),pointer               :: reflist
integer(kind=irg),INTENT(IN)            :: nn
integer(kind=irg),INTENT(IN)            :: numset
complex(kind=dbl),INTENT(INOUT)         :: Sgh(nn,nn,numset)
integer(kind=irg),INTENT(INOUT)         :: nat(numset)

integer(kind=irg)                       :: ip, ir, ic, kkk(3), ikk, n
real(kind=sgl)                          :: Znsq, DBWF, kkl
complex(kind=dbl)                       :: carg
real(kind=dbl)                          :: ctmp(192,3),arg, tpi
type(reflisttype),pointer               :: rltmpa, rltmpb

  tpi = 2.D0 * cPi

! for each special position we need to compute its contribution to the Sgh array
  do ip=1,cell % ATOM_ntype
    call CalcOrbit(cell,ip,n,ctmp)
    nat(ip) = cell%numat(ip)
! get Zn-squared for this special position, and include the site occupation parameter as well
    Znsq = float(cell%ATOM_type(ip))**2 * cell%ATOM_pos(ip,4)
! loop over all contributing reflections
! ir is the row index
    rltmpa => reflist%next    ! point to the front of the list
    do ir=1,nn
! ic is the column index
      rltmpb => reflist%next    ! point to the front of the list
      do ic=1,nn
        kkk = rltmpb%hkl - rltmpa%hkl
! We'll assume isotropic Debye-Waller factors for now ...
! That means we need the square of the length of s=  kk^2/4
        kkl = 0.25 * CalcLength(cell,float(kkk),'r')**2
! Debye-Waller exponential times Z^2
        DBWF = Znsq * exp(-cell%ATOM_pos(ip,5)*kkl)
! here is where we insert the proper weight factor, Z^2 exp[-M_{h-g}]
! and also the detector geometry...   For now, we do nothing with the detector
! geometry; the Rossouw et al 1994 paper lists a factor A that does not depend
! on anything in particular, so we assume it is 1. 
        do ikk=1,n
! get the argument of the complex exponential
!         arg = tpi*sum(kkk(1:3)*ctmp(ikk,1:3))
          arg = tpi*sum(kkk(1:3)*cell%apos(ip,ikk,1:3))
          carg = dcmplx(dcos(arg),dsin(arg))
! multiply with the prefactor and add
          Sgh(ir,ic,ip) = Sgh(ir,ic,ip) + carg * dcmplx(DBWF,0.D0)
        end do
        rltmpb => rltmpb%nexts  ! move to next column-entry
      end do
     rltmpa => rltmpa%nexts  ! move to next row-entry
   end do  
  end do
  
end subroutine CalcSgh

! ###################################################################
! 
! SUBROUTINE: CalcLgh
! 
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the Lgh matrix for EBSD, ECCI, ECP, etc simulations
!
!> @param DMat dynamical matrix
!> @param Lgh output array
!> @param thick integration thickness
!> @param kn normal wave vector component
!> @param nn number of strong beams
!> @param gzero index of incident beam
!> @param depthstep depth step size
!> @param lambdaE energy weight factors
!> @param izz number of energy weight factors
!
!> @date 10/13/98  MDG 1.0 original
!> @date 07/04/01  MDG 2.0 f90
!> @date 06/19/14  MDG 3.0 no globals
!> @date 06/23/14  MDG 4.0 moved to MBmodule
! ###################################################################
recursive subroutine CalcLgh(DMat,Lgh,thick,kn,nn,gzero,depthstep,lambdaE,izz)

use local
use io
use files
use constants
use error

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)        :: DMat(nn,nn)
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn)
real(kind=dbl),INTENT(IN)           :: thick
real(kind=dbl),INTENT(IN)           :: kn
integer(kind=irg),INTENT(IN)        :: nn
integer(kind=irg),INTENT(IN)        :: gzero
real(kind=dbl),INTENT(IN)           :: depthstep
real(kind=sgl),INTENT(IN)           :: lambdaE(izz)
integer(kind=irg),INTENT(IN)        :: izz

integer                             :: i,j, iz
complex(kind=dbl)                   :: CGinv(nn,nn), Minp(nn,nn), tmp3(nn,nn)

real(kind=dbl)                      :: tpi, dzt
complex(kind=dbl)                   :: Ijk(nn,nn), q, getMIWORK, qold

integer(kind=irg)                   :: INFO, LDA, LDVR, LDVL,  JPIV(nn), MILWORK
complex(kind=dbl)                   :: CGG(nn,nn), W(nn)
complex(kind=dbl),allocatable       :: MIWORK(:)

integer(kind=irg),parameter         :: LWMAX = 5000 
complex(kind=dbl)                   :: VL(nn,nn),  WORK(LWMAX)
real(kind=dbl)                      :: RWORK(2*nn)
character                           :: JOBVL, JOBVR
integer(kind=sgl)                   :: LWORK

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
! 
 Minp = DMat

! set some initial LAPACK variables 
 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0
 
! first initialize the parameters for the LAPACK ZGEEV, CGETRF, and CGETRI routines
 JOBVL = 'N'   ! do not compute the left eigenvectors
 JOBVR = 'V'   ! do compute the right eigenvectors
 LWORK = -1 ! so that we can ask the routine for the actually needed value

! call the routine to determine the optimal workspace size
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

! then call the eigenvalue solver
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  if (INFO.ne.0) call FatalError('Error in CalcLgh3: ','ZGEEV return not zero')

 CGinv = CGG

 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 MILWORK = -1
 call zgetri(nn,CGinv,LDA,JPIV,getMIWORK,MILWORK,INFO)
 MILWORK =  INT(real(getMIWORK))
 if (.not.allocated(MIWORK)) allocate(MIWORK(MILWORK))
 MIWORK = dcmplx(0.D0,0.D0)
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
 deallocate(MIWORK)

! in all the time that we've used these routines, we haven't
! had a single problem with the matrix inversion, so we don't
! really need to do this test:
!
! if ((cabs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) write (*,*) 'Error in matrix inversion; continuing'


! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)

! recall that alpha(1:nn) = CGinv(1:nn,gzero)

! first the Ijk matrix (this is Winkelmann's B^{ij}(t) matrix)
! combined with numerical integration over [0, z0] interval,
! taking into account depth profiles from Monte Carlo simulations ...
! the depth profile lambdaE must be added to the absorption 
! components of the Bloch wave eigenvalues.

tpi = 2.D0*cPi*depthstep
dzt = depthstep/thick
 do i=1,nn
  do j=1,nn
     q =  cmplx(0.D0,0.D0)
     qold = tpi * dcmplx(aimag(W(i))+aimag(W(j)),real(W(i))-real(W(j)))
     do iz = 1,izz
       q = q + dble(lambdaE(iz)) * cdexp( - qold * dble(iz) ) !MNS changed cexp to cdexp to be compatible with gfortran
     end do
     Ijk(i,j) = conjg(CGinv(i,gzero)) * q * CGinv(j,gzero)
  end do
 end do

Ijk = Ijk * dzt

! then the summations for Lgh and kin
tmp3 = matmul(CGG,transpose(Ijk)) 
Lgh = matmul(tmp3,transpose(conjg(CGG)))

end subroutine CalcLgh


!--------------------------------------------------------------------------
!
! SUBROUTINE: GetDynMat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the dynamical matrix, including Bethe potentials
!
!> @param cell unit cell pointer
!> @param listroot top of the main reflection list
!> @param listrootw top of the weak reflection list
!> @param Dyn dynamical scattering structure
!> @param nns number of strong reflections
!> @param nnw number of weak reflections
!
!> @date  04/22/14 MDG 1.0 new library version
!> @date  06/15/14 MDG 2.0 updated for removal of globals
!> @date  06/17/14 MDG 2.1 added listroot pointers etc to accommodate multiple threads
!> @date  06/18/14 MDG 2.2 corrected some pointer allocation errors in other routines; this one now works fine.
!--------------------------------------------------------------------------
recursive subroutine GetDynMat(cell, listroot, listrootw, rlp, DynMat, nns, nnw)

use local
use typedefs
use io
use crystal
use diffraction
use kvectors
use gvectors
use constants

IMPLICIT NONE

type(unitcell),pointer           :: cell
type(reflisttype),pointer        :: listroot
type(reflisttype),pointer        :: listrootw
type(gnode),INTENT(INOUT)        :: rlp
complex(kind=dbl),INTENT(INOUT)  :: DynMat(nns,nns)
integer(kind=irg),INTENT(IN)     :: nns
integer(kind=irg),INTENT(IN)     :: nnw

complex(kind=dbl)                :: czero, ughp, uhph, weaksum 
real(kind=dbl)                   :: weaksgsum
real(kind=sgl)                   :: Upz
integer(kind=sgl)                :: ir, ic, ll(3), istat, wc
type(reflisttype),pointer        :: rlr, rlc, rlw

czero = cmplx(0.0,0.0,dbl)      ! complex zero

nullify(rlr)
nullify(rlc)
nullify(rlw)

        DynMat = czero
        call CalcUcg(cell, rlp, (/0,0,0/) )
        Upz = rlp%Vpmod

        rlr => listroot%next
        ir = 1
        do
          if (.not.associated(rlr)) EXIT
          rlc => listroot%next
          ic = 1
          do
          if (.not.associated(rlc)) EXIT
          if (ic.ne.ir) then  ! not a diagonal entry
! here we need to do the Bethe corrections if necessary
            if (nnw.ne.0) then
              weaksum = czero
              rlw => listrootw
              do
               if (.not.associated(rlw)) EXIT
               ll = rlr%hkl - rlw%hkl
               ughp = cell%LUT(ll(1),ll(2),ll(3)) 
               ll = rlw%hkl - rlc%hkl
               uhph = cell%LUT(ll(1),ll(2),ll(3)) 
               weaksum = weaksum +  ughp * uhph *cmplx(1.D0/rlw%sg,0.0,dbl)
               rlw => rlw%nextw
              end do
!        ! and correct the dynamical matrix element to become a Bethe potential coefficient
              ll = rlr%hkl - rlc%hkl
              DynMat(ir,ic) = cell%LUT(ll(1),ll(2),ll(3))  - cmplx(0.5D0*cell%mLambda,0.0D0,dbl)*weaksum
             else
              ll = rlr%hkl - rlc%hkl
              DynMat(ir,ic) = cell%LUT(ll(1),ll(2),ll(3))
            end if
          else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! determine the total contribution of the weak beams
            if (nnw.ne.0) then
              weaksgsum = 0.D0
              rlw => listrootw
              do
               if (.not.associated(rlw)) EXIT
                ll = rlr%hkl - rlw%hkl
                ughp = cell%LUT(ll(1),ll(2),ll(3)) 
                weaksgsum = weaksgsum +  cdabs(ughp)**2/rlw%sg
                rlw => rlw%nextw
              end do
              weaksgsum = weaksgsum * cell%mLambda/2.D0
              DynMat(ir,ir) = cmplx(2.D0*rlr%sg/cell%mLambda-weaksgsum,Upz,dbl)
            else
              DynMat(ir,ir) = cmplx(2.D0*rlr%sg/cell%mLambda,Upz,dbl)
            end if           
        
           end if       
           rlc => rlc%nexts
           ic = ic + 1
          end do        
          rlr => rlr%nexts
          ir = ir+1
        end do

end subroutine GetDynMat


!--------------------------------------------------------------------------
!
! SUBROUTINE:CalcLghECP
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the Bloch wave Lgh matrix for ECP mode
!
!> @param nn number of strong beams
!> @param nt number of thickness values
!> @param thick array of thickness values
!> @param kn normal component of incident wave vector
!> @param gzero index of zero beam (should always be the first one; legacy parameter)
!> @param Lgh output array
!
!> @date 11/18/13  MDG 1.0 major rewrite from older ECP program; merged with ECPz
!--------------------------------------------------------------------------
recursive subroutine CalcLghECP(DMat,Lgh,nn,nt,thick,kn,gzero)

use local
use io
use files
use diffraction
use constants

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)        :: DMat(nn,nn)
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn,nt)
integer(kind=sgl),INTENT(IN)        :: nn
integer(kind=sgl),INTENT(IN)        :: nt
real(kind=sgl),INTENT(IN)           :: thick(nt)
real(kind=dbl),INTENT(IN)           :: kn
integer(kind=sgl),INTENT(IN)        :: gzero

integer                             :: i,j,it,ig,ih,IPIV(nn)
complex(kind=dbl),allocatable       :: CGinv(:,:), Minp(:,:),tmp3(:,:)
complex(kind=dbl)                   :: Ijk(nn,nn),q
complex(kind=dbl)                   :: CG(nn,nn), W(nn)


allocate(CGinv(nn,nn),Minp(nn,nn),tmp3(nn,nn))

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
! 
! then get the eigenvalues and eigenvectors
 Minp = DMat
 IPIV = 0

 call BWsolve(Minp,W,CG,CGinv,nn,IPIV)

! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)

! first do the Lgh matrices, looping over the thickness
do it=1,nt
! recall that alpha(1:nn) = CGinv(1:nn,gzero)
! first the Ijk matrix
 do i=1,nn
  do j=1,nn
    q = 2.0*cPi*thick(it)*cmplx(aimag(W(i))+aimag(W(j)),real(W(i))-real(W(j)))
    Ijk(i,j) = conjg(CGinv(i,gzero)) * (1.0-exp(-q))/q * CGinv(j,gzero)
  end do
 end do

! then the summations for Lgh
 do ih=1,nn
   do i=1,nn
      tmp3(ih,i) = sum(Ijk(i,1:nn)*CG(ih,1:nn))
   end do
 end do
 do ig=1,nn
  do ih=1,nn
     Lgh(ih,ig,it) = sum(conjg(CG(ig,1:nn))*tmp3(ih,1:nn))
  end do
 end do
end do ! thickness loop

deallocate(CGinv,Minp,tmp3)

end subroutine CalcLghECP




end module MBmodule
