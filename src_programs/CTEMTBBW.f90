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
! CTEMsoft2013:CTEMTBBW.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMTBBW 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Simple two-beam Bloch wave program
! 
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program CTEMTBBW

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

IMPLICIT NONE

integer(kind=irg)          		:: ns,g(3),k(3),fn(3), io_int(1)
real(kind=sgl)             		:: ktmax,io_real(1)
character(20)    			:: oname

 progname = 'CTEMTBBW.f90'
 progdesc = 'Two-beam bright field-dark field images, using Bloch waves'
 call CTEMsoft
                                                

! first get the crystal data and microscope voltage
 SG % SYM_reduce=.TRUE.
 call CrystalData
 call GetVoltage

! generate all atom positions
 call CalcPositions('v')

! get the reciprocal lattice vector
 mess = 'Diffraction vector :'; call Message("(A)")
 call GetIndex(g,'r')

! ask the user for the beam direction 
 mess = 'Beam direction :'; call Message("(A)")
 call GetIndex(k,'d')
 mess = 'Foil normal    :'; call Message("(A)")
 call GetIndex(fn,'d')

! some more parameters
 call ReadValue(' Enter maximum value of k_t in units of g: ', io_real, 1)
 ktmax = io_real(1)
 call ReadValue(' Number of orientations: ', io_int, 1)
 ns=io_int(1)

 call ReadValue(' Output file name : ',oname, "(A)")

! do the computation
 call CalcTBBW(g,float(k),float(fn),ns,ktmax,oname)

end program

!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcTBBW 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief This is where the actual computation is performed
! 
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine CalcTBBW(g,k,f,ns,ktmax,oname)

use local
use io
use crystal
use crystalvars
use diffraction
use constants
use dynamical

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: g(3)
real(kind=sgl),INTENT(INOUT)	:: k(3)
real(kind=sgl),INTENT(IN)	:: f(3)
integer(kind=irg),INTENT(IN)	:: ns
real(kind=sgl),INTENT(IN)	:: ktmax
character(20),INTENT(IN)		:: oname

real(kind=sgl)            		:: Vmod,Vphase,Vpmod,Vpphase,pre,upzero,find(3), &
				        kk,kt(3),kttb,kn,kz,io_real(1),pre2,dkt,gg,s
complex(kind=dbl)         		:: M(2,2),alph(2),CGinv(2,2),Mcp(2,2)
integer(kind=irg)         		:: ind(3),ivec(3),ik,izero, IPIV(2), io_int(2),i,j,nn

! pre converts from V to U
 pre = 2.0*sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18

! scaling factor for excitation error (2*k_0)
 pre2 = 2.0/sngl(cell%mLambda)

! normal aborption potential Uprime_0
 ind = [0,0,0]
 call CalcUcg(ind)
 Vmod = rlp%Vmod
 Vpmod = rlp%Vpmod
 Vphase = rlp%Vphase
 Vpphase = rlp%Vpphase
 upzero = pre*Vpmod

! tranmitted beam
 izero=1

! determine the dynamical matrix M (all but the diagonal)
! i is the row index
 do i=1,2
  ind = (i-1)*int(g) 
! j is the column index
  do j=1,2
   if (j.ne.i) then
    ivec = ind - (j-1)*int(g)
! use Weickenmeier-Kohl scattering parameters and form factors
    call CalcUcg(ivec)
    M(i,j) = -cPi*cmplx(-aimag(rlp%qg),real(rlp%qg))*cmplx(0.0,1.0/cPi/cell%mLambda)
   end if
  end do
 end do

!
! next we iterate over all incident beam directions, and for
! each direction we complete the M-matrix (diagonal).
!
 dkt = 2.0*ktmax/float(ns-1)
 io_real(1) = dkt
 call WriteValue(' beam tilt step size = ', io_real, 1, "(F8.4)")
 find = float(g)
 kk = CalcLength(k,'r')
 gg = CalcLength(find,'r')
 k = k/sngl(cell%mLambda)/kk
 kz = 1.0/cell%mLambda

! open the unformatted output file
 open (unit=15,file=oname,form='unformatted',status ='unknown')
 nn=2
 write (15) 'TB'
 write (15) cell % fname
 write (15) nn
 write (15) ns
 write (15) g
 write (15) k,kz

! loop over the beam directions
 allocate(W(2),CG(2,2))
 do ik = 1,ns
  if (mod(ik,25).eq.0) then
   io_int(1) = ik
   io_int(2) = ns
   call WriteValue(' ', io_int, 2,"('  -> completed column ',I4,' of ',I4)")
  endif

! rescale the wavevector and foil normal
  kt = k + dkt*(float(ik-ns/2)-0.5)*g
  kk = CalcLength(kt,'r')
  kt = kt/sngl(cell%mLambda)/kk

! then complete the diagonal of the M matrix
! i is the row index
  do i=1,nn
   ind = (i-1)*int(g) 
! get the excitation error
   find = float(ind)
   if (i.eq.1) then
    s = 0.0
   else
    s = Calcsg(find,kt,f)
   endif
! and multiply with 2k_0 and store, along with Uprime_0
   M(i,i) = cmplx(pre2*s,upzero)
  end do
!
! next, compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
!

! first, make a copy of M, since BWsolve destroys M
  Mcp = M

! then get the eigenvalues and eigenvectors
  call BWsolve(Mcp,W,CG,CGinv,nn,IPIV)

! the alpha coefficients are in the izero column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
  kttb = dkt*(float(ik-ns/2)-0.5)
  kn = -sqrt(kz**2-(kttb*gg)**2)
  W = W/cmplx(2.0*kn,0.0)
  do i=1,nn
   alph(i) = CGinv(i,izero)
  end do

! store eigenvalues in file along with excitation amplitudes, and eigenvector
! matrix (also the wave vector)
  write (15) kttb,kn
  write (15) W
  write (15) CG
  write (15) alph
 end do

! close the output file
 close(15, status='keep')
 mess = 'All data saved in file '//oname; call Message("(A)")
 
end subroutine CalcTBBW
