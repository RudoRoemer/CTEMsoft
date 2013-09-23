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
! CTEMsoft2013:CTEMTBSM.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMTBSM
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Two-beam scattering matrix program
!
!> @detail This is a very simplistic implementation to illustrate how the scattering matrix
!> approach can be implemented; this approach is also used in more sophisticated programs.
!
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program CTEMTBSM

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

IMPLICIT NONE

integer(kind=irg)     	:: g(3)
real(kind=sgl)        		:: wmax,io_real(1)

 progname = 'CTEMTBSM.f90'
 progdesc = 'Two-beam computations: comparing scattering matrix and analytical'
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

! some more parameters
 call ReadValue(' Enter maximum value of w: ', io_real, 1)
 wmax = io_real(1)
 call ReadValue(' Maximum foil thickness : ', io_real, 1)

! do the computation
 call CalcTBSM(g,wmax,io_real(1))

end program


!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcTBSM
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Actual two-beam scattering matrix computation
!
!> @param g reciprocal lattice vector
!> @param wmax maximum w parameter
!> @param tmax maximum foil thickness
!
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine CalcTBSM(g,wmax,tmax)

use local
use io
use crystal
use crystalvars
use diffraction
use constants
use postscript
use dynamical

integer(kind=irg),INTENT(IN)	:: g(3)
real(kind=sgl),INTENT(IN)	:: wmax
real(kind=sgl),INTENT(IN)	:: tmax

integer(kind=irg),parameter	:: ns=512, nt=512
integer(kind=irg)      		:: ind(3),io_int(1)
real(kind=sgl)         			:: qr,qi,bg,Vmod,Vphase,Vpmod,Vpphase,xg,xgp,xgpz, &
                          			   Ar(2,2), Ai(2,2), sg, dz, dsg, p(2),q(2), &
                          			   BF(512,512,2),DF(512,512,2),It,Is,io_real(1)

! normal aborption factor
 ind = [0,0,0]
 call CalcUcg(ind)
 xgpz= rlp%xgp

! extinction distance and anomalous absorption length
 call CalcUcg(g)
 xgp = rlp%xgp
 xg  = rlp%xg
 imanum=0
 dz = tmax/float(nt-1)

! compute the scattering matrix SM for each orientation
 dsg = 2.0*wmax/float(ns-1)
 mess = 'Starting scattering matrix computation'; call Message("(A)")
 io_int(1) = 8*ns
 call WriteValue(' Number of computations per * = ', io_int, 1, "(I4)")

! allocate arrays
 allocate(SMr(512,2,2))
 allocate(SMi(512,2,2))
 allocate(phir(512,2))
 allocate(phii(512,2))
 
T0 = SECNDS(0.0)
 do i=1,ns
  sg = (-wmax+dsg*float(i))/xg
  call TBCalcSM(Ar,Ai,sg,dz,xg,xgp,xgpz,bg)
  SMr(i,1:2,1:2)=Ar(1:2,1:2)
  SMi(i,1:2,1:2)=Ai(1:2,1:2)
 end do

! compute intensities for the first row
  do i=1,ns
   BF(i,1,1)=SMr(i,1,1)**2+SMi(i,1,1)**2
   DF(i,1,1)=SMr(i,2,1)**2+SMi(i,2,1)**2

! the first column of SM is the actual initial wavefunction at z=dz
   phir(i,1)=SMr(i,1,1)
   phir(i,2)=SMr(i,2,1)
   phii(i,1)=SMi(i,1,1)
   phii(i,2)=SMi(i,2,1)
  end do

! loop to maximum thickness
  do l=2,nt
   do k=1,ns
    do i=1,2
     p(i)=0.D0
     q(i)=0.D0
     do j=1,2
      p(i)=p(i)+SMr(k,i,j)*phir(k,j)-SMi(k,i,j)*phii(k,j)
      q(i)=q(i)+SMi(k,i,j)*phir(k,j)+SMr(k,i,j)*phii(k,j)
     end do
    end do
    phir(k,1:2)=p(1:2)
    phii(k,1:2)=q(1:2)
   end do
   if (mod(l,8).eq.0) then
     mess = '*'
     call Message("(1A,$)")
   end if
   do j=1,ns
    BF(j,l,1) = phir(j,1)**2+phii(j,1)**2
    DF(j,l,1) = phir(j,2)**2+phii(j,2)**2
   end do
  end do
  T1 = SECNDS(T0)
  mess = 'done'; call Message("(A)")
  io_real(1) = T1
  call WriteValue('  Total computation time [s] ', io_real, 1, "(F)")
  mess = 'Starting direct analytical solution'; call Message("(/A)")

! next, redo the computation, but this time use TBCalcInten
  T0 = SECNDS(0.0)
  do i=1,ns
   sg = (-wmax+dsg*float(i))/xg
   do j=1,nt
    t = float(j)*dz
    call TBCalcInten(It,Is,sg,t,xg,xgp,xgpz,bg)
    BF(i,j,2) = It
    DF(i,j,2) = Is
   end do
   if (mod(i,8).eq.0) then
     mess = '*'
     call Message("(1A,$)")
   end if
  end do
  T2 = SECNDS(T0)
  mess = 'done'; call Message("(A)")
  io_real(1) = T2
  call WriteValue('  Total computation time [s] ', io_real, 1, "(F)")

! compare the two computations
  bfdiff = 0.0
  dfdiff = 0.0
  do i=1,ns
   do j=1,nt
    bfdiff = bfdiff + (BF(i,j,1)-BF(i,j,2))**2
    dfdiff = dfdiff + (DF(i,j,1)-DF(i,j,2))**2
   end do
  end do
  bfdiff = bfdiff/float(ns)/float(nt)
  dfdiff = dfdiff/float(ns)/float(nt)
  io_real(1) = bfdiff
  call WriteValue(' Average difference in BF ', io_real, 1, "(F)")
  io_real(1) = dfdiff
  call WriteValue(' Average difference in DF ', io_real, 1, "(F)")
  mess = 'Images computed, preparing for PS output'; call Message("(A)")

! open PostScript file
  PS % pspage = 0
  call PS_openfile
  call PS_newpage(.TRUE.,'Two Beam Rocking Curves')
  call PS_setfont(PSfonts(2),0.16)
  call PS_text(2.5,8.6,'Scattering Matrix Solution')
  call PS_text(2.5,2.1,'Direct Analytical Solution')
  call PS_setfont(PSfonts(2),0.10)
  call PS_text(0.5,1.8,'Input file')
  call PS_text(2.6,1.8,cell % fname)
  call PS_text(0.5,1.6,'Active reflection')
  call PS_textint(2.5,1.6,' ',g(1))
  call PS_textint(2.6,1.6,' ',g(2))
  call PS_textint(2.7,1.6,' ',g(3))
  call PS_text(0.5,1.4,'Extinction distance  [nm]')
  call PS_textvar(2.5,1.4,' ',xg)
  call PS_text(0.5,1.2,'Anomalous absorption length [nm]')
  call PS_textvar(2.5,1.2,' ',xgp)
  call PS_text(0.5,1.0,'Normal absorption length [nm]')
  call PS_textvar(2.5,1.0,' ',xgpz)
  call PS_text(0.5,0.8,'Accelerating Voltage [V]')
  call PS_textvar(2.5,0.8,' ',sngl(mAccvol))
  call PS_text(0.5,0.6,'Maximum w')
  call PS_textvar(2.5,0.6,' ',wmax)
  call PS_text(0.5,0.4,'Maximum thickness [nm]')
  call PS_textvar(2.5,0.4,' ',tmax)

 imanum = 0	! must be defined, otherwise program might end in segmentation fault
! Bright Field image
  allocate(imaint(512,512))
  imaint(1:ns,1:nt)=int(255.0*BF(1:ns,1:nt,1))
  x0=0.2
  y0=5.5
  npx=ns
  npy=nt
  scl=3.0
  call PS_DumpImageDistort(x0,y0,npx,npy,scl,scl)

! Dark Field image
  imaint(1:ns,1:nt)=int(255.0*DF(1:ns,1:nt,1))
  x0=3.4
  y0=5.5
  npx=ns
  npy=nt
  scl=3.0
  call PS_DumpImageDistort(x0,y0,npx,npy,scl,scl)

! Bright Field image
  imaint(1:ns,1:nt)=int(255.0*BF(1:ns,1:nt,2))
  x0=0.2
  y0=2.3
  npx=ns
  npy=nt
  scl=3.0
  call PS_DumpImageDistort(x0,y0,npx,npy,scl,scl)

! Dark Field image
  imaint(1:ns,1:nt)=int(255.0*DF(1:ns,1:nt,2))
  x0=3.4
  y0=2.3
  npx=ns
  npy=nt
  scl=3.0
  call PS_DumpImageDistort(x0,y0,npx,npy,scl,scl)

! close Postscript file
  call PS_closefile

end subroutine CalcTBSM
       

