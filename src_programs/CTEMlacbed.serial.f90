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
! CTEMsoft2013:CTEMlacbed.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMlacbed 
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Zone axis LACBED
!
!> @todo implement full symmetry use; implement OpenMP; implement multiple reflection output
!> or, easier perhaps, selection of one reflection; a nice addition would be the ability to define a 
!> foil thickness profile that can be read from an external file. For each beam direction, the 
!> thickness would be slightly different, so that one can obtain more realistic LACBED patterns.
!> Usually, the beam illuminates a large area of the sample, sothe thickness can vary quite a 
!> bit; Joachim Mayer's [111] Si pattern at 100 kV is a nice example of that (page v of Spence&Zuo)
!
!> @date  11/29/01 MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!> @date 05/08/13 MDG 2.1 forked from mbcbed and adapted for large angle CBED patterns
!--------------------------------------------------------------------------
program CTEMlacbed

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

IMPLICIT NONE

real(kind=sgl)			:: io_real(1)

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData
 call GetVoltage

 call ReadValue(' Camera length L  [mm, real] ', io_real, 1)
 camlen = io_real(1)

! generate all atom positions
 call CalcPositions('v')

! generate a set of zone axis CBED patterns
 call LACBEDpattern

end program

!--------------------------------------------------------------------------
!
! SUBROUTINE:LACBEDpattern
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute a large angle zone axis convergent beam electron diffraction pattern
!
! 
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/13/13  MDG 3.0 OpenMP implementation added 
!--------------------------------------------------------------------------
subroutine LACBEDpattern

use local
use constants
use crystal
use crystalvars
use diffraction
use gvectors
use kvectors
use postscript, ONLY: GetIndex
use symmetry
use math
use dynamical
use io
use error
use files
use omp_lib

IMPLICIT NONE

real(kind=sgl)      			:: laL,kt,z0,thc,thb,hkl(3),ind(3),fn(3),ktmax, io_real(3), &
                       				   dom,glen,qr,qi,bg,xgpz,bragg,c(3),RR,gx(3),gy(3),gg(3), thetac, thick, &
                       				   sc, scmax, qx, qy, frac, dmin
integer(kind=irg)   			:: g(3),ira,dpcnt,ppi,ijmax,ga(3),gb(3),k(3),fcnt,ccnt,cnt, PX, PY, &
                       				  newcount,count_rate,count_max, io_int(6), i, j, isym, ir, order, nn, &
                       				  iorder, npx, npy, numt, im, numk, npix, ik, ip, jp, maxholz, numpix, istat
character(1)        			:: ans
character(2)        			:: srza
character(3)				:: method
character(fnlen)     			:: fname, tname


complex(kind=sgl)   		:: czero
logical             				:: overlap,np,first
real(kind=sgl),parameter     	:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
real(kind=sgl),allocatable    	:: disk(:,:,:), thickarray(:,:)
integer(kind=irg),allocatable 	:: diskoffset(:,:)
real(kind=sgl),allocatable    	:: inten(:,:)

! zone axis computation
 srza = 'ZA'

! get k and f
 mess = ' Enter wave vector direction'; call Message("(/A)")
 call GetIndex(k,'d')
 mess = ' Enter foil normal'; call Message("(/A)")
 call GetIndex(g,'d')
 DynFN = float(g)
 
! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do
! and the Bright Field symmetry
 call BFsymmetry(k,j,isym,ir)
 io_int(1:3) = k(1:3)
 call WriteValue('', io_int, 3, "(//,' ','[',3I2,'] has Bright Field symmetry ',$)")
 mess = PGTWD(isym)
 call Message("(A,$)")
 io_int(1) = ir 
 call WriteValue('; order = ',io_int, 1, "(I4,//)")


! determine the shortest reciprocal lattice points for this zone
 call ShortestG(k,ga,gb,isym)
 io_int(1:3)=ga(1:3)
 io_int(4:6)=gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! construct the list of all possible reflections
dmin = 0.025
method = 'ALL'
maxholz = 0
call Compute_ReflectionList(dmin,k,ga,gb,method,maxholz)

! determine all families of reciprocal lattice points
! that belong to the zone/row and rank them
!  call RankReflections(k,ga,gb,fcnt,srza,iorder)

! then select the reflections needed for the simulation
!  call SelectReflections(fcnt,nn,ccnt)

! enter range of incident beam directions
  mess = ' The program will use a symmetric cone of beam directions,'; call Message("(/A)")
  mess = ' centered on the incident beam direction entered above.'; call Message("(A)")
  bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5*1000.0
  io_real(1) = bragg
  call WriteValue(' The Bragg angle for the first reflection is equal to [mrad]: ',io_real, 1, "(F8.5)")

  call ReadValue(' Enter the beam convergence angle theta_c in mrad: ', io_real, 1)
  thetac = io_real(1)/1000.0
  
! convert to ktmax along ga
  ktmax = 0.5*thetac*1000.0/bragg

! compute number of pixels along diameter of central disk for given camera length
  RR = 300.0/25.4   ! dots per millimeter for 300 dots per inch
  npx = nint(RR*camlen*thetac)
  npy = npx
  io_int(1) = 2.0*npx
  call WriteValue(' Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")

  ijmax = float(npx)**2   ! truncation value for beam directions
! get number of thicknesses for which to compute the CBED pattern
  mess = ' This program computes LACBED patterns either for a series of thicknesses'; call Message("(A)")
  mess = ' or for a single variable thickness array; enter a negative thickness for second option'; call Message("(A)")
  call ReadValue(' Enter the first thickness [nm, R]', io_real, 1)
  if (io_real(1).gt.0.0) then
    thick = io_real(1)
    call ReadValue(' How many multiples of this thickness [I] ', io_int, 1) 
    numt = io_int(1)
  else
     call ReadValue(' Enter filename for thickness array : ', tname, '(A)')
     open(UNIT=dataunit,FILE=trim(tname),STATUS='old',FORM='unformatted')
     allocate(thickarray(900,900),stat=istat)
     read(UNIT=dataunit) numpix
     if (numpix.ne.900) call FatalError('LACBED ',' thickness array dimensions do not match required size')
     read(UNIT=dataunit) thickarray
     close(UNIT=dataunit,STATUS='keep')
     numt = 1
     thick = 0.0
  end if  
! determine all independent incident beam directions (use a linked list starting at khead)
! we're using the StandardConical method since we do want to make use of symmetry
! operations, but the illumination must be limited to a cone.  
!  call Calckvectors(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,'StandardConical')
  call Calckvectors(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,'Conical')
 

! allocate the disk variable which will hold the entire computed pattern
  npix = int(3.0*300.0)   ! 3 inches wide at 300 dpi
  allocate(disk(numt,npix,npix))
  disk=0.0

  sc = mLambda*camlen*RR
  scmax = 1.5*300.0 + npx
  PX = npix/2


! allocate the offset array
! to get this array, we need to do a mock initialization of the dynamical matrix in zone axis orientation
  call Compute_DynMat('BLOCHBETHE', khead%k, .TRUE.)
  frac = 0.05

  io_int(1)=numk
  call WriteValue(' Starting computation for # beam directions = ', io_int, 1, "(I8)")


! if we want to use OpenMP, then we need to first convert the wave vector information
! into regular allocatable arrays that can then become SHARED access to all threads.
!
! we need the following items from the khead linked pointer list:
! i, j, k

! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)



! point to the first beam direction
  ktmp => khead
! loop over all beam orientations, selecting them from the linked list
  do ik = 1,numk

! is this point inside the viewing square ?
      ip = PX  - ktmp%i
      jp = PX + ktmp%j

      if (((ip.ge.1).and.(ip.le.npix)).and.((jp.ge.1).and.(jp.le.npix))) then

  ! compute the dynamical matrix using Bloch waves with Bethe potentials 
        call Compute_DynMat('BLOCHBETHE', ktmp%k, .TRUE.)

! allocate the intensity array
        allocate(inten(numt,DynNbeams))
        inten = 0.0

! solve the dynamical eigenvalue equation
       if (allocated(thickarray)) then 
         call CalcBWint(DynNbeams,numt,thickarray(PX  - ktmp%i,PX + ktmp%j), inten)
      else
        call CalcBWint(DynNbeams,numt,thick,inten)
      end if

! copy in the correct locations
      do j=1,numt
         disk(j,ip,jp) = disk(j,ip,jp) + inten(j,1)
      end do

! and remove the intensity array
     deallocate(inten)
    
   end if

! select next beam direction
   if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

  end do



write (*,*) ' maximum intensity ', maxval(disk)
write (*,*) ' minimum # strong beams ',BetheParameter%minstrong
write (*,*) ' maximum # strong beams ',BetheParameter%maxstrong

! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  call WriteValue(' Total computation time [s] ' , io_real, 1, "(F)")

 loadingfile = .FALSE.
  call SafeOpenFile('d1','unformatted',fname)
  write (dataunit) numt,npix
  write (dataunit) disk
  call SafeCloseFile('d1','keep',fname)

end subroutine LACBEDpattern


! ###################################################################
! 
!  subroutine CalcBW
! 
!  Author: Marc De Graef
!  
!  Description: integrate the dynamical equations using the Bloch Wave
!  formalism.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   7/04/01 MDG 2.0 f90
! ###################################################################
subroutine CalcBWint(nn,nt,thick,inten)

use local
use io
use diffraction
use kvectors
use dynamical
use constants

IMPLICIT NONE

integer(kind=irg)                             	:: nn,i,j,ig,ih,IPIV(nn),nt
complex(kind=dbl)    				:: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), lalpha(nn)
real(kind=sgl)                             		:: thick,inten(nt,nn),th
! complex(kind=dbl)                			:: Ijk(nn,nn),Lgh(nn,nn),q

intent(IN)      :: nt, nn, thick
intent(OUT)     :: inten

! should we just return zero ?
if (nt.eq.1) then
  if (thick.eq.0.0) then 
    inten = 0.0
    return
  end if
end if

! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
 lW = cPi*Wloc/cmplx(ktmp%kn,0.0)
 do i=1,nn
  lalpha(i) = CGinv(i,1)
 end do

 do i=1,nt
  th = dble(thick*float(i))
  diag(1:nn)=exp(-th*imag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
  do j=1,nn
   inten(i,j) = abs(sum(lCG(j,1:nn)*diag(1:nn)))**2
  end do 
 end do
 
end subroutine CalcBWint
