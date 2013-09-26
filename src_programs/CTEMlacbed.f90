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
!> @todo implement full symmetry use; implement full Bloch wave output
!>
!> implement OpenMP multithreading for the actual computation part; requires modifications
!> in CTEMlib.a routines (mostly THREADPRIVATE commands in several modules)
!
!> @date 11/29/01 MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!> @date 05/08/13 MDG 2.1 forked from mbcbed and adapted for large angle CBED patterns
!> @date 05/14/13 MDG 2.2 replaced all IO by namelist file and added command line argument handling
!> @date 09/04/13 MDG 2.3 all command line argument handling now via files.f90 routine
!--------------------------------------------------------------------------
program CTEMlacbed

use local
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMlacbed.nml'
progname = 'CTEMlacbed.f90'
call Interpret_Program_Arguments(nmldeffile,1,(/ 10 /) )

! perform the zone axis computations
call LACBEDpattern(nmldeffile)

end program CTEMlacbed

!--------------------------------------------------------------------------
!
! SUBROUTINE:LACBEDpattern
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute a large angle zone axis convergent beam electron diffraction pattern
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!--------------------------------------------------------------------------
subroutine LACBEDpattern(nmlfile)

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

character(fnlen),INTENT(IN)	:: nmlfile

real(kind=sgl)      		:: laL,kt,z0,thc,thb,hkl(3),ind(3),ktmax, io_real(3), thickval(1), &
                       	   dom,glen,qr,qi,bg,xgpz,bragg,c(3),RR,gx(3),gy(3),gg(3), thetac, &
                       	   sc, scmax, qx, qy, frac, dmin, convergence, voltage, startthick, thickinc
integer(kind=irg)   		:: g(3),ira,dpcnt,ppi,ijmax,ga(3),gb(3),k(3),fcnt,ccnt,cnt,fn(3), PX, PY, numthick,  &
                      		   ii, newcount,count_rate,count_max, io_int(6), i, j, isym, ir, order, nn, skip, refpick, &
                      		  iorder, npx, npy, numt, im, numk, npix, ik, ip, jp, maxholz, numpix, istat, selection(3)
character(1)        		:: ans
character(3)			:: method
character(fnlen)     		:: outname, tname, xtalname


complex(kind=sgl)   		:: czero
logical             		:: overlap,np,first
real(kind=sgl),parameter     	:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
real(kind=sgl),allocatable    	:: disk(:,:,:), thickarray(:,:), thick(:)
integer(kind=irg),allocatable 	:: diskoffset(:,:)
real(kind=sgl),allocatable    	:: inten(:,:)


namelist /inputlist/ stdout, xtalname, voltage, camlen, k, fn, dmin, convergence, &
                              startthick, thickinc, numthick, tname, outname, selection

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		! initial value to check that the keyword is present in the nml file
voltage = 200000.0			! acceleration voltage [V]
camlen = 1000.0			! camera length [mm]
k = (/ 0, 0, 1 /)			! beam direction [direction indices]
fn = (/ 0, 0, 1 /)			! foil normal [direction indices]
dmin = 0.025				! smallest d-spacing to include in dynamical matrix [nm]
convergence = 20.0			! beam convergence angle [mrad]
startthick = 10.0			! starting thickness [nm]; if negative, then thickness profile file will be read
tname = 'thickness.data'		! thickness profile filename, used if startthick<0
thickinc = 10.0			! thickness increment
numthick = 10				! number of increments
selection = (/ 0, 0, 0/)		! default reflection to use for LACBED pattern
outname = 'lacbedout.data'	! output filename

! read the namelist file
open(UNIT=dataunit,FILE=nmlfile,DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=inputlist)
close(UNIT=dataunit,STATUS='keep')

if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMlacbed:',' structure file name is undefined in '//nmlfile)
end if

! print some information
 progname = 'CTEMlacbed.f90'
 progdesc = 'Large angle convergent beam pattern simulation'
 call CTEMsoft

 mess = 'Input parameter list: '; call Message("(A)")
 write (*,NML=inputlist)

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

! generate all atom positions
 call CalcPositions('v')
 
! set the foil normal 
 DynFN = float(fn)
 
! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do
! and the Bright Field symmetry, and print some information
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


! determine range of incident beam directions
  bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5*1000.0

!  call ReadValue(' Enter the beam convergence angle theta_c in mrad: ', io_real, 1)
  thetac = convergence/1000.0
  
! convert to ktmax along ga
  ktmax = 0.5*thetac*1000.0/bragg

! compute number of pixels along diameter of central disk for given camera length
  RR = 300.0/25.4   ! dots per millimeter for 300 dots per inch
  npx = nint(RR*camlen*thetac)
  npy = npx
  io_int(1) = 2.0*npx
  call WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")

  ijmax = float(npx)**2   ! truncation value for beam directions
! thickness information  [this needs to be changed for Bloch output]
  if (startthick.gt.0.0) then
    mess = ' This program will compute LACBED patterns for a series of thicknesses.'; call Message("(/A/)")
    numt = numthick
    allocate(thick(numt),stat=istat)
    thick = startthick + thickinc* (/ (float(i),i=0,numt-1) /)
  else
     mess = ' This program will compute an LACBED pattern for a specific thicknesses profile.'; call Message("(/A/)")
! this hasn't been programmed yet; will require a little work; make dims variable
     open(UNIT=dataunit,FILE=trim(tname),STATUS='old',FORM='unformatted')
     allocate(thickarray(900,900),stat=istat)
     read(UNIT=dataunit) numpix
     if (numpix.ne.900) call FatalError('LACBED ',' thickness array dimensions do not match required size')
     read(UNIT=dataunit) thickarray
     close(UNIT=dataunit,STATUS='keep')
     numt = 1
  end if  

! determine all independent incident beam directions (use a linked list starting at khead)
! we're using the Conical method since we do want to make use of symmetry
! operations, but the illumination must be limited to a cone.  
  call Calckvectors(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,'Conical')


! construct the list of all possible reflections
! this needs a little work because we really need to create a master list with ALL
! reflections that will be needed for at least one beam direction; furthermore, this 
! requires figuring out which reflections are strong or weak or absent.
method = 'ALL'
maxholz = 0	! not used
call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxholz)
 


! allocate the disk variable which will hold the entire computed pattern
  npix = int(3.0*300.0)   ! 3 inches wide at 300 dpi
  allocate(disk(numt,npix,npix))
  disk=0.0

  sc = mLambda*camlen*RR
  scmax = 1.5*300.0 + npx
  PX = npix/2

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters()

! allocate the offset array
! to get this array, we need to do a mock initialization of the dynamical matrix in zone axis orientation
  call Compute_DynMat('BLOCHBETHE', khead%k, .TRUE.)
  frac = 0.05

  io_int(1)=numk
  call WriteValue(' Starting computation for # beam directions = ', io_int, 1, "(I8)")

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
!write (*,*) ik, DynNbeams

! solve the dynamical eigenvalue equation
       if (allocated(thickarray)) then 
         thickval(1) = thickarray(PX  - ktmp%i,PX + ktmp%j)
         call CalcBWint(DynNbeams,numt,thickval, inten)
      else
        call CalcBWint(DynNbeams,numt,thick,inten)
      end if

! first determine which reflection we are going to plot the pattern for.
! normally, it's the bright field pattern, so refpick = 1, but if selection
! is not (/ 0,0,0 /), then we need to first figure out what the number
! is of the selected reflection in the main list, and then we need to see if
! that reflection was actually used in the Bloch wave computation...
	if (sum(abs(selection)).eq.0) then
	  refpick = 1
	else
	  refpick = -1
	  do ii = 1, BetheParameter%nns
	    if (sum(abs(selection(1:3)-BetheParameter%stronghkl(1:3,ii))).eq.0) then
	     refpick = ii
	    end if
	  end do
!	  if (refpick.gt.0) write (*,*) ik, BetheParameter%stronghkl(1:3,refpick), refpick, &
!	  BetheParameter%reflistindex(refpick)
	end if


! copy in the correct locations
      if (refpick.eq.-1) then
       disk(1:numt,ip,jp) = 0.0
     else
!      ii = BetheParameter%reflistindex(refpick)
      do j=1,numt
         disk(j,ip,jp) = disk(j,ip,jp) + inten(j,refpick)
      end do
     end if
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

!write (*,*) ' maximum intensity ', maxval(disk)
!write (*,*) ' minimum # strong beams ',BetheParameter%minstrong
!write (*,*) ' maximum # strong beams ',BetheParameter%maxstrong

! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F)")

  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
!  loadingfile = .FALSE.
!  call SafeOpenFile('d1','unformatted',outname)
  write (dataunit) numt,npix
  write (dataunit) disk
!  call SafeCloseFile('d1','keep',outname)
  close(UNIT=dataunit,STATUS='keep')
  
  mess = ' Data stored in '//outname; call Message("(/A/)") 
 
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
real(kind=sgl)                             		:: thick(nt),inten(nt,nn),th
! complex(kind=dbl)                			:: Ijk(nn,nn),Lgh(nn,nn),q

intent(IN)      :: nt, nn, thick
intent(OUT)     :: inten

! should we just return zero ?
if (nt.eq.1) then
  if (thick(1).eq.0.0) then 
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
  th =thick(i)
  diag(1:nn)=exp(-th*imag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
  do j=1,nn
   inten(i,j) = abs(sum(lCG(j,1:nn)*diag(1:nn)))**2
  end do 
 end do
 
end subroutine CalcBWint
