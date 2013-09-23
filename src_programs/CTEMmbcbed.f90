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
! CTEMsoft2013:CTEMmbcbed.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMmbcbed 
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Zone axis CBED
!
!> @todo OpenMP implementation; symmetry implementation
! 
!> @date  11/29/01 MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!> @date 05/14/13 MDG 2.1 replaced all IO by namelist file and added command line argument handling
!--------------------------------------------------------------------------
program CTEMmbcbed

use local
use io


IMPLICIT NONE

character(fnlen)			:: nmlfile

integer(kind=irg)			:: numarg		!< number of command line arguments
integer(kind=irg)			:: iargc		!< external function for command line
character(fnlen)    			:: arg		!< to be read from the command line

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------
!  Here is where the main program starts 
!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------
! process the command line argument (should be only one or none ...)
numarg = iargc()
if (numarg.gt.0) then ! there is an argument
        call getarg(1,arg)
        nmlfile = trim(arg)
        if (trim(nmlfile).eq.'-h') then
		mess = ' Program should be called as follows: '; call Message("(/A)")
		mess = '        CTEMmbcbed [nmlfile]'; call Message("(A)")
		mess = ' where nmlfile is an optional file name for the namelist file;'; call Message("(A)")
		mess = ' if absent, the default name ''CTEMmbcbed.nml'' will be used.'; call Message("(A/)")
		stop
	end if
else
	nmlfile = 'CTEMmbcbed.nml'    		! assign the default namelist file name
end if

! generate a set of zone axis CBED patterns
 call MBCBEDpage(nmlfile)

end program

!--------------------------------------------------------------------------
!
! SUBROUTINE:MBCBEDPage
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief draw zone axis convergent beam electron diffraction patterns
!
!> @param nmlfile namelist input file
! 
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!--------------------------------------------------------------------------
subroutine MBCBEDpage(nmlfile)

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

IMPLICIT NONE

real(kind=sgl)      			:: laL,kt,z0,thc,thb,hkl(3),ind(3),ktmax, io_real(3), voltage, convergence, &
                       				   dom,glen,qr,qi,bg,xgpz,bragg,c(3),RR,gx(3),gy(3),gg(3), thetac, startthick, thickinc, &
                       				   sc, scmax, PX, PY, qx, qy, frac, dmin
integer(kind=irg)   			:: g(3),ira,dpcnt,ppi,ijmax,ga(3),gb(3),k(3),fcnt,ccnt,cnt, skip, numthick, istat, &
                       				  newcount,count_rate,count_max, io_int(6), i, j, isym, ir, order, nn, fn(3), &
                       				  iorder, npx, npy, numt, im, numk, npix, ik, ip, jp, maxholz
character(1)        			:: ans
character(3)				:: method
character(fnlen)     			:: outname, nmlfile, xtalname


complex(kind=sgl)   		:: czero
logical             				:: overlap,np,first
real(kind=sgl),parameter     	:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
real(kind=sgl),allocatable    	:: disk(:,:,:), thick(:)
integer(kind=irg),allocatable 	:: diskoffset(:,:)
real(kind=sgl),allocatable    	:: inten(:,:)

namelist /inputlist/ stdout,xtalname, voltage, camlen, k, fn, dmin, convergence, startthick, thickinc, numthick, outname

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		! initial value to check that the keyword is present in the nml file
voltage = 200000.0			! acceleration voltage [V]
camlen = 1000.0			! camera length [mm]
k = (/ 0, 0, 1 /)			! beam direction [direction indices]
fn = (/ 0, 0, 1 /)			! foil normal [direction indices]
dmin = 0.025				! smallest d-spacing to include in dynamical matrix [nm]
convergence = 20.0			! beam convergence angle [mrad]
startthick = 10.0			! starting thickness [nm]
thickinc = 10.0			! thickness increment
numthick = 10				! number of increments
outname = 'mbcbedout.data'	! out put filename

! read the namelist file
open(UNIT=dataunit,FILE=nmlfile,DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=inputlist)
close(UNIT=dataunit,STATUS='keep')

if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMmbcbed:',' structure file name is undefined in '//nmlfile)
end if

! print some general program information (must be done here since the nmlist file 
! may have redefind the stdout entry
 progname = 'CTEMmbcbed.f90'
 progdesc = ' convergent beam pattern simulation'
 call CTEMsoft

mess = 'Input parameter list: '; call Message("(A)")
write (stdout,NML=inputlist)
mess=' '; call Message("(A/)")

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)
 skip = 3
 call CalcWaveLength(dble(voltage),skip)
 
 ! generate all atom positions
 call CalcPositions('v')
 
! get k and f
 DynFN = float(fn)
 
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
 call WriteValue('Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! construct the list of all possible reflections
method = 'ALL'
maxholz = 3
call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxholz)

! enter range of incident beam directions
  bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5

! call ReadValue(' Enter the beam convergence angle theta_c in mrad: ', io_real, 1)
  thetac = convergence/1000.0
  
! convert to ktmax along ga
  ktmax = 0.5*thetac*bragg

! compute number of pixels along diameter of central disk for given camera length
  RR = 300.0/25.4   ! dots per millimeter for 300 dots per inch; legacy code from when the output was in PostScript
  npx = nint(RR*camlen*thetac)
  npy = npx
  io_int(1) = 2.0*npx
  call WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")
  mess=' '; call Message("(A/)")
  
  ijmax = float(npx)**2   ! truncation value for beam directions
! get number of thicknesses for which to compute the CBED pattern
  numt = numthick
  allocate(thick(numt),stat=istat)
  thick = startthick + thickinc* (/ (float(i),i=0,numt-1) /)

! determine all independent incident beam directions (use a linked list starting at khead)
  call Calckvectors(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,'Conical')
 

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
  allocate(diskoffset(DynNbeamsLinked,3))
  diskoffset = 0

! compute the offset parameters for all diffraction disks (for the zone-axis orientation !!!)
! first normalize the zone axis in cartesian components; this is the z-axis
  call TransSpace(float(k),c,'d','c')
  call NormVec(c,'c')

! then make ga the x-axis
  call TransSpace(float(ga),gx,'r','c')
  call NormVec(gx,'c')

! compute the cross product between k and gx; this is the y-axis
  call CalcCross(c,gx,gy,'c','c',0)

! project every strong g-vector onto gx and gy to get the components
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    gg(1:3)=rltmpa%hkl
    call TransSpace(gg,c,'r','c')
    qx= CalcDot(c,gx,'c')*sc
    qy= CalcDot(c,gy,'c')*sc
    if ((abs(qx).lt.scmax).and.(abs(qy).lt.scmax)) then
      diskoffset(i,1) = 1
      diskoffset(i,2) = nint(qx)
      diskoffset(i,3) = nint(qy)
    else
      diskoffset(i,1) = 0
    end if
    rltmpa => rltmpa%next
  end do

  frac = 0.05

  mess = ' '; call Message("(A/)")
  io_int(1)=numk
  call WriteValue(' Starting computation for # beam directions = ', io_int, 1, "(I6)")

  call flush(stdout)
  
! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)

! point to the first beam direction
  ktmp => khead
! loop over all beam orientations, selecting them from the linked list
  do ik = 1,numk
  ! compute the dynamical matrix using Bloch waves with Bethe potentials 
   call Compute_DynMat('BLOCHBETHE', ktmp%k, .TRUE.)

! allocate the intensity array
  allocate(inten(numt,DynNbeams))
  inten = 0.0
   
! solve the dynamical eigenvalue equation
   call CalcBWint(DynNbeams,numt,thick,inten)

! and copy the intensities in the correct locations
   do i=1,DynNbeamsLinked
    if ( (diskoffset(i,1).eq.1).and.(BetheParameter%reflistindex(i).ne.0) ) then
      ip = PX + diskoffset(i,2) - ktmp%i
      jp = PX + diskoffset(i,3) + ktmp%j

! is this point inside the viewing square ?
      if (((ip.ge.1).and.(ip.le.npix)).and.((jp.ge.1).and.(jp.le.npix))) then
        do j=1,numt
         disk(j,ip,jp) = disk(j,ip,jp) + inten(j,BetheParameter%reflistindex(i))
        end do
      end if
    end if
   end do

  deallocate(inten)
  
! select next beam direction
   if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    call flush(stdout)
    frac = frac + 0.05
   end if  

  end do

! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F)")

  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
  write (dataunit) numt,npix
  write (dataunit) disk
  close(UNIT=dataunit,STATUS='keep')
  
  mess = ' Data stored in '//outname; call Message("(/A/)") 
 
end subroutine MBCBEDpage


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
  th = thick(i)
  diag(1:nn)=exp(-th*imag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
  do j=1,nn
   inten(i,j) = abs(sum(lCG(j,1:nn)*diag(1:nn)))**2
  end do 
 end do
 
end subroutine CalcBWint
