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
!> @date 09/25/13 MDG 2.2 improved command line argument handling
!> @date 09/26/13 MDG 2.3 corrected scaling error and added HOLZ
!--------------------------------------------------------------------------
program CTEMmbcbed

use local
use files
use io


IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMmbcbed.nml'
progname = 'CTEMmbcbed.f90'
call Interpret_Program_Arguments(nmldeffile,1,(/ 11 /) )

! generate a set of zone axis CBED patterns
call MBCBEDcomputation(nmldeffile)

end program

!--------------------------------------------------------------------------
!
! SUBROUTINE:MBCBEDcomputation
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief draw zone axis convergent beam electron diffraction patterns
!
!> @param nmlfile namelist input file
!
!> @todo implement symmetry
! 
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 09/25/13  MDG 2.2 minor program cleanup
!--------------------------------------------------------------------------
subroutine MBCBEDcomputation(nmlfile)

use local
use constants
use crystal
use crystalvars
use diffraction
use gvectors
use kvectors
use Lambert, ONLY: Apply2DPGSymmetry
use postscript, ONLY: GetIndex
use symmetry
use math
use dynamical
use io
use error
use files

IMPLICIT NONE

real(kind=sgl)      			:: ktmax, io_real(3), voltage, convergence, &
                       		bragg,c(3),RR,gx(3),gy(3),gg(3), thetac, startthick, thickinc, &
                       		sc, scmax, PX, qx, qy, frac, dmin, s, klaue(2)
integer(kind=irg)   			:: ijmax,ga(3),gb(3),k(3),cnt, skip, numthick, istat, dgn, badpoints, &
                       		newcount,count_rate,count_max, io_int(6), ii, i, j, isym, ir, fn(3), &
                       		npx, npy, numt, numk, npix, ik, ip, jp, maxholz, iequiv(2,12), nequiv, it
character(3)				:: method
character(fnlen)     			:: outname, nmlfile, xtalname


real(kind=sgl),parameter     	:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
real(kind=sgl),allocatable    	:: disk(:,:,:), thick(:)
integer(kind=irg),allocatable 	:: diskoffset(:,:)
real(kind=sgl),allocatable    	:: inten(:,:)
logical				:: usesym=.TRUE., bp

namelist /inputlist/ stdout,xtalname, voltage, camlen, k, fn, npix, dmin, convergence, & 
                     startthick, thickinc, numthick, outname, klaue

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		! initial value to check that the keyword is present in the nml file
voltage = 200000.0		! acceleration voltage [V]
camlen = 1000.0			! camera length [mm]
k = (/ 0, 0, 1 /)		! beam direction [direction indices]
fn = (/ 0, 0, 1 /)		! foil normal [direction indices]
klaue = (/ 0.0, 0.0 /)		! fractional coordinates of Laue circle center
npix = 900			! size of the (square) computed CBED pattern
dmin = 0.015			! smallest d-spacing to include in dynamical matrix [nm]
convergence = 20.0		! beam convergence angle [mrad]
startthick = 10.0		! starting thickness [nm]
thickinc = 10.0			! thickness increment
numthick = 10			! number of increments
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
 progdesc = ' Zone axis convergent beam pattern simulation'
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

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
 dgn = GetPatternSymmetry(k,j,.TRUE.)
 isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
 call ShortestG(k,ga,gb,isym)
 io_int(1:3)=ga(1:3)
 io_int(4:6)=gb(1:3)
 call WriteValue('Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! construct the list of all possible reflections
 method = 'ALL'
 maxholz = 0
 call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxholz)

! enter range of incident beam directions
  bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5

! call ReadValue(' Enter the beam convergence angle theta_c in mrad: ', io_real, 1)
  thetac = convergence/1000.0
  
! convert to ktmax along ga
  ktmax = 0.5*thetac/bragg

write (*,*) 'Bragg angle ', bragg
write (*,*) 'conv. angle ', thetac
write (*,*) 'ktmax ', ktmax

! compute number of pixels along diameter of central disk for given camera length
  RR = 300.0/25.4   ! dots per millimeter for 300 dots per inch; legacy code from when the output was in PostScript
  npx = int(RR*camlen*thetac)
  npy = npx
  io_int(1) = 2.0*npx
  call WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")
  mess=' '; call Message("(A/)")
  
! get number of thicknesses for which to compute the CBED pattern
  numt = numthick
  allocate(thick(numt),stat=istat)
  thick = startthick + thickinc* (/ (float(i),i=0,numt-1) /)

! if the Laue center is at the origin, then we can use symmetry groups to 
! speed up the simulation; otherwise we have to cover each incident wave
! vector separately.
  if (maxval(abs(klaue)).eq.0.0) then
    usesym=.TRUE.
  else
    usesym=.FALSE.
  end if

! determine all independent incident beam directions (use a linked list starting at khead)
  isym = WPPG(dgn)
  ijmax = float(npx)**2   ! truncation value for beam directions
  call CalckvectorsSymmetry(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue,.TRUE.)
  
! allocate the disk variable which will hold the entire computed pattern
  allocate(disk(numt,npix,npix))
  disk=0.0

  sc = mLambda*camlen*RR
  PX = npix/2
  scmax = PX + npx

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(.TRUE.)

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

! project every g-vector onto gx and gy to get the components
! and keep only the ones that will fall on the viewing region
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

  badpoints = 0

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

! allocate the intensity array to include both strong beams and weak beams (in that order)
  allocate(inten(numt,DynNbeams+BetheParameter%nnw))
  inten = 0.0
   
! solve the dynamical eigenvalue equation and return the intensities of ALL reflections,
! both strong and weak; the weak intensities should also be plotted at the correct locations....
! this uses a new version of the CalcBWint routine that implements both strong and weak beam intensities.
   call CalcBWint(DynNbeams,BetheParameter%nnw,numt,thick,inten)

! if the Bethe truncation parameters are not set correctly, then there may 
! be bad intensity points for some beams.  If this happens, we set the 
! corresponding intensities to zero and raise a flag that will cause the
! program to print a warning message at the end.  We also count the number of 
! bad points.  
  bp = .FALSE. 
  if (minval(inten).lt.0.0) then
    where (inten.lt.0.0) inten = 0.0
    bp = .TRUE.
  end if
  if (maxval(inten).gt.1.0) then
    where (inten.gt.1.0) inten = 0.0
    bp = .TRUE.
  end if
  if (bp) badpoints = badpoints+1
  
! we combine the reflistindex and weakreflistindex into one to make the CBED drawing.
BetheParameter%reflistindex = BetheParameter%reflistindex + BetheParameter%weakreflistindex
 
! and copy the strong intensities in the correct locations
   do i=1,DynNbeamsLinked
    if ( (diskoffset(i,1).eq.1).and.(BetheParameter%reflistindex(i).ne.0)) then
     if (usesym) then ! use 2D point group symmetry
      ip = diskoffset(i,2) - ktmp%i
      jp = diskoffset(i,3) + ktmp%j
      call Apply2DPGSymmetry(ip,jp,isym,iequiv,nequiv)
      iequiv = iequiv+PX

      do ii=1,nequiv
! is this point inside the viewing square ?
        ip = iequiv(1,ii)
        jp = iequiv(2,ii)
        if (((ip.ge.1).and.(ip.le.npix)).and.((jp.ge.1).and.(jp.le.npix))) then
         if ( (BetheParameter%reflistindex(i).eq.1) ) then
          disk(1:numt,ip,jp) = disk(1:numt,ip,jp) + inten(1:numt,BetheParameter%reflistindex(i))
 	 else
! this is a placeholder; it is not technically correct to do this, but the result looks quite reasonable
! it works fine when there is no overlap between diffraction disks, but when there is, then the result will
! be incorrect.  This should be noted in the manual. 
          do it=1,numt
           disk(it,ip,jp) = maxval( (/ inten(it,BetheParameter%reflistindex(i)), disk(it,ip,jp) /) )
          end do
         end if
        end if
      end do
     else  ! do not use symmetry
       ip = PX + diskoffset(i,2) - ktmp%i
       jp = PX + diskoffset(i,3) + ktmp%j
       if (((ip.ge.1).and.(ip.le.npix)).and.((jp.ge.1).and.(jp.le.npix))) then
          disk(1:numt,ip,jp) = disk(1:numt,ip,jp) + inten(1:numt,BetheParameter%reflistindex(i))
        end if
     end if
    end if
   end do
!if (maxval(disk).gt.1.0) then
!  write (*,*) 'disk larger than 1 at ik = ',ik, maxval(disk)
!  write (*,*) '   maxval(iten) = ',maxval(inten)
!  stop
!end if
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
 
  if (badpoints.gt.0) then 
    mess = '============================================================'
    call Message("(A)")
    mess = 'The program detected at least one point with invalid intensity (outside [0,1] interval)'
    call Message("(A/)")
    io_int(1) = badpoints
    io_int(2) = numk
    call WriteValue(' # bad points/ # total points = ', io_int, 2, "(I,'/',I)")
    mess = 'This may indicate that the Bethe cutoff parameters in BetheParameters.nml are too small' 
    call Message("(/A)")
    mess = 'Invalid intensities have been set to zero in the output file'
    call Message("(A)")
    mess = '============================================================'
    call Message("(A)/")    
  end if
 
 
 
end subroutine MBCBEDcomputation


!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcBWint
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute the scattered intensities for a range of thicknesses
!
!> @param nn number of strong beams
!> @param nw number of weak beams
!> @param nt number of thickness values
!> @param thick thickness array
!> @param inten output intensity list, including weak beams
!
!> @date  10/13/98 MDG 1.0 original
!> @date   7/04/01 MDG 2.0 f90
!> @date  04/29/13 MDG 3.0 inclusion of Bethe weak beams
!--------------------------------------------------------------------------
subroutine CalcBWint(nn,nw,nt,thick,inten)

use local
use io
use diffraction
use kvectors
use gvectors
use dynamical
use constants

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: nn			!< number of strong beams
integer(kind=irg),INTENT(IN)	:: nw			!< number of weak beams
integer(kind=irg),INTENT(IN)	:: nt			!< number of thickness values
real(kind=sgl),INTENT(IN)	:: thick(nt)		!< thickness array
real(kind=sgl),INTENT(INOUT)	:: inten(nt,nn+nw)	!< output intensities (both strong and weak)

integer(kind=irg)		:: i,j,IPIV(nn), ll(3), jp
complex(kind=dbl)		:: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), &
				lalpha(nn), delta(nn,nn), weak(nw,nn), Ucross(nw,nn), tmp(nw,nn), c
real(kind=sgl) 			:: th

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

! in preparation for the intensity computation, we need the prefactor array for the
! weak beam amplitude computation...
! we will also need a potential coefficient array for the cross coefficients, which we 
! compute in the same loop
 do j=1,nn   ! strong beam loop
   do jp=1,nw  ! weak beam loop
! prefactor value
     c = cmplx(2.D0*BetheParameter%weaksg(jp)/mLambda) - 2.D0*ktmp%kn*Wloc(j)
     weak(jp,j) = cmplx(-1.D0,0.D0)/c
! cross potential coefficient
     ll(1:3) = BetheParameter%weakhkl(1:3,jp) - BetheParameter%stronghkl(1:3,j)
     Ucross(jp,j) = LUT( ll(1),ll(2),ll(3) )
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
   inten(i,j) = cabs(sum(delta(j,1:nn)))**2
  end do 
! weak beams
  tmp = matmul(Ucross,delta)
  do jp=1,nw
   inten(i,nn+jp) = cabs( sum(weak(jp,1:nn)*tmp(jp,1:nn)) )**2
  end do  
 end do
   
end subroutine CalcBWint
