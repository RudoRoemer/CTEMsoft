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
! CTEMsoft2013:CTEMped.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMped 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Zone axis Precession Electron Diffraction
!
!> @todo implement full symmetry use; implement full Bloch wave output
!>
!
!> @date 11/29/01 MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!> @date 05/08/13 MDG 2.1 forked from mbcbed and adapted for large angle CBED patterns
!> @date 05/14/13 MDG 2.2 replaced all IO by namelist file and added command line argument handling
!> @date 09/04/13 MDG 2.3 all command line argument handling now via files.f90 routine
!--------------------------------------------------------------------------
program CTEMped

use local
use files
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMped.nml'
progname = 'CTEMped.f90'
call Interpret_Program_Arguments(nmldeffile,1,(/ 13 /) )

! perform the zone axis computations
call PEDpattern(nmldeffile)

end program CTEMped

!--------------------------------------------------------------------------
!
! SUBROUTINE:PEDpattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a zone axis precession electron diffraction pattern
!
!> @param nmlfile namelist file name
!
!> @date 03/05/14  MDG 1.0 original, based on CTEMlacbed program
!--------------------------------------------------------------------------
subroutine PEDpattern(nmlfile)

use local
use constants
use crystal
use crystalvars
use diffraction
use gvectors
use kvectors
use MBmodule
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

real(kind=sgl)      		:: ktmax, io_real(3), bragg, thetac, sc, minten, pxy(2), galen, prechalfwidth, DM(2,2), DD, X(2), &
                       	   frac, dmin, precangle, voltage, startthick, thickness, thick(1), klaue(2), thetam 
integer(kind=irg)   		:: ijmax,ga(3),gb(3),k(3),cnt,fn(3), PX, numthick, ss, icnt, pgnum, ih, nunique, famnum, &
                      		   newcount,count_rate,count_max, io_int(6), i, j, isym, ir, skip, ghkl(3), &
                      		   npx, npy, numt, numk, npix, ik, ip, jp, istat, dgn, nbeams, refcnt, &
                      		   ifamily, famhkl(3), inum, maxHOLZ, numksame, precazimuthal, precsample
character(3)			:: method
character(fnlen)     		:: outname, xtalname
character(5)                   :: filemode

real(kind=sgl),allocatable    	::  disk(:,:)
integer(kind=irg),allocatable 	:: familymult(:), familyhkl(:,:), whichHOLZ(:), gequiv(:,:)
real(kind=sgl),allocatable    	:: inten(:,:)
real(kind=dbl)			:: s(3)
logical,allocatable		:: ksame(:)


namelist /inputlist/ stdout, xtalname, voltage, k, fn, dmin, precangle, prechalfwidth, precsample, precazimuthal, &
                              thickness,  outname, npix, camlen, filemode

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		! initial value to check that the keyword is present in the nml file
stdout = 6			! standard output
voltage = 200000.0		! acceleration voltage [V]
k = (/ 0, 0, 1 /)		! beam direction [direction indices]
fn = (/ 0, 0, 1 /)		! foil normal [direction indices]
dmin = 0.025			! smallest d-spacing to include in dynamical matrix [nm]
precangle = 10.472		! beam precession angle [mrad]; default = 0.6 degrees
prechalfwidth = 0.25		! beam half width in the tilt direction [mrad]
precsample = 10		        ! number of samples (concentric circles) in beam half width (total = 2*precsample + 1)
precazimuthal = 360		! number of azimuthal samples for each precession circle
thickness = 10.0		! sample thickness [nm]
filemode = 'total'            ! 'total' mode or 'eachp'
npix = 256			! output arrays will have size npix x npix
outname = 'pedout.data'	! output filename
camlen = 1000.0			! camera length [mm]

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=inputlist)
close(UNIT=dataunit,STATUS='keep')

if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMped:',' structure file name is undefined in '//nmlfile)
end if

! print some information
 progname = 'CTEMped.f90'
 progdesc = 'Zone axis precession electron diffraction pattern simulation'
 call CTEMsoft

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

! generate all atom positions
 call CalcPositions('v')
 
! set the foil normal 
 DynFN = float(fn)
 numt = 1
 thick(1) = thickness
 
! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
 dgn = GetPatternSymmetry(k,j,.TRUE.)
 pgnum = j
 isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
 call ShortestG(k,ga,gb,isym)
 io_int(1:3)=ga(1:3)
 io_int(4:6)=gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! construct the list of all possible reflections
 method = 'ALL'
 thetac = precangle/1000.0
 call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxHOLZ,thetac)
 galen = CalcLength(float(ga),'r')

! determine the list of contributing wave vectors
  call CalckvectorsPrecession(dble(k),dble(ga),precangle,prechalfwidth,precsample,precazimuthal,numk)

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(.TRUE.)

! up to this point, everything is nearly identical to the mbcbed program,
! except that we do not use a camera length explicitly.  Now we
! need to do things a little differently.  First of all, we need a master 
! list for all the reflections that contribute in one way or another to the
! complete PED pattern.  We do this by going through the entire list of incident
! beam directions and flagging all reflections that contribute, either as weak
! or as strong reflections.
  mess = ' Pruning reflection list (this takes a while ...) '
  call Message("(A)")
  call Prune_ReflectionList(numk,nbeams)
  io_int(1) = nbeams
  call WriteValue('Number of contributing beams  : ', io_int, 1, '(I)')


! set the scale parameter for a default camera length of 1000 mm.
  sc = mLambda * 1000.0 * 300.0 / 25.4  ! the absolute value does not matter and is derived from legacy Postscript code
! The original code used 300 dpi (hence 300/25.4) which was convenient for Postscript output; in the current case, we
! do not actually use the true value, but in the IDL visualization program, we scale the user defined camera length by
! 1000.0, and use this ratio to scale the diskoffset coordinates.  So, there's no absolute length scale, only a relative scale.


! next we create the output array (not sure yet how to do this)
  allocate(disk(-npix:npix,-npix:npix))
  disk=0.0

  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")

! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)

! go through the linked reflection list and set all intensities to zero; we'll use the 
! xg entry in the rltmpa linked list to accumulate the intensities

if (filemode.eq.'total') then
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    rltmpa%xg = 0.D0
    rltmpa => rltmpa%next
  end do
end if


! the final bit of the program involves dumping all the results into a file,
! binary for now, but HDF5 in the future, for the IDL visualization program 
! to read.
  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
! write the program identifier
  write (dataunit) trim(progname)
! write the version number
  write (dataunit) scversion
! then the name of the crystal data file
  write (dataunit) xtalname
! the accelerating voltage [V]
  write (dataunit) voltage
! precangle angle [mrad]
  write (dataunit) precangle
! prechalfwidth angle [mrad]
  write (dataunit) prechalfwidth
! precsample
  write (dataunit) precsample
! precazimuthal
  write (dataunit) precazimuthal
! the zone axis indices
  write (dataunit) k
! the foil normal indices
  write (dataunit) fn
! number of k-values in disk
  write (dataunit) numk
! dmin value
  write (dataunit) dmin
! horizontal reciprocal lattice vector
  write (dataunit) ga  
! second reciprocal lattice vector
  write (dataunit) gb
! length horizontal reciprocal lattice vector (need for proper Laue center coordinate scaling)
  write (dataunit) galen
! length second reciprocal lattice vector (need for proper Laue center coordinate scaling)
  write (dataunit) CalcLength(float(gb),'r')
! angle between these two vectors (in radians)
  write (dataunit) CalcAngle(float(ga),float(gb),'r')
! thickness data
  write (dataunit) thickness
  

  
 DM(1,1) = CalcDot(float(gb),float(gb),'c')
 DM(1,2) = -CalcDot(float(ga),float(gb),'c')
 DM(2,1) = DM(1,2)
 DM(2,2) = CalcDot(float(ga),float(ga),'c')
 DD = 1.0/(DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1))

  


! here we go...

! point to the first beam direction
  ktmp => khead%next
! loop over all beam orientations, selecting them from the linked list
kvectorloop:  do ik = 1,numk

! compute the dynamical matrix using Bloch waves with Bethe potentials; note that the IgnoreFoilNormal flag
! has been set to .FALSE.; if it is set to .TRUE., the computation of the ZOLZ will still be mostly correct,
! but the excitation errors of the HOLZ reflections will be increasingly incorrect with HOLZ order.  This was
! useful during program testing but should probably be removed as an option altogether...
	call Compute_DynMat('BLOCHBETHE', ktmp%k, ktmp%kt, .FALSE.)

! allocate the intensity array to include both strong beams and weak beams (in that order)
	allocate(inten(numt,DynNbeams+BetheParameter%nnw))
  	inten = 0.0
 
! solve the dynamical eigenvalue equation and return the intensities of ALL reflections,
! both strong and weak; the weak intensities should also be plotted at the correct locations....
! this uses a new version of the CalcBWint routine that implements both strong and weak beam intensities.
   	call CalcBWint(DynNbeams,BetheParameter%nnw,numt,thick,inten)

! we combine the reflistindex and weakreflistindex into a single list so that we can match 
! each diffracted intensity with the correct diffraction disk in the disks array.
	BetheParameter%reflistindex = BetheParameter%reflistindex + BetheParameter%weakreflistindex

! ok, we have all the intensities.  Next we need to copy the relevant intensities into the 
! xg field of the rltmpa linked list; we simply add the intensities together
! of the disk array, one for each family.
	rltmpa => reflist%next
	do i=1,DynNbeamsLinked
	  if (BetheParameter%reflistindex(i).ne.0) then ! is this a reflection on the current list
	    rltmpa%xg = rltmpa%xg + inten(1,BetheParameter%reflistindex(i))
	  end if
	  rltmpa => rltmpa%next
	end do

! and remove the intensity array
     deallocate(inten)

! remove all the computed intensities for per pattern storage    
if (filemode.eq.'eachp') then

! first we count how many reflections have none-zero intensity  
  refcnt = 0
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    if (rltmpa%xg.ne.0.D0) refcnt = refcnt + 1
    rltmpa => rltmpa%next
  end do

! write refcnt
  write (dataunit) refcnt
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    if (rltmpa%xg.ne.0.D0) then
! decompose this point w.r.t ga and gb
      X(1) = CalcDot(float(rltmpa%hkl),float(ga),'c')
      X(2) = CalcDot(float(rltmpa%hkl),float(gb),'c')
      X = matmul(DM,X) * DD
      write (dataunit) rltmpa%hkl
      write (dataunit) rltmpa%xg
      write (dataunit) X
    end if
    rltmpa => rltmpa%next
  end do
  
! and reset the intensities for the next run
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    rltmpa%xg = 0.D0
    rltmpa => rltmpa%next
  end do
end if

! select next beam direction
   if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

  end do kvectorloop



! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F)")


! remove all the computed intensities for per pattern storage    
if (filemode.eq.'total') then

! first we count how many reflections have none-zero intensity  
  refcnt = 0
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    if (rltmpa%xg.ne.0.D0) refcnt = refcnt + 1
    rltmpa => rltmpa%next
  end do

! write refcnt
  write (dataunit) refcnt
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    if (rltmpa%xg.ne.0.D0) then
! decompose this point w.r.t ga and gb
      X(1) = CalcDot(float(rltmpa%hkl),float(ga),'c')
      X(2) = CalcDot(float(rltmpa%hkl),float(gb),'c')
      X = matmul(DM,X) * DD
      write (dataunit) rltmpa%hkl
      write (dataunit) rltmpa%xg
      write (dataunit) X
    end if
    rltmpa => rltmpa%next
  end do
end if
  
  close(UNIT=dataunit,STATUS='keep')

  mess = ' Data stored in '//outname; call Message("(/A/)") 

end subroutine PEDpattern




!--------------------------------------------------------------------------
!
! SUBROUTINE: CalckvectorsPrecession
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create a linked list of wave vectors for precession electron diffraction
!
!> @details This is a new version to test whether or not we can use the whole pattern
!> symmetry to determine the relevant list of incident wave vectors; this should be a 
!> general routine, so that we do not need to consider each symmetry case separately.
!> This will require a floating point version of the Apply2DPGSymmetry routine in symmetry.f90.
!
!> @param k central wave vector
!> @param ga reciprocal lattice vector normal to k
!> @param precangle main precession angle [mrad]
!> @param prechalfwidth precession beam half width [mrad]
!> @param precsample number of samples in half width 
!> @param precazimuthal number of samples around each precession circle
!> @param numk total number of wave vectors in list
!
!> @date  03/05/14 MDG 1.0 original
!--------------------------------------------------------------------------
recursive subroutine CalckvectorsPrecession(k,ga,precangle,prechalfwidth,precsample,precazimuthal,numk)

use local
use io
use error
use constants
use diffraction
use crystal
use crystalvars
use Lambert
use kvectors

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: k(3)		!< initial wave vector
real(kind=dbl),INTENT(IN)		:: ga(3)	!< "horizontal" reciprocal lattice vector
real(kind=sgl),INTENT(IN)		:: precangle	!< precession angle in [mrad]
real(kind=sgl),INTENT(IN)		:: prechalfwidth	!< halfwidth of tilted beam [mrad]
integer(kind=irg),INTENT(IN)		:: precsample		!< number of kvectors along beam tilt
integer(kind=irg),INTENT(IN)		:: precazimuthal	!< number of kvectors along circumference
integer(kind=irg),INTENT(OUT)		:: numk		!< total number of kvectors in linked list

integer(kind=irg)       		:: istat,i,j, iequiv(2,12), nequiv, jj, nx, ny, il, ith
real(kind=dbl)				:: gp, dgp, glen, gan(3), gperp(3), kstar(3), dth
real(kind=dbl),allocatable		:: gw(:), ct(:), st(:), th(:)
logical					:: hexgrid = .FALSE.
real(kind=sgl)				:: kt(3),kr(3)
real(kind=sgl)				:: ktlen

write (*,*) precangle, prechalfwidth, precsample, precazimuthal, mLambda

! compute geometrical factors 
 glen = CalcLength(ga,'r')              		! length of ga
 gan = ga/glen                                 	! normalized ga
 gp = 2.0*sin(precangle/1000.0)/mLambda		! precession angle converted to reciprocal length gp in units of glen
 dgp = 0.0
 if (precsample.gt.0) then
   dgp = 2.0*sin(0.001*(precangle-prechalfwidth))/mLambda/glen/float(precsample)	! half width step size converted to reciprocal length dgp in units of glen
 end if
 allocate(gw(2*precsample+1))				! sampling radii
 gw = gp + dgp * (/ (i,i=-precsample,precsample) /)	! sampling radii

! pre-compute cosines and sines
 allocate(ct(precazimuthal),st(precazimuthal), th(precazimuthal))
 dth = 2.D0*cPi / dble(precazimuthal)
 th = (/ (i-1,i=1,precazimuthal) /) * dth
 ct = cos(th)
 st = sin(th)
 
 call TransSpace(k,kstar,'d','r')       		! transform incident direction to reciprocal space
 call CalcCross(ga,kstar,gperp,'r','r',0)      	! compute g_perp = ga x k
 call NormVec(gperp,'r')                       	! normalize g_perp
 call NormVec(kstar,'r')                       	! normalize reciprocal beam vector

! allocate the head and tail of the linked list
 allocate(khead,stat=istat)   				! allocate new value
 if (istat.ne.0) call FatalError('CalckvectorsPrecession','unable to allocate khead pointer')
 ktail => khead                      			! tail points to new value
 nullify(ktail%next)                			! nullify next in new value
 numk = 0                          			! keep track of number of k-vectors so far

 
! next loop around each of the precession circles
 do il = 1,2*precsample+1  				! number of concentric circles
  do ith = 1,precazimuthal  				! number of points along each circle
! make a new one in the list, except for the first one
   if (numk.ne.0) then
     allocate(ktail%next,stat=istat)  			! allocate new value
     if (istat.ne.0) call FatalError('Add_knode:',' unable to allocate pointer')
     ktail => ktail%next               		! tail points to new value
     nullify(ktail%next)              			! nullify next in new value
   end if
! and populate the fields   
   kt = - gw(il)*ct(ith)*gan - gw(il)*st(ith)*gperp  	! tangential component of k
   ktail%kt = kt                    			! store tangential component of k
   ktlen = CalcLength(kt,'r')**2      			! squared length of tangential component

   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar 	! complete wave vector
   ktail%k = kr                     			! store in pointer list
   ktail%kn = CalcDot(ktail%k,kstar,'r')    		! normal component of k
   numk = numk + 1
  end do
 end do


end subroutine CalckvectorsPrecession

