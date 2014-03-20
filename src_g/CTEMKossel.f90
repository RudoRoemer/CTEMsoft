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
! CTEMsoft2013:CTEMKossel.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMKossel 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Kossel patterns, taken from master EBSD pattern simulation, with adjustment for 
!> Bloch wave part ...
!
!> @todo implement OpenMP multithreading for the actual computation part; requires modifications
!> in CTEMlib.a routines (mostly THREADPRIVATE commands in several modules)
!
!> @date 04/08/11 MDG 1.0 early version
!> @date 01/07/14 MDg 2.0 new version
!--------------------------------------------------------------------------
program CTEMKossel

use local
use files
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMKossel.nml'
progname = 'CTEMKossel.f90'
call Interpret_Program_Arguments(nmldeffile,1,(/ 12 /) )

! perform the zone axis computations
call Kosselpattern(nmldeffile)

end program CTEMKossel

!--------------------------------------------------------------------------
!
! SUBROUTINE:Kosselpattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief computation of a Kossel diffraction pattern
!
!> @param nmlfile namelist file name
!
!> @date 01/07/14  MDG 1.0 complete rewrite
!--------------------------------------------------------------------------
subroutine Kosselpattern(nmlfile)

use symmetryvars
use symmetry
use Lambert, ONLY: Apply2DPGSymmetry
use crystalvars
use crystal
use constants
use io
use error
use local
use files
use gvectors
use kvectors
use diffraction
use multibeams
use dynamical
use postscript
use timing
use MBmodule

IMPLICIT NONE

character(fnlen),INTENT(IN)	:: nmlfile

real(kind=sgl)      		:: ktmax, io_real(3), bragg, thetac, minten,  galen, delta, kstar(3), gperp(3), &
                       	   frac, dmin, convergence, voltage, startthick, thickinc, klaue(2), thetam 
integer(kind=irg)   		:: ijmax,ga(3),gb(3),k(3),cnt,fn(3), numthick, pgnum, &
                      		   newcount,count_rate,count_max, io_int(6), i, j, isym, skip, &
                      		   npx, npy, numt, numk, npix, ik, ip, jp, istat, dgn, &
                      		   maxHOLZ, numset, nn, gzero, ipx, ipy, ii, iequiv(2,12), nequiv
character(3)			:: method
character(fnlen)     		:: outname, xtalname

real(kind=sgl),allocatable    	:: thickarray(:)
real(kind=sgl),allocatable    	:: Iz(:), Izsum(:,:,:)
real(kind=dbl)			:: pre, tpi
complex(kind=dbl)              :: czero


namelist /Kossellist/ stdout, xtalname, voltage, k, fn, dmin, convergence, minten, &
                              startthick, thickinc, numthick, outname, npix, maxHOLZ

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		! initial value to check that the keyword is present in the nml file
stdout = 6			! standard output
voltage = 200000.0		! acceleration voltage [V]
k = (/ 0, 0, 1 /)		! beam direction [direction indices]
fn = (/ 0, 0, 1 /)		! foil normal [direction indices]
dmin = 0.025			! smallest d-spacing to include in dynamical matrix [nm]
convergence = 25.0		! beam convergence angle [mrad]
startthick = 10.0		! starting thickness [nm]
thickinc = 10.0			! thickness increment
numthick = 10			! number of increments
npix = 256			! output arrays will have size npix x npix
minten = 1.0E-5			! minimum intensity in diffraction disk to make it into the output file
outname = 'Kosselout.data'	! output filename

! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=Kossellist)
 close(UNIT=dataunit,STATUS='keep')

 if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMKossel:',' structure file name is undefined in '//nmlfile)
 end if

! print some information
 progname = 'CTEMKossel.f90'
 progdesc = 'Zone axis Kossel pattern simulation'
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
 
 ! for some of the 2D point groups, the standard orientation of the group according to ITC vol A
! may not be the orientation that we have here, so we need to determine by how much the 2D point
! group is rotated (CCW) with respect to the standard setting...
 call CheckPatternSymmetry(k,ga,isym,thetam)

! initialize the HOLZ geometry type
 call GetHOLZGeometry(float(ga),float(gb),k,fn) 

! construct the list of all possible reflections
 method = 'ALL'
 thetac = convergence/1000.0
 maxHOLZ = 3
 call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxHOLZ,thetac)
 galen = CalcLength(float(ga),'r')

! determine range of incident beam directions
  bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5
  
! convert to ktmax along ga
  ktmax = 0.5*thetac/bragg

! the number of pixels across the disk is equal to 2*npix + 1
  npx = npix
  npy = npx
  io_int(1) = 2.0*npx + 1
  call WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")
  mess=' '; call Message("(A/)")
  
! get number of thicknesses for which to compute the LACBED disk patterns
  numt = numthick
  allocate(thickarray(numt),stat=istat)
  thickarray = startthick + thickinc* (/ (float(i),i=0,numt-1) /)

! set parameters for wave vector computation
  klaue = (/ 0.0, 0.0 /)
  ijmax = float(npx)**2   ! truncation value for beam directions

! determine all independent incident beam directions (use a linked list starting at khead)
!  isym = WPPG(dgn)

! for now, the solution to the symmetry problem is to do the computation for the entire 
! illumination cone without application of symmetry.  Instead, we'll get the speed up by 
! going to multiple cores later on.
!  isym = 1
  call CalckvectorsSymmetry(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue,.FALSE.)

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(.TRUE.)


!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
! point to the first beam direction
  ktmp => khead
  czero = cmplx(0.D0,0.D0)
  pre = cmplx(0.D0,1.D0) * cPi
  call CalcUcg((/0,0,0/))   ! get the normal absorption parameter
  write (*,*) 'normal absorption length ',rlp%xgp
  DynUpz = rlp%Vpmod
  numset = cell % ATOM_ntype  ! number of special positions in the unit cell
  tpi = 2.D0*cPi
! allocate space for the results
  allocate(Iz(numthick),Izsum(2*npx+1,2*npy+1,numthick))
  
  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")
  gzero = 1
  
! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)

!  work through the beam direction list
  beamloop: do ik=1,numk

	ip = -ktmp%i
 	jp =  ktmp%j

! compute the dynamical matrix using Bloch waves with Bethe potentials; note that the IgnoreFoilNormal flag
! has been set to .FALSE.; if it is set to .TRUE., the computation of the ZOLZ will still be mostly correct,
! but the excitation errors of the HOLZ reflections will be increasingly incorrect with HOLZ order.  This was
! useful during program testing but should probably be removed as an option altogether...
	call Compute_DynMat('BLOCHBETHE', ktmp%k, ktmp%kt, .FALSE.)
        nn = DynNbeams

! solve the dynamical eigenvalue equation for this beam direction
        allocate(W(nn),CG(nn,nn),alpha(nn))
!write (*,*) ik,' calling CalcKint ',nn, numthick, ip, jp
        call CalcKint(nn,numthick,thickarray,Iz)
!write (*,*) ik,' returning from CalcKint '
        deallocate(W, CG, alpha)
        
        ipx = ktmp%i
        ipy = ktmp%j

! apply pattern symmetry 
       call Apply2DPGSymmetry(ipx,ipy,isym,iequiv,nequiv)

       do ii=1,nequiv
! is this point inside the viewing square ?
         ip = iequiv(1,ii) + npx + 1
         jp = iequiv(2,ii) + npy + 1
         if (((ip.ge.1).and.(ip.le.2*npix+1)).and.((jp.ge.1).and.(jp.le.2*npix+1))) then
           Izsum(ip,jp,1:numthick) = Iz(1:numthick)
         end if
       end do
      
! select next beam direction
        if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
	 if (float(ik)/float(numk) .gt. frac) then
	    io_int(1) = nint(100.0*frac) 
	    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
	    frac = frac + 0.05
	 end if  
	
    end do beamloop
    
! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F10.3)")
  
  
! store additional information for the IDL interface  
  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
! write the program identifier
  write (dataunit) trim(progname)
! write the version number
  write (dataunit) scversion
! first write the array dimensions
  write (dataunit) 2*npx+1,2*npy+1,numthick
! then the name of the crystal data file
  write (dataunit) xtalname
! the accelerating voltage [V]
  write (dataunit) voltage
! convergence angle [mrad]
  write (dataunit) convergence
! max kt value in units of ga
  write (dataunit) ktmax
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
! length horizontal reciprocal lattice vector (need for proper Laue center coordinate scaling)
  write (dataunit) galen
! we need to store the gperp vectors
  delta = 2.0*ktmax*galen/float(2*npix+1)        ! grid step size in nm-1 
  call TransSpace(float(k),kstar,'d','r')        ! transform incident direction to reciprocal space
  call CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
  call NormVec(gperp,'r')                        ! normalize g_perp
  write (dataunit) delta
  write (dataunit) gperp
! eight integers with the labels of various symmetry groups
  write (dataunit) (/ pgnum, PGLaue(pgnum), dgn, PDG(dgn), BFPG(dgn), WPPG(dgn), DFGN(dgn), DFSP(dgn) /)
! thickness data
  write (dataunit) startthick, thickinc
! and finaly the full pattern array
  write (dataunit) Izsum ! sr  
  close(UNIT=dataunit,STATUS='keep')
  
  mess = ' Data stored in '//trim(outname); call Message("(/A/)") 

end subroutine Kosselpattern

