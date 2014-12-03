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
! CTEMsoft2013:CTEMlacbed.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMlacbed 
!
!> @author Marc De Graef, Carnegie Mellon University
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
!> @date 09/24/13 MDG 2.4 replaced older output routines with Bloch parameter file
!--------------------------------------------------------------------------
program CTEMlacbed

use local
use files
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
!> @author Marc De Graef, Carnegie Mellon University
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
integer(kind=irg)   		:: g(3),ira,dpcnt,ppi,ijmax,ga(3),gb(3),k(3),fcnt,ccnt,cnt,fn(3), PX, PY, numthick, nstrong, ig, &
                      		   ii, newcount,count_rate,count_max, io_int(6), i, j, isym, ir, order, nn, skip, refpick, &
                      		  iorder, npx, npy, numt, im, numk, npix, ik, ip, jp, maxholz, numpix, istat, selection(3)
character(1)        		:: ans
character(3)			:: method
character(fnlen)     		:: outname, tname, xtalname


complex(kind=sgl)   		:: czero
logical             		:: overlap,np,first
real(kind=sgl),parameter     	:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
real(kind=sgl),allocatable    	:: disk(:,:,:), thickarray(:,:), thick(:)
integer(kind=irg),allocatable 	:: diskoffset(:,:), IPIV(:)
real(kind=sgl),allocatable    	:: inten(:,:), nab(:,:)
complex(kind=dbl),allocatable	:: CGinv(:,:), Wloc(:), lCG(:,:), lalpha(:)
integer(kind=irg),allocatable	:: strongreflections(:,:), StrongIndexList(:)


namelist /inputlist/ stdout, xtalname, voltage, camlen, k, fn, dmin, convergence, outname

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		! initial value to check that the keyword is present in the nml file
voltage = 200000.0		! acceleration voltage [V]
camlen = 1000.0			! camera length [mm]
k = (/ 0, 0, 1 /)		! beam direction [direction indices]
fn = (/ 0, 0, 1 /)		! foil normal [direction indices]
dmin = 0.025			! smallest d-spacing to include in dynamical matrix [nm]
convergence = 20.0		! beam convergence angle [mrad]
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

! some of this is legacy code and will need to be replaced by something a little more intelligent...
! the legacy version was intended for producing PostScript output, so the output was scaled in
! terms of inches/millimeters.  For the new IDL visualization interface, this is no longer needed.

! change from here ---------------

! determine range of incident beam directions
  bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5

! call ReadValue(' Enter the beam convergence angle theta_c in mrad: ', io_real, 1)
  thetac = convergence/1000.0
  
! convert to ktmax along ga
  ktmax = 0.5*thetac/bragg

! compute number of pixels along diameter of central disk for given camera length
  RR = 300.0/25.4   ! dots per millimeter for 300 dots per inch
  npx = nint(RR*camlen*thetac)
  npy = npx
  io_int(1) = 2.0*npx
  call WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")

  ijmax = float(npx)**2   ! truncation value for beam directions

! to here ---------------
! if we change this, then the Calckvectors routine may also need to be updated.

! determine all independent incident beam directions (use a linked list starting at khead)
! we should be using the StandardConical method since we do want to make use of symmetry
! operations, but for now we'll use the Conical method (no symmetry considerations).  
  call Calckvectors(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,'StandardConical')

! construct the list of all possible reflections.
! this needs a little work because we really need to create a master list with ALL
! reflections that will be needed for at least one beam direction; furthermore, this 
! requires figuring out which reflections are strong or weak or absent.
  method = 'ALL'
  maxholz = 0	! not used
  call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxholz)


! allocate the disk variable which will hold the entire computed pattern
!  npix = int(3.0*300.0)   ! 3 inches wide at 300 dpi
!  allocate(disk(numt,npix,npix))
!  disk=0.0

  sc = cell%mLambda*camlen*RR
!  scmax = 1.5*300.0 + npx
!  PX = npix/2

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters()

! allocate the offset array
! to get this array, we need to do a mock initialization of the dynamical matrix in zone axis orientation
  call Compute_DynMat('BLOCHBETHE', khead%k,khead%kt, .TRUE.) ! PGC added khead.kt
  frac = 0.05


! next, we need to get a list of all potential reflections for the given range of 
! wave vectors; this requires the basic steps of the dynamical matrix computations,
! without actually computing those matrices.  This list is then written to the 
! main data file.
  mess = 'Pruning reflection list (this takes a while ...) '
  call Message("(A)")
  call Prune_ReflectionList(numk,nstrong)
  io_int(1) = nstrong
  call WriteValue('Number of strong beams overall : ', io_int, 1, '(I8)') ! PGC I->I8

! compute the offset parameters for all diffraction disks (for the zone-axis orientation !!!)
! first normalize the zone axis in cartesian components; this is the z-axis
  call TransSpace(float(k),c,'d','c')
  call NormVec(c,'c')

! then make ga the x-axis
  call TransSpace(float(ga),gx,'r','c')
  call NormVec(gx,'c')

! compute the cross product between k and gx; this is the y-axis
  call CalcCross(c,gx,gy,'c','c',0)

! get the list of strong reflections and for each one determine the projections onto gx and gy  
! the IDL visualiation program will figure out which reflections are HOLZ and ZOLZ.
  allocate(strongreflections(3,nstrong),nab(2,nstrong))
  strongreflections(1:3,1) = (/ 0, 0, 0 /)
  nab(1:2,1) = (/ 0.0, 0.0 /)
  rltmpa => reflist%next
  rltmpa => rltmpa%next ! skip first reflection
  reflectionloop: do ig=2,DynNbeamsLinked
    if (rltmpa%famnum.ne.0) then
  	strongreflections(1:3,rltmpa%famnum) = rltmpa%hkl
        gg(1:3)=rltmpa%hkl
        call TransSpace(gg,c,'r','c')
        nab(1:2,rltmpa%famnum) = (/ CalcDot(c,gx,'c')*sc, CalcDot(c,gy,'c')*sc /)
    endif
! go to the next beam in the list
    rltmpa => rltmpa%next
  end do reflectionloop

! and write these arrays to the data file along with some other parameters
  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
! xtal file name
  write (dataunit) xtalname
! zone axis direction
  write (dataunit) k
! accelerating voltage
  write (dataunit) dble(voltage)
! number of incident beam directions
  write (dataunit) numk
! horizontal beam step number
  write (dataunit) npx
! max beam tilt angle [mrad]
  write (dataunit) convergence
! total number of strong reflections
  write (dataunit) nstrong
! Miller indices of strong reflections
  write (dataunit) strongreflections
! strong reflection projection components
  write (dataunit) nab
! keep this file open until the end of the program


  io_int(1)=numk
  call WriteValue(' Starting computation for # beam directions = ', io_int, 1, "(I8)")

! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)

! point to the first beam direction
  ktmp => khead
! loop over all beam orientations, selecting them from the linked list
  do ik = 1,numk
! compute the dynamical matrix using Bloch waves with Bethe potentials 
        call Compute_DynMat('BLOCHBETHE', ktmp%k,ktmp%kt, .TRUE.) ! PGC added ktmp%kt

! allocate a bunch of arrays for the Bloch wave eigenvalue computation
        allocate(CGinv(DynNbeams,DynNbeams), Wloc(DynNbeams), lCG(DynNbeams,DynNbeams), IPIV(DynNbeams))

! compute the eigenvalues and eigenvectors
 	IPIV = 0
 	call BWsolve(DynMat,Wloc,lCG,CGinv,DynNbeams,IPIV)

! and write the output to the data file
! write the number of strong reflections, and the ip,jp beam coordinate (these should really be k_t components)
	write (dataunit) BetheParameter%nns, ktmp%i, ktmp%j
! write the list of strong reflection IDs
	write (dataunit) BetheParameter%strongID
! write the scaled Bloch wave eigenvalue vector
	write (dataunit) Wloc/cmplx(ktmp%kn/cPi,0.0)
! write the Bloch wave eigenvector matrix
	write (dataunit) lCG
! write the Bloch wave excitation amplitude vector	
	write (dataunit) CGinv(1:DynNbeams,1)

       deallocate(CGinv, Wloc, lCG, IPIV, BetheParameter%strongID)

! intensity calculations must be done as follows for a thickness th
! diag=exp(-th*imag(lW))*cmplx(cos(th*real(lW)),sin(th*real(lW)))*lalpha   ! lalpha = CGinv(1:DynNbeams,1)
! inten = cabs(matmul(lCG,diag))**2



! first determine which reflection we are going to plot the pattern for.
! normally, it's the bright field pattern, so refpick = 1, but if selection
! is not (/ 0,0,0 /), then we need to first figure out what the number
! is of the selected reflection in the main list, and then we need to see if
! that reflection was actually used in the Bloch wave computation...
!	if (sum(abs(selection)).eq.0) then
!	  refpick = 1
!	else
!	  refpick = -1
!	  do ii = 1, BetheParameter%nns
!	    if (sum(abs(selection(1:3)-BetheParameter%stronghkl(1:3,ii))).eq.0) then
!	     refpick = ii
!	    end if
!	  end do
!!	  if (refpick.gt.0) write (*,*) ik, BetheParameter%stronghkl(1:3,refpick), refpick, &
!!	  BetheParameter%reflistindex(refpick)
!	end if


!! copy in the correct locations
!      if (refpick.eq.-1) then
!       disk(1:numt,ip,jp) = 0.0
!     else
!!      ii = BetheParameter%reflistindex(refpick)
!      do j=1,numt
!         disk(j,ip,jp) = disk(j,ip,jp) + inten(j,refpick)
!      end do
!     end if
!! and remove the intensity array
!     deallocate(inten)
!    
!   end if

! select next beam direction
   if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

  end do

! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F10.5)")

! close the output file
  close(UNIT=dataunit,STATUS='keep')  
  mess = ' Data stored in '//outname; call Message("(/A/)") 
 
end subroutine LACBEDpattern

