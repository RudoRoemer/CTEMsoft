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
!> @date 10/04/13  MDG 3.0 adaptation for new symmetry routines and output
!> @date 10/05/13  MDG 3.1 added output stuff for IDL visualization program
!> @date 10/07/13  MDG 3.2 corrected subtle error in reflection numbering
!> @date 10.08/13  MDG 3.3 added maxHOLZ output limitation; disk offset computation
!--------------------------------------------------------------------------
subroutine LACBEDpattern(nmlfile)

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

real(kind=sgl)      		:: ktmax, io_real(3), bragg, thetac, c(3), gx(3), gy(3), sc, qx, qy, &
                       	   frac, dmin, convergence, voltage, startthick, thickinc, klaue(2)
integer(kind=irg)   		:: ijmax,ga(3),gb(3),k(3),cnt,fn(3), PX, numthick, ss, icnt, pgnum, ih, &
                      		   newcount,count_rate,count_max, io_int(6), i, j, isym, ir, skip, &
                      		   npx, npy, numt, numk, npix, ik, ip, jp, maxholz, istat, dgn, nbeams, &
                      		   ifamily, isum, famhkl(3), inum, maxHOLZ
character(3)			:: method
character(fnlen)     		:: outname, xtalname

real(kind=sgl),allocatable    	:: diskoffset(:,:), disk(:,:,:,:), thick(:), familytwotheta(:), slice(:,:,:)
integer(kind=irg),allocatable 	:: familymult(:), familyhkl(:,:), whichHOLZ(:)
real(kind=sgl),allocatable    	:: inten(:,:)


namelist /inputlist/ stdout, xtalname, voltage, k, fn, dmin, convergence, &
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
maxHOLZ = 2			! maximum HOLZ layer index to be used for the output file; note that his number
				! does not affect the actual computations; it only determines which reflection 
				! families will end up in the output file
outname = 'lacbedout.data'	! output filename

camlen = 1000.0			! camera length [mm]   (this is not part of the namelist, but needed later)

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
 write (stdout,NML=inputlist)

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

! construct the list of all possible reflections
 method = 'ALL'
 maxholz = 0
 thetac = convergence/1000.0
 call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxholz,thetac)

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
  allocate(thick(numt),stat=istat)
  thick = startthick + thickinc* (/ (float(i),i=0,numt-1) /)

! set parameters for wave vector computation
  klaue = (/ 0.0, 0.0 /)
  ijmax = float(npx)**2   ! truncation value for beam directions

! determine all independent incident beam directions (use a linked list starting at khead)
!  isym = WPPG(dgn)

! for now, the solution to the symmetry problem is to do the computation for the entire 
! illumination cone without application of symmetry.  Instead, we'll get the spped up by 
! going to multiple cores later on.
  isym = 1
  call CalckvectorsSymmetry(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue,.TRUE.)

! set scaling parameters
  PX = npix/2

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(.TRUE.)

! up to this point, everything is nearly identical to the mbcbed program,
! except that we do not use a camera length explicitly.  Now we
! need to do things a little differently.  First of all, we need a master 
! list for all the reflections that contribute in one way or another to the
! complete LACBED pattern.  We do this by going through the entire list of incident
! beam directions and flagging all reflections that contribute, either as weak
! or as strong reflections.
  mess = ' Pruning reflection list (this takes a while ...) '
  call Message("(A)")
  call Prune_ReflectionList(numk,nbeams)
  io_int(1) = nbeams
  call WriteValue('Number of contributing beams  : ', io_int, 1, '(I)')

! since we're only going to store one diffraction disk per family of reflections,
! we need to decide which one we're going to keep.  For those, we'll need to 
! determine the 2D multiplicity and store that along with the Miller indices 
! of a representative family member.

! first normalize the zone axis in cartesian components; this is the z-axis
  call TransSpace(float(k),c,'d','c')
  call NormVec(c,'c')

! then make ga the x-axis
  call TransSpace(float(ga),gx,'r','c')
  call NormVec(gx,'c')

! compute the cross product between k and gx; this is the y-axis
  call CalcCross(c,gx,gy,'c','c',0)

! set the scale parameter for a default camera length of 1000 mm.
  sc = mLambda * 1000.0 * 300.0 / 25.4  ! the absolute value does not matter and is derived from legacy Postscript code
! The original code used 300 dpi (hence 300/25.4) which was convenient for Postscript output; in the current case, we
! do not actually use the true value, but in the IDL visualization program, we scale the user defined camera length by
! 1000.0, and use this ratio to scale the diskoffset coordinates.  So, there's no absolute length scale, only a relative scale.

! we do not know how many independent families there are yet, so we go through the 
! list once and then we'll allocate the proper arrays.
  ifamily = 1	! for the incident beam
  rltmpa => reflist%next
  rltmpa => rltmpa%next

outerloop: do while (associated(rltmpa%next))
    famhkl = rltmpa%famhkl
    ifamily = ifamily+1
    isum = 1
    rltmpa => rltmpa%next
    do while (sum(abs(rltmpa%famhkl-famhkl)).eq.0)
      rltmpa => rltmpa%next
      isum = isum+1
      if (.not.associated(rltmpa%next)) EXIT outerloop
    end do
end do outerloop

! ok, so there are ifamily families; next we need to store the corresponding
! hkl, and multiplicity, as well as the diffraction angle and the position of 
! the diffraction disk center for a standard camera length.
  allocate(familyhkl(3,ifamily), familymult(ifamily), familytwotheta(ifamily), diskoffset(2,ifamily))

! redo the above loop, but now fill in the data
  ifamily = 1	! for the incident beam
  familyhkl(1:3,ifamily) = (/ 0, 0, 0 /)
  familymult(ifamily) = 1
  diskoffset(1:2,ifamily) = (/ 0.0, 0.0 /)
  rltmpa => reflist%next
  rltmpa => rltmpa%next

outerloop2: do while (associated(rltmpa%next))
    famhkl = rltmpa%famhkl
    ifamily = ifamily+1
    familyhkl(1:3,ifamily) = famhkl(1:3)
    familytwotheta(ifamily) = CalcDiffAngle(famhkl(1),famhkl(2),famhkl(3))*1000.0
    familymult(ifamily) = 1
! get the disk offset parameters
    call TransSpace(float(famhkl),c,'r','c')
    qx = CalcDot(c,gx,'c')*sc
    qy = CalcDot(c,gy,'c')*sc
    diskoffset(1:2,ifamily) = (/ qx, qy /)
! and move on to the next one
    rltmpa => rltmpa%next
    do while (sum(abs(rltmpa%famhkl-famhkl)).eq.0)
      rltmpa => rltmpa%next
      familymult(ifamily) = familymult(ifamily) + 1
      if (.not.associated(rltmpa%next)) EXIT outerloop2
    end do
! print results for debugging
!write (*,*) ifamily, familyhkl(1,ifamily),familyhkl(2,ifamily),familyhkl(3,ifamily),familymult(ifamily),&
!familytwotheta(ifamily)
  end do outerloop2

! correct last counter
  familymult(ifamily) = familymult(ifamily)+1

! and print last entry (used for debugging)
!write (*,*) ifamily, familyhkl(1,ifamily),familyhkl(2,ifamily),familyhkl(3,ifamily),familymult(ifamily),&
!familytwotheta(ifamily)

  io_int(1) = ifamily
  call WriteValue('Number of unique families in output = ', io_int, 1, "(I5)")


! next we create the output array, which has one disk image for each 
! thickness and contributing family.  We make sure that each image is fully
! filled by a diffraction disk (i.e., no empty space along the main axes.
  allocate(disk(-npix:npix,-npix:npix,1:numt,1:ifamily))
  disk=0.0

  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")

! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)


! point to the first beam direction
  ktmp => khead
! loop over all beam orientations, selecting them from the linked list
kvectorloop:  do ik = 1,numk

	ip = -ktmp%i
 	jp =  ktmp%j

! compute the dynamical matrix using Bloch waves with Bethe potentials 
	call Compute_DynMat('BLOCHBETHE', ktmp%k, .TRUE.)

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

! ok, we have all the intensities.  Next we need to copy the relevant intensities into the slots 
! of the disk array, one for each family.
	rltmpa => reflist%next%next
	inum = 1
        disk(ip,jp,1:numt,inum) = inten(1:numt,1)
	do i=2,DynNbeamsLinked
	  if (BetheParameter%reflistindex(i).ne.0) then ! is this a reflection on the current list
! it is, so we need to determine which of the families corresponds to it
	    inum = -1
	    do ir=2,ifamily
	     ss = sum(abs(familyhkl(1:3,ir) - rltmpa%hkl(1:3)))
	     if (ss.eq.0) inum = ir
	    end do
  	    if (inum.ne.-1) then	
              disk(ip,jp,1:numt,inum) = inten(1:numt,BetheParameter%reflistindex(i))
	    end if
	  end if
	  rltmpa => rltmpa%next
	end do

! and remove the intensity array
     deallocate(inten)
    
! select next beam direction
   if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

  end do kvectorloop

! the following comment only applies if we use symmetry to determine the computational wedge.
! 
! next we need to apply the symmetry operators to ALL beams, including the special and 
! general diffraction symmetry for dark field disks... this is complicated, since the 
! reflections are not grouped by symmetrically equivalent classes.  There are two
! solutions; either we figure out on the spot which reflections are equivalent, which
! requires running through the entire reflection list several times, or we change the 
! ComputeReflections routine to list reflections consecutively by family.  Either way
! will work, but the latter one may result in a faster algorithm.  Whichever way we
! do this, we will need to make a copy of the entire disk array, so that we don't 
! apply the operators too many times...  To do so, we need to make sure that if one
! single family member has an entry in the disks array, then all equivalent family 
! members must also be in the array... (equivalent with respect to the zone axis WP symmetry)
! This will require a bit of thinking before implementation can begin.

! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F)")

! before we write the output file, we need to determine which reflection families 
! need to be written to the file; this is determined by the maxHOLZ parameter.  So,
! first we determine to which HOLZ layer each family belongs by using the zone 
! equation.   [check special case of hexagonal indices !!!!]
! also count the ones up to order maxHOLZ
  allocate(whichHOLZ(ifamily))
  icnt = 0
  do ir=1,ifamily
    whichHOLZ(ir) = iabs(k(1)*familyhkl(1,ir)+k(2)*familyhkl(2,ir)+k(3)*familyhkl(3,ir))
    if (whichHOLZ(ir).le.maxHOLZ) icnt = icnt+1
  end do  

! the final bit of the program involves dumping all the results into a file,
! binary for now, but HDF5 in the future, for the IDL visualization program 
! to read.
  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
! write the program identifier
  write (dataunit) trim(progname)
! write the version number
  write (dataunit) scversion
! first write the array dimensions
  write (dataunit) 2*npix+1,2*npix+1,numt,icnt
! then the name of the crystal data file
  write (dataunit) xtalname
! the accelerating voltage [V]
  write (dataunit) voltage
! convergence angle [mrad]
  write (dataunit) convergence
! the zone axis indices
  write (dataunit) k
! number of k-values in disk
  write (dataunit) numk
! horizontal reciprocal lattice vector
  write (dataunit) ga  
! maximum HOLZ layer in the output file
  write (dataunit) maxHOLZ
! eight integers with the labels of various symmetry groups
  write (dataunit) (/ pgnum, PGLaue(pgnum), dgn, PDG(dgn), BFPG(dgn), WPPG(dgn), DFGN(dgn), DFSP(dgn) /)
! thickness data
  write (dataunit) startthick, thickinc
! and from here one we write the individual diffraction disks with associated information  
! Miller indices, multiplicity, two-theta [mrad], and position for a reference camera length
! only for those families that belong the HOLZ layers <= maxHOLZ (to keep the file size down a bit)
  allocate(slice(-npix:npix,-npix:npix,1:numt))
  write (dataunit) familyhkl(1:3,1), familymult(1), familytwotheta(1), diskoffset(1:2,1), whichHOLZ(1)
  slice = disk(-npix:npix,-npix:npix,1:numt,1)
  write (dataunit) slice
! we'll write them in reverse order, so that the smaller Miller indices come first (at least for high symmetry structures);
! we'll also write them by HOLZ number
  do ih = 0,maxHOLZ
   do ir = ifamily,2,-1
    if (whichHOLZ(ir).eq.ih) then
      write (dataunit) familyhkl(1:3,ir), familymult(ir), familytwotheta(ir), diskoffset(1:2,ir), whichHOLZ(ir)
      slice = disk(-npix:npix,-npix:npix,1:numt,ir)
      write (dataunit) slice
    end if
   end do
  end do
  close(UNIT=dataunit,STATUS='keep')

  mess = ' Data stored in '//outname; call Message("(/A/)") 
  io_int(1) = maxHOLZ
  call WriteValue('Dats includes families of reflections in HOLZ layers 0 through ',io_int,1,"(I2)")
 
end subroutine LACBEDpattern




