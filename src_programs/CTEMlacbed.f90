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
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute a large angle zone axis convergent beam electron diffraction pattern
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 10/04/13  MDG 3.0 adaptation for new symmetry routines and output
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
                       	   sc, scmax, qx, qy, frac, dmin, convergence, voltage, startthick, thickinc, klaue(2)
integer(kind=irg)   		:: g(3),ira,dpcnt,ppi,ijmax,ga(3),gb(3),k(3),fcnt,ccnt,cnt,fn(3), PX, PY, numthick,  &
                      		   ii, newcount,count_rate,count_max, io_int(6), i, j, isym, ir, order, nn, skip, refpick, &
                      		  iorder, npx, npy, numt, im, numk, npix, ik, ip, jp, maxholz, numpix, istat, dgn, nbeams
character(1)        		:: ans
character(3)			:: method
character(fnlen)     		:: outname, tname, xtalname


complex(kind=sgl)   		:: czero
logical             		:: overlap,np,first
real(kind=sgl),parameter     	:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
real(kind=sgl),allocatable    	:: disk(:,:,:,:), thickarray(:,:), thick(:)
integer(kind=irg),allocatable 	:: diskoffset(:,:)
real(kind=sgl),allocatable    	:: inten(:,:)


namelist /inputlist/ stdout, xtalname, voltage, camlen, k, fn, dmin, convergence, &
                              startthick, thickinc, numthick, outname, npix

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		! initial value to check that the keyword is present in the nml file
stdout = 6			! standard output
voltage = 200000.0		! acceleration voltage [V]
camlen = 1000.0			! camera length [mm]
k = (/ 0, 0, 1 /)		! beam direction [direction indices]
fn = (/ 0, 0, 1 /)		! foil normal [direction indices]
dmin = 0.025			! smallest d-spacing to include in dynamical matrix [nm]
convergence = 20.0		! beam convergence angle [mrad]
startthick = 10.0		! starting thickness [nm]
thickinc = 10.0			! thickness increment
numthick = 10			! number of increments
npix = 512			! output arrays will have size npix x npix
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

write (*,*) 'ktmax = ',ktmax

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
! illumination cone without application of symmetry.  

  isym = 1
  call CalckvectorsSymmetry(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue,.TRUE.)

! set scaling parameters
  PX = npix/2

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters()

! up to this point, everything is nearly identical to the mbcbed program,
! except that we do not use a camera length explicitly.  Now we
! need to do things a little differently.  First of all, we need a master 
! list for all the reflections that contribute in one way or another to the
! complete LACBED pattern.  We do this by going through the entire list of incident
! beam diretions and flagging all reflections that contribute, either as weak
! or as strong reflections.
  mess = 'Pruning reflection list (this takes a while ...) '
  call Message("(A)")
  call Prune_ReflectionList(numk,nbeams)
  io_int(1) = nbeams
  call WriteValue('Number of contributing beams  : ', io_int, 1, '(I)')

! next we create the output array, which has one disk image for each 
! thickness and contributing beam.  We make sure that each image is fully
! filled by a diffraction disk (i.e., no empty space along the main axes.
  allocate(disk(-npix:npix,-npix:npix,1:numt,1:nbeams))
  disk=0.0

! allocate the offset array
! to get this array, we need to do a mock initialization of the dynamical matrix in zone axis orientation
!  call Compute_DynMat('BLOCHBETHE', khead%k, .TRUE.)

  io_int(1)=numk
  call WriteValue(' Starting computation for # beam directions = ', io_int, 1, "(I8)")

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

! ok, we have all the intensities.  now we need to copy them into the correct diffraction disks
! for now, we'll do this without application of the symmetry operators.  That works reasonably well,
! at least for non-three or sixfold symmetry, which we'll need to fix later.
	rltmpa => reflist%next
	do i=1,DynNbeamsLinked
	  if (BetheParameter%reflistindex(i).ne.0) then ! is this a reflection on the current list
! it is, so we need to determine which of the nbeams beams corresponds to it
	    ii = rltmpa%famnum
            disk(ip,jp,1:numt,ii) = inten(1:numt,BetheParameter%reflistindex(i))
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

! next we need to apply the symmetry operators to ALL beams, including the special and 
! general diffraction symmetry for dark field disks... this is complicated, since the 
! reflections are not grouped by symmetrically equivalent classes.  There are two
! solutions; either we figure out on the spot which reflections are equivalent, which
! requires running through the entire reflection list several times, or we change the 
! ComputeReflections routine to list reflections consecutively by family.  Either way
! will work, but the latter one may result in a faster algorithm.  Whichever way we
! do this, we will need to make a copy of the entire disk array, so that we don't 
! apply the operators too many times...  Either way, this routine will need to know
! ALL the details about the 31 diffraction groups.




!write (*,*) ' maximum intensity ', maxval(disk)
!write (*,*) ' minimum # strong beams ',BetheParameter%minstrong
!write (*,*) ' maximum # strong beams ',BetheParameter%maxstrong

! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F)")

! to avoid having more than 2Gb in one written section of the output file (IDL can not deal
! with sections that are larger) we will chop the array into smaller sections and write each
! one separately.  This will need to be recoded below.  In a later stage, we will replace
! all this with HDF5 output, so the problem will disappear.
  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
  write (dataunit) 2*npix+1,2*npix+1,numt,nbeams
  write (dataunit) disk
  close(UNIT=dataunit,STATUS='keep')


  
  mess = ' Data stored in '//outname; call Message("(/A/)") 
 
end subroutine LACBEDpattern




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
