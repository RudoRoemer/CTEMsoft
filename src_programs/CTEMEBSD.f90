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
! CTEMsoft2013:CTEMEBSD.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMEBSD
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief CTEMEBSD computes energy-weighted EBSD patterns
!
!> @date  08/01/12  MDG 1.0 EBSD extraction program for fundamental zone patterns
!> @date  08/17/12  MDG 1.1 generalized fundamental zone to other symmetries
!> @date  09/20/12  MDG 1.2 adapted for Lambert projection
!> @date  09/25/12  MDG 1.3 prepared for multithreaded version by separating computation steps
!> @date  12/11/12  MDG 2.0 new branch with energy-dependent Lambert projections (cubic only for now)
!> @date  02/26/14  MDG 3.0 incorporation into git and adapted to new libraries
!> @date  03/26/14  MDG 3.1 modification of file formats; made compatible with IDL visualization interface
! ###################################################################
! 

program CTEMEBSD

use local
use files
use io

IMPLICIT NONE

character(fnlen)	:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMEBSD.nml'
progname = 'CTEMEBSD.f90'
call Interpret_Program_Arguments(nmldeffile,1,(/ 21 /) )

! generate a set of EBSD patterns
 call ComputeEBSDPatterns(nmldeffile)

end program CTEMEBSD

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeEBSDPatterns
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute an energy-weighted EBSD pattern
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 08/01/13  MDG 3.0 complete rewrite, eliminated old Lambert projection
!> @date 09/25/13  MDG 3.1 replaced k-vector code by kvectors module
!> @date 02/26/14  MDG 4.0 new version
!> @date 03/26/14  MDG 4.1 adapted to new input and out file formats
!--------------------------------------------------------------------------
subroutine ComputeEBSDPatterns(nmlfile)


use local
use symmetryvars
use symmetry
use crystalvars
use crystal
use constants
use io
use files
use diffraction
use multibeams
use dynamical
use timing
use Lambert
use quaternions
use rotations
use noise

IMPLICIT NONE

character(fnlen),INTENT(IN)	:: nmlfile

! all geometrical parameters and filenames
real(kind=sgl)		:: L,L2,Ls,Lc	! distance between scintillator screen and interaction point [microns] (~working distance)
real(kind=sgl)		:: sig		! sample tile angle [degrees]
real(kind=sgl)		:: thetac	! detector tilt angle below horizontal [degrees]
real(kind=sgl)		:: delta	! scintillator step size [microns]
integer(kind=irg)	:: numsx 	! number of scintillator points along x
integer(kind=irg)	:: numsy	! number of scintillator points along y
real(kind=sgl)		:: xpc		! pattern center x [pixels]
real(kind=sgl)		:: ypc		! pattern center y [pixels]
character(fnlen)	:: eulerfile, FZfile, datafile, energyfile
real(kind=sgl)		:: energymin, energymax ! energy window for energy-filtered EBSD
integer(kind=irg)	:: numEbins, numzbins, nsx, nsy, nE, Emin, Emax, numelectrons  ! variables used in MC energy file
real(kind=dbl)		:: EkeV, Ehistmin, Ebinsize, depthmax, depthstep ! enery variables from MC program
real(kind=dbl)		:: beamcurrent, dwelltime, prefactor, MCsig, MComega
character(4)		:: MCmode	! Monte Carlo mode

! allocatable arrays
real(kind=sgl),allocatable		:: scin_x(:), scin_y(:) 		! scintillator coordinate ararays [microns]
real(kind=sgl),allocatable		:: rgx(:,:), rgy(:,:), rgz(:,:)  	! auxiliary arrays needed for interpolation
real(kind=sgl),allocatable		:: eulang(:,:)				! euler angle array
real(kind=sgl),allocatable		:: EBSDpattern(:,:)			! array with EBSD pattern
real(kind=sgl),allocatable 		:: sr(:,:,:), EkeVs(:)			! dynamical and kinematical parts of the FZ intensities
real(kind=sgl),allocatable		:: imagestack(:,:,:), z(:,:)		! used to store the computed patterns before writing to disk
integer(kind=irg),allocatable		:: accum_e(:,:,:)
real(kind=sgl),allocatable 		:: accum_e_detector(:,:,:)

! quaternion variables
real(kind=sgl),allocatable		:: quatang(:,:)
real(kind=sgl)				:: qq(4), qq1(4), qq2(4), qq3(4)

! various items
integer(kind=irg)	:: numeuler, numstacks	! number of Euler angle triplets in file
integer(kind=irg)	:: i, j, iang,k, io_int(6), num_el, MCnthreads, etotal		! various counters
integer(kind=irg)	:: istat		! status for allocate operations
integer(kind=irg)	:: nix, niy, npx, npy, offx, offy, imcnt	! various parameters
real(kind=sgl),parameter 	:: dtor = 0.0174533  ! convert from degrees to radians
real(kind=dbl),parameter	:: nAmpere = 6.241D+18   ! Coulomb per second
integer(kind=irg),parameter	:: storemax = 20	! number of EBSD patterns stored in one output block
real(kind=sgl)		:: alp, sa, ca		! angle and cosine and sine of alpha
real(kind=sgl)		:: dc(3), scl		! direction cosine array
real(kind=sgl)		:: sx, dx, dxm, dy, dym, rhos, io_real(6), x 	! various parameters
real(kind=sgl)		:: ixy(2)
character(3)		:: eulerconvention
character(5)		:: anglemode	! 'quats' or 'euler' for angular input
character(6)		:: sqorhe	! from Master file, square or hexagonal Lmabert projection
character(8)		:: Masterscversion, MCscversion
character(fnlen)	:: Masterprogname, Masterxtalname, Masterenergyfile, MCprogname, MCxtalname

! parameter for random number generator
integer, parameter 	:: K4B=selected_int_kind(9)      ! used by ran function in math.f90
integer(K4B) 		:: idum

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDdata / L, sig, thetac, delta, numsx, numsy, xpc, ypc, anglemode, eulerfile, eulerconvention, FZfile, &
			energyfile, datafile, beamcurrent, dwelltime, energymin, energymax

! spit out some information about the program 
progname = 'CTEMEBSD.f90'
progdesc = 'Dynamical EBSD patterns, using precomputed MC and master Lambert projections'
call CTEMsoft

! define reasonable default values for the namelist parameters
L		= 20000.0 	! [microns]
sig		= 70.0		! [degrees]
thetac		= 0.0		! [degrees]
delta		= 25.0		! [microns]
numsx		= 640		! [dimensionless]
numsy		= 480		! [dimensionless]
xpc		= 0.0		! [pixels]
ypc		= 0.0		! [pixels]
anglemode 	= 'euler'	! mode of the angular input data ('euler' or 'quats')
energymin	= 15.0		! minimum energy to consider
energymax	= 30.0		! maximum energy to consider
eulerfile	= 'euler.txt'	! filename
eulerconvention = 'tsl'	! convention for the first Euler angle ['tsl' or 'hkl']
FZfile		= 'FZ.data'	! filename
energyfile 	= 'energy.data' ! name of file that contains energy histograms for all scintillator pixels (output from MC program)
datafile	= 'EBSDout.data'	! output file name
beamcurrent 	= 14.513D-9	! beam current (actually emission current) in ampere
dwelltime 	= 100.0D-6	! in seconds
 
! then we read the rundata namelist, which may override some of these defaults  
! this could also be replaced (or duplicated) by the ability to directly call
! this program from IDL ?
 OPEN(UNIT=dataunit,FILE=trim(nmlfile),DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=EBSDdata)
 CLOSE(UNIT=dataunit)

! this routine should work as follows: the main routine sets up all relevant
! variables and geometry, so that the computation of a single EBSD pattern
! becomes just a function library call with the orientation as the only variable.
! That means that the complete detector geometry must be determined first; for
! each energy bin, we can then pre-compute the complete detector background.
! Computing the pattern is then only a matter of interpolating the master pattern
! with all its energy bins; then we combine the two arrays, collapse them in the 
! energy direction and apply the detector point spread function.  

! We should also include the ability to generate a dictionary; this means that 
! we should have the program sample the cubochoric space, then decide which
! points belong to the correct Rodrigues Fundamental Zone, and do the pattern
! computations for that complete set.  This may become part of the EBSD
! consulting project.



!====================================
! get the angular information, either in Euler angles or in quaternions, from a file
!====================================
if (anglemode.eq.'euler') then
! open the euler angle file and read the entries; there should be three reals on each line
! with the first line an integer indicating how many triplets there are in the file ... 
	open(unit=dataunit,file=trim(eulerfile),status='old',action='read')
	read(dataunit,*) numeuler  !MNS changed unit = dataunit to dataunit to be compatible with gfortran

	io_int(1) = numeuler
	call WriteValue('Number of euler angles = ',io_int,1)

	! allocate the euler angle array
	allocate(eulang(3,numeuler),stat=istat)
	! if istat.ne.0 then do some error handling ... 
	do i=1,numeuler
	  read(dataunit,*) eulang(1:3,i)
	end do
	close(unit=dataunit,status='keep')
	
	if (eulerconvention.eq.'hkl') then
	  mess = '  -> converting Euler angles to TSL representation'
	  call Message("(A/)")
	  eulang(1,1:numeuler) = eulang(1,1:numeuler) + 90.0
	end if
	
	! convert the euler angle triplets to quaternions
	allocate(quatang(4,numeuler),stat=istat)
	! if (istat.ne.0) then ...

	mess = '  -> converting Euler angles to quaternions'
	call Message("(A/)")

	do i=1,numeuler
	  quatang(1:4,i) = eu2qu(eulang(1:3,i)*dtor)
	end do
else
	open(unit=dataunit,file=trim(eulerfile),status='old',action='read')
	read(dataunit,*) numeuler
	io_int(1) = numeuler
	call WriteValue('Number of quaternions = ',io_int,1)
	allocate(quatang(4,numeuler),stat=istat)
	do i=1,numeuler
	  read(dataunit,*) quatang(1:4,i)
	end do
	close(unit=dataunit,status='keep')
end if
!====================================

!====================================
! for the creation of an EBSD dictionary, we'll need some special code
! here instead of the above reading in of Euler angles/quaternions; the
! program should create a uniform sampling in the appropriate fundamental
! zone (we should aim for Rodrigues) and then perform the dictionary
! computation for that representation and sampling.
!====================================


!====================================
! ------ read Monte Carlo data file
!====================================
! first, we need to load the data from the MC program.  This is an array of integers, which
! has an energy histogram for each sampled exit direction.  These are integer values to minimize
! storage, but they should all be divided by the total number of counts, which can be done 
! at the end of the computation.  We will then need to multiply by the beam current and by
! the dwell time to get units of electron counts.
!
! the datafile format is as follows (all in Lambert projections)
! write(dataunit) numEbins, numzbins, numsx, numsy, num_el, nthreads
! write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
! write(dataunit) accum_e
! write (dataunit) accum_z

mess = 'opening '//trim(energyfile)
call Message("(A)")

open(dataunit,file=trim(energyfile),status='unknown',form='unformatted')


! write the program identifier
 read (dataunit) MCprogname
! write the version number
 read (dataunit) MCscversion
! then the name of the crystal data file
 read (dataunit) MCxtalname
! energy information etc...



read(dataunit) numEbins, numzbins, nsx, nsy, num_el, MCnthreads
nsx = (nsx - 1)/2
nsy = (nsy - 1)/2

!io_int(1:6) = (/ numEbins, numzbins, nsx, nsy, num_el, MCnthreads /)
!call WriteValue(' NumEbins, numzbins, nsx, nsy, num_el, MCnthreads ',io_int,6,"(5I,',',I)")
etotal = num_el ! * MCnthreads

read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
!io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
!call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")

 read (dataunit) MCsig, MComega
 read (dataunit) MCmode

allocate(accum_e(numEbins,-nsx:nsx,-nsy:nsy),stat=istat)
read(dataunit) accum_e
num_el = sum(accum_e)

! we do not need the other array in this energyfile
! read(dataunit) accum_z    ! we only need this array for the depth integrations

close(dataunit,status='keep')

mess = ' -> completed reading '//trim(energyfile)
call Message("(A)")

! get the indices of the minimum and maximum energy
Emin = nint((energymin - Ehistmin)/Ebinsize) +1
if (Emin.lt.1)  Emin=1
if (Emin.gt.numEbins)  Emin=numEbins

Emax = nint((energymax - Ehistmin)/Ebinsize) +1
if (Emax.lt.1)  Emax=1
if (Emax.gt.numEbins)  Emax=numEbins
!====================================

write (*,*) 'shape(accum_e) = ',shape(accum_e)
write (*,*) 'energy range (Emin, Emax) = ', Emin, Emax, num_el

!====================================
! ----- Read energy-dispersed Lambert projections (master pattern)
! this has been updated on 3/26/14 to accommodate the new EBSDmaster file format
! but will need to be redone in HDF5 at a later time.
!====================================
open(unit=dataunit,file=FZfile,status='old',form='unformatted')
  read (dataunit) Masterprogname
! write the version number
  read (dataunit) Masterscversion
! then the name of the crystal data file
  read (dataunit) Masterxtalname
! then the name of the corresponding Monte Carlo data file
  read (dataunit) Masterenergyfile
! energy information and array size    
  read (dataunit) npx,npy,nE 
  if (numEbins.ne.nE) then
    write (*,*) 'Energy histogram and Lambert stack have different energy dimension; aborting program'
    write (*,*) 'energy histogram = ',shape(accum_e)
    write (*,*) 'Lambert stack = ', nE, npx, npy
    stop
  end if
  allocate(sr(-npx:npx,-npy:npy,nE),EkeVs(nE),stat=istat)
  read (dataunit) EkeVs
! is this a regular (square) or hexagonal projection ?
  read (dataunit) sqorhe
! and finally the results array
  read (dataunit) sr
close(unit=dataunit,status='keep')
mess = ' -> completed reading Master EBSD pattern file'
call Message("(A)")
!====================================


!====================================
! ------ generate the detector arrays
!====================================
! This needs to be done only once for a given detector geometry
allocate(scin_x(numsx),scin_y(numsy),stat=istat)
! if (istat.ne.0) then ...
scin_x = - ( xpc - ( 1.0 - numsx ) * 0.5 - (/ (i-1, i=1,numsx) /) ) * delta
scin_y = ( ypc - ( 1.0 - numsy ) * 0.5 - (/ (i-1, i=1,numsy) /) ) * delta

! auxiliary angle to rotate between reference frames
alp = 0.5 * cPi - (sig - thetac) * dtor
ca = cos(alp)
sa = sin(alp)

! we will need to incorporate a series of possible distortions 
! here as well, as described in Gert nolze's paper; for now we 
! just leave this place holder comment instead

! compute auxilliary interpolation arrays
allocate(rgx(numsx,numsy), rgy(numsx,numsy), rgz(numsx,numsy), stat=istat)
! if (istat.ne.0) then ...
L2 = L * L
do j=1,numsy
  sx = L2 + scin_y(j) * scin_y(j)
  Ls = ca * scin_y(j) + L*sa
  Lc = -sa * scin_y(j) + L*ca
  do i=1,numsx
   rhos = 1.0/sqrt(sx + scin_x(i)**2)
   rgx(i,j) = Ls * rhos
   rgy(i,j) = scin_x(i) * rhos
   rgz(i,j) = Lc * rhos
  end do
end do

! normalize the direction cosines.
allocate(z(numsx,numsy))
z = sqrt(rgx**2+rgy**2+rgz**2)
rgx = rgx/z
rgy = rgy/z
rgz = rgz/z
deallocate(z)
!====================================

!open(dataunit,file='test.data',status='unknown',form='unformatted')
!write (dataunit) rgx
!write (dataunit) rgy
!write (dataunit) rgz
!close(unit=dataunit,status='keep')

!====================================
! ------ create the equivalent detector energy array
!====================================
! from the Monte Carlo energy data, we need to extract the relevant
! entries for the detector geometry defined above.  Once that is 
! done, we can get rid of the larger energy arrays
!
! in the old version, we either computed the background model here, or 
! we would load a background pattern from file.  In this version, we are
! using the background that was computed by the MC program, and has 
! an energy histogram embedded in it, so we need to interpolate this 
! histogram to the pixels of the scintillator.  In other words, we need
! to initialize a new accum_e array for the detector by interpolating
! from the Lambert projection of the MC results.
!
  call InitLambertParameters
  allocate(accum_e_detector(numEbins,numsx,numsy), stat=istat)

write (*,*) 'shape(accum_e_detector) = ',shape(accum_e_detector), nsx

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
  scl = float(nsx) / LPs%sPio2

  do i=1,numsx
    do j=1,numsy
! do the coordinate transformation for this detector pixel
       dc = (/ rgx(i,j),rgy(i,j),rgz(i,j) /)
! make sure the third one is positive; if not, switch all 
       if (dc(3).lt.0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
	ixy = scl * LambertSphereToSquare( dc, istat )
	x = ixy(1)
	ixy(1) = ixy(2)
	ixy(2) = -x
! four-point interpolation (bi-quadratic)
        nix = int(nsx+ixy(1))-nsx
        niy = int(nsy+ixy(2))-nsy
        dx = ixy(1)-nix
        dy = ixy(2)-niy
        dxm = 1.0-dx
        dym = 1.0-dy
! interpolate the intensity 
        do k=Emin,Emax 
          accum_e_detector(k,i,j) =   accum_e(k,nix,niy) * dxm * dym + &
                          		accum_e(k,nix+1,niy) * dx * dym + &
					accum_e(k,nix,niy+1) * dxm * dy + &
                          		accum_e(k,nix+1,niy+1) * dx * dy
        end do
    end do
  end do 
  accum_e_detector = accum_e_detector * 0.25
! and finally, get rid of the original accum_e array which is no longer needed
  deallocate(accum_e)
  write (*,*) 'detector generation completed'
!====================================

num_el = nint(sum(accum_e_detector))

open(dataunit,file='test.data',status='unknown',form='unformatted')
write (dataunit) accum_e_detector
close(unit=dataunit,status='keep')

!

!====================================
! init a bunch of parameters
!====================================
! allocate the array that will hold the computed pattern
allocate(EBSDpattern(numsx,numsy),stat=istat)
! if (istat.ne.0) then ...
idum = -1		! to initialize the random number generator


! and allocate an array to store, say, 100 images; after every 100, all of them
! are stored in tiff format;  this is done to make sure that disk I/O doesn't
! become a bottleneck.
allocate(imagestack(numsx,numsy,storemax),stat=istat)
imcnt = 1

write (*,*) 'shape(imagestack) = ',shape(imagestack)
write (*,*) 'shape(Lambertstack) = ',shape(sr)

prefactor = 0.25D0 * nAmpere * beamcurrent * dwelltime / dble(num_el)
write (*,*) 'prefactor = ',prefactor

!====================================

!====================================
! ------ and open the output file for IDL visualization
!====================================
open(unit=dataunit,file=trim(datafile),status='unknown',form='unformatted',action='write')
! we need to write the imagestack dimensions, and also how many of those there are...
numstacks = numeuler/storemax + 1
write (dataunit) numsx, numsy, storemax, numstacks, numeuler

!====================================
! ------ start the actual image computation loop
!====================================
! start the timer
call Time_start

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
  scl = float(npx) / LPs%sPio2

do iang=1,numeuler
! convert the direction cosines to quaternions, include the 
! sample quaternion orientation, and then back to direction cosines...
! then convert these individually to the correct EBSD pattern location
        qq1 = conjg(quatang(1:4,iang))
        qq2 = quatang(1:4,iang)
        EBSDpattern = 0.0

	do i=1,numsx
	    do j=1,numsy
!  do the coordinate transformation for this euler agle
              qq = (/ 0.0, rgx(i,j),rgy(i,j),rgz(i,j) /)
              qq3 = quat_mult(qq2, quat_mult(qq,qq1) )
              dc(1:3) = (/ qq3(2), qq3(3), qq3(4) /) ! these are the direction cosines 
! make sure the third one is positive; if not, switch all 
              dc = dc/sqrt(sum(dc**2))
              if (dc(3).lt.0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
	      ixy = scl * LambertSphereToSquare( dc, istat )
! four-point interpolation (bi-quadratic)
              nix = int(npx+ixy(1))-npx
              niy = int(npy+ixy(2))-npy
              dx = ixy(1)-nix
              dy = ixy(2)-niy
              dxm = 1.0-dx
              dym = 1.0-dy
 ! interpolate the intensity 
              do k=Emin,Emax 
                EBSDpattern(i,j) = EBSDpattern(i,j) + accum_e_detector(k,i,j) * ( sr(nix,niy,k) * dxm * dym + &
                                           sr(nix+1,niy,k) * dx * dym + sr(nix,niy+1,k) * dxm * dy + &
                                           sr(nix+1,niy+1,k) * dx * dy )
              end do
          end do
       end do

! add sampling noise (Poisson noise in this case, so multiplicative, sort of)
!	do i=1,numsx
! 	  do j=1,numsy
!              EBSDpattern(i,j) = POIDEV(EBSDpattern(i,j),idum)
!	  end do
!	end do      

! we may need to deal with the energy sensitivity of the scintillator as well...


! all the following things should really be done in the IDL visualization program:
!
! - point spread function of camera
! - binning
! - brightness/contrast scaling
! - storage in image files of large HDF5 ? file
!  
! that means that at this point, we really only need to store all the patterns in a single
! file, at full resolution, as observed at the scintillator stage.
       
       

! next, apply the optical point spread function to this pattern
! [to optimize the speed, we may want to change the size of the array to a square shape]

! TO BE IMPLEMENTED


! and apply any camera binning
!	if (binning.ne.1) then 
!	  do i=1,numsx,binning
!	    do j=1,numsy,binning
!	        binned(i/binning+1,j/binning+1) = sum(EBSDpattern(i:i+binning-1,j:j+binning-1))
!	    end do
!	  end do  
!	else
!  	  binned = EBSDpattern
!	end if
!! and divide by binning^2
!	binned = binned * bindx
	

! this must also be removed; we're not going to apply contrast/brightness acling at this point;
! we'll allow the user to do this in the visualization code.
! and finally, before saving the patterns, apply the contrast function and scale between 0 and 255
!	EBSDmax = maxval(binned)
!	EBSDmin = minval(binned)
!	if (iang.eq.1) write (*,*) 'Limit Intensity Values : ',EBSDmin,EBSDmax
!	imagestack(1:binx,1:biny,imcnt) = int(255*(contrast * ( (binned - EBSDmin)/(EBSDmax-EBSDmin) - 1.0) + 1.0))

! the following is a "best fit" solution without any basis in physics and is used just as a place holder ... 
	imagestack(1:numsx,1:numsy,imcnt) = prefactor * EBSDpattern(1:numsx,1:numsy)
	imcnt = imcnt+1

! this will need to become an HDF5 formatted file with all the program output
! it should be readable in IDL as well as DREAM.3D.

	if ((imcnt.gt.storemax).or.(iang.eq.numeuler)) then 
	  write (*,*) 'storing EBSD patterns at angle ',iang,' of ',numeuler
	  write (dataunit) imagestack*sngl(prefactor)
	  imcnt = 1
	end if
! that's it for this pattern ... 
end do

write (dataunit) accum_e_detector

  close(unit=dataunit,status='keep')

 call Time_stop(numeuler)


end subroutine ComputeEBSDPatterns

