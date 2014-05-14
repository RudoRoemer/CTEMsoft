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
use quaternions
use rotations
use tiff_global
use tiff_f90
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
real(kind=sgl)		:: xscale, yscale  ! scale factors (only needed for current stereographic coordinates)
real(kind=sgl)		:: sigx, sigy 	! forward scattering Gaussian model FWHMs
real(kind=sgl)		:: posx, posy	! Gaussian peak position components
real(kind=sgl)		:: contrast	! contrast factor for camera, 1 = max contrast, 0 = no contrast
character(fnlen)	:: eulerfile, FZfile, prefix, energyfile
real(kind=sgl)		:: energymin, energymax ! energy window for energy-filtered EBSD
integer(kind=irg)	:: numEbins, numzbins, nsx, nsy, nE, Emin, Emax, numelectrons  ! variables used in MC energy file
integer(kind=irg)	:: num_el, etotal, MCnthreads !PGC
real(kind=dbl)		:: EkeV, Ehistmin, Ebinsize, depthmax, depthstep ! enery variables from MC program
real(kind=dbl)		:: beamcurrent, dwelltime, prefactor
character(fnlen)	:: outputmode !PGC

! allocatable arrays
real(kind=sgl),allocatable		:: scin_x(:), scin_y(:) 		! scintillator coordinate ararays [microns]
real(kind=sgl),allocatable		:: rgx(:,:), rgy(:,:), rgz(:,:)  ! auxiliary arrays needed for interpolation
real(kind=sgl),allocatable		:: eulang(:,:)			! euler angle array
real(kind=sgl),allocatable		:: EBSDpattern(:,:),binned(:,:),bgprofile(:,:)		! array with EBSD pattern
real(kind=sgl),allocatable 		:: sr(:,:,:), srcopy(:,:)		! dynamical and kinematical parts of the FZ intensities
integer(kind=irg),allocatable		:: imagestack(:,:,:)		! used to store the computed patterns before writing to disk
integer(kind=irg),allocatable		:: accum_e(:,:,:)

! quaternion variables
real(kind=sgl),allocatable		:: quatang(:,:)
real(kind=sgl)				:: qq(4), qq1(4), qq2(4), qq3(4)

! various items
integer(kind=irg)	:: numeuler	! number of Euler angle triplets in file
integer(kind=irg)	:: i, j, iang,k		! various counters
integer(kind=irg)	:: istat		! status for allocate operations
integer(kind=irg)	:: nix, niy, npx, npy, offx, offy, binx, biny, imcnt, storemax,firstfile, saveid		! various parameters
integer(kind=irg)	:: binning		! camera binning factor ([1,2,4,8] are supported)
real(kind=sgl),parameter 	:: dtor = 0.0174533  ! convert from degrees to radians
real(kind=dbl),parameter	:: nAmpere = 6.241D+18   ! Coulomb per second
real(kind=sgl)		:: alp, sa, ca		! angle and cosine and sine of alpha
real(kind=sgl)		:: dc(3)		! direction cosine array
real(kind=sgl)		:: ix, iy, sx, dx, dxm, dy, dym, dd, EBSDmax, EBSDmin, bindx, rhos, rx, ry	! various parameters
character(4)		:: filenumber	! used to number files (obviously)
character(3)		:: eulerconvention
character(5)		:: anglemode	! 'quats' or 'euler' for angular input
character(fnlen)	:: LambertPattern	! string for filename of Lambert projection (tiff format)
character(fnlen)	:: BackgroundPattern	! string for filename of background intensity pattern
character(12)		:: mapmode 	! 'PlainLambert' or 'RoscaLambert'  (mode of the input map)
real(kind=irg)		:: io_int(6)  !PGC added
real(kind=sgl)		:: io_real(5)  !PGC added

! parameter for random number generator
integer, parameter 	:: K4B=selected_int_kind(9)      ! used by ran function in math.f90
integer(K4B) 		:: idum

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDdata / L, sig, thetac, delta, numsx, numsy, xpc, ypc, eulerfile, FZfile, prefix, binning, &
                       storemax, eulerconvention, LambertPattern, &
			BackgroundPattern, anglemode, outputmode, saveid, mapmode, energymin, energymax, &
			energyfile, beamcurrent, dwelltime, numelectrons
 
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
binning		= 8		! [dimensionless]
storemax	= 100		! every storemax patterns, save them to disk
outputmode 	= 'tiff'	! output patterns in individual tiff files or single 'data' file
anglemode 	= 'euler'	! mode of the angular input data ('euler' or 'quats')
eulerfile	= 'euler.txt'	! filename
eulerconvention = 'tsl'	! convention for the first Euler angle ['tsl' or 'hkl']
FZfile		= 'FZ.data'	! filename
prefix		= 'EBSD_'	! prefix for tiff output filenames
LambertPattern = 'none'	! don't store a Lambert projectionl if different from 'none', then this is the filename
BackgroundPattern = 'none'	! don't load a preexisiting background profile file
saveid		= -1		! -1 if no single image intermediate data needs to be stored, positive with image number
energymin 	= 0.0		! minimum energy that contributes to the pattern(s)
energymax 	= 30.0		! maximum energy (both in keV)
energyfile 	= 'energy.data' ! name of file that contains energy histograms for all scintillator pixels (output from MC program)
beamcurrent 	= 14.513D-9	! beam current (actually emission current) in ampere
dwelltime 	= 100.0D-6	! in seconds
numelectrons 	= 1.0D9		! number of incident electrons in MonteCarlo simulation (to be removed)
 
! then we read the rundata namelist, which may override some of these defaults  
 OPEN(UNIT=dataunit,FILE='CTEMEBSD.nml',DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=EBSDdata)
 CLOSE(UNIT=dataunit)
! if (contrast.lt.0.01) contrast = 0.01
! if (contrast.gt.1.0) contrast = 1.0

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
	read(dataunit,*) numeuler ! PGC unit=dataunit -> dataunit

	io_int(1) = numeuler
	call WriteValue('Number of euler angles = ',io_int,1)

	! allocate the euler angle array
	allocate(eulang(3,numeuler),stat=istat)
	! if istat.ne.0 then do some error handling ... 
	do i=1,numeuler
	  read(dataunit,*) eulang(1:3,i)! PGC unit=dataunit -> dataunit
	end do
	close(unit=dataunit,status='keep')
	
	mess = '  -> converting Euler angles to TSL representation'
	call Message("(A/)")
	if (eulerconvention.eq.'hkl') eulang(1,1:numeuler) = eulang(1,1:numeuler) + 90.0
	
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
	read(dataunit,*) numeuler!PGC unit=dataunit -> dataunit
	io_int(1) = numeuler
	call WriteValue('Number of quaternions = ',io_int,1)
	allocate(quatang(4,numeuler),stat=istat)
	do i=1,numeuler
	  read(dataunit,*) quatang(1:4,i)!PGC unit=dataunit -> dataunit
	end do
	close(unit=dataunit,status='keep')
end if
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

read(dataunit) numEbins, numzbins, nsx, nsy, num_el, MCnthreads
nsx = (nsx - 1)/2
nsy = (nsy - 1)/2

io_int(1:6) = (/ numEbins, numzbins, nsx, nsy, num_el, MCnthreads /)
call WriteValue(' NumEbins, numzbins, nsx, nsy, num_el, MCnthreads ',io_int,6,"(5I,',',I)")
etotal = num_el * MCnthreads

read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")


allocate(accum_e(numEbins,-nsx:nsx,-nsy:nsy),stat=istat)
read(dataunit) accum_e

! we do not need the other array in this energyfile
! read(dataunit) accum_z    ! we only need this array for the depth integrations

close(dataunit,status='keep')

mess = ' -> completed reading '//trim(energyfile)
call Message("(A)")
!====================================










! get the indices of the minimum and maximum energy
Emin = nint((energymin - Ehistmin)/Ebinsize) +1
if (Emin.lt.1)  Emin=1
if (Emin.gt.numEbins)  Emin=numEbins

Emax = nint((energymax - Ehistmin)/Ebinsize) +1
if (Emax.lt.1)  Emax=1
if (Emax.gt.numEbins)  Emax=numEbins

write (*,*) 'energy range from bin ',Emin, ' to bin ',Emax



!write (*,*) 'Performing sleep test fo five seconds'
!call system('/bin/sleep 5')
!write (*,*) 'and we''re back...'



! ok, so we have the energy histogram.  Next we need to load the energy-dispersed Lambert
! projections.  There is no kinematical part for now, so we only read a single array
open(unit=dataunit,file=FZfile,status='old',form='unformatted')
  read (dataunit)  nE,npx,npy
! make sure that the number of energybins is the same as for the energy histogram
  if (numEbins.ne.nE) then
    write (*,*) 'Energy histogram and Lambert stack have different energy dimension; aborting program'
    write (*,*) 'energy histogram = ',shape(accum_e)
    write (*,*) 'Lambert stack = ',nE, npx, npy
    stop
  end if
  allocate(sr(nE,npx,npy),stat=istat)
  read (dataunit) sr
close(unit=dataunit,status='keep')


! next, allocate all necessary arrays and fill them up
allocate(scin_x(numsx),scin_y(numsy),stat=istat)
! if (istat.ne.0) then ...
do i=1,numsx
  scin_x(i) = xpc*delta - ( (i-1) - numsx*0.5 + 0.5)*delta 
end do
do i=1,numsy
  scin_y(i) = ypc*delta - ( (i-1) - numsy*0.5 + 0.5)*delta 
end do

! auxilliary angle to rotate between reference frames
alp = 0.5*cPi - (sig - thetac)*dtor
ca = cos(alp)
sa = sin(alp)

! we will need to incorporate a series of possible distortions 
! here, as described in Gert nolze's paper

! compute auxilliary interpolation arrays
allocate(rgx(numsx,numsy), rgy(numsx,numsy), rgz(numsx,numsy), stat=istat)
! if (istat.ne.0) then ...
L2 = L*L
do j=1,numsy
  sx = L2+scin_y(j)*scin_y(j)
  Ls = ca*scin_y(j) + L*sa
  Lc = -sa*scin_y(j)+L*ca
  do i=1,numsx
   rhos = 1.0/sqrt(sx + scin_x(i)**2)
   rgx(i,j) = Ls*rhos
   rgy(i,j) = scin_x(i)*rhos
   rgz(i,j) = Lc*rhos
  end do
end do

! allocate the array that will hold the computed pattern
allocate(EBSDpattern(numsx,numsy),stat=istat)
! if (istat.ne.0) then ...
EBSDpattern = 0.0


! scale factors for Lambert projection; we'll only allow the modified Lambert projection
  xscale = float(npx-1)/sqrt(2.0)/sqrt(cPi)
  yscale = xscale  


  offx = (npx-1)/2+1
  offy = offx

  write (*,*) 'precomputed Lambert projection dimensions = ',npx, npy

! in the old version, we either computed the background model here, or 
! we would load a background pattern from file.  In this version, we are
! using the background that was computed by the MC program, and has 
! an energy histogram embedded in it, so we need to interpolate this 
! histogram to the pixels of the scintillator.








idum = -1		! to initialize the random number generator

! set up the binning parameters
binx = numsx/binning
biny = numsy/binning
bindx = 1.0/float(binning)**2
allocate(binned(binx,biny),stat=istat)

! declare TIFF variables in TIFF_global if needed
if (outputmode.eq.'tiff') then
 TIFF_nx = binx
 TIFF_ny = biny
 allocate(TIFF_Image(0:TIFF_nx-1,0:TIFF_ny-1))
end if

! and allocate an array to store, say, 100 images; after every 100, all of them
! are stored in tiff format;  this is done to make sure that disk I/O doesn't
! become a bottleneck.
allocate(imagestack(binx,biny,storemax),stat=istat)
imcnt = 1

! start the timer  [something's not quite right here...]
call Time_start

firstfile = 1
if (outputmode.eq.'data') then
 open(unit=dataunit,file=trim(prefix)//'.data',status='replace',form='unformatted',action='write')
end if

! ok, we're done with the pre work.  Now we can compute the EBSD patterns
prefactor = 0.25D0 * nAmpere * beamcurrent * dwelltime / dble(numelectrons)
write (*,*) 'intensity prefactor = ',prefactor
write (*,*) 'max count = ',maxval(accum_e)

do iang=1,numeuler
! convert the direction cosines to quaternions, include the 
! sample quaternion orientation, and then back to direction cosines...
! then convert these individually to the correct EBSD pattern location
        qq1 = conjg(quatang(1:4,iang)) !PGC quatang(iang) -> quatang(1:4,iang)
        qq2 = quatang(1:4,iang) !PGC quatang(iang) -> quatang(1:4,iang)
	do i=1,numsx
	    do j=1,numsy
!  do the coordinate transformation for this euler agle
             ! qq = quaternion( 0.0, rgx(i,j),rgy(i,j),rgz(i,j) )
              qq = (/ 0.0, rgx(i,j),rgy(i,j),rgz(i,j) /) !PGC
              qq3 = qq2 * (qq * qq1)
              !dc(1:3) = (/ qq3%b, qq3%c, qq3%d /) ! these are the direction cosines 
              dc(1:3) = (/ qq3(2), qq3(3), qq3(4) /) ! PGC changed
! make sure the third one is positive; if not, switch all 
              dc = dc/sqrt(sum(dc**2))
              if (dc(3).lt.0.0) dc = -dc
! and finally do the interpolation part using Lambert equal area projection coordinates
	  	dd = sqrt(2.0*(1.0-dc(3)))
                if ( (abs(dc(2)).le.dc(1)) .or. (dc(1).le.-abs(dc(2))) ) then
 		  ix = ( abs(dc(1))/dc(1) ) * sqrt(cPi) * 0.5 * dd * xscale
		  iy = 2.0/sqrt(cPi) * atan2(dc(2),abs(dc(1))) * dd * yscale
	        else
		  ix = 2.0/sqrt(cPi) * atan2(dc(1),abs(dc(2))) * dd * xscale
		  iy = ( abs(dc(2))/dc(2) ) * sqrt(cPi) * 0.5 * dd * yscale
		end if
! four-point interpolation (bi-quadratic)
! needs to be adapted to the Lambert projection circle ... 
              nix = int(offx+ix)-offx
              niy = int(offy+iy)-offy
              dx = ix-nix
              dy = iy-niy
              dxm = 1.0-dx
              dym = 1.0-dy
              nix = nix+offx
              niy = niy+offy
 ! interpolate the intensity 
              do k=Emin,Emax 
                EBSDpattern(i,j) = EBSDpattern(i,j) + accum_e(k,i,j) * ( sr(k,nix,niy) * dxm * dym + &
                                           sr(k,nix+1,niy) * dx * dym + sr(k,nix,niy+1) * dxm * dy + &
                                           sr(k,nix+1,niy+1) * dx * dy )
              end do
              EBSDpattern(i,j) = prefactor * EBSDpattern(i,j)
          end do
       end do

! add sampling noise (Poisson noise in this case, so multiplicative, sort of)
	do i=1,numsx
 	  do j=1,numsy
              EBSDpattern(i,j) = POIDEV(EBSDpattern(i,j),idum)
	  end do
	end do      
       
! next, apply the optical point spread function to this pattern
! [to optimize the speed, we may want to change the size of the array to a square shape]

! TO BE IMPLEMENTED


! and apply any camera binning
	if (binning.ne.1) then 
	  do i=1,numsx,binning
	    do j=1,numsy,binning
	        binned(i/binning+1,j/binning+1) = sum(EBSDpattern(i:i+binning-1,j:j+binning-1))
	    end do
	  end do  
	else
  	  binned = EBSDpattern
	end if
! and divide by binning^2
	binned = binned * bindx
	

! this must also be removed; we're not going to apply contrast/brightness acling at this point;
! we'll allow the user to do this in the visualization code.
! and finally, before saving the patterns, apply the contrast function and scale between 0 and 255
	EBSDmax = maxval(binned)
	EBSDmin = minval(binned)
	if (iang.eq.1) write (*,*) 'Limit Intensity Values : ',EBSDmin,EBSDmax
!	imagestack(1:binx,1:biny,imcnt) = int(255*(contrast * ( (binned - EBSDmin)/(EBSDmax-EBSDmin) - 1.0) + 1.0))

! the following is a "best fit" solution without any basis in physics and is used just as a place holder ... 
	imagestack(1:binx,1:biny,imcnt) = int(35 + (242-35) * ((binned - EBSDmin)/(EBSDmax-EBSDmin))**contrast )

	imcnt = imcnt+1


! this will need to become an HDF5 formatted file with all the program output
! it should be reable in IDL as well as DREAM.3D.

	if ((imcnt.gt.storemax).or.(iang.eq.numeuler)) then 
	  write (*,*) 'storing EBSD patterns at angle ',iang,' of ',numeuler
	  if (outputmode.eq.'tiff') then
! and spit the images out in tiff format ... 
 	    do i=1,imcnt-1
     	   	  write (filenumber,"(I4.4)") iang-(imcnt-1)+i
	   	  TIFF_filename = trim(prefix)//filenumber//'.tiff'
! allocate memory for image
! fill the image with whatever data you have (between 0 and 255)
   	         TIFF_Image(0:TIFF_nx-1,0:TIFF_ny-1) = imagestack(1:TIFF_nx,1:TIFF_ny,i)
! create the file
 	 	  call TIFF_Write_File 
	     end do
	  else
	  do i=1,imcnt-1
	    write (dataunit) char(imagestack(1:binx,1:biny,i))
 	  end do
!          close(unit=dataunit,status='keep')
	  imcnt = 1
	 end if
	end if
! that's it for this pattern ... 
end do

if (outputmode.eq.'data') then 
  close(unit=dataunit,status='keep')
end if

 call Time_stop(numeuler)


end subroutine ComputeEBSDPatterns

