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
! CTEMsoft2013:CTEMMCTKD.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMMCTKD 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Monte Carlo BSE simulation with Lambert projection for data storage
!
!> @detail Monte Carlo Electron Trajectory Simulation for Transmission Kikuchi Diffraction
!>	This version uses the modified Lambert projection to store 
!>	the MC output data, so that we are not dependent on a
!>	particular detector geometry.  We store the energy and direction 
!>	cosines of a BSE electron along with depth information in the 
!>	Lambert projection array, which needs to be sufficiently fine in
!>	terms of sampling so that we can deal with a detector that's relatively
!>	far away.
!
!> @date  08/16/13 MDG  1.0 split off from MC.f90 to cover TKD patterns
!--------------------------------------------------------------------------

program CTEMMCTKD
use local
use io

IMPLICIT NONE

character(fnlen)			:: nmlfile

integer(kind=irg)			:: numarg	!< number of command line arguments
integer(kind=irg)			:: iargc	!< external function for command line
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
        nmlfile = arg
        if (trim(nmlfile).eq.'-h') then
		mess = ' Program should be called as follows: '; call Message("(/A)")
		mess = '        CTEMMCTKD [nmlfile]'; call Message("(A)")
		mess = ' where nmlfile is an optional file name for the namelist file;'; call Message("(A)")
		mess = ' if absent, the default name ''CTEMMCTKD.nml'' will be used.'; call Message("(A/)")
		stop
	end if
else
	nmlfile = 'CTEMMCTKD.nml'    		! assign the default namelist file name
end if

! perform a Monte Carlo simulation
 call DoMCTKDsimulation(nmlfile)
 
end program CTEMMCTKD 
 
 !--------------------------------------------------------------------------
!
! SUBROUTINE:DoMCTKDsimulation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Perform the MCTKD simulation
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 07/23/13  MDG 3.0 complete rewrite
!> @date 07/30/13  MDG 3.1 added Patrick's code for double sample tilt (sigma, omega)
!> @date 08/16/13  MDG 4.0 modified for TKD
!--------------------------------------------------------------------------
subroutine DoMCTKDsimulation(nmlfile)

use local
use crystal
use crystalvars
use symmetry
use error
use io
use files
use diffraction, only:CalcWaveLength
use rng
use Lambert
use omp_lib
 
IMPLICIT NONE

character(fnlen),INTENT(IN)	:: nmlfile	! namelist file

integer(kind=irg)	:: NUMTHREADS	! possible number of threads 

! all geometrical parameters for the scintillator setup
real(kind=dbl)		:: sig		! sample tile angle [degrees]
integer(kind=irg)	:: numsx 	! number of Lambert points along x
integer(kind=irg)	:: numsy	! number of Lambert points along y
integer,parameter	:: k12 = selected_int_kind(15)
integer(kind=k12)	:: num_el	! total number of electrons to try
integer(kind=irg)	:: nthreads ! number of threads requested
real(kind=dbl)		:: EkeV	! electron energy in keV
real(kind=dbl)		:: Ehistmin ! minimum energy for energy histogram (in keV)
real(kind=dbl)		:: Ebinsize  ! binsize in keV
integer(kind=irg)	:: numEbins, numzbins
real(kind=dbl)		:: depthmax ! maximum depth for which to keep track of exit energy statistics [in nm]
real(kind=dbl)		:: depthstep ! stepsize for depth-energy accumulator array [in nm]

! material parameters
character(fnlen)	:: xtalname	! crystal structure name (needed for density computation)
character(60)		:: dataname	! output file name
real(kind=dbl)		:: Ze		! average atomic number
real(kind=dbl)		:: density	! density in g/cm^3
real(kind=dbl)		:: at_wt	! average atomic weight in g/mole
real(kind=sgl)		:: avZ, avA, dens, dxyVEP

! variable passing array
real(kind=dbl)		:: varpas(13) 
integer(kind=irg)	:: i, TID
real(kind=sgl)		:: io_real(3)

! variables used for parallel random number generator (based on http://http://jblevins.org/log/openmp)
type(rng_t), allocatable :: rngs(:)
integer(kind=irg)	:: primeseed

! various allocatable arrays
integer(kind=irg),allocatable	:: accum_e(:,:,:,:), acc_e(:,:,:,:), accum_Ez(:,:), acc_Ez(:,:), &
				   accum_VEP(:,:,:), acc_VEP(:,:,:)

! various 
integer(kind=irg)	:: istat, skip, numVEP, nxyVEP

! filename stuff
character(5)		:: extension

! define the IO namelist to facilitate passing variables to the program.
namelist  / MCTKDdata / xtalname, sig, numsx, numsy, num_el, primeseed, EkeV, &
			dataname, nthreads, Ehistmin, Ebinsize, depthmax, depthstep 
 
! define reasonable default values for the namelist parameters
    xtalname = 'undefined'
    sig = 0.0
    numsx = 512
    numsy = 512
    primeseed	= 932117
    num_el = 10000000_k12  ! try 10 million electrons unless otherwise specified   
    EkeV	= 30.D0
    dataname  = 'sample'
    nthreads = 2
    Ehistmin = 5.0D0
    Ebinsize = 0.25D0  ! we'll keep this fixed for now
    depthmax = 100.D0
    depthstep = 1.D0
    
! then we read the MCTKDdata namelist, which may override some of these defaults  
 OPEN(UNIT=dataunit,FILE=trim(nmlfile),DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=MCTKDdata)
 CLOSE(UNIT=dataunit)

 if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMMCTKD:',' structure file name is undefined in '//nmlfile)
 end if

! print some information
 progname = 'CTEMMCTKD.f90'
 progdesc = 'Monte Carlo Electron Trajectory Simulation for TKD+Lambert'
 call CTEMsoft

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)
 skip = 3
 call CalcWaveLength(EkeV*1000.D0,skip)

! then get the density, average atomic number and average atomic weight
 call CalcDensity(dens, avZ, avA)
 density = dble(dens)
 Ze = dble(avZ)
 at_wt = dble(avA)
 io_real(1:3) = (/ sngl(density), sngl(Ze), sngl(at_wt) /)
 call WriteValue('Density, avZ, avA = ',io_real,3)


! ok, now we're ready to start the computation
!
! define the number of virtual exit planes and allocate exit plane accumulator
numVEP = nint(depthmax/depthstep)+1
dxyVEP = 1.0		! 1 nm step size for exit plane patterns
nxyVEP = 150		! number of xy steps in one direction
allocate(accum_VEP(-nxyVEP:nxyVEP,-nxyVEP:nxyVEP,numVEP),stat=istat)

! allocate the accumulator array for the energy histogram
numEbins =  int((EkeV-Ehistmin)/Ebinsize)+1
allocate(accum_Ez(numEbins,numVEP),stat=istat)

! and here we store the angular distribution of all scattered electrons
! by energy and exit depth
allocate(accum_e(numVEP,numEbins,numsx,numsy),stat=istat)

! now put most of these variables in an array to be passed to the single_run subroutine
varpas = (/ sig, dble(numsx), dble(numsy), dble(num_el), EkeV, dble(nxyVEP), &
	    Ze, density, at_wt, Ehistmin, Ebinsize, depthmax, depthstep /)

! set the number of OpenMP threads and allocate the corresponding number of random number streams
call OMP_SET_NUM_THREADS(nthreads)
allocate(rngs(nthreads),stat=istat)

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL  PRIVATE(i,TID,istat,acc_VEP,acc_Ez,acc_e) &
!$OMP& SHARED(NUMTHREADS,varpas,accum_e,accum_Ez,accum_VEP,num_el,numEbins,numVEP)

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()

if (TID.eq.0) write (*,*) 'Number of available threads = ',NUMTHREADS

! allocate the local accumulator arrays for number of electrons and energy
allocate(acc_e(numVEP,numEbins,numsx,numsy),acc_Ez(numEbins,numVEP),acc_VEP(-nxyVEP:nxyVEP,-nxyVEP:nxyVEP,numVEP),stat=istat)

! each thread gets to execute the entire program once
!$OMP DO SCHEDULE(STATIC,1)    

do i=1,NUMTHREADS
! get a unique seed for this thread  
  call rng_seed(rngs(i), 932117 + i)

! do the Monte Carlo run  
  call single_run(varpas,rngs(i),acc_e,acc_Ez,acc_VEP)

! make sure that only one thread copies its contents into the accumulator arrays at any given time
!$OMP CRITICAL
  accum_e = accum_e + acc_e
  accum_Ez = accum_Ez + acc_Ez
  accum_VEP = accum_VEP + acc_VEP
!$OMP END CRITICAL  

end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL

! and here we create the output file
extension = '.data'
write (*,*) ' '
write (*,*) ' All threads complete; saving data to file ',trim(dataname)

write (*,*) ' Total number of electrons generated = ',num_el*NUMTHREADS
write (*,*) ' Number of electrons on detector       = ',sum(accum_e)

open(dataunit,file=trim(dataname),status='unknown',form='unformatted')
write (dataunit) numEbins, numVEP, nxyVEP, numsx, numsy
write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
write (dataunit) accum_e
write (dataunit) accum_Ez
write (dataunit) accum_VEP 
close(dataunit,status='keep')


contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:single_run
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does a full simulation starting from a given random number seed
!
!> @param varpas variable list
!> @param rngt random number identifier
!> @param accum_e energy accumulator array
!> @param accum_z depth accumulator array
!> @param numEbins number of energy bins
!> @param numzbins number of depth bins
!
!> @date 11/**/12  PGC 1.0 IDL version
!> @date 12/04/12  MDG 1.1 conversion to Fortran-90
!> @date 12/05/12  MDG 1.2 created subroutine to run with OpenMP
!> @date 03/11/13  MDG 2.0 modified for Lambert projection
!> @date 07/23/13  MDG 3.0 complete rewrite, integration with CTEMsoft libraries
!> @date 07/30/13  MDG 3.1 added Patrick's code for tilted sample surface (sigma, omega)
!> @date 07/31/13  MDG 3.2 corrected off-by-one error in energy binning
!--------------------------------------------------------------------------
recursive subroutine  single_run(varpas,rngt,accum_e,accum_Ez,accum_VEP)

use local
use rng
use Lambert

! all geometrical parameters for the scintillator setup
real(kind=dbl)		:: sig		! sample tile angle [degrees]
integer(kind=irg)	:: numsx 	! number of scintillator points along x
integer(kind=irg)	:: numsy	! number of scintillator points along y
real(kind=dbl)		:: Ehistmin ! minimum energy for energy histogram (in keV)
real(kind=dbl)		:: Ebinsize ! binsize in keV
integer(kind=irg)	:: numEbins, numzbins
real(kind=dbl)		:: Emin	! Ehistmin - Ebinsize/2
integer(kind=irg)	:: iE		! energy bin counter
integer(kind=irg)	:: iz		! exit depth bin counter
real(kind=dbl)		:: depthmax ! maximum depth for which to keep track of exit energy statistics [in nm]
real(kind=dbl)		:: depthstep ! stepsize for depth-energy accumulator array [in nm]

! Monte Carlo related parameters
real(kind=dbl)		:: EkeV, Ec	! electron energy in keV
real(kind=dbl)		:: scaled = 1.0D8 ! cm to Angstrom scalefactor
real(kind=dbl)		:: min_energy = 1.D0  ! in keV
real(kind=dbl)		:: presig = 1.5273987D19   ! = 1/( 5.21D-21 * (4*cPi) )
real(kind=dbl)		:: x, y, z, xn, yn, zn 	! electron coordinates
real(kind=dbl)		:: cx, cy, cz, r1, r2, cxstart, czstart			! direction cosines
real(kind=dbl)		:: alpha, phi, psi		! angles
real(kind=dbl)		:: lambda, step1, step2, pre, sige, prealpha, predEds	! stepsize parameters
real(kind=dbl)		:: J, dEds, dE, dEsave			! energy-related variables
real(kind=dbl)		:: cphi, sphi, cpsi, spsi, tpi	! cosines and sines and such ... 
real(kind=dbl)		:: cxp, cyp, czp, dsq, rho, dx, dy, edis	! trajectory parameters

integer(kind=irg)	:: idx, idy	! scintillator coordinates

! derived quantities
real(kind=sgl)		:: alphad	! angle between sample normal and scintillator normal
real(kind=sgl)		:: ca, sa, tana, cota, La, Lsa, cai, sai, deltai ! angle functions

! auxiliary variables
integer(kind=irg)	:: num 	! number of scattering events to try
integer(kind=irg)	:: bsct	! back-scattered electron counter
integer,parameter	:: k12 = selected_int_kind(15)
integer(kind=k12)	:: num_el	! total number of electrons to try
integer(kind=k12)	:: el		! electron counter
integer(kind=irg)	:: traj, ipxy(2), z1, z2, ierr, idxy(2), iVEP	! trajectory counter

! material parameters
real(kind=dbl)		:: Ze		! average atomic number
real(kind=dbl)		:: density	! density in g/cm^3
real(kind=dbl)		:: at_wt	! average atomic weight in g/mole

real(kind=dbl)		:: delta, xyz(3), step, cxyz(3), cxyzp(3), dsqi, dd, xyzn(3), lam, dxy(2)

! variable passing arrays
real(kind=dbl),INTENT(IN)		:: varpas(13)
integer(kind=irg),INTENT(OUT)		:: accum_e(:,:,:,:), accum_Ez(:,:),accum_VEP(:,:,:) 
integer(kind=irg)			:: TID, nx, izmax

! parallel random number variable 
type(rng_t), intent(inout) 		:: rngt
real(kind=dbl)				:: rr	! random number

real(kind=dbl), parameter 		:: cPi=3.141592653589D0, cAvogadro = 6.0221415D+23, cDtoR = 0.017453293D0

! initalize the Lambert projection parameters
call InitLambertParameters

! get the current thread number
TID = OMP_GET_THREAD_NUM()

accum_e = 0
accum_Ez = 0
accum_VEP = 0

! initialize all the variables based on the varpas array
 sig = varpas(1)
 numsx = int(varpas(2),kind=irg)
 numsy = int(varpas(3),kind=irg)
 nx = (numsx-1)/2 
 num_el = int(varpas(4),kind=k12)
 EkeV = varpas(5)
 nxyVEP = int(varpas(6),kind=irg)
 Ze = varpas(7) 
 density = varpas(8) 
 at_wt = varpas(9)
 Ehistmin = varpas(10)
 Ebinsize = varpas(11)
 depthmax = varpas(12)
 depthstep = varpas(13)
 
 Emin = Ehistmin - Ebinsize/2.D0
 numVEP =  int(depthmax/depthstep)+1
 numEbins =  int((EkeV-Ehistmin)/Ebinsize)+1

!    deltai = 1.D0/delta
    
! prefactor for mean free path computation
    pre =  at_wt/cAvogadro/density
    prealpha = 3.4D-3 * Ze**(0.67) 
    J = (9.76D0 * Ze+58.5D0 / Ze**(0.19D0) )*1.0D-3 / 1.166
    predEds = -78500.0D0 * density * Ze / at_wt 
    
! parameter for the Lambert projection scaling
 delta = dble(nx) / LPs%sPio2

! initialize the electron counter to one
 num = 1 
 tpi = 2.D0 * cPi

 cxstart = dcos( (90.D0-sig) * cDtoR)		! direction cosines for beam on tilted sample
 czstart = -dsin( (90.D0-sig) * cDtoR)

izmax = int(depthmax/depthstep)+1

! and here is the main loop  
mainloop: do el = 1,num_el

! every now and then, print something to the screen 
  if (mod(el,5000000_k12).eq.0_k12) then
    write (*,*) TID,': Completed electron # ',el,'; detector hits = ',sum(accum_e), bsct
  end if
  Ec = EkeV   ! set the energy for this incident electron 

! these could in principle be sampled from an area corresponding to the beam size
    xyz = (/ 0.D0, 0.D0, 0.D0 /)  	! initial coordinates

 ! What we really need to have here is to calculate the distance travel after entering sample, i.e. step*cz and step*cx
    alpha=prealpha / Ec
    step = Ec*(Ec+1024.D0)/Ze/(Ec+511.D0)		! step is used here as a dummy variable
    sige =  presig * step * step * alpha * (1.D0+alpha)
    lambda = pre * sige
    rr = rng_uniform(rngt)
    step = -lambda * log(rr)
 
    cxyz = (/ cxstart, 0.D0, czstart /)		! direction cosines for beam on tilted sample
 
! advance the coordinates 
    xyz = xyz + step * scaled * cxyz
 
  traj = 0
  trajloop: do while (traj.lt.num)

! Subtract the energy that is lost for path length lambda.
    dEds = predEds * dlog( Ec / J+1.0D0 ) / Ec
    dE = dEds * step
    Ec = Ec+dE

! here we exit the trajloop (using the f90 EXIT command) if the energy becomes low enough
    if (Ec.lt.min_energy) EXIT trajloop
    
!  Find the angle the electron is deflected through by the scattering event.
    rr = rng_uniform(rngt)
    cphi = 1.D0-2.D0*alpha*rr/(1.D0+alpha-rr)
    sphi = dsqrt(1.D0-cphi*cphi)

! Find the azimuthal scattering angle psi
    rr = rng_uniform(rngt)
    psi = tpi * rr
    spsi = dsin(psi)
    cpsi = dcos(psi)

! From MCML paper START
    if (dabs(cxyz(3)).gt.0.99999D0) then
      cxyzp = (/ sphi * cpsi, sphi * spsi, (cxyz(3)/dabs(cxyz(3))) * cphi /)
    else 
      dsq = dsqrt(1.D0-cxyz(3)*cxyz(3))
      dsqi = 1.D0/dsq
      cxyzp = (/ sphi * (cxyz(1) * cxyz(3) * cpsi - cxyz(2) * spsi) * dsqi + cxyz(1) * cphi, &
      		sphi * (cxyz(2) * cxyz(3) * cpsi + cxyz(1) * spsi) * dsqi + cxyz(2) * cphi, &
        	-sphi * cpsi * dsq + cxyz(3) * cphi /)
    end if
!  From MCML paper END

! normalize the direction cosines
    dd = 1.D0/dsqrt(sum(cxyzp*cxyzp))
    cxyzp = cxyzp * dd
    
! get the step size between scattering events (inline code rather than function call)
    alpha = prealpha / Ec
    step = Ec*(Ec+1024.D0)/Ze/(Ec+511.D0)		! step is used here as a dummy variable
    sige =  presig * step * step * alpha * (1.D0+alpha)
    lambda = pre * sige

    rr = rng_uniform(rngt)
    step = - lambda * log(rr)

! apply the step to the next location    
    xyzn = xyz + step * scaled * cxyzp 

! next we need to analyze these numbers and determine which bins to add them to ...

! first of all, if the electron is backscattered, then we add 1 to the BSEcounter
! and move on to the next electron
    if (xyzn(3).ge.0.D0) then
        bsct = bsct + 1  	! we have a back-scattered electron
	CYCLE mainloop	   	! skip the rest and go to the next electron
    else
! if we are moving in the forward (negative z) direction, 
! then determine which virtual exit plane is crossed, else
! move on to the next 
     if (cxyzp(3).ge.0.D0) then
       CYCLE trajloop
     else			
! ok, it's going forward; next determine whether or not it has CROSSED one
! of the virtual exit planes AND is traveling in the general direction of the detector 
	z1 = int(-xyz(3)*0.1D0/depthstep)
	z2 = int(-xyzn(3)*0.1D0/depthstep)
!	write (*,*) 'forward electron ? ',z1,z2,xyz(3),xyzn(3)
	if (z2.gt.z1) then 		! the electron has crossed a virtual exit plane
! add the electron to the virtual exit plane spatial distribution map
	  lam = (-10.0D0*dble(depthstep)*dble(z2) - xyz(3)) / (xyzn(3) - xyz(3))
	  ipxy = nint(0.1D0*(xyz(1:2) + lam * (xyzn(1:2) - xyz(1:2))))	! we're taking 1 nm intervals
	  if (maxval(abs(ipxy)).le.nxyVEP) then 
	    ipxy = ipxy + nxyVEP + 1
!write (*,*) 'accum_VEP coordinates ',ipxy, z2, nxyVEP
	    accum_VEP(ipxy(1),ipxy(2),z2) =  accum_VEP(ipxy(1),ipxy(2),z2) + 1
	  end if
!write (*,*) 'direction cosines : ',cxyzp,z1,z2 
         if (cxyzp(1).lt.0.D0) then 	! it is traveling in the halfspace that is projected onto the detector
! Let's figure out where in the Lambert array this point should be projected ...   
! We know the direction cosines were normalized, so no reason to check the error flag
	    dxy = delta * LambertSphereToSquare( -cxyzp, ierr )       
	
! and get the nearest pixel [ take into account reversal of coordinate frame (x,y) -> (y,-x) ]
            idxy = (/ nint(dxy(2)), nint(-dxy(1)) /)
! then add it to the other accumulator arrays if    
!write (*,*) ' -> idxy = ',idxy, nx
            if (maxval(abs(idxy)).le.nx) then
! determine the energy bin
               iE = nint((Ec-Ehistmin)/Ebinsize)+1
! determine the depth bin, taking into account the dual storage in the modified Lambert array
	       iVEP = (z2-1)/2 + 1
	       if (mod(z2,2).eq.1) then 		! ! z2 is odd, so change sign of idxy(2)
		  idxy(2) = -idxy(2)
	       end if 	
	       idxy = idxy + nx + 1
! first add this electron to the correct exit distance vs. energy bin
	       edis = dabs(xyz(3)/cxyzp(3))   ! distance from last scattering point to surface along trajectory
	       iz = abs(nint(edis*0.1D0/depthstep)) +1
!write(*,*) ' idxy = ',idxy,edis*0.1D0,iVEP,z2
	       if ( (iz.gt.0).and.(iz.le.numVEP) ) then
                accum_Ez(iE,iz) = accum_Ez(iE,iz) + 1
              end if
! then add it to the modified Lambert accumulator array.
              accum_e(iVEP,iE,idxy(1),idxy(2)) = accum_e(iVEP,iE,idxy(1),idxy(2)) + 1
            end if
          end if
	  if (z2.gt.izmax) EXIT trajloop
        end if
     end if    
   end if  ! xyzn(3).lt.0.D0

! update the electron direction cosines, coordinates and the interaction event counter
   cxyz = cxyzp
   xyz = xyzn
   traj = traj + 1
  

  end do trajloop

end do mainloop

end subroutine single_run

end subroutine DoMCTKDsimulation
