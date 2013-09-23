!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MC.f90
! Copyright (c) 2012  Patrick Callahan/Marc De Graef
! ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -Fortran-90 Source Code-
!  
!  Description: Monte Carlo Electron Trajectory Simulation for EBSD
!                
!
!
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
! 11/**/12   PGC 1.0 IDL version
! 12/04/12  MDG 1.1 conversion to Fortran-90
! 12/06/12  MDG 1.2 conversion to OpenMP, with new random number generator
! 12/06/12  MDG 1.3 added energy histogram sampling
! 12/07/12  MDG 1.4 added energy vs. depth sampling
! ###################################################################
! 

program MC

use local
use rng
!use timing
use omp_lib

IMPLICIT NONE

integer(kind=irg)	:: NUMTHREADS	! possible number of threads 

! all geometrical parameters for the scintillator setup
real(kind=dbl)		:: L		! distance between scintillator screen and interaction point [microns] (~working distance)
real(kind=dbl)		:: sig		! sample tile angle [degrees]
real(kind=dbl)		:: omega	! sample rotation about RD angle [degrees] PGC
real(kind=dbl)		:: thetac	! detector tilt angle below horizontal [degrees]
real(kind=dbl)		:: delta	! scintillator step size [microns]
integer(kind=irg)	:: numsx 	! number of scintillator points along x
integer(kind=irg)	:: numsy	! number of scintillator points along y
real(kind=dbl)		:: xpc	! pattern center x [pixels]
real(kind=dbl)		:: ypc	! pattern center y [pixels]
integer(kind=irg)	:: num_el	! total number of electrons to try
integer(kind=irg)	:: nthreads ! number of threads requested
real(kind=dbl)		:: EkeV	! electron energy in keV
real(kind=dbl)		:: Ehistmin ! minimum energy for energy histogram (in keV)
real(kind=dbl)		:: Ebinsize  ! binsize in keV
integer(kind=irg)	:: numEbins, numzbins
real(kind=dbl)		:: depthmax ! maximum depth for which to keep track of exit energy statistics [in nm]
real(kind=dbl)		:: depthstep ! stepsize for depth-energy accumulator array [in nm]

! material parameters
character(30)		:: matl	! material descriptor
real(kind=dbl)		:: Ze		! average atomic number
real(kind=dbl)		:: density	! density in g/cm^3
real(kind=dbl)		:: at_wt	! average atomic weight in g/mole

! variable passing array
real(kind=dbl)		:: varpas(18) 
integer(kind=irg)	:: i, TID

! variables used for parallel random number generator (based on http://http://jblevins.org/log/openmp)
type(rng_t), allocatable :: rngs(:)
integer(kind=irg)	:: primeseed

! various allocatable arrays, energy histogram is first index, x,y on scintillator 2nd and 3rd indices
integer(kind=irg),allocatable	:: accum_e(:,:,:), acc_e(:,:,:), accum_z(:,:), acc_z(:,:)

! various 
integer(kind=irg)	:: istat

! filename stuff
character(5)		:: extension

! define the IO namelist to facilitate passing variables to the program.
namelist  / MCdata / L, sig, thetac, delta, numsx, numsy, xpc, ypc, num_el, primeseed, EkeV, &
				 matl, Ze, density, at_wt, nthreads, Ehistmin, depthmax, depthstep, omega 
 
! spit out some basic information
progname = 'MC.f90'
progdesc = 'Monte Carlo Electron Trajectory Simulation for EBSD'
call ProgramNote 

! define reasonable default values for the namelist parameters
    L = 15250.0
    sig = 70.0
    thetac = 0.0
    delta = 49.375
    numsx = 640
    numsy = 480
    xpc = 3.5679
    ypc =  113.449
    primeseed	= 932117
    num_el = 1000  ! try 10 million electrons unless otherwise specified   
    EkeV	= 30.D0
    matl	= 'Nickel'
    Ze		= 28.D0
    density	= 8.908D0
    at_wt	= 58.6934D0
    Ehistmin = EkeV - 10.0D0
    Ebinsize = 0.25D0  ! we'll keep this fixed for now
    omega = 2.0 !PGC
    
! then we read the MCdata namelist, which may override some of these defaults  
! OPEN(UNIT=dataunit,FILE='MC.nml',DELIM='APOSTROPHE')
 OPEN(UNIT=dataunit,FILE='MC.nml',DELIM='APOSTROPHE')

 READ(UNIT=dataunit,NML=MCdata)
 CLOSE(UNIT=dataunit)


! allocate the accumulator arrays for number of electrons and energy
numEbins =  int((EkeV-Ehistmin)/Ebinsize)+1
numzbins =  int(depthmax/depthstep)+1
allocate(accum_e(numEbins,numsx,numsy),accum_z(numEbins,numzbins),stat=istat)


! now put most of these variables in an array to be passed to the single_run subroutine
varpas = (/ L, sig, thetac, delta, dble(numsx), dble(numsy), xpc, ypc, dble(num_el), EkeV, &
		   Ze, density, at_wt, Ehistmin, Ebinsize, depthmax, depthstep, omega/)


! set the number of OpenMP threads and allocate the corresponding number of random number streams
call OMP_SET_NUM_THREADS(nthreads)
allocate(rngs(nthreads),stat=istat)

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL  PRIVATE(i,acc_e,acc_z,TID,istat) &
!$OMP& SHARED(NUMTHREADS,varpas,accum_e,accum_z,num_el,numEbins,numzbins)

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()

if (TID.eq.0) write (*,*) 'Number of available threads = ',NUMTHREADS

! allocate the accumulator arrays for number of electrons and energy
allocate(acc_e(numEbins,numsx,numsy),acc_z(numEbins,numzbins),stat=istat)

! each thread gets to execute the entire program once
!$OMP DO SCHEDULE(STATIC,1)    

! start the timer in thread 0 only
!if (TID.eq.0) call Time_start 

do i=1,NUMTHREADS
! get a unique seed for this thread  
  call rng_seed(rngs(i), 932117 + i)

! do the Monte Carlo run  
  call single_run(varpas,rngs(i),acc_e,acc_z)

! make sure that only one thread copies its contents into the accumulator arrays at any given time
!$OMP CRITICAL
  accum_e = accum_e + acc_e
  accum_z = accum_z + acc_z
!$OMP END CRITICAL  

end do
!$OMP END DO

! and report the timing parameters (after waiting for all threads to end)
!if (TID.eq.0) call Time_stop(num_el)

!$OMP END PARALLEL

! and here we create the output file
extension = '.data'
write (*,*) ' '
write (*,*) ' All threads complete; saving data to file ',trim(matl)//extension

write (*,*) ' Total number of electrons generated = ',num_el*NUMTHREADS
write (*,*) ' Number of electrons on detector       = ',sum(accum_e)

open(dataunit,file=trim(matl)//extension,status='unknown',form='unformatted')
write(dataunit) numEbins, numzbins, numsx, numsy
write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
write(dataunit) accum_e
write (dataunit) accum_z
close(dataunit,status='keep')


contains

! ###################################################################
! subroutine single_run
!
! does a full simulation starting from a given random number seed 
!
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
! 11/**/12   PGC 1.0 IDL version
! 12/04/12  MDG 1.1 conversion to Fortran-90
! 12/05/12  MDG 1.2 created subroutine to run with OpenMP
! ###################################################################
!
! at first, this didn't work at all, and the results were only correct for a single thread, but
! completely wrong for multiple threads.  After reading up on threads and such, I 
! figured out that it is necessary to tell the compiler to create memory allocations for 
! all variables for each thread separately...  This can be done by putting the "recursive"
! keyword in front of the subroutine; this forces allocation of separate memory for each
! invocation of the routine, i.e. for each thread, so that the threads become truly 
! independent of each other...
!
recursive subroutine  single_run(varpas,rng,accum_e,accum_z)

use local
use rng
!use timing

! all geometrical parameters for the scintillator setup
real(kind=dbl)		:: L		! distance between scintillator screen and interaction point [microns] (~working distance)
real(kind=dbl)		:: sig		! sample tile angle [degrees]
real(kind=dbl)		:: omega	! sample rotation about RD angle [degrees] PGC
real(kind=dbl)		:: thetac	! detector tilt angle below horizontal [degrees]
real(kind=dbl)		:: delta	! scintillator step size [microns]
integer(kind=irg)	:: numsx 	! number of scintillator points along x
integer(kind=irg)	:: numsy	! number of scintillator points along y
real(kind=dbl)		:: xpc	! pattern center x [pixels]
real(kind=dbl)		:: ypc	! pattern center y [pixels]
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
real(kind=dbl)		:: lambda, step, pre, sige, prealpha, predEds	! stepsize parameters
real(kind=dbl)		:: J, dEds, dE			! energy-related variables
real(kind=dbl)		:: cphi, sphi, cpsi, spsi, tpi	! cosines and sines and such ... 
real(kind=dbl)		:: cxp, cyp, czp, dsq, rho, dx, dy, edis	! trajectory parameters
real(kind=dbl)		:: cxpsave, cypsave, czpsave ! save these for calculating new basis PGC
real(kind=dbl)		:: nx, ny, nz, znmax ! PGC, define sample plane normal (From sig and omega)

integer(kind=irg)	:: idx, idy		! scintillator coordinates

! derived quantities
real(kind=sgl)		:: alphad	! angle between sample normal and scintillator normal
real(kind=sgl)		:: ca, sa, tana, cota, La, Lsa, cai, sai, deltai, cos70, sin70 ! angle functions

! auxiliary variables
integer(kind=irg)			:: num 	! number of scattering events to try
integer(kind=irg)			:: bsct	! back-scattered electron counter
integer(kind=irg)			:: num_el	! total number of electrons to try
integer(kind=irg)			:: el		! electron counter
integer(kind=irg)			:: traj	! trajectory counter

! material parameters
real(kind=dbl)				:: Ze		! average atomic number
real(kind=dbl)				:: density	! density in g/cm^3
real(kind=dbl)				:: at_wt	! average atomic weight in g/mole

! variable passing arrays
real(kind=dbl),INTENT(IN)		:: varpas(18)
integer(kind=irg),INTENT(OUT)	:: accum_e(:,:,:), accum_z(:,:)
integer(kind=irg)				:: TID

! parallel random number variable 
type(rng_t), intent(inout) 		:: rng
real(kind=dbl)				:: rr	! random number

real(kind=dbl), parameter :: cPi=3.141592653589D0, cAvogadro = 6.0221415D+23, cDtoR = 0.017453293D0


! get the current thread number
TID = OMP_GET_THREAD_NUM()

! initialize all the variables based on the varpas array
 L = varpas(1)
 sig = varpas(2)
 thetac = varpas(3)
 delta = varpas(4)
 numsx = int(varpas(5),kind=irg)
 numsy = int(varpas(6),kind=irg) 
 xpc = varpas(7) * delta
 ypc = varpas(8) * delta
 num_el = int(varpas(9),kind=irg)
 EkeV = varpas(10)
 Ze = varpas(11) 
 density = varpas(12) 
 at_wt = varpas(13)
 Ehistmin = varpas(14)
 Ebinsize = varpas(15)
 depthmax = varpas(16)
 depthstep = varpas(17)
 omega = varpas(18)
 
 Emin = Ehistmin - Ebinsize/2.D0
 numzbins =  int(depthmax/depthstep)+1
 numEbins =  int((EkeV-Ehistmin)/Ebinsize)+1

! compute angle quantities
    alphad = 90.D0-sig+thetac
    ca = dcos(alphad*cDtoR)
    sa = dsin(alphad*cDtoR)
    cai = 1.D0/ca
    sai = 1.D0/sa
    tana = dtan(alphad*cDtoR)
    cota = 1.D0/tana
    La = L * (tana+cota) 
    Lsa = L * sa
    deltai = 1.D0/delta
    cos70 = dcos(70.D0*cDtoR)
    sin70 = dsin(70.D0*cDtoR)
    
! prefactor for mean free path computation
    pre =  at_wt/cAvogadro/density
    prealpha = 3.4D-3 * Ze**(0.67) 
    J = (9.76D0 * Ze+58.5D0 / Ze**(0.19D0) )*1.0D-3 / 1.166
    predEds = -78500.0D0 * density * Ze / at_wt 
       
! Compute the plane normal for the sample PGC
	!Used Mathematica, got
	nx = -1*dcos(omega*cDtoR)*dsin(sig*cDtoR)
    ny = -1*dsin(omega*cDtoR)
    nz = dcos(omega*cDtoR)*dcos(sig*cDtoR)
    !write (*,*) nx,ny,nz
       
! initialize the electron counter to zero
 num = 15000  ! number of scattering events to try. Try 15000 ish to end from energy loss. 1 for ECCI
 tpi = 2.D0 * cPi

! cxstart = dcos( (90.D0-sig) * cDtoR)		! direction cosines for beam on tilted sample
! czstart = -dsin( (90.D0-sig) * cDtoR)
 cxstart = 0.D0		! direction cosines for beam on tilted sample
 czstart = -1.D0


! and here is the main loop
mainloop: do el = 1,num_el

! every now and then, print something to the screen
  if (mod(el,num_el/100).eq.0) then
    if (TID.eq.0) then
      write (*,*) ' '
      !call Time_remaining(el,num_el)
    end if
    write (*,*) TID,': Completed electron # ',el,'; bse hits = ',sum(accum_e),';num bsct =',bsct
  end if
  Ec = EkeV   ! set the energy for this incident electron

  x = 0.D0  	! initial coordinates
  y = 0.D0
! changed from Patrick's code so that the electron starts near the surface...
  !z = -1.D-1  !MDG
  z = 0.D0   !PGC
 
 ! What we really need to have here is to calculate the distance travel after entering sample, i.e. step*cz and step*cx
    alpha=prealpha / Ec
    sige =  presig * ( Ec*(Ec+1024.D0)/Ze/(Ec+511.D0) )**2 * (alpha*(1.D0+alpha))
    lambda = pre * sige
    rr =rng_uniform(rng)
    step = - lambda * log(rr)
 
    cx = cxstart		! direction cosines for beam on tilted sample
    cy = 0.D0
    cz = czstart
 
    step = step * scaled 
    x = x + step * cx
    y = y + step * cy
    z = z + step * cz


  traj = 0
  
  trajloop: do while (traj.lt.num)
    ! get the step size between scattering events (inline code rather than function call)
    alpha=prealpha / Ec
    sige =  presig * ( Ec*(Ec+1024.D0)/Ze/(Ec+511.D0) )**2 * (alpha*(1.D0+alpha))
    lambda = pre * sige

    rr =rng_uniform(rng)
    step = - lambda * log(rr)
    
    !  Find the angle the electron is deflected through by the scattering event.
    rr =rng_uniform(rng)
    phi = dacos(1.D0-2.D0*alpha*rr/(1.D0+alpha-rr))

    ! Find the azimuthal scattering angle psi
    rr =rng_uniform(rng)
    psi = tpi * rr

    ! Subtract the energy that is lost for path length lambda.
    dEds = predEds * alog( Ec / J+1.0D0 ) / Ec
    dE = dEds * step
    Ec = Ec+dE
    if (Ec.lt.min_energy) then
       traj=num   ! here we exit the trajloop (using the f90 EXIT command)
!       EXIT trajloop
    end if
    
! From MCML paper START
    sphi = dsin(phi)
    cphi = dcos(phi)
    spsi = dsin(psi)
    cpsi = dcos(psi)
    if (dabs(cz).gt.0.99999D0) then
      cxp = sphi * cpsi
      cyp = sphi * spsi
      czp = (cz/dabs(cz)) * cphi
    else 
      dsq = dsqrt(1.D0-cz**2)
      cxp = sphi * (cx * cz * cpsi - cy * spsi)/dsq+cx * cphi
      cyp = sphi * (cy * cz * cpsi + cx * spsi)/dsq+cy * cphi
      czp = -sphi * cpsi * dsq + cz * cphi
    end if
!  From MCML paper END

    step = step * scaled 
    xn = x + step * cxp 
    yn = y + step * cyp
    zn = z + step * czp


! HERE PGC
	znmax = -1*(nx*xn+ny*yn)/nz
	if (zn.gt.znmax) then
!	write (*,*) 'zn',zn,'znmax',znmax

!    if (zn.gt.0.D0) then
        bsct = bsct + 1  ! we have a back-scattered electron
        
        cxpsave=cxp
        cypsave=cyp
        czpsave=czp
        ! PGC Change to RD, TD, ND basis
        cxp=-1*cxpsave*cos70-czpsave*sin70
        cyp=-1*cypsave
        czp=-1*sin70*cxpsave+cos70*czpsave
        
        
        ! so let's see if it intersects the scintillator ... 
        rho = La / (czp*sai + cxp*cai)
        r1 = cxp * rho
        r2 = cyp * rho

	! convert to scintillator coordinates
        dx = (xpc - r2)*deltai 
        dy = (ypc - (r1 - Lsa)*cai )*deltai

	! and get the nearest scintillator pixel
        idx = nint(dx) + numsx/2
        idy = nint(dy) + numsy/2
        
        if (((idx.le.numsx).and.(idx.gt.0)).and.((idy.le.numsy).and.(idy.gt.0))) then
! If Ec larger than Emin, then we should count this electron
           if (Ec.gt.Emin) then     
             iE = int((Ec-Emin)/Ebinsize)+1
! first add this electron to the correct exit distance vs. energy bin
	     edis = dabs(z/czp)   ! distance from last scattering point to surface along trajectory
	     iz = int(edis*0.1D0/depthstep) +1
	     if ( (iz.gt.0).and.(iz.le.numzbins) ) accum_z(iE,iz) = accum_z(iE,iz) + 1
! then add it to the scintillator accumulator array.
             accum_e(iE,idx,idy) = accum_e(iE,idx,idy) + 1
           end if
        end if
        
        EXIT trajloop  ! and exit trajloop
    end if  ! zn.gt.0.D0

! update the electron coordinates and the interaction event counter
   cx = cxp
   cy = cyp
   cz = czp

   x = xn
   y = yn
   z = zn
   traj = traj + 1
  
  end do trajloop

end do mainloop

end subroutine single_run

end program MC
