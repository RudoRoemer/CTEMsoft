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
! CTEMsoft2013:CTEMSRdefect.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMSRdefect 
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Computes systematic row defect contrast for multiple defects
!                      using CTEM or STEM illumination and detector conditions 
!
!> @todo OpenMP implementation; overall verification of all coordinate transformations
! 
!> @date  01/14/10 MDG  1.0 original, based on SRdefect.f90 (collaboration with OSU) 
!> @date  05/14/13 MDG 2.1 replaced all IO by namelist file and added command line argument handling
!--------------------------------------------------------------------------
program CTEMSRdefect 
! 
use local

IMPLICIT NONE

character(fnlen)			:: nmlfile

integer(kind=irg)			:: numarg		!< number of command line arguments
integer(kind=irg)			:: iargc		!< external function for command line
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
        nmlfile = trim(arg)
        if (trim(nmlfile).eq.'-h') then
		mess = ' Program should be called as follows: '; call Message("(/A)")
		mess = '        CTEMSRdefect [nmlfile]'; call Message("(A)")
		mess = ' where nmlfile is an optional file name for the namelist file;'; call Message("(A)")
		mess = ' if absent, the default name ''CTEMSRdefect.nml'' will be used.'; call Message("(A/)")
		mess = ' To create templates of all possible input files, type CTEMSRdefect -init'; call Message("(A/)")
		stop ' end of program help information'
	end if
        if (trim(nmlfile).eq.'-init') then
! with this option the program creates template namelist files in the current folder so that the 
! user can edit them (file extension will be .nml.template)
		mess = ' Program will create template .nml files that can then be edited by the user.'; call Message("(A/)")
! to be written
		stop 'template files created'
	end if
else
	nmlfile = 'CTEMSRdefect.nml'    		! assign the default namelist file name
end if

! initialize all user-defined variables
call ComputeSRdefect(nmlfile)

end program CTEMSRdefect


!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeSRdefect
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute a systematic row defect image
!
!> @param nmlfile namelist file name
!
!> @todo check sign of kt for consistency 
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!--------------------------------------------------------------------------
subroutine ComputeSRdefect(nmlfile)

use local
use crystalvars 
use crystal
use symmetryvars
use symmetry
use postscript
use constants
use diffraction
use dynamical
use files
use io
use error
use foilmodule
use stacking_fault
use dislocation
use void
use inclusion
use defectmodule     
use TIFF_f90
use pgm
use timing
use STEMmodule
use quaternions

IMPLICIT NONE

character(fnlen),INTENT(IN)		:: nmlfile 
character(fnlen)				::  STEMnmlfile, foilnmlfile

integer(kind=irg)       :: ira,nn,izero,i,j,k,n,nsl,numi,npix,npiy,ii,jj,numvoids,numdisl,numsf, skip, &
                                  numinc,dinfo,t_interval,outputfirst,outputlast,DF_nums_new, io_int(3), grange, &
                                  DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, NTHR, jcnt, numCL, iCL, SR_g(3)
!                                  OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
real(kind=sgl)          :: ind(3),hkl(3),exerg,cosom,glen,exer,dgr,sl,thr,arg,zmax,thick,mi,ma, &
                                  att,xgp,DF_gf(3),DF_gd(3,maxdefects), io_real(1), voltage, kt
character(50)           :: fname,sgname,voidname,dislname(3*maxdefects),sfname(maxdefects),outputroot, &
                                  incname,dispfile
character(fnlen)        :: imname, xtalname
character(4)            :: outputformat, illumination_mode, dispmode
character(5)            :: fnumber(-10:10)
character(len=10)       :: ci

complex(kind=dbl),allocatable    :: DHWM(:,:),Afirst(:,:),DHWMvoid(:,:)
complex(kind=dbl),allocatable    :: q(:,:),qin(:,:),qout(:,:),r(:,:),amp(:),amp2(:)
complex(kind=dbl)                :: czero=cmplx(0.0,0.0,dbl),cone=cmplx(1.0,0.0,dbl)
real(kind=sgl),allocatable       :: weights(:), inten(:), STEMimages(:,:,:,:)
integer(kind=irg),allocatable    :: disparray(:,:,:)
real(kind=sgl),allocatable       :: BFweightsarray(:,:,:),ADFweightsarray(:,:,:)


namelist / rundata / DF_L, DF_npix, DF_npiy, DF_slice, Nmat, sgname, numvoids, incname, &
                                 voidname, numdisl, dislname, numsf, sfname, dinfo, outputformat, &
				 outputroot,t_interval,illumination_mode,outputfirst,outputlast, dispfile, &
				 dispmode, xtalname, voltage, SR_g, grange, STEMnmlfile, kt, exerg, foilnmlfile

! here we read the general simulation information from a namelist file SRdef_rundata.nml
! first we define the default values

! parameters specific to this run
 xtalname = 'undefined'		! initial value; MUST be present in nml file for program to execute
 voltage = 200000.0		! accelerating voltage
 SR_g = (/ 1,0,0 /)			! systematic row g-vector
 grange = 5				! maximum positive multiple of g-vector; total number will be 2*grange+1 
 kt = 0.5					! tangential component of wave vector in units of |G|  (for STEM)
 exerg = 0.0				! excitation error in nm^-1 for fundamental reflection (for CTEM)

! CTEM or STEM ?
 illumination_mode = 'CTEM'  	! default illumination mode (can be 'CTEM' or 'STEM')
 STEMnmlfile = 'STEM_rundata.nml'	! name of the STEM rundata namelist file
 foilnmlfile = 'SRdef_foildata.nml'	! name of the foil rundata namelist file
 
! column approximation parameters and image parameters 
 DF_L = 1.0             			! edge length of column in nanometers
 DF_npix = 256       			! number of image pixels along x
 DF_npiy = 256       			! number of image pixels along y 
 DF_slice = 1.0       			! slice thickness in nanometers
 Nmat = 10000       			! number of precomputed A matrices to be stored

 dinfo = 0               			! switch to make makedislocation verbose
 sgname = 'nofile'   			! if this variable is different from 'nofile', then an external sg array is read (to be implemented)

! defect parameters
 numdisl = 0           			! number of dislocation files
 numsf = 0             			! number of stacking fault files
 numinc = 0           			! number of inclusions
 numvoids = 0       			! number of voids
 voidname = 'none' 			! filename for void data
 dislname = ''         			! filenames for dislocation data
 sfname = ''            			! filenames for stacking fault data
 incname = 'none'   			! filename for inclusion data
 dispfile = 'none'     			! name of the displacement field output file (will be created if different from none)
 dispmode = 'not'  			! should a diplacement file be written ('new') or read ('old') or neither ('not')?

! output parameters
 outputformat = 'data' 		! format for output data, can be 'data', 'pgm', or 'tiff'
 outputroot = 'image'  		! default root name for output files
 outputfirst = 1       			! first image number to be written to file (will be set to first image in SR)
 outputlast = nn      		! last image number to be written to file (will be set to last image in SR)
 t_interval = 10       			! default timing interval (output every t_interval image columns)
 
! then we read the rundata namelist, which may override some of these defaults  
 OPEN(UNIT=dataunit,FILE=nmlfile,DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=rundata)
 CLOSE(UNIT=dataunit)

! make sure the xtalname variable has been properly defined
if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMSRdefect:',' structure file name is undefined in '//nmlfile)
end if

! display the standard program info
 progname = 'CTEMSRdefect.f90'
 progdesc = 'Systematic row CTEM/STEM defect image simulation'
 call CTEMsoft

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

 ! generate all atom positions
 call CalcPositions('v')
 
! use systematic row vector to compute G* as in eq. 8.28.
  DF_gf = float(SR_g)
  DF_gstar = DF_gf/CalcLength(DF_gf,'r')**2    ! define G* such that G.G* = 1
  call TransSpace(DF_gf,DF_gc,'r','c')                 ! convert to Cartesian reference frame

! compute total number of beams
  nn = 2*grange+1  ! total number of beams
  izero = (nn+1)/2   ! id number of g=0 beam 
  
! read the STEM parameters from a namelist file and initialize all STEM related arrays
if (illumination_mode.eq.'STEM') then
  call read_STEM_data(STEMnmlfile, nn, SR_g, kt, i)  
  mess = 'Found and read STEM namelist file'; call Message("(A)")
end if


! allocate and initialize DF_Sarray, theta, and DF_Svoid
  allocate(DF_Sarray(0:Nmat-1,nn,nn),theta(-nn:nn),DF_Svoid(nn,nn))
  DF_Sarray = czero; theta = czero

! allocate the various DHW Matrices
  allocate(DHWMz(nn,nn),DHWM(nn,nn),DHWMvoid(nn,nn))
  DHWMvoid = czero; DHWMz=czero; DHWM(nn,nn)=czero

  
! Compute the off-diagonal part of the complex DHW matrix (factor i is included)
! We can precompute those because they will not change at all during the run
!       (these lines implement the equations on page 476 of the CTEM book)
 ind = float(SR_g)
  do i=1,nn
   do j=1,nn
    hkl=(-grange+i-1)*ind-(-grange+j-1)*ind     ! difference vector
    if (i.ne.j) then
     call CalcUcg(int(hkl))                         ! compute the interaction parameters
     DHWMz(i,j) = cPi*cmplx(-aimag(rlp%qg),real(rlp%qg),dbl)  ! and initalize the off-diagonal matrix element
    else
     DHWMz(i,j) = czero                           ! for now at least; this will be filled in later
    endif
   end do
  end do
  mess = 'Reference Darwin-Howie-Whelan matrix initialized'; call Message("(A/)")
  
! display the diffraction information for the fundamental reflection of the systematic row
  call CalcUcg(SR_g)
  io_real(1) = rlp%xg
  call WriteValue('Extinction distance for g : ', io_real, 1, "(F10.5)")
  io_real(1) = rlp%xgp
  call WriteValue('Anomalous absorption length for g : ', io_real, 1, "(F10.5)")
  io_real(1) = rlp%xgp/rlp%xg
  call WriteValue('Absorption Ratio : ', io_real, 1, "(F10.5)")
! compute the normal absorption factor xgp
  call CalcUcg((/0,0,0/))
  xgp = aimag(rlp%qg)
  io_real(1) = 1.0/xgp
  call WriteValue('Normal absorption length : ', io_real, 1, "(F10.5/)")
  

! set up the excitation errors for CTEM illumination mode
! for STEM mode we'll use the precomputed values from read_STEM_data later on...
if (illumination_mode.eq.'CTEM') then  
! show the user what the +g, -g, and g-3g excitation errors are
!  glen = CalcLength(ind,'r')
!  exerg = -mLambda*0.5*glen**2
!  mess = ' Symmetric systematic row orientation for sg = '; oi_real(1)=exerg; call WriteReal(1,"(F10.5)")
!  exerg = -2.0*mLambda*glen**2/(2.0+mLambda**2*glen**2)
!  mess = ' Systematic row orientation (s_-g=0) for sg    = '; oi_real(1)=exerg; call WriteReal(1,"(F10.5)")
!  exerg = mLambda*glen**2/(1.0-1.5*mLambda**2*glen**2) 
!  mess = ' Weak beam row orientation (s_3g=0) for sg   = '; oi_real(1)=exerg; call WriteReal(1,"(F10.5/)")
!! read the orientation for the entire systematic row; in STEM mode, this will be the 
!! value of sg for the center of the diffraction disk
!  mess = 'Excitation error for reflection g [nm^-1] = '; call GetReal(1); exerg = io_real(1)

! fill the diagonal of the reference dynamical matrix and the void matrix
  glen = CalcLength(ind,'r')
  cosom = -(mLambda*glen**2+2.D0*exerg)/(2.D0*glen*(mLambda*exerg+1.D0))   ! cos(omega)
  do i=1,nn
   n = -grange+i-1
   exer = -0.5*n*glen*(2.D0*cosom+n*mLambda*glen)/(1.D0+n*mLambda*glen*cosom)   ! uses eq. (8.30)
   DHWMz(i,i)=2.0*cPi*cmplx(0.0,exer)    ! initialize the diagonal element of the dynamical matrix
   DHWMvoid(i,i) = DHWMz(i,i)
  end do
end if

! next, we read the foildata namelist from the foil namelist file
! this includes material property data, in this case the elastic moduli
call read_foil_data(foilnmlfile,DF_npix,DF_npiy,DF_L)
    
! define the foil thickness, attenuation, and number slices per column
thick = foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
att = exp(-2.0*cPi*thick*xgp)  ! this is the global attenuation factor; remove the factor 2.0 if amplitudes are needed
DF_nums = nint(thick/DF_slice)  ! this is the number of slices for this particular column


! next, deal with all the defects
!
! if there is a diplacement field file entered in the STEM_rundata.nml file,  
! then we simply read that file in; otherwise, we read all the defect descriptor files
if ((dispmode.eq.'new').or.(dispmode.eq.'not')) then

! is there a void data filename? If so, then read it  
   if (voidname.ne.'none') call read_void_data(numvoids,voidname,DF_L,DF_npix,DF_npiy,dinfo)

! read namelist files for all dislocations, if any
   if (numdisl.gt.0) call read_dislocation_data(dislname,numdisl,numsf,DF_npix,DF_npiy,DF_gf,DF_L,dinfo)

! read namelist files for all stacking faults, if any
   if (numsf.gt.0) call read_stacking_fault_data(numsf,numdisl,sfname,DF_L,DF_npix,DF_npiy,DF_g,dinfo)

! is there an inclusion data file? if so, then read it
   if (incname.ne.'none') call read_inclusion_data(numinc,incname,DF_L,DF_npix,DF_npiy,dinfo)

! transform the g-vector to the defect reference frames (needed for all dislocations in CalcR).
! this can only be done AFTER all dislocations and stacking faults have been created. [converted to quaternions 6/12/13]
   do i=1,numdisl
     DF_gd(0:2,i) = quat_rotate_vector(DL(i)%a_dc,dble(DF_gc))
   end do
   
! precompute ALL the defect columns and, if needed, store them in dispfile
! this portion should be carried out in multi-threaded mode as much as possible
  allocate(disparray(DF_nums,DF_npix,DF_npiy))
  disparray = 0
  mess = ' Starting Displacement Field Computation (multi-threaded)'; call Message("(A/)")

  call CalcRLocal(numvoids,numdisl,numsf,numinc,DF_nums,DF_npix,DF_npiy,t_interval,disparray)
  
! and, if needed, store the defect displacement field for re-runs
if (dispmode.ne.'not') then
  if ((dispfile.ne.'none').and.(dispmode.eq.'new')) then 
    mess = 'Displacement field data stored in file '//dispfile; call Message("(/A/)")
    open(unit=dataunit,file=dispfile,status='new',action='write',form='unformatted')
    write (dataunit) DF_nums,DF_npix,DF_npiy
    write (dataunit) disparray
    call SafeCloseFile('d1','keep',dispfile,.FALSE.)
  endif
  if ((dispfile.ne.'none').and.(dispmode.eq.'old')) then ! there is a pre-computed defect file, so let's load it
   allocate(disparray(DF_nums,DF_npix,DF_npiy))
   disparray = 0
   open(unit=dataunit,file=dispfile,status='old',action='read',form='unformatted')
   read (dataunit) DF_nums_new,DF_npix_new,DF_npiy_new
! check to make sure that these dimensions are the same as the ones used in the current run of the program
   if ((DF_nums_new.ne.DF_nums).or.(DF_npix_new.ne.DF_npix).or.(DF_npiy_new.ne.DF_npiy)) then
    io_int(1) = DF_nums_new; io_int(2) = DF_npix_new; io_int(3) = DF_npiy_new
    call WriteValue('The dimensions of the defect array in the file are : ', io_int, 3, "(3I5)")
    io_int(1) = DF_nums; io_int(2) = DF_npix; io_int(3) = DF_npiy
    call WriteValue('The dimensions in the SRdef_rundata file do not agree : ', io_int, 3, "(3I5)")
    mess = 'Terminating program run'; call Message("(A)")
    stop
  end if
! ok, we're good, so read the actual data...  
  read (dataunit) disparray
  close(unit=dataunit,status='keep')
 end if
end if
end if

! next we prepare for the actual image simulation
  npix = DF_npix
  npiy = DF_npiy
  if (illumination_mode.eq.'CTEM') then
     allocate(images(nn,npix,npiy))
     images = 0.0
  else   ! in STEM mode, there are only two signals, so two images ... 
     allocate(STEMimages(2,npix,npiy,STEM%numCL))
     STEMimages = 0.0
  end if
  
! allocate and initialize auxiliary variables 
  allocate(Afirst(nn,nn),q(nn,nn),qin(nn,nn),qout(nn,nn),r(nn,nn))
  Afirst = czero;  q = czero; qin = czero;  qout = czero;  r = czero
  
! next, set up the parameters for the scattering matrix formalism
  dgr = 2.D0*cPi/float(Nmat)       	! angular increment for theta values
  nsl = 4 						! this could also be taken to be 3, but computation will be less accurate
  sl = DF_slice/2**nsl             		! take a fraction of the slice thickness

! initialize the timer
if (illumination_mode.eq.'CTEM') then
  numstart = 1
  numstop = 1
  call Time_report(0.01*t_interval)  
else
  numstart = 2
  numstop = STEM%numberofsvalues-1    ! weight factors are zero for end-members ...
  call Time_report(1.0/float(numstop-numstart+1))
end if
call Time_start

! the following should not be needed, but let's make sure it works before we change it...

! to enable threaded excution, we must first make a local copy of 
! the weights arrays STEM%BFweightsarray and STEM%ADFweightsarray,
! since OpenMP does not allow for this type of variable constructs
allocate(BFweightsarray(nn,STEM%numberofsvalues,STEM%numCL),ADFweightsarray(nn,STEM%numberofsvalues,STEM%numCL))
BFweightsarray = STEM%BFweightsarray
ADFweightsarray = STEM%ADFweightsarray
numCL = STEM%numCL


mainloop: do isg = numstart,numstop   ! this is the main computational loop

! get the correct excitation errors for this beam orientation (in STEM mode)
if (illumination_mode.eq.'STEM') then  
! fill the diagonal of the reference dynamical matrix and the void matrix
  do i=1,nn
   DHWMz(i,i)=2.0*cPi*cmplx(0.0,STEM%sgarray(i,isg))    ! initialize the diagonal element of the dynamical matrix
   DHWMvoid(i,i) = DHWMz(i,i)
  end do
end if
 
  allocate(Az(nn,nn))
! compute the first slice scattering matrix from the exponential Taylor expansion  (eq. 5.27)
  Afirst = czero
  do i=1,nn
    Afirst(i,i) = cone    ! initialize Afirst to be the identity matrix
  end do
  thr  = 1.0d-10        ! convergence threshold for Taylor series truncation
! first the void matrix (essentially a vacuum propagator)
  Az = Afirst
  q = DHWMvoid*cmplx(sl,0.D0,dbl)
  qout = q
  zmax = 1.D0
  i = 2
  do while (zmax.gt.thr)     ! loop until convergence is obtained
   Az = Az + qout
   r = q/cmplx(dble(i),0.D0,dbl)
   qin = matmul(qout,r)
   zmax = maxval(abs(qout-qin))
   qout = qin
   i = i+1
  end do                             ! OK, the thinner slice has converged
  Az = Az + qout
! now multiply Az with itself nsl times to get to thickness DF_slice
  do i=1,nsl
    Az = matmul(Az,Az)
  end do
  DF_Svoid = Az  ! this is the vacuum propagator
! 
! main loop for the array of Nmat scattering matrices
  numi = 0
  do k=0,Nmat-1  
! initialize theta array
    do i=-nn,nn 
      arg = float(i)*float(k)*dgr
      theta(i) = cmplx(cos(arg),-sin(arg),dbl)
    end do
! then multiply DHWMz by appropriate theta values  
    do i=-grange,grange
      do j=-grange,grange
        DHWM(i+grange+1,j+grange+1) = DHWMz(i+grange+1,j+grange+1)*theta(i-j)
      end do
    end do
! and multiply with itself to get scattering matrix.
! use the series expansion for the scattering matrix (eq. 5.27)
    Az = Afirst
    q = DHWM*cmplx(sl,0.D0,dbl)
    qout = q
    zmax = 1.D0
    i = 2
    do while (zmax.gt.thr)
     Az = Az + qout
     r = q/cmplx(dble(i),0.D0,dbl)
     qin = matmul(qout,r)
     zmax = maxval(abs(qout-qin))
     qout = qin
     i = i+1
    end do
    Az = Az + qout
    numi = numi + i
! multiply Az with itself nsl times to get to thickness DF_slice
    do i=1,nsl
      Az = matmul(Az,Az)
    end do
! and store in the main array
      DF_Sarray(k,1:nn,1:nn) = Az(1:nn,1:nn)
   end do  ! main loop
   deallocate(Az)
   mess = ' Scattering matrices computed'; call Message("(/A)")

!------------------------------------------------!
! Finally, here it is: the actual (threaded!) image computation  !
!------------------------------------------------!
! NTHR = OMP_GET_NUM_THREADS()
!  call OMP_SET_NUM_THREADS(4)
!$OMP     PARALLEL PRIVATE(TID,i,j,k,ii,jj,iCL,amp,amp2,Az,inten,weights,jcnt) &
!$OMP&   SHARED(NTHR,npix,npiy,DF_nums,disparray,DF_Sarray,DF_Svoid,illumination_mode, &
!$OMP&   BFweightsarray,ADFweightsarray,att,images,STEMimages,t_interval,nn,numCL,izero)
!$OMP     

  NTHR = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM() 
  jcnt = 0

!$OMP DO SCHEDULE (GUIDED)
donpix: do i=1,npix
  donpiy:   do j=1,npiy
      allocate(Az(nn,nn),amp(nn),amp2(nn),weights(nn),inten(nn))   ! these are private variables, so each thread must allocate them !

! initialize the wave function for this pixel with (1.0,0.0) for the incident beam izero
      amp = cmplx(0.0,0.0,dbl)
      amp(izero) = cmplx(1.0,0.0,dbl)

! loop over the fixed thickness slices
      doslices: do k=1,DF_nums

! select the appropriate scattering matrix to propagate with (see section 8.3.3 in the book)
       if (disparray(k,i,j).eq.-10000) then  ! this is point inside a void
 	 Az = DF_Svoid    ! so we use the void propagator matrix
       else  ! it is not a void
	 Az = DF_Sarray(disparray(k,i,j),1:nn,1:nn)
       end if
! and multiply with this matrix
       amp2 = cmplx(0.0,0.0,dbl)
       do ii=1,nn
        do jj=1,nn
         amp2(ii) = amp2(ii) + Az(ii,jj) * amp(jj)
        end do
       end do
       amp = amp2
      end do doslices ! loop over slices
      
! compute the intensities and (for STEM mode) use the proper weight-factor
      if (illumination_mode.eq.'STEM') then
        do iCL=1,numCL
! BF detector
          inten(1:nn) = att*abs(amp(1:nn))**2
          weights(1:nn) = BFweightsarray(1:nn,isg,iCL)
          STEMimages(1,i,j,iCL) = STEMimages(1,i,j,iCL) + sum(inten*weights)
          weights(1:nn) = ADFweightsarray(1:nn,isg,iCL)
          STEMimages(2,i,j,iCL) = STEMimages(2,i,j,iCL) + sum(inten*weights)
        end do
      else
! compute the (attenuated) intensities for the CTEM images and store
        images(1:nn,i,j) = att*abs(amp(1:nn))**2
      end if
      deallocate(Az,amp,amp2,weights,inten)
    end do donpiy
    if ((TID.eq.0).and.(illumination_mode.eq.'CTEM')) then 
      jcnt = jcnt+1
      if (mod(jcnt,t_interval).eq.0) call Time_remaining(jcnt,npix/NTHR)
    endif
  end do donpix
!$OMP END DO 
!$OMP END PARALLEL
  
  
  if (illumination_mode.eq.'STEM') call Time_remaining(isg-numstart+1,numstop-numstart+1)
200 end do mainloop

    
! ok, so the main computation is complete; print some timing information
call Time_stop(npix*npiy)


if (outputformat.ne.'data') then
! we must normalize the image series to the intensity range from 0 to 255,
! making sure that individual images have the same relative scale
  if (illumination_mode.eq.'STEM') then
    mi = minval(STEMimages)
    ma = maxval(STEMimages)
    io_real(1)=mi; io_real(2)=ma;
    call WriteValue('Minimum and maximum image/detector intensities ', io_real, 2, "(F8.5,',',F8.5)")
    STEMimages = 255.0*(STEMimages-mi)/(ma-mi)
  else
    mi = minval(images)
    ma = maxval(images)
    io_real(1)=mi; io_real(2)=ma;
    call WriteValue('Minimum and maximum image/detector intensities ', io_real, 2, "(F8.5,',',F8.5)")
    images = 255.0*(images-mi)/(ma-mi)
  end if
! the fnumber variable is a temporary construction, and should be replaced by some better string handling
    fnumber = (/ '_-10g','_-09g','_-08g','_-07g','_-06g','_-05g','_-04g','_-03g','_-02g','_-01g','_-00g', &
                         '_+01g','_+02g','_+03g','_+04g','_+05g','_+06g','_+07g','_+08g','_+09g','_+10g' /)
end if


! and save the data in whatever the selected format is
mess(1:80) = ' '
select case (outputformat) 
case ('data')
      fname = trim(trim(outputroot)//'.data')
      open(unit=dataunit,file=trim(fname),status='unknown',action='write',form='unformatted')
      if (illumination_mode.eq.'CTEM') then
        write (dataunit) nn,npix,npiy
        write (dataunit) images
      else
        write (dataunit) 2,npix,npiy,numCL
        write (dataunit) STEMimages
      end if
      close(unit=dataunit,status='keep')
case ('tiff')
      TIFF_nx = npix
      TIFF_ny = npiy
      allocate(TIFF_Image(0:TIFF_nx-1,0:TIFF_ny-1))
      if (illumination_mode.eq.'CTEM') then
        do i=outputfirst,outputlast
          do j=1,npix
            TIFF_Image(j-1,0:npiy-1) = int(images(i,j,1:npiy))
          end do
          TIFF_filename = trim(outputroot)//fnumber(i-izero)//'.tiff'
          call TIFF_Write_File
        end do
      else
        do iCL=1,numCL
          do j=1,npix
            TIFF_Image(j-1,0:npiy-1) = int(STEMimages(1,j,1:npiy,iCL))
          end do
          write (ci,*) iCL
          TIFF_filename = trim(outputroot)//'_'//trim(adjustl(ci))//'_STEM_BF.tiff'
          call TIFF_Write_File
          do j=1,npix
            TIFF_Image(j-1,0:npiy-1) = int(STEMimages(2,j,1:npiy,iCL))
          end do
          TIFF_filename = trim(outputroot)//'_'//trim(adjustl(ci))//'_STEM_ADF.tiff'
          call TIFF_Write_File
        end do
      end if  
     deallocate(TIFF_Image)
case ('pgm')
      allocate(PGM_Image(npix,npiy))
      PGM_nx = npix
      PGM_ny = npiy
      if (illumination_mode.eq.'CTEM') then
        do i=outputfirst,outputlast
          do j=1,npix
            PGM_Image(j,1:npiy) = int(images(i,j,1:npiy))
          end do
          imname(1:100) = ' '
          imname =  trim(outputroot)//fnumber(i-izero)//'.pgm'
          call PGM_Write_File(imname)
        end do
      else
        do iCL=1,numCL
          do j=1,npix
            PGM_Image(j,1:npiy) = int(STEMimages(1,j,1:npiy,iCL))
          end do
          write (ci,*) iCL
          imname(1:100) = ' '
          imname =  trim(outputroot)//'_'//trim(adjustl(ci))//'_STEM_BF.pgm'
          call PGM_Write_File(imname)
          do j=1,npix
            PGM_Image(j,1:npiy) = int(STEMimages(2,j,1:npiy,iCL))
          end do
          imname(1:100) = ' '
          imname =  trim(outputroot)//'_'//trim(adjustl(ci))//'_STEM_ADF.pgm'
          call PGM_Write_File(imname)
        end do
      end if
      deallocate(PGM_Image)
end select  
  
end subroutine ComputeSRdefect


! we should not need to do this... all it takes is to change the settings in the defectmodule.f90 file 
! to enable OpenMP runs ... 


subroutine CalcRLocal(numvoids,numdisl,numsf,numinc,lDFnums,lDFnpix,lDFnpiy,t_interval,disparray)

! this routine returns the total displacement field (multithreaded with OPENMP)

use local
use constants
use crystal
use crystalvars
use dislocation
use foilmodule
use void
use stacking_fault
use inclusion
use defectmodule
use timing
use quaternions
use omp_lib


IMPLICIT NONE

integer(kind=irg)     :: i,j,k,ii,islice,numvoids,numdisl,numsf,numinc,lDFnums,lDFnpix,lDFnpiy,TID,NTHR,imat,t_interval
integer(kind=irg)     :: disparray(lDFnums,lDFnpix,lDFnpiy),jcnt
real(kind=dbl)        :: rx,ry,rz,dis,xpos,ypos,zpos,RR(3),sumR(3),thick,tmp(3),tmpf(3),u(3),zaamp,zaphase,zar,zai,zr(3),zi(3), &
                                 zt,fx,fy,fz,z0,gdotR,nunit(3)    !,&
!                                nu,x,y,z,zn,t,pre,r1,r2,r3,th,rn 
type(quaterniond)	:: afi, afc                         
real(kind=dbl)            :: lDFR(3)
complex(kind=dbl)     :: za(3)
complex(kind=sgl)     :: zero
logical               :: void
type (voidtype), allocatable  :: lvoids(:)
real(kind=sgl),allocatable :: lsg(:,:)
type (dislocationtype), allocatable  :: lDL(:)    
type (stackingfaulttype), allocatable  :: lSF(:)
type (inclusiontype), allocatable  :: linclusions(:)

! before we start the threads, we need to copy data from various modules
! into local variables that can then be accessed by the threads ...
 Nmat = 10000
 disparray = 0
 
! foil unit normal in microscope frame 
 nunit = quat_rotate_vector(conjg(foil%a_mc),dble(foil%Fn)) ! matmul(foil%Fn,transpose(foil%a_mc))
! foil normal parameters for zpos computation
 fx = -nunit(1)/nunit(3)
 fy = -nunit(2)/nunit(3)
 fz = foil%zb*0.5
! other parameters
 z0 = foil%z0
 thick = foil%zb
 zero = cmplx(0.0,0.0)
 afi = foil%a_fi
 afc = foil%a_fc
 if (allocated(voids)) then
   allocate(lvoids(numvoids))
   lvoids = voids
 endif
 allocate(lsg(foil%npix,foil%npiy))
 lsg = foil%sg
 if (allocated(DL)) then 
   allocate(lDL(numdisl+2*numsf))
   lDL = DL
 endif
 if (allocated(SF)) then 
   allocate(lSF(numsf))
   do ii=1,numsf 
     allocate(lSF(ii)%zpos(foil%npix,foil%npiy))
   end do
   lSF = SF
 endif
 if (allocated(inclusions)) then
   allocate(linclusions(numinc))
   linclusions = inclusions
 end if
 
! ok, we've copied all the necessary variables into local structures
! now we can perform the multi-threaded loop
call OMP_SET_NUM_THREADS(20)

! initiate multi-threaded segment
!$OMP   PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(TID,lDFR,gdotR,i,j,k,imat,zt,xpos,ypos,zpos,islice,dis,sumR,tmp,tmpf,ii,void,za,zar,zai,zaamp,zaphase,zr,zi,u,jcnt)
  NTHR = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

  if (TID.eq.0) then
! do time reporting only in the master thread
    write (*,*) 'Message from master thread ',TID,': splitting into ',NTHR,' threads '
    call Time_reset
    call Time_report(0.01*t_interval)
    call Time_start
    jcnt = 0
  end if

!$OMP barrier
!$OMP DO SCHEDULE (GUIDED)
  do i=1,lDFnpix  
    do j=1,lDFnpiy
! compute the displacement vectors lDFR for all points in this column
     
! scale the image coordinates with respect to the origin at the center of the image;
! this is where we need to include the zoom factor ...
      xpos = float(i-lDFnpix/2) * DF_L
      ypos = float(j-lDFnpiy/2) * DF_L
      zt =  (xpos*fx+ypos*fy+fz)
      
! loop over all slices (this is the main loop)
 sliceloop: do islice = 1,lDFnums 
      lDFR = 0.0

! zpos is the position down the column, starting at zt (in image coordinates)
      zpos = zt - float(islice)*DF_slice
    
! set the displacements to zero
       sumR = 0.0
    
! convert image point (xpos,ypos,zpos) to tmpf in the foil reference frame
!       tmpf = matmul(matmul( (/ xpos, ypos, zpos /),transpose(afi)),afc)
       tmpf = quat_rotate_vector( afi, dble( (/ xpos, ypos, zpos /) ))
    
! voids are easy to deal with; we simply return -10000 for each point tmpf that lies inside
! one of the voids; the calling routine then knows to use the void scattering matrix.
       if (numvoids.ne.0) then 
! are we inside a void ?
           void = .FALSE.
           voidloop: do ii=1,numvoids
! subtract the void position from the current slice position to get the relative position vector
            tmp = tmpf -  (/ lvoids(ii)%xpos, lvoids(ii)%ypos, lvoids(ii)%zpos /)
            dis = CalcLength(tmp,'c')
            if (dis.lt.lvoids(ii)%radius) then ! inside void
              void = .TRUE.
              exit voidloop
            end if
           end do voidloop
! skip the rest of the computation for this slice if we are inside a void
           if (void.eqv..TRUE.) then 
             lDFR(1) = -10000.0
             cycle sliceloop
           end if
        end if 

! ok, so we're not inside a void...
! first we take the foil shape into account using equations (8.28) and (8.29)
!      sumR = sumR + float(islice)*DF_slice*lsg(i,j)*DF_gstar

! let's put a few dislocations in ... (see section 8.4.2)

       do ii=1,numdisl
! convert the defect location from untilted image space to the tilted foil reference frame, and subtract it from the current
! column and slice position
         tmp =  (/ xpos, ypos, zpos /) - &
              quat_rotate_vector(conjg(afi), dble( (/ DF_L*lDL(ii)%id, DF_L*lDL(ii)%jd, lDL(ii)%zfrac*z0 /) ) )

! then convert the difference vector to the defect reference frame for this dislocation (we will only need the x and y coordinates)
         tmp = quat_rotate_vector(lDL(ii)%a_id,dble(tmp))
         
! check the z-coordinate; if it falls beyond the dislocation line that is inside the foil, then skip
! the displacement computation... the top and bottom coordinates of the dislocation intersections
! measured along the dislocation line were pre-computed when the dislocations were first read from
! the namelist files...
!         if (abs(tmp(3)).le.lDL(ii)%zu) then
! compute x1 + p_alpha x2  (eq. 8.38)
          za(1:3) = tmp(1) + lDL(ii)%pa(1:3)*tmp(2)
! compute the displacement vector u (eq. 8.38) [this expands the log of a complex number and takes the real part only] 
          if (tmp(1).gt.0.0) then
           do k=1,3
            zar =  real(za(k))
            zai = aimag(za(k))
            zaamp = abs(za(k))
            zaphase = abs(zai/zar)
            zr(k) = log(zaamp)
            zi(k) = atan(zaphase)
            if (zar.le.0.0) then
              if (zai.lt.0.0) zi(k) = -cPi+zi(k)
              if (zai.eq.0.0) zi(k) = cPi
              if (zai.gt.0.0) zi(k) = cPi-zi(k)
            else
              if (zai.lt.0.0) zi(k) = -zi(k)
            end if
           end do
          else
           do k=1,3
            zar =  real(za(k))
            zai = aimag(za(k))
            zaamp = abs(za(k))
            zaphase = abs(zai/zar)
            zr(k) = log(zaamp)
            zi(k) = atan(zaphase)
            if (zar.le.0.0) then
              if (zai.gt.0.0) zi(k) = cPi-zi(k)
              if (zai.eq.0.0) zi(k) = cPi
              if (zai.lt.0.0) zi(k) = cPi+zi(k)
            else
              if (zai.lt.0.0) zi(k) = 2.0*cPi-zi(k)
              if (zai.eq.0.0) zi(k) = 0.0
            end if
           end do  
          end if
          u = 2.0*real(matmul(lDL(ii)%dismat,cmplx(zr,zi)))
! transform displacement vector u to the Cartesian crystal reference frame
          sumR = sumR + quat_rotate_vector(lDL(ii)%a_dc,dble(u))
!         end if 
       end do


! stacking faults (this is easy because we've already done all the work in the stacking_fault module)
       do ii=1,numsf
         if ((zpos.lt.lSF(ii)%zpos(i,j)).and.(lSF(ii)%zpos(i,j).ne.-10000.0)) then 
           sumR = sumR + lSF(ii)%lpbc
         end if
       end do


! Mader's expression for the displacement field of a large inclusion
!   if (0.eq.1.) then 
!    nu = 0.25
!    ce = 0.005
!    rn = 25.0*DF_L
!    x = (float(i-DF_npix/2)-0.5)*DF_L
!    y = (float(j-DF_npiy/2)-0.5)*DF_L
!    z = float(k)*DF_slice
!    zn = 100.5*DF_slice
!    t = DF_slice * DF_nums
!    pre = (1.0+nu)/(3.0*(1.0-nu))*ce*rn**3
! 
!    r1 = sqrt(x**2+y**2+(z-zn)**2)
!    r2 = sqrt(x**2+y**2+(z+zn)**2)
!    r3 = sqrt(x**2+y**2+(2.0*t-z-zn)**2)
! 
!    if (((r1.eq.0.0).or.(r2.eq.0.0)).or.(r3.eq.0.0)) then
!      return
!    else
!     dis = (1.0/r1**3+(3.0-4.0*nu)/r2**3-6.0*z*(z+zn)/r2**5+(3.0-4.0*nu)/r3**3-6.0*(t-z)*(2.0*t-z-zn)/r3**5)
!     rx = x*dis
!     ry = y*dis
!     rz = (z-zn)/r1**3-(3.0-4.0*nu)*((z+zn)/r2**3+(2.0*t-z-zn)/r3**3)-6.0*z*(z+zn)**2/r2**5 + &
!          2.0*z/r2**3+6.0*(t-z)*(2.0*t-z-zn)**2/r3**5-2.0*(t-z)/r3**3
! 
!     sumR = pre*(/ rx, ry, rz /)
!     return
!    end if
!   end if

! then the coherent precipitates, using the model in section 8.4.1
        if (numinc.gt.0) then
         do ii=1,numinc
! subtract the inclusion position from the current slice position to get the relative position vector
           tmp = tmpf -  (/ linclusions(ii)%xpos, linclusions(ii)%ypos, linclusions(ii)%zpos /)
           dis = CalcLength(tmp,'c')
           if (dis.ge.linclusions(ii)%radius) then ! outside particle
             tmp = tmp*(linclusions(ii)%radius/dis)**3
           end if
           sumR = sumR + linclusions(ii)%C*tmp
         end do
        end if
 
! finally any displacement fields defined by the user routine UserDisp
! sumR = sumR + UserDisp()

! then convert to the dot-product 
! and select the appropriate scattering matrix to propagate with (see section 8.3.3 in the book)
       if (sumR(1).eq.-10000.0) then  ! this is point inside a void
 	disparray(islice,i,j) = -10000
       else  ! it is not a void, so use the full dot product g.R (both vectors must be Cartesian !)
         gdotR = Dot_Product(DF_gc,sumR)
         imat = nint(float(Nmat)*amod(sngl(gdotR+10000.D0),1.0))  ! select which pre-computed scattering matrix to use for this slice
	 if (imat.lt.0) imat = 0                          ! this shouldn't happen, but we'll check just in case ...
	 if (imat.ge.Nmat) imat = Nmat-1                  ! neither should this one ...
       end if
       disparray(islice,i,j) = imat

    end do sliceloop ! main loop over the slices

   end do  ! j-loop
  if (TID.eq.0) then 
    jcnt = jcnt+1
    if (mod(jcnt,t_interval).eq.0) call Time_remaining(jcnt,lDFnpix/NTHR)
  endif
end do  ! i-loop
!$OMP END DO 
 if (TID.eq.0) call Time_stop(DF_npix*DF_npiy)
!$OMP END PARALLEL
  

end subroutine CalcRLocal
       
