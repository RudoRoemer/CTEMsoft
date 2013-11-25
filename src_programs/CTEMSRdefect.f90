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
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Computes systematic row defect contrast for multiple defects
!                      using CTEM or STEM illumination and detector conditions 
!
!> @todo OpenMP implementation; overall verification of all coordinate transformations
! 
!> @date  01/14/10 MDG 1.0 original, based on SRdefect.f90 (collaboration with OSU) 
!> @date  05/14/13 MDG 2.0 replaced all IO by namelist file and added command line argument handling
!--------------------------------------------------------------------------
program CTEMSRdefect 
! 
use local
use files
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMSRdefect.nml'
progname = 'CTEMSRdefect.f90'
call Interpret_Program_Arguments(nmldeffile,7,(/ 2, 3, 4, 200, 201, 202, 203 /) )

! initialize all user-defined variables
call ComputeSRdefect(nmldeffile)

end program CTEMSRdefect


!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeSRdefect
!
!> @author Marc De Graef, Carnegie Mellon University
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
!> @date 11/22/13  MDG 3.0 added general output for STEM IDL routine, renamed old STEM mode to BFDF
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
use math
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

character(fnlen),INTENT(IN)      :: nmlfile 

integer(kind=irg)                :: ira,nn,izero,i,j,k,n,nsl,numi,npix,npiy,ii,jj,numvoids,numdisl,numYdisl,numsf, skip, &
                                  numinc,dinfo,t_interval,outputfirst,outputlast,DF_nums_new, io_int(3), grange, &
                                  DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, NTHR, jcnt, numCL, iCL, &
                                  SR_g(3), SETNTHR, imat
!                                  OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
real(kind=sgl)                   :: ind(3),hkl(3),exerg,cosom,glen,exer,dgr,sl,thr,arg,zmax,thick,mi,ma, &
                                  att,xgp,DF_gf(3),DF_gd(3,maxdefects), io_real(1), voltage, kt, frac, gdotR
character(fnlen)                 :: dataname,sgname,voidname,dislname(3*maxdefects),sfname(maxdefects), &
				     outputroot,incname,dispfile,tempname,xtalname,STEMname, foilnmlfile, STEMnmlfile
character(4)                     :: outputformat, illumination_mode, dispmode, progmode
character(5)                     :: fnumber(-10:10)
character(len=10)                :: ci

complex(kind=dbl),allocatable    :: DHWM(:,:),Afirst(:,:),DHWMvoid(:,:)
complex(kind=dbl),allocatable    :: q(:,:),qin(:,:),qout(:,:),r(:,:),amp(:),amp2(:)
complex(kind=dbl)                :: czero=cmplx(0.0,0.0,dbl),cone=cmplx(1.0,0.0,dbl)
real(kind=sgl),allocatable       :: weights(:), inten(:), STEMimages(:,:,:,:)
integer(kind=irg),allocatable    :: disparray(:,:,:)
real(kind=sgl),allocatable       :: BFweightsarray(:,:,:),ADFweightsarray(:,:,:)


namelist / SRdeflist / DF_L, DF_npix, DF_npiy, DF_slice, Nmat, sgname, numvoids, incname, &
                     voidname, numdisl, dislname, numsf, sfname, dinfo, outputformat, &
                     dataname,t_interval,dispfile, SETNTHR, &
		      dispmode, xtalname, voltage, SR_g, grange, STEMnmlfile, foilnmlfile

! here we read the general simulation information from a namelist file SRdef_rundata.nml
! first we define the default values

! parameters specific to this run
 xtalname = 'undefined'		! initial value; MUST be present in nml file for program to execute
 voltage = 200000.0		        ! accelerating voltage
 SR_g = (/ 1,0,0 /)			! systematic row g-vector
 grange = 5				! maximum positive multiple of g-vector; total number will be 2*grange+1 

! input files
 STEMnmlfile = 'STEM_rundata.nml'	! name of the STEM rundata namelist file
 foilnmlfile = 'SRdef_foildata.nml'	! name of the foil rundata namelist file
 
! column approximation parameters and image parameters 
 DF_L = 1.0             		! edge length of column in nanometers
 DF_npix = 256       			! number of image pixels along x
 DF_npiy = 256       			! number of image pixels along y 
 DF_slice = 1.0       			! slice thickness in nanometers

 dinfo = 0               		! switch to make makedislocation verbose
 sgname = 'nofile'   			! if this variable is different from 'nofile', then an external sg array is read (to be implemented)

! defect parameters
 numdisl = 0           		! number of dislocation files
 numsf = 0             		! number of stacking fault files
 numinc = 0           			! number of inclusions
 numvoids = 0       			! number of voids
 voidname = 'none' 			! filename for void data
 dislname = ''         		! filenames for dislocation data
 sfname = ''            		! filenames for stacking fault data
 incname = 'none'   			! filename for inclusion data
 dispfile = 'none'     		! name of the displacement field output file (will be created if different from none)
 dispmode = 'not'  			! should a diplacement file be written ('new') or read ('old') or neither ('not')?

! output parameters
 dataname = 'SRdefect.data'           ! default output file name
 t_interval = 10       		! default timing interval (output every t_interval image columns)
 SETNTHR = 6                          ! number of threads to use (OpenMP)

! other parameters that the user does not have access to
 Nmat = 3600       			! number of precomputed A matrices to be stored (every 0.1 degrees)
 progmode = 'STEM'
 frac = 0.05
 numYdisl = 0
 
! then we read the rundata namelist, which may override some of these defaults  
 OPEN(UNIT=dataunit,FILE=nmlfile,DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=SRdeflist)
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

! Weickenmeier-Kohl scattering parameters with absorption form factors
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

 ! generate all atom positions
 call CalcPositions('v')
 
! use systematic row vector to compute G* as in eq. 8.28.
  DF_gf = float(SR_g)
  DF_gstar = DF_gf/CalcLength(DF_gf,'r')**2    ! define G* such that G.G* = 1
  call TransSpace(DF_gf,DF_gc,'r','c')         ! convert to Cartesian reference frame

! compute total number of beams
  nn = 2*grange+1  ! total number of beams
  izero = (nn+1)/2 ! id number of g=0 beam 

! read the STEM parameters from a namelist file and initialize all STEM related arrays
  call read_STEM_data(STEMnmlfile, nn, SR_g, kt, i)  
  mess = 'Read STEM namelist file '//STEMnmlfile; call Message("(A)")

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
     call CalcUcg(int(hkl))                     ! compute the interaction parameters
     DHWMz(i,j) = cPi*cmplx(-aimag(rlp%qg),real(rlp%qg),dbl)  ! and initalize the off-diagonal matrix element
    else
     DHWMz(i,j) = czero                         ! for now at least; this will be filled in later
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
  
! next, we read the foildata namelist from the foil namelist file
! this includes material property data, in this case the elastic moduli
  call read_foil_data(foilnmlfile,DF_npix,DF_npiy,DF_L,dinfo)
    
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
   
! precompute ALL the defect columns and, if needed, store them in dispfile
! this portion should be carried out in multi-threaded mode as much as possible
  allocate(disparray(DF_nums,DF_npix,DF_npiy), DF_R(DF_nums,3))
  disparray = 0
  mess = ' Starting Displacement Field Computation (multi-threaded)'; call Message("(A/)")

  call CalcRLocal(numvoids,numdisl,numsf,numinc,DF_nums,DF_npix,DF_npiy,t_interval,disparray,SETNTHR)

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

! next the STEMimages array 
  npix = DF_npix
  npiy = DF_npiy
  allocate(STEMimages(nn,npix,npiy,STEM%numberofsvalues))
  STEMimages = 0.0

  
! allocate and initialize auxiliary variables 
  allocate(Afirst(nn,nn))
  Afirst = czero
  
! initialize the timer
  numstart = 1
  numstop = STEM%numberofsvalues    
  call Time_report(frac)
  call Time_start

! the following should not be needed, but let's make sure it works before we change it...

! before we start the main computational loop, we need to store the 
! necessary variables for the reconstruction of STEM BF/HAADF images
! using the post-processing IDL routine.  This file should contain all
! the information needed to recreate the full CBED patterns at each image
! pixel.  This is then followed by the actual data.

  mess ='Storing data for IDL visualization program in '//dataname
  call Message("(A/)")
  open(UNIT=dataunit,FILE=trim(dataname),STATUS='unknown',FORM='unformatted')  
  
! program mode
  write (dataunit) progmode

! filename of corresponding data set
  write (dataunit) dataname


write (*,*) 'foil normal test = ',foil%F

! various bits of useful information
  write (dataunit) xtalname						! crystal structure file names	
  write (dataunit) foil%F						! incident wave vector
  write (dataunit) CalcLength(float(SR_g),'r')                       ! reciprocal lattice vector length
  write (dataunit) CalcDiffAngle(SR_g(1),SR_g(2),SR_g(3))*0.5	! Bragg angle for the first reflection (whether allowed or not)
  write (dataunit) STEM%numberofsvalues				! number of pixels along disk radius
  write (dataunit) sngl(mLambda)					! wave length
  write (dataunit) STEM%beamconvergence				! beam divergence
  write (dataunit) DF_L						! pixel size
  write (dataunit) nn   ! number of reflections
! data array size
  write (dataunit) nn,npix,npiy,STEM%numberofsvalues


dgr = 2.D0*cPi/float(Nmat)       	! angular increment for theta values

mainloop: do isg = numstart,numstop   ! this is the main computational loop

! get the correct excitation errors for this beam orientation
! fill the diagonal of the reference dynamical matrix and the void matrix
  do i=1,nn
   DHWMz(i,i)=2.0*cPi*cmplx(0.0,STEM%sgarray(i,isg))    ! initialize the diagonal element of the dynamical matrix
   DHWMvoid(i,i) = DHWMz(i,i)
  end do

! compute the first slice scattering matrix from the exponential Taylor expansion  (eq. 5.27)
! for the void (i.e., sort of a vacuum propagator)
  Afirst = czero
  do i=1,nn
    Afirst(i,i) = cone    ! initialize Afirst to be the identity matrix
  end do
  call MatrixExponential(Afirst, DF_Svoid, dble(DF_slice), 'Pade', nn)  

! and for the "non-void" 
  allocate(Az(nn,nn))
! main loop for the array of Nmat scattering matrices
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

    call MatrixExponential(DHWM, Az, dble(DF_slice), 'Pade', nn)  

! and store in the main array
    DF_Sarray(k,1:nn,1:nn) = Az(1:nn,1:nn)
   end do  ! main loop
   deallocate(Az)
!   mess = ' Scattering matrices computed'; call Message("(/A)")

!------------------------------------------------!
! Finally, here it is: the actual (threaded!) image computation  !
!------------------------------------------------!
! NTHR = OMP_GET_NUM_THREADS()
!  call OMP_SET_NUM_THREADS(4)
!$OMP     PARALLEL PRIVATE(TID,i,j,k,ii,jj,iCL,amp,amp2,Az,inten,weights,jcnt) &
!$OMP&   SHARED(NTHR,npix,npiy,DF_nums,disparray,DF_Sarray,DF_Svoid,illumination_mode, &
!$OMP&   BFweightsarray,ADFweightsarray,att,images,STEMimages,t_interval,nn,numCL,izero,isg)
!$OMP     

  NTHR = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM() 
  jcnt = 0
  
!  write (*,*) 'max(DF_Sarray) ' ,maxval(cabs(DF_Sarray))
!  write (*,*) 'range disparray ',minval(disparray),maxval(disparray)

!$OMP DO SCHEDULE (GUIDED)
donpix: do i=1,npix
  donpiy:   do j=1,npiy
      allocate(Az(nn,nn),amp(nn),amp2(nn),weights(nn))   ! these are private variables, so each thread must allocate them !

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
! and multiply with this matrix; for some reason, the Absoft compiler has an issue with 
! this statement, and the matmul routine can not be used... probably a compiler bug...
!       amp2 = matmul(Az,amp)
! so we do the computation explicitly
       amp2 = cmplx(0.0,0.0,dbl)
       do ii=1,nn
        do jj=1,nn
         amp2(ii) = amp2(ii) + Az(ii,jj) * amp(jj)
        end do
       end do
       amp = amp2
      end do doslices ! loop over slices
      
! compute the (attenuated) intensities for the CTEM and STEM images and store
      STEMimages(1:nn,i,j,isg) = att*cabs(amp(1:nn))**2
 
     deallocate(Az,amp,amp2,weights)
    end do donpiy
  end do donpix
!$OMP END DO 
!$OMP END PARALLEL

   if (float(isg)/float(numstop) .gt. frac) then
     call Time_remaining(isg,numstop)
     frac = frac + 0.05
   end if  
   
200 end do mainloop

    
! ok, so the main computation is complete; print some timing information
call Time_stop(npix*npiy)

! and write all the data to the output file
write (dataunit) STEMimages
close(unit=dataunit,status='keep')
mess = 'Data stored in file '//trim(dataname)
call Message("(/A/)")

end subroutine ComputeSRdefect


! we should not need to do this... all it takes is to change the settings in the defectmodule.f90 file 
! to enable OpenMP runs ... 

!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcRLocal
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief displacement field parameters for a systematic row case
!
!> @param nmlfile namelist file name
!
!> @todo We should not have to use a separate routine for this; needs to be verified against 
!> the regular CalcR ... 
!
!> @date 11/22/13  MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine CalcRLocal(numvoids,numdisl,numsf,numinc,lDFnums,lDFnpix,lDFnpiy,t_interval,disparray,SETNTHR)

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
use rotations
use omp_lib


IMPLICIT NONE

integer(kind=irg)     :: i,j,k,ii,islice,numvoids,numdisl,numsf,numinc,lDFnums,lDFnpix,lDFnpiy,TID,NTHR,imat,t_interval
integer(kind=irg)     :: disparray(lDFnums,lDFnpix,lDFnpiy),jcnt, SETNTHR
real(kind=dbl)        :: rx,ry,rz,dis,xpos,ypos,zpos,RR(3),sumR(3),thick,tmp(3),tmp2(3),tmpf(3),u(3),zaamp,zaphase,zar, &
                         zai,zr(3),zi(3),zt,fx,fy,fz,z0,gdotR,nunit(3), a_fm(3,3)    !,&
!                        nu,x,y,z,zn,t,pre,r1,r2,r3,th,rn 
real(kind=dbl)        :: afi(4), afc(4)                         
real(kind=dbl)        :: lDFR(3)
complex(kind=dbl)     :: za(3)
complex(kind=sgl)     :: zero
logical               :: void
type (voidtype), allocatable          :: lvoids(:)
real(kind=sgl),allocatable            :: lsg(:,:)
type (dislocationtype), allocatable   :: lDL(:)    
type (stackingfaulttype), allocatable :: lSF(:)
type (inclusiontype), allocatable     :: linclusions(:)

! before we start the threads, we need to copy data from various modules
! into local variables that can then be accessed by the threads ...
 Nmat = 3600
 disparray = 0
 
! foil unit normal in microscope frame 
 a_fm = qu2om(foil%a_fm)
 fx = a_fm(3,1)
 fy = a_fm(3,2)
 fz = a_fm(3,3) 

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
call OMP_SET_NUM_THREADS(SETNTHR)

! initiate multi-threaded segment
!$OMP   PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(TID,lDFR,gdotR,i,j,k,imat,zt,xpos,ypos,zpos,islice,dis,sumR,tmp,tmp2,tmpf, &
!$OMP& ii,void,za,zar,zai,zaamp,zaphase,zr,zi,u,jcnt)
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
      zt = foil%zb*0.5 - (fx*xpos + fy*ypos)/fz
      
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


!------------------
!----CURVED FOIL---  We won't allow that in this program
!------------------
! first we take the foil shape into account using equations (8.28) and (8.29)
!      sumR = sumR + float(islice)*DF_slice*lsg(i,j)*DF_gstar

!-----------------
!--DISLOCATIONS--
!-----------------
! let's put a few dislocations in ... (see section 8.4.2)

       do ii=1,numdisl
! convert the defect location from untilted image space to the tilted foil reference frame, and subtract it from the current
! column and slice position
         tmp2 =  tmpf - dble( (/ DF_L*lDL(ii)%id, DF_L*lDL(ii)%jd, lDL(ii)%zfrac*z0 /) )

! then convert the difference vector to the defect reference frame for this dislocation (we will only need the x and y coordinates)
         tmp = quat_rotate_vector(lDL(ii)%a_df, tmp2)


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
          sumR = sumR + quat_rotate_vector(conjg(lDL(ii)%a_dc),dble(u))
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
       
