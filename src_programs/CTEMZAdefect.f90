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
! CTEMsoft2013:CTEMZAdefect.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMZAdefect 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief CTEMZAdefect computes zone axis defect contrast for multiple defects
!>                      using STEM illumination and detector conditions
!
!> @details produces two types of output: regular images, or a very large file with all
!> the CBED patterns, one for each image pixel.  An IDL program can then be used to
!> interactively change the detector settings (including camera length), as well as the 
!> segment combination for a hypothetical segmented detector.
! 
!> @date  03/07/10 MDG  1.0 original, based on STEMdefect.f90 (collaboration with OSU) 
!> @date  03/21/10 MDG  2.0 modified for STEM illumination, fast version with bi-variate interpolation
!> @date  06/19/13 MDG  3.0 conversion to new libraries 
!> @date  10/28/13 MDG  3.1 added Interpret_Program_Arguments line
!--------------------------------------------------------------------------
program CTEMZAdefect 

use local
use files
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMZAdefect.nml'
progname = 'CTEMZAdefect.f90'
call Interpret_Program_Arguments(nmldeffile,8,(/ 0, 1, 2, 3, 200, 201, 202, 203 /) )

! initialize all user-defined variables
call ComputeZAdefect(nmldeffile)

end program CTEMZAdefect


!--------------------------------------------------------------------------
!
! SUBROUTINE: ComputeZAdefect
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a zone axis STEM defect data set
!
!> @param nmlfile namelist file name
!
!> @todo check sign of kt for consistency 
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 11/13/13  MDG 2.2 implementation of Pade approximation for scattering matrix computation
!--------------------------------------------------------------------------
subroutine ComputeZAdefect(nmlfile)

use local
use error
use crystalvars 
use crystal
use symmetryvars
use symmetry
use postscript
use constants
use diffraction
use dynamical
use kvectors
use gvectors
use files
use io
use math
use foilmodule
use stacking_fault
use dislocation
use void
use inclusion
use defectmodule
use rotations
use TIFF_f90
use pgm
use timing
use STEMmodule

character(fnlen),INTENT(IN)		:: nmlfile

integer(kind=irg)    		:: nn,i,j,k,npix,npiy,ii,jj,numvoids,numdisl, numYdisl, &
					numsf,nCL,numinc,dinfo,t_interval, &
					DF_nums_new,DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, &
					NTHR, isym, ir, ga(3), gb(3),kk(3),ic,g,numd,ix,iy, &
					numk,ixp,iyp,SETNTHR, io_int(6), skip, gg(3), iSTEM
!                                  OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
integer(kind=irg),parameter 		:: numdd=180
real(kind=sgl)         		:: glen,exer,arg,thick, X(2), dmin, &
					lauec(2), g3(3), gdotR,att,xgp,DF_gf(3), &
					DM(2,2), DD, H,FNr(3),ll(3),lpg(3),gplen,LC3, c(3), gx(3), gy(3), &
					sgdenom, gac(3), gbc(3),zmax, beamdiv, ktmax, io_real(2), kt, qx, qy
character(fnlen)      			:: dataname,sgname,voidname,dislname(3*maxdefects),sfname(maxdefects), &
					incname,dispfile,xtalname,foilnmlfile, STEMnmlfile
character(4)            		:: dispmode, progmode
complex(kind=dbl),allocatable    	:: DHWM(:,:),DHWMvoid(:,:),DDD(:,:),Sarray(:,:,:,:)
complex(kind=dbl),allocatable    	:: amp(:),amp2(:),Azz(:,:)
complex(kind=dbl)                	:: czero,cone
complex(kind=dbl)                	:: para(0:numdd),dx,dy,dxm,dym
real(kind=sgl),allocatable       	:: inten(:), sgarray(:,:)
real(kind=sgl),allocatable    		:: disparray(:,:,:,:),imatvals(:,:), ZAimages(:,:,:,:)
integer(kind=irg),allocatable   	:: BFweightsarray(:,:,:),ADFweightsarray(:,:,:)
integer(kind=sgl),allocatable    	:: expval(:,:,:)

namelist / rundata / DF_L, DF_npix, DF_npiy, DF_slice, dmin, sgname, numvoids, incname, stdout, &
                                voidname, numdisl, dislname, numsf, sfname, dinfo, &
				 t_interval,progmode, dispfile, &
				 dispmode,SETNTHR,xtalname,voltage,kk, lauec, &
				 dataname, foilnmlfile, STEMnmlfile


! first we define the default values
czero=cmplx(0.0,0.0,dbl)
cone=cmplx(1.0,0.0,dbl)

! parameters specific to this run
 xtalname = 'undefined'		! initial value; MUST be present in nml file for program to execute
 voltage = 200000.0			! accelerating voltage
 kk = (/ 0, 0, 1 /)			! incident wave vector in crystal components (omitting wave length)
 lauec = (/ 0.0,0.0 /)			! Laue center coordinates (used for CTEM mode)
 dmin = 0.04			       ! smallest d-spacing to include in dynamical matrix [nm]

! CTEM or STEM ?
 progmode = 'CTEM'  		        ! default illumination mode (can be 'CTEM' or 'STEM')
 STEMnmlfile = 'STEM_rundata.nml'	! name of the STEM rundata namelist file
 foilnmlfile = 'FOIL_rundata.nml'	! name of the foil rundata namelist file
 
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
 dataname = 'ZAdefect.data'		! default outputfile name
 t_interval = 10       	        ! default timing interval (output every t_interval image columns)
! 
! then we read the actual rundata namelist, which may override some or all of these defaults  
 OPEN(UNIT=dataunit,FILE=nmlfile,DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=rundata)
 CLOSE(UNIT=dataunit)
  
! make sure the xtalname variable has been properly defined
if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMZAdefect:',' structure file name is undefined in '//nmlfile)
end if

! we got this far, so display the standard program info
 progname = 'CTEMZAdefect.f90'
 progdesc = 'Dynamical zone axis STEM illumination defect image simulation'
 call CTEMsoft
 
 numd = numdd
 
! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)

! Weickenmeier-Kohl scattering parameters with absorption form factors
 skip = 3	
 call CalcWaveLength(dble(voltage),skip)

 ! generate all atom positions
 call CalcPositions('v')

! determine the point group number and get the ZAP 2-D symmetry  NEEDS TO BE MODIFIED WITH NEW ROUTINES
 j=0
 do i=1,32
  if (SGPG(i).le.cell%SYM_SGnum) j=i
 end do
 call BFsymmetry(kk,j,isym,ir)
  
! determine and display the shortest reciprocal lattice vectors for this zone
 call ShortestG(kk,ga,gb,isym)
 io_int(1:3) = kk(1:3)
 call WriteValue('', io_int, 3,  "(//,' ','[',3I2,'] has Bright Field symmetry ',$)")
 mess = PGTWD(isym)
 call Message("(A,$)")
 io_int(1) = ir
 call WriteValue(' order = ', io_int, 1, "(I4/)")
 mess = 'Reciprocal lattice vectors : '; 
 io_int(1:3) = ga(1:3)
 io_int(4:6) = gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6, "('(',3I3,') and (',3I3,')',/)")

! determine the cartesian components of ga
 DF_gf = float(ga)
 DF_gstar = DF_gf/CalcLength(DF_gf,'r')**2    ! define G* such that G.G* = 1
 call TransSpace(DF_gf,DF_gc,'r','c')                 ! convert to Cartesian reference frame
  
! we'll need to compute the list of wavevectors at this point
! read the STEM parameters from a namelist file and initialize all STEM related arrays
if (progmode.ne.'CTEM') then
  beamdiv = 0.0
  kt = 0.0
  call read_STEM_data(STEMnmlfile,'ZA',nn,ga,kt,numk,beamdiv)          ! first we need the beam divergence angle from this file ...
  ktmax = 2.0*sin(beamdiv/2000.)/mLambda/CalcLength(float(ga),'r')     ! ktmax in units of |ga|
! this next line needs to be verified !!!
  ijmax = float(STEM%numberofsvalues)**2
  call Calckvectors(dble(kk),dble(ga),dble(ktmax),STEM%numberofsvalues,STEM%numberofsvalues,numk,isym,ijmax,'Conical')    ! here we figure out how many beams there are
  call Compute_ReflectionList(dmin,kk,ga,gb,'ALL',.FALSE.,0,beamdiv)
else  ! progmode = CTEM
  ijmax = 0.0
  call Calckvectors(dble(kk),dble(ga),dble(ktmax),STEM%numberofsvalues,STEM%numberofsvalues,numk,isym,ijmax,'Conical')    ! here we figure out how many beams there are
  call Compute_ReflectionList(dmin,kk,ga,gb,'ALL',.FALSE.,0)
end if

! for now, we do not consider weak beams at all
  if (BetheParameter%cutoff.eq.0.0) call Set_Bethe_Parameters
  BetheParameter%weakcutoff = BetheParameter%cutoff

! next, we read the foildata namelist from the SRdef_foildata.nml file
! [yes, we're using the same file as for the systematic row case]
! this includes material property data, in this case the elastic moduli,
! and the foil normal, which we will need in the next step
  call read_foil_data(foilnmlfile,DF_npix,DF_npiy,DF_L,dinfo)
  DynFN = foil%F
  
! then we need to prune the reflection list to have only reflections that will actually occur in the 
! computation
  mess = ' Pruning reflection list (this takes a while ...) '
  call Message("(A)")
  call Prune_ReflectionList(numk,nbeams)
  io_int(1) = nbeams
  call WriteValue('Number of contributing beams  : ', io_int, 1, '(I)')
  nn = nbeams

  if (progmode.ne.'CTEM') then
    call read_STEM_data(STEMnmlfile,'ZA',nn,ga,kt,numk)        ! and then we go back and initialize the rest of the STEM parameters
  end if

! ideally, we should use Bethe potentials to reduce the size of the dynamical matrix;
! while the theory has been worked out to do this, it would require tremendous changes
! to the program starting about here; given the time limitations, there won't be any 
! chance to do this before the end of the year (2013), so we'll leave that for some other time...
  
! allocate and initialize DF_Sarray, theta, and DF_Svoid
  allocate(theta(-nn:nn),DF_Svoid(nn,nn))
  DF_Sarray = czero; theta = czero

! allocate the various DHW Matrices
  allocate(DHWMz(nn,nn),DHWM(nn,nn),DHWMvoid(nn,nn))
  DHWMvoid = czero; DHWMz=czero; DHWM(nn,nn)=czero
  
 ! Compute the off-diagonal part of the complex DHW matrix (factor i is included)
 ! We can precompute those because they will not change at all during the run
 !       (these lines implement the equations on page 476 of the CTEM book)
 ! In this program, the reflections are stored using linked lists, which does not lend itself to
 ! OpenMP acceleration;  so, we'll have to re-write this at some point in the future...
 !
 ! this is also where we compute the decomposition of the reflection indices w.r.t. ga and gb,
 ! and store them in the nab(1:2) field of the linked list; these are used to compute the 
 ! defect contributions to the dynamical matrix in the displace routine of the MEmath module
 !
 DM(1,1) = CalcDot(float(gb),float(gb),'c')
 DM(1,2) = -CalcDot(float(ga),float(gb),'c')
 DM(2,1) = DM(1,2)
 DM(2,2) = CalcDot(float(ga),float(ga),'c')
 DD = DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1)

 rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
  do ir=1,nn
   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
   do ic=1,nn
    if (ic.ne.ir) then  ! exclude the diagonal
! compute Fourier coefficient of electrostatic lattice potential 
     call CalcUcg(rltmpa%hkl - rltmpb%hkl)
     DHWMz(ir,ic) = cPi*cmplx(-aimag(rlp%qg),real(rlp%qg),dbl)  ! and initialize the off-diagonal matrix element (including i)
    end if
    rltmpb => rltmpb%next  ! move to next column-entry
   end do
! decompose this point w.r.t ga and gb
   X(1) = CalcDot(float(rltmpa%hkl),float(ga),'c')
   X(2) = CalcDot(float(rltmpa%hkl),float(gb),'c')
   X = matmul(DM,X)/DD
   rltmpa%nab(1:2) = int(X(1:2))
   rltmpa => rltmpa%next   ! move to next row-entry
  end do
  mess = 'Reference Darwin-Howie-Whelan matrix initialized'; call Message("(A/)")

! compute the normal absorption factor xgp
  call CalcUcg((/0,0,0/))
  xgp = aimag(rlp%qg)
  io_real(1) = 1.0/xgp
  call WriteValue('Normal absorption length : ', io_real, 1, "(F10.5/)")

! Set up the excitation errors for CTEM illumination mode;
! distance between consecutive HOLZ layers in nm-1
if (progmode.eq.'CTEM') then  
  H = 1.0/CalcLength(float(kk),'d')
! g3 basis vector, properly scaled
  call CalcCross(float(ga),float(gb),g3,'r','r',1)
  call NormVec(g3,'r')
  g3 = H * g3
! unit foil normal in reciprocal space  
  call TransSpace(sngl(foil%F),FNr,'d','r')
  call NormVec(FNr,'r')
 
! parallel illumination, so the excitation errors need to be computed only once
! fill the diagonal of the reference dynamical matrix and the void matrix
  rltmpa => reflist%next   ! point to the front of the list
! ir is the row index
  do ir=1,nn
   glen = CalcLength(float(rltmpa%hkl),'r')
   if (glen.eq.0.0) then
    DHWMz(ir,ir) = czero
   else  ! compute the excitation error
    ll = lauec(1)*ga + lauec(2)*gb   ! Laue center vector
    lpg = ll + rltmpa%hkl                ! Laue + g
    gplen = CalcLength(lpg,'r')
    LC3 = sqrt(1.0-mLambda**2*(CalcLength(ll,'r')**2))   ! to ensure proper normalization of wave vector
    if (gplen.eq.0.0) then
           exer = -mLambda*CalcDot(float(rltmpa%hkl),ll+lpg,'r')/2.0*LC3*cos(CalcAngle(dble(kk),foil%F,'d'))	 
    else
	   sgdenom = 2.0*LC3*cos(CalcAngle(dble(kk),foil%F,'d'))-2.0*mLambda*gplen*cos(CalcAngle(lpg,FNr,'r'))
           exer = -(mLambda*CalcDot(float(rltmpa%hkl),ll+lpg,'r')-2.0*LC3*CalcDot(g3,lpg,'r'))/sgdenom
    end if
    rltmpa%sg = exer
    DHWMz(ir,ir) = cmplx(0.0,2.D0*cPi*exer,dbl)
    call CalcUcg(rltmpa%hkl)
   end if 
   DHWMvoid(ir,ir) = DHWMz(ir,ir)
   rltmpa => rltmpa%next   ! move to next row-entry
  end do
end if

! define the foil thickness, attenuation, and number slices per column
thick = foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
att = exp(-2.0*cPi*thick*xgp)  ! this is the global intensity attenuation factor; remove the factor 2.0 if amplitudes are needed
DF_nums = nint(thick/DF_slice)  ! this is the number of slices for this particular column


! next, deal with all the defects  This is the same as for the systematic row, except that in the zone
! axis we have two fundamental vectors and we can not use the same integer-based approach;  therefore,
! we have to store two floats for each slice in each column instead of a single integer.
!
! if there is a diplacement field file entered in the STEM_rundata.nml file,  
! then we simply read that file in; otherwise, we read all the defect descriptor files
if ((dispmode.eq.'new').or.(dispmode.eq.'not')) then

! is there a void data filename? If so, then read it  
   if (voidname.ne.'none') call read_void_data(numvoids,voidname,DF_L,DF_npix,DF_npiy,dinfo)

! read namelist files for all dislocations, if any
   if (numdisl.gt.0) call read_dislocation_data(dislname,numdisl,numsf,DF_npix,DF_npiy,DF_gf,DF_L,dinfo)

! read namelist files for all stacking faults, if any
   numYdisl = 0
   if (numsf.gt.0) call read_stacking_fault_data(numsf,numdisl,numYdisl,sfname,DF_L,DF_npix,DF_npiy,DF_g,dinfo)

! is there an inclusion data file? if so, then read it
   if (incname.ne.'none') call read_inclusion_data(numinc,incname,DF_L,DF_npix,DF_npiy,dinfo)

   
! the following will also have to be changed for the ZA case; it might not be necessary (or possible) 
! to pre-compute all the possible scattering matrices; perhaps we need to precompute an array of 
! 180x180 scattering matrices (about 40,000 of them) and then use bi-linear interpolation to select
! the right one for each slice.  Or, perhaps better in the long run, we simply compute the one we
! need when we need it...

! precompute ALL the defect columns and, if needed, store them in dispfile
! this portion should be carried out in multi-threaded mode as much as possible
  allocate(disparray(2,DF_nums,DF_npix,DF_npiy),imatvals(2,DF_nums))
  disparray = 0.0; imatvals = 0
  mess = 'displacement field computation '; call Message("(A)")
  
! initiate multi-threaded segment
!!$OMP     PARALLEL PRIVATE(TID,DF_R,imatvals,gdotR,i,j,k,imat) &
!!$OMP&   SHARED(NTHR,DF_npix,DF_npiy,DF_nums,numvoids,numdisl,numsf,numinc,disparray,t_interval)
!  NTHR = OMP_GET_NUM_THREADS()
!  TID = OMP_GET_THREAD_NUM()
!  write (*,*) TID,': entering parallel region'
 TID = 0
 NTHR = 1
  if (TID.eq.0) then
! do time reporting only in the master thread
    call Time_report(t_interval*float(NTHR)/float(DF_npix))
    call Time_start
  end if
  allocate(DF_R(DF_nums,3))     ! each thread has its own DF_R array 

 call TransSpace(float(ga),gac,'r','c')
 call TransSpace(float(gb),gbc,'r','c')
 
!!$OMP DO SCHEDULE (GUIDED)
!  write(*,*) TID,': starting Do Schedule'
  do i=1,DF_npix  
    do j=1,DF_npiy
      DF_R = 0.0
! compute the displacement vectors DF_R for all points in the column
      call CalcR(i,j,numvoids,numdisl,numYdisl,numsf,numinc)
! loop over the fixed thickness slices
      do k=1,DF_nums
! then convert to the dot-product 
       if (DF_R(k,1).eq.-10000.0) then  ! this is point inside a void
 	imatvals(1:2,k) = -10000
       else  ! it is not a void, so use the full dot product g.R (all vectors must be Cartesian !)
! use gac and gbc to get the two dot products and store both of them as integers mapped onto the 0..numd range
         gdotR = Dot_Product(gac,DF_R(k,1:3))
         imatvals(1,k) = numd*amod(gdotR+1000.0,1.0)
         gdotR = Dot_Product(gbc,DF_R(k,1:3))
         imatvals(2,k) = numd*amod(gdotR+1000.0,1.0)
       end if
     end do ! k loop
     disparray(1:2,1:DF_nums,i,j) = imatvals(1:2,1:DF_nums)
   end do
  if ((mod(i,t_interval).eq.0).and.(TID.eq.0)) call Time_remaining(i,DF_npix)
end do
!!$OMP END DO 
 if (TID.eq.0) call Time_stop(DF_npix*DF_npiy)
!!$OMP END PARALLEL
end if

! and, if needed, store the defect displacement field for re-runs
if (dispmode.ne.'not') then 
  if (dispmode.eq.'new') then 
    mess = 'Displacement field data will be stored in '//dispfile; call Message("(/A/)")
    open(unit=dataunit,file=dispfile,status='new',action='write',form='unformatted')
    i = 2
    write (dataunit) i,DF_nums,DF_npix,DF_npiy
    write (dataunit) disparray
    call SafeCloseFile('d1','keep',dispfile)
  else  ! there is a pre-computed defect file, so let's load it
   allocate(disparray(2,DF_nums,DF_npix,DF_npiy))
   disparray = 0
   open(unit=dataunit,file=dispfile,status='old',action='read',form='unformatted')
   read (dataunit) i,DF_nums_new,DF_npix_new,DF_npiy_new
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
  write (*,*) 'read displacement array from file ',dispfile
 end if
end if
io_real(1) = minval(disparray)
io_real(2) = maxval(disparray)
call WriteValue('disparray bounds: ', io_real, 2, "(2(F10.5,' '))")

! ok, all the set up is now complete;
! next, we prepare for the actual image simulation
! there are three types of simulation: CTEM, BF/HAADF, or STEM with storage of CBED patterns
  npix = DF_npix
  npiy = DF_npiy
  if (progmode.eq.'CTEM') then
    allocate(ZAimages(nn,npix,npiy,1))
  end if
  if (progmode.eq.'BFDF') then
      allocate(ZAimages(2,npix,npiy,STEM%numCL))
  end if
  if (progmode.eq.'STEM') then
! for STEM mode; we'll store data in blocks of 5 k-vectors to avoid issues with
! files that have too large a record size for IDL to read
      allocate(ZAimages(npix,npiy,nn,5))
  end if
  ZAimages = 0.0

! the following needs to be modified

! to enable threaded excution, we must first make a local copy of 
! the weights arrays STEM%ZABFweightsarray and STEM%ZAADFweightsarray,
! since OpenMP does not allow for this type of variable constructs
if (progmode.ne.'CTEM') then  
  allocate(BFweightsarray(nn,STEM%numk,STEM%numCL),ADFweightsarray(nn,STEM%numk,STEM%numCL))
  BFweightsarray = STEM%ZABFweightsarray
  ADFweightsarray = STEM%ZAADFweightsarray
! do the same with the excitation error array
  allocate(sgarray(nn,STEM%numk))
  sgarray = STEM%sgarray
end if

! loop over all reflections to get the appropriate powers
allocate(expval(2,nn,nn))
 rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
do ir=1,nn
   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
   do ic=1,nn
     if (ic.ne.ir) then  ! exclude the diagonal
       expval(1,ir,ic) = rltmpa%nab(1)-rltmpb%nab(1) 
       expval(2,ir,ic) = rltmpa%nab(2)-rltmpb%nab(2)
     end if
     rltmpb => rltmpb%next  ! move to next column-entry
  end do
  rltmpa => rltmpa%next   ! move to next row-entry
end do

! initialize the timer
if (progmode.eq.'CTEM') then
  numstart = 1
  numstop = 1
  call Time_report(float(t_interval)/float(npix))
else
  numstart = 1
  numstop = STEM%numk   
  io_int(1) = STEM%numk
  call WriteValue('STEM number of beam directions =  ', io_int, 1, "(I5)")
  call Time_report(1.0/float(numstop-numstart+1))
end if

nCL = STEM%numCL     ! set the number of camera length values

! define the numd complex defect parameters
  do i=0,numd
    arg = 2.D0*cPi*float(i)/dble(numd)
    para(i) = cmplx(cos(arg),sin(arg),dbl)
  end do

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

! various bits of useful information
  write (dataunit) xtalname					! crystal structure file names	
  write (dataunit) kk				               ! incident wave vector
  write (dataunit) CalcDiffAngle(ga(1),ga(2),ga(3))*0.5	! Bragg angle for the first reflection (whether allowed or not)
  write (dataunit) STEM%numberofsvalues			! number of pixels along disk radius
  write (dataunit) sngl(mLambda)				! wave length
  write (dataunit) beamdiv					! beam divergence
  write (dataunit) DF_L				        ! pixel size
	
! number of reflections, and associated information (hkl, ...)
  call TransSpace(float(kk),c,'d','c')
  call NormVec(c,'c')
! then make ga the x-axis
  call TransSpace(float(ga),gx,'r','c')
  call NormVec(gx,'c')
! compute the cross product between k and gx; this is the y-axis
  call CalcCross(c,gx,gy,'c','c',0)
  
  write (dataunit) nn   ! number of reflections
  rltmpa => reflist%next
  do ic = 1,nn
    gg=rltmpa%hkl
    call TransSpace(float(gg),c,'r','c')
    qx= CalcDot(c,gx,'c')
    qy= CalcDot(c,gy,'c')
    write (dataunit) rltmpa%hkl(1:3)
    write (dataunit) qx,qy
    rltmpa => rltmpa%next   ! move to next entry
  end do

! number of wave vectors, tangential components, etc...
  write (dataunit) numk   ! number of wave vectors
  ktmp => khead
  do ic=1,numk
    write (dataunit) ktmp%i, ktmp%j
    ktmp => ktmp%next
  end do

! for progmode=BFDF, we need to store the list of camera lengths
  if (progmode.eq.'BFDF') then
    write (dataunit) STEM%numCL
    write (dataunit) STEM%CLarray        ! note that this array has 20 entries !
  end if

if (progmode.eq.'STEM') then
! then prepare to store the individual blocks of data   
  write (dataunit) npix,npiy,nn,numk
  iSTEM = 0
end if

allocate(Sarray(nn,nn,0:numd,0:numd))

call Time_start
!--------------------------------------------------------------
!--------------------------------------------------------------
!--------------------------------------------------------------
mainloop: do isg = numstart,numstop   ! this is the main computational loop
 iSTEM = iSTEM+1
!--------------------------------------------------------------
! here we precompute an array of scattering matrices that can 
! then be used, either directly, or via bi-linear interpolation,
! by the image computation portion of this program.
!
! For starters, we'll subdivide the range of possible alpha values
! in 180 segments (2 degrees each), with a copy of the last one
! (i.e., 181x181 = 32761 entries or 240 Mb for 31 beams)
!
! this part is essentially the same as the threaded section of the older
! ZAdefect.all.f90 program (which was a test program).
!

! get the correct excitation errors for this beam orientation (in STEM mode);
! no need to do anything in CTEM mode
if (progmode.ne.'CTEM') then  
! fill the diagonal of the reference dynamical matrix and the void matrix
  forall (i=1:nn)
   DHWMz(i,i)=2.0*cPi*cmplx(0.0,sgarray(i,isg))    ! initialize the diagonal elements of the dynamical matrix
   DHWMvoid(i,i) = DHWMz(i,i)
  end forall
end if


NTHR = SETNTHR

!$OMP  PARALLEL NUM_THREADS(NTHR) PRIVATE(TID,i,j,k,ii,jj,ic,ir,g,Azz,DDD,zmax) &
!$OMP& SHARED(NTHR,Sarray,thr,para,t_interval,nn,DHWMz,DF_slice,expval,cone,czero,numd)    
TID = OMP_GET_THREAD_NUM() 

allocate(Azz(nn,nn), DDD(nn,nn))   ! these are private variables, so each thread must allocate them !

!$OMP DO SCHEDULE(STATIC)
do j=0,numd
 do i=0,numd 
! loop over all reflections in the array DD using the information in expval
! ir is the row index
  do ir=1,nn
! ic is the column index
   do ic=1,nn
    if (ic.ne.ir) then  ! exclude the diagonal
     DDD(ir,ic) = DHWMz(ir,ic) * para(i)**expval(1,ir,ic) * para(j)**expval(2,ir,ic)
    else
     DDD(ir,ic) = DHWMz(ir,ic) 
    end if
   end do
  end do
    
  call MatrixExponential(DDD, Azz, dble(DF_slice), 'Pade', nn)  
    
   Sarray(1:nn,1:nn,i,j) = Azz(1:nn,1:nn)
 end do
 if ((TID.eq.0).and.(mod(j,10).eq.0)) then
    mess = '.'; call Message("(A1,$)")
 end if
end do
!$OMP END DO

deallocate(Azz, DDD)
!$OMP END PARALLEL

mess = ' 181x181 scattering matrices precomputed '; call Message("(A,' ',$)")

!----------------------------------------------------!
! Finally, here it is: the actual image computation  !
!----------------------------------------------------!

NTHR = SETNTHR
!$OMP  PARALLEL NUM_THREADS(NTHR) PRIVATE(TID,i,j,k,ii,Azz,amp,amp2,ix,iy,dx,dy,dxm,dym,inten,ixp,iyp) &
!$OMP& SHARED(NTHR,Sarray,DF_Svoid,disparray,t_interval,nn,cone,czero,images,npix,npiy,att,DF_nums)    
TID = OMP_GET_THREAD_NUM() 
  
 allocate(Azz(nn,nn),amp(nn),amp2(nn),inten(nn))
 
!$OMP DO SCHEDULE (STATIC)
 donpix: do i=1,npix
 if ((TID.eq.0).and.(mod(i,10).eq.0)) then
   mess = '.'
   call Message("(A1,$)")
 end if
 donpiy:   do j=1,npiy
! initialize the wave function for this pixel with (1.0,0.0) for the incident beam
    amp = czero
    amp(1) = cone
    doslices: do k=1,DF_nums    ! loop over the fixed thickness slices
! compute the appropriate scattering matrix to propagate with (see section 8.3.3 in the book)
       if (disparray(1,k,i,j).eq.-10000) then  ! this is point inside a void
 	 Azz = DF_Svoid    ! so we use the void propagator matrix
       else  ! it is not a void
! in this version, we use complex bi-variate interpolation to select the appropriate Azz 
! matrix from the pre-computed Sarray.
         ix = int(disparray(1,k,i,j))
         ixp = ix+1
         if (ix.eq.numd) ixp=1
         iy = int(disparray(2,k,i,j))
         iyp = iy+1
         if (iy.eq.numd) iyp=1
         dx = cmplx(amod(disparray(1,k,i,j),1.0),0.0)
         dy = cmplx(amod(disparray(2,k,i,j),1.0),0.0)
         dxm = cone-dx
         dym = cone-dy
	 Azz = dxm*dym*Sarray(1:nn,1:nn,ix,iy)+dx*dym*Sarray(1:nn,1:nn,ixp,iy)+ &
                    dxm*dy*Sarray(1:nn,1:nn,ix,iyp)+dx*dy*Sarray(1:nn,1:nn,ixp,iyp)
       end if
! and multiply with this matrix
       amp2 = matmul(Azz, amp)
!       do ii=1,nn
!           amp2(ii) = sum(Azz(ii,1:nn) * amp(1:nn))
!       end do
       amp = amp2
    end do doslices ! loop over slices 
      
! compute the intensities and (for STEM mode) use the proper weight-factor
    inten(1:nn) = att*cdabs(amp(1:nn))**2
    if (progmode.eq.'STEM') then
        ZAimages(i,j,1:nn,iSTEM) = inten(1:nn)
    end if
    if (progmode.eq.'BFDF') then
        do ii=1,nn
         do jj=1,nCL
          if (BFweightsarray(ii,isg,jj).eq.1) ZAimages(1,i,j,jj) = ZAimages(1,i,j,jj) + inten(ii)
          if (ADFweightsarray(ii,isg,jj).eq.1) ZAimages(2,i,j,jj) = ZAimages(2,i,j,jj) + inten(ii)
         end do
        end do
      end if
    if (progmode.eq.'CTEM') then
      ZAimages(1:nn,i,j,1) = inten(1:nn)
    end if

    end do donpiy
end do donpix
!$OMP END DO
  deallocate(Azz,amp,amp2,inten)
!$OMP END PARALLEL

  
  if (progmode.ne.'CTEM') call Time_remaining(isg-numstart+1,numstop-numstart+1)

  if (progmode.eq.'STEM') then
    if (mod(isg,5).eq.0) then
!      mess = 'Storing block of 5 k-vector results'; call Message("(A)")
      write (dataunit) ZAimages
      iSTEM = 0
      ZAimages = 0.0
    end if
  end if

200 end do mainloop

deallocate(Sarray)

! ok, so the main computation is complete; print some timing information
call Time_stop(npix*npiy)

! and save the data in the appropriate format
if (progmode.eq.'CTEM') then
  write (dataunit) nn,DF_npix,DF_npiy,1
  write (dataunit) ZAimages
end if 

if (progmode.eq.'BFDF') then
  write (dataunit) 2,DF_npix,DF_npiy,STEM%numCL
  write (dataunit) ZAimages
end if 

! flush the final ZAimages array if necessary (depends on value of iSTEM)
if ((progmode.eq.'STEM').AND.(iSTEM.ne.0)) then
  write (dataunit) ZAimages
end if 

! in all cases we need to close the file here...
close(UNIT=dataunit,STATUS='keep')

mess = 'Data stored in file '//trim(dataname)
call Message("(/A/)")

end subroutine ComputeZAdefect

