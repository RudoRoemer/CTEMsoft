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
! CTEMsoft2013:CTEMECCI.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMECCI
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief tryECCI computes ECCI defect contrast for multiple defects
!
!> @details This is a test version based on the CTEMZAdefect program
! 
!> @date  03/07/10 MDG  1.0 original, based on STEMdefect.f90 (collaboration with OSU) 
!> @date  03/21/10 MDG  2.0 modified for STEM illumination, fast version with bi-variate interpolation
!> @date  06/19/13 MDG  3.0 conversion to new libraries 
!> @date  10/28/13 MDG  3.1 added Interpret_Program_Arguments line
!> @date  12/03/13 MDG  4.0 new start 
!> @date  12/08/13 MDG  4.1 added trace (line scan) mode
!> @date  02/10/14 MDG  4.2 added apbs
!--------------------------------------------------------------------------
program CTEMECCI

use local
use files
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMECCI.nml'
progname = 'CTEMECCI.f90'
call Interpret_Program_Arguments(nmldeffile,8,(/ 0, 3, 41, 200, 201, 202, 203, 204 /) )

! initialize all user-defined variables
call ComputeECCI(nmldeffile)

end program CTEMECCI


!--------------------------------------------------------------------------
!
! SUBROUTINE: ComputeECCI
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a zone axis ECCI defect data set
!
!> @param nmlfile namelist file name
!
!> @todo check sign of kt for consistency 
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 11/13/13  MDG 2.2 implementation of Pade approximation for scattering matrix computation
!> @date 12/04/13  MDG 3.0 new implementation of the thickness integration, nmore along the lines of CTEMECP 
!> @date 12/07/13  MDG 3.1 added line scan mode (called "trace")
!> @date 02/10/14  MDG 3.2 added apbs
!> @date 02/24/14  MDG 3.3 removal of double-counted phase factor
!> @date 03/05/14  MDG 3.4 correction of integration to double integration for Lgh array
!> @date 03/12/14  MDG 3.5 conversion of Sgh and Lgh arrays to diagonal only
!--------------------------------------------------------------------------
subroutine ComputeECCI(nmlfile)

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
use YSHmodule
use void
use inclusion
use apb
use defectmodule
use rotations
use timing
use STEMmodule

IMPLICIT NONE

character(fnlen),INTENT(IN)		:: nmlfile

integer(kind=irg)    		        :: nn,i,j,k,npix,npiy,ii,jj,numvoids,numdisl, numset, &
					numYdisl,numsf,numinc,numapb,dinfo,t_interval,nat(100), &
					DF_nums_new,DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, &
					NTHR, isym, ir, ga(3), gb(3),kk(3),ic,g,numd,ix,iy,nkt,nbeams, ik, ig, &
					numk,ixp,iyp,SETNTHR, io_int(6), skip, gg(3), iSTEM, nktstep
!                                  OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
integer(kind=irg),parameter 		:: numdd=360 ! 180
real(kind=sgl)         		:: thick, X(2), dmin, dkt, bragg, thetac, kstar(3), gperp(3), &
					lauec(2),lauec2(2),gdotR,DF_gf(3), tpi, &
					DM(2,2), DD, c(3), gx(3), gy(3), &
					gac(3), gbc(3),zmax, ktmax, io_real(2), voltage, ijmax
real(kind=dbl)                        :: arg, glen
character(fnlen)      			:: dataname,sgname,voidname,dislname(3*maxdefects),sfname(maxdefects),ECPname, &
					incname,dispfile,xtalname,foilnmlfile, STEMnmlfile,dislYname(3*maxdefects), apbname
character(4)            		:: dispmode, summode
character(5)                          :: progmode
complex(kind=dbl),allocatable    	:: DHWM(:,:),DHWMvoid(:,:),DDD(:,:),Sarray(:,:,:,:)
complex(kind=dbl),allocatable    	:: amp(:),amp2(:),Azz(:,:)
complex(kind=dbl)                	:: czero,cone
complex(kind=dbl)                	:: para(0:numdd),dx,dy,dxm,dym, xgp
real(kind=sgl),allocatable       	:: sgarray(:,:)
real(kind=sgl),allocatable    		:: disparray(:,:,:,:),imatvals(:,:), ECCIimages(:,:)
integer(kind=sgl),allocatable    	:: expval(:,:,:)
complex(kind=dbl),allocatable         :: Lgh(:),Sgh(:) ! Lgh(:,:),Sgh(:,:)
logical	                               :: ECCI, NANCHK

namelist / ECCIlist / DF_L, DF_npix, DF_npiy, DF_slice, dmin, sgname, numvoids, incname, stdout, &
                                voidname, numdisl, dislname, numYdisl, dislYname, numsf, sfname, dinfo, &
				 t_interval,progmode, dispfile, ktmax, dkt, ECPname, summode, lauec, lauec2, &
				 dispmode,SETNTHR,xtalname,voltage,kk, lauec, nktstep, &
				 dataname, foilnmlfile, STEMnmlfile, apbname

ECCI = .TRUE.

! first we define the default values
czero=dcmplx(0.D0,0.D0)
cone=dcmplx(1.D0,0.D0)
tpi = 2.0*sngl(cPi)

! parameters specific to this run
 xtalname = 'undefined'		! initial value; MUST be present in nml file for program to execute
 voltage = 20000.0			! accelerating voltage
 kk = (/ 0, 0, 1 /)			! incident wave vector in crystal components (omitting wave length)
 lauec = (/ 0.0,0.0 /)			! Laue center coordinates (used for single image mode, start point in trace modea)
 lauec2 = (/ 0.0,0.0 /)		! Laue center 2 coordinates (used for trace mode, end point)
 nktstep = 10                         ! number of steps in line scan mode
 dmin = 0.04			        ! smallest d-spacing to include in dynamical matrix [nm]
 ktmax = 5.0
 dkt = 0.5
 progmode = 'array'                   ! 'array' for array of images, 'trace' for line scan
 summode = 'diag'                     ! 'full' for complete summation, 'diag' for diagonal only

! CTEM or STEM ?
 progmode = 'CTEM'  		        ! default illumination mode (can be 'CTEM' or 'STEM')
 STEMnmlfile = 'STEM_rundata.nml'	! name of the STEM rundata namelist file
 foilnmlfile = 'FOIL_rundata.nml'	! name of the foil rundata namelist file
 ECPname = 'undefined'                ! name of the corresponding ECP file (must exist!)
 
! column approximation parameters and image parameters 
 DF_L = 1.0             		! edge length of column in nanometers
 DF_npix = 256       			! number of image pixels along x
 DF_npiy = 256       			! number of image pixels along y 
 DF_slice = 1.0       			! slice thickness in nanometers

 dinfo = 0               		! switch to make makedislocation verbose
 sgname = 'nofile'   			! if this variable is different from 'nofile', then an external sg array is read (to be implemented)

! defect parameters
 numdisl = 0           	! number of dislocation files
 numYdisl = 0           	! number of Yoffe dislocation files
 numsf = 0             	! number of stacking fault files
 numinc = 0           		! number of inclusions
 numvoids = 0       		! number of voids
 voidname = 'none' 		! filename for void data
 dislname = ''         	! filenames for dislocation data
 dislYname = ''         	! filenames for Yoffe dislocation data
 sfname = ''            	! filenames for stacking fault data
 incname = 'none'   		! filename for inclusion data
 apbname = 'none'              ! filename for (circular) apbs
 dispfile = 'none'     	! name of the displacement field output file (will be created if different from none)
 dispmode = 'not'  		! should a diplacement file be written ('new') or read ('old') or neither ('not')?

! output parameters
 dataname = 'ECCIdefect.data'	! default outputfile name
 t_interval = 10       	! default timing interval (output every t_interval image columns)
 
! then we read the actual rundata namelist, which may override some or all of these defaults  
 OPEN(UNIT=dataunit,FILE=trim(nmlfile),DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=ECCIlist)
 CLOSE(UNIT=dataunit)
  
! make sure the xtalname variable has been properly defined
if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMECCI:',' structure file name is undefined in '//nmlfile)
end if

! make sure the ECPname variable has been properly defined
if (trim(ECPname).eq.'undefined') then
  call FatalError('CTEMECCI:',' ECP pattern file name is undefined in '//nmlfile)
end if

! we got this far, so display the standard program info
 progname = 'CTEMECCI.f90'
 progdesc = 'Dynamical zone axis ECCI defect image simulation'
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
 call TransSpace(DF_gf,DF_gc,'r','c')         ! convert to Cartesian reference frame

! we'll need to compute the list of wavevectors at this point
  nkt = nint( ktmax / dkt)
  ijmax = nkt**2
  bragg = CalcDiffAngle(ga(1),ga(2),ga(3)) * 0.5
  if (ktmax.eq.0.0) then
    thetac = bragg
  else
    thetac = (ktmax * 2.0 * bragg) * 0.5
  end if
! here we figure out how many beams there are
  if (progmode.eq.'array') then 
!    call Calckvectors(dble(kk),dble(ga),dble(ktmax),nkt,nkt,numk,isym,ijmax,'Conical')   
    call Calckvectorcone(kk,ga,lauec(1),lauec(2),ktmax,nkt,numk)
    io_int(1)=numk 
    call WriteValue('Total number of independent incident beam directions inside cone = ', io_int, 1,"(I8)")
  end if
  if (progmode.eq.'trace') then  
! single image mode, potentially with a narrow range of incident beam directions
    call Calckvectortrace(kk,ga,lauec(1),lauec(2),lauec2(1),lauec2(2),1.0,nktstep,numk)
    io_int(1)=numk 
    call WriteValue('Total number of independent incident beam directions along trace = ', io_int, 1,"(I8)")
  end if

! and determine the overall reflection list
  call Compute_ReflectionList(dmin,kk,ga,gb,'ALL',.FALSE.,0,thetac)

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

! print the list for debugging purposes...
 rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
  do ir=1,nn
    write (*,*) ir,': ',rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3)
    rltmpa => rltmpa%next
  end do


! ideally, we should use Bethe potentials to reduce the size of the dynamical matrix;
! while the theory has been worked out to do this, it would require tremendous changes
! to the program starting about here; given the time limitations, there won't be any 
! chance to do this before the end of the year (2013), so we'll leave that for some other time...
  
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

 cone = dcmplx(0.D0,cPi)
 write (*,*) 'complex iPi = ',cone

 rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
  do ir=1,nn
   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
   do ic=1,nn
    if (ic.ne.ir) then  ! exclude the diagonal
! compute Fourier coefficient of electrostatic lattice potential 
     call CalcUcg(rltmpa%hkl - rltmpb%hkl)
     DHWMz(ir,ic) = rlp%qg * cone
!cmplx(- cPi * aimag(rlp%qg), cPi * real(rlp%qg),dbl)  ! and initialize the off-diagonal matrix element (including i Pi)
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

! compute the normal absorption factor xgp (which equals rlp%qg with g=0)
  rlp%qg = cmplx(0.D0,0.D0)
  call CalcUcg((/0,0,0/))
  xgp = cmplx(-cPi/rlp%xgp,0.0) ! cPi * rlp%qg * cone
  io_real(1) = rlp%xgp
  call WriteValue('Normal absorption length : ', io_real, 1, "(F10.5/)")
  write (*,*) 'rlp%qg = ',rlp%qg
write (*,*) 'i Pi / q_0 = ',xgp

! define the foil thickness, attenuation, and number slices per column
  thick = foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
  DF_nums = nint(thick/DF_slice)  ! this is the number of slices for this particular column
write (*,*) 'foil thickness', foil%zb, thick, DF_nums

! next, deal with all the defects  This is the same as for the systematic row, except that in the zone
! axis case we have two fundamental vectors and we can not use the same integer-based approach;  therefore,
! we have to store two floats for each slice in each column instead of a single integer.
!


! if there is a diplacement field file entered in the STEM_rundata.nml file,  
! then we simply read that file in; otherwise, we read all the defect descriptor files
if ((dispmode.eq.'new').or.(dispmode.eq.'not')) then

! is there a void data filename? If so, then read it  
   if (voidname.ne.'none') call read_void_data(numvoids,voidname,DF_L,DF_npix,DF_npiy,dinfo)

! read namelist files for all dislocations, if any
   if (numdisl.gt.0) call read_dislocation_data(dislname,numdisl,numsf,DF_npix,DF_npiy,DF_gf,DF_L,dinfo)

! read namelist files for all Yoffe dislocations, if any
   if (numYdisl.gt.0) call read_YSH_dislocation_data(dislYname,numYdisl,DF_npix,DF_npiy,DF_gf,DF_L,dinfo)

! read namelist files for all stacking faults, if any
   if (numsf.gt.0) call read_stacking_fault_data(numsf,numdisl,numYdisl,sfname,DF_L,DF_npix,DF_npiy,DF_g,dinfo,ECCI)

! is there an inclusion data file? if so, then read it
   if (incname.ne.'none') call read_inclusion_data(numinc,incname,DF_L,DF_npix,DF_npiy,dinfo)

! is there an apb file ?
   if (apbname.ne.'none') call read_apb_data(numapb,apbname,DF_L,DF_npix,DF_npiy,dinfo)
   
! the following will also have to be changed for the ZA case; it might not be necessary (or possible) 
! to pre-compute all the possible scattering matrices; perhaps we need to precompute an array of 
! 180x180 scattering matrices (about 40,000 of them) and then use bi-linear interpolation to select
! the right one for each slice.  Or, perhaps better in the long run, we simply compute the one we
! need when we need it...

! precompute ALL the defect columns and, if needed, store them in dispfile
! this portion should be carried out in multi-threaded mode as much as possible
  allocate(disparray(2,DF_nums,DF_npix,DF_npiy),imatvals(2,DF_nums))
  disparray = 0.0; imatvals = 0
  
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
      call CalcR(i,j,numvoids,numdisl,numYdisl,numsf,numinc,numapb)
! loop over the fixed thickness slices
      do k=1,DF_nums
! then convert to the dot-product 
       if (DF_R(k,1).eq.-10000.0) then  ! this is point inside a void
 	imatvals(1:2,k) = -10000
       else  ! it is not a void, so use the full dot product g.R (all vectors must be Cartesian !)
! use gac and gbc to get the two dot products and store both of them as integers mapped onto the 0..numd range
         gdotR = Dot_Product(gac,DF_R(k,1:3))
! due to the logarithmic singularity at a dislocation core, it is possible for gdotR to be NaN; we need
! to intercept these cases, and replace the value of gdotR by 0.0
         if (NANCHK(gdotR)) gdotR = 0.0         
         imatvals(1,k) = numd*amod(gdotR+1000.0,1.0)
         gdotR = Dot_Product(gbc,DF_R(k,1:3))
         if (NANCHK(gdotR)) gdotR = 0.0         
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
  
  allocate(ECCIimages(npix,npiy))
  ECCIimages = 0.0

! compute the excitation error array
  allocate(sgarray(nn,numk))
! loop over the wave vector linked list
  ktmp => khead
  beamloopCL: do ik=1,numk
!    ll = ktmp%kt        ! this is the tangential component of the wave vector
! and loop over all reflections
    rltmpa => reflist%next
    reflectionloopCL: do ig=1,nn
      gg = float(rltmpa%hkl)
!      glen = CalcLength(dble(gg),'r')
!      lpg = ll + gg                ! Laue + g
!      gplen = CalcLength(lpg,'r')
!      kpg = 2000.0*asin(0.50*sngl(mLambda)*gplen)    ! 2theta in mrad
      sgarray(ig,ik) = Calcsg(float(gg),sngl(ktmp%k),DynFN)
 ! and we move to the next reflection in the list
      rltmpa => rltmpa%next
    end do reflectionloopCL  
    ktmp => ktmp%next
  end do beamloopCL

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

! define the numd complex defect parameters
  do i=0,numd
    arg = 2.D0*cPi*dble(i)/dble(numd)
    para(i) = dcmplx(dcos(arg),-dsin(arg))
  end do

! determine the Sgh array, which is sort of a glorified structure factor...
! since all incident beam orientations use the same number of reflections,
! we can compute Sgh here instead of inside the main loop.  This will need
! to be modified when we start using Bethe potentials.

  mess = 'Computing Sgh array'; call Message("(A)")

  numset = cell % ATOM_ntype  ! number of special positions in the unit cell

  allocate(Sgh(nn))
  call CalcSgh(nn,Sgh,nat)
  


! store necessary data in data file
  mess ='Storing data for IDL visualization program in '//dataname
  call Message("(A/)")
  open(UNIT=dataunit,FILE=trim(dataname),STATUS='unknown',FORM='unformatted')  
  
! program mode
  write (dataunit) 'ECCI'

! filename of corresponding data set
  write (dataunit) dataname

! filename of corresponding ECP pattern information
  write (dataunit) ECPname

! program mode
  write (dataunit) progmode

! summation mode 
  write (dataunit) summode

! various bits of useful information
  write (dataunit) xtalname					! crystal structure file names	
  write (dataunit) kk						! incident wave vector
  write (dataunit) CalcDiffAngle(ga(1),ga(2),ga(3))*0.5      ! Bragg angle for the first reflection (whether allowed or not)
  if (progmode.eq.'array') then 
    write (dataunit) ktmax		                       ! max beam convergence in units of |g_a|
  else
    write (dataunit) 1.0		                       ! max beam convergence in units of |g_a|
  end if

  write (dataunit) dkt				               ! angular step size along disk radius
  write (dataunit) voltage				        ! microscope voltage
  write (dataunit) thetac					! beam divergence
  write (dataunit) DF_L					! pixel size

! number of reflections, and associated information (hkl, ...)
  call TransSpace(float(kk),c,'d','c')
  call NormVec(c,'c')
! then make ga the x-axis
  call TransSpace(float(ga),gx,'r','c')
  call NormVec(gx,'c')
! compute the cross product between k and gx; this is the y-axis
  call CalcCross(c,gx,gy,'c','c',0)
  
  write (dataunit) nn   ! number of reflections


! this needs to be fixed !!!!!

  call TransSpace(float(kk),kstar,'d','r')        ! transform incident direction to reciprocal space
  call CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
  call NormVec(gperp,'r')                        ! normalize g_perp

 DM(1,1) = CalcDot(gperp,gperp,'c')
 DM(1,2) = -CalcDot(float(ga),gperp,'c')
 DM(2,1) = DM(1,2)
 DM(2,2) = CalcDot(float(ga),float(ga),'c')
 DM = transpose(DM)
 DD = DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1)
 glen = CalcLength(float(ga),'r')



! number of wave vectors, tangential components, etc...
  write (dataunit) numk   ! number of wave vectors
  ktmp => khead
!  if (progmode.eq.'trace') then  
    do ic=1,numk
      X(1) = -CalcDot(sngl(ktmp%kt),float(ga),'c') ! / glen
      X(2) = -CalcDot(sngl(ktmp%kt),gperp,'c') * glen
!      write (*,*) X(1), X(2), ktmp%kt
      write (dataunit) X(1), X(2)
      ktmp => ktmp%next
    end do
!  end if

  
  write (dataunit) DF_npix,DF_npiy
  
  allocate(Sarray(nn,nn,0:numd,0:numd))

! initialize the timer
  numstart = 1
  numstop = numk   
  io_int(1) = numk
  call WriteValue('ECCI: number of beam directions =  ', io_int, 1, "(I5)")
  call Time_report(1.0/float(numstop-numstart+1))
  call Time_start
  
  cone = dcmplx(1.D0,0.D0)
  write (*,*) 'cone = ',cone
  write (*,*) 'czero = ',czero
  
  
!--------------------------------------------------------------
!--------------------------------------------------------------
!--------------------------------------------------------------
mainloop: do isg = numstart,numstop   ! this is the main computational loop
 iSTEM = iSTEM+1
 ECCIimages = 0.0

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
! fill the diagonal of the reference dynamical matrix and the void matrix
! make sure that normal absorption is properly taken into account...
  forall (i=1:nn)
   DHWMz(i,i)= dcmplx(0.D0,2.D0*cPi*sgarray(i,isg)) + xgp ! xgp already has i Pi in it.
   DHWMvoid(i,i) = DHWMz(i,i)
  end forall

NTHR = SETNTHR

!$OMP  PARALLEL NUM_THREADS(NTHR) DEFAULT(SHARED) PRIVATE(TID,i,j,k,ii,jj,ic,ir,g,Azz,DDD,zmax)
TID = OMP_GET_THREAD_NUM() 
!TID = 0
allocate(Azz(nn,nn), DDD(nn,nn))   ! these are private variables, so each thread must allocate them !

!$OMP DO SCHEDULE(STATIC)
do j=0,numd
 do i=0,numd 
! loop over all reflections in the array DD using the information in expval
! ic is the column index
 do ic=1,nn
! ir is the row index
    do ir=1,nn
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

mess = ' Scattering matrices precomputed '; call Message("(A,' ',$)")

!----------------------------------------------------!
! Finally, here it is: the actual image computation  !
! This was modified from the regular (S)TEM mode to  !
! reflect the depth integration, which requires a    !
! summation of the product of Sgh and Lgh.           !
!----------------------------------------------------!

 NTHR = SETNTHR
!$OMP  PARALLEL NUM_THREADS(NTHR) DEFAULT(SHARED) PRIVATE(TID,i,j,k,ii,Azz,amp,amp2,ix,iy,dx,dy,dxm,dym,ixp,iyp,Lgh,ir,ic)
 TID = OMP_GET_THREAD_NUM() 
! TID = 0 
 allocate(Azz(nn,nn),amp(nn),amp2(nn),Lgh(nn))
 
!$OMP DO SCHEDULE (STATIC)
 donpix: do i=1,npix
 if ((TID.eq.0).and.(mod(i,10).eq.0)) then
   mess = '.'
   call Message("(A1,$)")
 end if
 donpiy:   do j=1,npiy
! initialize the wave function for this pixel with (1.0,0.0) for the incident beam
    Lgh = czero
    amp = czero
    amp(1) = cone
!    if ((TID.eq.0).and.(isg.eq.1).and.(i.eq.1).and.(j.eq.1)) then
!      write (*,*) 'Storing data for comparison with ECP program; dimensions : ',nn,DF_nums
!      open(unit=dataunit,file='ECCIcheck.data',status='unknown',action='write',form='unformatted')
!      write (dataunit) nn,DF_nums
!    end if

    doslices: do k=1,DF_nums    ! loop over the fixed thickness slices
! compute the appropriate scattering matrix to propagate with (see section 8.3.3 in the book)
       if (disparray(1,k,i,j).eq.-10000) then  ! this is point inside a void
 	 Azz = DHWMvoid    ! so we use the void propagator matrix
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

       amp2 = matmul(Azz,amp)

       if (k.eq.1) then 
	  Lgh = abs(amp2)**2
       else
	  Lgh = Lgh + abs(amp2)**2
       end if


       amp = amp2
       
!    if ((isg.eq.1).and.(i.eq.1).and.(j.eq.1)) then
!      write (*,*) k,maxval(cabs(amp)**2)
!!         write (dataunit) Lgh
!    end if
       
    end do doslices ! loop over slices 

! then we need to multiply Sgh and Lgh, sum, and take the real part which will
! produce the desired BSE intensity
    ECCIimages(i,j) = sngl(real(sum( Sgh * Lgh )))
 
    end do donpiy
end do donpix
!$OMP END DO
  deallocate(Azz,amp,amp2,Lgh)
!$OMP END PARALLEL
  
  ECCIimages = ECCIimages / float(DF_nums) / float(sum(nat))
  call Time_remaining(isg-numstart+1,numstop-numstart+1)

  write (dataunit) ECCIimages

200 end do mainloop

close(UNIT=dataunit,STATUS='keep')

! ok, so the main computation is complete; print some timing information
call Time_stop(npix*npiy)

mess = 'Data stored in file '//trim(dataname)
call Message("(/A/)")

end subroutine ComputeECCI


!--------------------------------------------------------------------------
!
! SUBROUTINE: Calckvectorcone
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a set of incident beam directions for single image ECCI mode
!
!> @param k incident wave vector (zone axis)
!> @param ga principal g vector
!> @param ktx tangential x component
!> @param kty tangential y component
!> @param ktrad cone opening angle
!> @param ktstep number of steps along cone radius
!> @param numk resulting number of incident beam directions
!
!> @date 11/29/01  MDG 1.0 original
!> @date 12/05/13  MDG 2.0 adapted for ECCI simulations 
!--------------------------------------------------------------------------
subroutine Calckvectorcone(k,ga,ktx,kty,ktrad,ktstep,numk)

use io
use error
use diffraction
use crystal
use kvectors
use dynamical

IMPLICIT NONE

integer(kind=irg),INTENT(IN)        :: k(3)
integer(kind=irg),INTENT(IN)        :: ga(3)
real(kind=sgl),INTENT(IN)           :: ktx
real(kind=sgl),INTENT(IN)           :: kty
real(kind=sgl),INTENT(IN)           :: ktrad
integer(kind=irg),INTENT(IN)        :: ktstep
integer(kind=irg),INTENT(OUT)       :: numk

integer                             :: istat,imin,imax,jmin,jmax,ijmax,i,j,ic,jc,ki
real                                :: kr(3),glen,delta,kstar(3),kt(3),gan(3),gperp(3),ktlen, dkt

! compute geometrical factors 
 glen = CalcLength(float(ga),'r')              ! length of ga
 gan = ga/glen                                 ! normalized ga
! delta = 2.0*ktrad*glen/float(2*ktstep+1)     ! grid step size in nm-1 
 delta = ktrad*glen/float(ktstep)              ! grid step size in nm-1 
 dkt = ktrad/float(ktstep)
!write (*,*) ktrad, ktstep, glen, delta
 call TransSpace(float(k),kstar,'d','r')       ! transform incident direction to reciprocal space
 call CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
 call NormVec(gperp,'r')                       ! normalize g_perp
 call NormVec(kstar,'r')                       ! normalize reciprocal beam vector

!  kt = -(klaue(1)+float(i)*delta)*gan - (klaue(2)+float(j)*delta)*gperp  ! tangential component of k

! deal only with the incident beam (parallel illumination)
if (ktstep.eq.0) then
 if (.not.associated(khead)) then     ! allocate the head and ktail of the linked list
   allocate(khead,stat=istat)         ! allocate new value
   if (istat.ne.0) call FatalError('Calckvectorcone: unable to allocate head pointer',' ')
   ktail => khead                      ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   numk = 1                          ! keep track of number of k-vectors so far
 ! this should be the center vector of the illumination cone !!!
   kt = - glen * (ktx*gan + kty * gperp)
!   kt = - glen * (-ktx*gan + kty * gperp)
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktx**2+kty**2)         ! squared length of tangential component
   
!   write (*,*) 0, 0, kt, ktlen
   
   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end if
else
! next, put the center of the cone in units of (i,j) (original ECP "screen" coordinates)
  ic = int(ktx*glen/delta)
  jc = int(kty*glen/delta)
  ki = ktstep
!write (*,*) ktx, kty, ktx*glen/delta, kty*glen/delta 

 if (.not.associated(khead)) then     ! allocate the head and ktail of the linked list
   allocate(khead,stat=istat)         ! allocate new value
   if (istat.ne.0) call FatalError('Calckvectorcone: unable to allocate head pointer',' ')
   ktail => khead                      ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   numk = 1                          ! keep track of number of k-vectors so far
 ! this should be the center vector of the illumination cone !!!
   ktail%i = ic                            ! i-index of beam
   ktail%j = jc                            ! j-index of beam
!   kt = float(ktail%i)*delta*gan - float(ktail%j)*delta*gperp  ! tangential component of k
   kt = -float(ktail%i)*delta*gan - float(ktail%j)*delta*gperp  ! tangential component of k

!   kt = delta * (ktx * gan - kty * gperp ) ! tangential component of k
! write (*,*) ic,jc,kt
   ktail%kt = kt                           ! store tangential component of k
   ktlen = delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component

!   write (*,*) ic, jc, kt, ktlen

   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 else
   call FatalError('Calckvectorcone: pointer head already allocated',' ')
 end if

! the following lines are quite different if symmetry is taken into account;
! check the MBsym.f90 program to determine how that can be done.
  imin =  -ki; imax = ki; jmin = -ki; jmax = ki; 
  ijmax = ki**2
! now do the real work
  do i=imin,imax
   do j=jmin,jmax
    if (.not.((i.eq.0).and.(j.eq.0))) then  ! the point (0,0) has already been taken care of
     if ((i**2+j**2).le.ijmax) then   ! is point inside the incident cone ?
      allocate(ktail%next,stat=istat)  ! allocate new value
      if (istat.ne.0) call FatalError('Calckvectorcone: unable to allocate pointer',' ')
      ktail => ktail%next               ! ktail points to new value
      nullify(ktail%next)              ! nullify next in new value
      numk = numk + 1                 ! keep track of number of k-vectors so far
      ktail%i = ic+i                   ! i-index of beam
      ktail%j = jc+j                   ! j-index of beam
!      kt = float(ktail%i)*delta*gan - float(ktail%j)*delta*gperp  ! tangential component of k
      kt = - float(ktail%i)*delta*gan - float(ktail%j)*delta*gperp  ! tangential component of k
!     kt = delta * ((ktx + float(i)*dkt) * gan - (kty + float(j)*dkt) * gperp ) ! tangential component of k
! write (*,*) ic+i,jc+j,kt
      ktail%kt = kt                    ! store tangential component of k
      ktlen = delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component
 
!   write (*,*) i, j, kt, ktlen

     kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
     end if
    end if
   end do
  end do
end if

end subroutine Calckvectorcone




!--------------------------------------------------------------------------
!
! SUBROUTINE: Calckvectortrace
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a set of incident beam directions for line scan ECCI mode
!
!> @param k incident wave vector (zone axis)
!> @param ga principal g vector
!> @param ktx tangential x component of trace start point
!> @param kty tangential y component
!> @param kt2x tangential x component of trace end point
!> @param kt2y tangential y component
!> @param ktrad cone opening angle
!> @param ktstep number of steps along cone radius
!> @param numk resulting number of incident beam directions
!
!> @date 12/08/13  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine Calckvectortrace(k,ga,ktx,kty,kt2x,kt2y,ktrad,ktstep,numk)

use io
use error
use diffraction
use crystal
use kvectors
use dynamical

IMPLICIT NONE

integer(kind=irg),INTENT(IN)        :: k(3)
integer(kind=irg),INTENT(IN)        :: ga(3)
real(kind=sgl),INTENT(IN)           :: ktx
real(kind=sgl),INTENT(IN)           :: kty
real(kind=sgl),INTENT(IN)           :: kt2x
real(kind=sgl),INTENT(IN)           :: kt2y
real(kind=sgl),INTENT(IN)           :: ktrad
integer(kind=irg),INTENT(IN)        :: ktstep
integer(kind=irg),INTENT(OUT)       :: numk

integer                             :: istat,j
real                                :: kr(3),glen,delta,kstar(3),kt(3),gan(3),gperp(3),ktlen, dktx, dkty

!write (*,*) 'ktx,kty: ', ktx, kty, kt2x, kt2y

! compute geometrical factors 
 glen = CalcLength(float(ga),'r')              ! length of ga
 gan = ga/glen                                 ! normalized ga
 delta = 2.0*ktrad*glen/float(2*ktstep+1)      ! grid step size in nm-1 
 call TransSpace(float(k),kstar,'d','r')       ! transform incident direction to reciprocal space
 call CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
 call NormVec(gperp,'r')                       ! normalize g_perp
 call NormVec(kstar,'r')                       ! normalize reciprocal beam vector

!write (*,*) 'Calckvectortrace: ',gan, gperp, glen, kstar 

! kt = -(klaue(1)+float(i)*delta)*gan - (klaue(2)+float(j)*delta)*gperp  ! tangential component of k

 j = 0
 if (.not.associated(khead)) then     ! allocate the head and ktail of the linked list
   allocate(khead,stat=istat)         ! allocate new value
   if (istat.ne.0) call FatalError('Calckvectortrace: unable to allocate head pointer',' ')
   ktail => khead                     ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   numk = 1                           ! keep track of number of k-vectors so far
! this should be the starting point of the line trace
!   kt = - glen * ( - ktx*gan + kty * gperp)
   kt = - glen * ( ktx*gan + kty * gperp)
!   write (*,*) j, kt
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktx**2+kty**2)         ! squared length of tangential component
   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end if
! write (*,*) 'kt start = ',kt
 
 dktx = (kt2x - ktx)/float(ktstep-1)
 dkty = (kt2y - kty)/float(ktstep-1)

! write (*,*) 'stepsizes :', dktx, dkty
 do j=1,ktstep-1
      allocate(ktail%next,stat=istat)  ! allocate new value
      if (istat.ne.0) call FatalError('Calckvectortrace: unable to allocate pointer',' ')
      ktail => ktail%next              ! ktail points to new value
      nullify(ktail%next)              ! nullify next in new value
      numk = numk + 1                  ! keep track of number of k-vectors so far
!      kt = - glen * (-(ktx+float(j)*dktx)*gan + (kty+float(j)*dkty) * gperp) ! tangential component of k
      kt = - glen * ( (ktx+float(j)*dktx)*gan + (kty+float(j)*dkty) * gperp) ! tangential component of k
!   write (*,*) j, kt
      ktail%kt = kt                    ! store tangential component of k
      ktlen = delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component
      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end do

end subroutine Calckvectortrace


!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcSgh
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute structure factor-like array for ECCI and ECP simulations
!
!> @param nn dimension of array
!> @param Sgh output array
!> @param nat normalization array
!
!> @date 03/05/14  MDG 1.0 original (used to be in-line in ECP and ECCI programs)
!> @date 03/11/14  MDG 1.1 converted to diagonal Sgh array only
!--------------------------------------------------------------------------
subroutine CalcSgh(nn,Sgh,nat)

use local
use crystalvars
use crystal
use gvectors
use constants
use symmetry

IMPLICIT NONE

integer(kind=irg),INTENT(IN)		:: nn
complex(kind=dbl),INTENT(INOUT)	:: Sgh(nn,nn)
integer(kind=irg),INTENT(INOUT)	:: nat(100)

integer(kind=irg)			:: ip, ir, ic, kkk(3), ikk, n, numset
real(kind=sgl)				:: Znsq, DBWF, kkl
complex(kind=dbl)			:: carg
real(kind=dbl)                        :: ctmp(192,3),arg, tpi

  tpi = 2.D0 * cPi
  Sgh = dcmplx(0.D0,0.D0)
  numset = cell % ATOM_ntype  ! number of special positions in the unit cell

! for each special position we need to compute its contribution to the Sgh array
  do ip=1,numset
    call CalcOrbit(ip,n,ctmp)
    nat(ip) = n
! get Zn-squared for this special position, and include the site occupation parameter as well
    Znsq = float(cell%ATOM_type(ip))**2 * cell%ATOM_pos(ip,4)
! loop over all contributing reflections
! ir is the row index
    rltmpa => reflist%next    ! point to the front of the list
    do ir=1,nn
! ic is the column index
      rltmpb => reflist%next    ! point to the front of the list
      do ic=1,nn
        kkk = rltmpb%hkl - rltmpa%hkl
! We'll assume isotropic Debye-Waller factors for now ...
! That means we need the square of the length of s=  kk^2/4
        kkl = 0.25 * CalcLength(float(kkk),'r')**2
! Debye-Waller exponential times Z^2
        DBWF = Znsq * exp(-cell%ATOM_pos(ip,5)*kkl)
! here is where we should insert the proper weight factor, Z^2 exp[-M_{h-g}]
! and also the detector geometry...   For now, we do nothing with the detector
! geometry; the Rossouw et al 1994 paper lists a factor A that does not depend
! on anything in particular, so we assume it is 1. 
        do ikk=1,n
! get the argument of the complex exponential
          arg = tpi*sum(kkk(1:3)*ctmp(ikk,1:3))
          carg = dcmplx(dcos(arg),dsin(arg))
! multiply with the prefactor and add
          Sgh(ir,ic) = Sgh(ir,ic) + carg * dcmplx(DBWF,0.D0)
        end do
        rltmpb => rltmpb%next  ! move to next column-entry
      end do
     rltmpa => rltmpa%next  ! move to next row-entry
   end do  
  end do
  
end subroutine CalcSgh




!C***********************************************************************
!C
!C                        naninfchk.f
!C
!C  	*****************************************************************
!C 	* 								*
!C	* 	Absoft Corporation 					* 
!C 	*	2781 Bond Street					*
!C	*	Rochester Hills, MI  48309				*
!C	*								*
!C	*	This file contains example code for demonstration	*
!C	*	purposes only.  Absoft makes no warranty of the	* 
!C	*	suitability of this code for any purpose.		*
!C	*								*
!C	*	In no event shall Absoft be liable for any incidental,*
!C	*	indirect, special, or consequential damages arising	*
!C	*	out of the use of this code.				*
!C	*								*
!C	***************************************************************** 
!C
!C Routines to test real and double values against NaN and INF
!C
!C            NANCHK(X) - tests REAL*4 value X against NaN
!!C            DNANCHK(X) - tests REAL*8 value X against NaN
!!C            INFCHK(X) - tests REAL*4 value X against INF
!!C            DINFCHK(X) - test REAL*8 value X against INF
!C
!C For little endian machines (Intel x86), compile with
!C
!C      f77 -c -DBYTE_SWAPPED=1 naninfchk.f
!C	or
!C      f90 -c -DBYTE_SWAPPED=1 naninfchk.f -YBOZTYPE=INT
!C
!C For big endian machines (PowerPC), compile with
!C
!C      f77 -c naninfchk.f
!C	or
!C      f90 -c naninfchk.f -YBOZTYPE=INT
!C
!C***********************************************************************
RECURSIVE LOGICAL FUNCTION NANCHK(X)
IMPLICIT NONE
REAL,INTENT(IN)      :: X
REAL                 :: Y
INTEGER              :: I
EQUIVALENCE(Y,I)

Y = X
NANCHK = isnan(Y) !((I .AND. z'7f800000') .EQ. z'7f800000') .AND.((I .AND. z'007fffff') .NE. z'00000000')

RETURN
END





