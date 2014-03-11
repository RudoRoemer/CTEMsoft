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
! CTEMsoft2013:CTEMECP.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMECP 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Zone axis electron channeling patterns
!
!> @date 03/18/10 MDG 1.0 f90
!> @date 08/09/10 MDG 2.0 corrected weight factors and g-vector ordering problem
!> @date 11/18/13 MDG 3.0 major rewrite with new libraries 
!--------------------------------------------------------------------------
program tryECCI2

use local
use files
use io

IMPLICIT NONE

character(fnlen)	:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMECCI.nml'
progname = 'CTEMECCI.f90'
call Interpret_Program_Arguments(nmldeffile,8,(/ 0, 3, 41, 200, 201, 202, 203, 204 /) )

! initialize all user-defined variables
call ComputeECCI(nmldeffile)

end program tryECCI2


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
!> @date 03/05/14  MDG 3.4 rewrite, starting from ECP program
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

character(fnlen),INTENT(IN)	       :: nmlfile

integer(kind=irg)    		        :: nn,i,j,k,npix,npiy,ii,jj,numvoids,numdisl, numset, dgn, maxHOLZ, pgnum, ia, ib, &
					numYdisl,numsf,numinc,numapb,dinfo,t_interval,nat(100),kkk(3), &
					DF_nums_new,DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, &
					NTHR, isym, ir, ga(3), gb(3),kk(3),ic,g,ix,iy,nkt,nbeams, ik, ig, &
					numk,ixp,iyp,SETNTHR, io_int(6), skip, gg(3), iSTEM, nktstep, ip, n, ikk
!                                  OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
integer(kind=irg),parameter 		:: numd = 3600
real(kind=sgl)         		:: thick, X(2), dmin, dkt, bragg, thetac, frac, &
					lauec(2),lauec2(2),gdotR,DF_gf(3),Znsq,tpi, &
					DM(2,2), DD, ll(3), c(3), gx(3), gy(3), DBWF, &
					gac(3), gbc(3),zmax, ktmax, io_real(2), kkl, voltage, ijmax
real(kind=dbl)                        :: ctmp(192,3),arg
character(fnlen)      			:: dataname,sgname,voidname,dislname(3*maxdefects),sfname(maxdefects),ECPname, &
					incname,dispfile,xtalname,foilnmlfile, STEMnmlfile,dislYname(3*maxdefects), apbname
character(4)            		:: dispmode, summode
character(3)				:: method
character(5)                          :: progmode
complex(kind=dbl),allocatable    	:: DHWM(:,:),DHWMvoid(:,:),DDD(:,:),Sarray(:,:,:,:)
complex(kind=dbl),allocatable    	:: amp(:),amp2(:),Azz(:,:)
complex(kind=dbl)                	:: czero,cone,carg
complex(kind=dbl)                	:: para(0:numd),dx,dy,dxm,dym, xgp
real(kind=sgl),allocatable       	:: sgarray(:,:)
real(kind=sgl),allocatable    		:: disparray(:,:,:,:),imatvals(:,:), ECCIimages(:,:)
integer(kind=sgl),allocatable    	:: expval(:,:,:)
complex(kind=dbl),allocatable         :: Lgh(:,:),Sgh(:,:)
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
 summode = 'full'                     ! 'full' for complete summation, 'diag' for diagonal only

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
 
 frac = 0.05
 

!! set the input parameters to default values (except for xtalname, which must be present)
!xtalname = 'undefined'		        ! initial value to check that the keyword is present in the nml file
!stdout = 6			        ! standard output
!voltage = 30000.0		        ! acceleration voltage [V]
!k = (/ 0, 0, 1 /)		        ! beam direction [direction indices]
!fn = (/ 0, 0, 1 /)		        ! foil normal [direction indices]
!dmin = 0.025			        ! smallest d-spacing to include in dynamical matrix [nm]
!distort = .FALSE.                      ! distort the input unit cell ?  
!abcdist = (/ 0.4, 0.4, 0.4/)           ! distorted a, b, c [nm]
!albegadist = (/ 90.0, 90.0, 90.0 /)    ! distorted angles [degrees]
!ktmax = 0.0                            ! beam convergence in units of |g_a|
!thetac = 0.0                           ! beam convergence in mrad (either ktmax or thetac must be given)
!startthick = 2.0		        ! starting thickness [nm]
!thickinc = 2.0			        ! thickness increment
!numthick = 10			        ! number of increments
!npix = 256			        ! output arrays will have size npix x npix
!outname = 'ecp.data'        	        ! output filename
!compmode = 'Bloch'                     ! 'Blochwv' or 'ScatMat' solution mode (Bloch is default)
!zintstep = 1.0                        ! integration step size for ScatMat mode

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)

! initialize the wave length and lattice potential computations
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

! generate all atom positions
 call CalcPositions('v')

! next, we read the foildata namelist from the SRdef_foildata.nml file
! [yes, we're using the same file as for the systematic row case]
! this includes material property data, in this case the elastic moduli,
! and the foil normal, which we will need in the next step
  call read_foil_data(foilnmlfile,DF_npix,DF_npiy,DF_L,dinfo)
  DynFN = foil%F
  call NormVec(DynFN,'r')

! define the foil thickness, attenuation, and number slices per column
  thick = foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
  DF_nums = nint(thick/DF_slice)  ! this is the number of slices for this particular column
write (*,*) 'foil thickness', foil%zb, thick, DF_nums


! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
 dgn = GetPatternSymmetry(kk,j,.TRUE.)
 pgnum = j
 isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
 call ShortestG(kk,ga,gb,isym)
 io_int(1:3)=ga(1:3)
 io_int(4:6)=gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

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

  ktmp => khead
  do ik=1,numk
    write (*,*) ktmp%k, ktmp%kt
    ktmp => ktmp%next
  end do

! construct the list of all possible reflections
  method = 'ALL'
  maxHOLZ = 0  ! not used for method = 'ALL'
  call Compute_ReflectionList(dmin,kk,ga,gb,method,.FALSE.,maxHOLZ)

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(.TRUE.)

! and prune this list, based on the list of k-vectors
! this reduces the lengh of the list to only those reflections that are important 
! for the current geometry.
  call Prune_ReflectionList(numk,nbeams)
  write (*,*) 'Pruned reflections list has # beams ',nbeams

! compute the projections of these potential reflections onto the zone axis basis vectors
 DM(1,1) = CalcDot(float(gb),float(gb),'c')
 DM(1,2) = -CalcDot(float(ga),float(gb),'c')
 DM(2,1) = DM(1,2)
 DM(2,2) = CalcDot(float(ga),float(ga),'c')
 DD = 1.0/(DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1))

 rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
  do ir=1,nbeams
! decompose this point w.r.t ga and gb
   X(1) = CalcDot(float(rltmpa%hkl),float(ga),'c')
   X(2) = CalcDot(float(rltmpa%hkl),float(gb),'c')
   X = matmul(DM,X) * DD
   rltmpa%nab(1:2) = int(X(1:2))
   rltmpa => rltmpa%next   ! move to next reflection entry
  end do


  nat = 0

! deal with all the potential defects
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

  
!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
! for each incident wave vector, we consider all the image columns
! point to the first beam direction
  ktmp => khead
  czero = cmplx(0.0,0.0,dbl)
  call CalcUcg((/0,0,0/))   ! get the normal absorption parameter
  DynUpz = rlp%Vpmod
  numset = cell % ATOM_ntype  ! number of special positions in the unit cell
  tpi = 2.0*cPi
! allocate space for the results
  npix = DF_npix
  npiy = DF_npiy
  
  allocate(ECCIimages(npix,npiy))
  ECCIimages = 0.0

  allocate(DF_R(DF_nums,3))     ! each thread has its own DF_R array 
  call TransSpace(float(ga),gac,'r','c')
  call TransSpace(float(gb),gbc,'r','c')

  do i=0,numd
    arg = 2.D0*cPi*dble(i)/dble(numd)
    para(i) = dcmplx(dcos(arg),-dsin(arg))
  end do



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

! number of wave vectors, tangential components, etc...
  write (dataunit) numk   ! number of wave vectors
  ktmp => khead
  do ic=1,numk
    write (dataunit) sngl(ktmp%kt(1)),sngl(ktmp%kt(2))
    ktmp => ktmp%next
  end do
  
  write (dataunit) DF_npix,DF_npiy

  call Time_report(1.0/float(numk))
  call Time_start
  
!  work through the beam direction list
  ktmp => khead
  beamloop: do ik=1,numk

  write (*,*) 'starting wave vector ',ktmp%k

! compute the dynamical matrix using Bloch waves with Bethe potentials; note that the IgnoreFoilNormal flag
! has been set to .FALSE.; if it is set to .TRUE., the computation of the ZOLZ will still be mostly correct,
! but the excitation errors of the HOLZ reflections will be increasingly incorrect with HOLZ order.  This was
! useful during program testing but should probably be removed as an option altogether...
	call Compute_DynMat('BLOCHBETHE', ktmp%k, ktmp%kt, .FALSE.)
        nn = DynNbeams
        
! this is the dynamical matrix for the perfect crystal case, which we 
! need to convert to the structure matrix A by multiplication by i pi lambda
	DynMat = DynMat * dcmplx(0.0, cPi * mLambda)
write (*,*) 'structure matrix computed ', shape(DynMat), nn
 
! next, we need to compute the Sgh matrix
        allocate(Sgh(nn,nn), DDD(nn,nn))        
	call CalcSgh(nn, Sgh, nat)
write (*,*) 'Sgh matrix computed'

! loop over all reflections to get the appropriate powers for the phase shifts
	allocate(expval(2,nn,nn), amp(nn), Azz(nn,nn), Lgh(nn,nn), amp2(nn))
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
write (*,*) 'expval matrix computed'

! then we need a loop over all the image columns; for each column we
! need to compute the displacements; then we loop over the slices, and
! for each slice we compute the defect phase shift array which we then
! multiply by the structure matrix.  Computation of the column displacements
! is done by the CalcR routine.

! single core for now; OpenMP later
write (*,*) 'starting image loop ',shape(DDD), shape(DF_R)

	donpix: do i=1,DF_npix  
	  donpiy: do j=1,DF_npiy
            DF_R = 0.0
! compute the displacement vectors DF_R for all points in the column
            call CalcR(i,j,numvoids,numdisl,numYdisl,numsf,numinc,numapb)

! initialize the wave function
   	    Lgh = czero
    	    amp = czero
    	    amp(1) = cone

! loop over the fixed thickness slices
            do k=1,DF_nums

! create the defect phase shift matrix DDD for this slice
	      DDD = dcmplx(0.D0,0.D0)
	      forall (ir=1:nn) 
		DDD(ir,ir) = dcmplx(1.0D0,0.0D0)	! the diagonal must have 1 in all entries, also covers void case
	      end forall

! then convert displacements to dot-products
              if (DF_R(k,1).ne.-10000.0) then  ! this is point inside a void
! use gac and gbc to get the two dot products and convert both of them to integers mapped onto the 0..numd range
! Due to the logarithmic singularity at a dislocation core, it is possible for gdotR to be NaN; we need
! to intercept these cases, and replace those values of gdotR by 0.0
                gdotR =  Dot_Product(gac,DF_R(k,1:3))
                if (NANCHK(gdotR)) gdotR = 0.0         
		 ia = nint( numd * amod(gdotR+1000.0,1.0) )

                gdotR = Dot_Product(gbc,DF_R(k,1:3))
                if (NANCHK(gdotR)) gdotR = 0.0         
		 ib = nint( numd * amod(gdotR+1000.0,1.0) )
! use ia and ib as indices to the para array and compute the DDD defect phase factor matrix 
	        do ic=1,nn
   	          do ir=1,nn
	           if (ic.ne.ir) then  ! exclude the diagonal
	             DDD(ir,ic) = para(ia)**expval(1,ir,ic) * para(ib)**expval(2,ir,ic)
                   end if
                 end do
                end do
              end if

! multiply DDD by the structure matrix DynMat and exponentiate
	      DDD = DDD * DynMat
  	      call MatrixExponential(DDD, Azz, dble(DF_slice), 'Pade', nn)  

! now do the usual thing and propagate the wave function to the next slice
	      amp2 = matmul(Azz, amp)

!!	      if (k.eq.1) then 
!        	Lgh = spread(amp2(1:nn),dim=2,ncopies=nn) * spread(conjg(amp2(1:nn)),dim=1,ncopies=nn)
!   	      else
!       		Lgh = Lgh + spread(amp2(1:nn),dim=2,ncopies=nn) * spread(conjg(amp2(1:nn)),dim=1,ncopies=nn)
!   	      end if
	      if (k.eq.1) then 
        	forall (ir=1:nn)
		  Lgh(ir,ir) = abs(amp2(ir))**2
		end forall
   	      else
        	forall (ir=1:nn)
		  Lgh(ir,ir) = Lgh(ir,ir) + abs(amp2(ir))**2
		end forall
   	      end if
   	      
   	      amp = amp2

	   end do  ! end of column integration

    	   ECCIimages(i,j) = sngl(real(sum( Sgh * Lgh )))
 
    	 end do donpiy
    	 if (mod(i,10).eq.0) write (*,*) 'line ',i
	end do donpix


  	ECCIimages = ECCIimages / float(DF_nums) / float(sum(nat))
  	write (dataunit) ECCIimages


       deallocate(Lgh, Sgh, amp, DDD, Azz, amp2, expval)
! select next beam direction
       ktmp => ktmp%next
          
! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
     call Time_remaining(ik,numk)
     frac = frac + 0.05
   end if  

  end do beamloop

close(UNIT=dataunit,STATUS='keep')

! ok, so the main computation is complete; print some timing information
call Time_stop(npix*npiy)

mess = 'Data stored in file '//trim(dataname)
call Message("(/A/)")


end subroutine ComputeECCI




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
!        kkk = rltmpb%hkl - rltmpa%hkl
        kkk = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)

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
!   kt = - glen * (ktx*gan + kty * gperp)
   kt = - glen * (-ktx*gan + kty * gperp)
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
   kt = float(ktail%i)*delta*gan - float(ktail%j)*delta*gperp  ! tangential component of k

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
      kt = float(ktail%i)*delta*gan - float(ktail%j)*delta*gperp  ! tangential component of k
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

! compute geometrical factors 
 glen = CalcLength(float(ga),'r')              ! length of ga
 gan = ga/glen                                 ! normalized ga
 delta = 2.0*ktrad*glen/float(2*ktstep+1)      ! grid step size in nm-1 
 call TransSpace(float(k),kstar,'d','r')       ! transform incident direction to reciprocal space
 call CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
 call NormVec(gperp,'r')                       ! normalize g_perp
 call NormVec(kstar,'r')                       ! normalize reciprocal beam vector


! kt = -(klaue(1)+float(i)*delta)*gan - (klaue(2)+float(j)*delta)*gperp  ! tangential component of k

 j = 0
 if (.not.associated(khead)) then     ! allocate the head and ktail of the linked list
   allocate(khead,stat=istat)         ! allocate new value
   if (istat.ne.0) call FatalError('Calckvectortrace: unable to allocate head pointer',' ')
   ktail => khead                     ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   numk = 1                           ! keep track of number of k-vectors so far
! this should be the starting point of the line trace
   kt = - glen * ( - ktx*gan + kty * gperp)
!   write (*,*) j, kt
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktx**2+kty**2)         ! squared length of tangential component
   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end if
 dktx = (kt2x - ktx)/float(ktstep-1)
 dkty = (kt2y - kty)/float(ktstep-1)

 do j=1,ktstep-1
      allocate(ktail%next,stat=istat)  ! allocate new value
      if (istat.ne.0) call FatalError('Calckvectortrace: unable to allocate pointer',' ')
      ktail => ktail%next              ! ktail points to new value
      nullify(ktail%next)              ! nullify next in new value
      numk = numk + 1                  ! keep track of number of k-vectors so far
      kt = - glen * (-(ktx+float(j)*dktx)*gan + (kty+float(j)*dkty) * gperp) ! tangential component of k
!   write (*,*) j, kt
      ktail%kt = kt                    ! store tangential component of k
      ktlen = delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component
      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end do

end subroutine Calckvectortrace



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
NANCHK = ((I .AND. z'7f80 0000') .EQ. z'7f80 0000') .AND.((I .AND. z'007f ffff') .NE. z'0000 0000')

RETURN
END

