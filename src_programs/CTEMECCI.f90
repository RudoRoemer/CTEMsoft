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
!> @brief Zone axis ECCI defect image computations
!
!> @date   03/07/10 MDG  1.0 original, based on STEMdefect.f90 (collaboration with OSU) 
!> @date   03/21/10 MDG  2.0 modified for STEM illumination, fast version with bi-variate interpolation
!> @date   10/15/10 MDG  3.0 modified for ECCI, using ECPz concepts
!> @date   06/02/11 MDG  3.1 verified integration direction after dislocation changes to main CTEMlib routines
!> @date   08/04/12 MDG  4.0 complete rework based on Stefan Zaeferrer's comments; defect counts only
!>                           for incident beam, not for outgoing beam
!> @date   11/19/13 MDG  5.0 complete rewrite based on new ctemlib.a library routines
!--------------------------------------------------------------------------

! is this module still needed ?
module Svars

use local

IMPLICIT NONE

complex(kind=dbl),allocatable    :: DHWM(:,:),Afirst(:,:),DHWMvoid(:,:),DDD(:,:),Sarray(:,:,:,:)
complex(kind=dbl),allocatable    :: q(:,:),qin(:,:),qout(:,:),r(:,:),amp(:),amp2(:),Azz(:,:),para(:)
complex(kind=dbl)                :: czero,cone
integer(kind=irg),allocatable    :: expval(:,:,:)
integer(kind=irg)                :: nsl
real(kind=sgl)                   :: zmax, sl
real(kind=dbl)                   :: thr
real(kind=sgl),allocatable       :: ECCIimage(:,:,:)

end module Svars



program CTEMECCI 

use local
use files
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMECCI.nml'
progname = 'CTEMECCI.f90'
call Interpret_Program_Arguments(nmldeffile,7,(/ 3, 41, 200, 201, 202, 203, 204 /) )

! perform the zone axis computations
call ECCIsimulation(nmldeffile)

end program CTEMECCI

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECCIsimulation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute ECCI images for a given range of beam directions
!
!> @note This is somewhat similar to a zone axis defect computation
!
!> @param nmlfile namelist file name
!
!> @date 11/19/13  MDG 1.0 major rewrite from older ZAECCI program
!--------------------------------------------------------------------------
subroutine ECCIsimulation(nmlfile) 

use local
use crystalvars 
use crystal
use symmetryvars
use symmetry
use constants
use diffraction
use dynamical
use gvectors
use kvectors
use error
use files
use io
use foilmodule
use stacking_fault
use dislocation
use YSHModule
use void
use inclusion
use defectmodule
use timing
use STEMmodule
use Svars

IMPLICIT NONE

character(fnlen),INTENT(IN)	 :: nmlfile

integer(kind=irg)              :: nn,i,j,ik,npix,npiy,numvoids,numdisl,numYdisl,numsf, &
                                  numinc,dinfo,t_interval,DF_nums_new,io_int(6),npx, npy, nkt, ijmax, &
                                  DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, NTHR, isym, ir, ga(3), gb(3), nbeams, &
                                  ic,numd,ix,iy, numk, k(3), fn(3), istart,ig, gg(3), ixp, iyp, skip, dgn, pgnum, SETNTHR

integer(kind=irg),parameter    :: numdd=180
real(kind=sgl)                 :: arg,thick, X(2),gdotR,xgp,DF_gf(3), DM(2,2), DD, voltage, c(3), gx(3), gy(3), dmin,&
                                  gac(3), gbc(3), io_real(3), ktmax, dkt, bragg, thetac, galen, qx, qy

character(fnlen)               :: sgname, voidname, dislname(3*maxdefects), sfname(maxdefects), dataname, &
                                  incname, dispfile,  dislYname(3*maxdefects), xtalname, foilnmlfile
character(4)                   :: dispmode, method

complex(kind=dbl)              :: dx,dy,dxm,dym
real(kind=sgl),allocatable     :: sgarray(:,:),thickarray(:)
real(kind=sgl),allocatable     :: disparray(:,:,:,:),imatvals(:,:)
real(kind=dbl),allocatable     :: att(:)
real(kind=dbl)                 :: g3n(3),g3(3),H,FNg(3),ll(3),lpg(3),gplen,LC3,sgdenom,exer

namelist / ECCIlist / stdout, DF_L, DF_npix, DF_npiy, DF_slice, k, dkt, ktmax, foilnmlfile, sgname, numvoids, &
                      voidname, numinc, incname, numdisl, dislname, numYdisl, dislYname, numsf, sfname, dinfo, &
                      dispfile, dispmode, SETNTHR, dataname, xtalname, dmin, voltage


! here we read the general simulation information from a namelist file
! first we define the default values
 xtalname = 'undefined' ! input crystal structure file
 stdout = 6          ! output channel (screen)
 SETNTHR = 1          ! number of threads for OpenMP sections
 DF_L = 1.0          ! edge length of column in nanometers
 DF_npix = 256       ! number of image pixels along x
 DF_npiy = 256       ! number of image pixels along y 
 DF_slice = 1.0      ! slice thickness in nanometers
 voltage = 30000.0   ! microscope voltage
 k = (/ 0, 0, 1 /)   ! zone axis orientation
 dkt = 0.1           ! step size along kt   
 ktmax = 5.0         ! max tangential wave vector in units of |g_a|
 dmin = 0.04         ! minimum d-spacing
 dinfo = 0           ! switch to make makedislocation verbose
 foilnmlfile = 'FOIL_rundata.nml'	! name of the foil rundata namelist file
 sgname = 'nofile'   ! if this variable is different from 'nofile', then an external sg array is read (to be implemented)
 numdisl = 0         ! number of dislocation files
 numYdisl = 0        ! number of relaxed dislocations
 numsf = 0           ! number of stacking fault files
 numinc = 0          ! number of inclusions
 numvoids = 0        ! number of voids
 voidname = 'none'   ! filename for void data
 dislname = ''       ! filenames for dislocation data
 sfname = ''         ! filenames for stacking fault data
 incname = 'none'    ! filename for inclusion data
 dataname = 'test.data'! default root name for output files
 dispfile = 'none'   ! name of the displacement field output file (will be created if different from none)
 dispmode = 'not'    ! should a diplacement file be written ('new') or read ('old') or neither ('not')?
 t_interval = 5
 
! then we read the rundata namelist, which may override some of these defaults  
 OPEN(UNIT=dataunit,FILE=trim(nmlfile),DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=ECCIlist)
 CLOSE(UNIT=dataunit)

if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMECCI:',' structure file name is undefined in '//nmlfile)
end if

! display the standard program info
 progname = 'CTEMECCI.f90'
 progdesc = 'Near zone axis ECCI defect image simulation (version 3)'
 call CTEMsoft

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 hexset = .FALSE.
 call CrystalData(xtalname)

! initialize the wave length and lattice potential computations
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

! generate all atom positions
 call CalcPositions('v')

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
 dgn = GetPatternSymmetry(k,j,.TRUE.)
 pgnum = j
 isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
 call ShortestG(k,ga,gb,isym)
 io_int(1:3)=ga(1:3)
 io_int(4:6)=gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")
 
 ! initialize the HOLZ geometry type
 call GetHOLZGeometry(float(ga),float(gb),k,fn) 

! determine range of incident beam directions
 bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5
 thetac = (ktmax * 2.0 * bragg) ! *1000.0
 
! the number of pixels across the illumination disk is equal to 2*npix + 1
 npx = DF_npix
 npy = DF_npiy
 io_int(1:2) = (/ 2*npx + 1, 2*npy+1 /)
 call WriteValue('Output image size = ', io_int, 2, "(I6i,I6)")
 mess=' '; call Message("(A/)")

 isym = 1
 nkt = nint(ktmax / dkt)        ! number of pixels across half the illumination disk
 ijmax = nkt**2
 call Calckvectors(dble(k),dble(ga),dble(ktmax),nkt,nkt,numk,isym,ijmax,'Conical')    ! here we figure out how many beams there are
 io_int(1)=numk
 call WriteValue('Total # incident beam directions = ', io_int, 1, "(I8)")

! determine the reflection list
 method = 'ALL'
 call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,0,thetac)
 galen = CalcLength(float(ga),'r')

 numd = numdd
 czero=cmplx(0.0,0.0,dbl)
 cone=cmplx(1.0,0.0,dbl)
 
! force dynamical matrix routine to read new Bethe parameters from file
! in the present implementation, we're not splitting off weak reflections at all
 call Set_Bethe_Parameters(.TRUE.)
 BetheParameter%weakcutoff = BetheParameter%cutoff

! next, we read the foildata namelist from the foil nml file
! this includes material property data, in this case the elastic moduli,
! and the foil normal
 call read_foil_data(foilnmlfile,DF_npix,DF_npiy,DF_L,dinfo)
 fn = foil%F
 
! transform the foil normal
 call TransSpace(float(fn),DynFN,'d','r')
 call NormVec(DynFN,'r')

! then we need to prune the reflection list to have only reflections that will actually occur in the computation
 mess = ' Pruning reflection list (this takes a while ...) '
 call Message("(A)")
 write (*,*) numk, nbeams, 'entering pruning routine'
 call Prune_ReflectionList(numk,nbeams)
 io_int(1) = nbeams
 call WriteValue('Maximum number of contributing beams  : ', io_int, 1, '(I)')
 nn = nbeams

! in this version of the program, we take all beams in to account all the time;
! when there's time, this needs to be modified substantially to include the Bethe
! potential approximation !

! reset the incident beam direction to the actual beam direction
 ktmp => khead
 call TransSpace(ktmp%k,foil%B,'r','d')

! define the foil thickness, attenuation, and number slices per column
thick = foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
DF_slice = 1.0
DF_nums = nint(thick/DF_slice)  ! this is the number of slices for each column

! in the first version of this program, we would do a Kossel depth profile computation
! at this point.  After a long conversation with Stefan Zaeferrer during the 2012 M&M
! conference in Phoenix, I decided to change the code so that the elastic scattering
! phenomena affect only the incident electrons, not the outgoing ones (since those
! are incoherently detected on a large area detector).  In the first version, we did not
! take into account scattering by the defect(s) for the incident electron beam, but only 
! for the exiting electrons, which is the equivalent of using a point detector in an 
! EBSD experiment.  So, this second version of the program takes the incident beam
! and integrates it through the defect(s) layer by layer, and adds together the 
! intensities for all depths upto some maximum depth of several tens of nanometers.
! This maximum depth could still be computed using a Kossel style computation, but 
! we are no longer using the depth profile as weight factors, as we did in the first
! version of this program.

! To really do things right, we need to perform a Monte Carlo simulation to determine 
! the incoherent background of the ECCI pattern, and then the coherent contrast information
! on top of that background...  that will be for another day...

! this is taken from the systematic row program, so we need to verify that it is still correct!
! in this program, we're not allowing for a bent foil of any kind, but we still need to set 
! these parameters so that the CalcR routine will know not to include it.
 DF_gf = float(ga)
 DF_gstar = DF_gf/CalcLength(DF_gf,'r')**2            ! define G* such that G.G* = 1
 call TransSpace(DF_gf,DF_gc,'r','c')                 ! convert to Cartesian reference frame
 foil%sg = 0.0
 DF_gstar = 0.0

! allocate and initialize DF_Sarray, and DF_Svoid   TO BE CHECKED
 allocate(DF_Svoid(nn,nn))
 DF_Sarray = czero
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
    if (rltmpa%num.eq.0) then 
      do 
       rltmpa => rltmpa%next    ! keep going until the first contributing reflection
       if (rltmpa%num.ne.0) exit
      end do
    end if
   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
   do ic=1,nn
      if (rltmpb%num.eq.0) then 
       do 
        rltmpb => rltmpb%next    ! keep going until the first contributing reflection
        if (rltmpb%num.ne.0) exit
       end do
      end if
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

  
! compute normal absorption array for all thicknesses
  allocate (att(1:DF_nums),thickarray(1:DF_nums))
  do i=1,DF_nums
    thickarray(i) = dble(i)*DF_slice
  end do
  att = exp(-2.D0*cPi*thickarray*xgp)
  


! THIS PART NEEDS TO MOVE TO LATER IN THE PROGRAM !!!!

! loop over all reflections to get the appropriate expval powers
allocate(expval(2,nn,nn))
expval = 0.0
 rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
  do ir=1,nn
    if (rltmpa%num.eq.0) then 
      do 
       rltmpa => rltmpa%next    ! keep going until the first contributing reflection
       if (rltmpa%num.ne.0) exit
      end do
    end if   
   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
   do ic=1,nn
      if (rltmpb%num.eq.0) then 
       do 
        rltmpb => rltmpb%next    ! keep going until the first contributing reflection
        if (rltmpb%num.ne.0) exit
       end do
      end if
      if (ic.ne.ir) expval(1:2,ir,ic) = rltmpa%nab(1:2)-rltmpb%nab(1:2)
      rltmpb => rltmpb%next  ! move to next column-entry
     end do
     rltmpa => rltmpa%next   ! move to next row-entry
  end do



! next, deal with all the defects; this is the same sequence in several separate programs,
! so we really ought to make this its own module, or include it all in the defectmodule.... TO BE COMPLETED
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
   if (numsf.gt.0) call read_stacking_fault_data(numsf,numdisl,sfname,DF_L,DF_npix,DF_npiy,DF_g,dinfo)

! is there an inclusion data file? if so, then read it
   if (incname.ne.'none') call read_inclusion_data(numinc,incname,DF_L,DF_npix,DF_npiy,dinfo)

! transform the g-vector to the defect reference frames (needed for all dislocations in CalcR).
! this can only be done AFTER all dislocations and stacking faults have been created.  
!
!! This may need to be looked at for the zone axis case!  Which DF_gc should we be using here?   CHECK THIS !!!
!!
!   do i=1,numdisl
!     DF_gd(0:2,i) = matmul(DL(i)%a_dc,DF_gc)
!   end do
   
! precompute ALL the defect columns and, if needed, store them in dispfile
! this portion should be carried out in multi-threaded mode as much as possible
  allocate(disparray(2,DF_nums,DF_npix,DF_npiy),imatvals(2,DF_nums))
  disparray = 0.0; imatvals = 0
  mess = 'displacement field computation'; call Message("(A)")
  
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
    call Time_report(10*t_interval*float(NTHR)/float(DF_npix))
    call Time_start
!    write (*,*) TID,': completed timing initialization'
  end if
  allocate(DF_R(DF_nums,3))     ! each thread has its own DF_R array 

 call TransSpace(float(ga),gac,'r','c')
 call TransSpace(float(gb),gbc,'r','c')
 
!!$OMP DO SCHEDULE (GUIDED)
  do i=1,DF_npix  
    do j=1,DF_npiy
      DF_R = 0.0
! compute the displacement vectors DF_R for all points in the column
      call CalcR(i,j,numvoids,numdisl,numYdisl,numsf,numinc)
! loop over the fixed thickness slices
      do ik=1,DF_nums
! then convert to the dot-product 
       if (DF_R(ik,1).eq.-10000.0) then  ! this is point inside a void
 	imatvals(1:2,ik) = -10000
       else  ! it is not a void, so use the full dot product g.R (all vectors must be Cartesian !)
! use gac and gbc to get the two dot products and store both of them as integers mapped onto the 0..180 or 0..360 range
         gdotR = Dot_Product(gac,DF_R(ik,1:3))
         imatvals(1,ik) = numd*amod(gdotR+1000.0,1.0)
         gdotR = Dot_Product(gbc,DF_R(ik,1:3))
         imatvals(2,ik) = numd*amod(gdotR+1000.0,1.0)
       end if
     end do ! k loop
     disparray(1:2,1:DF_nums,i,j) = imatvals(1:2,1:DF_nums)
   end do
  if ((mod(i,10*t_interval).eq.0).and.(TID.eq.0)) call Time_remaining(i,DF_npix)
end do
!!$OMP END DO 
 if (TID.eq.0) call Time_stop(DF_npix*DF_npiy)
!!$OMP END PARALLEL
 
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

end if

! ok, all the set up is now complete;
! next, we prepare for the actual image simulation
  npix = DF_npix
  npiy = DF_npiy
  allocate(ECCIimage(npix,npiy,numk))
  ECCIimage = 0.0
  

 ! determine all the excitation errors for all reflections and beam directions
 mess = 'Computing excitation error array' ; call Message("(A)")
! distance between consecutive HOLZ layers in nm-1
 H = 1.0/CalcLength(float(k),'d')
! determine g3 basis vector (vector normal to ZOLZ, with length H)
 call CalcCross(dble(ga),dble(gb),g3,'r','r',1)
 call NormVec(g3,'r')
 g3n = g3
 g3 = H * g3
 FNg = (/ CalcDot(dble(DynFN),dble(ga),'r'), CalcDot(dble(DynFN),dble(gb),'r'), CalcDot(dble(DynFN),dble(g3),'r') /)
 allocate(sgarray(1:nn,1:numk))
 ktmp => khead
! open(unit=20,file=glistname,status='unknown',action='write',form='formatted')  
 do ik=1,numk
  ! get the Laue center (x and y components of the wave vector in units of ga and gb)
    ll = ktmp%kt  !  (1) * float(ga) + tmp%kt(2) * float(gb)
    rltmpa => reflist%next
    reflectionloop: do ig=1,nn
      if (rltmpa%num.eq.0) then 
        do 
         rltmpa => rltmpa%next    ! keep going until the first contributing reflection
         if (rltmpa%num.ne.0) exit
        end do
      end if   
!    if (ik.eq.1) write (20,*) ig,rltmpa%hkl,CalcDiffAngle(rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3))
      gg = rltmpa%hkl
      lpg = ll+gg
      gplen = CalcLength(lpg,'r')
      LC3 = sqrt(1.D0-mLambda**2*CalcLength(ll,'r')**2)
      if (gplen.eq.0.D0) then
        exer = -mLambda*CalcDot(dble(gg),2.D0*ll+gg,'r')/(2.D0*LC3*CalcDot(g3n,dble(DynFN),'r'))  
        sgdenom = 0.0
      else
        sgdenom = 2.D0*LC3*CalcDot(g3n,dble(DynFN),'r')-2.D0*mLambda*CalcDot(lpg,dble(DynFN),'r')
        exer = -(mLambda*CalcDot(dble(gg),2.D0*ll+gg,'r')-2.D0*LC3*CalcDot(g3n,lpg,'r'))/sgdenom
      end if
      sgarray(ig,ik) = exer
      rltmpa => rltmpa%next
    end do reflectionloop
    ktmp => ktmp%next
end do
!close(unit=20,status='keep')


! define the numd complex defect parameters
  allocate(para(0:numd))
  do i=0,numd
    arg = 2.D0*cPi*dble(i)/dble(numd)
    para(i) = cmplx(cos(arg),-sin(arg),dbl)
  end do

! store necessary data in data file
  mess ='Storing data for IDL visualization program in '//dataname
  call Message("(A/)")
  open(UNIT=dataunit,FILE=trim(dataname),STATUS='unknown',FORM='unformatted')  
  
! program mode
  write (dataunit) 'ECCI'

! filename of corresponding data set
  write (dataunit) dataname

! various bits of useful information
  write (dataunit) xtalname						! crystal structure file names	
  write (dataunit) k							! incident wave vector
  write (dataunit) CalcDiffAngle(ga(1),ga(2),ga(3))*0.5		! Bragg angle for the first reflection (whether allowed or not)
  write (dataunit) nkt				                       ! number of pixels along disk radius
  write (dataunit) sngl(mLambda)					! wave length
  write (dataunit) thetac						! beam divergence
  write (dataunit) DF_L						! pixel size
	
! number of reflections, and associated information (hkl, ...)
  call TransSpace(float(k),c,'d','c')
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

  write (dataunit) DF_npix,DF_npiy,numk


! initialize the timer
  numstart = 1
  numstop = numk
  call Time_report(float(t_interval)/float(numk))
  call Time_start

write (*,*) 'Shape(ECCIimage) = ',shape(ECCIimage)

!--------------------------------------------------------------
!--------------------------------------------------------------
!--------------------------------------------------------------
mainloop: do isg = numstart,numstop   ! this is the main computational loop over the wave vectors

! the original lines for the computation of the Sarray are now in the CalcSEntry subroutine ...
! we do need to reset the Sarray to -10000.0
  allocate(Sarray(nn,nn,0:numd,0:numd))
  Sarray(1,1,0:numd,0:numd) = -10000.0

  forall (i=1:nn)
   DHWMz(i,i)=2.0*cPi*cmplx(0.0,sgarray(i,isg))    ! initialize the diagonal elements of the dynamical matrix
   DHWMvoid(i,i) = DHWMz(i,i)
  end forall
  
  allocate(Afirst(nn,nn))
  Afirst = czero
  forall (i=1:nn)
    Afirst(i,i) = cone
  end forall
  
!---------------------------------------------------------!
! Finally, here it is: the actual ECCI image computation  !
!---------------------------------------------------------!
  
! do we have to start at the bottom slice ?
 istart = 1

! allocate auxiliary arrays
 allocate(Azz(nn,nn),amp(nn),amp2(nn))

! loop over the image pixels
 donpix: do i=1,npix
 donpiy:   do j=1,npiy
    amp = czero
    amp(1) = cone
doslices: do ik=1,DF_nums    ! loop over the fixed thickness slices
! compute the appropriate scattering matrix to propagate with (see section 8.3.3 in the book)
       if (disparray(1,ik,i,j).eq.-10000) then  ! this is point inside a void
 	 Azz(:,:) = DF_Svoid(:,:)    ! so we use the void propagator matrix
       else  ! it is not a void
! in this version, we use complex bi-variate interpolation to select the appropriate Azz 
! matrix from the Sarray; if the required S matrix is not yet known, then we compute it first
         ix = int(disparray(1,ik,i,j))
         iy = int(disparray(2,ik,i,j))
         dx = cmplx(amod(disparray(1,ik,i,j),1.0),0.0)
         dy = cmplx(amod(disparray(2,ik,i,j),1.0),0.0)
         dxm = cone-dx
         dym = cone-dy
         if (minval( (/ dble(Sarray(1,1,ix,iy)),dble(Sarray(1,1,ix+1,iy)),dble(Sarray(1,1,ix,iy+1)), &
                dble(Sarray(1,1,ix+1,iy+1)) /) ).eq.-10000.0) call CalcSEntry(nn,ix,iy,DF_slice)
              ixp = ix+1
              iyp = iy+1
	      Azz = dxm*dym*Sarray(1:nn,1:nn,ix,iy)+dx*dym*Sarray(1:nn,1:nn,ixp,iy)+ &
                        dxm*dy*Sarray(1:nn,1:nn,ix,iyp)+dx*dy*Sarray(1:nn,1:nn,ixp,iyp)
         end if
! and multiply with this matrix, to integrate over the thickness
         amp2 = matmul(Azz,amp)
! then compute all the intensities at this depth, and add them all together ... 
         ECCIimage(i,j,isg) = ECCIimage(i,j,isg) + att(ik) * sum(cabs(amp2)**2)
! update the amplitudes to the next slice
         amp = amp2
      end do doslices ! loop over slices 
      
    end do donpiy
    if (mod(i,32).eq.0) write (*,*) isg,i,' done'
  end do donpix
  deallocate(Afirst,Azz,amp,amp2,Sarray)
  call Time_remaining(isg-numstart+1,numstop-numstart+1)
  
! do an intermediate save of the resulting image intensities  
!if (intname.ne.'none') then
!  write (*,*) 'intermediate image stored in intermediate.data'
!  open(unit=20,file=intname,status='unknown',action='write',form='unformatted')  
!  write (20) npix,npiy 
!  write (20) ECCIimage
!  close(unit=20,status='keep')
!end if
  
end do mainloop

! and save the data 
mess = 'data written to file '//trim(dataname)
call Message("(/A/)")
write (dataunit) ECCIimage
close(UNIT=dataunit,STATUS='keep')
  
    
! ok, so the main computation is complete; print some timing information
call Time_stop(npix*npiy)

end subroutine ECCIsimulation


!
!
!! ###################################################################
!! 
!!  subroutine Calckvectorcone
!!
!!  Author: Marc De Graef
!!  
!!  Description: computes the independent incident beam directions
!!  for the multi-beam case and returns them as a linked list in the
!!  global variables head and tail
!! 
!!  History
!! 
!!  modified by  rev reason
!!  -------- --- --- -----------
!!   6/4/01  MDG 1.0 original
!! ###################################################################
!subroutine Calckvectorcone(k,ga,ktx,kty,ktrad,ktstep,numk)
!
!use io
!use error
!use diffraction
!use crystal
!use dynamical
!
!IMPLICIT NONE
!
!integer              :: numk,k(3),ga(3),istat,imin,imax,jmin,jmax,ijmax,i,j,ktstep,ic,jc,ki
!real                 :: ktrad,ktx,kty,kr(3),glen,delta,kstar(3),kt(3),gan(3),gperp(3),ktlen
!
!intent(IN)           :: ktx,kty,k,ga,ktrad,ktstep
!intent(OUT)          :: numk
!
!! compute geometrical factors 
! glen = CalcLength(float(ga),'r')              ! length of ga
! gan = ga/glen                                 ! normalized ga
! delta = 2.0*ktrad*glen/float(2*ktstep+1)             ! grid step size in nm-1 
! call TransSpace(float(k),kstar,'d','r')       ! transform incident direction to reciprocal space
! call CalcCross(float(ga),kstar,gperp,'r','r',0)      ! compute g_perp = ga x k
! call NormVec(gperp,'r')                       ! normalize g_perp
! call NormVec(kstar,'r')                       ! normalize reciprocal beam vector
! write (*,*) 'Calckvectors : ',k,ga,glen,ktrad,ktstep,ktx,kty,delta
! write (*,*) 'Calckvectors : ',kstar,gperp
!
! open(unit=20,file='ECPzkt.txt',status='unknown',form='formatted')
!
!! deal only with the incident beam (parallel illumination)
!if (ktstep.eq.0) then
! if (.not.associated(head)) then     ! allocate the head and tail of the linked list
!   allocate(head,stat=istat)         ! allocate new value
!   if (istat.ne.0) call FatalError('Calckvectorcone: unable to allocate head pointer',' ')
!   tail => head                      ! tail points to new value
!   nullify(tail%next)                ! nullify next in new value
!   numk = 1                          ! keep track of number of k-vectors so far
! ! this should be the center vector of the illumination cone !!!
!   kt = - glen * (ktx*gan + kty * gperp)
!   tail%kt = kt                           ! store tangential component of k
!   ktlen = glen**2*(ktx**2+kty**2)         ! squared length of tangential component
!   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
!   tail%k = kr                            ! store in pointer list
!   tail%kn = CalcDotd(tail%k,dble(kstar),'r')    ! normal component of k
!  write(unit=20,*) tail%kt,tail%k,CalcLengthd(tail%k,'r')
! end if
!else
!! next, put the center of the cone in units of (i,j) (original ECP "screen" coordinates)
!  ic = int(ktx*glen/delta)
!  jc = int(kty*glen/delta)
!  ki = ktstep
!  write(*,*) ic,jc,ki,1.0/mLambda
!
!
! if (.not.associated(head)) then     ! allocate the head and tail of the linked list
!   allocate(head,stat=istat)         ! allocate new value
!   if (istat.ne.0) call FatalError('Calckvectorcone: unable to allocate head pointer',' ')
!   tail => head                      ! tail points to new value
!   nullify(tail%next)                ! nullify next in new value
!   numk = 1                          ! keep track of number of k-vectors so far
! ! this should be the center vector of the illumination cone !!!
!   tail%i = ic                            ! i-index of beam
!   tail%j = jc                            ! j-index of beam
!   kt = -float(tail%i)*delta*gan - float(tail%j)*delta*gperp  ! tangential component of k
!   tail%kt = kt                           ! store tangential component of k
!   ktlen = delta**2*(tail%i**2+tail%j**2)         ! squared length of tangential component
!   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
!   tail%k = kr                            ! store in pointer list
!   tail%kn = CalcDotd(tail%k,dble(kstar),'r')    ! normal component of k
!  write(unit=20,*) tail%i,tail%j,tail%kt,tail%k,CalcLengthd(tail%k,'r')
! else
!   call FatalError('Calckvectorcone: pointer head already allocated',' ')
! end if
!
!! the following lines are quite different if symmetry is taken into account;
!! check the MBsym.f90 program to determine how that can be done.
!  imin =  -ki; imax = ki; jmin = -ki; jmax = ki; 
!  ijmax = ki**2
!! now do the real work
!  do i=imin,imax
!   do j=jmin,jmax
!    if (.not.((i.eq.0).and.(j.eq.0))) then  ! the point (0,0) has already been taken care of
!     if ((i**2+j**2).le.ijmax) then   ! is point inside the incident cone ?
!      allocate(tail%next,stat=istat)  ! allocate new value
!      if (istat.ne.0) call FatalError('Calckvectors: unable to allocate pointer',' ')
!      tail => tail%next               ! tail points to new value
!      nullify(tail%next)              ! nullify next in new value
!      numk = numk + 1                 ! keep track of number of k-vectors so far
!      tail%i = ic+i                   ! i-index of beam
!      tail%j = jc+j                   ! j-index of beam
!      kt = -float(tail%i)*delta*gan - float(tail%j)*delta*gperp  ! tangential component of k
!      tail%kt = kt                    ! store tangential component of k
!      ktlen = delta**2*(tail%i**2+tail%j**2)         ! squared length of tangential component
!      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
!      tail%k = kr                     ! store in pointer list
!      tail%kn = CalcDotd(tail%k,dble(kstar),'r')    ! normal component of k
!  write(unit=20,*) tail%i,tail%j,tail%kt,tail%k,CalcLengthd(tail%k,'r')
!!     write (*,*) i,j,tail%kt,tail%kn
!     end if
!    end if
!   end do
!  end do
!end if
!close(unit=20,status='keep')
!  
!end subroutine
!
!




!--------------------------------------------------------------------------
!
! SUBROUTINE:CalcSEntry
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute entries in the scattering matrix array
!
!> @note Here we precompute four entries in the array of scattering matrices that can 
! then be used, either directly, or via bi-linear interpolation,
! by the image computation portion of this program.
!
! For starters, we'll subdivide the range of possible alpha values
! in 180 segments (2 degrees each), with a copy of the last one
! (i.e., 181x181 = 32761 entries or 240 Mb for 31 beams)
!
! this part is essentially the same as the threaded section of the older
! ZAdefect.all.f90 program (which was a test program).

!> @param nn number of beams
!> @param ini x-parameter
!> @param inj y-parameter
!
!> @date 11/19/13  MDG 1.0 rewrite including MatrixExponential call
!--------------------------------------------------------------------------
recursive subroutine CalcSEntry(nn,ini,inj,DF_slice)

use local
use dynamical
use Svars
use math

IMPLICIT NONE

integer(kind=irg),INTENT(IN)        :: nn
integer(kind=irg),INTENT(IN)        :: ini
integer(kind=irg),INTENT(IN)        :: inj
real(kind=sgl),INTENT(IN)           :: DF_slice

integer(kind=irg)     :: i,j,ir,ic

allocate(DDD(nn,nn))

! loop over the four points closest to the one requested
do i=ini,ini+1
 do j=inj,inj+1
  if (Sarray(1,1,i,j).eq.-10000.0) then 
! loop over all reflections in the array DDD using the information in expval
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
  
  end if
 end do
end do

deallocate(DDD)

end subroutine CalcSEntry
