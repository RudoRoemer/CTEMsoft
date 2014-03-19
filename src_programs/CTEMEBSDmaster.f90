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
! CTEMsoft2013:CTEMEBSDmaster.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMEBSDmaster
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief CTEMEBSDzE computes the energy-dependent master EBSD pattern for a given structure
!
!> @todo implement full symmetry use; implement multiple reflection output
!> or, easier perhaps, selection of one reflection;
!>
!> implement OpenMP multithreading for the actual computation part; requires modifications
!> in CTEMlib.a routines (mostly THREADPRIVATE commands in several modules)
!
!> @date  03/08/12  MDG 1.0 EBSD program for fundamental zone patterns
!> @date  08/17/12  MDG 1.1 added generalized fundamental zone for all crystal symmetries
!> @date  08/20/12  MDG 1.2 modifid FZ to use the Lambert projection of a square sampling grid
!> @date  09/06/12  MDG 1.3 added support for second setting of Laue group -3m
!> @date  11/21/12  MDG 2.0 added full Lambert projection support for both square and hexagonal grids
!>				   the older code is still available for now, but will be removed after validation
!>				   of the newer code.
!> @date  12/04/12  MDG 2.1 added support for equilateral triangle mapping; needs to be validated.
!>				   also modified the structure of the output file, so that EBSD.f90 will
!>				   know which of the inverse mapping methods it should use.
!> @date  12/10/12  MDG 3.0 expanded EBSDFZ program to include energy-dependencies from Monte Carlo
!> @date  12/12/12  MDG 3.1 test to do an actual numerical integration for the I_jk integrals, using MC profiles
!> @date  08/01/13  MDG 4.0 complete rewrite with Lambert format for MC output and new ctemlib.a routines 
!>                          also, the earlier versions would do only one energy value, whereas this new 
!>                          implementation does the complete energy-dependent master pattern
!> @date  01/27/14  MDG 4.1 continued rewrite, fixed problem with kvector list, replaced gvector routines
!>		   	     with updated routines; changed program name to CTEMEBSDmaster.
!--------------------------------------------------------------------------
program CTEMEBSDmaster

use local
use files
use io

IMPLICIT NONE

character(fnlen)	:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMEBSDmaster.nml'
progname = 'CTEMEBSDmaster.f90'
call Interpret_Program_Arguments(nmldeffile,1,(/ 21 /) )

! generate a set of master EBSD patterns
 call ComputeMasterPattern(nmldeffile)

end program CTEMEBSDmaster

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeMasterPattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute an EBSD master pattern as a function of energy
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 08/01/13  MDG 3.0 complete rewrite, eliminated old Lambert projection
!> @date 09/25/13  MDG 3.1 replaced k-vector code by kvectors module
!--------------------------------------------------------------------------
subroutine ComputeMasterPattern(nmlfile)

use symmetryvars
use symmetry
use crystalvars
use crystal
use constants
use error
use gvectors
use kvectors
use io
use local
use files
use diffraction
use multibeams
use dynamical
use timing
use Lambert

IMPLICIT NONE

character(fnlen),INTENT(IN)	:: nmlfile


real(kind=dbl)   	:: ctmp(192,3),arg
integer(kind=irg)      :: isym,i,j,ik,npx,npy,ipx,ipy,debug,iE,izz, izzmax, iequiv(2,12), nequiv, num_el, MCnthreads, & ! counters
                    	numk, & ! number of independent incident beam directions
                    	ir,nat(100),kk(3), npyhex, skip, ijmax, Esel, one, &
                    	numset,n,ix,iy,iz, io_int(6),  &
                    	istat,gzero,ic,ip,ikk     ! counters
real(kind=dbl)         :: tpi,Znsq, kkl, DBWF,  totZnsq, kin, s, sthr, Bethe_store !
real(kind=sgl) 		:: dmin, io_real(5), selE
real(kind=sgl),allocatable 	:: sr(:,:,:), srkin(:,:), srhex(:,:,:), srkinhex(:,:), EkeVs(:) ! results
complex(kind=dbl)  		:: czero
complex(kind=dbl),allocatable 	:: Lgh(:,:), Sgh(:,:)
logical 		:: usehex, switchmirror
character(fnlen) 	:: xtalname, outname
! the following will need to be moved elsewhere at some point...
integer(kind=irg),parameter	:: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees

! Monte Carlo derived quantities
integer(kind=irg)   	:: numEbins, numzbins, nsx, nsy    ! variables used in MC energy file
real(kind=dbl)     	:: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, etotal ! enery variables from MC program
integer(kind=irg),allocatable :: accum_e(:,:,:), accum_z(:,:,:,:), thick(:)
real(kind=sgl),allocatable :: lambdaE(:,:)
character(80) 		:: energyfile


namelist /EBSDmastervars/ xtalname,dmin,npx,outname,energyfile, Esel

! note that the CTEMMC.f90 program creates a statistical output file that 
! must be read by the present program, so that things like energy etc are 
! known to the program...  The content of the MC output file is as follows:
! write(dataunit) numEbins, numzbins, numsx, numsy
! write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
! write(dataunit) accum_e
! write (dataunit) accum_z
! So that means that the beam voltage range is already known; we do not need to
! ask for this again...

! read namelist file with all input parameters
! here are the default values in case some parameters are absent from the file
xtalname 	= 'undefined'	! program must read existing structure file name
dmin 		= 0.015   	! minimum d-spacing in order to be admitted to the list (must become user entry)
npx		= 500		! Nx pixels (total = 2Nx+1)
outname		= 'EBSDmasterout.data'	! default filename for final output
energyfile	= 'z0E.data' 	! default filename for z_0(E_e) data from CTEMMC Monte Carlo simulations
Esel           = -1            ! selected energy value for single energy run

! then we read the rundata namelist, which may override some of these defaults  
OPEN(UNIT=dataunit,FILE=nmlfile,DELIM='APOSTROPHE')
READ(UNIT=dataunit,NML=EBSDmastervars)
CLOSE(UNIT=dataunit)

if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMEBSDmaster:',' structure file name is undefined in '//nmlfile)
end if

! print some information
progname = 'CTEMEBSDmaster.f90'
progdesc = 'EBSD Energy-dependent Master Pattern Simulation'
call CTEMsoft

! square master pattern
npy = npx
sthr = 10.D0

! ---------- start of symmetry and crystallography section
! first get the crystal data and some other crystallographic and symmetry information 
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)

! generate all atom positions
 call CalcPositions('v')

! the following line needs to be verified ... 
!  hexset = .TRUE.    ! if hexagonal structure this switch selects between three and four index notation (4 if true)

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do
 isym = j
 
! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
! dgn = GetPatternSymmetry( (/ 0, 0, 1 /),j,.TRUE.)
! pgnum = j
! isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers
 
! and convert this to the corresponding  Laue point group number since diffraction 
! patterns are always centrosymmetric (hence, there are only 11 different cases for
! the symmetry of the incident beam).   
  isym = PGLaueinv(isym)  

! If the Laue group is # 7, then we need to determine the orientation of the mirror plane.
! The second orientation of the mirror plane is represented by "Laue group" # 12 in this program.
 switchmirror = .FALSE.
 if (isym.eq.7) then
  do i=1,11
    if (cell%SYM_SGnum.eq.LaueTest(i)) switchmirror = .TRUE.
  end do
 end if
 if (switchmirror) then
  isym = 12
  mess = ' Switching computational wedge to second setting for this space group'; call Message("(A)")
 end if
 write (*,*) ' Laue group # ',isym, PGTHD(j)

! if this point group is trigonal or hexagonal, we need to switch usehex to .TRUE. so that
! the program will use the hexagonal sampling method
usehex = .FALSE.
if (((isym.ge.6).and.(isym.le.9)).or.(isym.eq.12)) usehex = .TRUE.
! ---------- end of symmetry and crystallography section

! ---------- read Monte Carlo output file and extract necessary parameters
! first, we need to load the data from the MC program.  This is an array of integers, which
! has an energy histogram for each scintillator pixel.  These are integer values to minimize
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

! MCnthreads = 8
io_int(1:6) = (/ numEbins, numzbins, nsx, nsy, num_el, MCnthreads /)
call WriteValue(' NumEbins, numzbins, nsx, nsy, num_el, MCnthreads ',io_int,6,"(5I,',',I)")
etotal = num_el ! * MCnthreads

read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")

! allocate(accum_e(numEbins,-nx:nx,-nx:nx),accum_z(numEbins,numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)

allocate(accum_e(numEbins,-nsx:nsx,-nsy:nsy),accum_z(numEbins,numzbins,-nsx/10:nsx/10,-nsy/10:nsy/10),stat=istat)
read(dataunit) accum_e
! actually, we do not yet need the accum_e array; that is for the actual EBSD image computation program
! but we need to skip it in this unformatted file...
deallocate(accum_e)

read(dataunit) accum_z    ! we only need this array for the depth integrations

close(dataunit,status='keep')
mess = ' -> completed reading '//trim(energyfile)
call Message("(A)")

! this is where we determine the value for the thickness integration limit for the CalcLgh3 routine...
allocate(EkeVs(numEbins),thick(numEbins))

do i=1,numEbins
  EkeVs(i) = Ehistmin + float(i-1)*Ebinsize
end do

! then, for each energy determine the 95% histogram thickness
izzmax = 0
do iE = 1,numEbins
 do ix=-nsx/10,nsx/10
  do iy=-nsy/10,nsy/10
   istat = sum(accum_z(iE,:,ix,iy))
   izz = 1
   do while (sum(accum_z(iE,1:izz,ix,iy)).lt.(0.95*istat)) 
    izz = izz+1
   end do
   if (izz.gt.izzmax) izzmax = izz
  end do
 end do
 thick(iE) = dble(izzmax) * depthstep
end do

izz = nint(maxval(thick)/depthstep)
allocate(lambdaE(1:numEbins,1:izz),stat=istat)
do iE=1,numEbins
 do iz=1,izz
  lambdaE(iE,iz) = float(sum(accum_z(iE,iz,-nsx/10:nsx/10,-nsy/10:nsy/10)))/etotal
 end do
end do

! and get rid of the accum_z array
deallocate(accum_z)
! ---------- end of 'read Monte Carlo output file and extract necessary parameters' section

! ---------- allocate memory for the master pattern
! we need to sample the stereographic projection Northern hemisphere or a portion
! thereoff, depending on the order of the Laue group.  There are 11 Laue groups, 
! which leads to 9 different shapes for the stereographic asymmetric unit for the 
! independent incident beam directions.  
! allocate space for the results (needs to be altered for general symmetry case)
if (Esel.eq.-1) then
  allocate(sr(-npx:npx,-npy:npy,1:numEbins),stat=istat)
else
  allocate(sr(-npx:npx,-npy:npy,1),stat=istat)
end if 

! in the trigonal/hexagonal case, we need intermediate storage arrays
  if (usehex) then
   npyhex = nint(2.0*float(npy)/sqrt(3.0))
   allocate(srhex(-npx:npx,-npyhex:npyhex,1:numEbins),stat=istat)
  end if

! set various arrays to zero
   sr = 0.0
   if (usehex) then
     srhex = 0.0
   end if
! ---------- end allocate memory for the master pattern

! force dynamical matrix routine to read new Bethe parameters from file
! this will all be changed with the new version of the Bethe potentials
   call Set_Bethe_Parameters()
   BetheParameter%weakcutoff = 30.0
   BetheParameter%cutoff = 60.0
   
   numset = cell % ATOM_ntype  
   gzero = 1  ! index of incident beam
   debug = 0  ! no longer used

! ---------- from here on, we need to repeat the entire computation for each energy value
energyloop: do iE=numEbins,1,-1
! is this a single-energy run ?
   if (Esel.ne.-1) then
     if (Esel.ne.iE) CYCLE energyloop
   end if
   
! print a message to indicate where we are in the computation
   io_int(1)=iE
   mess = 'Starting computation for energy bin'
   call Message("(/A$)")
   call WriteValue(' ',io_int,1,"(I4$)")
   io_real(1) = EkeVs(iE)
   call WriteValue('; energy [keV] = ',io_real,1,"(F6.2/)")
   selE = EkeVs(iE)

   nat = 0
!   LUT = dcmplx(0.D0,0.D0)

! set the accelerating voltage
   skip = 3
   call CalcWaveLength(dble(EkeVs(iE)*1000.0),skip)

! ---------- create the master reflection list
! Then we must determine the masterlist of reflections (also a linked list);
! This list basically samples a large reciprocal space volume; it does not 
! distinguish between zero and higher order Laue zones, since that 
! distinction becomes meaningless when we consider the complete 
! fundamental zone.  
 
  call Compute_ReflectionList(dmin,(/0,0,1/),(/1,0,0/),(/0,1,0/),'ALL',.FALSE.,0)
! ---------- end of "create the master reflection list"

! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation

   if (usehex) then
    call Calckvectors( (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,npx,npyhex,numk,isym,ijmax,'RoscaLambert',usehex)
   else 
! Calckvectors(k,ga,ktmax,npx,npy,numk,isym,ijmax,mapmode,usehex)
    call Calckvectors( (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,npx,npy,numk,isym,ijmax,'RoscaLambert',usehex)
   end if
   io_int(1)=numk
   call WriteValue('# independent beam directions to be considered = ', io_int, 1, "(I8)")

! point to the first beam direction
   ktmp => khead
! loop over all beam orientations, selecting them from the linked list
   beamloop:do ik = 1,numk

! compute the dynamical matrix using Bloch waves with Bethe potentials 
     DynFN = ktmp%k
     call Compute_DynMat('BLOCHBETHE', ktmp%k, ktmp%kt, .FALSE.)



! then we need to initialize the Sgh and Lgh arrays
     if (allocated(Sgh)) deallocate(Sgh)
     if (allocated(Lgh)) deallocate(Lgh)

     allocate(Sgh(BetheParameter%nns,BetheParameter%nns),Lgh(BetheParameter%nns,BetheParameter%nns))
     Sgh = czero

! for each special position in the asymmetric unit ...
     totZnsq = 0.D0
     do ip=1,numset
       call CalcOrbit(ip,n,ctmp)  ! get all equivalent points
       nat(ip) = n
! get Zn-squared for this special position and keep track of the total value
       Znsq = float(cell%ATOM_type(ip))**2
       totZnsq = totZnsq + Znsq
! ir is the row index
       do ir=1,BetheParameter%nns
! ic is the column index
         do ic=1,BetheParameter%nns
               kk = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
! We'll assume isotropic Debye-Waller factors for now ...
! That means we need the square of the length of s=  kk^2/4
               kkl = 0.25 * CalcLength(float(kk),'r')**2
               do ikk=1,n
! get the argument of the complex exponential
                 arg = tpi*sum(kk(1:3)*ctmp(ikk,1:3))
! Debye-Waller exponential
                 DBWF = exp(-cell%ATOM_pos(ip,5)*kkl)
!  multiply with the prefactor and add to the structure matrix Sgh
                 Sgh(ir,ic) = Sgh(ir,ic) + cmplx(Znsq * DBWF,0.0) * cmplx(cos(arg),sin(arg))
               end do
           end do
       end do  
     end do

! for now, we're disabling the kinematical part
! solve the dynamical eigenvalue equation for this beam direction  Lgh,thick,kn,nn,gzero,kin,debug
     call CalcLgh3(DynMat,Lgh,dble(thick(iE)),ktmp%kn,BetheParameter%nns,gzero,kin,debug,depthstep,lambdaE(iE,1:izzmax),izzmax)

! dynamical contribution
     s = real(sum(Lgh*Sgh))/float(sum(nat))




! we need to make sure that the Bethe potentials did not accidentally cause a problem (divergence);
! so, if s is very large, then there is likely a problem and we should rerun this incident beam direction
! without the Bethe approximation; there shouldn't be too many of these cases...
     if (s.gt.sthr) then ! we redo this without the Bethe potentials

       write (*,*) 'Warning: caught divergence; redoing computation with full DynMat for point ', ik
 
! temporarily switch off the Bethe potential strong threshold
       Bethe_store = BetheParameter%weakcutoff
       BetheParameter%weakcutoff = BetheParameter%cutoff
              
! compute the dynamical matrix using Bloch waves with Bethe potentials 
       DynFN = ktmp%k
       call Compute_DynMat('BLOCHBETHE', ktmp%k, ktmp%kt, .FALSE.)

! then we need to initialize the Sgh and Lgh arrays
       if (allocated(Sgh)) deallocate(Sgh)
       if (allocated(Lgh)) deallocate(Lgh)

       allocate(Sgh(BetheParameter%nns,BetheParameter%nns),Lgh(BetheParameter%nns,BetheParameter%nns))
       Sgh = czero

! for each special position in the asymmetric unit ...
       totZnsq = 0.D0
       do ip=1,numset
         call CalcOrbit(ip,n,ctmp)  ! get all equivalent points
         nat(ip) = n
! get Zn-squared for this special position and keep track of the total value
         Znsq = float(cell%ATOM_type(ip))**2
         totZnsq = totZnsq + Znsq
! ir is the row index
         do ir=1,BetheParameter%nns
! ic is the column index
           do ic=1,BetheParameter%nns
                 kk = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
! We'll assume isotropic Debye-Waller factors for now ...
! That means we need the square of the length of s=  kk^2/4
                 kkl = 0.25 * CalcLength(float(kk),'r')**2
                 do ikk=1,n
! get the argument of the complex exponential
                   arg = tpi*sum(kk(1:3)*ctmp(ikk,1:3))
! Debye-Waller exponential
                   DBWF = exp(-cell%ATOM_pos(ip,5)*kkl)
!  multiply with the prefactor and add to the structure matrix Sgh
                   Sgh(ir,ic) = Sgh(ir,ic) + cmplx(Znsq * DBWF,0.0) * cmplx(cos(arg),sin(arg))
                 end do
             end do
         end do  
       end do

! for now, we're disabling the kinematical part
! solve the dynamical eigenvalue equation for this beam direction  Lgh,thick,kn,nn,gzero,kin,debug
       call CalcLgh3(DynMat,Lgh,dble(thick(iE)),ktmp%kn,BetheParameter%nns,gzero,kin,debug,depthstep,lambdaE(iE,1:izzmax),izzmax)

! dynamical contribution (we'll assume that this time it is ok
       s = real(sum(Lgh*Sgh))/float(sum(nat))

! and reset the weak cutoff
       BetheParameter%weakcutoff = Bethe_store
     end if


! and store the resulting values
     ipx = ktmp%i
     ipy = ktmp%j
     call Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)
     if (usehex) then
       do ix=1,nequiv
         srhex(iequiv(1,ix),iequiv(2,ix),iE) = s
	end do
     else
        if (Esel.eq.-1) then
         do ix=1,nequiv
           sr(iequiv(1,ix),iequiv(2,ix),iE) = s
	  end do
	else
         do ix=1,nequiv
           sr(iequiv(1,ix),iequiv(2,ix),1) = s
	  end do
        endif
     end if
  
     if (mod(ik,2500).eq.0) write (*,*) 'completed beam direction ',ik
! select next beam direction
     ktmp => ktmp%next
    end do beamloop

! stop the clock and report the total time     
!   call Time_stop(numk)

! write (*,*) 'Some statistics :'
! write (*,*) 'Average number of strong beams : ',float(BetheParameter%totstrong)/float(numk)
! write (*,*) '          (min,max) : ',BetheParameter%minstrong,BetheParameter%maxstrong
! write (*,*) 'Average number of weak beams : ',float(BetheParameter%totweak)/float(numk)
! write (*,*) '          (min,max) : ',BetheParameter%minweak,BetheParameter%maxweak


! Finally, if this was sampled on a hexagonal array, we need to do barycentric interpolation
! to the standard square array for final output and use of the subsequent program.
! [this interpolation scheme must be verified; it is possible that there is an off-by-one error somewhere ...]
!if (usehex) then
!  delta = dsqrt(2.D0)/dble(npx)
!  srt = 2.D0/dsqrt(3.D0)
!! copy the central row without modifications
!  sr(-npx:npx,0) = srhex(-npx:npx,0)
!  srkin(-npx:npx,0) = srkinhex(-npx:npx,0)
!! we'll go through the array with pairs of horizontal rows at a time
!  do j=1,npy-1
!! determine which way the triangle is oriented for this row of the square array
!    jh = floor(j*srt)
!    if (mod(jh,2).eq.0) then ! even numbers mean triangle points down
!      h = delta/srt - (j*delta - float(jh)*delta/srt)
!      lambda = 0.5D0 - h/delta/dsqrt(3.D0)
!      omtl = 1.D0-2.D0*lambda
!      do i=-npx+1,npx-1  ! perform the barycentric interpolation
!! positive row, pay attention to hexagonal coordinate transformation !
!	sr(i,j) = ( srhex(i-1,jh+1) + srhex(i,jh+1) )*lambda + omtl * srhex(i,jh)
!	srkin(i,j) = ( srkinhex(i-1,jh+1) + srkinhex(i,jh+1) )*lambda + omtl * srkinhex(i,jh)
!! negative row
!	sr(i,-j) = ( srhex(i-1,-jh-1) + srhex(i,-jh-1) )*lambda + omtl * srhex(i,-jh)
!	srkin(i,-j) = ( srkinhex(i-1,-jh-1) + srkinhex(i,-jh-1) )*lambda + omtl * srkinhex(i,-jh)
!      end do
!    else
!      h = j*delta - float(jh)*delta/srt
!      lambda = 0.5D0 - h/delta/dsqrt(3.D0)
!      omtl = 1.D0-2.D0*lambda
!      do i=-npx+1,npx-1  ! perform the barycentric interpolation
!! positive row, pay attention to hexagonal coordinate transformation !
!	sr(i,j) = ( srhex(i-1,jh) + srhex(i,jh) )*lambda + omtl * srhex(i,jh+1)
!	srkin(i,j) = ( srkinhex(i-1,jh) + srkinhex(i,jh) )*lambda + omtl * srkinhex(i,jh+1)
!! negative row
!	sr(i,-j) = ( srhex(i-1,-jh) + srhex(i,-jh) )*lambda + omtl * srhex(i,-jh-1)
!	srkin(i,-j) = ( srkinhex(i-1,-jh) + srkinhex(i,-jh) )*lambda + omtl * srkinhex(i,-jh-1)
!      end do
!    end if
!  end do
!end if

! since these computations can take a long time, here we store 
! all the output at the end of each pass through the energyloop.

  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form = 'unformatted')
! write the program identifier
  write (dataunit) trim(progname)
! write the version number
  write (dataunit) scversion
! then the name of the crystal data file
  write (dataunit) xtalname
! then the name of the corresponding Monte Carlo data file
  write (dataunit) energyfile
! energy information and array size    
  if (Esel.eq.-1) then
    write (dataunit) npx,npy,numEbins 
    write (dataunit) EkeVs
  else
    one = 1
    write (dataunit) npx,npy,one 
    write (dataunit) selE
  end if
! is this a regular (square) or hexagonal projection ?
  if (usehex) then 
    write (dataunit) 'hexago'
  else
    write (dataunit) 'square'
  end if
! and finally the results array
  write (dataunit) sr
  close(unit=dataunit,status='keep')

 if ((Esel.eq.-1).and.(iE.ne.1)) then 
  mess = 'Intermediate data stored in file '//trim(outname)
  call Message("(A/)")
 end if

 if ((Esel.eq.-1).and.(iE.eq.1)) then 
  mess = 'Final data stored in file '//trim(outname)
  call Message("(A/)")
 end if

end do energyloop

if (Esel.ne.-1) then
  mess = 'Final data stored in file '//trim(outname)
  call Message("(A/)")
end if

end subroutine ComputeMasterPattern

! ###################################################################
! 
!  subroutine CalcLgh3
! 
!  Author: Marc De Graef
!  
!  Description: integrate the dynamical equations using the Bloch Wave
!  formalism.
! 
!  History
! 
!> @todo this may need to be rewritten to use the original LAPACK routines
!> instead of the eispack routines
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   7/04/01 MDG 2.0 f90
!  12/12/12 MDG 3.0 attempt at actual numerical integration for the matrix elements I_jk
! ###################################################################
subroutine CalcLgh3(DMat,Lgh,thick,kn,nn,gzero,kin,debug,depthstep,lambdaE,izz)

use local
use io
use files
use constants
use error

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)        :: DMat(nn,nn)
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn)
real(kind=dbl),INTENT(IN)           :: thick
real(kind=dbl),INTENT(IN)           :: kn
integer(kind=irg),INTENT(IN)        :: nn
integer(kind=irg),INTENT(IN)        :: gzero
integer(kind=irg),INTENT(IN)        :: debug
real(kind=dbl),INTENT(OUT)          :: kin
real(kind=dbl),INTENT(IN)           :: depthstep
real(kind=sgl),INTENT(IN)	      :: lambdaE(izz)
integer(kind=irg),INTENT(IN)        :: izz

integer         		:: i,j, iz
complex(kind=dbl) 		:: CGinv(nn,nn), Minp(nn,nn), tmp3(nn,nn)

real(kind=dbl)  		:: tpi, dzt
complex(kind=dbl) 		:: Ijk(nn,nn), q, getMIWORK, qold

integer(kind=irg)    		:: INFO, LDA, LDVR, LDVL,  JPIV(nn), MILWORK
complex(kind=dbl)    		:: CGG(nn,nn), W(nn)
complex(kind=dbl),allocatable 	:: MIWORK(:)

integer(kind=irg),parameter   	:: LWMAX = 5000 
complex(kind=dbl)    		:: VL(nn,nn),  WORK(LWMAX)
real(kind=dbl)       		:: RWORK(2*nn)
character            		:: JOBVL, JOBVR
integer(kind=sgl)              :: LWORK


! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
! 
 Minp = DMat

! set some initial LAPACK variables 
 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0
 
! first initialize the parameters for the LAPACK ZGEEV, CGETRF, and CGETRI routines
 JOBVL = 'N'   ! do not compute the left eigenvectors
 JOBVR = 'V'   ! do compute the right eigenvectors
 LWORK = -1 ! so that we can ask the routine for the actually needed value

! call the routine to determine the optimal workspace size
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

! then call the eigenvalue solver
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  if (INFO.ne.0) call FatalError('Error in CalcLgh3: ','ZGEEV return not zero')

!if (debug.eq.1) then 
!  write (*,*) 'inside CalcLgh3'
!  write (*,*) 'thick, depthstep ', thick, depthstep
!  write (*,*) 'kn, nn ', kn, nn
!  write (*,*) 'lambdaE', lambdaE
!  write (*,*) 'gzero = ', gzero
!  write (*,*) 'max(DMat) ', maxval(cabs(DMat))
!end if

!! then get the eigenvalues and eigenvectors
! 
! IPIV = 0
! W = cmplx(0.D0,0.D0)
! CGG = cmplx(0.D0,0.D0)
! CGinv = cmplx(0.D0,0.D0)
!
! LDA = nn
! LDVL = nn
! LDVR = nn
! INFO = 0
!! the eispack routines use the real and imaginary parts of arrays as separate arrays instead
!! of as complex variable ...  
!!if (debug.eq.1) write (*,*) 'splitting real and imaginary'
!
!ar = dble(Minp)
!ai = aimag(Minp)
!
!! first balance the matrix
!!if (debug.eq.1) write (*,*) 'cbal'
!
!call cbal(nn, ar, ai, low, igh, scle)
!! transform the upper Hessenberg form
!allocate(ortr(igh), orti(igh))
!!if (debug.eq.1) write (*,*) 'corth'
!call corth(nn, low, igh, ar, ai, ortr, orti)
!! get eigenvalues and eigenvectors
!!if (debug.eq.1) write (*,*) 'comqr2'
!call comqr2(nn, low, igh, ortr, orti, ar, ai, wr, wi, vr, vi, ierr )
!if ( ierr.ne.0) then 
!  write (*,*) 'Error in comqr2 eispack routine'
!  stop
!end if
!deallocate(ortr, orti)
!! undo the cbal transformation
!!if (debug.eq.1) write (*,*) 'cbabk2'
!call cbabk2(nn, low, igh, scle, nn, vr, vi)
!
!! return to complex variables
!W = cmplx(wr, wi)
!CGG = cmplx(vr, vi)
! 

 CGinv = CGG

!if (debug.eq.1) write (*,*) 'inverting matrix'
 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 MILWORK = -1
 call zgetri(nn,CGinv,LDA,JPIV,getMIWORK,MILWORK,INFO)
 MILWORK =  INT(real(getMIWORK))
 if (.not.allocated(MIWORK)) allocate(MIWORK(MILWORK))
 MIWORK = dcmplx(0.D0,0.D0)
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
 deallocate(MIWORK)

! in all the time that we've used these routines, we haven't
! had a single problem with the matrix inversion, so we don't
! really need to do this test:
!
! if ((cabs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) write (*,*) 'Error in matrix inversion; continuing'


! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)
 
! if (debug.eq.1) then
!   write (*,*) 'W = ',W
! end if
 
! recall that alpha(1:nn) = CGinv(1:nn,gzero)


! first the Ijk matrix (this is Winkelmann's B^{ij}(t) matrix)
! combined with numerical integration over [0, z0] interval,
! taking into account depth profiles from Monte Carlo simulations ...
! the depth profile lambdaE must be added to the absorption 
! components of the Bloch wave eigenvalues.

tpi = 2.D0*cPi*depthstep
dzt = depthstep/thick
 do i=1,nn
  do j=1,nn
     q =  cmplx(0.D0,0.D0)
     qold = tpi * dcmplx(aimag(W(i))+aimag(W(j)),real(W(i))-real(W(j)))
     do iz = 1,izz
       q = q + dble(lambdaE(iz)) * cexp( - qold * dble(iz) )
     end do
     Ijk(i,j) = conjg(CGinv(i,gzero)) * q * CGinv(j,gzero)
  end do
 end do

Ijk = Ijk * dzt

!if (debug.eq.1) then 
!  write (*,*) 'tpi, dzt ', tpi, dzt
!  write (*,*) 'maxval(Ijk) ', maxval(cabs(Ijk))
!end if


!if (debug.eq.1) write (*,*) ' -> done'

!if (debug.eq.1) write (*,*) 'matmul operations'

! then the summations for Lgh and kin
tmp3 = matmul(CGG,transpose(Ijk)) 
Lgh = matmul(tmp3,transpose(conjg(CGG)))

! there might be a problem with the Absoft implementation of the 
! matmul routine...  So let's do this multiplication explicitly...
!do i=1,nn
!  do j=1,nn
!     tmp3(i,j) = sum( CGG(i,1:nn) * Ijk(j,1:nn) )
!  end do
!end do
!
!! we no longer need CGinv, so we'll use the array to store the conjugate of CGG
!CGinv = conjg(CGG)
!do i=1,nn
!  do j=1,nn
!    Lgh(i,j) = sum( tmp3(i,1:nn) * CGinv(j,1:nn) )
!  end do
!end do
!

! we'll approximate this by the sum of the diagonal entries 
! tmp3 = matmul(transpose(CGG),conjg(CGG)) 
! kin = 1.D0-real(sum(Ijk * tmp3))
!if (debug.eq.1) write (*,*) ' -> done'

! there may be an issue with this part of the routine ... 
!do i=1,nn
!  do j=1,nn
!    tmp3(i,j) = sum( CGG(1:nn,i) * CGinv(1:nn,j) )
!  end do
!end do
!kin = 1.D0-real(sum(Ijk * tmp3))
!


!deallocate(CGinv,Minp,diag,tmp3)

end subroutine CalcLgh3

     
