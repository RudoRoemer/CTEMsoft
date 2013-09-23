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
! CTEMsoft2013:CTEMEBSDzE.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMEBSDzE 
!
!> @author Marc De Graef, Carnegie Melon University
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
!--------------------------------------------------------------------------
program CTEMEBSDzE

use local
use io

IMPLICIT NONE

character(fnlen)			:: nmlfile

integer(kind=irg)			:: numarg		!< number of command line arguments
integer(kind=irg)			:: iargc		!< external function for command line
character(fnlen)    			:: arg			!< to be read from the command line

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
		mess = '        CTEMEBSDzE [nmlfile]'; call Message("(A)")
		mess = ' where nmlfile is an optional file name for the namelist file;'; call Message("(A)")
		mess = ' if absent, the default name ''CTEMEBSDzE.nml'' will be used.'; call Message("(A/)")
		stop
	end if
else
	nmlfile = 'CTEMEBSDzE.nml'    		! assign the default namelist file name
end if


! generate a set of zone axis CBED patterns
 call ComputeMasterPattern(nmlfile)

end program CTEMEBSDzE

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeMasterPattern
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute an EBSD master pattern as a function of energy
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 08/01/13  MDG 3.0 complete rewrite, eliminated old Lambert projection
!--------------------------------------------------------------------------
subroutine ComputeMasterPattern(nmlfile)

use symmetryvars
use symmetry
use crystalvars
use crystal
use constants
use io
use local
use files
use diffraction
use multibeams
use dynamical
use timing
use Lambert

IMPLICIT NONE

character(100)    	:: fname   ! = output filename
real(kind=dbl)   	:: ctmp(192,3),arg
integer(kind=irg)      :: isym,i,j,ik,npx,npy,ipx,ipy,debug,iE,izz, & ! counters
                    	numk, intfreq, & ! number of independent incident beam directions
                    	nn,nns,nnw, & ! total number of reflections in computation, and thicknesses
                    	ir,nat(100),kk(3),numrefs, iw, npyhex,jh, &
                    	numset,n,ig,gg(3),ix,iy,iz, imh, imk, iml, voltage,  &
                    	istat,gzero,ic,ip,ikk,iweak,istrong, minweak,maxweak,totweak, minstrong, maxstrong, totstrong     ! counters
real(kind=dbl)         :: pre, sgp,delta,srt, lambda, omtl, h,  & 
                    	tpi,thick,ll(3),lpg(3),gplen,cutoff, Znsq, kkl, weakcutoff, sgcutoff,&
                    	DBWF,  totZnsq, kin, ave  !
real(kind=sgl) 		:: dmin, dhkl 
real(kind=sgl),allocatable 	:: sr(:,:), srkin(:,:), weaksg(:), strongsg(:), srhex(:,:), srkinhex(:,:) ! results
integer(kind=irg),allocatable 	:: weaklist(:), stronglist(:), weakhkl(:,:), stronghkl(:,:)
complex(kind=dbl)  		:: czero, weaksum, ughp, uhph
complex(kind=dbl),allocatable 	:: Lgh(:,:), Sgh(:,:), LUT(:,:,:), DMat(:,:)
logical 		:: usehex, switchmirror
character(fnlen) 	:: xtalname
character(100) 		:: intfname
character(12) 		:: mapmode    ! 'PlainLambert' or 'RoscaLambert'
integer(kind=irg),parameter	:: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees

! Monte Carlo derived quantities
integer(kind=irg)   	:: numEbins, numzbins, nsx, nsy, nE, Emin, Emax  ! variables used in MC energy file
real(kind=dbl)     	:: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, etotal ! enery variables from MC program
integer(kind=irg),allocatable :: accum_e(:,:,:), accum_z(:,:), lambdaE(:)
character(80) 		:: energyfile


namelist /EBSDzEvars/ xtalname,voltage,dmin, weakcutoff, cutoff, sgcutoff,npx,npy, &
                      thick,intfreq,intfname,fname,debug,mapmode,energyfile, etotal

! read namelist file with all input parameters
! here are the default values in case some parameters are absent from the file
voltage		= 30000  	! default accelerating voltage
dmin 		= 0.015   	! minimum d-spacing in order to be admitted to the list (must become user entry)
weakcutoff 	= 30.0     	! dimensionless cutoff parameter (smaller=strong, larger=weak)
cutoff 		= 80.0		! overall cutoff parameter (dimensionless)
sgcutoff 	= 0.1		! cutoff for excitation error to avoid small denominators in Bethe perturbation computation
npx		= 500		! Nx pixels (total = 2Nx+1)
npy		= 500		! Ny pixels (total = 2Ny+1)
!thick		= 100.0		! integration thickness [nm]  This is replaced by the z_0(E_e) dependence.
intfreq		= 0		! frequency of intermediate data storage (every intfreq wave vector directions)
intfname	= 'none'	! full pathname to intermediate storage file
fname		= 'EBSDzEout.data'	! default filename for final output
mapmode		= 'RoscaLambert'	! use square grid -> sphere (PlainLambert) or square grid -> circle -> sphere (RoscaLambert) mapping
debug		= 0
energyfile	= 'z0E.data' 	! default filename for z_0(E_e) data from Monte Carlo simulations
etotal		= 1.0D9 	! total number of electrons from MC simulation.

! then we read the rundata namelist, which may override some of these defaults  
OPEN(UNIT=dataunit,FILE=nmlfile,DELIM='APOSTROPHE')
READ(UNIT=dataunit,NML=EBSDzEvars)
CLOSE(UNIT=dataunit)

if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMEBSDzE:',' structure file name is undefined in '//nmlfile)
end if

! print some information
progname = 'CTEMEBSDzE.f90'
progdesc = 'EBSD Energy-dependent Master Pattern Simulation'
call CTEMsoft

! first get the crystal data and microscope voltage (that's the initial voltage)
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

! generate all atom positions
 call CalcPositions('v')
 
!  hexset = .TRUE.    ! if hexagonal structure this switch selects between three and four index notation (4 if true)

! determine the crystal point group
  j=0
  do i=1,32
   if (SGPG(i).le.cell % SYM_SGnum) j=i
  end do
  isym = j   ! point group number corresponding to the crystal structure's space group
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


! first, we need to load the data from the MC program.  This is an array on integers, which
! has an energy histogram for each scintillator pixel.  These are integer values to minimize
! storage, but they should all be divided by the total number of counts, which can be done 
! at the end of the computation.  We will then need to multiply by the beam current and by
! the dwell time to get units of electron counts.
!
write (*,*) 'opening ',trim(energyfile)
open(dataunit,file=trim(energyfile),status='unknown',form='unformatted')
read(dataunit) numEbins, numzbins, nsx, nsy
write (*,*) numEbins, numzbins, nsx, nsy
read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
write (*,*) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
allocate(accum_e(numEbins,nsx,nsy),accum_z(numEbins,numzbins),stat=istat)
read(dataunit) accum_e
deallocate(accum_e)
read(dataunit) accum_z    ! we only need this array for the depth integrations
close(dataunit,status='keep')

! this is where we determine which line we need from the array accum_z, and
! also we set the value for the thickness integration limit for the CalcLgh3 routine...

! first get the correct energy index
iE = nint((voltage/1000. - Ehistmin)/Ebinsize) +1
! then, for this energy determine the 95% histogram thickness
istat = sum(accum_z(iE,:))
izz = 1
do while (sum(accum_z(iE,1:izz)).lt.(0.95*istat)) 
  izz = izz+1
end do
thick = dble(izz) * depthstep
allocate(lambdaE(1:izz),stat=istat)
lambdaE(1:izz) = accum_z(iE,1:izz)/etotal
deallocate(accum_z)
write (*,*) 'Integration depth = ',thick
write (*,*) 'energy index        = ',iE
write (*,*) 'depth index	     = ',izz
write (*,*) lambdaE


! we need to sample the stereographic projection Northern hemisphere or a portion
! thereoff, depending on the order of the Laue group.  There are 11 Laue groups, 
! which leads to 9 different shapes for the stereographic asymmetric unit for the 
! independent incident beam directions.  
! allocate space for the results (needs to be altered for general symmetry case)
   allocate(sr(-npx:npx,-npy:npy),srkin(-npx:npx,-npy:npy),stat=istat)
   sr = 0.0
   srkin = 0.0

! in the trigonal/hexagonal case, we need intermediate storage arrays
  if (usehex) then
   npyhex = nint(2.0*float(npy)/sqrt(3.0))
   allocate(srhex(-npx:npx,-npyhex:npyhex),srkinhex(-npx:npx,-npyhex:npyhex),stat=istat)
   srhex = 0.0
   srkinhex = 0.0
  end if

  nat = 0

! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation
if (usehex) then
  call Calckvectors(npx,npyhex,numk,isym,usehex,mapmode)
else 
  call Calckvectors(npx,npy,numk,isym,usehex,mapmode)
end if
mess = '# independent beam directions to be considered = '; oi_int(1)=numk; call WriteInt(1,"(I8)")

! Then we must determine the masterlist of reflections (also a linked list);
! This list basically samples a large reciprocal space volume; it does not 
! distinguish between zero and higher order Laue zones, since that 
! distinction becomes meaningless when we consider the complete 
! fundamental zone.  

! The master list is easily created by brute force
 imh = 1
 do 
   imh = imh + 1
   dhkl = 1.0/CalcLength(  (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
 end do
 imk = 1
 do 
   imk = imk + 1
   dhkl = 1.0/CalcLength( (/0.0_sgl,float(imk),0.0_sgl/), 'r')
  if (dhkl.lt.dmin) EXIT
 end do
 iml = 1
 do 
   iml = iml + 1
   dhkl = 1.0/CalcLength( (/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.dmin) EXIT
 end do
 write (*,*) 'Range of reflections along a*, b* and c* = ',imh,imk,iml

  numrefs = 1
  gzero = 1
  kk = (/ 0,0,0 /)
  call AddReflection( kk )   ! this guarantees that 000 is always the first reflection
! the LUT array stores all the Fourier coefficients, so that we only need to compute them once...
  allocate(LUT(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
  LUT = dcmplx(0.D0,0.D0)
  call CalcUcg(kk)   
  DynUpz = rlp%Vpmod
  write (*,*) 'Normal absorption length = ',rlp%xgp
! and add this reflection to the look-up table
  LUT(kk(1),kk(2),kk(3)) = rlp%Ucg
! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections !!!
    do ix=-2*imh,2*imh
      do iy=-2*imk,2*imk
       do iz=-2*iml,2*iml
        gg = (/ ix, iy, iz /)
        if ((IsGAllowed(gg)).AND.(sum(abs(gg)).ne.0)) then
           if ((abs(ix).le.imh).and.(abs(iy).le.imk).and.(abs(iz).le.iml)) then 
             numrefs = numrefs + 1
             call AddReflection(gg)
           end if
           call CalcUcg(gg)
           LUT(gg(1),gg(2),gg(3)) = rlp%Ucg
        end if
       end do
      end do
    end do
  mess = 'Length of the master list of reflections : '; oi_int(1) = numrefs; call WriteInt(1,"(I8,/)")

   
! next, we start the major loop over all incident beam directions...
! for each value, we must compute how many reflections must be taken 
! into account (we'll keep track of min and max number), and then
! we need to initialize the dynamical matrix and solve the equations;
! this is to be repeated for all indicent beam directions.  To reduce the
! size of the dynamical matrix, we use Bethe potentials.
  
! define the cutoff for the product of excitation error and extinction
! distance; only reflections for which this product is smaller than
! cutoff can contribute to the dynamical matrix.  If the product is 
! smaller than weakcutoff, then the beam is strong, if it is between
! weakcutoff and cutoff, then it is a weak beam that needs to be 
! taken into account via the Bethe potential mechanism.  These 
! parameters are set by trial and error, and they should probably be
! made user-adjustable.

  allocate(weaklist(numrefs),stronglist(numrefs))
  weaklist = 0
  stronglist = 0
  minweak = numrefs
  maxweak = 0
  totweak = 0
  minstrong = numrefs
  maxstrong = 0
  totstrong = 0

!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
  mess = 'Starting main computational loop'; call Message("(A)")

! point to the first beam direction and init complex zero and a complex prefactor
  tmp => head
  czero = cmplx(0.0,0.0,dbl)
  pre = cmplx(0.0,cPi,dbl)

  numset = cell % ATOM_ntype  ! number of special positions in the unit cell
  tpi = 2.D0*cPi
  
! we should be able to do this loop as a multithreaded loop ... 

!open(unit=25,file='pg02.txt',status='unknown',form='formatted')

!  work through the beam direction list
  beamloop: do ik=1,numk
   if (mod(ik,1000).eq.0) write (*,*) ik,' out of',numk,' beam orientations completed'
! time the computation
   if (ik.eq.1) call Time_start
if (debug.eq.1) write (*,*) 'starting loop index ',ik

!write (*,*) 'at top of beamloop'

! pick the wave vector from the pointer list
    ll = tmp%k 

if (debug.eq.1) write (*,*) 'starting excitation errors '

! first, for this beam direction, determine the excitation errors of 
! all the reflections in the master list, and count the ones that are
! needed for the dynamical matrix (weak as well as strong)
! Also, we put the foil normal parallel to the incident beam direction for each incident beam direction.
! That means that the typical non-uniform intensity distribution of an individual EBSD
! pattern is not incorporated in this approach; it needs to be computed separately
! by means of a Monte Carlo simulation of the BSE yield vs. scintillator position and energy.
    rltmpa => reflist%next
    nn = 1
    stronglist = 0
    weaklist = 0
    stronglist(1) = 1   ! make sure that the origin is always a strong reflection...
    weaklist(1) = 0
    rltmpa%sg = 0.D0    
    rltmpa => rltmpa%next
    reflectionloop: do ig=2,numrefs
      gg = rltmpa%hkl        ! this is the reciprocal lattice parameter 
      lpg = ll+gg                ! k0 + g (vectors)
      gplen = CalcLength(lpg,'r')  ! |k0+g|
! we're taking the foil normal to be parallel to the incident beam direction at each point of
! the standard stereographic triangle, so cos(alpha) = 1 in eqn. 5.11 of CTEM
      rltmpa%sg = (1.0/mLambda**2 - gplen**2)*0.5/gplen
! use the reflection num entry to indicate whether or not this
! reflection should be used for the dynamical matrix
! there are two criteria: one is that |sg|*xi must be smaller than some cutoff.
! the other is that the excitation error itself must be smaller than a cutoff in
! order for the reflectio to be a strong reflection.  So if a reflection is considered
! weak according to the first criterion, it could still be made strong is the 
! excitation error is really small (to avoid a blowup of the perturbation algorithm).
        sgp = abs(rltmpa%sg) * rltmpa%xg
        if (sgp.le.cutoff) then
          nn = nn+1
! is this a weak or a strong reflection (in terms of Bethe potentials)? 
          if ((sgp.le.weakcutoff).or.(abs(rltmpa%sg).le.sgcutoff)) then 
            stronglist(ig) = 1
          else
            weaklist(ig) = 1
          end if
        end if
      rltmpa => rltmpa%next
    end do reflectionloop
if (debug.eq.1) write (*,*) 'ending excitation errors ',nn

! if we don't have any beams in this list (unlikely, but possible if the cutoff and 
! weakcutoff parameters have unreasonable values) then we skip the rest of this
! incident beam computation (and we report it to the user) 
 if (nn.eq.0) then
   mess = 'no beams found for the following parameters:'; call Message("(A)")
   write (*,*) ik,tmp%i,tmp%j,ll,'  -> number of beams = ',nn
   mess =  ' -> check cutoff and weakcutoff parameters for reasonableness'; call Message("(A)")
   mess =  ' [skipping this beam direction] '; call Message("(A)")
   tmp => tmp%next
   cycle beamloop
end if

! next, we define nns to be the number of strong beams, and nnw the number 
! of weak beams.
 nns = sum(stronglist)
 nnw = sum(weaklist)
if (debug.eq.1) write (*,*) 'strong/weak ',nns,nnw

! We may want to keep track of the total and average numbers of strong and weak beams  
 totweak = totweak + nnw
 totstrong = totstrong + nns
 if (nnw.lt.minweak) minweak=nnw
 if (nnw.gt.maxweak) maxweak=nnw
 if (nns.lt.minstrong) minstrong=nns
 if (nns.gt.maxstrong) maxstrong=nns

if (debug.eq.1) write (*,*) 'allocating weak/strong arrays '

! allocate arrays for weak and strong beam information
if (allocated(weakhkl)) deallocate(weakhkl)
if (allocated(weaksg)) deallocate(weaksg)
if (allocated(stronghkl)) deallocate(stronghkl)
if (allocated(strongsg)) deallocate(strongsg)
allocate(weakhkl(3,nnw),weaksg(nnw))
allocate(stronghkl(3,nns),strongsg(nns))

if (debug.eq.1) write (*,*) 'converting linked lists'

! here's where we extract the relevant information from the linked list (much faster
! than traversing the list each time...)
rltmpa => reflist%next    ! reset the a list
iweak = 0
istrong = 0
do ir=1,numrefs
     if (weaklist(ir).eq.1) then
        iweak = iweak+1
        weakhkl(1:3,iweak) = rltmpa%hkl(1:3)
        weaksg(iweak) = rltmpa%sg
     end if
     if (stronglist(ir).eq.1) then
        istrong = istrong+1
        stronghkl(1:3,istrong) = rltmpa%hkl(1:3)
        strongsg(istrong) = rltmpa%sg
     end if
   rltmpa => rltmpa%next
end do

! deallocate any previous Bloch wave arrays and reallocate them
! with the current number of beams
    if (allocated(DMat)) deallocate(DMat)

! initialize the dynamical matrix
    allocate(DMat(nns,nns),stat=istat) 
    DMat = czero

if (debug.eq.1) write (*,*) 'init DMat'
    
! ir is the row index
    do ir=1,nns
! ic is the column index
         do ic=1,nns
! compute the Bethe Fourier coefficient of the electrostatic lattice potential 
              if (ic.ne.ir) then  ! not a diagonal entry
                kk = stronghkl(1:3,ir) - stronghkl(1:3,ic)
                DMat(ir,ic) = LUT(kk(1),kk(2),kk(3)) 
! and subtract from this the total contribution of the weak beams
                weaksum = czero
                do iw=1,nnw
                      kk = stronghkl(1:3,ir) - weakhkl(1:3,iw)
                      ughp = LUT(kk(1),kk(2),kk(3)) 
                      kk = weakhkl(1:3,iw) - stronghkl(1:3,ic)
                      uhph = LUT(kk(1),kk(2),kk(3)) 
                      weaksum = weaksum +  ughp * uhph *cmplx(1.D0/weaksg(iw),0.0,dbl)
                 end do
! and correct the dynamical matrix element to become a Bethe potential coefficient
                 DMat(ir,ic) = DMat(ir,ic) - cmplx(0.5D0*mLambda,0.0D0,dbl)*weaksum
               else  ! it is a diagonal entry, so we need the excitation error and the absorption length
                 DMat(ir,ir) = cmplx(2.D0*strongsg(ir)/mLambda,DynUpz,dbl)
               end if
         end do
      end do
! that should do it for the initialization of the dynamical matrix, except for the first element
DMat(1,1) =  cmplx(0.0,DynUpz,dbl)

if (debug.eq.1) write (*,*) 'DMat ready'

! then we need to initialize the Sgh array;
  if (allocated(Sgh)) deallocate(Sgh)
  if (allocated(Lgh)) deallocate(Lgh)

  allocate(Sgh(nns,nns),Lgh(nns,nns))
  Sgh = czero

if (debug.eq.1) write (*,*) 'init position part'

! for each special position in the asymmetric unit ...
  totZnsq = 0.D0
  do ip=1,numset
    call CalcOrbit(ip,n,ctmp)  ! get all equivalent points
    nat(ip) = n
! get Zn-squared for this special position and keep track of the total value
    Znsq = float(cell%ATOM_type(ip))**2
    totZnsq = totZnsq + Znsq
! ir is the row index
    do ir=1,nns
! ic is the column index
         do ic=1,nns
               kk = stronghkl(1:3,ir) - stronghkl(1:3,ic)
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

if (ik.eq.-1) then
  open(unit=20,file='array.data',status='unknown',form='unformatted')
  write (20) nns
  write (20) DMat
  close(unit=20,status='keep')
end if

if (debug.eq.1) write (*,*) 'Solving eigenvalue problem'  

! set up the z integration arrays  This might actually be incorrect, but we need to try it !!!


! solve the dynamical eigenvalue equation for this beam direction  Lgh,thick,kn,nn,gzero,kin,debug
  call CalcLgh3(DMat,Lgh,thick,tmp%kn,nns,gzero,kin,debug,depthstep,lambdaE,izz)
! call CalcLgh2(DMat,Lgh,thick,tmp%kn,nns,gzero,kin,debug)
!  call CalcLgh(DMat,nns,thick,tmp%kn,gzero,Lgh,kin,debug)
if (debug.eq.1) write (*,*) 'Return fromCalcLgh ',nnw,nns,sum(Lgh),kin,sum(Sgh),tmp%i,tmp%j

! and store the resulting values
   ipx = tmp%i
   ipy = tmp%j
 if (usehex) then
! dynamical contribution
   srhex(ipx,ipy) = real(sum(Lgh*Sgh))/float(sum(nat))
! kinematical contribution
   srkinhex(ipx,ipy) = kin * totZnsq/float(sum(nat))
   if (srhex(ipx,ipy).gt.500.0) then
     write (*,*) ' abnormal value detected; weak excitation errors are: '
     write (*,*) weaksg
   end if
   if (debug.eq.1) write (*,*) 'results ', srhex(ipx,ipy), srkinhex(ipx,ipy), intfreq, ik
else
! dynamical contribution
   sr(ipx,ipy) = real(sum(Lgh*Sgh))/float(sum(nat))
! kinematical contribution
   srkin(ipx,ipy) = kin * totZnsq/float(sum(nat))
!   if (sr(ipx,ipy).gt.5000.0) then
!     write (*,*) ' abnormal value detected; weak excitation errors are: '
!     write (*,*) weaksg
!     write (*,*) 'This might indicate that sgcutoff is too small, leading to '
!     write (*,*) 'problems with the Bethe potential perturbation model ... '
!   end if
   if (debug.eq.1) write (*,*) 'results ', sr(ipx,ipy), srkin(ipx,ipy), intfreq, ik
end if


! save intermediate data if requested
 if (intfreq.ne.0) then
  if (mod(ik,intfreq).eq.0) then
     mess= 'intermediate data saved in '//intfname; call Message("(A)")
     open(unit=dataunit,file=intfname,status='unknown',form='unformatted')
   if (usehex) then
     write (dataunit) 2*npx+1,2*npyhex+1
     write (dataunit) srhex
     write (dataunit) srkinhex  
   else
     write (dataunit) 2*npx+1,2*npy+1
     write (dataunit) sr
     write (dataunit) srkin  
   end if
     close(unit=dataunit, status='keep')
     write (*,*) 'Average number of strong beams : ',float(totstrong)/float(ik)
     write (*,*) '          (min,max) : ',minstrong,maxstrong
     write (*,*) 'Average number of weak beams : ',float(totweak)/float(ik)
     write (*,*) '          (min,max) : ',minweak,maxweak
  end if  
end if
  
! select next beam direction
!write (*,*) 'switching to next beam direction'
   tmp => tmp%next
if (debug.eq.1) write (*,*) 'reached end of loop'
end do beamloop

!close(unit=25,status='keep')

tmp=>head
if (usehex) then 
   call ApplySymmetry(numk,isym,srhex,srkinhex,npx,npyhex)
else
   call ApplySymmetry(numk,isym,sr,srkin,npx,npy)
end if

! stop the clock and report the total time     
  call Time_stop(numk)

 write (*,*) 'Some statistics :'
 write (*,*) 'Average number of strong beams : ',float(totstrong)/float(numk)
 write (*,*) '          (min,max) : ',minstrong,maxstrong
 write (*,*) 'Average number of weak beams : ',float(totweak)/float(numk)
 write (*,*) '          (min,max) : ',minweak,maxweak

! clean up the kinematical data (remove spikes); right now (8/29/12), I'm not sure why 
! these spikes occur.  It must be because of something that's ill-defined 
! in the perturbation approach of the Bethe potentials, but it appears to
! only affect the kinematical part, not the dynamical part and that's rather odd.
! For now, we look for these spikes and remove them by placing the average
! value of the nearest neighbours in that location ... 
if (usehex) then ! the hexagonal array has different dimensions 
  do i=-npx+1,npx-1
   do j=-npyhex+1,npyhex-1
     ave = 0.25*(srkinhex(i-1,j) + srkinhex(i+1,j) + srkinhex(i,j-1) + srkinhex(i,j+1) ) 
     if (abs(srkinhex(i,j)).gt.ave*3.0) srkinhex(i,j)=ave
     if (abs(srkinhex(i,j)).lt.ave*0.333) srkinhex(i,j)=ave
   end do
 end do
else
  do i=-npx+1,npx-1
   do j=-npy+1,npy-1
     ave = 0.25*(srkin(i-1,j) + srkin(i+1,j) + srkin(i,j-1) + srkin(i,j+1) ) 
     if (abs(srkin(i,j)).gt.ave*3.0) srkin(i,j)=ave
     if (abs(srkin(i,j)).lt.ave*0.333) srkin(i,j)=ave
   end do
 end do
end if
!
mess = 'removed spikes from kinematical array';  call Message("(A)")

! Finally, if this was sampled on a hexagonal array, we need to do barycentric interpolation
! to the standard square array for final output and use of the subsequent program.
! [this interpolation scheme must be verified; it is possible that there is an off-by-one error somewhere ...]
if (usehex) then
  delta = dsqrt(2.D0)/dble(npx)
  srt = 2.D0/dsqrt(3.D0)
! copy the central row without modifications
  sr(-npx:npx,0) = srhex(-npx:npx,0)
  srkin(-npx:npx,0) = srkinhex(-npx:npx,0)
! we'll go through the array with pairs of horizontal rows at a time
  do j=1,npy-1
! determine which way the triangle is oriented for this row of the square array
    jh = floor(j*srt)
    if (mod(jh,2).eq.0) then ! even numbers mean triangle points down
      h = delta/srt - (j*delta - float(jh)*delta/srt)
      lambda = 0.5D0 - h/delta/dsqrt(3.D0)
      omtl = 1.D0-2.D0*lambda
      do i=-npx+1,npx-1  ! perform the barycentric interpolation
! positive row, pay attention to hexagonal coordinate transformation !
	sr(i,j) = ( srhex(i-1,jh+1) + srhex(i,jh+1) )*lambda + omtl * srhex(i,jh)
	srkin(i,j) = ( srkinhex(i-1,jh+1) + srkinhex(i,jh+1) )*lambda + omtl * srkinhex(i,jh)
! negative row
	sr(i,-j) = ( srhex(i-1,-jh-1) + srhex(i,-jh-1) )*lambda + omtl * srhex(i,-jh)
	srkin(i,-j) = ( srkinhex(i-1,-jh-1) + srkinhex(i,-jh-1) )*lambda + omtl * srkinhex(i,-jh)
      end do
    else
      h = j*delta - float(jh)*delta/srt
      lambda = 0.5D0 - h/delta/dsqrt(3.D0)
      omtl = 1.D0-2.D0*lambda
      do i=-npx+1,npx-1  ! perform the barycentric interpolation
! positive row, pay attention to hexagonal coordinate transformation !
	sr(i,j) = ( srhex(i-1,jh) + srhex(i,jh) )*lambda + omtl * srhex(i,jh+1)
	srkin(i,j) = ( srkinhex(i-1,jh) + srkinhex(i,jh) )*lambda + omtl * srkinhex(i,jh+1)
! negative row
	sr(i,-j) = ( srhex(i-1,-jh) + srhex(i,-jh) )*lambda + omtl * srhex(i,-jh-1)
	srkin(i,-j) = ( srkinhex(i-1,-jh) + srkinhex(i,-jh) )*lambda + omtl * srkinhex(i,-jh-1)
      end do
    end if
  end do
end if
 
! save the results  
! this file format was revised to account for the different sampling/mapping 
! modes; the file now contains 4 integers, the first two provide the dimensions
! of the array, the the third provides the Laue group, and the last one the 
! particular Lambert map mode (square, hexagonal or triangular).
! in the next version of this program, we'll also have an extra dimension for
! the simulation energy, so that we can interpolate both by location and 
! energy, to get a more accurate distribution on the scintillator.
  open(unit=dataunit,file=trim(fname),status='unknown',action='write',form = 'unformatted')
  write (dataunit)  2*npx+1,2*npy+1 ! , isym, placeholder
  write (dataunit) sr
  write (dataunit) srkin 
  close(unit=dataunit,status='keep')

end subroutine ComputeMasterPattern


! ###################################################################
! 
!  subroutine Calckvectors
!
!  Author: Marc De Graef
!  
!  Description: computes the independent incident beam directions
!  for the multi-beam case and returns them as a linked list in the
!  global variables head and tail; to do this we use the Lambert 
!  equal area projection of a square grid on a plane tangent to the
!  projection sphere.  The square grid allows for the consideration
!  of 7 of the 11 Laue symmetries; a hexagonal grid covers the 
!  other 4.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
!   8/20/12 MDG 1.1 k-vectors determined from Lambert equal are projection 
!  11/13/12 MDG 1.2 added RoscaLambert and PlainLambert distinction for sphere mapping
! ###################################################################
subroutine Calckvectors(npx,npy,numk,isym,usehex,mapmode)

use local
use io
use error
use constants
use diffraction
use crystal
use crystalvars
use dynamical

IMPLICIT NONE

integer(kind=irg)       	:: npx,npy,numk,istat,i,j,isym, istart,iend,jstart,jend
real(kind=dbl)                 :: delta,kstar(3)
character(*)			:: mapmode
logical				:: usehex

intent(IN)           :: npx,npy,isym,usehex,mapmode
intent(OUT)        :: numk


 if (.not.associated(head)) then     ! allocate the head and tail of the linked list
   allocate(head,stat=istat)         ! allocate new value
   if (istat.ne.0) call FatalError('Calckvectors: unable to allocate head pointer',' ')
   tail => head                      ! tail points to new value
   nullify(tail%next)                ! nullify next in new value
   numk = 1                          ! keep track of number of k-vectors so far
   tail%i = 0                             ! i-index of beam
   tail%j = 0                             ! j-index of beam
   kstar = (/ 0.0, 0.0, 1.0 /)
   call NormVecd(kstar,'c')                ! normalize incident direction
   kstar = kstar/mLambda                 ! divide by wavelength
! and transform to reciprocal crystal space using the structure matrix
   tail%k = matmul(transpose(cell%dsm),kstar)
   tail%kn = 1.0/mLambda
 else
   call FatalError('Calckvectors: pointer head already allocated',' ')
 end if
!
if (mapmode.eq.'PlainLambert') then 
   delta = dsqrt(2.D0)/dble(npx)
else 
   if (usehex) then
!      delta =  0.6733868435D0/dble(npx)
      delta =  2.D0*dsqrt(cPi)/3.0D0**0.75D0/dble(npx)
   else
      delta = dsqrt(cPi*0.5D0)/dble(npx)
   end if
end if

select case (isym)
   case (1)  ! triclinic symmetry
   	istart = -npx
	iend = npx
   	jstart = -npy
	jend = npy
	  do j=jstart,jend
	    do i=istart,iend   ! 
		call AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)
	    end do
	  end do
  case (2)   !  monoclinic symmetry
  	istart = -npx
	iend = npx
   	jstart = 0
	jend = npy
	  do j=jstart,jend
	   do i=istart,iend   ! 
		call AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)
	   end do
	  end do
  case (3,4,10)  ! orthorhombic mmm, tetragonal 4/m, cubic m-3
  	istart = 0
	iend = npx
   	jstart = 0
	jend = npy
	  do j=jstart,jend
	   do i=istart,iend   ! 
		call AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)
	   end do
	  end do
  case (5,11)  ! tetragonal 4/mmm, cubic m-3m
    	istart = 0
	iend = npx
   	jstart = 0
	jend = npy
	  do i=istart,iend
	   do j=jstart,i   ! 
		call AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)
	   end do
	  end do
  case (6,7)   ! npy is now npyhex !
   	istart = 0
	iend = npy
   	jstart = 0
	jend = npy
	  do j=jstart,jend
	    do i=istart,iend   ! 
		call AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)
	    end do
	  end do
  case (8)   ! npy is now npyhex !
   	istart = 0
	iend = npy
   	jstart = 0
	jend = npy
	  do j=jstart,jend
	    do i=j,iend   ! 
		call AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)
	    end do
	  end do
  case (9)   ! npy is now npyhex !
   	istart = 0
	iend = npx   ! was npy
   	jstart = 0
	jend = npx   ! was npy
	  do j=jstart,jend
	    do i=2*j,iend   ! 
		call AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)
	    end do
	  end do
  case (12)   ! npy is now npyhex !  This is the second setting of Laue group -3m.
   	istart = 0
	iend = npy
   	jstart = 0
	jend = npy/2
	  do j=jstart,jend
	    do i=2*j,iend   ! 
		call AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)
		call AddkVector(npx,npy,numk,delta,i-j,-j,usehex,mapmode)
	    end do
	  end do
end select

end subroutine

function GetSextant(x,y) result(res)
! ###################################################################
! 
!  function GetSextant
!
!  Author: Marc De Graef
!  
!  Description: determine to which sextant a point in hexagonal coordinates belongs
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   11/21/12 MDG 1.0 original 
! ###################################################################
use local

IMPLICIT NONE

real(kind=dbl),parameter	:: srt = 1.732050808   ! sqrt(3.D0)
integer(kind=irg)		:: res
real(kind=dbl),INTENT(IN):: x, y 
real(kind=dbl)			:: xx

xx = dabs(x*srt)    	! |x| sqrt(3)

if (y.ge.0) then
  if (y.ge.xx) then
	res = 0
  else
	if (x.gt.0.D0) then
	  res = 1
	else
	  res = 5
	end if
  end if
else
  if (dabs(y).ge.xx) then
	res = 3
  else
	if (x.gt.0.D0) then
	  res = 2
	else
	  res = 4
	end if
  end if
end if

end function GetSextant


! ###################################################################
! 
!  subroutine AddkVector
!
!  Author: Marc De Graef
!  
!  Description: add a k-vector for square or hexagonal grid sampling mode
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   8/20/12 MDG 1.0 original 
! ###################################################################
subroutine AddkVector(npx,npy,numk,delta,i,j,usehex,mapmode)

use local
use io
use constants
use error
use diffraction
use crystal
use crystalvars
use dynamical

IMPLICIT NONE

integer(kind=irg)         	:: npx,npy,numk,istat,i,j, ks, GetSextant
real(kind=dbl)                 :: delta, kstar(3), x, y, rr, q, iPi, XX, YY, xp, yp
logical 				:: usehex, goahead
real(kind=dbl),parameter	:: srt = 0.86602540D0   ! sqrt(3.D0)/2.D0
real(kind=dbl),parameter	:: isrt = 0.577350269D0   ! 1.D0/sqrt(3.D0)
!real(kind=dbl),parameter:: alpha = 1.346773687D0   !  sqrt(pi)/3^(1/4)
real(kind=dbl),parameter:: rtt = 1.7320508076D0   !  sqrt(3)
real(kind=dbl),parameter:: prea = 0.525037568D0   !  3^(1/4)/sqrt(2pi)
real(kind=dbl),parameter:: preb = 1.050075136D0   !  3^(1/4)sqrt(2/pi)
real(kind=dbl),parameter:: prec = 0.90689968D0   !  pi/2sqrt(3)
real(kind=dbl),parameter:: pred = 2.09439510D0   !  2pi/3

character(*)			:: mapmode

intent(IN)           :: npx,npy,i,j,delta,usehex,mapmode
intent(OUT)        :: numk

iPi = 1.D0/cPi  ! inverse of pi
goahead = .FALSE.

if (usehex) then
  x = (i - j*0.5)*delta
  y = j*delta*srt
else
  x = i*delta
  y = j*delta
end if
rr = x*x+y*y
	

if (mapmode.eq.'PlainLambert') then
	if ((.not.((i.eq.0).and.(j.eq.0))).and.(rr.le.2.D0)) then  ! skip (0,0) and stay inside circle mapped onto northern hemisphere
	     q = 0.5D0*dsqrt(4.D0-rr)
	     kstar = (/ x*q, y*q, 1.D0-0.5D0*rr /)
	     goahead = .TRUE.
	end if    

else ! mapmode = 'RoscaLambert'
! this is the more correct mapping from a uniform square grid to a sphere via an equal area map
! to the 2D circle first; it is computationally slightly longer due to the trigonometric function calls
! but it is mathematically more correct.  We do distinguish here between the square and hexagon
! projections; Note that the hexagonal case must be mirrored x <-> y with respect to the 
! analytical derivation due to a rotation of the hexagon cell.

      if (usehex) then  ! we're projecting from a hexagonal array
	if ( .not.((i.eq.0).and.(j.eq.0)) ) then  ! skip (0,0) 
! decide which sextant the point (i,j) is located in.
          ks = GetSextant(x,y)
	  select case (ks)
	  case (0,3)
	  	XX = preb*y*dcos(x*prec/y)
	  	YY = preb*y*dsin(x*prec/y)
	  case (1,4)
	  	xp = y+rtt*x
	  	yp = y*pred/xp
	  	XX = prea*xp*dsin(yp)
	  	YY = prea*xp*dcos(yp)
	  case (2,5)
	  	xp = y-rtt*x
	  	yp = y*pred/xp
	  	XX = prea*xp*dsin(yp)
	  	YY = -prea*xp*dcos(yp)	  
	  end select
	  q = XX**2+YY**2
	  kstar = (/ 0.5D0*XX*dsqrt(4.D0-q), 0.5D0*YY*dsqrt(4.D0-q),1.D0-0.5D0*q /)
          goahead = .TRUE.
	end if
      else   ! we're projecting from a square array
	if ( .not.((i.eq.0).and.(j.eq.0)) ) then  ! skip (0,0) 
! decide which equation to use  [ (8) or (9) from Rosca's paper, with r=1 ]
	     if (dabs(x).le.dabs(y)) then
                 q = 2.D0*y*iPi*dsqrt(cPi-y*y)
  	          kstar = (/ q*dsin(x*cPi*0.25D0/y), q*dcos(x*cPi*0.25D0/y), 1.D0-2.D0*y*y*iPi /)  
	     else
                  q = 2.D0*x*iPi*dsqrt(cPi-x*x)
  	          kstar = (/ q*dcos(y*cPi*0.25D0/x), q*dsin(y*cPi*0.25D0/x), 1.D0-2.D0*x*x*iPi /)  
    	     end if
             goahead = .TRUE.
         end if
      end if
end if 
 
 if (goahead.eq..TRUE.) then 
	     allocate(tail%next,stat=istat)  ! allocate new value
	     if (istat.ne.0) call FatalError('Calckvectors: unable to allocate pointer',' ')
	     tail => tail%next          ! tail points to new value
	     nullify(tail%next)          ! nullify next in new value
	     numk = numk + 1       ! keep track of number of k-vectors so far
	     if (usehex) then ! transform the hex coordinates to square-array coordinates
	       tail%i = i - j/2+mod(j,2)/2                 ! i-index of beam
	       tail%j = j                      ! j-index of beam
	     else  ! leave the square coordinates unchanged
	       tail%i = i                      ! i-index of beam
	       tail%j = j                      ! j-index of beam
	     end if 
	     call NormVec(kstar,'c')                ! normalize incident direction in cartesian space
	     kstar = kstar/mLambda                 ! divide by wavelength
	! and transform to reciprocal crystal space using the structure matrix
	     tail%k = matmul(transpose(cell%dsm),kstar)
	     tail%kn = 1.0/mLambda
end if

end subroutine AddkVector


! ###################################################################
! 
!  subroutine ApplySymmetry
! 
!  Author: Marc De Graef
!  
!  Description: Apply Laue group symmetry to input arrays
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   7/04/01 MDG 2.0 f90
!  12/05/12 MDG 2.1 adapted as a subroutine
! ###################################################################
subroutine ApplySymmetry(numk,isym,sr,srkin,npx,npy)

use local
use io
use diffraction
use multibeams
use dynamical

IMPLICIT NONE
 
integer(kind=irg)	:: ik, ipx, ipy, fly, ix, iy
real(kind=sgl)		:: my

integer(kind=irg),INTENT(IN)	:: numk, isym, npx, npy
real(kind=sgl),INTENT(INOUT) :: sr(-npx:npx,-npy:npy), srkin(-npx:npx,-npy:npy)


! depending on the symmetry group, this result needs to be copied to multiple locations ...
mess = ' Applying Laue point group symmetry to EBSD intensity array '; call Message("(A)")

symmetryloop: do ik=1,numk
   ipx = tmp%i
   ipy = tmp%j

select case (isym)
  case (1)
    ! do nothing (triclinic)
  case (2)
     sr(ipx,-ipy) = sr(ipx,ipy)
     srkin(ipx,-ipy) = srkin(ipx,ipy)
  case (4)
     sr(-ipy,ipx) = sr(ipx,ipy)
     sr(ipy,-ipx) = sr(ipx,ipy)
     sr(-ipx,-ipy) = sr(ipx,ipy)
     srkin(-ipy,ipx) = srkin(ipx,ipy)
     srkin(ipy,-ipx) = srkin(ipx,ipy)
     srkin(-ipx,-ipy) = srkin(ipx,ipy)     
  case (3,10)
     sr(-ipx,ipy) = sr(ipx,ipy)
     sr(ipx,-ipy) = sr(ipx,ipy)
     sr(-ipx,-ipy) = sr(ipx,ipy)
     srkin(-ipx,ipy) = srkin(ipx,ipy)
     srkin(ipx,-ipy) = srkin(ipx,ipy)
     srkin(-ipx,-ipy) = srkin(ipx,ipy)     
  case (5,11)
     sr(-ipx,ipy) = sr(ipx,ipy)
     sr(ipx,-ipy) = sr(ipx,ipy)
     sr(-ipx,-ipy) = sr(ipx,ipy)
     sr(ipy,ipx) = sr(ipx,ipy)
     sr(-ipy,ipx) = sr(ipx,ipy)
     sr(ipy,-ipx) = sr(ipx,ipy)
     sr(-ipy,-ipx) = sr(ipx,ipy)
     srkin(-ipx,ipy) = srkin(ipx,ipy)
     srkin(ipx,-ipy) = srkin(ipx,ipy)
     srkin(-ipx,-ipy) = srkin(ipx,ipy)
     srkin(ipy,ipx) = srkin(ipx,ipy)
     srkin(-ipy,ipx) = srkin(ipx,ipy)
     srkin(ipy,-ipx) = srkin(ipx,ipy)
     srkin(-ipy,-ipx) = srkin(ipx,ipy)
  case (6,7)  ! the transformations were computed and verified using Mathematica on 8/28/12
  ! from here on, sr and srkin are really srhex and srkinhex !
    fly = floor( float(ipy)/2.0)

    ix = ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0 )
    iy = ipx-ipy+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-fly-floor( float(ipx-fly)/2.0 )
    iy = -ipx-fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

  case (13)  ! the transformations were computed and verified using Mathematica on 9/14/12
  ! mirror is parallel to the a_1 axis
    my = float(mod(ipy,2))/2.0
    fly = floor( float(ipy)/2.0)

    ix = ipy-floor(float(ipx+fly)/2.0)
    iy = ipx+ipy/2-my
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0)
    iy = ipx-ipy+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)
    
    ix = ipx
    iy = -ipy
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-fly-floor( float(ipx-fly)/2.0)
    iy = -ipx-fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = -fly-floor(float(ipx+ipy-fly)/2.0)
    iy = -ipx+ipy/2+my
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)
    
   case (8)  ! the transformations were computed and verified using Mathematica on 9/10/12
    my = float(mod(ipy,2))
    fly = floor( float(ipy)/2.0)

    ix = ipx-fly-floor( float(ipx+2*ipy-3*fly)/2.0)
    iy = ipx+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0)
    iy = ipx-ipy+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)
    
    ix = -ipx+my
    iy = -ipy
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-fly-floor( float(ipx-fly)/2.0)
    iy = -ipx-fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipx+ipy-floor( float(ipx+ipy-fly)/2.0)
    iy = -ipx+ipy/2+my/2
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)
    
  case (9)  ! the transformations were computed and verified using Mathematica on 8/28/12
    my = float(mod(ipy,2))/2.0
    fly = floor( float(ipy)/2.0)

! {-Floor[ipy/2] - Floor[1/2 (ipx + ipy - Floor[ipy/2])],  1/2 (-2 ipx + ipy + Mod[ipy, 2])}
    ix= fly + floor( float(ipx+ipy -fly)/2.0 )
    iy = ipx-ipy/2-my
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipx-fly-floor( float(ipx+2*ipy-3*fly)/2.0)
    iy = ipx+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    sr(ipx,-ipy) = sr(ipx,ipy)
    srkin(ipx,-ipy) = srkin(ipx,ipy)

    ix = ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0)
    iy = ipx-ipy+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)
    
    ix = ipx-floor( float(ipx+2*ipy-fly)/2.0)
    iy = -ipx-fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = -ipx+2*my
    iy = -ipy
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = -fly-floor( float(ipx+ipy-fly)/2.0)
    iy = -ipx+ipy/2+my
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-fly-floor( float(ipx-fly)/2.0)
    iy = -ipx-fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = -ipx+2*my
    iy = ipy
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipx+ipy-floor( float(ipx+ipy-fly)/2.0)
    iy = (-2*ipx+ipy+2*my)/2
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-2*fly -floor( float(ipx-3*fly)/2.0)
    iy = ipx+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

! finally, take the special case of the -3m Laue group with mirror at 30 degrees from a_1.... 
 case (12)  ! the transformations were computed and verified using Mathematica on 9/6/12
    my = float(mod(ipy,2))
    fly = floor( float(ipy)/2.0)

    ix = ipx+fly-floor( float(ipx-ipy+fly)/2.0)
    iy = ipx-ipy+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0)
    iy = ipx-ipy+fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)
    
    ix = ipx-floor( float(ipx+2*ipy-fly)/2.0)
    iy = -ipx-fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = ipy-fly-floor( float(ipx-fly)/2.0)
    iy = -ipx-fly
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

    ix = -ipx+my
    iy = ipy
    sr(ix,iy) = sr(ipx,ipy)
    srkin(ix,iy) = srkin(ipx,ipy)

end select
  
! select next beam direction
   tmp => tmp%next
 end do symmetryloop
 
 
end subroutine ApplySymmetry

! ###################################################################
! 
!  subroutine CalcLgh
! 
!  Author: Marc De Graef
!  
!  Description: integrate the dynamical equations using the Bloch Wave
!  formalism.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   7/04/01 MDG 2.0 f90
! ###################################################################
subroutine CalcLgh(DMat,nn,thick,kn,gzero,Lgh,kin,debug)

use local
use io
use files
! use diffraction
! use dynamical
use constants

IMPLICIT NONE

integer         :: nn,i,j,IPIV(nn),gzero,debug
complex(kind=dbl) :: CGinv(nn,nn), Minp(nn,nn), tmp3(nn,nn),DMat(nn,nn)

real(kind=dbl)  :: kn,thick,kin
complex(kind=dbl) :: Ijk(nn,nn),Lgh(nn,nn),q,getMIWORK

integer    :: INFO, LDA, LDVR, LDVL, LWORK, JPIV(nn),MILWORK
complex(kind=dbl)    :: VL(nn,nn), CGG(nn,nn), W(nn)
real(kind=dbl)       :: RWORK(2*nn)
character            :: JOBVL, JOBVR
complex(kind=dbl),allocatable :: MIWORK(:), WORK(:)

intent(IN)      :: nn,thick, kn, gzero,DMat,debug
intent(OUT)     :: Lgh,kin

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
! 
! then get the eigenvalues and eigenvectors
 Minp = DMat
 IPIV = 0
 W = cmplx(0.D0,0.D0)
 CGG = cmplx(0.D0,0.D0)
 CGinv = cmplx(0.D0,0.D0)

 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0
 JOBVL = 'N'   ! do not compute the left eigenvectors
 JOBVR = 'V'   ! compute the right eigenvectors
! call the routine to determine the optimal workspace size
!if (debug.eq.1) write (*,*) 'sum(DynMat) = ',sum(DMat)
!if (debug.eq.1) write (*,*) 'calling zgeev'
 LWORK = -1
 allocate(WORK(2))
 WORK = cmplx(0.D0,0.D0)
 call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
 LWORK =  INT( WORK( 1 ) )  
! if (debug.eq.1) write (*,*) 'LWORK = ',LWORK
 deallocate(WORK)
! then call the eigenvalue solver
 allocate(WORK(LWORK))
 Minp = DMat
 call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
 deallocate(WORK)
!if (debug.eq.1) write (*,*) ' zgeev done '
!if (debug.eq.1) write (*,*) maxval(cabs(W))
! make a new copy of CG for the matrix inversion routines
 CGinv = CGG

!if (debug.eq.1) write (*,*) 'inverting matrix'
 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 MILWORK = -1
 call zgetri(nn,CGinv,LDA,JPIV,getMIWORK,MILWORK,INFO)
 MILWORK =  INT(real(getMIWORK))
 if (.not.allocated(MIWORK)) allocate(MIWORK(MILWORK))
 MIWORK = dcmplx(0.D0,0.D0)
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)

 if ((cabs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) write (*,*) 'Error in matrix inversion; continuing'
 deallocate(MIWORK)
!if (debug.eq.1) write (*,*) ' -> done'

! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)
! recall that alpha(1:nn) = CGinv(1:nn,gzero)

!if (debug.eq.1) write (*,*) 'computing B matrix'

! first the Ijk matrix (this is Winkelmann's B^{ij}(t) matrix)
 do i=1,nn
  do j=1,nn
   q = 2.0*cPi*thick*cmplx(aimag(W(i))+aimag(W(j)),real(W(i))-real(W(j)))
    Ijk(i,j) = conjg(CGinv(i,gzero)) * (1.0-exp(-q))/q * CGinv(j,gzero)
  end do
 end do
!if (debug.eq.1) write (*,*) ' -> done'
!
!if (debug.eq.1) write (*,*) 'matmul operations'
!
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

! we no longer need CGinv, so we'll use the array to store the conjugate of CGG
!CGinv = conjg(CGG)
!do i=1,nn
!  do j=1,nn
!    Lgh(i,j) = sum( tmp3(i,1:nn) * CGinv(j,1:nn) )
!  end do
!end do


! we'll approximate this by the sum of the diagonal entries 
 tmp3 = matmul(transpose(CGG),conjg(CGG)) 
 kin = 1.D0-real(sum(Ijk * tmp3))
!if (debug.eq.1) write (*,*) ' -> done'

!do i=1,nn
!  do j=1,nn
!    tmp3(i,j) = sum( CGG(1:nn,i) * CGinv(1:nn,j) )
!  end do
!end do
!kin = 1.D0-real(sum(Ijk * tmp3))
!


!deallocate(CGinv,Minp,diag,tmp3)

end subroutine


! ###################################################################
! 
!  subroutine CalcLgh2
! 
!  Author: Marc De Graef
!  
!  Description: integrate the dynamical equations using the Bloch Wave
!  formalism.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   7/04/01 MDG 2.0 f90
! ###################################################################
subroutine CalcLgh2(DMat,Lgh,thick,kn,nn,gzero,kin,debug)

use local
use io
use files
! use diffraction
! use dynamical
use constants

IMPLICIT NONE

integer         :: nn,i,j,IPIV(nn),gzero,debug, low, igh, ierr
complex(kind=dbl) :: CGinv(nn,nn), Minp(nn,nn), tmp3(nn,nn),DMat(nn,nn)

real(kind=dbl)  :: kn,thick,kin
complex(kind=dbl) :: Ijk(nn,nn),Lgh(nn,nn),q,getMIWORK

integer(kind=irg)    :: INFO, LDA, LDVR, LDVL,  JPIV(nn),MILWORK
complex(kind=dbl)    ::  CGG(nn,nn), W(nn)
complex(kind=dbl),allocatable :: MIWORK(:)
real(kind=dbl),allocatable :: ortr(:), orti(:)

! Uncomment the following two lines if the EISPACK cg routine is needed
real(kind=dbl)       :: ar(nn,nn), ai(nn,nn), wr(nn), wi(nn), vr(nn,nn), vi(nn,nn), scle(nn)

intent(IN)      :: nn,thick, kn, gzero,DMat,debug
intent(OUT)     :: Lgh,kin

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
! 
! then get the eigenvalues and eigenvectors
!if (debug.eq.1) write (*,*) 'entered CalcLgh2'
 Minp = DMat
!if (debug.eq.1) write (*,*) 'copied DMat'
 
 IPIV = 0
 W = cmplx(0.D0,0.D0)
 CGG = cmplx(0.D0,0.D0)
 CGinv = cmplx(0.D0,0.D0)

 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0
! the eispack routines use the real and imaginary parts of arrays as separate arrays instead
! of as complex variable ...  
!if (debug.eq.1) write (*,*) 'splitting real and imaginary'

ar = dble(Minp)
ai = aimag(Minp)

! first balance the matrix
!if (debug.eq.1) write (*,*) 'cbal'

call cbal(nn, ar, ai, low, igh, scle)
! transform the upper Hessenberg form
allocate(ortr(igh), orti(igh))
!if (debug.eq.1) write (*,*) 'corth'
call corth(nn, low, igh, ar, ai, ortr, orti)
! get eigenvalues and eigenvectors
!if (debug.eq.1) write (*,*) 'comqr2'
call comqr2(nn, low, igh, ortr, orti, ar, ai, wr, wi, vr, vi, ierr )
if ( ierr.ne.0) then 
  write (*,*) 'Error in comqr2 eispack routine'
  stop
end if
deallocate(ortr, orti)
! undo the cbal transformation
!if (debug.eq.1) write (*,*) 'cbabk2'
call cbabk2(nn, low, igh, scle, nn, vr, vi)

! return to complex variables
W = cmplx(wr, wi)
CGG = cmplx(vr, vi)
 

 CGinv = CGG

!if (debug.eq.1) write (*,*) 'inverting matrix'
 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 MILWORK = -1
 call zgetri(nn,CGinv,LDA,JPIV,getMIWORK,MILWORK,INFO)
 MILWORK =  INT(real(getMIWORK))
 if (.not.allocated(MIWORK)) allocate(MIWORK(MILWORK))
 MIWORK = dcmplx(0.D0,0.D0)
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)

 if ((cabs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) write (*,*) 'Error in matrix inversion; continuing'
 deallocate(MIWORK)
!if (debug.eq.1) write (*,*) ' -> done'

! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)
! recall that alpha(1:nn) = CGinv(1:nn,gzero)

!if (debug.eq.1) write (*,*) 'computing B matrix'

! first the Ijk matrix (this is Winkelmann's B^{ij}(t) matrix)  [old code]
 do i=1,nn
  do j=1,nn
   q = 2.0*cPi*thick*cmplx(aimag(W(i))+aimag(W(j)),real(W(i))-real(W(j)))
    Ijk(i,j) = conjg(CGinv(i,gzero)) * (1.0-exp(-q))/q * CGinv(j,gzero)
  end do
end do

!if (debug.eq.1) write (*,*) ' -> done'

!if (debug.eq.1) write (*,*) 'matmul operations'

! then the summations for Lgh and kin
!tmp3 = matmul(CGG,transpose(Ijk)) 
!Lgh = matmul(tmp3,transpose(conjg(CGG)))

! there might be a problem with the Absoft implementation of the 
! matmul routine...  So let's do this multiplication explicitly...
do i=1,nn
  do j=1,nn
     tmp3(i,j) = sum( CGG(i,1:nn) * Ijk(j,1:nn) )
  end do
end do

! we no longer need CGinv, so we'll use the array to store the conjugate of CGG
CGinv = conjg(CGG)
do i=1,nn
  do j=1,nn
    Lgh(i,j) = sum( tmp3(i,1:nn) * CGinv(j,1:nn) )
  end do
end do


! we'll approximate this by the sum of the diagonal entries 
! tmp3 = matmul(transpose(CGG),conjg(CGG)) 
! kin = 1.D0-real(sum(Ijk * tmp3))
!if (debug.eq.1) write (*,*) ' -> done'

do i=1,nn
  do j=1,nn
    tmp3(i,j) = sum( CGG(1:nn,i) * CGinv(1:nn,j) )
  end do
end do
kin = 1.D0-real(sum(Ijk * tmp3))
!


!deallocate(CGinv,Minp,diag,tmp3)

end subroutine CalcLgh2

     
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
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   7/04/01 MDG 2.0 f90
! 12/12/12 MDg 3.0 attempt at actual numerical integration for the matrix elements I_jk
! ###################################################################
subroutine CalcLgh3(DMat,Lgh,thick,kn,nn,gzero,kin,debug,depthstep,lambdaE,izz)

use local
use io
use files
! use diffraction
! use dynamical
use constants

IMPLICIT NONE

integer         :: nn,i,j,IPIV(nn),gzero,debug, low, igh, ierr, izz, iz
complex(kind=dbl) :: CGinv(nn,nn), Minp(nn,nn), tmp3(nn,nn),DMat(nn,nn)

real(kind=dbl)  :: kn,thick,kin, dz, tpi, dzt, depthstep
complex(kind=dbl) :: Ijk(nn,nn), Lgh(nn,nn), q, getMIWORK, qold

integer(kind=irg)    :: INFO, LDA, LDVR, LDVL,  JPIV(nn), MILWORK, lambdaE(izz)
complex(kind=dbl)    ::  CGG(nn,nn), W(nn)
complex(kind=dbl),allocatable :: MIWORK(:)
real(kind=dbl),allocatable :: ortr(:), orti(:)

! Uncomment the following two lines if the EISPACK cg routine is needed
real(kind=dbl)       :: ar(nn,nn), ai(nn,nn), wr(nn), wi(nn), vr(nn,nn), vi(nn,nn), scle(nn)

intent(IN)      :: nn,thick, kn, gzero,DMat,debug
intent(OUT)     :: Lgh,kin

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
! 
! then get the eigenvalues and eigenvectors
!if (debug.eq.1) write (*,*) 'entered CalcLgh2'
 Minp = DMat
!if (debug.eq.1) write (*,*) 'copied DMat'
 
 IPIV = 0
 W = cmplx(0.D0,0.D0)
 CGG = cmplx(0.D0,0.D0)
 CGinv = cmplx(0.D0,0.D0)

 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0
! the eispack routines use the real and imaginary parts of arrays as separate arrays instead
! of as complex variable ...  
!if (debug.eq.1) write (*,*) 'splitting real and imaginary'

ar = dble(Minp)
ai = aimag(Minp)

! first balance the matrix
!if (debug.eq.1) write (*,*) 'cbal'

call cbal(nn, ar, ai, low, igh, scle)
! transform the upper Hessenberg form
allocate(ortr(igh), orti(igh))
!if (debug.eq.1) write (*,*) 'corth'
call corth(nn, low, igh, ar, ai, ortr, orti)
! get eigenvalues and eigenvectors
!if (debug.eq.1) write (*,*) 'comqr2'
call comqr2(nn, low, igh, ortr, orti, ar, ai, wr, wi, vr, vi, ierr )
if ( ierr.ne.0) then 
  write (*,*) 'Error in comqr2 eispack routine'
  stop
end if
deallocate(ortr, orti)
! undo the cbal transformation
!if (debug.eq.1) write (*,*) 'cbabk2'
call cbabk2(nn, low, igh, scle, nn, vr, vi)

! return to complex variables
W = cmplx(wr, wi)
CGG = cmplx(vr, vi)
 

 CGinv = CGG

!if (debug.eq.1) write (*,*) 'inverting matrix'
 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 MILWORK = -1
 call zgetri(nn,CGinv,LDA,JPIV,getMIWORK,MILWORK,INFO)
 MILWORK =  INT(real(getMIWORK))
 if (.not.allocated(MIWORK)) allocate(MIWORK(MILWORK))
 MIWORK = dcmplx(0.D0,0.D0)
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)

 if ((cabs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) write (*,*) 'Error in matrix inversion; continuing'
 deallocate(MIWORK)
!if (debug.eq.1) write (*,*) ' -> done'

! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)
! recall that alpha(1:nn) = CGinv(1:nn,gzero)

!if (debug.eq.1) write (*,*) 'computing B matrix'

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

!if (debug.eq.1) write (*,*) ' -> done'

!if (debug.eq.1) write (*,*) 'matmul operations'

! then the summations for Lgh and kin
!tmp3 = matmul(CGG,transpose(Ijk)) 
!Lgh = matmul(tmp3,transpose(conjg(CGG)))

! there might be a problem with the Absoft implementation of the 
! matmul routine...  So let's do this multiplication explicitly...
do i=1,nn
  do j=1,nn
     tmp3(i,j) = sum( CGG(i,1:nn) * Ijk(j,1:nn) )
  end do
end do

! we no longer need CGinv, so we'll use the array to store the conjugate of CGG
CGinv = conjg(CGG)
do i=1,nn
  do j=1,nn
    Lgh(i,j) = sum( tmp3(i,1:nn) * CGinv(j,1:nn) )
  end do
end do


! we'll approximate this by the sum of the diagonal entries 
! tmp3 = matmul(transpose(CGG),conjg(CGG)) 
! kin = 1.D0-real(sum(Ijk * tmp3))
!if (debug.eq.1) write (*,*) ' -> done'

! there may be an issue with this part of the routine ... 
do i=1,nn
  do j=1,nn
    tmp3(i,j) = sum( CGG(1:nn,i) * CGinv(1:nn,j) )
  end do
end do
kin = 1.D0-real(sum(Ijk * tmp3))
!


!deallocate(CGinv,Minp,diag,tmp3)

end subroutine CalcLgh3

     
