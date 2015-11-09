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
! EMsoft:EMECPSMMultiLayer.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMECPSMMultiLayer
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Zone axis electron channeling patterns for two layer structure using scattering matrix
!
!> @date 11/25/14 SS 1.0 f90
!> @date 08/31/15 SS 1.1 Revisions+OpenMP support
!--------------------------------------------------------------------------
program EMECPSMMultiLayer

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)            :: nmldeffile
character(fnlen)            :: progdesc
character(fnlen)            :: progname
type(ECPNameListType)       :: ecpnl

! deal with the command line arguments, if any
nmldeffile = 'EMECPMultiLayer.nml'
progname = 'EMECPSMMultiLayer.f90'
progdesc = 'Electron channeling pattern computation for two layer structure using Scattering Matrix approach'

! set all global variables and print information
call EMsoft(progname, progdesc)

! deal with command line arguments if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

!deal with the namelist file
call GetECPNameList(nmldeffile,ecpnl)

! perform the zone axis computations
call ECPMultiLayerSM(ecpnl, progname)

end program EMECPSMMultiLayer

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECPMultiLayerSM
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief compute a zone axis channeling pattern for two layer structure
!
!> @note This program computes the channeling patterns for a thin film of top of
!> a substrate. The contribution of the substrate becomes important at small film
!> thicknesses and we need to take that into account. This program uses the
!> scattering matrix approach. Later on the program will be integrated with the
!> bloch wave approach program to have one single program with a handle to switch.
!> The equations being used for this computations were derived assuming that each
!> diffracted beam is an incident beam to the substrate.
!
!> @param ecpnl namelist file
!
!> @date 11/25/14  SS 1.0 new routine
!> @date 08/31/15  SS 1.1 revisions+OpenMP support
!--------------------------------------------------------------------------

subroutine ECPMultiLayerSM(ecpnl, progname)


use local
use typedefs
use NameListTypedefs
use crystal
use constants
use symmetry
use Lambert
use initializers
use constants
use gvectors
use kvectors
use error
use io
use files
use diffraction
use rotations
use MBmodule
use math
use omp_lib

IMPLICIT NONE

type(ECPNameListType),INTENT(IN)        :: ecpnl
character(fnlen),INTENT(IN)             :: progname

integer(kind=irg)                       :: TID, nthreads
integer(kind=irg)                       :: numEbins,numzbins,nx,ny,totnum_el
!integer(kind=8)                        :: totnum_el this will be done later
real(kind=dbl)                          :: EkeV,Ehistmin,Ebinsize,depthmax,depthstep,sig,omega
real(kind=sgl)                          :: totnum_el_det
integer(kind=irg)                       :: iz,izz,ii,jj,kk,ll ! all loop variables
integer(kind=irg)                       :: ipx,ipy,iequiv(2,12),nequiv,gzero
integer(kind=irg),allocatable           :: accum_z(:,:,:,:), accum_e(:,:,:)
integer(kind=irg)                       :: istat,filmthickness,substhickness
real(kind=sgl),allocatable              :: lambdaZ(:),thick(:)
real(kind=dbl)                          :: io_real(5)
integer(kind=irg)                       :: io_int_sgl(1)

character(fnlen)                        :: oldprogname
character(fnlen)                        :: xtalname
character(8)                            :: MCscversion
character(4)                            :: MCmode

logical                                 :: verbose

type(unitcell),pointer                  :: cell_film
type(unitcell),pointer                  :: cell_subs
type(orientation),pointer               :: orel
real(kind=sgl)                          :: TTinv(3,3),kg(3),kg1(3)
real(kind=dbl)                          :: eWavelength_subs,eWavelength_film
type(gnode)                             :: rlp_film
type(gnode)                             :: rlp_subs
type(DynType)                           :: Dyn_film
type(DynType)                           :: Dyn_subs
integer(kind=irg)                       :: pgnum_film,pgnum_subs,isym_film,isym_subs,numset_film,numset_subs
type(reflisttype),pointer               :: reflist_film,reflisttmp,firstw_film,rltmpb
type(reflisttype),pointer               :: reflist_subs,firstw_subs,rltmpa
type(refliststrongsubstype),pointer     :: refliststrong_subs,refliststrong_subs_tmp
real(kind=dbl)                          :: kn,K0c(3),FNc(3)

type(kvectorlist),pointer               :: khead,ktmp
integer(kind=irg)                       :: numk,nref_film,nns_film,nnw_film
integer(kind=irg)                       :: nref_subs,nns_subs,nnw_subs,numg
integer(kind=irg),allocatable           :: kij(:,:)
real(kind=sgl),allocatable              :: klist(:,:)
real(kind=sgl)                          :: rotmat(3,3),FN(3),k0(3),dmin,k1(3),natsubs
type(BetheParameterType)                :: BetheParameters
real(kind=sgl),allocatable              :: sr_film(:,:),sr_subs(:,:)
complex(kind=dbl),allocatable           :: DynMat_film(:,:),Dynmat_subs(:,:),mexp(:,:),ampl2(:),ampl(:)
complex(kind=dbl),allocatable           :: Sghfilm(:,:),Lghfilm(:,:),Sghfilmtmp(:,:,:),Lghfilmtmp(:,:,:)
complex(kind=dbl),allocatable           :: Sghsubs(:,:),Lghsubs(:,:),Sghsubstmp(:,:,:),Lghsubstmp(:,:,:)
complex(kind=dbl),allocatable           :: S0_subs(:)
complex(kind=dbl)                       :: czero
integer(kind=irg),allocatable           :: nat_film(:),nat_subs(:)
real(kind=dbl)                          :: xgp,subscontribution
real(kind=dbl),allocatable              :: sigmagg(:,:)

call Message('opening '//trim(ecpnl%energyfile), frm = "(A)" )

open(dataunit,file=trim(ecpnl%energyfile),status='unknown',form='unformatted')

! lines from CTEMMCCL.f90... these are the things we need to read in...
! write (dataunit) progname
!! write the version number
! write (dataunit) MCscversion
!! then the name of the crystal data file
! write (dataunit) xtalname
!! energy information etc...
! write (dataunit) numEbins, numzbins, numsx, numsy, , totnum_el
! write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
! write (dataunit) sig, omega
! write (dataunit) MCmode
!! and here are the actual results
! write (dataunit) accum_e
! write (dataunit) accum_z

read (dataunit) oldprogname
read (dataunit) MCscversion
read (dataunit) xtalname

read(dataunit) numEbins, numzbins, nx, ny, totnum_el
nx = (nx - 1)/2
ny = (ny - 1)/2

read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")

read (dataunit) sig, omega
read (dataunit) MCmode

allocate(accum_e(numEbins,-nx:nx,-nx:nx),accum_z(numEbins,numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)

read(dataunit) accum_e
! actually, we do not yet need the accum_e array for ECP. This will be removed with an updated version of the MC code
! but we need to skip it in this unformatted file so that we can read the accum_z array ...
deallocate(accum_e)

read(dataunit) accum_z    ! we only need this array for the depth integrations
close(dataunit,status='keep')
call Message(' -> completed reading '//trim(ecpnl%energyfile), frm = "(A)")

!=============================================
! completed reading monte carlo file
!=============================================

! calculate the weight factors from the Monte Carlo input file
allocate(lambdaZ(numzbins),thick(numzbins),stat=istat)
totnum_el_det = float(sum(accum_z(:,:,:,:)))

do ii = 1,numzbins
    lambdaZ(ii) = float(sum(accum_z(:,ii,:,:)))/totnum_el_det
    thick(ii) = depthstep*(ii-1)
end do
call Message(' -> completed calculating lambda values from file '//trim(ecpnl%energyfile), frm = "(A)")
deallocate(accum_z)
filmthickness = nint(ecpnl%filmthickness)
substhickness = numzbins - filmthickness

if(ecpnl%filmthickness .gt. depthmax) call FatalError('ECPpattern','Lambda values not available for all film thicknesses')

allocate(cell_film)

verbose = .TRUE.
call Initialize_Cell(cell_film,Dyn_film,rlp_film,ecpnl%xtalname,ecpnl%dmin,sngl(1000.0*EkeV),verbose)
numset_film = cell_film%ATOM_ntype
allocate(nat_film(numset_film),stat=istat)

allocate(cell_subs)

call Initialize_Cell(cell_subs,Dyn_subs,rlp_subs,ecpnl%xtalname2,ecpnl%dmin,sngl(1000.0*EkeV),verbose)
numset_subs = cell_subs%ATOM_ntype
allocate(nat_subs(numset_subs),stat=istat)

jj=0
do ii=1,32
    if (SGPG(ii).le.cell_film%SYM_SGnum) jj=ii
end do
isym_film = jj
pgnum_film = jj

jj=0
do ii=1,32
    if (SGPG(ii).le.cell_subs%SYM_SGnum) jj=ii
end do
isym_subs = jj
pgnum_subs = jj

!=============================================
! completed initializing the two crystals
!=============================================

allocate(orel)
orel%gA = float(ecpnl%gF)
orel%tA = float(ecpnl%tF)
orel%gB = float(ecpnl%gS)
orel%tB = float(ecpnl%tS)

! check for orthonormality using zone equation
if (sum(orel%tA*orel%gA) .ne. 0.0) call FatalError('ECPpatternfos','Plane does not contain direction (Film)')
if (sum(orel%tB*orel%gB) .ne. 0.0) call FatalError('ECPpatternfos','Plane does not contain direction (Substrate)')
! and compute the (inverse) transformation matrix
TTinv = ComputeOR(orel, cell_film, cell_subs, 'BA')

rotmat(:,1) = (/1.0,0.0,0.0/)
rotmat(:,1) = rotmat(:,1)/NORM2(rotmat(:,1))

rotmat(:,2) = (/0.0,1.0,0.0/)
rotmat(:,2) = rotmat(:,2)/NORM2(rotmat(:,2))

rotmat(:,3) = (/0.0,0.0,1.0/)
rotmat(:,3) = rotmat(:,3)/NORM2(rotmat(:,3))

nullify(khead)
nullify(ktmp)

call CalckvectorsECP(khead,cell_film,rotmat,ecpnl%thetac,ecpnl%npix,ecpnl%npix,numk)

io_int_sgl(1) = numk
call WriteValue('# independent beam directions to be considered = ', io_int_sgl, 1, "(I8)")

allocate(kij(2,numk),klist(3,numk),stat=istat)

ktmp => khead
!kij(1:2,1) = (/ ktmp%i, ktmp%j /)
do ii = 1,numk
    kij(1:2,ii) = (/ ktmp%i, ktmp%j /)
    klist(1:3,ii) = ktmp%k
    ktmp => ktmp%next
end do

!==============================================================
! completed initializing the list of incident wavevectors
!==============================================================
call Set_Bethe_Parameters(BetheParameters,.TRUE.)
dmin = ecpnl%dmin

allocate(sr_film(-ecpnl%npix:ecpnl%npix,-ecpnl%npix:ecpnl%npix),stat=istat)
sr_film = 0.0
allocate(sr_subs(-ecpnl%npix:ecpnl%npix,-ecpnl%npix:ecpnl%npix),stat=istat)
sr_subs = 0.0

!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
! point to the first beam direction
 czero = dcmplx(0.D0,0.D0)

!$OMP PARALLEL PRIVATE(TID,nthreads,k0,FN,kn,nref_film,firstw_film,reflist_film) &
!$OMP& PRIVATE(nns_film,nnw_film,DynMat_film,Sghfilmtmp,Lghfilmtmp,Sghfilm,Lghfilm,S0_subs,refliststrong_subs) &
!$OMP& PRIVATE(jj,rltmpa,kg,kg1,reflist_subs,nref_subs,firstw_subs,nns_subs,nnw_subs,ipx,ipy) &
!$OMP& PRIVATE(DynMat_subs,Sghsubstmp,Lghsubstmp,Sghsubs,Lghsubs,sigmagg,refliststrong_subs_tmp) &
!$OMP& SHARED(lambdaZ,sr_film,filmthickness,ecpnl,ii,cell_film,cell_subs,BetheParameters,dmin,rlp_film) &
!$OMP& SHARED(numset_film,numset_subs,nat_film,nat_subs,substhickness,sr_subs,czero,numzbins,thick,TTinv,rlp_subs)

TID = OMP_GET_THREAD_NUM()
nthreads = OMP_GET_NUM_THREADS()
if(TID .eq. 0) write(*,*)'Number of threads = ',nthreads 

!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC)
beamloop: do ii = 1,numk

     k0 = klist(1:3,ii) 
     FN = ecpnl%gF
     call NormVec(cell_film,FN,'r')
     kn = CalcDot(cell_film,FN,k0,'r')
     nullify(reflist_film)
     call Initialize_ReflectionList(cell_film, reflist_film, BetheParameters, FN, k0, dmin, nref_film)

! determine strong and weak reflections
     call Apply_BethePotentials(cell_film, reflist_film, firstw_film, BetheParameters, nref_film, nns_film, nnw_film)

     if (allocated(DynMat_film)) deallocate(DynMat_film)
     if (allocated(Sghfilmtmp)) deallocate(Sghfilmtmp)
     if (allocated(Lghfilmtmp)) deallocate(Lghfilmtmp)
     if (allocated(Sghfilm)) deallocate(Sghfilm)
     if (allocated(Lghfilm)) deallocate(Lghfilm)
     if (allocated(S0_subs)) deallocate(S0_subs)

     allocate(DynMat_film(nns_film,nns_film),stat=istat)
     call GetDynMat(cell_film, reflist_film, firstw_film, rlp_film, DynMat_film, nns_film, nnw_film)
    
     allocate(Sghfilmtmp(nns_film,nns_film,numset_film), Sghfilm(nns_film,nns_film),&
            Lghfilmtmp(nns_film,nns_film,filmthickness),Lghfilm(nns_film,nns_film),stat=istat)
     allocate(S0_subs(1:nns_film),stat =istat)



     Sghfilmtmp = czero
     Sghfilm = czero
     Lghfilmtmp = czero
     Lghfilm = czero
     nat_film = 0
     call CalcSgh(cell_film,reflist_film,nns_film,numset_film,Sghfilmtmp,nat_film)
     Sghfilm = sum(Sghfilmtmp,3)

!   call CalcLgh(DynMat_film,Lghfilm,dble(filmthickness),kn,nns_film,gzero,1.D0,lambdaZ,numzbins)

     call CalcLghfilm(cell_film%mLambda,nns_film,numzbins,thick,ktmp%kn,Lghfilmtmp,DynMat_film,lambdaZ,gzero,filmthickness,S0_subs)
     Lghfilm = Lghfilmtmp(:,:,filmthickness)/dcmplx(ecpnl%filmthickness,0.D0)
     ipx = kij(1,ii)
     ipy = kij(2,ii)

     sr_film(ipx,ipy) = real(sum(Lghfilm(1:nns_film,1:nns_film)*&
                Sghfilm(1:nns_film,1:nns_film)))/2.0

!print*,'Film contribution = ',sr_film(ipx,ipy)
!print*,S0_subs
!stop
!if(TID .eq. 3) print*,'Film contribution = ',sr_film(ipx,ipy) 

!==================================================================================
! film contribution done
!==================================================================================

! THIS IS THE FULL APPROACH IN WHICH THE DIFFRACTED BEAMS FROM DIFFERENT INCIDENT BEAM ON SUBSTRATE INTERACT


! this routine gives the list of reflections and corresponding dynamical matrices for every incident direction on the substrate
! the members of the linked list are 
! kg:       the incident beam vector converted to the substrate reference frame
! DynMat:   Dynamical matrix for that particular beam direction
! hlist:    list of reflections (hkl)
! nns:      number of strong beams

     FN = ecpnl%gS
     call NormVec(cell_subs,FN,'r')

     !if(allocated(refliststrong_subs)) deallocate(refliststrong_subs)
    
     call GetStrongBeamsSubs(cell_film,cell_subs,reflist_film,refliststrong_subs,&
         k0,dble(FN),nns_film,dmin,TTinv,rlp_subs,1.D0)

!if(TID .eq. 0) then    
!     refliststrong_subs_tmp => refliststrong_subs
!     do jj = 1,nns_film
        !print*,'No. of strong beams for last incident beam = ',refliststrong_subs_tmp%nns
!        call Initialize_ReflectionList(cell_subs,reflist_subs,BetheParameters,sngl(FN),&
!        sngl(refliststrong_subs_tmp%kg),dmin,nref_subs)

!        call Apply_BethePotentials(cell_subs, reflist_subs, firstw_subs, BetheParameters, nref_subs, nns_subs, nnw_subs)
!        if(allocated(Sghsubs)) deallocate(Sghsubs)
        
!        allocate(Sghsubs(nns_subs,nns_subs))
!        Sghsubs = czero

!        call CalcSgh(cell_subs,reflist_subs,nns_subs,numset_subs,Sghsubs,nat_subs)

!        if(allocated(DynMat_subs)) deallocate(DynMat_subs)
!        if(allocated(mexp)) deallocate(mexp)
!        allocate(DynMat_subs(nns_subs,nns_subs),mexp(nns_subs,nns_subs))
        !print*,'No. of strong beams for this beam inline = ',nns_subs
        !print*,''
        
!        call GetDynMat(cell_subs,reflist_subs,firstw_subs,rlp_subs,DynMat_subs,nns_subs,nnw_subs)
        
!        DynMat_subs = DynMat_subs*cmplx(0.D0,cPi*cell_subs%mLambda)
!        call MatrixExponential(DynMat_subs,mexp,1.D0,'Pade',nns_subs)
!        if(allocated(Lghsubs)) deallocate(Lghsubs)
!        if(allocated(ampl)) deallocate(ampl)
!        if(allocated(ampl2)) deallocate(ampl2)

!        allocate(Lghsubs(nns_subs,nns_subs),ampl(nns_subs),ampl2(nns_subs))
!        ampl = czero
!        ampl2 = czero
!        Lghsubs = czero
!        ampl(1) = S0_subs(jj)

!        do kk = 1,substhickness
!            ampl2 = matmul(mexp,ampl)
!            Lghsubs = Lghsubs + spread(ampl2(1:nns_subs),dim=2,ncopies=nns_subs)*&
!            spread(conjg(ampl2(1:nns_subs)),dim=1,ncopies=nns_subs)*lambdaZ(kk+filmthickness)
!            ampl = ampl2
!        end do

       
!       print*,Lghsubs(1,1:5)
!subscontribution = subscontribution + real(sum(Lghsubs*Sghsubs))/4.0

!if(real(sum(Lghsubs*Sghsubs))/float(substhickness) .gt. 100) then
!print*,maxval(abs(Lghsubs)),maxval(abs(Sghsubs))
!end if

!print*,'Inline =',real(sum(Lghsubs*Sghsubs))/float(substhickness)        
!print*,'Inline = ', mexp(1,1:5)
!print*,'subroutine = ',refliststrong_subs_tmp%DynMat(1,1:5)
!print*,''

!        refliststrong_subs_tmp => refliststrong_subs_tmp%next
!     end do
!stop
!end if

     if (allocated(sigmagg)) deallocate(sigmagg)
     allocate(sigmagg(1:nns_film,1:nns_film),stat=istat)
     if (istat .ne. 0) call FatalError("STOP:"," cannot allocate pointer")
     sigmagg = 0.D0

     call CalcSigmaggSubstrate(cell_subs,nns_film,refliststrong_subs,S0_subs,sigmagg,&
         filmthickness,substhickness,lambdaZ,thick,numzbins)

     sr_subs(ipx,ipy) = sum(sigmagg)

!if(TID .eq. 0) then
!print*,'Substrate contribution = ',sr_subs(ipx,ipy) 
!print*,'Symmetric check',sigmagg(11,17),sigmagg(17,11)
!print*,''
!subscontribution = 0.0
!do ll = 1,nns_film
!do kk = 1,nns_film
!if(ll .eq. kk) subscontribution = subscontribution  + sigmagg(ll,kk)
!end do
!end do
!print*,'Diagonal term = ',subscontribution
!print*,'Off diagonal term = ',sum(sigmagg) - subscontribution
!end if

     if (mod(ii,100) .eq. 0) then
        io_int_sgl = ii
        call WriteValue('Completed beam # ', io_int_sgl, 1, "(I8)")
     end if

 
end do beamloop

!$OMP END DO
!$OMP END PARALLEL

open(unit=13,file='test_film_10.txt',action='write')
open(unit=14,file='test_subs_10.txt',action='write')


do ii = -ecpnl%npix,ecpnl%npix
   do jj = -ecpnl%npix,ecpnl%npix
      write(13,'(F15.6)',advance='no')sr_film(ii,jj)
      write(14,'(F15.6)',advance='no')sr_subs(ii,jj) 
   end do
   write(13,*)''
   write(14,*)''
end do

close(13)
close(14)

end subroutine ECPMultiLayerSM

!--------------------------------------------------------------------------
!
! SUBROUTINE:CalcLghfilm
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief integrate the Bloch wave function over the foil thickness
!
!> @param nn number of strong beams
!> @param nt number of thickness values
!> @param thick array of thickness values
!> @param kn normal component of incident wave vector
!> @param gzero index of zero beam (should always be the first one; legacy parameter)
!> @param Lgh output array
!> @param Kossel Kossel intensity array
!
!> @date 11/18/13  MDG 1.0 major rewrite from older ECP program; merged with ECPz
!> @date 03/03/14  MDG 2.0 version that works with the scattering matrix...
!--------------------------------------------------------------------------
recursive subroutine CalcLghfilm(mLambda,nn,nt,thick,kn,Lgh,DynMat,lambdaZ,gzero,filmthickness,S0) !,Kossel)

use local
use io
use files
use diffraction
use MBmodule
use constants
use math

IMPLICIT NONE

real(kind=dbl),INTENT(IN)           :: mLambda
integer(kind=irg),INTENT(IN)        :: nn
integer(kind=irg),INTENT(IN)        :: nt
real(kind=sgl),INTENT(IN)           :: thick(nt)
real(kind=dbl),INTENT(IN)           :: kn
complex(kind=dbl),INTENT(IN)        :: DynMat(nn,nn)
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn,filmthickness)
real(kind=sgl),INTENT(IN)           :: lambdaZ(nt)
integer(kind=irg),INTENT(IN)        :: filmthickness
complex(kind=dbl),INTENT(OUT)       :: S0(nn)
integer(kind=irg),INTENT(IN)        :: gzero
!real(kind=dbl),INTENT(OUT)          :: Kossel(nt)

integer                             :: i
complex(kind=dbl),allocatable       :: Minp(:,:),Azz(:,:),ampl(:),ampl2(:)
real(kind=sgl)                      :: dthick

allocate(Minp(nn,nn),Azz(nn,nn),ampl(nn),ampl2(nn))

dthick = thick(2)-thick(1)
Minp = DynMat * dcmplx(0.D0,cPi * mLambda)
call MatrixExponential(Minp, Azz, dble(dthick), 'Pade', nn)
ampl = dcmplx(0.D0,0.D0)
ampl(1) = dcmplx(1.0D0,0.D0)
Lgh = dcmplx(0.D0,0.D0)
S0(1:nn) = ampl(1:nn)
! add some thickness handling here !!!
! thge integration uses a small step size, but there might
! be only a small number of thicknesses for which output is

do i=1,filmthickness
    ampl2 = matmul(Azz,ampl)
    if (i .eq. 1) then
        Lgh(1:nn,1:nn,i) = spread(ampl2(1:nn),dim=2,ncopies=nn)*spread(conjg(ampl2(1:nn)),dim=1,ncopies=nn)*lambdaZ(i)

    else
        Lgh(1:nn,1:nn,i) = Lgh(1:nn,1:nn,i-1)+spread(ampl2(1:nn),dim=2,ncopies=nn)*spread(conjg(ampl2(1:nn)),dim=1,ncopies=nn)*&
                   lambdaZ(i)
    end if
!    do ih=1,nn
!      do ig=1,nn
!        Lgh(ig,ih,i) = Lgh(ig,ih,i) + ampl2(ih) * conjg(ampl2(ig))
!      end do
!    end do
    ampl = ampl2
    S0(1:nn) = ampl(1:nn)
end do


deallocate(Minp,Azz,ampl,ampl2)

end subroutine CalcLghfilm

