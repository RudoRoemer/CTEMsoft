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
! CTEMsoft2013:CTEMKossel.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMlacbed 
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Kossel patterns, taken from master EBSD pattern simulation, with adjustment for 
!> Bloch wave part ...
!
!> @todo implement OpenMP multithreading for the actual computation part; requires modifications
!> in CTEMlib.a routines (mostly THREADPRIVATE commands in several modules)
!
!> @date 04/08/11 MDG 1.0 early version
!--------------------------------------------------------------------------
program CTEMKossel

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
use postscript
use timing

IMPLICIT NONE

character(2)     :: srza    ! = SR for systematic row and ZA for zone axis orientation
character(4),parameter :: methods(3) = (/'SMTE','RK90','BWAV'/)
character(20)    :: fname, gname   ! = output filename
real(kind=dbl)   :: ctmp(192,3),arg
integer(kind=irg)          :: ik,ih,npx,npy,ipx,ipy, & ! counters
                    k(3), & ! incident beam direction
                    numk,ktstep, & ! number of independent incident beam directions
                    ga(3),& ! first reciprocal lattice vector
                    gb(3),& ! second reciprocal lattice vector (not used for systematic row)
                    fn(3),& ! foil normal
                    nn,nthick,& ! total number of reflections in computation, and thicknesses
                    isym,&  ! number of 2D point group symmetry
                    i,j,ir,nat(100),kk(3),gshort(3),NN,maxholz,minholz,numrefs,N, &
                    numset,n,ig,gg(3),ix,iy,im, &
                    istat,gzero,ic,ip,ikk     ! counters
real(kind=dbl)             :: ktx,kty,ktrad,thickinc,ktmax,& ! maximum tangential component of wave vector
                    pre, pp, &   ! prefactors 
                    exer,&  ! excitation error
                    tpi,thick,g3n(3),g3(3),H,FNg(3),gp(3),RHOLZ,ll(3),lpg(3),gplen,LC3,sgdenom,cutoff, Znsq, kkl, &
                    DBWF, maxangle, braggangle  ! 
real,allocatable :: sr(:)
real(kind=dbl),allocatable :: Iz(:),Izsum(:,:,:),thickarray(:),percent(:,:) ! thickness array, results
integer(kind=irg),allocatable :: rlist(:)
complex(kind=dbl)  :: czero,cc
complex(kind=dbl),allocatable :: UghTable(:,:)


! display the standard program info
 progname = 'ECPz.f90'
 progdesc = 'Dynamical orientation-averaged depth profiles with HOLZ'
 call CTEMsoft

! set up section (*.xtal, atom positions, accelerating voltage and such)
 SG % SYM_reduce = .FALSE.
 call CrystalData
  
 call GetVoltage
! generate all atom positions
 call CalcPositions('v')

! zone axis computation
 srza = 'ZA'
! get k and f
 mess = ' Enter nearest high-symmetry zone axis direction'; call Message("(/A)")
 call GetIndex(k,'d')
 mess = ' Enter foil normal'; call Message("(/A)")
 call GetIndex(fn,'d')
 call TransSpace(float(fn),DynFN,'d','r')
 call NormVec(DynFN,'r')

write (*,*) 'Dynamical foil normal = ',DynFN

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do
 call BFsymmetry(k,j,isym,ir)

! determine the shortest reciprocal lattice points for this zone
 call ShortestG(k,ga,gb,isym)
 write (*,"(//,' ','[',3I2,'] has Bright Field symmetry ',A,'; order = ',I4,//)") k(1:3),PGTWD(isym),ir
 mess = 'Reciprocal lattice vectors : '; 
 oi_int(1:3)=ga(1:3)
 oi_int(4:6)=gb(1:3)
 call WriteInt(6,"('(',3I3,') and (',3I3,')',/)")

! enter range of incident beam directions, semi-number of pixels in image, and foil thickness
  mess = ' The program will use a symmetric range of incident beam directions,'; call Message("(/A)")
  mess = ' centered on the incident beam direction entered above.'; call Message("(A)")
  mess = ' Do you want to specify maximum k_t/|g| (1) or max angle (2) ? '; call GetInt(1)
  if (io_int(1).eq.1) then
    mess = ' Enter the maximum (positive) value ktmax of the tangential component of k'; call Message("(A)")
    mess = ' along the vector '; oi_int(1:3)=ga(1:3); call WriteInt(3,"(' (',3I3,')',$)")
    mess = ' in units of the length of this vector : '; call GetReal(1)
    ktmax = io_real(1)
  else
    mess = ' Enter the maximum beam tilt angle (in mrad) '; call Message("(A)")
    mess = ' along the vector '; oi_int(1:3)=ga(1:3); call WriteInt(3,"(' (',3I3,')')")
    mess = ' maximum angle : '; call GetReal(1)
    maxangle = io_real(1)
! next, convert this to a ktmax value by taking the ratio of maxangle and the Bragg angle for ga
    braggangle = 1000.0 * CalcDiffAngle(ga(1),ga(2),ga(3))/2.0
    ktmax = 0.5 * maxangle/braggangle
    mess = 'Program will use the following value for ktmax : '; oi_real(1) = ktmax; call WriteReal(1,"(F8.3)")
  end if
  mess = 'How many image pixels are needed for the beam tilt range [0,ktmax] ? '; call GetInt(1)
  npx = io_int(1)
  npy = npx
  mess = 'Max foil thickness and increment size [nm] '; call GetReal(2)
  thick = io_real(1)
  thickinc = dble(io_real(2))
  nthick = int(thick/thickinc)
  allocate(thickarray(nthick))
  do ih=1,nthick
    thickarray(ih) = dble(ih)*thickinc
  end do
  nat = 0

! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation
  call Calckvectors(k,ga,ktmax,npx,npy,numk)
  mess = 'Total number of independent incident beam directions inside cone = '; oi_int(1)=numk; call WriteInt(1,"(I8)")

! Then we must determine the masterlist of reflections (also a linked list);
! this list contains all reflections that are within some distance from the 
! Ewald sphere, starting with the ZOLZ upto, for now, the FOLZ (can be extended 
! to higher order if needed).  This list contains the zone-axis excitation errors
! and the extinction distances as well as the Miller indices.

! distance between consecutive HOLZ layers in nm-1
  H = 1.0/CalcLength(float(k),'d')

! determine g3 basis vector (vector normal to ZOLZ, with length H)
  call CalcCrossd(dble(ga),dble(gb),g3,'r','r',1)
  call NormVecd(g3,'r')
  g3n = g3
  g3 = H * g3
  FNg = (/ CalcDotd(dble(DynFN),dble(ga),'r'), CalcDotd(dble(DynFN),dble(gb),'r'), CalcDotd(dble(DynFN),g3,'r') /)
  
  write (*,*) 'FNg = ',FNg
    
  mess = 'basis vectors for this computation: '
  oi_real(1:3) = ga(1:3)
  oi_real(4:6) = gb(1:3)
  oi_real(7:9) = g3(1:3)  
  call WriteReal(9,"(/'ga = ',3f10.5,/'gb = ',3f10.5,/'g3 = ',3f10.5,/)")
  mess = 'reciprocal interplanar spacing H = '; oi_real(1) = H; call WriteReal(1,"(F10.4,' nm^-1'/)")

  call ShortestGFOLZ(k,ga,gb,gshort,gp,NN)
  mess = 'shortest vector to FOLZ = '; oi_int(1:3) = gshort(1:3); call WriteInt(3,"('(',3I3,')',/)")

! The master list is most easily created by brute force; we'll compute the 
! radius of the FOLZ ring, scale it by the length of ga or gb, turn that into an integer
! and make the range twice as large...  All reflections from FOLZ and ZOLZ within
! this range will become part of the list
  maxholz = 2
  minholz = 0
  i = maxholz
  numrefs = 1
  RHOLZ = sqrt(2.0*H*float(i)/mLambda - (float(i)*H)**2)
!  write (*,*) 'radius of HOLZ ring is ',RHOLZ
  im = max( nint(RHOLZ/CalcLengthd(dble(ga),'r')),nint(RHOLZ/CalcLengthd(dble(gb),'r')) )
!  write (*,*) 'range of reflections is ',im
  gzero = 1
  call AddReflection( (/0,0,0/) )   ! this guarantees that 000 is always the first reflection
  N = 0
!  do N=minholz,maxholz
    do ix=-im,im
      do iy=-im,im
        gg = ix*ga + iy*gb + N*gshort
        if ((IsGAllowed(gg)).AND.(CalcLengthd(dble(gg),'r').ne.0.D0)) then
           numrefs = numrefs + 1
           call AddReflection(gg)
        end if
      end do
    end do
!  end do
  mess = 'Length of the master list of reflections : '; oi_int(1) = numrefs; call WriteInt(1,"(I5,/)")
  rltmpa => reflist%next
  
! next, we start the major loop over all incident beam directions...
! for each value, we must compute how many reflections must be taken 
! into account (we'll keep track of min and max number), and then
! we need to initialize the dynamical matrix and solve the equations;
! this is to be repeated for all indicent beam directions.
  
! define the cutoff for the product of excitation error and extinction
! distance; only reflections for which this product is smaller than
! cutoff can contribute to the dynamical matrix.
  cutoff = 80.0
  


!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
! point to the first beam direction
  tmp => head
  czero = cmplx(0.0,0.0,dbl)
  pre = cmplx(0.0,cPi,dbl)
  call CalcUcg((/0,0,0/))   ! get the normal absorption parameter
  write (*,*) 'normal absorption length ',rlp%xgp
  DynUpz = rlp%Vpmod
  numset = cell % ATOM_ntype  ! number of special positions in the unit cell
  tpi = 2.D0*cPi
! allocate space for the results
  allocate(sr(nthick),rlist(numrefs),Iz(nthick),Izsum(2*npx+1,2*npy+1,nthick),percent(2*npx+1,2*npy+1))
  
! runtime analysis with Shark showed that the following loop spends nearly
! 75% of its time recomputing the potential Fourier coefficients.
! Therefore, we need to compute these ahead of time and simply use
! a look up table !!!!
!  write (*,*) 'Pre-computing all potential coefficients...',numrefs**2
  allocate(UghTable(numrefs,numrefs))
  cc = cmplx(-10000.0,0.0,dbl)
  UghTable = cc
!  rltmpa => reflist%next    ! reset the a list
! ir is the row index
!    do ir=1,numrefs
!     rltmpb => reflist%next   ! reset the b list
! ic is the column index
!     do ic=1,numrefs
!      if (ic.ne.ir) then  ! not a diagonal entry
! compute Fourier coefficient of electrostatic lattice potential 
!        call CalcUcg(rltmpa%hkl - rltmpb%hkl)
!        UghTable(ir,ic) = rlp%Ucg
!      end if
!      rltmpb => rltmpb%next
!     end do
!     rltmpa => rltmpa%next
!    end do
!   write (*,*) 'Done.  Starting main loop...'



!  work through the beam direction list
  beamloop: do ik=1,numk
  if (mod(ik,500).eq.0) write (*,*) ik,' out of',numk,' beam orientations'
! time the computation
   if (ik.eq.1) call Time_start

! get the Laue center (x and y components of the wave vector in units of ga and gb)
    ll = tmp%kt  !  (1) * float(ga) + tmp%kt(2) * float(gb)

! first, for this beam direction, determine the excitation errors of 
! all the reflections in the master list, and count the ones that are
! needed for the dynamical matrix (only if ki = 1)
    rltmpa => reflist%next
    nn = 0
    rlist = 0
    reflectionloop: do ig=1,numrefs
      gg = rltmpa%hkl
      lpg = ll+gg
      gplen = CalcLengthd(lpg,'r')
      LC3 = sqrt(1.D0-mLambda**2*CalcLengthd(ll,'r')**2)
      if (gplen.eq.0.D0) then
        exer = -mLambda*CalcDotd(dble(gg),2.D0*ll+gg,'r')/(2.D0*LC3*CalcDotd(g3n,dble(DynFN),'r'))  
      else
        sgdenom = 2.D0*LC3*CalcDotd(g3n,dble(DynFN),'r')-2.D0*mLambda*CalcDotd(lpg,dble(DynFN),'r')
        exer = -(mLambda*CalcDotd(dble(gg),2.D0*ll+gg,'r')-2.D0*LC3*CalcDotd(g3n,lpg,'r'))/sgdenom
      end if
      rltmpa%sg = exer  
! use the reflection num entry to indicate whether or not this
! reflection should be used for the dynamical matrix
        rltmpa%num = 0
        if (abs(exer) * rltmpa%xg .le. cutoff) then
          nn = nn+1
          rlist(nn) = ig
          rltmpa%num = nn
!          write (*,*) nn,rltmpa%hkl,rltmpa%sg,CalcDot(float(gg),(/1.0,1.0,1.0/),'r'),abs(exer) * rltmpa%xg
        end if
      rltmpa => rltmpa%next
    end do reflectionloop

    
! deallocate any previous Bloch wave arrays and reallocate them
! with the current number of beams
    if (allocated(DynMat)) deallocate(DynMat)
    if (allocated(W)) deallocate(W)
    if (allocated(CG)) deallocate(CG)
    if (allocated(alpha)) deallocate(alpha)

! initialize the dynamical matrix
    allocate(DynMat(nn,nn),stat=istat); DynMat = czero
    rltmpa => reflist%next    ! reset the a list
! ir is the row index
    do ir=1,nn
     if (rltmpa%num.eq.0) then 
      do 
       rltmpa => rltmpa%next    ! keep going until the first contributing reflection
       if (rltmpa%num.ne.0) exit
      end do
     end if
     rltmpb => reflist%next   ! reset the b list
! ic is the column index
     do ic=1,nn
      if (rltmpb%num.eq.0) then 
       do 
        rltmpb => rltmpb%next    ! keep going until the first contributing reflection
        if (rltmpb%num.ne.0) exit
       end do
      end if
      if (ic.ne.ir) then  ! not a diagonal entry
! compute Fourier coefficient of electrostatic lattice potential 
!          call CalcUcg(rltmpa%hkl - rltmpb%hkl)
!          DynMat(ir,ic) = rlp%Ucg
           if (UghTable(rlist(ir),rlist(ic)).eq.cc) then
             call CalcUcg(rltmpa%hkl - rltmpb%hkl)
             UghTable(rlist(ir),rlist(ic)) = rlp%Ucg
           end if
           DynMat(ir,ic) = UghTable(rlist(ir),rlist(ic))
      else  ! it is a diagonal entry, so we need the excitation error and the absorption length
        DynMat(ir,ir) = cmplx(2.D0*rltmpa%sg/mLambda,DynUpz,dbl)
      end if
      rltmpb => rltmpb%next
     end do
     rltmpa => rltmpa%next
    end do
! that should do it for the initialization of the dynamical matrix
DynMat(1,1) =  cmplx(0.0,DynUpz,dbl)

  
! solve the dynamical eigenvalue equation for this beam direction
  allocate(W(nn),CG(nn,nn),alpha(nn))
!  if (tmp%i.eq.50) then
!    if ((tmp%j.gt.-145).and.(tmp%j.lt.-140)) then
!      write (*,*) tmp%i,tmp%j
!      pp=1
!    end if
!  end if

  call CalcLghz(nn,thickarray,nthick,tmp%kn,gzero,Iz,pp)
! and add the resulting matrix to the existing one...

   ipx = tmp%i + npx + 1
   ipy = tmp%j + npy + 1
   Izsum(ipx,ipy,1:nthick) = Iz(1:nthick)
   percent(ipx,ipy) = pp
  
! select next beam direction
   tmp => tmp%next
  end do beamloop
    
! stop the clock and report the total time     
  call Time_stop(numk)
  
! save the Izsum array  
  call SafeOpenFile('d1','unformatted',fname)
  write (dataunit) 2*npx+1,2*npy+1,nthick
  write (dataunit) Izsum ! sr  
  call SafeCloseFile('d1','keep',fname)
! and store the 1% thickness array  
!  gname = '1ptest.data'
!  call SafeOpenFile('d1','unformatted',gname)
!  write (dataunit) 2*npx+1,2*npy+1,nthick
!  write (dataunit) percent 
!  call SafeCloseFile('d1','keep',gname)
end program


! ###################################################################
! 
!  subroutine Calckvectors
!
!  Author: Marc De Graef
!  
!  Description: computes the independent incident beam directions
!  for the multi-beam case and returns them as a linked list in the
!  global variables head and tail
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
! ###################################################################
subroutine Calckvectors(k,ga,ktmax,npx,npy,numk)

use local
use io
use error
use diffraction
use crystal
use dynamical

IMPLICIT NONE

integer(kind=irg)              :: npx,npy,numk,k(3),ga(3),istat,imin,imax,jmin,jmax,i,j
real(kind=dbl)                 :: ktmax,kr(3),glen,delta,kstar(3),kt(3),gan(3),gperp(3),ktlen

intent(IN)           :: npx,npy,k,ga,ktmax
intent(OUT)          :: numk

! compute geometrical factors 
 glen = CalcLengthd(dble(ga),'r')              ! length of ga
 gan = ga/glen                                 ! normalized ga
 delta = 2.D0*ktmax*glen/(2.D0*dble(npx)+1.D0)   ! grid step size in nm-1 
 call TransSpaced(dble(k),kstar,'d','r')       ! transform incident direction to reciprocal space
 call CalcCrossd(dble(ga),kstar,gperp,'r','r',0)      ! compute g_perp = ga x k
 call NormVecd(gperp,'r')                       ! normalize g_perp
 call NormVecd(kstar,'r')                       ! normalize reciprocal beam vector
write (*,*) 'Calckvectors : ',k,ga,ktmax,npx,npy,delta
write (*,*) 'Calckvectors : ',kstar,gperp

open(unit=20,file='ECPkt.txt',status='unknown',form='formatted')

 if (.not.associated(head)) then     ! allocate the head and tail of the linked list
   allocate(head,stat=istat)         ! allocate new value
   if (istat.ne.0) call FatalError('Calckvectors: unable to allocate head pointer',' ')
   tail => head                      ! tail points to new value
   nullify(tail%next)                ! nullify next in new value
   numk = 1                          ! keep track of number of k-vectors so far
   tail%i = 0                             ! i-index of beam
   tail%j = 0                             ! j-index of beam
   tail%kt = (/0.0,0.0,0.0/)              ! tangential component of k
   call NormVecd(kstar,'r')                ! normalize incident direction
   tail%k = kstar/mLambda                 ! divide by wavelength and assign
   tail%kn = CalcDotd(tail%k,kstar,'r')    ! normal component of k
     write(unit=20,*) tail%i,tail%j,tail%kt,tail%k,tail%kn,CalcLengthd(tail%k,'r')
 else
   call FatalError('Calckvectors: pointer head already allocated',' ')
 end if
! the following lines are quite different if symmetry is taken into account;
! check the MBsym.f90 program to determine how that can be done.
  imin =  -npx; imax = npx; jmin = -npy; jmax = npy; 
! now do the real work
  do i=imin,imax
   do j=jmin,jmax
    if (.not.((i.eq.0).and.(j.eq.0))) then  ! the point (0,0) has already been taken care of
     allocate(tail%next,stat=istat)  ! allocate new value
     if (istat.ne.0) call FatalError('Calckvectors: unable to allocate pointer',' ')
     tail => tail%next               ! tail points to new value
     nullify(tail%next)              ! nullify next in new value
     numk = numk + 1                 ! keep track of number of k-vectors so far
     tail%i = i                      ! i-index of beam
     tail%j = j                      ! j-index of beam
     kt = -dble(i)*delta*gan - dble(j)*delta*gperp  ! tangential component of k
     tail%kt = kt                    ! store tangential component of k
     ktlen = delta**2*(i*i+j*j)      ! squared length of tangential component
     kr = kt + sqrt(1.D0/mLambda**2 - ktlen)*kstar ! complete wave vector
     tail%k = kr                     ! store in pointer list
     tail%kn = CalcDotd(tail%k,kstar,'r')    ! normal component of k
       write(unit=20,*) tail%i,tail%j,tail%kt,tail%k,tail%kn,CalcLengthd(tail%k,'r')
!     write (*,*) i,j,tail%kt,tail%kn
    end if
   end do
  end do
  close(unit=20,status='keep')
  
end subroutine



! ###################################################################
! 
!  subroutine CalcLghz
! 
!  Author: Marc De Graef
!  
!  Description: Bloch wave portion for the orientation-averaged depth
!  profile.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   8/24/20 MDG 1.0 original
! ###################################################################
subroutine CalcLghz(nn,thickarray,nthick,kn,gzero,Iz,pp)

use local
use io
use files
use diffraction
use dynamical
use constants

IMPLICIT NONE

integer         :: nn,i,j,k,kk,ig,ih,IPIV(nn),gzero,nthick
complex(kind=dbl),allocatable :: CGinv(:,:), Minp(:,:),diag(:),tmp3(:,:)
real(kind=dbl) :: thickarray(nthick),Iz(nthick),s,t,pp
real(kind=dbl) :: q,kn


intent(IN)      :: nn,thickarray,nthick,kn,gzero

allocate(CGinv(nn,nn),Minp(nn,nn),diag(nn),tmp3(nn,nn))

! get the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 W = cmplx(0.D0,0.D0)
 CG = cmplx(0.D0,0.D0)
 CGinv = cmplx(0.D0,0.D0)
 call BWsolve(Minp,W,CG,CGinv,nn,IPIV,.FALSE.)    ! .FALSE. employs the Schur vector route, omitting it uses direct eigenvectors
 
! scale the eigenvalues by the wave number
 W = W/cmplx(2.D0*kn,0.D0,dbl)
! recall that alpha(1:nn) = CGinv(1:nn,gzero)

! make sure the alpha excitation coefficients are normalized 
s = sum(cabs(CGinv(1:nn,gzero))**2)
if (s.gt.1.D0) then
  s = cmplx(1.D0/sqrt(s),0.0,dbl)
  CGinv(1:nn,gzero) = CGinv(1:nn,gzero)*s
endif 

! compute the thickness array 
Iz = 0.D0
 do j=1,nn
    q = -4.D0*cPi*aimag(W(j))
    s = cabs(CGinv(j,gzero))**2
    do k=1,nthick
      t = q*thickarray(k)
      if (abs(t).lt.30.D0) Iz(k) = Iz(k) +  s * exp(t)
    end do
 end do

 
! and finally compute at which thickness only 1% of the total Bloch wave intensity remains
! (this is used in the ZAECCI program as an estimate of where to begin the integrations)

! first find the rough position from the Iz array
k = 1
do 
  if (Iz(k).lt.0.01) exit 
  k=k+1
end do

! then use linear interpolation
pp = thickarray(k-1) + (Iz(k-1)-0.01D0)*(thickarray(k)-thickarray(k-1))/(Iz(k-1)-Iz(k))

deallocate(CGinv,Minp,diag,tmp3)
end subroutine



! ###################################################################
!
!  subroutine ShortestGFOLZ
!
!  Author: Marc De Graef
!
!  Description: find the G vector (displacement FOLZ w.r.t. ZOLZ, Chapter 3)
!
!  History
!
!  modified by  rev reason
!  -------
!  01/29/02 MDG 1.0 original
! ###################################################################
subroutine ShortestGFOLZ(k,ga,gb,gshort,gp)

use local
use io
use crystal
use crystalvars
use error

IMPLICIT NONE

real(kind=sgl)               :: gmin,gam11,gam12,gam22,gp(3),g1(3),g2(3),g3(3),glen
integer(kind=irg),parameter  :: inm = 8
integer(kind=irg)            :: ih,ik,il,NN,ga(3),gb(3),gshort(3),k(3)

! look for the shortest reflection satisfying hu+kv+lw = 1
! This could be replaced by code from Jackson's paper (1987),
! but it does essentially the same thing.
 gmin = 100.0
 NN=1
 g1 = float(ga)
 g2 = float(gb)
 do while((gmin.eq.100.0).and.(NN.lt.4))
  do ih=-inm,inm
   do ik=-inm,inm
    do il=-inm,inm
! does this reflection lie in the plane NN ?
     if ((ih*k(1)+ik*k(2)+il*k(3)).eq.NN) then
      glen = CalcLength(float((/ih,ik,il/)),'r')
      if (glen.lt.gmin) then
       gmin = glen
       gshort =  (/ ih,ik,il /) 
      end if
     end if
    end do
   end do
  end do
  if (gmin.eq.100.0) then 
    mess = 'Could not find any reflections with hu+kv+lw = '; oi_int(1)=NN; call WriteInt(1,"(I2)")
    NN = NN+1
  end if
 end do
 if (gmin.eq.100.0) then ! for some reason there is no reflection with N<=3 ...
  call FatalError('HOLZ: could not find any reflections with hu+kv+lw<=3 ...',' ')
 end if
 g3 = float(gshort)
! projected components of G
 gam11 = CalcDot(g1,g1,'r')
 gam12 = CalcDot(g1,g2,'r')
 gam22 = CalcDot(g2,g2,'r')
 gmin = 1.0/(gam11*gam22-gam12**2)
 gp(1) = (CalcDot(g3,g1,'r')*gam22-CalcDot(g3,g2,'r')*gam12)*gmin
 gp(2) = (CalcDot(g3,g2,'r')*gam11-CalcDot(g3,g1,'r')*gam12)*gmin
end subroutine ShortestGFOLZ


