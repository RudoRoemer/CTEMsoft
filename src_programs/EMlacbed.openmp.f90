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
! EMsoft:EMlacbed.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMlacbed 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Zone axis LACBED
!
!> @todo implement full symmetry use; implement OpenMP; implement multiple reflection output
!> or, easier perhaps, selection of one reflection; a nice addition would be the ability to define a 
!> foil thickness profile that can be read from an external file. For each beam direction, the 
!> thickness would be slightly different, so that one can obtain more realistic LACBED patterns.
!> Usually, the beam illuminates a large area of the sample, sothe thickness can vary quite a 
!> bit; Joachim Mayer's [111] Si pattern at 100 kV is a nice example of that (page v of Spence&Zuo)
!
!> @date  11/29/01 MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!> @date 05/08/13 MDG 2.1 forked from mbcbed and adapted for large angle CBED patterns
!--------------------------------------------------------------------------
program EMlacbed

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

IMPLICIT NONE

real(kind=sgl)			:: io_real(1)

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData
 call GetVoltage

 call ReadValue(' Camera length L  [mm, real] ', io_real, 1)
 camlen = io_real(1)

! generate all atom positions
 call CalcPositions('v')

! generate a set of zone axis CBED patterns
 call LACBEDpattern

end program

!--------------------------------------------------------------------------
!
! SUBROUTINE:LACBEDpattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a large angle zone axis convergent beam electron diffraction pattern
!
! 
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/13/13  MDG 3.0 OpenMP implementation added 
!--------------------------------------------------------------------------
subroutine LACBEDpattern

use local
use constants
use crystal
use crystalvars
use diffraction
use gvectors
use kvectors
use postscript, ONLY: GetIndex
use symmetry
use math
use dynamical
use io
use error
use files
use omp_lib

IMPLICIT NONE

real(kind=sgl)      			:: laL,kt,z0,thc,thb,hkl(3),ind(3),fn(3),ktmax, io_real(3), &
                       				   dom,glen,qr,qi,bg,xgpz,bragg,c(3),RR,gx(3),gy(3),gg(3), thetac, thick, &
                       				   sc, scmax, qx, qy, frac, dmin
integer(kind=irg)   			:: g(3),ira,dpcnt,ppi,ijmax,ga(3),gb(3),k(3),fcnt,ccnt,cnt, PX, PY, &
                       				  newcount,count_rate,count_max, io_int(6), i, j, isym, ir, order, nn, &
                       				  iorder, npx, npy, numt, im, numk, npix, ik, ip, jp, maxholz, numpix, istat, nthreads, TID, &
                       				  NUMTHREADS,DynN
character(1)        			:: ans
character(2)        			:: srza
character(3)				:: method
character(fnlen)     			:: fname, tname


complex(kind=sgl)   		:: czero
logical             				:: overlap,np,first, ConvertList
real(kind=sgl),parameter     	:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
real(kind=sgl),allocatable    	:: disk(:,:,:), thickarray(:,:)
real(kind=dbl),allocatable	:: ktmpk(:,:), ktmpkn(:)
integer(kind=irg),allocatable 	:: diskoffset(:,:), ktmpi(:), ktmpj(:)
real(kind=sgl),allocatable    	:: inten(:,:)

!$OMP THREADPRIVATE(DynMat,BPstronglist,BPweaklist,BPstronghkl,BPweakhkl,BPstrongsg,BPweaksg, &
!$OMP& BPreflistindex,BPnns,BPnnw)

! zone axis computation
 srza = 'ZA'

! get k and f
 mess = ' Enter wave vector direction'; call Message("(/A)")
 call GetIndex(k,'d')
 mess = ' Enter foil normal'; call Message("(/A)")
 call GetIndex(g,'d')
 DynFN = float(g)
 
! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do
! and the Bright Field symmetry
 call BFsymmetry(k,j,isym,ir)
 io_int(1:3) = k(1:3)
 call WriteValue('', io_int, 3, "(//,' ','[',3I2,'] has Bright Field symmetry ',$)")
 mess = PGTWD(isym)
 call Message("(A,$)")
 io_int(1) = ir 
 call WriteValue('; order = ',io_int, 1, "(I4,//)")


! determine the shortest reciprocal lattice points for this zone
 call ShortestG(k,ga,gb,isym)
 io_int(1:3)=ga(1:3)
 io_int(4:6)=gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! construct the list of all possible reflections and convert it to regular allocatable arrays
dmin = 0.025
method = 'ALL'
maxholz = 0
ConvertList = .TRUE.
call Compute_ReflectionList(dmin,k,ga,gb,method,ConvertList,maxholz)

! enter range of incident beam directions
  mess = ' The program will use a symmetric cone of beam directions,'; call Message("(/A)")
  mess = ' centered on the incident beam direction entered above.'; call Message("(A)")
  bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5*1000.0
  io_real(1) = bragg
  call WriteValue(' The Bragg angle for the first reflection is equal to [mrad]: ',io_real, 1, "(F8.5)")

  call ReadValue(' Enter the beam convergence angle theta_c in mrad: ', io_real, 1)
  thetac = io_real(1)/1000.0
  
! convert to ktmax along ga
  ktmax = 0.5*thetac*1000.0/bragg

! compute number of pixels along diameter of central disk for given camera length
  RR = 300.0/25.4   ! dots per millimeter for 300 dots per inch
  npx = nint(RR*camlen*thetac)
  npy = npx
  io_int(1) = 2.0*npx
  call WriteValue(' Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")

  ijmax = float(npx)**2   ! truncation value for beam directions
! get number of thicknesses for which to compute the CBED pattern
  mess = ' This program computes LACBED patterns either for a series of thicknesses'; call Message("(A)")
  mess = ' or for a single variable thickness array; enter a negative thickness for second option'; call Message("(A)")
  call ReadValue(' Enter the first thickness [nm, R]', io_real, 1)
  if (io_real(1).gt.0.0) then
    thick = io_real(1)
    call ReadValue(' How many multiples of this thickness [I] ', io_int, 1) 
    numt = io_int(1)
  else
     call ReadValue(' Enter filename for thickness array : ', tname, '(A)')
     open(UNIT=dataunit,FILE=trim(tname),STATUS='old',FORM='unformatted')
     allocate(thickarray(900,900),stat=istat)
     read(UNIT=dataunit) numpix
     if (numpix.ne.900) call FatalError('LACBED ',' thickness array dimensions do not match required size')
     read(UNIT=dataunit) thickarray
     close(UNIT=dataunit,STATUS='keep')
     numt = 1
     thick = 0.0
  end if  
! determine all independent incident beam directions (use a linked list starting at khead)
! we're using the StandardConical method since we do want to make use of symmetry
! operations, but the illumination must be limited to a cone.  
!  call Calckvectors(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,'StandardConical')
  call Calckvectors(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,'Conical')
 

! allocate the disk variable which will hold the entire computed pattern
  npix = int(3.0*300.0)   ! 3 inches wide at 300 dpi
  allocate(disk(numt,npix,npix))
  disk=0.0

  sc = cell%mLambda*camlen*RR
  scmax = 1.5*300.0 + npx
  PX = npix/2


! allocate the offset array
! to get this array, we need to do a mock initialization of the dynamical matrix in zone axis orientation
!  call Compute_DynMat('BLOCHBETHE', khead%k, .TRUE.)
  frac = 0.05

  io_int(1)=numk
  call WriteValue(' Starting computation for # beam directions = ', io_int, 1, "(I8)")

! if we want to use OpenMP, then we need to first convert the wave vector information
! into regular allocatable arrays that can then become SHARED access to all threads.
!
! we need the following items from the khead linked pointer list:
! i, j, k(3)
allocate(ktmpi(numk), ktmpj(numk), ktmpkn(numk), ktmpk(3,numk),stat=istat)
! point to the first beam direction
  ktmp => khead
! loop over list and then deallocate the linked list
  do ik=1,numk
    ktmpi(ik) = ktmp%i
    ktmpj(ik) = ktmp%j
    ktmpk(1:3,ik) = ktmp%k
    ktmpkn(ik) = ktmp%kn
    ktmp => ktmp%next
  end do
  call Delete_kvectorlist()

! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)

if (BetheParameter%weakcutoff.eq.0.0) call Set_Bethe_Parameters

nthreads = 2

call OMP_SET_NUM_THREADS(nthreads)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,ip,jp,DynN,inten,j,ik )

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()

if (TID.eq.0) write (*,*) 'Number of allocated threads = ',NUMTHREADS

!$OMP DO SCHEDULE(DYNAMIC,32)
! loop over all beam orientations
  do ik = 1,numk

! is this point inside the viewing square ?
      ip = PX  - ktmpi(ik)
      jp = PX + ktmpj(ik)

      if (((ip.ge.1).and.(ip.le.npix)).and.((jp.ge.1).and.(jp.le.npix))) then

  ! compute the dynamical matrix using Bloch waves with Bethe potentials 
        call LocalCompute_DynMat('BLOCHBETHE', ktmpk(1:3,ik), .TRUE.,DynN,TID)

! allocate the intensity array
        allocate(inten(numt,DynN))
        inten = 0.0

! solve the dynamical eigenvalue equation
       if (allocated(thickarray)) then 
         call CalcBWint(DynN, ktmpkn(ik), numt, thickarray(ip,jp), inten)
      else
        call CalcBWint(DynN, ktmpkn(ik),numt,thick,inten)
      end if

! copy in the correct locations
      do j=1,numt
         disk(j,ip,jp) = disk(j,ip,jp) + inten(j,1)
      end do

! and remove the intensity array
     deallocate(inten)
    
   end if

! update computation progress [will need to be rewritten for OpenMP version]
   if ((TID.eq.0).and.(float(ik)/float(numk) .gt. frac)) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

  end do
!$OMP END DO

!$OMP BARRIER
!$OMP END PARALLEL


write (*,*) ' maximum intensity ', maxval(disk)

! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  call WriteValue(' Total computation time [s] ' , io_real, 1, "(F)")

 loadingfile = .FALSE.
  call SafeOpenFile('d1','unformatted',fname)
  write (dataunit) numt,npix
  write (dataunit) disk
  call SafeCloseFile('d1','keep',fname)

end subroutine LACBEDpattern


!--------------------------------------------------------------------------
!
! SUBROUTINE: LocalCompute_DynMat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the entire dynamical matrix, including HOLZ and Bethe potentials
!
!> @details This is a complicated routine because it forms the core of the dynamical
!> scattering simulations.  The routine must be capable of setting up the dynamical 
!> matrix for systematic row and zone axis, with or without HOLZ reflections, and 
!> with Bethe potentials for the Bloch wave case.  The routine must also be able to
!> decide which reflections are weak, and which are strong (again for the Bloch wave
!> case, but potentially also for other cases? Further research needed...).
!
!> @param calcmode string that describes the particular matrix mode
!
!> @date   05/06/13 MDG 1.0 original
!--------------------------------------------------------------------------
recursive subroutine LocalCompute_DynMat(calcmode,kk,IgnoreFoilNormal,DynN,TID)

use local
use dynamical
use error
use constants
use crystal
use diffraction
use gvectors
use io

IMPLICIT NONE

character(*),INTENT(IN)			:: calcmode		!< computation mode
real(kind=dbl),INTENT(IN)		:: kk(3)			!< incident wave vector
logical,INTENT(IN)				:: IgnoreFoilNormal	!< how to deal with the foil normal
integer(kind=irg),INTENT(OUT)	:: DynN		!< number of beams for OpenMP implementation
integer(kind=irg),INTENT(IN)		:: TID

complex(kind=dbl)  			:: czero,pre, weaksum, ughp, uhph
integer(kind=irg) 		 		:: istat,ir,ic,nn, iweak, istrong, iw, ig, ll(3), gh(3), nnn
real(kind=sgl)     				:: glen,exer,gg(3), kpg(3), gplen, sgp
 
! has the list of reflections been allocated ?
if (.not.associated(reflist)) call FatalError('Compute_DynMat',' reflection list has not been allocated')

! if the dynamical matrix has already been allocated, deallocate it first
! this is partially so that no program will allocate DynMat itself; it must be done
! via this routine only.
if (allocated(DynMat)) deallocate(DynMat)

! initialize some parameters
czero = cmplx(0.0,0.0,dbl)	! complex zero
pre = cmplx(0.0,cPi,dbl)		! i times pi


! we don't know yet how many strong reflections there are so we'll need to determine this first
! this number depends on some externally supplied parameters, which we will get from a namelist
! file (which should be read only once by the Set_Bethe_Parameters routine), or from default values
! if there is no namelist file in the folder.
!	if (BetheParameter%weakcutoff.eq.0.0) call Set_Bethe_Parameters


! reset the value of DynNbeams in case it was modified in a previous call 
  	DynNbeams = DynNbeamsLinked
  
! first, for the input beam direction, determine the excitation errors of 
! all the reflections in the master list, and count the ones that are
! needed for the dynamical matrix (weak as well as strong)
        if (.not.allocated(BPweaklist)) allocate(BPweaklist(DynNbeams))
        if (.not.allocated(BPstronglist)) allocate(BPstronglist(DynNbeams))
        if (.not.allocated(BPreflistindex)) allocate(BPreflistindex(DynNbeams))

	BPweaklist = 0
	BPstronglist = 0
	BPreflistindex = 0

! deal with the transmitted beam first
    nn = 1
    nnn = 1
    BPstronglist(nn) = 1   ! make sure that the origin is always a strong reflection...
    BPweaklist(nn) = 0
    BPreflistindex(nn) = 1

    Reflist_sg(nn) = 0.D0    

reflectionloop: do ig=2,DynNbeamsLinked
      gg(1:3) = Reflist_hkl(1:3,ig)      

! deal with the foil normal; if IgnoreFoilNormal is .TRUE., then assume it is parallel to the beam direction
     if (IgnoreFoilNormal) then 
! we're taking the foil normal to be parallel to the incident beam direction at each point of
! the standard stereographic triangle, so cos(alpha) = 1 in eqn. 5.11 of EM
        kpg = kk+gg                		! k0 + g (vectors)
        gplen = CalcLength(kpg,'r')  	! |k0+g|
        Reflist_sg(ig) = (1.0/cell%mLambda**2 - gplen**2)*0.5/gplen
     else
	Reflist_sg(ig) = Calcsg(gg,sngl(kk),DynFN)
     end if

! use the reflection num entry to indicate whether or not this
! reflection should be used for the dynamical matrix
! there are two criteria: one is that |sg|*xi must be smaller than some cutoff.
! the other is that the excitation error itself must be smaller than a cutoff in
! order for the reflection to be a strong reflection.  So, if a reflection is considered
! weak according to the first criterion, it could still be made strong if the 
! excitation error is really small (to avoid a blowup of the Bethe perturbation algorithm).
        sgp = abs(Reflist_sg(ig)) * Reflist_xg(ig)
        if (sgp.le.BetheParameter%cutoff) then
          nn = nn+1
! is this a weak or a strong reflection (in terms of Bethe potentials)? 
          if ((sgp.le.BetheParameter%weakcutoff).or.(abs(Reflist_sg(ig)).le.BetheParameter%sgcutoff)) then 
            nnn = nnn+1
	    BPstronglist(ig) = 1
            BPreflistindex(ig) = nnn
          else
            BPweaklist(ig) = 1
          end if
!          write (*,*) TID,nn,nnn
        end if
    end do reflectionloop
    
! if we don't have any beams in this list (unlikely, but possible if the cutoff and 
! weakcutoff parameters have unreasonable values) then we skip the rest of this
! incident beam computation (and we report it to the user) 
	 if (nn.eq.0) then
	   mess = ' no beams found for the following parameters:'; call Message("(A)")
	   write (*,*) ' wave vector = ',kk,'  -> number of beams = ',nn
	   mess =  '   -> check cutoff and weakcutoff parameters for reasonableness'; call Message("(A)")
	  call FatalError('Compute_DynMat','No beams in list')
	end if

! next, we define nns to be the number of strong beams, and nnw the number 
! of weak beams.
	 BPnns = sum(BPstronglist)
	 BPnnw = sum(BPweaklist)

! We may want to keep track of the total and average numbers of strong and weak beams  
!	 BetheParameter%totweak = BetheParameter%totweak + BetheParameter%nnw
!	 BetheParameter%totstrong = BetheParameter%totstrong + BetheParameter%nns
!	 if (BetheParameter%nnw.lt.BetheParameter%minweak) BetheParameter%minweak=BetheParameter%nnw
!	 if (BetheParameter%nnw.gt.BetheParameter%maxweak) BetheParameter%maxweak=BetheParameter%nnw
!	 if (BetheParameter%nns.lt.BetheParameter%minstrong) BetheParameter%minstrong=BetheParameter%nns
!	 if (BetheParameter%nns.gt.BetheParameter%maxstrong) BetheParameter%maxstrong=BetheParameter%nns

! allocate arrays for weak and strong beam information
	if (allocated(BPweakhkl)) deallocate(BPweakhkl)
	if (allocated(BPweaksg)) deallocate(BPweaksg)
	if (allocated(BPstronghkl)) deallocate(BPstronghkl)
	if (allocated(BPstrongsg)) deallocate(BPstrongsg)
	allocate(BPweakhkl(3,BPnnw),BPweaksg(BPnnw))
	allocate(BPstronghkl(3,BPnns),BPstrongsg(BPnns))

! here's where we extract the relevant information from the linked list (much faster
! than traversing the list each time...)
	iweak = 0
	istrong = 0
	do ir=1,DynNbeamsLinked
	     if (BPweaklist(ir).eq.1) then
	        iweak = iweak+1
	        BPweakhkl(1:3,iweak) = Reflist_hkl(1:3,ir)
	        BPweaksg(iweak) = Reflist_sg(ir)
	     end if
	     if (BPstronglist(ir).eq.1) then
	        istrong = istrong+1
	        BPstronghkl(1:3,istrong) = Reflist_hkl(1:3,ir)
	        BPstrongsg(istrong) = Reflist_sg(ir)
	     end if
	end do


! now we are ready to create the dynamical matrix
	DynN = BPnns


! allocate DynMat if it hasn't already been allocated
	  if (allocated(DynMat)) deallocate(DynMat)
	  allocate(DynMat(DynN,DynN),stat=istat)
	  DynMat = czero

! get the absorption coefficient
	  call CalcUcg((/0,0,0/))
	  DynUpz = rlp%Vpmod

! ir is the row index
       do ir=1,BPnns
! ic is the column index
          do ic=1,BPnns
! compute the Bethe Fourier coefficient of the electrostatic lattice potential 
              if (ic.ne.ir) then  ! not a diagonal entry
                ll = BPstronghkl(1:3,ir) - BPstronghkl(1:3,ic)
                DynMat(ir,ic) = LUT(ll(1),ll(2),ll(3)) 
! and subtract from this the total contribution of the weak beams
                weaksum = czero
                do iw=1,BPnnw
                      ll = BPstronghkl(1:3,ir) - BPweakhkl(1:3,iw)
                      ughp = LUT(ll(1),ll(2),ll(3)) 
                      ll = BPweakhkl(1:3,iw) - BPstronghkl(1:3,ic)
                      uhph = LUT(ll(1),ll(2),ll(3)) 
                      weaksum = weaksum +  ughp * uhph *cmplx(1.D0/BPweaksg(iw),0.0,dbl)
                 end do
! and correct the dynamical matrix element to become a Bethe potential coefficient
                 DynMat(ir,ic) = DynMat(ir,ic) - cmplx(0.5D0*cell%mLambda,0.0D0,dbl)*weaksum
               else  ! it is a diagonal entry, so we need the excitation error and the absorption length
                 DynMat(ir,ir) = cmplx(2.D0*BPstrongsg(ir)/cell%mLambda,DynUpz,dbl)
               end if
          end do
        end do
! that should do it for the initialization of the dynamical matrix, except for the first element
       DynMat(1,1) =  cmplx(0.0,DynUpz,dbl)

end subroutine LocalCompute_DynMat


! ###################################################################
! 
!  subroutine CalcBW
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
recursive subroutine CalcBWint(nn,kn,nt,thick,inten)

use local
use io
use diffraction
use kvectors
use dynamical
use constants

IMPLICIT NONE

integer(kind=irg)                             	:: nn,i,j,ig,ih,IPIV(nn),nt
complex(kind=dbl)    				:: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), lalpha(nn)
real(kind=sgl)                             		:: thick,inten(nt,nn),th, kn
! complex(kind=dbl)                			:: Ijk(nn,nn),Lgh(nn,nn),q

intent(IN)      :: nt, nn, thick, kn
intent(OUT)     :: inten

! should we just return zero ?
if (nt.eq.1) then
  if (thick.eq.0.0) then 
    inten = 0.0
    return
  end if
end if

! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
 lW = cPi*Wloc/cmplx(kn,0.0)
 do i=1,nn
  lalpha(i) = CGinv(i,1)
 end do 
 
 do i=1,nt
  th = dble(thick*float(i))
  diag(1:nn)=exp(-th*imag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
  do j=1,nn
   inten(i,j) = abs(sum(lCG(j,1:nn)*diag(1:nn)))**2
  end do 
 end do
 
end subroutine CalcBWint
