
program tester

use local
use constants
use initializers
use crystalvars
use crystal
use files
use gvectors
use dynamical
use diffraction

IMPLICIT NONE

real(kind=sgl)		        :: voltage, dmin, k(3)
character(fnlen)		:: xtalname
type(unitcell)	                :: cellA
logical	                        :: verbose
integer(kind=irg)              :: i,ii, istat, numt, j1, j2, c1min, c1max, c2min, c2max, c1step, c2step, &
                                  nref, ir
type(reflisttype),pointer	:: rlr, rlw
complex(kind=dbl)		:: czero, pre
real(kind=sgl)			:: kn
real(kind=dbl),allocatable 	:: thick(:), inten2(:,:), intenref(:,:)
integer(kind=sgl),allocatable	:: hkllist(:,:), rlist(:)


verbose = .FALSE.
numt = 250
c1min = 4
c1max = 40
c1step = 2
c2min = 30
c2max = 60
c2step = 2

allocate(thick(numt), stat=istat)
thick = dble( (/ (i-1,i=1,numt) /) ) * 0.5D0

! initialize some parameters
czero = cmplx(0.0,0.0,dbl)	! complex zero
pre = cmplx(0.0,cPi,dbl)		! i times pi

xtalname = 'Cu.xtal'
dmin = 0.01
voltage = 200000.0
call Initialize_Cell(xtalname, dmin, voltage, verbose)
call CopyFromCell(cellA)

open(unit=dataunit,file='testrun.data',status='unknown',form='unformatted')

do ii=0,0
  DynFN = (/ sin(float(ii)*0.02), 0.0, cos(float(ii)*0.02) /) 
  k = DynFN
  call NormVec(k,'r')
  k = k/sngl(mLambda)
  write (*,*) 'wave vector = ',k
  
! this needs to be modified
  kn = k(3)

! and generate the reflectionlist
  call Initialize_ReflectionList(k,dmin,verbose)

        write (*,*) 'Total number of reflections in list : ', cell%DynNbeams
        
        
! first, we need to do the full computation, using a whole lot of reflections
! and we'll take that as the "correct" result.  Then we compare similar computations
! with fewer beams and various Bethe cutoff ranges to figure out the range of 
! parameters for which the result remains accurate to 10^(-5)
  BetheParameter%c2 = 20000.0
  BetheParameter%c1 = 20000.0 
  call Apply_BethePotentials()

! generate the dynamical matrix
  call GetDynMat()
  write (*,*) 'reference dynmat = ',shape(DynMat)

! next, we do a Bloch wave computation and generate the resulting scattered intensities
! for a series of 250 thickness values (increments of 1 nm)
  if (allocated(intenref)) deallocate(intenref)
  allocate(intenref(numt,cell%nns))
  intenref = 0.0
 
! solve the dynamical eigenvalue equation and return the intensities of ALL reflections,
! both strong and weak; the weak intensities should also be plotted at the correct locations....
! this uses a new version of the CalcBWint routine that implements both strong and weak beam intensities.
  call CalcBWint(cell%nns,cell%nnw,numt,thick,kn,intenref)

! store the data to a file
  write (dataunit) cell%nns, numt
  allocate(hkllist(3,cell%nns),stat=istat)
  rlr => cell%reflist%next
  ir = 1
! go through this list and assign each a sequential number in the famnum field
  do
    if (.not.associated(rlr)) EXIT
    hkllist(1:3,ir) = rlr%hkl(1:3)
    rlr%famnum = ir
    rlr => rlr%nexts
    ir = ir + 1
  end do
  write (dataunit) hkllist
  write (dataunit) thick
  write (dataunit) intenref
  write (dataunit) c1min, c1max, c1step, c2min, c2max, c2step
  deallocate(intenref)

! ok, so we have the reference intensity, computed with lots of beams and without 
! the Bethe approximation.  Next, we need to scan over a range of Bethe truncation
! parameters and compare the intensities.  One issue will be that the reflections
! are in a different order each time, so we'll need a translation routine to convert
! the inten array into the same order as the intenref array...

  do j1=c1min,c1max,c1step
    BetheParameter%c1 = float(j1)  
    do j2=c2min,c2max,c2step
      BetheParameter%c2 = float(j2)  
      if (j2.lt.j1) then
! write negative numbers so that the IDL routine "knows" that this combination was omitted
        write (dataunit) -j1, -j2, -j1, -j2
        EXIT
      end if
      
! nullify all the nexts and nextw pointers
      call ResetPointers()

! create the new lists      
      call Apply_BethePotentials()
      if (allocated(rlist)) deallocate(rlist)
      allocate(rlist(cell%nns+cell%nnw),stat=istat)

! dynamical simulation
      call GetDynMat()
      write (*,*) 'ended GetDynMat ', float(j1), float(j2), shape(DynMat), cell%nns, cell%nnw
      if (allocated(inten2)) deallocate(inten2)
      allocate(inten2(numt,cell%nns+cell%nnw))
      inten2 = 0.0

      call CalcBWint(cell%nns,cell%nnw,numt,thick,kn,inten2)

! make the reference list in the correct order
      rlist = 0
! first the strong reflections
      rlr => cell%reflist%next
      ir = 1
      do
        if (.not.associated(rlr)) EXIT
        rlist(ir) = rlr%famnum
        rlr => rlr%nexts
        ir = ir + 1
      end do
! then the weak reflections      
      rlw => cell%firstw
      do
        if (.not.associated(rlw)) EXIT
        rlist(ir) = rlw%famnum
        rlw => rlw%nextw
        ir = ir + 1
      end do

! and write the result to the data file
      write (dataunit) j1, j2, cell%nns, cell%nnw
      write (dataunit) rlist
      write (dataunit) inten2
    end do
    write (*,*) 'completed c1 = ',j1
  end do

end do

close(unit=dataunit,status='keep')


end program tester








subroutine ResetPointers()

use local
use crystalvars

IMPLICIT NONE

type(reflisttype),pointer	:: rlr

rlr => cell%reflist%next
do
  if (.not.associated(rlr)) EXIT
  nullify(rlr%nexts)
  nullify(rlr%nextw)
  rlr => rlr%next
end do

end subroutine ResetPointers





!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcBWint
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the scattered intensities for a range of thicknesses
!
!> @param nn number of strong beams
!> @param nw number of weak beams
!> @param nt number of thickness values
!> @param thick thickness array
!> @param kn normal component of wave vector
!> @param inten output intensity list, including weak beams
!
!> @date  10/13/98 MDG 1.0 original
!> @date   7/04/01 MDG 2.0 f90
!> @date  04/29/13 MDG 3.0 inclusion of Bethe weak beams
!> @date  04/21/14 MDG 4.0 new library version
!--------------------------------------------------------------------------
subroutine CalcBWint(nn,nw,nt,thick,kn,inten)

use local
use io
use diffraction
use crystalvars
use crystal
use kvectors
use gvectors
use dynamical
use constants

IMPLICIT NONE

integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
integer(kind=irg),INTENT(IN)    :: nw                   !< number of weak beams
integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
real(kind=dbl),INTENT(IN)       :: thick(nt)            !< thickness array
real(kind=sgl),INTENT(IN)       :: kn			   !< normal component of wave vector
real(kind=dbl),INTENT(INOUT)    :: inten(nt,nn+nw)      !< output intensities (both strong and weak)

integer(kind=irg)               :: i,j,IPIV(nn), ll(3), jp
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), &
                                lalpha(nn), delta(nn,nn), weak(nw,nn), Ucross(nw,nn), tmp(nw,nn), c
real(kind=sgl)                  :: th
type(reflisttype),pointer	:: rlr, rlw


! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
 lW = cPi*Wloc/cmplx(kn,0.0)
 lalpha(1:nn) = CGinv(1:nn,1)

! in preparation for the intensity computation, we need the prefactor array for the
! weak beam amplitude computation...
! we will also need a potential coefficient array for the cross coefficients, which we 
! compute in the same loop
if (nw.ne.0) then
 	rlr => cell%reflist%next
 	j = 1
	do
	  if (.not.associated(rlr)) EXIT
 	  rlw => cell%firstw
 	  jp = 1
	  do
	    if (.not.associated(rlw)) EXIT
! prefactor value
            c = cmplx(2.D0*rlw%sg/mLambda) - 2.D0*kn*Wloc(j)
            weak(jp,j) = cmplx(-1.D0,0.D0)/c
! cross potential coefficient
            ll(1:3) = rlw%hkl - rlr%hkl
            Ucross(jp,j) = cell%LUT( ll(1),ll(2),ll(3) )
 	    rlw => rlw%nextw
	    jp = jp + 1
 	  end do
 	  rlr => rlr%nexts
          j = j + 1
        end do
end if

! compute the strong beam intensities, stored in the first nn slots of inten 
! we can also compute the weak beams, since they make use of the same diag(1:nn) expression
! as the strong beams, plus a few other factors (excitation error, wave length, Fourier coefficients)
 do i=1,nt
  th = sngl(thick(i))
  diag=exp(-th*aimag(lW))*cmplx(cos(th*real(lW)),sin(th*real(lW)))*lalpha

! the delta array is common to the strong and weak beam intensity computation, so we compute it first
  do j=1,nn
!   delta(j,1:nn) = lCG(j,1:nn)*diag(1:nn)
   delta(1:nn,j) = lCG(1:nn,j)*diag(j)
  end do
! strong beams
  do j=1,nn
   inten(i,j) = cabs(sum(delta(j,1:nn)))**2
  end do 

  if (nw.ne.0) then  
! weak beams
   tmp = matmul(Ucross,delta)
   do jp=1,nw
    inten(i,nn+jp) = cabs( sum(weak(jp,1:nn)*tmp(jp,1:nn)) )**2
   end do  
  end if
  
 end do
   
end subroutine CalcBWint


!--------------------------------------------------------------------------
!
! SUBROUTINE: GetDynMat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the dynamical matrix, including Bethe potentials
!
!> @param 
!
!> @date  04/22/14 MDG 1.0 new library version
!--------------------------------------------------------------------------
subroutine GetDynMat()

use local
use io
use crystalvars
use crystal
use diffraction
use kvectors
use gvectors
use dynamical
use constants

IMPLICIT NONE


complex(kind=dbl)                :: czero, ughp, uhph, weaksum 
real(kind=dbl)                   :: weaksgsum
integer(kind=sgl)                :: ir, ic, ll(3), istat
type(reflisttype),pointer	   :: rlr, rlc, rlw


czero = cmplx(0.0,0.0,dbl)	! complex zero

! create the dynamical matrix for his reference case
! get the absorption coefficient
	if (allocated(DynMat)) deallocate(DynMat)
	allocate(DynMat(cell%nns,cell%nns),stat=istat)
	DynMat = czero
write (*,*) ' allocated ',cell%nns, cell%nnw
	call CalcUcg( (/0,0,0/) )
	DynUpz = rlp%Vpmod

 	rlr => cell%reflist%next
 	ir = 1
	do
	  if (.not.associated(rlr)) EXIT
   	  rlc => cell%reflist%next
	  ic = 1
	  do
	  if (.not.associated(rlc)) EXIT
          if (ic.ne.ir) then  ! not a diagonal entry
! here we need to do the Bethe corrections if necessary
            if (cell%nnw.ne.0) then
              weaksum = czero
 	      rlw => cell%firstw
	      do
	        if (.not.associated(rlw)) EXIT
               ll = rlr%hkl - rlw%hkl
               ughp = cell%LUT(ll(1),ll(2),ll(3)) 
               ll = rlw%hkl - rlc%hkl
               uhph = cell%LUT(ll(1),ll(2),ll(3)) 
               weaksum = weaksum +  ughp * uhph *cmplx(1.D0/rlw%sg,0.0,dbl)
 	       rlw => rlw%nextw
              end do
!        ! and correct the dynamical matrix element to become a Bethe potential coefficient
              ll = rlr%hkl - rlc%hkl
              DynMat(ir,ic) = cell%LUT(ll(1),ll(2),ll(3))  - cmplx(0.5D0*mLambda,0.0D0,dbl)*weaksum
	     else
              ll = rlr%hkl - rlc%hkl
              DynMat(ir,ic) = cell%LUT(ll(1),ll(2),ll(3))
            end if
          else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! determine the total contribution of the weak beams
            if (cell%nnw.ne.0) then
              weaksgsum = 0.D0
 	      rlw => cell%firstw
	      do
	       if (.not.associated(rlw)) EXIT
                ll = rlr%hkl - rlw%hkl
                ughp = cell%LUT(ll(1),ll(2),ll(3)) 
                weaksgsum = weaksgsum +  cabs(ughp)**2/rlw%sg
 	        rlw => rlw%nextw
              end do
              weaksgsum = weaksgsum * mLambda/2.D0
              DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda-weaksgsum,DynUpz,dbl)
            else
              DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda,DynUpz,dbl)
	    end if	     
	
	   end if	
	   rlc => rlc%nexts
	   ic = ic + 1
	  end do	
	  rlr => rlr%nexts
	  ir = ir+1
	end do

end subroutine GetDynMat
