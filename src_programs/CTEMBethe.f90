
program CTEMBethe

use local
use typedefs
use constants
use initializers
use crystal
use files
use error
use io
use gvectors
use diffraction
use MBmodule

IMPLICIT NONE

type(unitcell),pointer          :: cell
type(gnode)                     :: rlp
type(BetheParameterType)        :: BetheParameter
type(DynType)                   :: Dyn
type(reflisttype),pointer       :: reflist, firstw

complex(kind=dbl),allocatable   :: DynMat(:,:)

integer(kind=irg)               :: numarg       !< number of command line arguments
integer(kind=irg)               :: iargc        !< external function for command line
character(fnlen)                :: arg          !< to be read from the command line
real(kind=sgl)                  :: voltage, dmin, k(3), c3
character(fnlen)                :: xtalname, nmlfile
logical                         :: verbose
integer(kind=irg)               :: i,ii, istat, numt, j1, j2, c1min, c1max, c2min, c2max, c1step, c2step, &
                                   nref, ir, nns, nnw, io_int(1)
type(reflisttype),pointer       :: rlr, rlw
complex(kind=dbl)               :: czero, pre
real(kind=sgl)                  :: kn
real(kind=dbl),allocatable      :: thick(:), inten2(:,:), intenref(:,:)
integer(kind=sgl),allocatable   :: hkllist(:,:), rlist(:)

namelist /parameters/ numt, c1min, c1max, c1step, c2min, c2max, c2step, c3, xtalname, dmin, voltage

interface
        recursive subroutine CalcBWintlocal(cell,DynMat,reflist,firstw,nn,nw,nt,thick,kn,inten)

        use local
        use io
        use diffraction
        use crystal
        use kvectors
        use gvectors
        use constants
        
        IMPLICIT NONE
        
        type(unitcell),pointer          :: cell
        complex(kind=dbl),INTENT(IN)    :: DynMat(nn,nn)
        type(reflisttype),pointer       :: reflist, firstw
        integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
        integer(kind=irg),INTENT(IN)    :: nw                   !< number of weak beams
        integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
        real(kind=dbl),INTENT(IN)       :: thick(nt)            !< thickness array
        real(kind=sgl),INTENT(IN)       :: kn                   !< normal component of wave vector
        real(kind=dbl),INTENT(INOUT)    :: inten(nt,nn+nw)      !< output intensities (both strong and weak)
        end subroutine CalcBWintlocal

        subroutine ResetPointers(reflist)

                use local
        use typedefs

        IMPLICIT NONE

        type(reflisttype),pointer       :: reflist, rlr
        end subroutine ResetPointers

end interface

numarg = iargc()
nmlfile = 'CTEMBethe.nml'
i = 1
if (numarg.gt.1) then 
    call getarg(i,arg)
    nmlfile = arg
end if

verbose = .FALSE.

! set default values for all namelist variables
xtalname = 'undefined'
numt = 250
c1min = 4
c1max = 40
c1step = 2
c2min = 20
c2max = 60
c2step = 2
c3 = 500.0
dmin = 0.05
voltage = 200000.0

! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=parameters)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMBethe:',' structure file name is undefined in '//nmlfile)
 end if


allocate(thick(numt), stat=istat)
thick = dble( (/ (i-1,i=1,numt) /) ) * 0.5D0

! initialize some parameters
czero = cmplx(0.0,0.0,dbl)      ! complex zero
pre = cmplx(0.0,cPi,dbl)        ! i times pi

xtalname = 'Cu.xtal'
dmin = 0.01
voltage = 200000.0

allocate(cell)
call Initialize_Cell(cell,Dyn,rlp,xtalname, dmin, voltage, verbose)

open(unit=dataunit,file='testrun.data',status='unknown',form='unformatted')

nullify(reflist)
nullify(firstw)
nns = 0
nnw = 0

do ii=0,0
  Dyn%FN = (/ sin(float(ii)*0.02), 0.0, cos(float(ii)*0.02) /) 
  k = Dyn%FN
  call NormVec(cell,k,'r')
  k = k/sngl(cell%mLambda)
  write (*,*) 'wave vector = ',k
  
! this needs to be modified
  kn = k(3)

  BetheParameter%c3 = c3

! and generate the reflectionlist
  call Initialize_ReflectionList(cell, reflist, BetheParameter, Dyn%FN, k, dmin, nref, verbose)

  io_int(1) = nref
  call WriteValue('Total number of reflections in list : ', io_int, 1, "(I6)")
        
! first, we need to do the full computation, using a whole lot of reflections
! and we'll take that as the "correct" result.  Then we compare similar computations
! with fewer beams and various Bethe cutoff ranges to figure out the range of 
! parameters for which the result remains accurate to 10^(-5)
  BetheParameter%c2 = 20000.0
  BetheParameter%c1 = 20000.0 
  call Apply_BethePotentials(cell, reflist, firstw, BetheParameter, nref, nns, nnw)

! generate the dynamical matrix
  allocate(DynMat(nns,nns))
  call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)
  write (*,*) 'reference dynmat = ',shape(DynMat)

! next, we do a Bloch wave computation and generate the resulting scattered intensities
! for a series of 250 thickness values (increments of 1 nm)
  if (allocated(intenref)) deallocate(intenref)
  allocate(intenref(numt,nns))
  intenref = 0.0
 
! solve the dynamical eigenvalue equation and return the intensities of ALL reflections,
! both strong and weak; the weak intensities should also be plotted at the correct locations....
! this uses a new version of the CalcBWint routine that implements both strong and weak beam intensities.
  call CalcBWintlocal(cell,DynMat,reflist,firstw,nns,nnw,numt,thick,kn,intenref)
  deallocate(DynMat)

! store the data to a file
  write (dataunit) nns, numt
  allocate(hkllist(3,nns),stat=istat)
  rlr => reflist%next
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
innerloop:    do j2=c2min,c2max,c2step
      BetheParameter%c2 = float(j2)  
      if (j2.lt.j1) then
! write negative numbers so that the IDL routine "knows" that this combination was omitted
        write (dataunit) -j1, -j2, -j1, -j2
        CYCLE innerloop
      end if
      
! nullify all the nexts and nextw pointers
      call ResetPointers(reflist)
      nullify(firstw)

! create the new lists      
      call Apply_BethePotentials(cell, reflist, firstw, BetheParameter, nref, nns, nnw)
      if (allocated(rlist)) deallocate(rlist)
      allocate(rlist(nns+nnw),stat=istat)

! dynamical simulation
      allocate(DynMat(nns,nns))
      call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)
      write (*,*) 'ended GetDynMat ', float(j1), float(j2), shape(DynMat), nns, nnw
      if (allocated(inten2)) deallocate(inten2)
      allocate(inten2(numt,nns+nnw))
      inten2 = 0.0

      call CalcBWintlocal(cell,DynMat,reflist,firstw,nns,nnw,numt,thick,kn,inten2)
      deallocate(DynMat)

! make the reference list in the correct order
      rlist = 0
! first the strong reflections
      rlr => reflist%next
      ir = 1
      do
        if (.not.associated(rlr)) EXIT
        rlist(ir) = rlr%famnum
        rlr => rlr%nexts
        ir = ir + 1
      end do
! then the weak reflections      
      rlw => firstw
      do
        if (.not.associated(rlw)) EXIT
        rlist(ir) = rlw%famnum
        rlw => rlw%nextw
        ir = ir + 1
      end do

! and write the result to the data file
      write (dataunit) j1, j2, nns, nnw
      write (dataunit) rlist
      write (dataunit) inten2
    end do innerloop
    write (*,*) 'completed c1 = ',j1
  end do

end do

close(unit=dataunit,status='keep')


end program CTEMBethe


subroutine ResetPointers(reflist)

use local
use typedefs

IMPLICIT NONE

type(reflisttype),pointer       :: reflist, rlr

rlr => reflist%next
do
  if (.not.associated(rlr)) EXIT
  nullify(rlr%nexts)
  nullify(rlr%nextw)
  rlr => rlr%next
end do

end subroutine ResetPointers




!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcBWintlocal
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
!> @date  06/24/14 MDG 4.1 updated for removal of all global variables
!--------------------------------------------------------------------------
recursive subroutine CalcBWintlocal(cell,DynMat,reflist,firstw,nn,nw,nt,thick,kn,inten)

use local
use io
use diffraction
use crystal
use kvectors
use gvectors
use constants

IMPLICIT NONE

type(unitcell),pointer          :: cell
complex(kind=dbl),INTENT(IN)    :: DynMat(nn,nn)
type(reflisttype),pointer       :: reflist, firstw
integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
integer(kind=irg),INTENT(IN)    :: nw                   !< number of weak beams
integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
real(kind=dbl),INTENT(IN)       :: thick(nt)            !< thickness array
real(kind=sgl),INTENT(IN)       :: kn                   !< normal component of wave vector
real(kind=dbl),INTENT(INOUT)    :: inten(nt,nn+nw)      !< output intensities (both strong and weak)

integer(kind=irg)               :: i,j,IPIV(nn), ll(3), jp
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), &
                                lalpha(nn), delta(nn,nn), weak(nw,nn), Ucross(nw,nn), tmp(nw,nn), c
real(kind=sgl)                  :: th
type(reflisttype),pointer       :: rlr, rlw


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
        rlr => reflist%next
        j = 1
        do
          if (.not.associated(rlr)) EXIT
          rlw => firstw
          jp = 1
          do
            if (.not.associated(rlw)) EXIT
! prefactor value
            c = cmplx(2.D0*rlw%sg/cell%mLambda) - 2.D0*kn*Wloc(j)
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
   inten(i,j) = cdabs(sum(delta(j,1:nn)))**2
  end do 

  if (nw.ne.0) then  
! weak beams
   tmp = matmul(Ucross,delta)
   do jp=1,nw
    inten(i,nn+jp) = cdabs( sum(weak(jp,1:nn)*tmp(jp,1:nn)) )**2
   end do  
  end if
  
 end do
   
end subroutine CalcBWintlocal


