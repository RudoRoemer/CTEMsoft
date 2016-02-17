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
! EMsoft:EMdymodHDF.f90
!--------------------------------------------------------------------------
!
! MODULE: EMdymodHDF
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines that can be called by external code
!
!> @date  10/16/15 MDG 1.0 original
!> @date  01/11/15 MDG 2.0 split from EMdymod (everything that depends on HDF library)
!--------------------------------------------------------------------------
module EMdymodHDF

use EMdymod

!abstract interface
!        subroutine func (nipar, nfpar, ninit, ipar, fpar, initmeanval, expt, n, x, f, fname)  !! calfun interface

!            use local

!            implicit none

!            integer(8),intent(in)                :: ipar(nipar)
!            real(sgl),intent(inout)              :: fpar(nfpar)
!            real(sgl),intent(in)                 :: initmeanval(ninit)
!            integer(irg),intent(in)              :: nipar
!            integer(irg),intent(in)              :: nfpar
!            integer(irg),intent(in)              :: ninit
!            real(kind=sgl),intent(in)            :: expt(ipar(2)*ipar(3))
!            integer(irg),intent(in)              :: n
!            real(dbl),dimension(:),intent(in)    :: x
!            real(dbl),intent(out)                :: f
!            character(fnlen),intent(in),optional :: fname(2)

!    end subroutine func
!end interface

contains


!--------------------------------------------------------------------------
!
! SUBROUTINE:SinglePEDPattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief This function can be called as a standalone function to compute a kinematical PED pattern
!
!> @etails The main purpose of this routine and its accompanying wrapper routine is to
!> provide a way for an external program to compute a kinematical PED pattern.  The idea is that 
!> all the necessary arrays and variables are passed in by reference as arguments, without
!> the need for the routine to fetch any other data from files etc...  The initial goal is
!> to have a function that can be called with the CALL_EXTERNAL mechanism in IDL, but 
!> in the long run this will also be the approach for calling the routine from C/C++, which
!> is an essential part of integration with DREAM.3D.  
!>
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param cpar array with string parameters
!> @param PEDPattern output array
!> @param quats array of quaternions
!
!> @date 11/15/15 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine SinglePEDPattern(ipar, fpar, cpar, PEDpattern, quats)

! the input parameters are all part of a ipar and fpar input arrays instead of the usual namelist structures.
! The following is the mapping:
!
! ipar(1) = 1 if GetVectorsCone detector arrays need to be computed, 0 if not (arrays will have save status)
! ipar(2) = pednl%npix
! ipar(3) = numquats
! ipar(4) = 
! ipar(5) = 
! ipar(6) = 

! fpar(1) = pednl%dmin
! fpar(2) = pednl%voltage
! fpar(3) = pednl%rnmpp
! fpar(4) = pednl%thickness

! cpar(1) = pednl%xtalname


use local
use typedefs
use dictmod
use crystal
use initializersHDF
use gvectors
use io
use diffraction
use symmetry
use quaternions
use NameListTypedefs
use constants
use rotations
use so3
use math
use,INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

integer(c_size_t),PARAMETER             :: nipar=6
integer(c_size_t),PARAMETER             :: nfpar=4
integer(c_size_t),PARAMETER             :: ncpar=1
integer(c_size_t),PARAMETER             :: nq=4
integer(c_size_t),INTENT(IN)            :: ipar(nipar)
real(kind=sgl),INTENT(IN)               :: fpar(nfpar)
character(fnlen),INTENT(IN)             :: cpar(ncpar)
real(kind=sgl),INTENT(OUT)              :: PEDpattern(ipar(2),ipar(2),ipar(4))
real(kind=sgl),INTENT(IN)               :: quats(nq,ipar(4))

type(unitcell),pointer,SAVE             :: cell
type(DynType),save                      :: Dyn
type(gnode),save                        :: rlp
type(reflisttype),pointer,SAVE          :: reflist, nexts, rltmpa
real(kind=sgl),SAVE                     :: Igmax, xgmin, sgmax, rnmpp
real(kind=sgl),allocatable,SAVE         :: xx(:,:), yy(:,:), dot(:,:)
integer(kind=irg),SAVE                  :: ww, tdp, nsize, nref

integer(kind=irg)                       :: pgnum, ival, i, j
real(kind=sgl)                          :: la, dval, dmin, glen, gmax, om(3,3), k(3), FN(3), Ig, & 
                                           maxint, w, ku(3), kp(3), dx, dy, eu(3), tstart, tstop, x, y, ma, mi
integer(kind=irg)                       :: gp(3), imh, imk, iml, gg(3), ix, iy, iz, sx, sy, ip
real(kind=sgl),allocatable              :: line(:), pedp(:,:)


! should we initialize the necessary arrays ? [this should be done the first time we call this routine]
if (ipar(1).gt.0) then
  sgmax = 0.50

!=============================================
! crystallography section
  nullify(cell)
  allocate(cell)
  call Initialize_Cell(cell,Dyn,rlp,cpar(1), fpar(1), fpar(2))

!=============================================
! generation of all potential reflections inside a reciprocal space sphere
! computed from the camera length and the detector size ...

! first set the maximum |g| value that can possibly give rise to a diffracted beam on the detector (diagonal)
  gmax = sqrt(2.0) * float(ipar(2)) * fpar(3)

! this code is taken from the Initialize_ReflectionList routine, but we do not
! need everything from that routine; first get the size of the lookup table
  gp = shape(cell%LUT)
  imh = (gp(1)-1)/4
  imk = (gp(2)-1)/4
  iml = (gp(3)-1)/4

! initialize the reflection list
  nullify(reflist)
  nullify(rltmpa)
  nref = 0
 
! transmitted beam always has excitation error zero
  gg = (/ 0,0,0 /)
  call AddReflection(rltmpa, reflist, cell, nref, gg)   ! this guarantees that 000 is always the first reflection
  rltmpa%xg = 0.0
  xgmin = 100000.0
  Igmax = 0.0

! now compute |U_g|^2 for all allowed reflections and place the values in the linked reflist; 
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
        if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
         gg = (/ ix, iy, iz /)
         glen = CalcLength(cell, float(gg), 'r' )

! find all reflections, ignoring double diffraction spots
         if ((IsGAllowed(cell,gg)).and.(glen.le.gmax)) then ! allowed by the lattice centering, if any
            call AddReflection(rltmpa, reflist, cell, nref, gg )
! we'll use the sangle field of the rltail structure to store |Ug|^2; we will also need the extinction distance
            rltmpa%sangle = cdabs(cell%LUT(ix, iy, iz))**2
            if (rltmpa%sangle.gt.Igmax) Igmax = rltmpa%sangle
            rltmpa%xg = 1.0/(cdabs(cell%LUT(ix,iy,iz))*cell%mLambda)
            if (rltmpa%xg.lt.xgmin) xgmin = rltmpa%xg
         end if ! IsGAllowed
        end if
       end do izl
      end do iyl
    end do ixl

!=============================================
! create the coordinate arrays for the Gaussian peaks that will represent the diffraction spots
  rnmpp = 1.0/fpar(3)
  ww = 4
  tdp = 2*ww+1
  nsize = ipar(2)/2 + ww 
  allocate(xx(-ww:ww,-ww:ww), yy(-ww:ww,-ww:ww), line(-ww:ww), dot(-ww:ww,-ww:ww))
  line = (/ (float(i),i=-ww,ww) /) * rnmpp
  xx = spread(line,dim=1,ncopies=tdp)
  yy = transpose(xx)
  deallocate(line)
end if

! here we start the actual pattern calculation ... 
allocate(pedp(-nsize:nsize,-nsize:nsize))

do ip=1,ipar(3)   ! loop over all requested orientations
  pedp = 0.0
! convert quaternion to orientation matrix
  om = qu2om(quats(1:4,ip)) 
! determine the appropriate wave vector and foil normal
  k = (/ 0.0, 0.0, 1.0 /)
  ku = matmul(om,k)
  FN = ku
  k = ku/sngl(cell%mLambda)

! first we go through the entire reflection list and compute the excitation errors
! those points that satisfy the cutoff are linked via the nexts pointers
    rltmpa => reflist%next
    nexts => rltmpa
    do j=1,nref
      gg = rltmpa%hkl
      rltmpa%sg = Calcsg(cell,float(gg),k,FN)
! should we consider this point any further ? If so, add it to the strong reflection linked list
      if (abs(rltmpa%sg).le.sgmax) then 
        nexts%nexts => rltmpa
        nexts => rltmpa
      end if
      rltmpa => rltmpa%next
    end do

! then, for each point in the nexts list, we compute the components of k' = k+g+s
! and place them in the proper reference frame; we skip the incident beam since it is 
! meaningless in the kinematical approximation
  nexts => reflist%next%nexts
  om = transpose(om)
  do 
! determine the vector k'
    kp = k + float(nexts%hkl) + nexts%sg*ku
    kp = matmul(om,kp)

! get the intensity for each point
    w = sngl(cPi)*nexts%sg*fpar(4)
    if (abs(w).lt.1.0e-6) then
      Ig = nexts%sangle  ! * (sngl(cPi)*pednl%thickness/nexts%xg)**2
    else 
      Ig = nexts%sangle * (sin(w)/w)**2 ! * (sngl(cPi)*pednl%thickness/nexts%xg)**2
    end if

! determine the spot coordinates on the detector
    x = rnmpp * kp(1)
    y = rnmpp * kp(2)

! and plot that spot as a small Gaussian in the pedpattern array, assuming it falls on the detector.
    if ((abs(x).le.nsize-ww).and.(abs(y).le.nsize-ww)) then
      sx = nint(x)
      sy = nint(y)
      dx = x-sx
      dy = y-sy
      dot = (Ig/Igmax)**0.2 * exp(-((xx-dx)**2+(yy-dy)**2)*0.003)
      pedp(sx-ww:sx+ww,sy-ww:sy+ww) = pedp(sx-ww:sx+ww,sy-ww:sy+ww) + dot(-ww:ww,-ww:ww)
    end if

! and repeat this until the end of the list
    if (.not. associated(nexts%nexts)) EXIT
    nexts => nexts%nexts
  end do

! put the pedp pattern in the output array
  PEDpattern(1:ipar(2),1:ipar(2),ip) = pedp

! reset the nexts linked list and start over
  nexts => reflist%next
  rltmpa => nexts%nexts
  do 
    nullify(nexts%nexts)
    if (.not. associated(rltmpa%nexts)) EXIT
    nexts => rltmpa
    rltmpa => rltmpa%nexts
  end do
end do

deallocate(pedp)

end subroutine SinglePEDPattern

!--------------------------------------------------------------------------
!
! SUBROUTINE:SinglePEDPatternWrapper
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief wrapper routine for SinglePEDPattern; nearly identical to ECP case
!>
!> see https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/Gk0xxVFbW8E
!
!> @param argc number of argument
!> @param argv pointers to subroutine parameters
!
!> @date 11/15/15 MDG 1.0 first version
!--------------------------------------------------------------------------
function SinglePEDPatternWrapper(argc, argv) bind(c, name='SinglePEDPatternWrapper') 

use,INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

INTEGER(c_size_t), VALUE, INTENT(IN)            :: argc 
type(c_ptr), dimension(argc), INTENT(INOUT)     :: argv
REAL(c_float)                                   :: SinglePEDPatternWrapper

! wrapper function dependent declarations; they are all pointers 
! since we pass everything by reference from IDL 
integer(c_size_t)                               :: nipar, nfpar, ncpar, nq
integer(c_size_t),dimension(:), pointer         :: ipar
real(c_float), dimension(:), pointer            :: fpar
character(c_char), dimension(:), pointer        :: cpar
real(c_float), dimension(:,:), pointer          :: quats
real(c_float), dimension(:,:,:), pointer        :: PEDPattern

! the following line just helps in identifying the correct order of the subroutine arguments...
!                              1      2     3       4        5
!subroutine SinglePEDpattern(ipar, fpar, cpar, PEDPattern, quats)
!
! transform the C pointers above to fortran pointers, and use them in the regular function call
nipar = 6
nfpar = 1
nq = 4
call c_f_pointer(argv(1),ipar,(/nipar/)) 
call c_f_pointer(argv(2),fpar,(/nfpar/)) 
call c_f_pointer(argv(3),cpar,(/ncpar/)) 
call c_f_pointer(argv(4),PEDpattern,(/ipar(2),ipar(2),ipar(4)/))
call c_f_pointer(argv(5),quats,(/nq,ipar(4)/))

call SinglePEDPattern(ipar, fpar, cpar, PEDpattern, quats)

SinglePEDPatternWrapper = 1._c_float
end function SinglePEDPatternWrapper



!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDcalfun
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief This function is used by bobyqa to fit an EBSD pattern
!
!> @etails The main purpose of this routine is to calculte the difference of 1 with the dot
!> product of an experimental pattern with the given set of detector parameters. This is used
!> by bobyqa module to fit an EBSD pattern.
!>
!> This routine will first compute the detector arrays rgx etc. if necessary, and then perform
!> the usual interpolation from the square Lambert projection. The pattern will be a basic pattern,
!> without any intensity scaling or binning etc; the calling program should take care of those 
!> operations.
!
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param initmeanval mean value of search space
!> @param EBSDpattern output array
!> @param quats quaternion input array
!> @param accum_e array with Monte Carlo histogram
!> @param mLPNH Northern hemisphere master pattern
!> @param mLPSH Southern hemisphere master pattern
!
!> @date 12/12/15 SS 1.0 original
!--------------------------------------------------------------------------

recursive subroutine EBSDcalfun(nipar, nfpar, ninit, ipar, fpar, initmeanval, expt, n, x, f, fname, verbose)

! the input parameters are all part of a ipar and fpar input arrays instead of the usual namelist structures.
! The following is the mapping:
!
! ipar(1) = 1 if hdf5 file needs to be read, 0 if not (arrays will have save status)
! ipar(2) = detnumsx
! ipar(3) = detnumsy
! ipar(4) = detnumEbins
! ipar(5) = mcnsx
! ipar(6) = mpnpx
! ipar(7) = numset
! ipar(8) = numquats
! ipar(9) = 0/1 ;0 for no mask, 1 for mask
! ipar(10) = binning

! fpar(1) = enl%xpc
! fpar(2) = enl%ypc
! fpar(3) = enl%delta
! fpar(4) = enl%MCsig
! fpar(5) = enl%omega
! fpar(6) = enl%thetac
! fpar(7) = enl%L
! fpar(8) = enl%beamcurrent
! fpar(9) = enl%dwelltime
! fpar(10) = enl%gammavalue
! fpar(11) = maskradius

! initmeanval(1) = L
! initmeanval(2) = phi1
! initmeanval(3) = phi
! initmeanval(4) = phi2

            
use local
use rotations
use constants
use HDF5
use HDFsupport
use distortion
use filters
use,INTRINSIC :: ISO_C_BINDING
           
implicit none
            
integer(c_size_t),intent(in)            :: ipar(nipar)
real(sgl),intent(inout)                 :: fpar(nfpar)
real(sgl),intent(in)                    :: initmeanval(ninit)
integer(irg),intent(in)                 :: nipar
integer(irg),intent(in)                 :: nfpar
integer(irg),intent(in)                 :: ninit
real(sgl),intent(in)                    :: expt(ipar(2)*ipar(3)/ipar(10)/ipar(10))
integer(irg),intent(in)                 :: n
real(dbl),dimension(:),intent(in)       :: x
real(dbl),intent(out)                   :: f
character(fnlen),intent(in),optional    :: fname(2)
logical,intent(in),optional             :: verbose

integer(kind=irg)                       :: nnx, nny, binx, biny
complex(dbl)                            :: D
real(kind=sgl)                          :: quats(4,1), bindx, ma, mi
integer(kind=irg),allocatable           :: accum_e(:,:,:)
real(kind=sgl),allocatable              :: EBSDpattern(:,:,:), binned(:,:)
real(kind=sgl),allocatable              :: EBSDpatternintd(:,:)
integer(kind=irg),allocatable           :: EBSDpatterninteger(:,:), EBSDpatternad(:,:)

real(kind=sgl),allocatable,save         :: accum_e_real(:,:,:)
real(kind=sgl),allocatable,save         :: mLPNH(:,:,:,:)
real(kind=sgl),allocatable,save         :: mLPSH(:,:,:,:)

! variables that must be saved for the next time this function is called
real(kind=sgl)                          :: prefactor

! other variables
real(kind=sgl),parameter                :: dtor = 0.0174533  ! convert from degrees to radians
real(kind=sgl)                          :: ixy(2), eu(3)
real(kind=sgl), allocatable             :: EBSDvector(:), EBSDflip(:,:), mask(:,:)
integer(kind=irg)                       :: i, j, istat

logical                                 :: stat, readonly
integer(kind=irg)                       :: hdferr, nlines
integer(HSIZE_T)                        :: dims(1), dims4(4), dims3(3)
character(fnlen)                        :: groupname, dataset, masterfile

type(HDFobjectStackType),pointer        :: HDF_head


fpar(1) = sngl(X(1))*10.0 - 5.0 + initmeanval(5) ! xpc +/- 5 pixels
fpar(2) = sngl(X(2))*10.0 - 5.0 + initmeanval(6) ! ypc +/- 5 pixels
fpar(5) = sngl(X(3))*10.0 - 5.0 ! omega 0 +/- 5 degrees
fpar(7) = sngl(X(4))*2000.0 -1000.0 + initmeanval(1) ! mean +/- 1000 microns

eu = (/X(5)*10.0 - 5.0 + initmeanval(2), X(6)*10.0 - 5.0 + initmeanval(3), X(7)*10.0 - 5.0 + initmeanval(4)/)*dtor ! mean +/- 10 degrees
quats(1:4,1) = eu2qu(eu)

D = dcmplx(X(8)*0.000002D0 - 0.000001D0 + dble(initmeanval(5)), X(9)*0.000002D0 - 0.000001D0 + dble(initmeanval(6)))

binx = ipar(2)/ipar(10)
biny = ipar(3)/ipar(10)
bindx = 1.0/float(ipar(10)**2)
if(present(fname)) then
   if (ipar(1) .eq. 1) then
! open the fortran HDF interface
     call h5open_f(hdferr)

     nullify(HDF_head, HDF_head)

! is this a propoer HDF5 file ?
     call h5fis_hdf5_f(trim(fname(1)), stat, hdferr)

     if (stat) then 
! open the master file 
       readonly = .TRUE.
       hdferr =  HDF_openFile(fname(1), HDF_head, readonly)

       groupname = 'EMData'
       hdferr = HDF_openGroup(groupname, HDF_head)

       dataset = 'mLPNH'
       mLPNH = HDF_readDatasetFloatArray4D(dataset, dims4, HDF_head)

       dataset = 'mLPSH'
       mLPSH = HDF_readDatasetFloatArray4D(dataset, dims4, HDF_head)

       call HDF_pop(HDF_head)

     end if

! is this a propoer HDF5 file ?
     call h5fis_hdf5_f(trim(fname(2)), stat, hdferr)

     if (stat) then 

! open the MC file using the default properties.
       readonly = .TRUE.
       hdferr =  HDF_openFile(fname(2), HDF_head, readonly)

! open the Data group
       groupname = 'EMData'
       hdferr = HDF_openGroup(groupname, HDF_head)

       dataset = 'accum_e'
       accum_e = HDF_readDatasetIntegerArray3D(dataset, dims3, HDF_head)
       nnx = (dims3(2)-1)/2
       allocate(accum_e_real(1:dims3(1),-nnx:nnx,-nnx:nnx))
! and close everything
       call HDF_pop(HDF_head,.TRUE.)

     end if

! close the fortran HDF interface
     call h5close_f(hdferr)

     accum_e_real = float(accum_e)
     deallocate(accum_e)

   end if

end if

if (present(verbose)) then
    if(verbose) then    
        print*,'xpc, ypc, omega, L, eu = ',fpar(1),fpar(2),fpar(5),fpar(7), X(5)*10.0 - 5.0 + initmeanval(2), &
        X(6)*10.0 - 5.0 + initmeanval(3), X(7)*10.0 - 5.0 + initmeanval(4)
    end if
end if

allocate(EBSDvector(binx*biny),mask(binx,biny))
allocate(EBSDpattern(ipar(2),ipar(3),1))
allocate(binned(binx,biny))
allocate(EBSDpatternintd(ipar(2),ipar(3)),EBSDpatterninteger(ipar(2),ipar(3)), EBSDpatternad(ipar(2),ipar(3)))
binned = 0.0
EBSDpatternintd = 0.0
EBSDpatterninteger = 0
EBSDpatternad = 0

mask = 1.0
do i = 1,binx
    do j = 1,biny
        if(((float(i)-ceiling(float(binx)/2.0))**2 + (float(j)-ceiling(float(biny)/2.0))**2) .gt. fpar(11)**2) then
            mask(i,j) = 0.0
        end if
    end do
end do

call getEBSDPatterns(ipar, fpar, EBSDpattern, quats, accum_e_real, mLPNH, mLPSH)

nnx = ipar(2)
nny = ipar(3)

ma = maxval(EBSDpattern)
mi = minval(EBSDpattern)

EBSDpatternintd = ((EBSDPattern(:,:,1) - mi)/ (ma-mi))
EBSDpatterninteger = nint(EBSDpatternintd*255.0)
EBSDpatternad =  adhisteq(10,nnx,nny,EBSDpatterninteger)
EBSDPattern(:,:,1) = float(EBSDpatternad)


if (ipar(10) .ne. 1) then
    do i=1,binx
        do j=1,biny
            binned(i,j) = sum(EBSDpattern((i-1)*ipar(10)+1:i*ipar(10),(j-1)*ipar(10)+1:j*ipar(10),1))
        end do
    end do 
    binned = binned * bindx 
else
    binned(1:binx,1:biny) = EBSDpattern(1:binx,1:biny,1)
end if

if (ipar(9) .eq. 1) then
    binned(1:binx,1:biny) = binned(1:binx,1:biny)*mask(1:binx,1:biny)
end if

do i=1,binx
    do j=1,biny
        EBSDvector((i-1)*biny+j) = binned(j,i)
    end do
end do


EBSDvector = EBSDvector**fpar(10)
EBSDvector = EBSDvector/NORM2(EBSDvector)

F = 1.0 - DOT_PRODUCT(EBSDvector,expt)

end subroutine EBSDcalfun

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECPcalfun
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief This function is used by bobyqa to fit an EBSD pattern
!
!> @etails The main purpose of this routine is to calculte the difference of 1 with the dot
!> product of an experimental pattern with the given set of detector parameters. This is used
!> by bobyqa module to fit an EBSD pattern.
!>
!> This routine will first compute the detector arrays rgx etc. if necessary, and then perform
!> the usual interpolation from the square Lambert projection. The pattern will be a basic pattern,
!> without any intensity scaling or binning etc; the calling program should take care of those 
!> operations.
!
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param initmeanval array with mean value of search space
!> @param ECPattern output array
!> @param quats array of quaternions
!> @param accum_e array with Monte Carlo histogram
!> @param mLPNH Northern hemisphere master pattern
!> @param mLPSH Southern hemisphere master pattern
!
!> @date 12/12/15  SS 1.0 original
!--------------------------------------------------------------------------

recursive subroutine ECPcalfun (nipar, nfpar, ninit, ipar, fpar, initmeanval, expt, n, x, f, fname, verbose)

! the input parameters are all part of a ipar and fpar input arrays instead of the usual namelist structures.
! The following is the mapping:
!
! ipar(1) = 1 if hdf5 file needs to be read, 0 if not (arrays will have save status)
! ipar(2) = detnumsx
! ipar(3) = detnumsy
! ipar(4) = numangle
! ipar(5) = mcnsx
! ipar(6) = numset
! ipar(7) = mpnpx
! ipar(8) = numquats
! ipar(9) = 0/1 ;0 for no mask, 1 for mask

! fpar(1) = ecpnl%thetac
! fpar(2) = ecpnl%sampletilt
! fpar(3) = ecpnl%workingdistance
! fpar(4) = ecpnl%Rin
! fpar(5) = ecpnl%Rout
! fpar(6) = ecpnl%sigstart
! fpar(7) = ecpnl%sigend
! fpar(8) = ecpnl%sigstep
! fpar(9) = ecpnl%gammavalue
! fpar(10) = maskradius

! initmeanval(1) = thetac
! initmeanval(2) = sampletilt
! initmeanval(3) = working distance 
! initmeanval(4) = phi1
! initmeanval(5) = phi
! initmeanval(6) = phi2

use local
use rotations
use constants
use HDF5
use HDFsupport 
use distortion 
use,INTRINSIC :: ISO_C_BINDING
use filters

IMPLICIT NONE

integer(c_size_t),intent(in)            :: ipar(nipar)
real(sgl),intent(inout)                 :: fpar(nfpar)
real(sgl),intent(in)                    :: initmeanval(ninit)
integer(irg),intent(in)                 :: nipar
integer(irg),intent(in)                 :: nfpar
integer(irg),intent(in)                 :: ninit
real(sgl),intent(in)                    :: expt(ipar(2)*ipar(3))
integer(irg),intent(in)                 :: n
real(dbl),dimension(:),intent(in)       :: x
real(dbl),intent(out)                   :: f
character(fnlen),intent(in),optional    :: fname(2)
logical,intent(in),optional             :: verbose

integer(kind=irg)                       :: nnx, nny
complex(dbl)                            :: D
real(kind=sgl)                          :: quats(4,1), ma, mi
integer(kind=irg),allocatable           :: accum_e(:,:,:)
real(kind=sgl),allocatable              :: ECPpattern(:,:,:)
real(kind=sgl),allocatable,save         :: accum_e_real(:,:,:)

real(kind=sgl),allocatable,save         :: mLPNH(:,:,:)
real(kind=sgl),allocatable,save         :: mLPSH(:,:,:)

real(kind=sgl),allocatable              :: binned(:,:)
real(kind=sgl),allocatable              :: ECPpatternintd(:,:)
integer(kind=irg),allocatable           :: ECPpatterninteger(:,:), ECPpatternad(:,:)
real(kind=sgl),allocatable              :: ECPvector(:), ECPvectorcpy(:), ECPtmp(:,:)
real(kind=sgl),allocatable              :: mask(:,:)

integer(kind=irg)                       :: istat, i, j
real(kind=sgl),parameter                :: dtor = 0.0174533  ! convert from degrees to radians
real(kind=sgl)                          :: eu(3)
logical                                 :: stat, readonly
integer(kind=irg)                       :: hdferr, nlines
integer(HSIZE_T)                        :: dims(1), dims4(4), dims3(3)
character(fnlen)                        :: groupname, dataset, masterfile

type(HDFobjectStackType),pointer        :: HDF_head


fpar(1) = sngl(X(1))*4.0 - 2.0 + initmeanval(1) ! thetac mean +/- 2 degrees degrees only
fpar(2) = sngl(X(2))*2.0 - 1.0 + initmeanval(2) ! mean +/- 2 degrees sampletilt
fpar(3) = sngl(X(3))*2.0 - 1.0 + initmeanval(3) ! workingdistance +/ 2 mm

eu = (/X(4)*10.0 - 5.0 + initmeanval(4), X(5)*10.0 - 5.0 + initmeanval(5), X(6)*10.0 - 5.0 + initmeanval(6)/)*cPi/180.0 ! mean +/- 5 degrees
quats(1:4,1) = eu2qu(eu)
!D = dcmplx(X(7)*0.000001D0 - 0.0000005D0 + dble(initmeanval(7)), X(8)*0.000001D0 - 0.0000005D0 + dble(initmeanval(8)))

! read all the files 
if(present(fname)) then
    
    if (ipar(1) .ge. 1) then
! open the fortran HDF interface
        call h5open_f(hdferr)

        nullify(HDF_head, HDF_head)

! is this a propoer HDF5 file ?
        call h5fis_hdf5_f(trim(fname(1)), stat, hdferr)

        if (stat) then 
! open the master file 
          readonly = .TRUE.
          hdferr =  HDF_openFile(fname(1), HDF_head, readonly)

          groupname = 'EMData'
          hdferr = HDF_openGroup(groupname, HDF_head)

          dataset = 'mLPNH'
          mLPNH = HDF_readDatasetFloatArray3D(dataset, dims3, HDF_head)

          dataset = 'mLPSH'
          mLPSH = HDF_readDatasetFloatArray3D(dataset, dims3, HDF_head)

          call HDF_pop(HDF_head)

        end if

! is this a propoer HDF5 file ?
        call h5fis_hdf5_f(trim(fname(2)), stat, hdferr)

        if (stat) then 

! open the MC file using the default properties.
          readonly = .TRUE.
          hdferr =  HDF_openFile(fname(2), HDF_head, readonly)

! open the Data group
          groupname = 'EMData'
          hdferr = HDF_openGroup(groupname, HDF_head)

          dataset = 'accum_e'
          accum_e = HDF_readDatasetIntegerArray3D(dataset, dims3, HDF_head)
          allocate(accum_e_real(1:dims3(1),1:dims3(2),1:dims3(3)))
! and close everything
          call HDF_pop(HDF_head,.TRUE.)

        end if

! close the fortran HDF interface
        call h5close_f(hdferr)

        accum_e_real = float(accum_e)
        deallocate(accum_e)

    end if
end if

allocate(ECPvector(ipar(2)*ipar(3)),mask(ipar(2),ipar(3)))
allocate(ECPpattern(ipar(2),ipar(3),ipar(8)))
allocate(ECPpatternintd(ipar(2),ipar(3)),ECPpatterninteger(ipar(2),ipar(3)), ECPpatternad(ipar(2),ipar(3)))
ECPpatternintd = 0.0
ECPpatterninteger = 0
ECPpatternad = 0

if (present(verbose)) then
    if(verbose) then    
        print*,'thetac, sampletilt, WD, eu = ',sngl(X(1))*4.0 - 2.0 + initmeanval(1),sngl(X(2))*2.0 - 1.0 + initmeanval(2), &
        sngl(X(3))*2.0 - 1.0 + initmeanval(3),eu*180.0/cPi
    end if
end if

mask = 1.0
do i = 1,ipar(2)
    do j = 1,ipar(3)
        if(((float(i)-ceiling(float(ipar(2))/2.0))**2 + (float(j)-ceiling(float(ipar(3))/2.0))**2) .gt. fpar(10)**2) then
            mask(i,j) = 0.0
        end if
    end do
end do

call getECPatterns(ipar, fpar, ECPpattern, quats, accum_e_real, mLPNH, mLPSH)

nnx = ipar(2)
nny = ipar(3)
!call BarrelDistortion(D, ECPpattern(:,:,1), nnx, nny) 

ma = maxval(ECPpattern)
mi = minval(ECPpattern)

ECPpatternintd = ((ECPPattern(:,:,1) - mi)/ (ma-mi))
ECPpatterninteger = nint(ECPpatternintd*255.0)
ECPpatternad =  adhisteq(10,nnx,nny,ECPpatterninteger)
ECPPattern(:,:,1) = float(ECPpatternad)

allocate(ECPvectorcpy(ipar(2)*ipar(3)),stat=istat)

if (ipar(9) .eq. 1) then
    do i = 1,ipar(8)
        ECPpattern(:,:,i) = ECPpattern(:,:,i)*mask
    end do
end if

do i=1,ipar(2)
    do j=1,ipar(3)
        ECPvectorcpy((i-1)*ipar(3)+j) = ECPpattern(j,i,1)
    end do
end do

ECPvector = 0.0

do i = 1,ipar(2)
    ECPvector((i-1)*ipar(3)+1:i*ipar(3)) = ECPvectorcpy((ipar(2)-i)*ipar(3)+1:(ipar(2)-i+1)*ipar(3))
end do

ECPvector = ECPvector**fpar(9)
ECPvector = ECPvector/NORM2(ECPvector)

F = 1.0 - DOT_PRODUCT(ECPvector,expt)

end subroutine ECPcalfun


!--------------------------------------------------------------------------
!
! SUBROUTINE:efitECPWrapper
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief wrapper routine for fitting ECP pattern
!>
!> see https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/Gk0xxVFbW8E
!>
!
!> @param argc number of argument
!> @param argv pointers to subroutine parameters
!
!> @date 12/15/15  SS 1.0 original
!--------------------------------------------------------------------------
function efitECPWrapper(argc, argv) bind(c, name='efitECPWrapper') 

use,INTRINSIC :: ISO_C_BINDING
use bobyqa_module
use local

IMPLICIT NONE

INTEGER(c_size_t), VALUE, INTENT(IN)            :: argc 
type(c_ptr), dimension(argc), INTENT(INOUT)     :: argv
REAL(c_float)                                   :: efitECPWrapper

! wrapper function dependent declarations; they are all pointers 
! since we pass everything by reference from IDL 
integer(4)                                      :: nipar, nfpar, ninit, n, iprint, maxfun, npt
real(c_double), dimension(:), pointer           :: rhobeg, rhoend
integer(c_size_t),dimension(:), pointer         :: ipar
real(c_float), dimension(:), pointer            :: fpar
character(c_char), dimension(:), pointer        :: fname
real(c_float), dimension(:), pointer            :: initmeanval
real(c_float), dimension(:), pointer            :: expt
real(c_double), dimension(:), pointer           :: X
real(c_double), dimension(:), pointer           :: XL
real(c_double), dimension(:), pointer           :: XU


! the following line just helps in identifying the correct order of the subroutine arguments...
!                                        1      2     3       4           5           6   7   8   9       10
!subroutine BOBYQA(nipar, nfpar, ninit, ipar, fpar, fname, initmeanval, expt, N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT, MAXFUN, ECPCALFUN)
!
! transform the C pointers above to fortran pointers, and use them in the regular function call
nipar = 9
nfpar = 10
ninit = 6
n = 6
iprint = 2
maxfun = 10000
npt = n + 6

call c_f_pointer(argv(1),ipar,(/nipar/)) 
call c_f_pointer(argv(2),fpar,(/nfpar/)) 
call c_f_pointer(argv(3),fname,(/2/))
call c_f_pointer(argv(4),initmeanval,(/n/))
call c_f_pointer(argv(5),expt,(/ipar(2)*ipar(3)/))
call c_f_pointer(argv(6),X,(/n/))
call c_f_pointer(argv(7),XL,(/n/))
call c_f_pointer(argv(8),XU,(/n/))
call c_f_pointer(argv(9),RHOBEG,(/1/))
call c_f_pointer(argv(10),RHOEND,(/1/))

call BOBYQA(NIPAR, NFPAR, NINIT, IPAR, FPAR, FNAME, INITMEANVAL, EXPT, N, NPT, X, XL, XU, RHOBEG(1), RHOEND(1),&
     IPRINT, MAXFUN, ECPCALFUN)

efitECPWrapper = 1._c_float

end function efitECPWrapper

!--------------------------------------------------------------------------
!
! SUBROUTINE:efitEBSDWrapper
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief wrapper routine for fitting EBSD pattern
!>
!> see https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/Gk0xxVFbW8E
!>
!
!> @param argc number of argument
!> @param argv pointers to subroutine parameters
!
!> @date 12/15/15  SS 1.0 original
!--------------------------------------------------------------------------
function efitEBSDWrapper(argc, argv) bind(c, name='efitEBSDWrapper') 

use,INTRINSIC :: ISO_C_BINDING
use bobyqa_module
use local

IMPLICIT NONE

INTEGER(c_size_t), VALUE, INTENT(IN)            :: argc 
type(c_ptr), dimension(argc), INTENT(INOUT)     :: argv
REAL(c_float)                                   :: efitEBSDWrapper

! wrapper function dependent declarations; they are all pointers 
! since we pass everything by reference from IDL 
integer(4)                                      :: nipar, nfpar, ninit, n, iprint, maxfun, npt
real(c_double), dimension(:), pointer           :: rhobeg, rhoend
integer(c_size_t),dimension(:), pointer         :: ipar
real(c_float), dimension(:), pointer            :: fpar
character(c_char), dimension(:), pointer        :: fname
real(c_float), dimension(:), pointer            :: initmeanval
real(c_float), dimension(:), pointer            :: expt
real(c_double), dimension(:), pointer           :: X
real(c_double), dimension(:), pointer           :: XL
real(c_double), dimension(:), pointer           :: XU



! the following line just helps in identifying the correct order of the subroutine arguments...
!                                        1      2     3       4           5           6   7   8   9        10
!subroutine BOBYQA(nipar, nfpar, ninit, ipar, fpar, fname, initmeanval, expt, N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT, MAXFUN, EBSDCALFUN)
!
! transform the C pointers above to fortran pointers, and use them in the regular function call
nipar = 9
nfpar = 11
ninit = 4
n = 7
iprint = 2
maxfun = 10000
npt = n + 6

call c_f_pointer(argv(1),ipar,(/nipar/)) 
call c_f_pointer(argv(2),fpar,(/nfpar/)) 
call c_f_pointer(argv(3),fname,(/2/))
call c_f_pointer(argv(4),initmeanval,(/n/))
call c_f_pointer(argv(5),expt,(/ipar(2)*ipar(3)/))
call c_f_pointer(argv(6),X,(/n/))
call c_f_pointer(argv(7),XL,(/n/))
call c_f_pointer(argv(8),XU,(/n/))
call c_f_pointer(argv(9),RHOBEG,(/1/))
call c_f_pointer(argv(10),RHOEND,(/1/))

!call BOBYQA(NIPAR, NFPAR, NINIT, IPAR, FPAR, FNAME, INITMEANVAL, EXPT, N, NPT, X, XL, XU, RHOBEG(1), RHOEND(1),&
!     IPRINT, MAXFUN, EBSDCALFUN)

efitEBSDWrapper = 1._c_float

end function efitEBSDWrapper

end module EMdymodHDF
