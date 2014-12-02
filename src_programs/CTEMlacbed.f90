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
! CTEMsoft2013:CTEMLACBED.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMLACBED 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Zone axis LACBED
!
!> @todo implement full symmetry use; implement full Bloch wave output
!>
!> implement OpenMP multithreading for the actual computation part; requires modifications
!> in CTEMlib.a routines (mostly THREADPRIVATE commands in several modules)
!
!> @date 11/29/01 MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!> @date 05/08/13 MDG 2.1 forked from mbcbed and adapted for large angle CBED patterns
!> @date 05/14/13 MDG 2.2 replaced all IO by namelist file and added command line argument handling
!> @date 09/04/13 MDG 2.3 all command line argument handling now via files.f90 routine
!> @date 07/01/14 MDG 3.0 removal of all globals; separate handling of namelist; uppercased program name
!--------------------------------------------------------------------------
program CTEMLACBED

use local
use files
use NameListTypedefs
use NameListHandlers
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(KosselNameListType)                :: lacbednl

nmldeffile = 'CTEMLACBED.nml'
progname = 'CTEMLACBED.f90'
progdesc = 'Large angle convergent beam pattern simulation'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 10 /), progname)

! deal with the namelist stuff
call GetLACBEDNameList(nmldeffile,lacbednl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the zone axis computations
call LACBEDpattern(lacbednl, progname)

end program CTEMLACBED

!--------------------------------------------------------------------------
!
! SUBROUTINE:LACBEDpattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a large angle zone axis convergent beam electron diffraction pattern
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 10/04/13  MDG 3.0 adaptation for new symmetry routines and output
!> @date 10/05/13  MDG 3.1 added output stuff for IDL visualization program
!> @date 10/07/13  MDG 3.2 corrected subtle error in reflection numbering
!> @date 10/08/13  MDG 3.3 added maxHOLZ output limitation; disk offset computation
!> @date 10/15/13  MDG 3.4 modified handling of Whole Pattern symmetry; add minten output limitation
!> @date 10/16/13  MDG 3.5 added correct HOLZ coordinate handling
!> @date 07/01/14  MDG 4.0 removal of globals and namelist handling
!--------------------------------------------------------------------------
subroutine LACBEDpattern(lacbenl, progname)

use local
use typedefs
use NameListTypedefs
use constants
use crystal
use diffraction
use gvectors
use kvectors
use MBmodule
use postscript, ONLY: GetIndex
use symmetry
use math
use io
use error
use files
use omp_lib

IMPLICIT NONE

type(LACBEDNameListType),INTENT(IN)     :: lacbednl
character(fnlen),INTENT(IN)             :: progname

real(kind=sgl)                  :: ktmax, io_real(3), bragg, thetac, sc, pxy(2), galen, &
                                   frac,  klaue(2), thetam 
integer(kind=irg)               :: ijmax,ga(3),gb(3),cnt, PX, numthick, ss, icnt, pgnum, ih, nunique, famnum, &
                                   newcount,count_rate,count_max, io_int(6), i, j, isym, ir, skip, ghkl(3), &
                                   npx, npy, numt, numk, npix, ik, ip, jp, istat, dgn, nbeams, &
                                   ifamily, famhkl(3), inum, numksame
character(3)                    :: method
integer(kind=irg)               :: itmp(48,3)                   !< array used for family computations etc

real(kind=sgl),allocatable      :: diskoffset(:,:), disk(:,:,:,:), thick(:), familytwotheta(:), slice(:,:,:)
integer(kind=irg),allocatable   :: familymult(:), familyhkl(:,:), whichHOLZ(:), gequiv(:,:)
real(kind=sgl),allocatable      :: inten(:,:)
real(kind=dbl)                  :: s(3)
logical,allocatable             :: ksame(:)
complex(kind=dbl),allocatable   :: DynMat(:,:)


type(unitcell),pointer          :: cell
type(gnode)                     :: rlp
type(DynType)                   :: Dyn
type(kvectorlist),pointer       :: khead, ktmp
type(symdata2D)                 :: TDPG
type(BetheParameterType)        :: BetheParameters
type(reflisttype),pointer       :: reflist, firstw, rltmp, rltmpa, rltmpb


! camlen = 1000.0                 ! camera length [mm]   (this is not part of the namelist, but needed later)

  nullify(cell)
  nullify(khead)
  nullify(ktmp)

  allocate(cell)

  verbose = .TRUE.
  call Initialize_Cell(cell,Dyn,rlp,lacbednl%xtalname, lacbednl%dmin, lacbednl%voltage, verbose)

! set the foil normal 
  Dyn%FN = float(lacbednl%fn)
 
! determine the point group number
  j=0
  do i=1,32
   if (SGPG(i).le.cell % SYM_SGnum) j=i
  end do

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
  dgn = GetPatternSymmetry(lacbednl%k,j,.TRUE.)
  pgnum = j
  isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
  call ShortestG(lacbednl%k,ga,gb,isym)
  io_int(1:3)=ga(1:3)
  io_int(4:6)=gb(1:3)
  call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! for some of the 2D point groups, the standard orientation of the group according to ITC vol A
! may not be the orientation that we have here, so we need to determine by how much the 2D point
! group is rotated (CCW) with respect to the standard setting...
  call CheckPatternSymmetry(lacbednl%k,ga,isym,thetam)

! initialize the HOLZ geometry type
  call GetHOLZGeometry(cell,float(ga),float(gb),lacbednl%k,lacbednl%fn) 

! construct the list of all possible reflections
! method = 'ALL'
! thetac = convergence/1000.0
! call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxHOLZ,thetac)

  galen = CalcLength(cell,float(ga),'r')

! determine range of incident beam directions
  bragg = CalcDiffAngle(cell,ga(1),ga(2),ga(3))*0.5
  
! convert to ktmax along ga
  ktmax = 0.5*thetac/bragg

! the number of pixels across the disk is equal to 2*npix + 1
  npx = lacbednl%npix
  npy = npx
  io_int(1) = 2.0*npx + 1
  call WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4/)")
  
! get number of thicknesses for which to compute the LACBED disk patterns
  numt = lacbednl%numthick
  allocate(thick(numt),stat=istat)
  thick = lacbednl%startthick + lacbednl%thickinc* (/ (float(i),i=0,numt-1) /)

! set parameters for wave vector computation
  klaue = (/ 0.0, 0.0 /)
  ijmax = float(npx)**2   ! truncation value for beam directions

! determine all independent incident beam directions (use a linked list starting at khead)
!  isym = WPPG(dgn)

! for now, the solution to the symmetry problem is to do the computation for the entire 
! illumination cone without application of symmetry.  Instead, we'll get the speed up by 
! going to multiple cores later on.
  isym = 1
  call CalckvectorsSymmetry(khead,cell,TDPG,dble(lacbednl%k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue,.FALSE.)

! set scaling parameters
  PX = npix/2

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(BetheParameters,.TRUE.)

! construct the list of all possible reflections ...
  call Initialize_ReflectionList(cell, reflist, BetheParameters, float(lacbednl%fn), &
                                 float(lacbednl%k), lacbednl%dmin, nref, verbose)

! up to this point, everything is nearly identical to the mbcbed program,
! except that we do not use a camera length explicitly.  Now we
! need to do things a little differently.  First of all, we need a master 
! list for all the reflections that contribute in one way or another to the
! complete LACBED pattern.  We do this by going through the entire list of incident
! beam directions and flagging all reflections that contribute, either as weak
! or as strong reflections.
  mess = ' Pruning reflection list (this takes a while ...) '
  call Message("(A)")
  call Prune_ReflectionList(numk,nbeams)
  io_int(1) = nbeams
  call WriteValue('Number of contributing beams  : ', io_int, 1, '(I8)')

! since we're only going to store one diffraction disk per family of reflections,
! we need to decide which one we're going to keep.  For those, we'll need to 
! determine the 2D multiplicity and store that along with the Miller indices 
! of a representative family member.

! set the scale parameter for a default camera length of 1000 mm.
  sc = cell%mLambda * 1000.0 * 300.0 / 25.4  ! the absolute value does not matter and is derived from legacy Postscript code
! The original code used 300 dpi (hence 300/25.4) which was convenient for Postscript output; in the current case, we
! do not actually use the true value, but in the IDL visualization program, we scale the user defined camera length by
! 1000.0, and use this ratio to scale the diskoffset coordinates.  So, there's no absolute length scale, only a relative scale.

! Next, we need to go through the pruned reflection list and identify for each 
! reflection the multiplicity with respect to the family generated by the
! 2D whole pattern point group [see manual for IDL visualization program
! for an explicit example].  In other words, we need to tag each reflection
! that will need to end up in the final output file for the LACBED program.
! For instance, for the [112] Copper zone axis orientation, (1 1 -1) and
! (-1 -1 1) do not belong to the same family (since they lie on the mirror
! plane that makes up the 2D Whole Pattern group m.  Therefore, these must 
! be tagged as separate families with appropriate multiplicities, otherwise
! the visualization program will only know of one reflection (1 1 -1) and
! will not be able to generate the other one (-1 -1 1) {which need not be 
! identical to begin with, due to the m symmetry].
! At the same time we need to determine how many independent families there 
! are, as well as their multiplicities.
 
! to do all this requires knowledge of the subset of 3D symmetry operators that keeps 
! the incident beam direction invariant, so we determine this subset first.  This 
! code comes from the CalcStar routine, but we only need a portion of it.
  allocate(ksame(cell%SG%SYM_NUMpt))
  ksame = .FALSE.
  numksame = 1
  ksame(1) = .TRUE.
  
! get all the symmetry operation IDs that leave the zone axis invariant (skip the identity operation)
  do i=2,cell%SG%SYM_NUMpt 
    s = matmul(cell%SG%SYM_direc(i,1:3,1:3),dble(k)) 
    if (sum(abs(s-dble(k))).lt.1.0D-10) then 
      ksame(i) = .TRUE.
      numksame = numksame+1
    end if
  end do

! and output the number of 3D symmetry operators that leave k invariant; they form
! the subset that we are after... This should really coincide with the 2D whole pattern symmetry group,
! but with 3D operators instead of 2D operators; so now we have the actual symmetry matrices.
  io_int(1) = numksame
  call WriteValue('Number of 3D symmetry operators that leave k invariant : ',io_int, 1, "(I3)")
  io_int(1) = PGTWDorder(WPPG(dgn))
  call WriteValue('Order of Whole Pattern point group : ',io_int, 1, "(I3)")
  allocate(gequiv(numksame,3))
  
! to be safe, we need to reset the family number for each of the reflections in the current list
! to zero, to reflect the fact that we haven't considered that reflection yet in terms of 2D symmetry.
  rltmpa => reflist%next
  do while (associated(rltmpa))
    rltmpa%famnum = 0
    if (.not.associated(rltmpa%next)) EXIT
    rltmpa => rltmpa%next
  end do

! next we go through the entire list again, this time with a double loop, and we
! determine for each hkl all the 2D equivalent ones and tag them as belonging to
! the same family in terms of 2D symmetry.  
  rltmpa => reflist%next
  rltmpa%famnum = 1             ! reflection at origin is always its own family
  ifamily = 1
  rltmpa => rltmpa%next
whileloop1: do while (associated(rltmpa))
! only look at points that haven't been tagged yet
   if (rltmpa%famnum.eq.0) then
    ghkl = rltmpa%hkl           ! get the reflection
    call Calc2DFamily(ghkl,ksame,numksame,nunique,itmp)
! we add this point to a new famnum value
    ifamily = ifamily + 1
    rltmpa%famnum = ifamily
    rltmpa%famhkl = ghkl
    famhkl = ghkl
    if (nunique.gt.1) then
! the order is larger than 1, so we need to go through the list and tag all the equivalent ones;
! we need only consider those that have the same original famhkl.
      do i=2,nunique
        rltmpb => rltmpa%next
whileloop2: do while (associated(rltmpb))
! look for this reflection on the list
          if (sum(abs(itmp(i,1:3) - rltmpb%hkl(1:3))).eq.0) then  ! found it!
           rltmpb%famnum = ifamily
           rltmpb%famhkl = rltmpa%famhkl
           EXIT whileloop2
          end if
          if (.not.associated(rltmpb%next)) EXIT whileloop2
          rltmpb => rltmpb%next
        end do whileloop2
      end do
    end if
   end if
   if (.not.associated(rltmpa%next)) EXIT whileloop1
   rltmpa => rltmpa%next
  end do whileloop1

! ok, so there are ifamily families; next we need to store the corresponding
! hkl, and multiplicity, as well as the diffraction angle and the position of 
! the diffraction disk center for a standard camera length.
  allocate(familyhkl(3,ifamily), familymult(ifamily), familytwotheta(ifamily), diskoffset(2,ifamily))

! redo the above loop, sort of, but now fill in the actual data
! we no longer need to keep the famnum entries in the linked list, so we
! can reset those to zero to keep track of points already visited.
  ifamily = 1   ! for the incident beam
  familyhkl(1:3,ifamily) = (/ 0, 0, 0 /)
  familymult(ifamily) = 1
  diskoffset(1:2,ifamily) = (/ 0.0, 0.0 /)
  familytwotheta(1) = 0.0
  rltmpa => reflist%next%next

outerloop2: do while (associated(rltmpa))
  if (rltmpa%famnum.ne.0) then 
    ifamily = ifamily+1
    famnum = rltmpa%famnum
    rltmpa%famnum = 0
    famhkl = rltmpa%famhkl
    familyhkl(1:3,ifamily) = famhkl(1:3)
    familytwotheta(ifamily) = CalcDiffAngle(famhkl(1),famhkl(2),famhkl(3))*1000.0
    familymult(ifamily) = 1
! get the disk offset parameters
    pxy = sc * GetHOLZcoordinates(cell,float(famhkl), (/ 0.0, 0.0, 0.0 /), sngl(cell%mLambda))
    diskoffset(1:2,ifamily) = pxy
  
! and remove the equivalent reflections from the list
    rltmpb => rltmpa%next
whileloop3: do while (associated(rltmpb))
      if (rltmpb%famnum.eq.famnum) then
         familymult(ifamily) = familymult(ifamily) + 1
         rltmpb%famnum = 0
      end if
      if (.not.associated(rltmpb%next)) EXIT whileloop3
      rltmpb => rltmpb%next
    end do whileloop3
! print results for debugging
!write (*,*) ifamily, familyhkl(1,ifamily),familyhkl(2,ifamily),familyhkl(3,ifamily),familymult(ifamily),&
!familytwotheta(ifamily)
   end if
   if (.not.associated(rltmpa%next)) EXIT outerloop2
   rltmpa => rltmpa%next
  end do outerloop2

  io_int(1) = ifamily
  call WriteValue('Maximum number of unique families in output = ', io_int, 1, "(I5)")

! next we create the output array, which has one disk image for each 
! thickness and contributing family.  We make sure that each image is fully
! filled by a diffraction disk (i.e., no empty space along the main axes.
  allocate(disk(-npix:npix,-npix:npix,1:numt,1:ifamily))
  disk=0.0

  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")

! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)


! point to the first beam direction
  ktmp => khead
! loop over all beam orientations, selecting them from the linked list
kvectorloop:  do ik = 1,numk

        ip = -ktmp%i
        jp =  ktmp%j

! determine strong and weak reflections
        call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)

! generate the dynamical matrix
        allocate(DynMat(nns,nns))
        call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)

! compute the dynamical matrix using Bloch waves with Bethe potentials; note that the IgnoreFoilNormal flag
! has been set to .FALSE.; if it is set to .TRUE., the computation of the ZOLZ will still be mostly correct,
! but the excitation errors of the HOLZ reflections will be increasingly incorrect with HOLZ order.  This was
! useful during program testing but should probably be removed as an option altogether...
!       call Compute_DynMat('BLOCHBETHE', ktmp%k, ktmp%kt, .FALSE.)

! allocate the intensity array to include both strong beams and weak beams (in that order)
        allocate(inten(numt,nns+nnw))
        inten = 0.0
 
! solve the dynamical eigenvalue equation and return the intensities of ALL reflections,
! both strong and weak; the weak intensities should also be plotted at the correct locations....
! this uses a new version of the CalcBWint routine that implements both strong and weak beam intensities.
        call CalcBWint(nns,nnw,numt,thick,inten)

! we combine the reflistindex and weakreflistindex into a single list so that we can match 
! each diffracted intensity with the correct diffraction disk in the disks array.
        BetheParameter%reflistindex = BetheParameter%reflistindex + BetheParameter%weakreflistindex

! ok, we have all the intensities.  Next we need to copy the relevant intensities into the slots 
! of the disk array, one for each family.
        rltmpa => reflist%next%next
        inum = 1
        disk(ip,jp,1:numt,inum) = inten(1:numt,1)
        do i=2,DynNbeamsLinked
          if (BetheParameter%reflistindex(i).ne.0) then ! is this a reflection on the current list
! it is, so we need to determine which of the families corresponds to it
            inum = -1
            do ir=2,ifamily
             ss = sum(abs(familyhkl(1:3,ir) - rltmpa%hkl(1:3)))
             if (ss.eq.0) inum = ir
            end do
            if (inum.ne.-1) then        
              disk(ip,jp,1:numt,inum) = inten(1:numt,BetheParameter%reflistindex(i))
            end if
          end if
          rltmpa => rltmpa%next
        end do

! and remove the intensity array
     deallocate(inten)
    
! select next beam direction
   if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

  end do kvectorloop

! the following comment only applies if we use symmetry to determine the computational wedge.
! 
! next we need to apply the symmetry operators to ALL beams, including the special and 
! general diffraction symmetry for dark field disks... this is complicated, since the 
! reflections are not grouped by symmetrically equivalent classes.  There are two
! solutions; either we figure out on the spot which reflections are equivalent, which
! requires running through the entire reflection list several times, or we change the 
! ComputeReflections routine to list reflections consecutively by family.  Either way
! will work, but the latter one may result in a faster algorithm.  Whichever way we
! do this, we will need to make a copy of the entire disk array, so that we don't 
! apply the operators too many times...  To do so, we need to make sure that if one
! single family member has an entry in the disks array, then all equivalent family 
! members must also be in the array... (equivalent with respect to the zone axis WP symmetry)
! This will require a bit of thinking before implementation can begin.

! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F10.4)")

! before we write the output file, we need to determine which reflection families 
! need to be written to the file; this is determined by the maxHOLZ parameter.  So,
! first we determine to which HOLZ layer each family belongs by using the zone 
! equation.   [check special case of hexagonal indices !!!!]
! Along the way, we count the ones up to order maxHOLZ.
! Also, to reduce the size of the output file a bit, we ignore those reflections that
! have a maximum intensity less than minten for the initial thickness value. On the
! other hand, if a reflection is due to double diffraction, then we include it, always.
  allocate(whichHOLZ(ifamily))
  icnt = 0
  do ir=1,ifamily
    whichHOLZ(ir) = iabs(k(1)*familyhkl(1,ir)+k(2)*familyhkl(2,ir)+k(3)*familyhkl(3,ir))
    if (whichHOLZ(ir).le.maxHOLZ) then 
      if ( (maxval(disk(:,:,:,ir)).ge.minten) .or. ( dbdiff(familyhkl(1,ir),familyhkl(2,ir),familyhkl(3,ir))) ) then
        icnt = icnt+1
      else  ! just change the HOLZ value to some large value to make sure it does not get written to the file
        whichHOLZ(ir) = 10
      end if
    end if
  end do    

! the final bit of the program involves dumping all the results into a file,
! binary for now, but HDF5 in the future, for the IDL visualization program 
! to read.
  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
! write the program identifier
  write (dataunit) trim(progname)
! write the version number
  write (dataunit) scversion
! first write the array dimensions
  write (dataunit) 2*npix+1,2*npix+1,numt,icnt
! then the name of the crystal data file
  write (dataunit) xtalname
! the accelerating voltage [V]
  write (dataunit) voltage
! convergence angle [mrad]
  write (dataunit) convergence
! the zone axis indices
  write (dataunit) k
! the foil normal indices
  write (dataunit) fn
! number of k-values in disk
  write (dataunit) numk
! dmin value
  write (dataunit) dmin
! horizontal reciprocal lattice vector
  write (dataunit) ga  
! length horizontal reciprocal lattice vector (need for proper Laue center coordinate scaling)
  write (dataunit) galen
! maximum HOLZ layer in the output file
  write (dataunit) maxHOLZ
! intensity cutoff
  write (dataunit) minten
! eight integers with the labels of various symmetry groups
  write (dataunit) (/ pgnum, PGLaue(pgnum), dgn, PDG(dgn), BFPG(dgn), WPPG(dgn), DFGN(dgn), DFSP(dgn) /)
! write the 2D point group rotation angle
  write (dataunit) thetam
! thickness data
  write (dataunit) startthick, thickinc
! and from here one we write the individual diffraction disks with associated information  
! Miller indices, multiplicity, two-theta [mrad], and position for a reference camera length
! only for those families that belong the HOLZ layers <= maxHOLZ (to keep the file size down a bit)
  allocate(slice(-npix:npix,-npix:npix,1:numt))
  write (dataunit) familyhkl(1:3,1), familymult(1), familytwotheta(1), diskoffset(1:2,1), whichHOLZ(1)
  slice = disk(-npix:npix,-npix:npix,1:numt,1)
  write (dataunit) slice
! we'll write them in reverse order, so that the smaller Miller indices come first (at least for high symmetry structures);
! we'll also write them by HOLZ number
  do ih = 0,maxHOLZ
   do ir = ifamily,2,-1
    if (whichHOLZ(ir).eq.ih) then
      write (dataunit) familyhkl(1:3,ir), familymult(ir), familytwotheta(ir), diskoffset(1:2,ir), whichHOLZ(ir)
      slice = disk(-npix:npix,-npix:npix,1:numt,ir)
      write (dataunit) slice
    end if
   end do
  end do
  close(UNIT=dataunit,STATUS='keep')

  mess = ' Data stored in '//outname; call Message("(/A/)") 
  io_int(1) = maxHOLZ
  call WriteValue('Data includes families of reflections in HOLZ layers 0 through ',io_int,1,"(I2)")
  io_int(1) = icnt
  call WriteValue('Total number of independent families above intensity threshold saved to file : ',io_int,1,"(I4/)")
 
 write (*,*) 'Some statistics :'
 write (*,*) 'Average number of strong beams : ',float(BetheParameter%totstrong)/float(numk)
 write (*,*) '          (min,max) : ',BetheParameter%minstrong,BetheParameter%maxstrong
 write (*,*) 'Average number of weak beams : ',float(BetheParameter%totweak)/float(numk)
 write (*,*) '          (min,max) : ',BetheParameter%minweak,BetheParameter%maxweak


end subroutine LACBEDpattern




