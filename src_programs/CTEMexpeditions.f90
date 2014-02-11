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
! CTEMsoft2013:CTEMexpeditions.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMexpeditions 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief CTEMexpeditions computes zone axis defect contrast for a simply inclusion
!
!> @details test program used to generate trialdata sets for the NSF Expeditions pre-proposal
!> The sample measures 128x128 pixels, with DF_L = 2 nm, a step size of 0.5 nm along z for 
!> a total of 128 steps (64 nm thick foil), with a spherical inclusion located at the center
!> of the computational cell.  We generate the structure matrix A for a 9 beam case, identify
!> the beams for the [001] zone axis orientation, print, the 1/q factors, and finally output
!> the complete A-matrix for each slice and every image column in text format, one text file
!> per column, using the column position in the file name.  No images are computed by this 
!> program, only the 9x9 A matrices are output. 
! 
!> @date  03/07/10 MDG  1.0 original, based on STEMdefect.f90 (collaboration with OSU) 
!> @date   1/16/14 MDG  1.1 forked for NSF Expeditions program
!--------------------------------------------------------------------------
program CTEMexpeditions 

use local
use files
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMexpeditions.nml'
progname = 'CTEMexpeditions.f90'
call Interpret_Program_Arguments(nmldeffile,8,(/ 0, 1, 2, 3, 200, 201, 202, 203 /) )

! initialize all user-defined variables
call ComputeExpeditions(nmldeffile)

end program CTEMexpeditions


!--------------------------------------------------------------------------
!
! SUBROUTINE: ComputeExpeditions
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a near zone axis A matrix data set
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date  1/16/14  MDG 1.1 adapted for NSF Expeditions test run
!--------------------------------------------------------------------------
subroutine ComputeExpeditions(nmlfile)

use local
use error
use crystalvars 
use crystal
use symmetryvars
use symmetry
use postscript
use constants
use diffraction
use dynamical
use kvectors
use gvectors
use files
use io
use math
use foilmodule
use stacking_fault
use dislocation
use void
use inclusion
use defectmodule
use rotations
use TIFF_f90
use pgm
use timing
use STEMmodule

character(fnlen),INTENT(IN)		:: nmlfile

integer(kind=irg)    		:: nn,i,j,k,npix,npiy,ii,jj,numvoids,numdisl, numYdisl, &
					numsf,nCL,numinc,dinfo,t_interval, &
					DF_nums_new,DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, &
					NTHR, isym, ir, ga(3), gb(3),kk(3),ic,g,numd,ix,iy,iz, &
					numk,ixp,iyp,SETNTHR, io_int(6), skip, gg(3), iSTEM, gmgp(3)
!                                  OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
integer(kind=irg),parameter 		:: numdd=180
real(kind=sgl)         		:: glen,exer,arg,thick, X(2), dmin, kr(3), RR(3), ggp(3), al, &
					lauec(2), g3(3), gdotR,att,xgp,DF_gf(3), &
					DM(2,2), DD, H,FNr(3),ll(3),lpg(3),gplen,LC3, c(3), gx(3), gy(3), &
					sgdenom, gac(3), gbc(3),zmax, beamdiv, ktmax, io_real(2), kt, qx, qy
character(fnlen)      			:: dataname,sgname,voidname,dislname(3*maxdefects),sfname(maxdefects), &
					incname,dispfile,xtalname,foilnmlfile, STEMnmlfile
character(4)            		:: dispmode, progmode
character(12)                         :: fid
complex(kind=dbl),allocatable    	:: DHWM(:,:),DHWMvoid(:,:),DDD(:,:),Sarray(:,:,:,:)
complex(kind=dbl),allocatable    	:: amp(:),amp2(:),Azz(:,:)
complex(kind=dbl)                	:: czero,cone
complex(kind=dbl)                	:: para(0:numdd),dx,dy,dxm,dym
real(kind=sgl),allocatable       	:: inten(:), sgarray(:,:), gmat(:,:,:)
real(kind=sgl),allocatable    		:: disparray(:,:,:,:),imatvals(:,:), ZAimages(:,:,:,:)
integer(kind=irg),allocatable   	:: BFweightsarray(:,:,:),ADFweightsarray(:,:,:)
integer(kind=sgl),allocatable    	:: expval(:,:,:)

namelist / rundata / DF_L, DF_npix, DF_npiy, DF_slice, dmin, sgname, numvoids, incname, numinc, stdout, &
                                voidname, numdisl, dislname, numsf, sfname, dinfo, &
				 t_interval,progmode, dispfile, &
				 dispmode,SETNTHR,xtalname,voltage,kk, lauec, &
				 dataname, foilnmlfile, STEMnmlfile


! first we define the default values
czero=cmplx(0.0,0.0,dbl)
cone=cmplx(1.0,0.0,dbl)

! parameters specific to this run
 xtalname = 'undefined'		! initial value; MUST be present in nml file for program to execute
 voltage = 200000.0			! accelerating voltage
 kk = (/ 0, 0, 1 /)			! incident wave vector in crystal components (omitting wave length)
 lauec = (/ 0.0,0.0 /)			! Laue center coordinates (used for CTEM mode)
 dmin = 0.04			       ! smallest d-spacing to include in dynamical matrix [nm]

! CTEM or STEM ?
 progmode = 'CTEM'  		        ! default illumination mode (can be 'CTEM' or 'STEM')
 STEMnmlfile = 'STEM_rundata.nml'	! name of the STEM rundata namelist file
 foilnmlfile = 'FOIL_rundata.nml'	! name of the foil rundata namelist file
 
! column approximation parameters and image parameters 
 DF_L = 1.0             		! edge length of column in nanometers
 DF_npix = 256       			! number of image pixels along x
 DF_npiy = 256       			! number of image pixels along y 
 DF_slice = 1.0       			! slice thickness in nanometers

 dinfo = 0               		! switch to make makedislocation verbose
 sgname = 'nofile'   			! if this variable is different from 'nofile', then an external sg array is read (to be implemented)

! defect parameters
 numdisl = 0           		! number of dislocation files
 numsf = 0             		! number of stacking fault files
 numinc = 0           			! number of inclusions
 numvoids = 0       			! number of voids
 voidname = 'none' 			! filename for void data
 dislname = ''         		! filenames for dislocation data
 sfname = ''            		! filenames for stacking fault data
 incname = 'none'   			! filename for inclusion data
 dispfile = 'none'     		! name of the displacement field output file (will be created if different from none)
 dispmode = 'not'  			! should a diplacement file be written ('new') or read ('old') or neither ('not')?

! output parameters
 dataname = 'ZAdefect.data'		! default outputfile name
 t_interval = 10       	        ! default timing interval (output every t_interval image columns)
! 
! then we read the actual rundata namelist, which may override some or all of these defaults  
 OPEN(UNIT=dataunit,FILE=nmlfile,DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=rundata)
 CLOSE(UNIT=dataunit)
  
! make sure the xtalname variable has been properly defined
if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMExpeditions:',' structure file name is undefined in '//nmlfile)
end if

! we got this far, so display the standard program info
 progname = 'CTEMExpeditions.f90'
 progdesc = 'Dynamical zone axis CTEM structure matrix computation'
 call CTEMsoft
 
 numd = numdd
 
! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)

! Weickenmeier-Kohl scattering parameters with absorption form factors
 skip = 3	
 call CalcWaveLength(dble(voltage),skip)

 ! generate all atom positions
 call CalcPositions('v')

! determine the point group number and get the ZAP 2-D symmetry  NEEDS TO BE MODIFIED WITH NEW ROUTINES
 j=0
 do i=1,32
  if (SGPG(i).le.cell%SYM_SGnum) j=i
 end do
 call BFsymmetry(kk,j,isym,ir)
  
! determine and display the shortest reciprocal lattice vectors for this zone
 call ShortestG(kk,ga,gb,isym)
 io_int(1:3) = kk(1:3)
 call WriteValue('', io_int, 3,  "(//,' ','[',3I2,'] has Bright Field symmetry ',$)")
 mess = PGTWD(isym)
 call Message("(A,$)")
 io_int(1) = ir
 call WriteValue(' order = ', io_int, 1, "(I4/)")
 mess = 'Reciprocal lattice vectors : '; 
 io_int(1:3) = ga(1:3)
 io_int(4:6) = gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6, "('(',3I3,') and (',3I3,')',/)")

! determine the cartesian components of ga
 DF_gf = float(ga)
 DF_gstar = DF_gf/CalcLength(DF_gf,'r')**2    ! define G* such that G.G* = 1
 call TransSpace(DF_gf,DF_gc,'r','c')                 ! convert to Cartesian reference frame
  
! we'll need to compute the list of wavevectors at this point
  ijmax = 0.0
  call Calckvectors(dble(kk),dble(ga),dble(ktmax),STEM%numberofsvalues,STEM%numberofsvalues,numk,isym,ijmax,'Conical')    ! here we figure out how many beams there are
  call Compute_ReflectionList(dmin,kk,ga,gb,'ALL',.FALSE.,0,0.003)

! for now, we do not consider weak beams at all
  if (BetheParameter%cutoff.eq.0.0) call Set_Bethe_Parameters
  BetheParameter%weakcutoff = BetheParameter%cutoff

! next, we read the foildata namelist from the SRdef_foildata.nml file
! [yes, we're using the same file as for the systematic row case]
! this includes material property data, in this case the elastic moduli,
! and the foil normal, which we will need in the next step
  call read_foil_data(foilnmlfile,DF_npix,DF_npiy,DF_L,dinfo)
  DynFN = foil%F
  
! then we need to prune the reflection list to have only reflections that will actually occur in the 
! computation
  mess = ' Pruning reflection list (this takes a while ...) '
  call Message("(A)")
  call Prune_ReflectionList(numk,nbeams)
  io_int(1) = nbeams
  call WriteValue('Number of contributing beams  : ', io_int, 1, '(I)')
  nn = nbeams
  
! allocate the various Dynamical matrix
  allocate(DynMat(nn,nn))
  DynMat=czero
  
! define the foil thickness, attenuation, and number slices per column
thick = foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
DF_nums = nint(thick/DF_slice)  ! this is the number of slices for this particular column

! is there an inclusion data file? if so, then read it
  if (incname.ne.'none') call read_inclusion_data(numinc,incname,DF_L,DF_npix,DF_npiy,dinfo)

  allocate(disparray(3,DF_nums,DF_npix,DF_npiy))
  disparray = 0.0
  mess = 'displacement field computation '; call Message("(A)")
  
  allocate(DF_R(3,DF_nums))     ! each thread has its own DF_R array 
 
  do i=1,DF_npix  
    do j=1,DF_npiy
      DF_R = 0.0
      call CalcRinclusion(i,j)
      disparray(1:3,1:DF_nums,i,j) = DF_R(1:3,1:DF_nums)
   end do
  end do

! ok, now we have the displacement field, verified in IDL by examining the 
! output in the following file:
!open(unit=dataunit,file=trim(dataname)//'test.data',status='unknown',form='unformatted')
!write(dataunit) shape(disparray)
!write(dataunit) disparray
!close(unit=dataunit,status='keep')
!

! next step is to compute the dynamical matrix without defect phase shifts, which 
! is a simple 7x7 matrix;.
! ir is the row index

gmgp = (/ 0, 0, 0/)
call CalcUcg(gmgp)
xgp = rlp%xgp

rltmpa => reflist%next
ktmp => khead

kr = ktmp%k

! gmat = g-g'
allocate(gmat(3,nn,nn))
gmat = 0.0

write (*,*) 'incident beam direction ', kk
write (*,*) 'Foil normal ', DynFN
write (*,*) 'wave vector ',kr

do ir=1,nn
! ic is the column index
  rltmpb => reflist%next
  do ic=1,nn
! compute the Fourier coefficient of the electrostatic lattice potential 
    if (ic.ne.ir) then  ! not a diagonal entry
      gmgp = rltmpa%hkl - rltmpb%hkl
      gmat(1:3,ir,ic) = float(gmgp(1:3))
      call CalcUcg(gmgp)
      DynMat(ir,ic) = cPi * rlp%qg
    else  ! it is a diagonal entry, so we need the excitation error and the absorption length
      sg = Calcsg(dble(rltmpa%hkl),dble(kr),dble(DynFN))	
write (*,*) rltmpa%hkl, sg
      DynMat(ir,ir) = cmplx(2.D0*cPi*sg,cPi/xgp,dbl)
      write (*,*) DynMat(ir,ir)
    end if
    rltmpb => rltmpb%next
  end do
  rltmpa => rltmpa%next
end do

do i=1,7 
  write (*,*) cabs(DynMat(i,1:7))**2
end do

fid = 'Dynarray.txt'
open(unit=dataunit,file=trim(dataname)//fid,status='unknown',form='formatted')
do i=1,nn
  do j=1,nn
    write (dataunit,"(ES14.4,',',ES14.4)") real(DynMat(i,j)),aimag(DynMat(i,j))
  end do
end do
close(unit=dataunit,status='keep')

! and here is the final computation of 2pi(g-g').R(x,y,z)
  do i=1,DF_npix  
    do j=1,DF_npiy
! for each image column
! open a file for each image column
      write (fid,"(A,I3.3,A,I3.3,A4)") '_',i,'_',j,'.txt'
      open(unit=dataunit,file=trim(dataname)//fid,status='unknown',form='formatted')
      do iz = 1,DF_nums
        write (dataunit,"(I4)") iz 
! and for each slice in a given column
        RR(1:3) = disparray(1:3,iz,i,j)
! we get the displacement vector and compute 2pi(g-g').R for all g-g'
        do ig=1,nn
          do igp=1,nn
            ggp(1:3) = gmat(1:3,ig,igp)
            al = 2.0 * cPi * dot_product(RR,ggp)
            write (dataunit,"(F12.8)") al
          end do
        end do        
      end do
      close(unit=dataunit,status='keep')
    end do
  end do



end subroutine ComputeExpeditions




subroutine CalcRinclusion(i,j)

use local
use constants
use crystal
use crystalvars
use dislocation
use defectmodule
use YSHModule
use foilmodule
use void
use stacking_fault
use inclusion
use quaternions
use rotations

IMPLICIT NONE

integer(kind=irg),INTENT(IN)    	:: i,j
integer(kind=irg)			:: k, islice, ii
real(kind=dbl)        			:: dis,xpos,ypos,zpos,sumR(3),thick,tmp(3),tmp2(3), &
					   tmpf(3),u(3),zaamp,zaphase,zar,zai,zr(3),zi(3), &
                                 	   zt,fx,fy,fz,a_fm(3,3)    !,&
!                                nu,x,y,z,zn,t,pre,r1,r2,r3,th,rn 
                         			 
complex(kind=dbl)     			:: za(3)
complex(kind=sgl)     			:: zero
logical               			:: void

! scale the image coordinates with respect to the origin at the center of the image
 xpos = float(i-DF_npix/2)*DF_L
 ypos = float(j-DF_npiy/2)*DF_L

! determine the starting point of the z-integration for the tilted foil
! this depends on the foil normal components which give the equation
! of the top foil plane as F . r = z0/2, from which we get zt...
 a_fm = qu2om(foil%a_fm)
 fx = a_fm(3,1)
 fy = a_fm(3,2)
 fz = a_fm(3,3) 
 zt = foil%zb*0.5 - (fx*xpos + fy*ypos)/fz
 
! initialize some other variables
 thick = foil%zb
 zero = cmplx(0.0,0.0)

! loop over all slices (this is the main loop)
 sliceloop: do islice = 1,DF_nums 

! zpos is the position down the column, starting at zt (in image coordinates)
    zpos = zt - float(islice)*DF_slice
    
! set the displacements to zero
    sumR = 0.0
        
! set the position in the foil reference frame
    tmpf = quat_rotate_vector( foil%a_fi, dble( (/ xpos, ypos, zpos /)) )

!--------------------
!--SMALL INCLUSIONS--
!--------------------
! then the coherent precipitates, using the model in section 8.4.1
ii = 1
! subtract the inclusion position from the current slice position to get the relative position vector
     tmp = tmpf - (/ inclusions(ii)%xpos, inclusions(ii)%ypos, inclusions(ii)%zpos /)
     dis = CalcLength(tmp,'c')
     if (dis.ge.inclusions(ii)%radius) then ! outside particle
       tmp = tmp*(inclusions(ii)%radius/dis)**3
     end if
     sumR = sumR + inclusions(ii)%C*tmp

   DF_R(1:3,islice) = sumR(1:3)
  end do sliceloop ! main loop over the slices

end subroutine CalcRinclusion


