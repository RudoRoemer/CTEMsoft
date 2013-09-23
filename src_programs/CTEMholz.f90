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
! CTEMsoft2013:CTEMholz.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMholz.f90 
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Generate kinematical HOLZ patterns
!
!> @details  This program generates a multipage PostScript file with HOLZ patterns for the standard lattice parameters
!> and also for slightly modified parameters.   All computations are kinematical only.  A few important global
!> variables are defined in the HOLZvars module.
!
!> @date 1/5/99   MDG 1.0 original
!> @date 4/8/13  MDG 2.0 revision
!--------------------------------------------------------------------------
module HOLZvars

use local

IMPLICIT NONE

type holzreflection
  integer(kind=irg)             		:: hkl(3),n1,n2,N
  logical                       			:: draw,dbdiff
  real(kind=sgl)                		:: hlphi,hlx(2),hly(2),sg,Ig,pxy(2)
  type(holzreflection),pointer        	:: next
end type holzreflection

type(holzreflection),pointer    		:: top,temp,bot
real(kind=sgl)                  		:: g1(3),g2(3),g3(3),H,FNg(3),FNr(3),gshort(3),gp(3),LC1,LC2,thickness,rectangle(2), &
                                   			   PX,PY,thetac,laL,Gmax,Imax,gtoc(2,2),glen,phi,CBEDrad,CBEDsc
integer(kind=irg)               		:: uvw(3),FN(3)

end module HOLZvars


program CTEMholz

use local
use crystalvars
use crystal
use symmetryvars
use symmetry
use graphics
use files
use postscript
use io
use diffraction
use HOLZvars

IMPLICIT NONE

real(kind=sgl)			:: io_real(1)

 progname = 'CTEMholz.f90'
 progdesc = 'Kinematical HOLZ pattern and HOLZ line simulations'
 call CTEMsoft

 SG % SYM_reduce=.TRUE.

! read crystal information, microscope voltage, and camera length
 call CrystalData
 call GetVoltage
 call ReadValue(' Camera length L  [mm, real] ', io_real, 1)
 camlen = io_real(1)

! generate all atom positions in the fundamental unit cell
 call CalcPositions('v')

! open PostScript file
 loadingfile = .FALSE.
 call PS_openfile
 PS%pspage = 0

! generate a set of HOLZ patterns
 call HOLZPage

! close Postscript file
 call PS_closefile

end program

!--------------------------------------------------------------------------
!
! SUBROUTINE: HOLZPage
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief draw zone axis HOLZ diffraction pattern and line pattern
!
! 
!> @date 1/16/02   MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine HOLZPage

use local
use postscript
use crystal
use crystalvars
use symmetry
use symmetryvars
use math
use io
use constants
use diffraction
use dynamical
use HOLZvars

IMPLICIT NONE

type (unitcell)                 		:: savecell

logical                         		:: again, first, newzone, nexttop
real(kind=sgl)                  	:: negative(2),twopi,thr, oi_real(9), &
                                   		  RHOLZmax,RHOLZ(20),xo,yo,sc,pos(2),dy
real(kind=sgl),parameter        	:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/), &
                                   		   yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/), &
                                   		   eps = 1.0E-3
integer(kind=irg)               	:: g1i(3),g2i(3),i,numHOLZ, io_int(1)
real(kind=sgl),parameter        :: le=3.25,he=2.9375

! dimensions of a TEM negative in inches:  3.9375 x 3.1875
 negative = (/ 3.9375, 3.1875 /)
 rectangle = negative*0.5
 SG % SYM_reduce=.TRUE.
 thr = 1.E-4 
 twopi = 2.0*cPi
 Imax = 0.0

! camera length
 laL = sngl(mLambda) * camlen
 oi_real(1) = sngl(mLambda)
 call WriteValue(' wavelength [nm] = ', oi_real, 1, "(F10.6)")
 oi_real(1) = camlen
 call WriteValue(' L         [mm] = ', oi_real, 1, "(f10.2)")
 oi_real(1) = laL
 call WriteValue(' camera length lambda*L [mm nm] = ', oi_real, 1, "(f10.5)")

! what portion of reciprocal space is covered by the negative along horizontal direction ?
 Gmax = sqrt(rectangle(1)**2+rectangle(2)**2)*25.4/laL

! this is also the maximum allowable HOLZ radius
 RHOLZmax = Gmax

! next loop over multiple zone axis orientations or different lattice parameters
 again = .TRUE.
 first = .TRUE.
 do while (again)

! new zone axis or modify lattice parameters ?
  if (first.eqv..TRUE.) then
   newzone = .TRUE.

! get the zone axis
   call GetIndex(uvw,'d')
   call ReadValue(' Enter foil thickness [nm, R] (to get relrods of proper length) : ', oi_real, 1)
   thickness = oi_real(1)
   first = .FALSE.

! allocate the linked list to store all reflections
   if (.not.associated(top)) then 
     allocate(top)
     bot => top
     nullify(bot%next)
   end if
  else
! either get a new zone or change the lattice parameters and keep the zone
   call ReadValue(' New zone (1) or change lattice parameters for present zone (2) ', io_int, 1)
   if (io_int(1).eq.1) then
    newzone = .TRUE.
    call GetIndex(uvw,'d')
    cell = savecell
! deallocate the previous linked list and allocate a new one
    temp => top%next
    do while (associated(temp%next))
     deallocate(top)
     top => temp
     temp => top%next
    end do
    deallocate(top)
    allocate(top)
    bot => top
    nullify(bot%next)
   else
! show the current lattice parameters and save the parameters
     newzone = .FALSE.
   end if
  end if

! it is a new zone, so draw the diffraction pattern on the top half of the page
  if (newzone.eqv..TRUE.) then
! get the basis vectors g1 and g2
    call ShortestG(uvw,g1i,g2i,i)
    g1 = float(g1i); g2 = float(g2i)
    
! get the beam divergence angle to determine the diameter of the central disk
    oi_real(1) = 500.0*minval( (/ CalcDiffAngle(g1i(1),g1i(2),g1i(3)), &
                               CalcDiffAngle(g2i(1),g2i(2),g2i(3)) /) )
    call WriteValue(' Maximum disk diameter without overlap [mrad]= ', oi_real, 1, "(f10.4)")
    call ReadValue(' Enter the beam divergence angle [mrad, R] ', oi_real, 1); 
    thetac = oi_real(1)*0.001
    
! distance between consecutive HOLZ layers in nm-1
    H = 1.0/CalcLength(float(uvw),'d')

! determine g3 basis vector
    call CalcCross(g1,g2,g3,'r','r',1)
    call NormVec(g3,'r')
    g3 = H * g3

! get foil normal
    mess = 'Enter Foil Normal F [real space indices]'; call Message("(A)")
    call GetIndex(FN,'d')
! compute components of FN with respect to g1, g2, g3
    call TransSpace(float(FN),FNr,'d','r')
    call NormVec(FNr,'r')
    FNg = (/ CalcDot(FNr,g1,'r'), CalcDot(FNr,g2,'r'), CalcDot(FNr,g3,'r') /)

! determine shortest vector of FOLZ layer
    call ShortestGFOLZ
    oi_real(1:3) = gp(1)*g1(1:3)+gp(2)*g2(1:3)
    call WriteValue(' HOLZ shift vector = ', oi_real, 3, "(3f9.4)") 

! get Laue center in terms of g1 and g2
    oi_real(1:3) = g1(1:3)
    oi_real(4:6) = g2(1:3)
    oi_real(7:9) = g3(1:3)
    call WriteValue( 'The new basis vectors for this zone axis are ', &
          oi_real, 9, "(/'g1 = ',3f10.5,/'g2 = ',3f10.5,/'g3 = ',3f10.5,/)")
    oi_real(1) = H
    call WriteValue(' reciprocal interplanar spacing H = ', oi_real, 1, "(F10.4,' nm^-1'/)")
    call ReadValue(' Enter the coordinates of the Laue center with respect to g1 and g2', oi_real, 2)
    LC1 = oi_real(1)
    LC2 = oi_real(2)

! compute how many HOLZ layers need to be drawn.
! this follows from the camera length and the size of the micrograph
    i=1
    numHOLZ = 0
    do while(i.lt.20)
     RHOLZ(i) = sqrt(2.0*H*float(i)/mLambda - (float(i)*H)**2)    
     if (RHOLZ(i).lt.RHOLZmax) numHOLZ = numHOLZ+1
     i=i+1
    end do
    
! print the  number and radii of the HOLZ rings in the field of view
    mess = 'Number and radii of possible HOLZ rings inside field of view'; call Message("(A)")
    oi_real(1) = RHOLZmax
    call WriteValue(' RHOLZmax = ', oi_real, 1, "(F10.5)")
    do i=1,numHOLZ
      oi_real(1)=float(i); oi_real(2)=RHOLZ(i)
      call WriteValue(' Ring ', oi_real, 2, "(F5.0,3x,F10.5)")
    end do

! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
    if (newzone) then
     call PS_newpage(.FALSE.,'Kinematical HOLZ Diffraction Patterns')
     call PS_text(5.25,-0.05,'scale bar in reciprocal nm')
     call PS_textvar(5.25,PS % psfigheight+0.02,'Camera Constant [nm mm]',laL)
     call PS_setfont(PSfonts(2),0.15)
     call PS_text(-0.25,PS % psfigheight+0.02,'Structure File : '//cell % fname)
    end if

! draw frame and related stuff
    xo = 2.25
    yo = 5.00
    PX = xo + rectangle(1)
    PY = yo + rectangle(2)
    CBEDrad = 1.5
    CBEDsc = 1.3
    call PS_setlinewidth(0.012)
    call PS_balloon(xo,yo,negative(1),negative(2),0.0312)

! zone axis
    call PS_setfont(PSfonts(2),0.12)
    call PS_text(xo+0.05,yo+negative(2)+0.12,'Zone axis ')
    call PrintIndices('d',uvw(1),uvw(2),uvw(3),xo+0.6,yo+negative(2)+0.12)

! add other data lines to the upper left
    call PS_setfont(PSfonts(2),0.15)
    call PS_textvar(-0.25,PS%psfigheight-0.18,'Acc. Voltage [kV] ',sngl(mAccvol)*0.001)
    call PS_text(-0.25,PS%psfigheight-0.38,'Foil normal ')
    call PrintIndices('d',FN(1),FN(2),FN(3),-0.25+1.5,PS%psfigheight-0.38)
    call PS_textvar(-0.25,PS%psfigheight-0.58,'Foil thickness [nm] ',thickness)
    call PS_text(-0.25,PS%psfigheight-0.78,'Laue center ')
    call PS_textvar(-0.25+1.1,PS%psfigheight-0.78,'',LC1)
    call PS_textvar(-0.25+1.6,PS%psfigheight-0.78,'',LC2)

! HOLZ ring radii
    call PS_text(xo-1.5,PS%psfigheight-1.45,'HOLZ radii [nm-1] ')
    do i=1,numHOLZ
        call PS_textint(xo-1.5,PS%psfigheight-1.5-float(i)*0.14,'',i)
        call PS_textvar(xo-1.3,PS%psfigheight-1.5-float(i)*0.14,'',RHOLZ(i))
    end do

! CBED 000 disk text
    call PS_setfont(PSfonts(2),0.12)
    call PS_textvar(-0.25,0.5,'Convergence angle [mrad] ',thetac*1000.0)

! lattice parameters
    call PS_setfont(PSfonts(4),0.14)
    call PS_text(-0.25,2.00,'a :')
    call PS_text(-0.25,1.84,'b :') 
    call PS_text(-0.25,1.68,'c :')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,2.00,cell % a
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.84,cell % b
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.68,cell % c
    call PS_setfont(PSfonts(1),0.14)
    call PS_text(-0.25,1.52,'a :')
    call PS_text(-0.25,1.36,'b :')
    call PS_text(-0.25,1.20,'g :')
    call PS_setfont(PSfonts(4),0.14)
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.52,cell % alpha
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.36,cell % beta 
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.20,cell % gamma

! scale bar (sc is the conversion factor from nm-1 to inches)
    sc = laL/25.4
    call PS_setlinewidth(0.020)
    call PS_line(xo+0.05,yo+0.06,xo+0.05+5.0*sc,yo+0.06)
    call PS_setfont(PSfonts(2),0.15)
    call PS_text(xo+0.05+2.5*sc,yo+0.10,'5 ')

! plot origin of reciprocal space 
    call PS_filledcircle(PX,PY,0.03,0.0)

! set clip path
    call PS_closepath
    call PS_gsave
    call PS_move(xo,yo)
    call PS_draw(xo,yo+negative(2))
    call PS_draw(xo+negative(1),yo+negative(2))
    call PS_draw(xo+negative(1),yo)
    call PS_clippath
    call PS_newpath

! now we are ready to compute the HOLZ layers and store them in the linked list
! first the ZOLZ layer
    call CalcHOLZ(0)

! then the HOLZ layers 
    do i=1,numHOLZ
      call PS_setlinewidth(0.005)
      call PS_circle(PX,PY,RHOLZ(i)*laL/25.4)
      call CalcHOLZ(i)
    end do

! once that is done, we can determine the intensities to be drawn and 
! draw all reflections.
    call PlotHOLZ

! and eliminate the current clippath
    call PS_closepath
    call PS_grestore

! draw the vectors g1 and g2, and the projection gp of G
    call PS_setlinewidth(0.02)
    call PS_filledcircle(0.5,PY-2.5,0.015,0.0)

! g1
    pos = matmul(gtoc,(/1.0,0.0/) )
    glen = 2.0*sqrt(pos(1)**2+pos(2)**2)
    pos = pos/glen
    call PS_line(0.5,PY-2.5,0.5+pos(1),PY-2.5+pos(2))
    call PS_filledcircle(0.5+pos(1),PY-2.5+pos(2),0.04,0.0)
    call PrintIndices('r',int(g1(1)),int(g1(2)),int(g1(3)),0.5+pos(1)+0.1,PY-2.5+pos(2))

! g2
    pos = matmul(gtoc,(/0.0,1.0/) )
    pos = pos/glen
    call PS_line(0.5,PY-2.5,0.5+pos(1),PY-2.5+pos(2))
    call PS_filledcircle(0.5+pos(1),PY-2.5+pos(2),0.04,0.0)
    call PrintIndices('r',int(g2(1)),int(g2(2)),int(g2(3)),0.5+pos(1)+0.1,PY-2.5+pos(2))

! and then the projection of G onto g1,g2
    pos = matmul(gtoc,(/gp(1),gp(2)/) )
    pos = pos/glen
    call PS_setlinewidth(0.02)
    call PS_line(0.45+pos(1),PY-2.5+pos(2),0.55+pos(1),PY-2.5+pos(2))
    call PS_line(0.5+pos(1),PY-2.45+pos(2),0.5+pos(1),PY-2.55+pos(2))
    call PS_text(-0.5,3.2,'Basis vectors g1, g2,')
    call PS_text(-0.5,3.0,'and projection of G (cross)') 
! draw fixed radius circle for bright field CBED disk  
    call PS_setlinewidth(0.025)
    call PS_circle(PX,yo-2.5,CBEDrad)
! indicate center of pattern
    call PS_setlinewidth(0.01)
    call PS_line(PX-0.05,yo-2.5,PX+0.05,yo-2.5)
    call PS_line(PX,yo-2.45,PX,yo-2.55)

    call PlotHOLZlines(0.0)
    nexttop = .TRUE.
  else  ! this is not a new zone axis
! let the user define new lattice parameters 
    savecell = cell
    oi_real(1) = cell%a; oi_real(2) = cell%b; oi_real(3) = cell%c
    call WriteValue(' Current lattice parameters [nm] ', oi_real, 3, "(/'a = ',f7.5,', b = ',f7.5,', c = ',f7.5)")
    oi_real(1) = cell%alpha; oi_real(2) = cell%beta; oi_real(3) = cell%gamma
    call WriteValue('', oi_real, 3, "(/'alpha = ',f7.2,', beta = ',f7.2,', gamma = ',f7.2)")

! ask for the new parameters (all must be entered) and recompute metric information
    call ReadValue(' Enter new lattice parameters a, b, and c [nm] ', oi_real, 3)
    cell%a = oi_real(1); cell%b = oi_real(2); cell%c = oi_real(3)
   call ReadValue(' Enter new angles alpha, beta, and gamma [degrees] ', oi_real, 3)
    cell%alpha = oi_real(1); cell%beta = oi_real(2); cell%gamma = oi_real(3)
    call CalcMatrices

! redo the geometry
    call ReCalcHOLZ

! move to top or bottom for next drawing ?
    if (nexttop.eqv..TRUE.) then
     call PS_newpage(.FALSE.,'Kinematical HOLZ Diffraction Patterns')
     dy = 4.25
     nexttop=.FALSE.
    else
     dy = -0.25
     nexttop=.TRUE.
    end if

! CBED 000 disk text
    call PS_setfont(PSfonts(2),0.12)
    call PS_textvar(-0.25,0.5+dy,'Convergence angle [mrad] ',thetac*1000.0)

! lattice parameters
    call PS_setfont(PSfonts(4),0.14)
    call PS_text(-0.25,2.00+dy,'a :')
    call PS_text(-0.25,1.84+dy,'b :') 
    call PS_text(-0.25,1.68+dy,'c :')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,2.00+dy,cell % a
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.84+dy,cell % b
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.68+dy,cell % c
    call PS_setfont(PSfonts(1),0.14)
    call PS_text(-0.25,1.52+dy,'a :')
    call PS_text(-0.25,1.36+dy,'b :')
    call PS_text(-0.25,1.20+dy,'g :')
    call PS_setfont(PSfonts(4),0.14)
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.52+dy,cell % alpha
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.36+dy,cell % beta 
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.20+dy,cell % gamma

! draw fixed radius circle for bright field CBED disk  
    call PS_setlinewidth(0.025)
    call PS_circle(PX,yo-2.5+dy,CBEDrad)

! indicate center of pattern
    call PS_setlinewidth(0.01)
    call PS_line(PX-0.05,yo-2.5+dy,PX+0.05,yo-2.5+dy)
    call PS_line(PX,yo-2.45+dy,PX,yo-2.55+dy)

    call PlotHOLZlines(dy)
  end if ! newzone .eq. .TRUE.

  call ReadValue(' Another pattern ? [1/0] ', io_int, 1)
  if (io_int(1).ne.1) again=.FALSE.
 end do  ! end of main loop

end subroutine

!--------------------------------------------------------------------------
! 
! SUBROUTINE:ShortestGFOLZ
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief find the shortest G vector
!
!> @details  see chapter 3
!
! 
!> @date 1/29/02   MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine ShortestGFOLZ

use local
use io
use postscript
use crystal
use crystalvars
use error
use HOLZvars

IMPLICIT NONE

real(kind=sgl)               			:: gmin,gam11,gam12,gam22
integer(kind=irg),parameter  		:: inm = 8
integer(kind=irg)            			:: ih,ik,il,NN, oi_int(1)

! look for the shortest reflection satisfying hu+kv+lw = 1
! This could be replaced by code from Jackson's paper (1987),
! but it does essentially the same thing.
 gmin = 100.0
 NN=1
 do while((gmin.eq.100.0).and.(NN.lt.4))
  do ih=-inm,inm
   do ik=-inm,inm
    do il=-inm,inm
! does this reflection lie in the plane NN ?
     if ((ih*uvw(1)+ik*uvw(2)+il*uvw(3)).eq.NN) then
      glen = CalcLength(float((/ih,ik,il/)),'r')
      if (glen.lt.gmin) then
       gmin = glen
       gshort = float( (/ ih,ik,il /) )
      end if
     end if
    end do
   end do
  end do
  oi_int(1) = NN
  call WriteValue(' Could not find any reflections with hu+kv+lw = ', oi_int, 1, "(I2)")
  NN = NN+1
 end do
 if (gmin.eq.100.0) then ! for some reason there is no reflection with N<=3 ...
  call FatalError('ShortestGFOLZ: ',' could not find any reflections with hu+kv+lw<=3 ...')
 end if

! projected components of G
 gam11 = CalcDot(g1,g1,'r')
 gam12 = CalcDot(g1,g2,'r')
 gam22 = CalcDot(g2,g2,'r')
 gmin = 1.0/(gam11*gam22-gam12**2)
 gp(1) = (CalcDot(gshort,g1,'r')*gam22-CalcDot(gshort,g2,'r')*gam12)*gmin
 gp(2) = (CalcDot(gshort,g2,'r')*gam11-CalcDot(gshort,g1,'r')*gam12)*gmin

end subroutine ShortestGFOLZ


!--------------------------------------------------------------------------
! 
! SUBROUTINE:CalcHOLZ
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute the HOLZ reflections for zone N 
!
!> @param N HOLZ number
!
! 
!> @date 1/29/02   MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine CalcHOLZ(N)

use local
use io
use postscript
use crystal
use crystalvars
use symmetry
use error
use diffraction
use dynamical
use constants
use HOLZvars

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: N 

integer(kind=irg)                	:: inmhkl(2),hc,i,j,nref,istat, oi_int(1)
real(kind=sgl)                   	:: correction,gg(3),Ig,smax,kk(3),gxy(2),pxy(2),exer,sgdenom,x,tgm,qx,qy,y,det,LC3, &
                                    		   ll(3),lpg(3),gplen
logical                          		:: a,dbdiff
character(1)                     	:: q

 mess = 'Computing HOLZ reflection data'; call Message("(/A)")

! set the index boundaries
 inmhkl(1) = int(1.1*Gmax/CalcLength(g1,'r'))
 inmhkl(2) = int(1.1*Gmax/CalcLength(g2,'r'))

! we will take g1 to be the x-axis of the pattern
! so the transformation matrix from g1,g2 to Cartesian is given by
 if (N.eq.0) then
  phi = CalcAngle(g1,g2,'r')
  glen = CalcLength(g2,'r')
  gtoc(1,1) = CalcLength(g1,'r')
  gtoc(1,2) = glen*cos(phi)
  gtoc(2,1) = 0.0
  gtoc(2,2) = glen*sin(phi)
  gtoc = gtoc*laL/25.4  ! nm^-1 to inches via the camera length
end if

! loop over all possible reflections in this layer (2D loop !!!)
 smax = 20.0/thickness  ! the number 10 is arbitrary and could be changed at will
 nref=0
 do i=-inmhkl(1),inmhkl(1)
  do j=-inmhkl(2),inmhkl(2)
! make sure this reflection is close enough to the Ewald sphere to give some intensity;
! use crystal thickness to determine this according to the argument of the sinc function.
      gg = i*g1 + j*g2 + N*gshort

! Do this only for those reflections that are allowed by
! the lattice centering !
      a = IsGAllowed(int(gg))
      if (a) then

! compute excitation error, including Laue center, foil normal, and HOLZ reflection.
       glen = CalcLength(gg,'r')
       if (glen.ne.0.0) then
         ll = LC1*g1 + LC2*g2
         lpg = ll + gg
         gplen = CalcLength(lpg,'r')
         LC3 = sqrt(1.0-mLambda**2*CalcLength(ll,'r')**2)
	 if (gplen.eq.0.0) then
           exer = -mLambda*CalcDot(gg,2.0*ll+gg,'r')/(2.0*LC3*cos(CalcAngle(float(uvw),float(FN),'d'))	 )
	 else
	   sgdenom = 2.0*LC3*cos(CalcAngle(float(uvw),float(FN),'d'))- &
	           2.0*mLambda*gplen*cos(CalcAngle(lpg,FNr,'r'))
           exer = -(mLambda*CalcDot(gg,2.0*ll+gg,'r')-2.0*LC3*gplen*cos(CalcAngle(g3,lpg,'r')))/sgdenom
	 end if
       else
         exer = 10000.0
       end if

! exclude the 000 reflection
       if (abs(exer).le.smax) then

! OK, it is close enough.  Does it have any intensity ?
        call CalcUcg(int(gg))

! yes, it does.  get the scaled intensity using the sinc function
        Ig = rlp%Vmod**2
	if (Ig.lt.1.0e-16) then 
	  dbdiff = .TRUE.
	else
	  dbdiff = .FALSE.
	end if
        if (abs(exer).ge.1.0e-5) then
	 Ig = Ig * (sin(cPi*exer*thickness)/(cPi*exer*thickness))**2
        end if 

! store maximum intensity
        if (Ig.gt.Imax) Imax = Ig

! next, determine the drawing coordinates, first in terms of g1 and g2
        correction = 1.0/(1.0-mLambda*H*(float(N)+exer*FNg(3)))
        gxy = (/ (i+N*gp(1)+exer*FNg(1)), (j+N*gp(2)+exer*FNg(2))  /) * correction

! convert to Cartesian drawing coordinates
        pxy = matmul(gtoc,gxy)

! and add the point to the linked list 
        allocate(bot%next,stat=istat)
        if (istat.ne.0) then
	  call PlotHOLZ
	  call PS_closefile
	  call FatalError('CalcHOLZ: ',' unable to allocate memory for linked list')
	end if
        bot => bot%next
        nullify(bot%next)
        bot%hkl = gg
        bot%n1  = i
        bot%n2  = j
        bot%N   = N
        bot%sg  = exer
        bot%Ig  = Ig
        bot%pxy = pxy
	bot%dbdiff = dbdiff
	nref = nref+1

! would this point contribute to the HOLZ line drawing in the central disk ?
        phi = asin(mLambda*glen*0.5) - asin(N*H/glen)
        bot%draw = .FALSE.
	if (abs(phi).le.thetac) then
         x = phi/thetac  *  CBEDrad
	 if (pxy(1).ne.0.0) then
  	   tgm = pxy(2)/pxy(1)
	   y = atan2(pxy(2),pxy(1))
	   qx = x*cos(y)
	   qy = x*sin(y)
	   det = 1.0-(1.0+tgm**2)*(1.0-(tgm*CBEDrad*CBEDsc/(qx+tgm*qy))**2)
	   if (det.gt.0.0) then  ! there is an intersection for this line so it should be drawn
	     bot%draw = .TRUE.
	     bot%hlx(1) = (qx+tgm*qy)*(1.0-sqrt(det))/(1.0+tgm**2)
	     bot%hly(1) = qy-(bot%hlx(1)-qx)/tgm
	     bot%hlx(2) = (qx+tgm*qy)*(1.0+sqrt(det))/(1.0+tgm**2)
	     bot%hly(2) = qy-(bot%hlx(2)-qx)/tgm
	   end if
	 else  ! parallel to the y-axis (easy to deal with)
	     bot%draw = .TRUE.
	     bot%hlx(1) = qx
	     bot%hly(1) = sqrt((CBEDrad*CBEDsc)**2-qx**2)
	     bot%hlx(2) = qx
	     bot%hly(2) = -bot%hly(1)	   
	 end if
	end if

       end if
    end if
   end do
  end do

  oi_int(1) = nref
  Call WriteValue(' number of reflections to be drawn : ', oi_int, 1, "(I6)")

end subroutine CalcHOLZ

!--------------------------------------------------------------------------
! 
! SUBROUTINE:ReCalcHOLZ
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute the HOLZ reflections for zone N for new lattice parameters
!
!> @date 2/03/02   MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine ReCalcHOLZ

use local
use io
use postscript
use crystal
use crystalvars
use symmetry
use error
use diffraction
use dynamical
use constants
use HOLZvars

IMPLICIT NONE


real(kind=sgl)                   :: correction,gg(3),gxy(2),pxy(2),exer,sgdenom,x,tgm,qx,qy, &
                                    	     y,det,LC3,ll(3),lpg(3),gplen


 mess = 'Computing HOLZ reflection data'; call Message("(/A/)")

    temp => top%next
    do while (associated(temp))
       gg = temp%hkl

! compute excitation error
       glen = CalcLength(gg,'r')
       ll = LC1*g1 + LC2*g2
       lpg = ll + gg
       gplen = CalcLength(lpg,'r')
       LC3 = sqrt(1.0-mLambda**2*CalcLength(ll,'r')**2)
       if (gplen.eq.0.0) then
         exer = -mLambda*CalcDot(gg,2.0*ll+gg,'r')/2.0*LC3*cos(CalcAngle(float(uvw),float(FN),'d'))	 
       else
         sgdenom = 2.0*LC3*cos(CalcAngle(float(uvw),float(FN),'d'))- &
	         2.0*mLambda*gplen*cos(CalcAngle(lpg,FNr,'r'))
         exer = -(mLambda*CalcDot(gg,2.0*ll+gg,'r')-2.0*LC3*gplen*cos(CalcAngle(g3,lpg,'r')))/sgdenom
       end if

! next, determine the drawing coordinates, first in terms of g1 and g2
       correction = 1.0/(1.0-mLambda*H*(float(temp%N)+exer*FNg(3)))
       gxy = (/ (temp%n1+temp%N*gp(1)+exer*FNg(1)), (temp%n2+temp%N*gp(2)+exer*FNg(2))  /) * correction

! convert to Cartesian drawing coordinates
       pxy = matmul(gtoc,gxy)
       phi = asin(mLambda*glen*0.5) - asin(temp%N*H/glen)
       temp % draw = .FALSE.
       if (abs(phi).le.thetac) then
        x = phi/thetac  *  CBEDrad
        if (pxy(1).ne.0.0) then
         tgm = pxy(2)/pxy(1)
	 y = atan2(pxy(2),pxy(1))
	 qx = x*cos(y)
         qy = x*sin(y)
	 det = 1.0-(1.0+tgm**2)*(1.0-(tgm*CBEDrad*CBEDsc/(qx+tgm*qy))**2)
	 if (det.gt.0.0) then  ! there is an intersection for this line so it should be drawn
	  temp%draw = .TRUE.
          temp%hlx(1) = (qx+tgm*qy)*(1.0-sqrt(det))/(1.0+tgm**2)
	  temp%hly(1) = qy-(temp%hlx(1)-qx)/tgm
	  temp%hlx(2) = (qx+tgm*qy)*(1.0+sqrt(det))/(1.0+tgm**2)
	  temp%hly(2) = qy-(temp%hlx(2)-qx)/tgm
         end if
	else  ! parallel to the y-axis (easy to deal with)
	  temp%draw = .TRUE.
	  temp%hlx(1) = qx
	  temp%hly(1) = sqrt((CBEDrad*CBEDsc)**2-qx**2)
	  temp%hlx(2) = qx
	  temp%hly(2) = -temp%hly(1)	   
	end if
       end if

! move to the next reflection
       temp=>temp%next
     end do

end subroutine ReCalcHOLZ

!--------------------------------------------------------------------------
! 
! SUBROUTINE:PlotHOLZ
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief draw a single HOLZ zone axis diffraction pattern
!
!> @date 1/22/02   MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine PlotHOLZ

use local
use io
use postscript
use HOLZvars

IMPLICIT NONE

real(kind=sgl)         :: V,qx,qy

 mess = 'Plotting HOLZ reflections'; call Message("(/A/)")

! point to the top of the linked list
 temp => top%next

! move through the entire list and draw all reflections with exponentially scaled intensity;
 do while (associated(temp))
  V=0.05*(temp%Ig/Imax)**0.1
  qx = temp%pxy(1)
  qy = temp%pxy(2)

! make sure the point is inside the rectangle
  if ((abs(qx).lt.rectangle(1)).and.(abs(qy).lt.rectangle(2))) then
   if (temp%dbdiff) then  ! potential double diffraction spot
    call PS_setlinewidth(0.005)
    call PS_square(PX+qx,PY+qy,0.04)
   else ! regular reflection
    call PS_filledcircle(PX+qx,PY+qy,V,0.0)
   end if
  end if

! move to the next reflection
  temp=>temp%next
 end do

end subroutine PlotHOLZ


!--------------------------------------------------------------------------
! 
! SUBROUTINE:PlotHOLZlines
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief draw a single HOLZ line convergent beam disk
!
!> @param dy offset parameter
!
!> @date 1/22/02   MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------

subroutine PlotHOLZlines(dy)

use local
use io
use postscript
use constants
use HOLZvars

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: dy

real(kind=sgl)         			:: V,V2,qx,qy,CB
character(12)          			:: txt
integer(kind=irg)  			:: i,nref

 mess = 'Plotting HOLZ lines and labels'; call Message("(/A/)")

! point to the top of the linked list
 temp => top%next
 nref = 0
 open(unit=30,file='temp.txt',status='unknown',form='formatted')
 call PS_setfont(PSfonts(4),0.08)

! move through the entire list and draw the corresponding HOLZ line across the 000 diffraction disk
 do while (associated(temp))
 
  if (temp%draw.eqv..TRUE.) then ! draw the HOLZ line on the second drawing 
   V=0.075*(temp%Ig/Imax)**0.1
   V2=1.0-4.0*(temp%Ig/Imax)**0.1
   call PS_setlinewidth(V)

! make sure that all the lines will actually fit inside the region of interest
   if ((maxval(abs(temp%hlx)).le.CBEDrad*CBEDsc).and.(maxval(abs(temp%hly)).le.CBEDrad*CBEDsc)) then
    call PS_line_gray(PX+temp%hlx(1),2.5+temp%hly(1)+dy,PX+temp%hlx(2),2.5+temp%hly(2)+dy,V2)

! add indices along continuation of lines
    qx= PX+1.02*temp%hlx(2)
    qy= 2.5+1.02*temp%hly(2)+dy
    V = 180.0*atan2(temp%hly(2)-temp%hly(1),temp%hlx(2)-temp%hlx(1))/cPi
    write (30,"(1x,I3,1x,I3,1x,I3,1x,3(f10.5,1x))") (temp%hkl(i),i=1,3),qx,qy,V
!    write (*,"(1x,I3,1x,I3,1x,I3,1x,3(f10.5,1x))") (temp%hkl(i),i=1,3),qx,qy,V
    nref = nref+1
   end if
  end if

! move to the next reflection
  temp=>temp%next
 end do

! add indices along continuation of lines
 close(unit=30,status='keep')

 CB = CBEDrad**2
 open(unit=30,file='temp.txt',status='old',form='formatted')
 do i=1,nref
   read (30,"(A12,1x,3(f10.5,1x))") txt,qx,qy,V
   if (((qx-PX)**2+(qy-(2.5+dy))**2).gt.CB) then
    call PS_move(qx,qy)  ! just outside clipping ring
    write (psunit,*) V,' rotate'
    write (psunit,"(1x,'( ',A12,' ) show')") txt
    write (psunit,*) -V,' rotate'
   end if
 end do
 close(unit=30,status='delete')

end subroutine PlotHOLZlines
