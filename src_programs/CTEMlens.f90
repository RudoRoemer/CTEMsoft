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
!
! MODULE: lens 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief A few auxiliary variables and functions used by CTEMlens program
! 
!> @date   6/10/98  MDG 1.0 original
!> @date   2/11/99  MDG 1.1 replaced Numerical Recipes routine by rksuite.f
!> @date   5/23/01  MDG 2.0 f90, including rksuite_90 implementation
!> @date  4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
module lens
! defines the field profile and also the derivatives for rksuite_90

use local

! type declaration
type glaser
 real(kind=sgl)    	:: pre2,bmax(2),aa(2),z(2)
 integer(kind=irg) 	:: nof
end type

type (glaser) 		:: GL

contains

!--------------------------------------------------------------------------
!
! FUNCTION: bz 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief axial field value
!
!> @param t position along lens axis
! 
!> @date   6/10/98  MDG 1.0 original
!> @date   2/11/99  MDG 1.1 replaced Numerical Recipes routine by rksuite.f
!> @date   5/23/01  MDG 2.0 f90, including rksuite_90 implementation
!> @date  4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
function bz(t) result(res)

IMPLICIT NONE

real(kind=dbl), intent(in) 		:: t
real(kind=dbl)             		:: res
integer(kind=irg)			:: i

 res = 0.0_dbl
 do i=1,GL % nof
  res = res + GL % bmax(i)/(1.0_dbl+((t-GL % z(i))/GL % aa(i))**2)
 end do

end function bz

!--------------------------------------------------------------------------
!
! FUNCTION: f 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief axial field
!
!> @param t position along lens axis
!> @param y 
! 
!> @date   6/10/98  MDG 1.0 original
!> @date   2/11/99  MDG 1.1 replaced Numerical Recipes routine by rksuite.f
!> @date   5/23/01  MDG 2.0 f90, including rksuite_90 implementation
!> @date  4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
function f(t,y)

real(kind=dbl), intent(in) :: t
real(kind=dbl), dimension(:), intent(in) :: y
real(kind=dbl), dimension(size(y)) :: f

 f(:) = (/ y(2), -GL % pre2*bz(t)**2*y(1) /)
end function f

end module lens




!--------------------------------------------------------------------------
! CTEMsoft2013:CTEMlens.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMlens 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Electron trajectory simulation for simple field
!
!> @details   Simulation of electron trajectory through
!> a lens with a predefined magnetic induction B(z)
!> Output in PostScript.  This incorporates the lens2.f
!> program.
!>
!> Note: internal units for the integration are meters
!> and Tesla.  All quantities are rescaled to a standard
!> range for PostScript output.  
! 
!> @date   6/10/98  MDG 1.0 original
!> @date   2/11/99  MDG 1.1 replaced Numerical Recipes routine by rksuite.f
!> @date   5/23/01  MDG 2.0 f90, including rksuite_90 implementation
!> @date  4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program CTEMlens

use local
use io
use postscript
use constants
use math
use lens
use rksuite
use crystalvars
use crystal
use diffraction
use graphics
! IMPLICIT NONE

integer(kind=irg)            			:: np,p(3),io_int(1)
integer(kind=irg),parameter  		:: nmax = 50
real(kind=sgl)               			:: u(3),v(3),qi(4),qo(4),cc(4),n(3),M(4,4),lw(4),phitot, &
							   foctot,r1,r2,y1,y2,ddz,x1,x2,fysc, io_real(1), oi_real(1)
real(kind=dbl)               			:: ysc
real(kind=dbl),allocatable   		:: r(:),phi(:),bfield(:),rx(:),ry(:),dz(:),db(:),z(:),yy(:)

! rksuite_90 variables
integer(kind=irg)            			:: nout, totf, l
real(kind=dbl)      				:: t_start=0.0_wp, tol=5.0e-5_wp
real(kind=dbl)      				:: t_end, t_got, t_want, pi, t_inc
real(kind=dbl), dimension(2) 		:: y_start, thres = 1.0e-10_wp, y_got, yderiv_got, y_maxvals

type(rk_comm_real_1d) :: comm

        
 progname = 'CTEMlens.f90'
 progdesc = 'Computation of particle trajectories in round magnetic lens'
 call CTEMsoft
 
! generate magnetic field parameters
!
! Note: all internal computations are carried out in 
! units of meters and Tesla
!
! number of sampling points
 call ReadValue(' Number of sampling points : ', io_int, 1)
 np = io_int(1)
 fnp=float(np)
 allocate(r(np))
 allocate(yy(np))
 allocate(phi(np))
 allocate(bfield(np))
 allocate(rx(np))
 allocate(ry(np))
 allocate(dz(np))
 allocate(db(np))
 allocate(z(np))
! how many fields ?
 call ReadValue(' Number of magnetic fields along axis (1 or 2) : ', io_int, 1)
 GL % nof = io_int(1)

 do i=1,GL % nof
   io_int(1) = i
   call WriteValue(' Field # ', io_int, 1, "(I3)")

! maximum strength of B field
   call ReadValue('     Maximum axial field strength [T]   : ', io_real, 1)
   GL % bmax(i) = io_real(1)

! extent of the field in mm
   call ReadValue('    Lateral extent parameter a  [mm]   : ', io_real, 1)
   GL % aa(i) = 0.001D0 * io_real(1)

! location of field maximum along axis
   call ReadValue('    Location of field maximum   [mm]   : ', io_real, 1)
   GL % z(i) = 0.001 * io_real(1)

   mess = '--------'; call Message("(A)")
 end do

! range of z
 call ReadValue('  Total length of optical axis [mm]  : ', io_real, 1)
 zmax = 0.001 * io_real(1)
 zmin=-zmax*0.5
 zmax=zmax*0.5
 zabs=np/2
 ddz=(zmax-zmin)/float(np-1)
! equidistant z array
 do i=1,np
  z(i)=zmin+float(i)*ddz
  bfield(i)=bz(z(i))
 end do 

 1040   format(1x,i3,f12.6)
!
! get acceleration voltage and other constants
 call GetVoltage
!
! initialize a Cartesian reference frame
 cell % xtal_system=1
 cell % a=1.0
 cell % b=1.0
 cell % c=1.0
 cell % alpha=90.0
 cell % beta=90.0
 cell % gamma=90.0
 call CalcMatrices
! solve for phi(z), initial value pi/2
 eta = sngl(dsqrt(cCharge*0.5D0/cRestmass))
 pre1=0.5*eta/sqrt(sngl(mPsihat))
 phi(1)=cPi*0.5
 focal=bfield(1)**2*ddz
 do i=2,np
  phi(i)=phi(i-1)+pre1*bfield(i)*ddz
  focal=focal+bfield(i)**2*ddz
 end do
 focal=1.0/focal/pre1**2
! compute the rotation and focal length at the end of 
! the trajectory, using Simpson''s rule.
 phitot=0.0
 foctot=0.0
 do i=1,np/2-1
  phitot=phitot+2.0*bfield(1+2*i)
  foctot=foctot+2.0*(bfield(1+2*i))**2
 end do
 do i=1,np/2
  phitot=phitot+4.0*bfield(1+(2*i-1))
  foctot=foctot+4.0*(bfield(1+(2*i-1)))**2
 end do
 phitot=pre1*ddz*(phitot+bfield(1)+bfield(np))/3.0
 foctot=ddz*pre1**2*(foctot+bfield(1)**2+bfield(np)**2)/3.0
 foctot=1.0/foctot 
 mess = 'rotation angle [rad]'; call Message("(A)")
 oi_real(1)=phitot
 call WriteValue('  Simpson integration ', oi_real, 1, "(F10.5)")
 oi_real(1)=phi(np)-phi(1)
 call WriteValue('  Naive integration   ', oi_real, 1, "(F10.5)")
 if (GL % nof.eq.1) then
   bb=pre1*GL %bmax(1)*cPi*GL % aa(1) 
   oi_real(1)=bb
   call WriteValue('  Analytical solution ', oi_real, 1, "(F10.5)")
 end if
 mess = 'focal length   [m]'; call Message("(A)")
 oi_real(1)=foctot
 call WriteValue('  Simpson integration ', oi_real, 1, "(F10.5)")
 oi_real(1)=focal 
 call WriteValue('  Naive integration   ',oi_real, 1, "(F10.5)")
 if (GL % nof.eq.1) then 
  bb=pre1**2*GL % bmax(1)**2*cPi*GL % aa(1)*0.5
  oi_real(1)=1.0/bb
  call WriteValue('  Analytical solution ', oi_real, 1, "(F10.5)")
 end if
 mess = 'done solving for phi(z)'; call Message("(A)")
!
! solve for r(z)
!
 GL % pre2=0.25D0*eta**2/mPsihat
 nstep=np-1
 t_start = z(1)
 t_end = z(np)
 y_start = (/ 1.0_wp, 0.0_wp /)
 x2 = t_end
 call setup(comm,t_start,y_start,t_end,tol,thres)
 nout = np; t_inc = (t_end-t_start)/nout
 do l = 1, nout
  t_want = t_end + (l-nout)*t_inc
  call range_integrate(comm,f,t_want,t_got,y_got,yderiv_got)
  yy(l) = y_got(1)
 end do
 call statistics(comm,y_maxvals=y_maxvals,total_f_calls=totf)
 mess = 'done solving for phi(z)'; call Message("(A)")
 io_int(1) = totf
 call WriteValue(' Number of Runge-Kutta integration calls : ', io_int, 1, "(I8)")
!
! and the rest is simply graphics output...
!
! open PostScript outputfile
! and set up PostScript fonts
 PS % psname='lens.ps'
 call PS_openfile(.TRUE.)
!
! generate the page-layout
!
 call PS_newpage(.TRUE.,'Electron Trajectory Computation')
!
! put the text on the top of the page
!
 sp=0.25
 eps=0.15
 call PS_setfont(PSfonts(4),0.14)
! first column
 q=sngl(mAccvol)
 call PS_textvar(sp,PS % psfigheight-0.75,'Voltage [V]',q)
 q=sngl(mPsihat)
 call PS_textvar(sp,PS % psfigheight-0.75-0.18,'Relativistic voltage [V]',q)
 call PS_textvar(sp,PS % psfigheight-0.75-0.36,'Wavelength [pm]',sngl(1000.0*mLambda))
! second column
 call PS_textvar(0.5*PS % psfigwidth,PS % psfigheight-0.75,'Total rotation [rad]',sngl(phi(np)-phi(1)))
 call PS_textvar(0.5*PS % psfigwidth,PS % psfigheight-0.75-0.18,'Focal length  [mm]',1000.0*focal)
!
! define the rectangles
 r3x=sp
 r3y=0.5
 r3w=PS % psfigwidth-2.0*sp
 r3h=0.75*r3w
 r1x=sp
 r1y=2.0*r3y+r3h
 r2w=2.0
 r2h=2.0
 r2x=PS % psfigwidth-sp-r2w
 r2y=r1y
 r1w=PS % psfigwidth-3.0*sp-r2w
 r1h=r2h
!
! Fill rectangle 1
!
 call PS_setlinewidth(0.004)
 call PS_drawrect(r1x,r1y,r1x+r1w,r1y+r1h) 
 call PS_setlinewidth(0.01)
 call PS_line(r1x+eps,r1y+0.5*r1h,r1x+r1w-eps,r1y+0.5*r1h)
 call PS_line(r1x+eps,r1y+eps,r1x+eps,r1y+r1h-eps)
 call PS_setfont(PSfonts(2),0.15)
 call PS_text(r1x+1.2*eps,r1y+r1h-1.3*eps,'B(z), r(z)')
 call PS_setfont(PSfonts(2),0.17)
 call PS_text(r1x+2.0*eps,r1y-0.20,'Radial Plot')
! draw the magnetic field profile
! first normalize to the proper scale
 fysc=0.4*(r1h-2.0*eps)
 ysc=fysc
 call PS_line(r1x+0.7*eps,r1y+0.5*r1h+fysc,r1x+1.3*eps,r1y+0.5*r1h+fysc)
 call PS_line(r1x+0.7*eps,r1y+0.5*r1h-fysc,r1x+1.3*eps,r1y+0.5*r1h-fysc)
 xsc=r1w-2.0*eps
! place the B-max value at the top of the page 
 bmax = maxval(GL % bmax)
 call PS_setfont(PSfonts(4),0.14)
 call PS_textvar(sp,PS % psfigheight-0.75-0.54,'B-max  [T]',bmax)
 dx=r1x+eps
 dy=r1y+0.5*r1h
 do i=1,np
  dz(i)=xsc*(z(i)-zmin)/(zmax-zmin)
  db(i)=ysc*bfield(i)/bmax+dy
 end do 
! label the z-axis
 call PS_setlinewidth(0.01)
 call PS_line(dx-0.2*eps,dy-0.2*eps,dx+0.2*eps,dy+0.2*eps)
 call PS_line(dx-0.2*eps,dy+0.2*eps,dx+0.2*eps,dy-0.2*eps)
 call PS_setfont(PSfonts(2),0.08)
 call PS_textvar(dx-0.4*eps,dy-0.5*eps,' ',1000.0*zmin)
 call PS_line(dx+r1w-2.2*eps,dy-0.2*eps,dx+r1w-1.8*eps,dy+0.2*eps)
 call PS_line(dx+r1w-2.2*eps,dy+0.2*eps,dx+r1w-1.8*eps,dy-0.2*eps)
 call PS_textvar(dx+r1w-4.2*eps,dy-0.5*eps,' ',1000.0*zmax)
 call PS_line(dx+sngl(dz(zabs))-0.2*eps,dy-0.2*eps,dx+sngl(dz(zabs))+0.2*eps,dy+0.2*eps)
 call PS_line(dx+sngl(dz(zabs))-0.2*eps,dy+0.2*eps,dx+sngl(dz(zabs))+0.2*eps,dy-0.2*eps)
 call PS_textvar(dx+sngl(dz(zabs))-0.8*eps,dy-0.5*eps,' ',0.0)
!
 call PS_setlinewidth(0.02)
 do i=1,np-2
  call PS_line(dx+sngl(dz(i)),sngl(db(i)),dx+sngl(dz(i+1)),sngl(db(i+1)))
 end do 
 call PS_setlinewidth(0.03)
 do i=1,np-2
  r1=yy(i)*ysc+dy
  r2=yy(i+1)*ysc+dy
  if (abs(yy(i+1)).le.1.0) then
   call PS_line(dx+sngl(dz(i)),r1,dx+sngl(dz(i+1)),r2)
  end if
 end do
!
! Fill rectangle 2
!
 call PS_setlinewidth(0.004)
 call PS_drawrect(r2x,r2y,r2x+r2w,r2y+r2h) 
 fysc=0.4*(r2h-2.0*eps)
 ysc=fysc
 call PS_setlinewidth(0.01)
 call PS_line(r2x+eps,r2y+0.5*r2h,r2x+r2w-eps,r2y+0.5*r2h)
 call PS_line(r2x+0.5*r2w,r2y+eps,r2x+0.5*r2w,r2y+r2h-eps)
 call PS_setlinewidth(0.007)
 call PS_circle(r2x+0.5*r2w,r2y+0.5*r2h,fysc)
 call PS_setfont(PSfonts(2),0.15)
 call PS_text(r2x+r2w-0.7*eps,r2y+0.5*r2h-0.3*eps,'x')
 call PS_text(r2x+0.49*r2w,r2y+r2h-0.7*eps,'y')
 call PS_setfont(PSfonts(2),0.17)
 call PS_text(r2x+2.0*eps,r2y-0.20,'Azimuthal Plot')
! plot the rotational trajectory
 dx=r2x+0.5*r2w
 dy=r2y+0.5*r2h
 call PS_setlinewidth(0.03)
 do i=1,np-2
  x=ysc*yy(i)*cos(phi(i))+dx
  y=ysc*yy(i)*sin(phi(i))+dy
  rx(i)=x-dx
  ry(i)=y-dy
  r1=ysc*yy(i+1)*cos(phi(i+1))+dx
  r2=ysc*yy(i+1)*sin(phi(i+1))+dy
  if (abs(yy(i+1)).le.1.0) call PS_line(x,y,r1,r2)
 end do
 rx(np)=x2-dx
 ry(np)=y2-dy
!
! Perspective drawing in rectangle 3
! define the clipping region to prevent the trajectories
! from going of the page
!
 call PS_setlinewidth(0.004)
 call PS_closepath
 call PS_move(r3x,r3y)
 call PS_draw(r3x,r3y+r3h)
 call PS_draw(r3x+r3w,r3y+r3h)
 call PS_draw(r3x+r3w,r3y)
 call PS_clippath
 call PS_setfont(PSfonts(2),0.17)
 call PS_text(r3x+2.0*eps,r3y-0.20,'Perspective Plot')
! define some scaling and positioning parameters
 sc=1.6
 dx=r3x+0.15*r3w
 dy=r3y+0.55*r3h
! compute the viewing transformation matrix
! based on pages 405-412 in Computer Graphics: Systems and Concepts
! The vector p(3) defines the viewing point and direction (origin)
 p(1)=4
 p(2)=4
 p(3)=2
 call ComputeViewTrans(p,M,1.0)
! first draw the reference frame
 call get2d(0.0_wp,dz(1),0.0_wp,M,r1,y1,sc)
 call get2d(0.0_wp,dz(np),0.0_wp,M,r2,y2,sc)
 call PS_line(dx+r1,dy+y1,dx+r2,dy+y2)
 call PS_setfont(PSfonts(2),0.15)
 call PS_text(dx+r2+0.1*eps,dy+y2-0.1*eps,'z')
 call get2d(ysc,0.0_wp,0.0_wp,M,r2,y2,sc)
 call PS_line(dx+r1,dy+y1,dx+r2,dy+y2)
 call PS_text(dx+r2-0.1*eps,dy+y2+0.3*eps,'x')
 call get2d(0.0_wp,0.0_wp,ysc,M,r2,y2,sc)
 call PS_line(dx+r1,dy+y1,dx+r2,dy+y2)
 call PS_text(dx+r2-0.4*eps,dy+y2,'y')
! then the beam trajectories
 cc(1)=0.0
 cc(2)=cPi*0.5
 cc(3)=cPi
 cc(4)=cPi*1.5
 lw(1)=0.04
 lw(2)=0.02
 lw(3)=0.02
 lw(4)=0.02
 do j=1,4
  call PS_setlinewidth(lw(j))
  do i=1,np-2
   call get2d(ysc*yy(i)*cos(cc(j)+phi(i)),dz(i),ysc*yy(i)*sin(cc(j)+phi(i)),M,r1,y1,sc)
   call get2d(ysc*yy(i+1)*cos(cc(j)+phi(i+1)),dz(i),ysc*yy(i+1)*sin(cc(j)+phi(i+1)),M,r2,y2,sc)
   call PS_line(dx+r1,dy+y1,dx+r2,dy+y2)
   if (mod(i-1,50).eq.0) then
    call PS_setlinewidth(0.004)
    call get2d(0.0_wp,dz(i),0.0_wp,M,r2,y2,sc)
    call PS_line(dx+r1,dy+y1,dx+r2,dy+y2)
    call PS_setlinewidth(lw(j))
   endif
  end do
 end do
!
 call PS_closefile
 call system(psviewer//'lens.ps')
end program
! ###################################################################
!
!  subroutine get2d
!
!  Author: Marc De Graef
!
!  Description: convert 3D to 2D coordinates for given viewing
!               transformation (normal coordinates)
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/23/01 MDG 2.0 f90
! ###################################################################
subroutine get2d(a,b,c,M,x,y,sc)

use local
use math

real(kind=dbl) :: a,b,c
real(kind=sgl) :: p(4),q(4),M(4,4),x,y,sc
! for the point (a,b,c), compute the 2-D projected
! coordinates using the viewing transformation matrix M
 p(1)=a
 p(2)=b
 p(3)=c  
 p(4)=1.0
 q = matmul(p,M)
 x=sc*q(1)/q(4) 
 y=sc*q(2)/q(4) 
end subroutine
