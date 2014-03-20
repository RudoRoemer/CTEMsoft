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
! CTEMsoft2013:CTEMSRBW.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMSRBW 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Systematic row Bloch wave program
!
!> @todo program needs to be debugged thoroughly !
! 
!> @date   4/24/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!> @date  03/04/14 PGC 3.0.1 gfortran compatibility stuff
!--------------------------------------------------------------------------
program CTEMSRBW

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

IMPLICIT NONE

integer(kind=irg)          		:: ira,ns,g(3),k(3),fn(3),io_int(1)
real(kind=sgl)             		:: ktmax,io_real(1)

 progname = 'CTEMSRBW.f90'
 progdesc = 'Systematic row convergent beam pattern (Bloch waves)'
 call CTEMsoft

! first get the crystal data and microscope voltage
 SG % SYM_reduce=.TRUE.
 call CrystalData
 call GetVoltage

! generate all atom positions
 call CalcPositions('v')

! get the reciprocal lattice vector
 mess = 'Diffraction vector :'; call Message("(A)")
 call GetIndex(g,'r')

! ask the user for the beam direction 
 mess = 'Beam direction :'; call Message("(A)")
 call GetIndex(k,'d')
 mess = 'Foil normal    :'; call Message("(A)")
 call GetIndex(fn,'d')

! some more parameters
 call ReadValue(' Maximum multiple of g to consider : ', io_int, 1)
 ira=io_int(1)
 call ReadValue(' Enter maximum value of k_t in units of g: ', io_real, 1)
 ktmax = io_real(1)
 call ReadValue(' Number of orientations: ', io_int, 1)
 ns=io_int(1)

! do the computation
 call CalcSRBW(g,float(k),float(fn),ira,ns,ktmax)

end program CTEMSRBW

!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcSRBW 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief performs the actual Bloch wave calculation
! 
!> @param g reciprocal lattice vector
!> @param k wave vector
!> @param f foil normal
!> @param ira index range
!> @param ns number of orientations
!> @param ktmax maximum tangential wave vector component
!
!> @date   4/24/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine CalcSRBW(g,k,f,ira,ns,ktmax)

use local
use io
use crystal
use crystalvars
use diffraction
use constants
use dynamical

IMPLICIT NONE

integer(kind=irg),INTENT(IN)		:: g(3)
real(kind=sgl),INTENT(INOUT)		:: k(3)
real(kind=sgl),INTENT(IN)		:: f(3)
integer(kind=irg),INTENT(IN)		:: ira
integer(kind=irg),INTENT(IN)		:: ns
real(kind=sgl),INTENT(IN)		:: ktmax
	
real(kind=sgl)             			:: pre,pre2,glen,upzero,gg,kk,kt(3),kz,kn,kttb,dkt,frac,exer, io_real(1)
integer(kind=irg)			        :: i,j,ind(3),ivec(3),ik,ir,jr,nn,IPIV(2*ira+1),nmin,nmax,izero, &
							   icount,count_rate,count_max,newcount
complex(kind=dbl)          		:: M(2*ira+1,2*ira+1),alph(2*ira+1)
complex(kind=dbl),allocatable 	:: CGinv(:,:),Mcp(:,:)

 allocate(CG(2*ira+1,2*ira+1),W(2*ira+1),Mcp(2*ira+1,2*ira+1),CGinv(2*ira+1,2*ira+1))
 pre = 2.0*sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18

! scaling factor for excitation error (2*k_0)
 pre2 = 2.0/sngl(mLambda)

! min and max of systematic row
 nmin = -ira
 nmax = ira

! normal aborption potential Uprime_0
 ind = (/0,0,0/)
 call CalcUcg(ind)
 upzero = pre*rlp%Vpmod

! tranmitted beam
 izero=1-nmin

! determine the dynamical matrix M (all but the diagonal)
 nn = 2*ira+1
 M(1:nn,1:nn) = cmplx(0.0,0.0)
 glen = ira*CalcLength(float(g),'r')

! i is the row index
 do i=nmin,nmax
  ir = i-nmin+1
  ind = i*g 
 write (*,*) ind
! j is the column index
  do j=nmin,nmax
   if (j.ne.i) then
    jr = j-nmin+1
    ivec = ind - j*g
! use Weickenmeier-Kohl scattering parameters and form factors
    call CalcUcg(ivec)
    M(ir,jr) = rlp%Ucg
   end if
  end do
 end do
!
! next we iterate over all incident beam directions, and for
! each direction we complete the M-matrix (diagonal).
!
! time the computation
 call system_clock(icount,count_rate,count_max)
 dkt = 2.0*ktmax/float(ns-1)
 io_real(1) = dkt
 call WriteValue(' beam tilt step size = ', io_real, 1, "(F8.4)")
 kk = CalcLength(k,'r')
 gg = CalcLength(float(g),'r')
 k =  k/sngl(mLambda)/kk
 kz = 1.0/mLambda

! open the unformatted output file
 open (unit=15,form='unformatted',status ='unknown')
 write (15) nn
 write (15) ns
 write (15) g
 write (15) k,kz

! loop over the beam directions
 frac=0.05
 do ik = 1,ns
  if (float(ik)/float(ns) .gt. frac) then
    write (*,"(1x,I3,' percent completed')") int(100.0*frac)
    frac = frac + 0.05
  end if

! rescale the wavevector and foil normal
  kt = k + dkt*(float(ik-ns/2))*g
  kk = CalcLength(kt,'r')
  kt = kt/sngl(mLambda)/kk

! then complete the diagonal of the M matrix
! i is the row index
  do i=nmin,nmax
   ir = i-nmin+1
   ind = i*g 
! get the excitation error
   if (i.eq.0) then
    exer = 0.0
   else
    exer = Calcsg(float(ind),kt,f)
   endif
! and multiply with 2k_0 and store, along with Uprime_0
   M(ir,ir) = cmplx(pre2*exer,upzero)
  end do

!
! next, compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
!

! first, make a copy of M, since BWsolve would destroy M
  Mcp = M

! then get the eigenvalues and eigenvectors
  call BWsolve(Mcp,W,CG,CGinv,nn,IPIV)

! the alpha coefficients are in the izero column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
  kttb = dkt*(float(ik-ns/2))
  kn = sqrt(kz**2-(kttb*gg)**2)
  W = W/cmplx(2.0*kn,0.0)
  do i=1,nn
   alph(i) = CGinv(i,izero)
  end do

! store eigenvalues in file along with excitation amplitudes, and eigenvector
! matrix (also the wave vector)
  write (15) kttb,kn
  write (15) W
  write (15) CG
  write (15) alph
 end do

! close the output file
 close(15, status='keep')

 call system_clock(newcount,count_rate,count_max)
 io_real(1)=float(newcount-icount)/float(count_rate)
 call WriteValue('Total computation time [s] ', io_real, 1, "(F)")

 mess = 'All data saved in temporary file fort.15'; call Message("(A)")
 
 call BWshow

end subroutine CalcSRBW


!--------------------------------------------------------------------------
!
! SUBROUTINE: BWshow 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create some output for this computation, using the temporary file fort.15
!
!> @date   4/24/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine BWshow

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript
use error

IMPLICIT NONE

integer(kind=irg)		:: ns,nn, io_int(1)
real(kind=sgl)           	:: io_real(2)

 open (unit=15,form='unformatted',status='old')
 read (15) nn
 read (15) ns
 mess = 'Systematic Row data set'; call Message("(A)")
 io_int(1) = nn
 call WriteValue(' Number of beams        = ', io_int, 1, "(I3)")
 io_int(1) = ns
 call WriteValue(' Number of orientations = ', io_int, 1, "(I3)")
 close (unit=15,status='keep')

! PostScript output file
! call ReadValue(' PostScript output filename : ', PS % psname, "(A)")
 call PS_openfile
 PS % pspage=0

 call ReadValue(' Minimum and maximum foil thickness : ', io_real, 2)
 call ReadValue(' Number of thicknesses : ',io_int, 1)

! draw grayscale BF and DF images and any other line drawings         
 call BWtoPS(nn,ns,io_int(1),io_real(1),io_real(2))

 call PS_closefile

end subroutine BWshow

!--------------------------------------------------------------------------
!
! SUBROUTINE: BWtoPS 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create PostScript output
!
!> @param nn number of beams
!> @param ns number of orientations
!> @param nt number of thicknesses
!> @param tmin minimum thickness
!> @param tmax maximum thickness
!
!> @date   4/24/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine BWtoPS(nn,ns,nt,tmin,tmax)

use local 
use io
use constants
use postscript
use graphics

IMPLICIT NONE


integer(kind=irg),INTENT(IN)		:: nn
integer(kind=irg),INTENT(IN)		:: ns
integer(kind=irg),INTENT(IN)		:: nt
real(kind=sgl),INTENT(IN)		:: tmin
real(kind=sgl),INTENT(IN)		:: tmax

real(kind=sgl)                 			:: images(ns,nt,nn),y(ns),workx(ns),worky(ns,nn),xmin,xmax,ymin,ymax,imax,imin, &
							   scl, lp, kn
integer(kind=irg)              		:: l,i,ii,io_int(1),npx,npy,ira,id,irmin,irmax, nf, isel,j
logical              				:: npg
character(12),parameter 			:: yt(5) = ['gamma^(j)   ', &
                                    			'k^(j)_z-k_0 ', &
                                    			'alpha^(j)   ', &
                                    			'q^(j)       ', &
                                    			'            ']
character(40),parameter 			:: gt(5) = ['Bloch wave eigenvalues                  ', &
                                   			'Bloch wave eigenvalues                  ', &
                                    			'Bloch wave excitation amplitudes        ', &
                                   			'Bloch wave absorption parameters        ', &
                                    			'Intensity profiles                      ']
real(kind=sgl),parameter          	:: xo(4)=[0.5,4.25,0.5,4.25],yo(4)=[5.5,5.5,1.0,1.0]
real(kind=sgl),parameter          	:: xx(11) = [0.25,2.35,4.45,0.25,2.35,4.45,0.25,2.35,4.45,0.25,2.35], &
                           				  yy(11) = [6.6,6.6,6.6,4.5,4.5,4.5,2.4,2.4,2.4,0.3,0.3]
complex(kind=sgl),allocatable     	:: CCGG(:,:)

 PS % pspage = 0
 call PS_newpage(.FALSE.,'BF/DF Images')

! next compute the beam intensities for a wedge shaped sample
 mess = 'Computing bright field and dark field images'; call Message("(A)")
 call BWtoI(nn,ns,nt,tmin,tmax,images)

 imanum = 0
! print the images on the top half of the page
 allocate(imaint(ns,nt))
 npx=ns
 npy=nt
 scl = 2.0
 imax=maxval(images)
 imin=minval(images)
 ira = (nn-1)/2
 id  = 6-ira-1

! is it s 2-beam case ?
 if (nn.eq.2) then 
  irmin = 1
  irmax = 2
  id = 0
  lp = 0
! or a multibeam case ?
 else
  irmin = -ira
  irmax = ira
  id = 6-ira-1
  lp = 6
 end if

 do l=irmin,irmax
  j = l+lp
  if ((j.gt.0).and.(j.lt.12)) then
   do i=1,ns
    do ii=1,nt
     imaint(i,ii)=int(255.0*images(i,ii,j-id))
    end do
   end do
   io_int(1) = j-id
   call WriteValue('  dumping image # ', io_int, 1, "(I3)")
   call PS_DumpImageDistort(xx(j),yy(j),npx,npy,scl,scl)
  end if
 end do

! next, ask if any other curves are needed
 nf = 0
 isel = 1
open(unit=25,file='sss.dat',status='unknown',form='unformatted')
 do while (isel.ne.0)
  mess = 'Output of various curves'; call Message("(//A)")
  mess = '0. Quit program' ; call Message("(A)")
  mess = '1. Bloch wave eigenvalues (gamma^(j))'; call Message("(A)") 
  mess = '2. Bloch wave eigenvalues (k^(j)_z-k_0)' ; call Message("(A)")
  mess = '3. Bloch wave excitation amplitudes ' ; call Message("(A)")
  mess = '4. Bloch wave absorption parameters ' ; call Message("(A)")
  mess = '5. Bloch wave amplitude profiles ' ; call Message("(A)")
  mess = '6. Store beam intensity vs. thickness' ; call Message("(A)")
  call ReadValue('   Enter your selection :', io_int, 1)
  isel = io_int(1)      

  if (isel.ne.0) then
   if (isel.ne.5) then 
     call ExtractBWdata(workx,worky,kn,isel,ns,nn)
   else
     allocate(CCGG(nn,nn))
     call ExtractBWdata(workx,worky,kn,isel,ns,nn,CCGG)
   endif
!
   nf = nf+1
   npg = .FALSE.
   if (nf.eq.5) nf=1
   if (nf.eq.1) npg = .TRUE.
   select case (isel) 

    case(1,2,3,4);
     xmin = minval(workx)
     xmax = maxval(workx)
     ymin = minval(worky)
     ymax = maxval(worky)

! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
     if (npg.eqv..TRUE.) then
      call PS_newpage(.FALSE.,'Bloch Wave Results')
     endif

! and call the axis routine to create the plot
! define the drawing location and size
     AX % axw = 4.2
     AX % xll = xo(nf)
     AX % yll = yo(nf)
     y(1:ns) = worky(1:ns,1)
     call axis(ns,workx,y,xmin,xmax,ymin,ymax,.FALSE.,.FALSE.,'lin','lin','CON',1, &
               'BOT','LEF',.FALSE.,.TRUE.,gt(isel),'kt/g',yt(isel))
! superimpose more curves
     do j=2,nn
      y(1:ns) = worky(1:ns,j)
      call axis(ns,workx,y,xmin,xmax,ymin,ymax,.FALSE.,.FALSE.,'lin','lin','CON',1, &
                'BOT','LEF',.TRUE.,.FALSE.,' ',' ',' ')
     end do

    case(5);
     write (25) CCGG

    case(6);
        open (dataunit,file='BWshow.out',status='unknown',form='unformatted')
        write (dataunit) images
        close (dataunit,status='keep')

    case default;

   end select

  end if

 end do

close (unit=25,status='keep')

end subroutine BWtoPS


!--------------------------------------------------------------------------
!
! SUBROUTINE: ExtractBWdata 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Produces postscript output
!
!> @param workx wave vector components
!> @param worky various quantities
!> @param kn normal component of wave vector
!> @param isel selector
!> @param ns number of orientations
!> @param nn number of beams
!> @param CCGG optional image array
!
!> @date   4/24/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine ExtractBWdata(workx,worky,kn,isel,ns,nn,CCGG)

use io

IMPLICIT NONE

real(kind=sgl),INTENT(OUT)	:: workx(ns)
real(kind=sgl),INTENT(OUT)	:: worky(ns,nn)
real(kind=sgl),INTENT(OUT)	:: kn
integer(kind=irg),INTENT(IN)	:: isel
integer(kind=irg),INTENT(IN)	:: ns
integer(kind=irg),INTENT(IN)	:: nn
complex(kind=sgl),optional 	:: CCGG(nn,nn)

real(kind=sgl)             		:: kttb, k(3),kzero
integer(kind=sgl)          		:: g(3),ipix, io_int(1),i,ik
complex(kind=sgl)          	:: W(nn), alph(nn), CG(nn,nn)

 open (unit=15,form='unformatted',status = 'old')
 if (isel.eq.5) then
   call ReadValue(' Enter pixel value for CCGG output ', io_int, 1)
   ipix = io_int(1)
 end if
 
! Systematic Row
 read (15) i
 read (15) i
 read (15) g
 read (15) k,kzero

 do ik = 1,ns
  read (15) kttb,kn
  read (15) W
  read (15) CG
  read (15) alph
   workx(ik) = kttb 

   if (isel.eq.1) worky(ik,1:nn) = real(W(1:nn))

   if (isel.eq.2) worky(ik,1:nn) = (kn+real(W(1:nn))-kzero)

   if (isel.eq.3) worky(ik,1:nn) = abs(alph(1:nn))**2

   if (isel.eq.4) worky(ik,1:nn) = imag(W(1:nn))
  
   if ((isel.eq.5).and.(ik.eq.ipix)) CCGG = CG

 end do

 close(unit=15,status='keep')

end subroutine ExtractBWdata

!--------------------------------------------------------------------------
!
! SUBROUTINE: BWtoI 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Convert Bloch wave parameters to intensities
!
!> @param ns number of orientations
!> @param nn number of beams
!> @param nt number of thicknesses
!> @param tmin minimum thickness
!> @param tmax maximum thickness
!> @param oname filename
!> @param images output images
! 
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine BWtoI(nn,ns,nt,tmin,tmax,images)

use local
use constants

integer(kind=irg),INTENT(IN)	:: nn
integer(kind=irg),INTENT(IN)	:: ns
integer(kind=irg),INTENT(IN)	:: nt
real(kind=sgl),INTENT(IN)	:: tmin
real(kind=sgl),INTENT(IN)	:: tmax
real(kind=sgl),INTENT(OUT) 	:: images(ns,nt,nn)

real(kind=sgl)				:: kt(3),kttb,kn,kzero,Wr(nn),Wi(nn)
integer(kind=irg)  			:: g(3),i
complex(kind=dbl)  		:: W(nn), alph(nn), CG(nn,nn), amp, diag(nn), q(nn)
character(15)      			:: fname
character(2)       			:: dtype

 open (unit=15,form='unformatted',status = 'old')
  read (15) i
  read (15) j
  read (15) g
  read (15) kt,kzero
  dz = (tmax-tmin)/float(nt-1)
write (*,*) dz, tmin, tmax, g, kt, kzero

  open(unit=25,file='ttt.dat',status='unknown',form='unformatted')
  do ik = 1,ns
   read (15) kttb,kn
   read (15) W
   W = W*2.0*cPi
   Wr = real(W) 
   Wi = imag(W)
   read (15) CG
   read (15) alph
   do i=1,nt
    z = tmin + dz*float(i-1)
    diag=exp(-z*Wi)*cmplx(cos(z*Wr),sin(z*Wr))*alph
    do j=1,nn
     images(ik,i,j) = abs(sum(CG(j,1:nn)*diag(1:nn)))**2 ! PGC cabs -> abs
    end do 
    if (i.eq.256) write (25) (images(ik,i,j),j=1,nn)
   end do
  end do
 close (unit=25,status='keep')

close (unit=15,status='keep')

end subroutine BWtoI




