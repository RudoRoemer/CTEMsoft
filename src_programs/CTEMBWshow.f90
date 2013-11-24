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
! CTEMsoft2013:CTEMBWshow.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMBWshow 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Program reads Bloch wave input file from various programs and displays output
!
! 
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program CTEMBWshow

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript
use error
use dynamical

IMPLICIT NONE

integer(kind=irg)         	:: na,nb,nt,ns,nn,n1,io_int(1)
real(kind=sgl)            	:: io_real(2)
character(fnlen)             	:: oname
character(2)              	:: dtype
character(15)      		:: fname

 progname = 'CTEMBWshow.f90'
 progdesc = 'Display program for various Bloch wave data files'
 call CTEMsoft
                                                                               
! get the input filename
 call ReadValue(' Bloch wave input filename : ', oname, "(A)")

! open file, and determine what kind of dataset it is
 open (unit=15,file=trim(oname),form='unformatted',status='old')
 read (15) dtype
 read (15) fname
 select case (dtype)

! Two Beam
  case('TB'); read (15) nn
              read (15) ns
              mess = 'Two Beam data set'; call Message("(A)")
              io_int(1) = ns
              call WriteValue('    Number of orientations = ', io_int, 1, "(I3)")
        
! Systematic Row
  case('SR'); read (15) nn
              read (15) ns
              mess = 'Systematic Row data set'; call Message("(A)")
	      io_int(1) = nn
              call WriteValue('    Number of beams        = ', io_int, 1, "(I3)")
	      io_int(1) = n1
              call WriteValue('    Number of orientations = ', io_int, 1, "(I3)")
        
! Zone Axis
  case('ZA'); read (15) na,nb
              read (15) ns,nt
              nn = (2*na+1)*(2*nb+1)
              mess = 'Zone Axis data set'; call Message("(A)")
              
  case default; call FatalError('CTEMBWshow','unknown data type')

 end select
 close (unit=15,status='keep')

! PostScript output file
 call PS_openfile
 PS % pspage=0
!
 if ((dtype.eq.'TB').or.(dtype.eq.'SR')) then
  call ReadValue(' Minimum and maximum foil thickness : ', io_real, 2)
  call ReadValue(' Number of thicknesses : ', io_int, 1)

! first, draw grayscale BF and DF images          
  call BWtoPS(nn,ns,io_int(1),io_real(1),io_real(2),oname)
 endif

 call PS_closefile
end program


!--------------------------------------------------------------------------
!
! SUBROUTINE: BWtoPS 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Produces postscript output
!
!> @param nn Number of beams
!> @param ns Number of orientations
!> @param nt Number of thicknesses
!> @param tmin minimum thickness
!> @param tmax maximum thickness
!> @param oname output filename
! 
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine BWtoPS(nn,ns,nt,tmin,tmax,oname)

use local 
use io
use constants
use postscript
use graphics

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: nn
integer(kind=irg),INTENT(IN)	:: ns
integer(kind=irg),INTENT(IN)	:: nt
real(kind=sgl),INTENT(IN)	:: tmin
real(kind=sgl),INTENT(IN)	:: tmax
character(*),INTENT(IN)		:: oname

real(kind=sgl)       			:: images(ns,nt,nn),y(ns),workx(ns),worky(ns,nn), &
						   xmin,xmax,ymin,ymax,imax,imin,scl,kn
integer(kind=irg)    			:: npx,npy,ira,id,irmin,irmax,lp,l,j,i,ii,nf,io_int(1),isel
logical              			:: npg
character(12),parameter 		:: yt(5) = ['gamma^(j)   ', &
                                  		  'k^(j)_z-k_0 ', &
                                   		  'alpha^(j)   ', &
                                    		  'q^(j)       ', &
                                    		  '            ']
character(40),parameter 		:: gt(5) = ['Bloch wave eigenvalues                  ', &
                                    		'Bloch wave eigenvalues                  ', &
                                    		'Bloch wave excitation amplitudes        ', &
                                    		'Bloch wave absorption parameters        ', &
                                    		'Intensity profiles                      ']
real(kind=sgl),parameter  	:: xo(4)=[0.5,4.25,0.5,4.25],yo(4)=[5.5,5.5,1.0,1.0]
real(kind=sgl),parameter   	:: xx(11) = [0.25,2.35,4.45,0.25,2.35,4.45,0.25,2.35,4.45,0.25,2.35], &
                           			   yy(11) = [6.6,6.6,6.6,4.5,4.5,4.5,2.4,2.4,2.4,0.3,0.3]

 PS % pspage = 0
 call PS_newpage(.FALSE.,'BF/DF Images')

! next compute the beam intensities for a wedge shaped sample
 mess = 'Computing bright field and dark field images'; call Message("(A)")
 call BWtoI(nn,ns,nt,tmin,tmax,oname,images)

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
   call WriteValue('    dumping image # ', io_int, 1, "(I3)")
   call PS_DumpImage(xx(j),yy(j),npx,npy,scl)
  end if
 end do

! next, ask if any other curves are needed
 nf = 0
 isel = 1
 do while (isel.ne.0)
  mess = 'Output of various curves'; call Message("(A)")
  mess = '0. Quit program' ; call Message("(A)")
  mess = '1. Bloch wave eigenvalues (gamma^(j))'; call Message("(A)") 
  mess = '2. Bloch wave eigenvalues (k^(j)_z-k_0)' ; call Message("(A)")
  mess = '3. Bloch wave excitation amplitudes ' ; call Message("(A)")
  mess = '4. Bloch wave absorption parameters ' ; call Message("(A)")
  mess = '5. Store beam intensity vs. thickness' ; call Message("(A)")
  call ReadValue(' Enter your selection :', io_int, 1)
  isel = io_int(1)      

  if (isel.ne.0) then
   call ExtractBWdata(workx,worky,kn,isel,ns,nn,oname)
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
        open (dataunit,file='BWshow.out',status='unknown',form='unformatted')
        write (dataunit) images
        close (dataunit,status='keep')

    case default;

   end select

   mess  = 'done'; call Message("(A)")
  end if
 end do

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
!> @param oname filename
! 
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine ExtractBWdata(workx,worky,kn,isel,ns,nn,oname)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(OUT)	:: workx(ns)
real(kind=sgl),INTENT(OUT)	:: worky(ns,nn)
real(kind=sgl),INTENT(OUT)	:: kn
integer(kind=irg),INTENT(IN)	:: isel
integer(kind=irg),INTENT(IN)	:: ns
integer(kind=irg),INTENT(IN)	:: nn
character(*),INTENT(IN) 		:: oname

real(kind=sgl)   			:: kttb, kt(3), k(3), kzero
integer(kind=irg)			:: g(3),i,ik
complex(kind=dbl)			:: W(nn), alph(nn), CG(nn,nn), amp, diag(nn)
character(2)     				:: dtype
character(15)				:: fname

 open (unit=15,file=trim(oname),form='unformatted',status = 'old')
 read (15) dtype
 read (15) fname

! Two Beam or Systematic Row
 if ((dtype.eq.'TB').or.(dtype.eq.'SR')) then 
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

   if (isel.eq.2) worky(ik,1:nn) = -(kn+real(W(1:nn))+kzero)

   if (isel.eq.3) worky(ik,1:nn) = abs(alph(1:nn))**2

   if (isel.eq.4) worky(ik,1:nn) = imag(W(1:nn))

  end do

 end if

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
subroutine BWtoI(nn,ns,nt,tmin,tmax,oname,images)

use local
use constants

integer(kind=irg),INTENT(IN)	:: nn
integer(kind=irg),INTENT(IN)	:: ns
integer(kind=irg),INTENT(IN)	:: nt
real(kind=sgl),INTENT(IN)	:: tmin
real(kind=sgl),INTENT(IN)	:: tmax
character(*),INTENT(IN)		:: oname
real(kind=sgl),INTENT(OUT) 	:: images(ns,nt,nn)

real(kind=sgl)				:: kt(3),kttb,kn,kzero,Wr(nn),Wi(nn)
integer(kind=irg)  			:: g(3)
complex(kind=dbl)  		:: W(nn), alph(nn), CG(nn,nn), amp, diag(nn), q(nn)
character(15)      			:: fname
character(2)       			:: dtype

 open (unit=15,file=oname,form='unformatted',status = 'old')
 read (15) dtype
 read (15) fname

! Two Beam
 if ((dtype.eq.'TB').or.(dtype.eq.'SR')) then 
  read (15) i
  read (15) i
  read (15) g
  read (15) kt,kzero
  dz = (tmax-tmin)/float(nt)

  do ik = 1,ns
   read (15) kttb,kn
   read (15) W
   read (15) CG
   read (15) alph

   do i=1,nt
    z = dz*float(i)
    arg = 2.0*sngl(cPi)*z
    Wr = arg*real(W)
    Wi = arg*imag(W)
    q = cmplx(cos(Wr),sin(Wr))
    diag=exp(Wi)*cmplx(cos(Wr),sin(Wr))*alph

    do j=1,nn
     amp = sum(CG(j,1:nn)*diag(1:nn))
     images(ik,i,j) = abs(amp)**2
    end do 

   end do

  end do

 end if

 close (unit=15,status='keep')

end subroutine BWtoI






