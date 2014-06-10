! ###################################################################
! Copyright (c) 2014, Marc De Graef/Carnegie Mellon University
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
! CTEMsoft:constants.f90
!--------------------------------------------------------------------------
!
! MODULE: constants
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief 
!> Definition of constants and constant arrays used by other routines
!
!> @details  
!> physical and mathematical constants used by various programs; periodic
!>  table information; atomic weights; 
! 
!> @date 1/5/99   MDG 1.0 original
!> @date 5/18/01  MDG 2.0 f90
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 03/19/13 MDG 3.0 added atomic weights 
!> @date 01/10/14 MDG 4.0 new version
!> @date 04/29/14 MDG 4.1 constants updated from NIST physical constants tables
!--------------------------------------------------------------------------

module constants

use local

IMPLICIT NONE

! various physical constants
!> cPi		= pi [dimensionless]
!> cLight	= velocity of light [m/s]
!> cPlanck	= Planck''s constant [Js]
!> cBoltzmann	= Boltmann constant [J/K]
!> cPermea	= permeability of vacuum [dimensionless]
!> cPermit	= permittivity of vacuum [F/m]
!> cCharge	= electron charge [C]
!> cRestmass	= electron rest mass [kg]
!> cMoment	= electron magnetic moment [J/T]
!> cJ2eV	= Joules per eV
!> cAvogadro	= Avogadro's constant [mol^-1]

real(kind=dbl), parameter :: cPi=3.141592653589793238D0, cLight = 299792458.D0, &
                             cPlanck = 6.62606957D-34, cBoltzmann = 1.3806488D-23,  &
                             cPermea = 1.2566370614D-6, cPermit = 8.854187817D-12, &
                             cCharge = 1.602176565D-19, cRestmass = 9.10938291D-31, &
                             cMoment = 9.28476430D-24, cJ2eV = 1.602176565D-19, &
			      cAvogadro = 6.02214129D23

!> element symbols (we'll do 1-98 for all parameter lists)
character(2), parameter :: ATOM_sym(98)=(/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
                                          'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca', &
                                          'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                                          'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr', &
                                          'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                                          'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                                          'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                                          'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', &
                                          'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                                          'Pa',' U','Np','Pu','Am','Cm','Bk','Cf'/)

!> Shannon-Prewitt ionic radii in nanometer
real(kind=sgl), parameter :: ATOM_SPradii(98)=(/0.010,0.010,0.680,0.300,0.160,0.150,0.148,0.146,0.133,0.500, &
                                                0.098,0.065,0.450,0.380,0.340,0.190,0.181,0.500,0.133,0.940, &
                                                0.068,0.060,0.740,0.690,0.670,0.640,0.630,0.620,0.720,0.740, &
                                                0.113,0.073,0.580,0.202,0.196,0.500,0.148,0.110,0.880,0.770, &
                                                0.067,0.068,0.500,0.500,0.500,0.860,0.126,0.970,0.132,0.930, &
                                                0.076,0.222,0.219,0.500,0.167,0.129,0.104,0.111,0.500,0.108, &
                                                0.050,0.104,0.500,0.970,0.500,0.990,0.500,0.960,0.500,0.940, &
                                                0.050,0.050,0.680,0.600,0.520,0.500,0.500,0.500,0.137,0.112, &
                                                0.140,0.132,0.740,0.230,0.227,0.500,0.175,0.137,0.111,0.990, &
                                                0.090,0.083,0.500,0.108,0.500,0.500,0.500,0.500/)

!> atomic (metallic) radii in nanometer (0.100 if not known/applicable)
real(kind=sgl), parameter :: ATOM_MTradii(98)=(/0.100,0.100,0.156,0.112,0.100,0.100,0.100,0.100,0.100,0.100, &
                                                0.191,0.160,0.142,0.100,0.100,0.100,0.100,0.100,0.238,0.196, &
                                                0.160,0.146,0.135,0.128,0.136,0.127,0.125,0.124,0.128,0.137, &
                                                0.135,0.139,0.125,0.116,0.100,0.100,0.253,0.215,0.181,0.160, &
                                                0.147,0.140,0.135,0.133,0.134,0.137,0.144,0.152,0.167,0.158, &
                                                0.161,0.143,0.100,0.100,0.270,0.224,0.187,0.182,0.182,0.181, &
                                                0.100,0.100,0.204,0.178,0.177,0.175,0.176,0.173,0.174,0.193, &
                                                0.173,0.158,0.147,0.141,0.137,0.135,0.135,0.138,0.144,0.155, &
                                                0.171,0.174,0.182,0.168,0.100,0.100,0.100,0.100,0.100,0.180, &
                                                0.163,0.154,0.150,0.164,0.100,0.100,0.100,0.100/)

!> atom colors for PostScript drawings
character(3), parameter :: ATOM_color(98)=(/'blu','grn','blu','blu','red','bro','blu','red','grn','grn', &
                                            'blu','pnk','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','grn'/)

!> atomic weights for things like density computations (from NIST elemental data base)
real(kind=sgl),parameter	:: ATOM_weights(98) = (/1.00794, 4.002602, 6.941, 9.012182, 10.811, &
							12.0107, 14.0067, 15.9994, 18.9984032, 20.1797, &
							22.98976928, 24.3050, 26.9815386, 28.0855, 30.973762, &
							32.065, 35.453, 39.948, 39.0983, 40.078, &
							44.955912, 47.867, 50.9415, 51.9961, 54.938045, &
							55.845, 58.933195, 58.6934, 63.546, 65.38, &
							69.723, 72.64, 74.92160, 78.96, 79.904, &
							83.798, 85.4678, 87.62, 88.90585, 91.224, &
							92.90638, 95.96, 98.9062, 101.07, 102.90550, &
							106.42, 107.8682, 112.411, 114.818, 118.710, &
							121.760, 127.60, 126.90447, 131.293, 132.9054519, &
							137.327, 138.90547, 140.116, 140.90765, 144.242, &
							145.0, 150.36, 151.964, 157.25, 158.92535, &
							162.500, 164.93032, 167.259, 168.93421, 173.054, &
							174.9668, 178.49, 180.94788, 183.84, 186.207, &
							190.23, 192.217, 195.084, 196.966569, 200.59, &
							204.3833, 207.2, 208.98040, 209.0, 210.0, &
							222.0, 223.0, 226.0, 227.0, 232.03806, &
							231.03588, 238.02891, 237.0, 244.0, 243.0, &
							247.0, 251.0, 252.0 /)


! these are a bunch of constants used for Lambert and related projections; they are all in double precision
type LambertParametersType
	real(kind=dbl)		:: Pi=3.141592653589793D0    	!  pi
	real(kind=dbl)		:: iPi=0.318309886183791D0    	!  1/pi
	real(kind=dbl)		:: sPi=1.772453850905516D0	!  sqrt(pi)
	real(kind=dbl)		:: sPio2=1.253314137315500D0	!  sqrt(pi/2)
	real(kind=dbl)		:: sPi2=0.886226925452758D0	!  sqrt(pi)/2
	real(kind=dbl)		:: srt=0.866025403784439D0    	!  sqrt(3)/2
	real(kind=dbl)		:: isrt=0.577350269189626D0    !  1/sqrt(3)
	real(kind=dbl)		:: alpha=1.346773687088598D0   !  sqrt(pi)/3^(1/4)
	real(kind=dbl)		:: rtt=1.732050807568877D0  	!  sqrt(3)
	real(kind=dbl)		:: prea=0.525037567904332D0    !  3^(1/4)/sqrt(2pi)
	real(kind=dbl)		:: preb=1.050075135808664D0   	!  3^(1/4)sqrt(2/pi)
	real(kind=dbl)		:: prec=0.906899682117109D0    !  pi/2sqrt(3)
	real(kind=dbl)		:: pred=2.094395102393195D0   	!  2pi/3
	real(kind=dbl)		:: pree=0.759835685651593D0   	!  3^(-1/4)
	real(kind=dbl)		:: pref=1.381976597885342D0	!  sqrt(6/pi)
! the following constants are used for the cube to quaternion hemisphere mapping
	real(kind=dbl)		:: a=1.925749019958253D0	! pi^(5/6)/6^(1/6)
	real(kind=dbl)		:: ap=2.145029397111025D0	! pi^(2/3)
	real(kind=dbl)		:: sc=0.897772786961286D0	! a/ap
	real(kind=dbl)		:: beta=0.962874509979126D0	! pi^(5/6)/6^(1/6)/2
	real(kind=dbl)		:: R1=1.330670039491469D0	! (3pi/4)^(1/3)
	real(kind=dbl)		:: r2=1.414213562373095D0	! sqrt(2)
	real(kind=dbl)		:: r22=0.707106781186547D0	! 1/sqrt(2)
	real(kind=dbl)		:: pi12=0.261799387799149D0	! pi/12
	real(kind=dbl)		:: pi8=0.392699081698724D0	! pi/8
	real(kind=dbl)		:: prek=1.643456402972504D0	! R1 2^(1/4)/beta
	real(kind=dbl)		:: r24=4.898979485566356D0	! sqrt(24)
        real(kind=dbl)	        :: rvmax2=13131225398.930372238D0      ! square of max rodrigues vector length  
	real(kind=dbl)		:: tfit(7)=(/ -0.5000096149170321D0, -0.02486606148871731D0, &
                                             -0.004549381779362819D0, 0.0005118668366387526D0, &
                                             -0.0016500827333575548D0, 0.0007593352203388718D0, &
                                             -0.0002040422502566876D0 /)	! fit parameters
	real(kind=dbl)              :: BP(6)= (/ 0.D0, 1.D0, 0.577350269189626D0, 0.414213562373095D0, 0.D0,  &
                                              0.267949192431123D0 /)       ! used for Fundamental Zone determination
end type LambertParametersType

type(LambertParametersType)        :: LPs

end module
