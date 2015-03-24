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
! EMsoft:constants.f90
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
!> @date 07/06/14 MDG 4.2 added omegamax to Lambert constant type
!> @date 08/11/14 MDG 4.3 added infty for handling of +Infinity in rotations module
!> @date 08/11/14 MDG 4.4 added epsijk option to package
!> @date 09/30/14 MDG 4.5 added some additional comments about epsijk
!> @date 10/02/14 MDG 4.6 removed omegamax again, since we now have properly dealt with 180 degree rotations
!> @date 03/11/15 MDG 4.7 added some additional comments about epsijk 
!--------------------------------------------------------------------------

module constants

use local

IMPLICIT NONE

! ****************************************************
! ****************************************************
! ****************************************************
! used to change the sign of the permutation symbol from Adam Morawiec's book to
! the convention used for the EMsoft package.  If you want to use Adam's convention,
! both of these parameters should be set to +1; -1 will change the sign everywhere
! for all representations that involve the unit vector.  The quaternion product is 
! also redefined to include the epsijk parameter.  Doing so guarantees that the 
! quat_Lp operator ALWAYS returns an active result, regardless of the choice of epsijk;
! quat_LPstar ALWAYS returns a passive result.

! Uncomment these for an alternative way of doing things
real(kind=sgl), parameter :: epsijk = -1.0
real(kind=dbl), parameter :: epsijkd = -1.D0

! uncomment these for the Morawiec version.
!real(kind=sgl), parameter :: epsijk = 1.0
!real(kind=dbl), parameter :: epsijkd = 1.D0

! In the first case, epsijk=-1, the rotation 120@[111] will result in 
! an axis angle pair of [111], 2pi/3.  In the second case, the axis-angle 
! pair will be -[111], 2pi/3.  In all cases, the rotations are interpreted
! in the passive sense.  The case epsijk=+1 corresponds to the mathematically 
! consistent case, using the standard definition for the quaternion product; in
! the other case, epsijk=-1, one must redefine the quaternion product in order
! to produce consistent results.  This takes a lengthy explanation ... see the
! rotations tutorial paper for an in-depth explanation.  These changes propagate
! to a number of files, notably quaternions.f90, and everywhere else that quaternions
! and rotations in general are used.
! ****************************************************
! ****************************************************
! ****************************************************

! the rotations.f90 routines need to have access to the value +Infinity
! which is defined here (using the LaTeX name infty)
INTEGER,private :: inf
REAL,public :: infty
EQUIVALENCE (inf,infty) ! stores two variable at the same address
DATA inf/z'7f800000'/ !Hex for +Infinity

INTEGER(kind=8),private :: infd
REAL(kind=dbl),public :: inftyd
EQUIVALENCE (infd,inftyd) ! stores two variable at the same address
DATA infd/z'7FF0000000000000'/ !Hex for +Infinity




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
real(kind=sgl),parameter        :: ATOM_weights(98) = (/1.00794, 4.002602, 6.941, 9.012182, 10.811, &
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
        real(kind=dbl)          :: Pi=3.141592653589793D0       !  pi
        real(kind=dbl)          :: iPi=0.318309886183791D0      !  1/pi
        real(kind=dbl)          :: sPi=1.772453850905516D0      !  sqrt(pi)
        real(kind=dbl)          :: sPio2=1.253314137315500D0    !  sqrt(pi/2)
        real(kind=dbl)          :: sPi2=0.886226925452758D0     !  sqrt(pi)/2
        real(kind=dbl)          :: srt=0.866025403784439D0      !  sqrt(3)/2
        real(kind=dbl)          :: isrt=0.577350269189626D0    !  1/sqrt(3)
        real(kind=dbl)          :: alpha=1.346773687088598D0   !  sqrt(pi)/3^(1/4)
        real(kind=dbl)          :: rtt=1.732050807568877D0      !  sqrt(3)
        real(kind=dbl)          :: prea=0.525037567904332D0    !  3^(1/4)/sqrt(2pi)
        real(kind=dbl)          :: preb=1.050075135808664D0     !  3^(1/4)sqrt(2/pi)
        real(kind=dbl)          :: prec=0.906899682117109D0    !  pi/2sqrt(3)
        real(kind=dbl)          :: pred=2.094395102393195D0     !  2pi/3
        real(kind=dbl)          :: pree=0.759835685651593D0     !  3^(-1/4)
        real(kind=dbl)          :: pref=1.381976597885342D0     !  sqrt(6/pi)
! the following constants are used for the cube to quaternion hemisphere mapping
        real(kind=dbl)          :: a=1.925749019958253D0        ! pi^(5/6)/6^(1/6)
        real(kind=dbl)          :: ap=2.145029397111025D0       ! pi^(2/3)
        real(kind=dbl)          :: sc=0.897772786961286D0       ! a/ap
        real(kind=dbl)          :: beta=0.962874509979126D0     ! pi^(5/6)/6^(1/6)/2
        real(kind=dbl)          :: R1=1.330670039491469D0       ! (3pi/4)^(1/3)
        real(kind=dbl)          :: r2=1.414213562373095D0       ! sqrt(2)
        real(kind=dbl)          :: r22=0.707106781186547D0      ! 1/sqrt(2)
        real(kind=dbl)          :: pi12=0.261799387799149D0     ! pi/12
        real(kind=dbl)          :: pi8=0.392699081698724D0      ! pi/8
        real(kind=dbl)          :: prek=1.643456402972504D0     ! R1 2^(1/4)/beta
        real(kind=dbl)          :: r24=4.898979485566356D0      ! sqrt(24)
        real(kind=dbl)          :: tfit(16) = (/1.0000000000018852D0, -0.5000000002194847D0, & 
                                             -0.024999992127593126D0, - 0.003928701544781374D0, & 
                                             -0.0008152701535450438D0, - 0.0002009500426119712D0, & 
                                             -0.00002397986776071756D0, - 0.00008202868926605841D0, & 
                                             +0.00012448715042090092D0, - 0.0001749114214822577D0, & 
                                             +0.0001703481934140054D0, - 0.00012062065004116828D0, & 
                                             +0.000059719705868660826D0, - 0.00001980756723965647D0, & 
                                             +0.000003953714684212874D0, - 0.00000036555001439719544D0 /)
	real(kind=dbl)              :: BP(6)= (/ 0.D0, 1.D0, 0.577350269189626D0, 0.414213562373095D0, 0.D0,  &
                                              0.267949192431123D0 /)       ! used for Fundamental Zone determination in so3 module
end type LambertParametersType

type(LambertParametersType)        :: LPs




! Following numbers are coefficients used to calculate the exponential of a matrix


! The following two arrays are used to determine the FZtype (FZtarray) and primary rotation axis order (FZoarray)
! for each of the 32 crystallographic point group symmetries (in the order of the International Tables)
!
! 1 (C1), -1 (Ci), [triclinic]
! 2 (C2), m (Cs), 2/m (C2h), [monoclinic]
! 222 (D2), mm2 (C2v), mmm (D2h), [orthorhombic]
! 4 (C4), -4 (S4), 4/m (C4h), 422 (D4), 4mm (C4v), -42m (D2d), 4/mmm (D4h), [tetragonal]
! 3 (C3), -3 (C3i), 32 (D3), 3m (C3v), -3m (D3d), [trigonal]
! 6 (C6), -6 (C3h), 6/m (C6h), 622 (D6), 6mm (C6v), -6m2 (D3h), 6/mmm (D6h), [hexagonal]
! 23 (T), m3 (Th), 432 (O), -43m (Td), m-3m (Oh) [cubic]
!
! FZtype
! 0        no symmetry at all
! 1        cyclic symmetry
! 2        dihedral symmetry
! 3        tetrahedral symmetry
! 4        octahedral symmetry
!
! these parameters are used in the so3 module
!
integer(kind=irg),parameter     :: FZtarray(32) = (/ 0,0,1,1,1,2,2,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,2,2,2,2,3,3,4,3,4 /)
integer(kind=irg),parameter     :: FZoarray(32) = (/ 0,0,2,2,2,2,2,2,4,4,4,4,4,4,4,3,3,3,3,3,6,6,6,6,6,6,6,0,0,0,0,0 /)


end module
