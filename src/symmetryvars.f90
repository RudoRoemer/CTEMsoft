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
! CTEMsoft2013:symmetryvars.f90
!--------------------------------------------------------------------------
!
! MODULE: symmetryvars
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief all symmetry and related variable definitions
!
!> @details  contains a list of space group names and generator strings,
!> numbering according to International Tables for Crystallography
!> [hexagonal groups in rhombohedral setting
!> are described in locations 231-237]
!> 
!> The space group information is encoded in the following way:
!> From the International Crystallographic Tables one finds that
!> the point symmetry parts of the space group generators are 
!> represented by 14 out of 71 different matrices.  When expressed 
!> in the crystallographic reference frame, those matrices contain 
!> only the entries 1, 0, and -1.  
!> 
!> The 14 matrices are represented by the lower case 
!> letters a through n. (see program and Appendix A3, page 666, for details)
!> The translational parts of the space group generators are limited
!> to the following set (encoded with upper-case letters):
!> 
!> A             1/6
!> B             1/4
!> C             1/3
!> D             1/2
!> E             2/3
!> F             3/4
!> G             5/6
!> O             0
!> X             -3/8
!> Y             -1/4
!> Z             -1/8
!>
!> In encoding the space groups we have selected the [unique axis b, 
!> cell choice 1] settings for the monoclinic point groups.
!> The first regular space group with a multiple origin choice is 
!> No. 48 (Pnnn).  
!>
!> <b> WARNING FOR  THE USER: </b> This is a very tricky file!  Make sure that you 
!> understand all the details before you attempt to change anything!!!
! 
!> @note  9/25/2011: corrected the generators for space groups 39, 90, 107, 108, 131, 214, and 222.\n
!> [Thank you Matthew O'Brien for pointing out the errors !]\n
!> [Thank you also to Marco Schowalter for pointing out an error in space group 222; corrected on 9/24/2012]

!> @date 1/5/99   MDG 1.0 original
!> @date   5/19/01  MDG 2.0 f90
!> @date  11/27/01  MDG 2.1 added kind support
!> @date   03/19/13 MDG 3.0 updated file
!--------------------------------------------------------------------------

module symmetryvars

use local

!>  SYM_SGname all space group names
! TRICLINIC SPACE GROUPS
character(11),parameter :: SYM_SGname(237)= (/" P  1      " ," P -1      ", & ! MONOCLINIC SPACE GROUPS
        " P 2       " ," P 21      " ," C 2       " ," P m       ", &
        " P c       " ," C m       " ," C c       " ," P 2/m     ", &
        " P 21/m    " ," C 2/m     " ," P 2/c     " ," P 21/c    ", &
        " C 2/c     ", &                                              ! ORTHORHOMBIC SPACE GROUPS
        " P 2 2 2   " ," P 2 2 21  " ," P 21 21 2 " ," P 21 21 21", &
        " C 2 2 21  " ," C 2 2 2   " ," F 2 2 2   " ," I 2 2 2   ", &
        " I 21 21 21" ," P m m 2   " ," P m c 21  " ," P c c 2   ", &
        " P m a 2   " ," P c a 21  " ," P n c 2   " ," P m n 21  ", &
        " P b a 2   " ," P n a 21  " ," P n n 2   " ," C m m 2   ", &
        " C m c 21  " ," C c c 2   " ," A m m 2   " ," A b m 2   ", &
        " A m a 2   " ," A b a 2   " ," F m m 2   " ," F d d 2   ", &
        " I m m 2   " ," I b a 2   " ," I m a 2   " ," P m m m   ", &
        " P n n n   " ," P c c m   " ," P b a n   " ," P m m a   ", &
        " P n n a   " ," P m n a   " ," P c c a   " ," P b a m   ", &
        " P c c n   " ," P b c m   " ," P n n m   " ," P m m n   ", &
        " P b c n   " ," P b c a   " ," P n m a   " ," C m c m   ", &
        " C m c a   " ," C m m m   " ," C c c m   " ," C m m a   ", &
        " C c c a   " ," F m m m   " ," F d d d   " ," I m m m   ", &
        " I b a m   " ," I b c a   " ," I m m a   ", &                ! TETRAGONAL SPACE GROUPS  
        " P 4       " ," P 41      " ," P 42      " ," P 43      ", &
        " I 4       " ," I 41      " ," P -4      " ," I -4      ", &
        " P 4/m     " ," P 42/m    " ," P 4/n     " ," P 42/n    ", &
        " I 4/m     " ," I 41/a    " ," P 4 2 2   " ," P 4 21 2  ", &
        " P 41 2 2  " ," P 41 21 2 " ," P 42 2 2  " ," P 42 21 2 ", &
        " P 43 2 2  " ," P 43 21 2 " ," I 4 2 2   " ," I 41 2 2  ", &
        " P 4 m m   " ," P 4 b m   " ," P 42 c m  " ," P 42 n m  ", &
        " P 4 c c   " ," P 4 n c   " ," P 42 m c  " ," P 42 b c  ", &
        " I 4 m m   " ," I 4 c m   " ," I 41 m d  " ," I 41 c d  ", &
        " P -4 2 m  " ," P -4 2 c  " ," P -4 21 m " ," P -4 21 c ", &
        " P -4 m 2  " ," P -4 c 2  " ," P -4 b 2  " ," P -4 n 2  ", &
        " I -4 m 2  " ," I -4 c 2  " ," I -4 2 m  " ," I -4 2 d  ", &
        " P 4/m m m " ," P 4/m c c " ," P 4/n b m " ," P 4/n n c ", &
        " P 4/m b m " ," P 4/m n c " ," P 4/n m m " ," P 4/n c c ", &
        " P 42/m m c" ," P 42/m c m" ," P 42/n b c" ," P 42/n n m", &
        " P 42/m b c" ," P 42/m n m" ," P 42/n m c" ," P 42/n c m", &
        " I 4/m m m " ," I 4/m c m " ," I 41/a m d" ," I 41/a c d", & ! RHOMBOHEDRAL SPACE GROUPS  
        " P 3       " ," P 31      " ," P 32      " ," R 3       ", &
        " P -3      " ," R -3      " ," P 3 1 2   " ," P 3 2 1   ", &
        " P 31 1 2  " ," P 31 2 1  " ," P 32 1 2  " ," P 32 2 1  ", &
        " R 3 2     " ," P 3 m 1   " ," P 3 1 m   " ," P 3 c 1   ", &
        " P 3 1 c   " ," R 3 m     " ," R 3 c     " ," P -3 1 m  ", &
        " P -3 1 c  " ," P -3 m 1  " ," P -3 c 1  " ," R -3 m    ", &
        " R -3 c    ", &                                              ! HEXAGONAL SPACE GROUPS   
        " P 6       " ," P 61      " ," P 65      " ," P 62      ", &
        " P 64      " ," P 63      " ," P -6      " ," P 6/m     ", &
        " P 63/m    " ," P 6 2 2   " ," P 61 2 2  " ," P 65 2 2  ", &
        " P 62 2 2  " ," P 64 2 2  " ," P 63 2 2  " ," P 6 m m   ", &
        " P 6 c c   " ," P 63 c m  " ," P 63 m c  " ," P -6 m 2  ", &
        " P -6 c 2  " ," P -6 2 m  " ," P -6 2 c  " ," P 6/m m m ", &
        " P 6/m c c " ," P 63/m c m" ," P 63/m m c", &                ! CUBIC SPACE GROUPS
        " P 2 3     " ," F 2 3     " ," I 2 3     " ," P 21 3    ", &
        " I 21 3    " ," P m 3     " ," P n 3     " ," F m 3     ", &
        " F d 3     " ," I m 3     " ," P a 3     " ," I a 3     ", &
        " P 4 3 2   " ," P 42 3 2  " ," F 4 3 2   " ," F 41 3 2  ", &
        " I 4 3 2   " ," P 43 3 2  " ," P 41 3 2  " ," I 41 3 2  ", &
        " P -4 3 m  " ," F -4 3 m  " ," I -4 3 m  " ," P -4 3 n  ", &
        " F -4 3 c  " ," I -4 3 d  " ," P m 3 m   " ," P n 3 n   ", &
        " P m 3 n   " ," P n 3 m   " ," F m 3 m   " ," F m 3 c   ", &
        " F d 3 m   " ," F d 3 c   " ," I m 3 m   " ," I a 3 d   ", & ! TRIGONAL GROUPS RHOMBOHEDRAL SETTING
        " R 3   |146" ," R -3  |148" ," R 3 2 |155" ," R 3 m |160", &
        " R 3 c |161" ," R -3 m|166" ," R -3 c|167"/)


!>  SYM_GL	encoded generator strings
character(40),parameter :: SYM_GL(237)= (/  &
"000                                     ","100                                     ","01cOOO0                                 ", &
"01cODO0                                 ","02aDDOcOOO0                             ","01jOOO0                                 ", &
"01jOOD0                                 ","02aDDOjOOO0                             ","02aDDOjOOD0                             ", &
"11cOOO0                                 ","11cODO0                                 ","12aDDOcOOO0                             ", &
"11cOOD0                                 ","11cODD0                                 ","12aDDOcOOD0                             ", &
"02bOOOcOOO0                             ","02bOODcOOD0                             ","02bOOOcDDO0                             ", &
"02bDODcODD0                             ","03aDDObOODcOOD0                         ","03aDDObOOOcOOO0                         ", &
"04aODDaDODbOOOcOOO0                     ","03aDDDbOOOcOOO0                         ","03aDDDbDODcODD0                         ", &
"02bOOOjOOO0                             ","02bOODjOOD0                             ","02bOOOjOOD0                             ", &
"02bOOOjDOO0                             ","02bOODjDOO0                             ","02bOOOjODD0                             ", &
"02bDODjDOD0                             ","02bOOOjDDO0                             ","02bOODjDDO0                             ", &
"02bOOOjDDD0                             ","03aDDObOOOjOOO0                         ","03aDDObOODjOOD0                         ", &
"03aDDObOOOjOOD0                         ","03aODDbOOOjOOO0                         ","03aODDbOOOjODO0                         ", &
"03aODDbOOOjDOO0                         ","03aODDbOOOjDDO0                         ","04aODDaDODbOOOjOOO0                     ", &
"04aODDaDODbOOOjBBB0                     ","03aDDDbOOOjOOO0                         ","03aDDDbOOOjDDO0                         ", &
"03aDDDbOOOjDOO0                         ","12bOOOcOOO0                             ","03bOOOcOOOhDDD1BBB                      ", &
"12bOOOcOOD0                             ","03bOOOcOOOhDDO1BBO                      ","12bDOOcOOO0                             ", &
"12bDOOcDDD0                             ","12bDODcDOD0                             ","12bDOOcOOD0                             ", &
"12bOOOcDDO0                             ","12bDDOcODD0                             ","12bOODcODD0                             ", &
"12bOOOcDDD0                             ","03bOOOcDDOhDDO1BBO                      ","12bDDDcOOD0                             ", &
"12bDODcODD0                             ","12bDODcODO0                             ","13aDDObOODcOOD0                         ", &
"13aDDObODDcODD0                         ","13aDDObOOOcOOO0                         ","13aDDObOOOcOOD0                         ", &
"13aDDObODOcODO0                         ","04aDDObDDOcOOOhODD1OBB                  ","14aODDaDODbOOOcOOO0                     ", &
"05aODDaDODbOOOcOOOhBBB1ZZZ              ","13aDDDbOOOcOOO0                         ","13aDDDbOOOcDDO0                         ", &
"13aDDDbDODcODD0                         ","13aDDDbODOcODO0                         ","02bOOOgOOO0                             ", &
"02bOODgOOB0                             ","02bOOOgOOD0                             ","02bOODgOOF0                             ", &
"03aDDDbOOOgOOO0                         ","03aDDDbDDDgODB0                         ","02bOOOmOOO0                             ", &
"03aDDDbOOOmOOO0                         ","12bOOOgOOO0                             ","12bOOOgOOD0                             ", &
"03bOOOgDDOhDDO1YBO                      ","03bOOOgDDDhDDD1YYY                      ","13aDDDbOOOgOOO0                         ", &
"04aDDDbDDDgODBhODB1OYZ                  ","03bOOOgOOOcOOO0                         ","03bOOOgDDOcDDO0                         ", &
"03bOODgOOBcOOO0                         ","03bOODgDDBcDDB0                         ","03bOOOgOODcOOO0                         ", &
"03bOOOgDDDcDDD0                         ","03bOODgOOFcOOO0                         ","03bOODgDDFcDDF0                         ", &
"04aDDDbOOOgOOOcOOO0                     ","04aDDDbDDDgODBcDOF0                     ","03bOOOgOOOjOOO0                         ", &
"03bOOOgOOOjDDO0                         ","03bOOOgOODjOOD0                         ","03bOOOgDDDjDDD0                         ", &
"03bOOOgOOOjOOD0                         ","03bOOOgOOOjDDD0                         ","03bOOOgOODjOOO0                         ", &
"03bOOOgOODjDDO0                         ","04aDDDbOOOgOOOjOOO0                     ","04aDDDbOOOgOOOjOOD0                     ", &
"04aDDDbDDDgODBjOOO0                     ","04aDDDbDDDgODBjOOD0                     ","03bOOOmOOOcOOO0                         ", &
"03bOOOmOOOcOOD0                         ","03bOOOmOOOcDDO0                         ","03bOOOmOOOcDDD0                         ", &
"03bOOOmOOOjOOO0                         ","03bOOOmOOOjOOD0                         ","03bOOOmOOOjDDO0                         ", &
"03bOOOmOOOjDDD0                         ","04aDDDbOOOmOOOjOOO0                     ","04aDDDbOOOmOOOjOOD0                     ", &
"04aDDDbOOOmOOOcOOO0                     ","04aDDDbOOOmOOOcDOF0                     ","13bOOOgOOOcOOO0                         ", &
"13bOOOgOOOcOOD0                         ","04bOOOgOOOcOOOhDDO1YYO                  ","04bOOOgOOOcOOOhDDD1YYY                  ", &
"13bOOOgOOOcDDO0                         ","13bOOOgOOOcDDD0                         ","04bOOOgDDOcDDOhDDO1YBO                  ", &
"04bOOOgDDOcDDDhDDO1YBO                  ","13bOOOgOODcOOO0                         ","13bOOOgOODcOOD0                         ", &
"04bOOOgDDDcOODhDDD1YBY                  ","04bOOOgDDDcOOOhDDD1YBY                  ","13bOOOgOODcDDO0                         ", &
"13bOOOgDDDcDDD0                         ","04bOOOgDDDcDDDhDDD1YBY                  ","04bOOOgDDDcDDOhDDD1YBY                  ", &
"14aDDDbOOOgOOOcOOO0                     ","14aDDDbOOOgOOOcOOD0                     ","05aDDDbDDDgODBcDOFhODB1OBZ              ", &
"05aDDDbDDDgODBcDOBhODB1OBZ              ","01nOOO0                                 ","01nOOC0                                 ", &
"01nOOE0                                 ","02aECCnOOO0                             ","11nOOO0                                 ", &
"12aECCnOOO0                             ","02nOOOfOOO0                             ","02nOOOeOOO0                             ", &
"02nOOCfOOE0                             ","02nOOCeOOO0                             ","02nOOEfOOC0                             ", &
"02nOOEeOOO0                             ","03aECCnOOOeOOO0                         ","02nOOOkOOO0                             ", &
"02nOOOlOOO0                             ","02nOOOkOOD0                             ","02nOOOlOOD0                             ", &
"03aECCnOOOkOOO0                         ","03aECCnOOOkOOD0                         ","12nOOOfOOO0                             ", &
"12nOOOfOOD0                             ","12nOOOeOOO0                             ","12nOOOeOOD0                             ", &
"13aECCnOOOeOOO0                         ","13aECCnOOOeOOD0                         ","02nOOObOOO0                             ", &
"02nOOCbOOD0                             ","02nOOEbOOD0                             ","02nOOEbOOO0                             ", &
"02nOOCbOOO0                             ","02nOOObOOD0                             ","02nOOOiOOO0                             ", &
"12nOOObOOO0                             ","12nOOObOOD0                             ","03nOOObOOOeOOO0                         ", &
"03nOOCbOODeOOC0                         ","03nOOEbOODeOOE0                         ","03nOOEbOOOeOOE0                         ", &
"03nOOCbOOOeOOC0                         ","03nOOObOODeOOO0                         ","03nOOObOOOkOOO0                         ", &
"03nOOObOOOkOOD0                         ","03nOOObOODkOOD0                         ","03nOOObOODkOOO0                         ", &
"03nOOOiOOOkOOO0                         ","03nOOOiOODkOOD0                         ","03nOOOiOOOeOOO0                         ", &
"03nOOOiOODeOOO0                         ","13nOOObOOOeOOO0                         ","13nOOObOOOeOOD0                         ", &
"13nOOObOODeOOD0                         ","13nOOObOODeOOO0                         ","03bOOOcOOOdOOO0                         ", &
"05aODDaDODbOOOcOOOdOOO0                 ","04aDDDbOOOcOOOdOOO0                     ","03bDODcODDdOOO0                         ", &
"04aDDDbDODcODDdOOO0                     ","13bOOOcOOOdOOO0                         ","04bOOOcOOOdOOOhDDD1YYY                  ", &
"15aODDaDODbOOOcOOOdOOO0                 ","06aODDaDODbOOOcOOOdOOOhBBB1ZZZ          ","14aDDDbOOOcOOOdOOO0                     ", &
"13bDODcODDdOOO0                         ","14aDDDbDODcODDdOOO0                     ","04bOOOcOOOdOOOeOOO0                     ", &
"04bOOOcOOOdOOOeDDD0                     ","06aODDaDODbOOOcOOOdOOOeOOO0             ","06aODDaDODbODDcDDOdOOOeFBF0             ", &
"05aDDDbOOOcOOOdOOOeOOO0                 ","04bDODcODDdOOOeBFF0                     ","04bDODcODDdOOOeFBB0                     ", &
"05aDDDbDODcODDdOOOeFBB0                 ","04bOOOcOOOdOOOlOOO0                     ","06aODDaDODbOOOcOOOdOOOlOOO0             ", &
"05aDDDbOOOcOOOdOOOlOOO0                 ","04bOOOcOOOdOOOlDDD0                     ","06aODDaDODbOOOcOOOdOOOlDDD0             ", &
"05aDDDbDODcODDdOOOlBBB0                 ","14bOOOcOOOdOOOeOOO0                     ","05bOOOcOOOdOOOeOOOhDDD1YYY              ", &
"14bOOOcOOOdOOOeDDD0                     ","05bOOOcOOOdOOOeDDDhDDD1YYY              ","16aODDaDODbOOOcOOOdOOOeOOO0             ", &
"16aODDaDODbOOOcOOOdOOOeDDD0             ","07aODDaDODbODDcDDOdOOOeFBFhBBB1ZZZ      ","07aODDaDODbODDcDDOdOOOeFBFhFFF1XXX      ", &
"15aDDDbOOOcOOOdOOOeOOO0                 ","15aDDDbDODcODDdOOOeFBB0                 ","01dOOO0                                 ", &
"11dOOO0                                 ","02dOOOfOOO0                             ","02dOOOlOOO0                             ", &
"02dOOOlDDD0                             ","12dOOOfOOO0                             ","12dOOOfDDD0                             "/) 



!>  SGPG contains the first space group # for a given point group
integer(kind=irg),parameter :: SGPG(32) =(/1,2,3,6,10,16,25,47,75,81,83,89,99,111,123,143, &
                                          147,149,156,162,168,174,175,177,183,187,191,195, &
                                          200,207,215,221/)

!>  SGsym contains the numbers of all the symmorphic space groups
integer(kind=irg),parameter :: SGsym(73) =(/1,2,3,5,6,8,10,12,16,21,22,23,25,35,38,42,44,47, &
                                            65,69,71,75,79,81,82,83,87,89,97,99,107,111,115, &
                                            119,121,123,139,143,146,147,148,149,150,155,156, &
                                            157,160,162,164,165,168,174,175,177,183,187,189, &
                                            191,195,196,197,200,202,204,207,209,211,215,216, &
                                            217,221,225,229/)

! these parameters implement the diffraction group
! formalism described in the BESR paper.

!> 10 2D point group symbols in International Tables order
character(10),parameter  :: PGTWD(0:11) = (/ ' none     ','    1     ','    2     ','    m     ','  2mm     ','    4     ', &
                                             '  4mm     ','    3     ','   3m1    ','    6     ','  6mm     ','   31m    ' /)

!> 10 2D point group orders in International Tables order
integer(kind=irg),parameter       :: PGTWDorder(0:11) = (/0,1,2,2,4,4,8,3,6,6,12,6/)

!> inverse table for 2D point groups; this essentially implements the inverse of Table 4 in BESR paper for the Bright Field symmetry.
integer(kind=irg),parameter       :: PGTWDinverse(12,11) = reshape((/ & 
                                   1,0,0,0,0,0,0,0,0,0,0,0,  1,2,0,0,0,0,0,0,0,0,0,0, &
                                   1,3,0,4,0,0,0,0,0,0,0,0,  1,3,0,5,0,0,0,0,0,0,0,0, &
                                   1,3,0,4,0,0,0,6,0,0,0,0,  1,0,7,0,0,0,0,0,0,0,0,0, &
                                   1,2,0,0,0,8,0,0,0,0,0,0,  1,3,0,0,0,9,0,0,0,0,0,0, &
                                   1,3,0,4,0,0,0,0,0,0,0,10, 1,3,7,4,0,0,0,0,0,0,0,0, &
                                   1,3,0,4,0,8,0,6,0,0,0,0 /), (/ 12,11 /))


!> 32 3D point group symbols in International Tables order
character(5),parameter  :: PGTHD(32) =(/'    1','   -1','    2','    m','  2/m','  222', &
                                        '  mm2','  mmm','    4','   -4','  4/m','  422', &
                                        '  4mm',' -42m','4/mmm','    3','   -3','   32', &
                                        '   3m','  -3m','    6','   -6','  6/m','  622', &
                                        '  6mm',' -6m2','6/mmm','   23','   m3','  432', &
                                        ' -43m',' m-3m'/)

!> 3D point groups : Laue group number
integer(kind=irg),parameter       :: PGLaue(32) =(/2,2,5,5,5,8,8,8,11,11,11,15,15,15,15,17,17, &
                                                  20,20,20,23,23,23,27,27,27,27,29,29,32,32,32/)
!> 3D point groups : inverted Laue group number
integer(kind=irg),parameter       :: PGLaueinv(32) = (/1,1,2,2,2,3,3,3,4,4,4,5,5,5,5,6,6, &
                                                       7,7,7,8,8,8,9,9,9,9,10,10,11,11,11/)

!> 31 diffraction group symbols in BESR order
character(5),parameter  :: DG(31) =(/'    1','   1R','    2','   2R','  21R','   mR', &
                                     '    m','  m1R','2mRmR','  2mm','2RmmR','2mm1R', &
                                     '    4','   4R','  41R','4mRmR','  4mm','4RmmR', &
                                     '4mm1R','    3','   6R','  3mR','   3m','6RmmR', &
                                     '    6','  31R','  61R','6mRmR','  6mm',' 3m1R', &
                                     '6mm1R'/)

!> 31 diffraction group orders in BESR order
integer(kind=irg),parameter       :: DGorder(31) =(/1, 2, 2, 2, 4, 2, 2, 4, 4, 4, 4, 8, &
                                          4, 4, 8, 8, 8, 8,16, 3, 6, 6, 6,12, &
                                          6, 6,12,12,12,12,24/)

!> Bright Field planar point group for 31 diffraction groups (Table 2, column 2, BESR, with change in row ordering)
integer(kind=irg),parameter       :: BFPG(31) =(/1,2,2,1,2,3,3,4,4,4,3,4,5,5,5,6,6,6,6,7,7,8,8,8,9,9,9,10,10,10,10/)

!> Whole Pattern planar point group for 31 diffraction groups (Table 2, column 3, BESR, with change in row ordering)
integer(kind=irg),parameter       :: WPPG(31) =(/1,1,2,1,2,1,3,3,2,4,3,4,5,2,5,5,6,4,6,7,7,7,8,8,9,7,9,9,10,8,10/)

!> Dark Field planar point group for 31 diffraction groups (Table 2, column 4, BESR, with change in row ordering)
integer(kind=irg),parameter       :: DFGN(31) = (/1,2,1,1,2,1,1,2,1,1,1,2,1,1,2,1,1,1,2,1,1,1,1,1,1,2,2,1,1,2,2/)

!> Dark Field planar point group for 31 diffraction groups (Table 2, column 5, BESR, with change in row ordering)
integer(kind=irg),parameter       :: DFSP(31) = (/0,0,0,0,0,3,3,4,3,3,3,4,0,0,0,3,3,3,4,0,0,3,3,3,0,0,0,3,3,4,4/)

!> 10 projection diffraction groups in BESR order (Table 2, column 8, BESR, with change in row ordering)
integer(kind=irg),parameter       :: PDG(31) = (/2,2,5,5,5,8,8,8,12,12,12,12,15,15,15,19,19,19,19,26,27,30,30, &
                                                 31,27,26,27,31,31,30,31/)

!> short hand for .FALSE. logical parameter
logical,parameter :: F=.FALSE.

!> short hand for .TRUE. logical parameter
logical,parameter :: T=.TRUE.

!> Table 3 from BESR paper
logical,parameter       :: DGPG(32,31) = reshape((/ &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,T, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,T,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,T,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T, &
     F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F, &
     F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F, &
     F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F, &
     F,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,T,F,T,F,F,T, &
     F,F,F,F,T,F,F,T,F,F,T,F,F,F,T,F,F,F,F,T,F,F,T,F,F,F,T,F,T,F,F,T, &
     F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F, &
     F,F,F,F,F,T,F,F,F,F,F,T,F,T,F,F,F,F,F,F,F,F,F,T,F,F,F,T,F,T,F,F, &
     F,F,F,F,F,F,T,F,F,F,F,F,T,T,F,F,F,F,F,F,F,F,F,F,T,T,F,F,F,F,T,F, &
     F,F,F,T,F,F,T,F,F,F,F,F,T,T,F,F,F,F,T,F,F,T,F,F,T,T,F,F,F,F,T,F, &
     F,F,T,F,F,T,T,F,T,T,F,T,T,T,F,F,F,T,F,F,T,F,F,T,T,T,F,T,F,T,T,F, &
     F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F, &
     F,T,F,F,T,F,F,T,F,F,T,F,F,F,T,F,T,F,F,T,F,F,T,F,F,F,T,F,T,F,F,T, &
     F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F, &
     F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F, &
     T,F,T,T,F,T,T,F,T,T,F,T,T,T,F,T,F,T,T,F,T,T,F,T,T,T,F,T,F,T,T,F/), (/32,31/))


! these lines are from an older version; not sure if they are still needed...
! SYM_SGtworig	= point symmetries for two origin choices
! SYM_NUMtworig	= number of space groups with two origin settings


! declare user-defined types
!> symdata type definition (for 3D symmetry operations with space groups)
type symdata
  integer(kind=irg) 	:: SYM_GENnum			!< number of generator matrices
  integer(kind=irg) 	:: SYM_MATnum			!< number of non-zero symmetry matrices
  integer(kind=irg) 	:: SYM_NUMpt			!< number of point group operators
  logical          	:: SYM_reduce			!< switch to enable/disable reduction to fundamental cell
  logical          	:: SYM_trigonal			!< switch for hexagonal vs. rhombohedral settings
  logical          	:: SYM_second			!< switch for second setting of spacegroup (if any)
  logical          	:: SYM_centrosym		!< switch for presence of centrosymmetry
  real(kind=dbl)    	:: SYM_data(192,4,4)		!< all symmetry matrices for a given spacegroup
  real(kind=dbl)    	:: SYM_direc(48,3,3)		!< direct space point group matrices
  real(kind=dbl)    	:: SYM_recip(48,3,3)		!< reciprocal space point group matrices
  real(kind=dbl)    	:: SYM_c(4,4)			!< dummy 4x4 matrix used for various computations
  character(11)    	:: SYM_name			!< current space group name
end type

! declare global variables
!> @param SG entire space group structure
type (symdata)   	:: SG
	
! arrays used by CalcFamily, CalcPositions and related routines
integer(kind=irg)	:: itmp(48,3)			!< auxiliary array used for family computations etc
integer(kind=irg)	:: numat(maxpasym)		!< number of atoms of each type in the asymmetric unit

! atom coordinates
real(kind=dbl),allocatable :: apos(:,:,:)		!< atom position array

! for many diffraction calculations we need the 2D planar point groups; 
! the maximum order of such a group is 12, and there are only 10 of them, with
! two settings for one of them (3m1 and 31m).
type symdata2D
  integer(kind=irg)	:: SYM_pgnum			!< 2D point group number
  integer(kind=irg)	:: SYM_MATnum			!< number of non-zero symmetry matrices (order)
  integer(kind=irg)	:: SYM_direc(12,2,2)		!< point group matrices (filled in by Generate2DSymmetry)
end type

! define the 2D point group variable
type (symdata2D)	:: TDPG


end module symmetryvars

