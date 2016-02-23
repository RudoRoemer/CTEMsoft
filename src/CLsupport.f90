! ###################################################################
! Copyright (c) 2013-2016, Marc De Graef Research Group/Carnegie Mellon University
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
module CLsupport

use local
use clfortran

IMPLICIT NONE


contains


!--------------------------------------------------------------------------
!
! SUBROUTINE:CLquery_platform_info
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief display information about an OpenCL device (based on (clfortran's query_platforms_devices.f90)
!
!> @param platform_id id number of platform
!
!> @date 02/18/16  MDG 1.0 modification of clfortran's original routine
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! original Copyright information (clfortran's query_platforms_devices.f90)
! -----------------------------------------------------------------------------
!
! Copyright (C) 2013-2014 Company for Advanced Supercomputing Solutions LTD
! Bosmat 2a St.
! Shoham
! Israel 60850
! http://www.cass-hpc.com
!
! Author: Mordechai Butrashvily <support@cass-hpc.com>
!
! -----------------------------------------------------------------------------
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
subroutine CLquery_platform_info(platform_id)

use ISO_C_BINDING

IMPLICIT NONE

! Input variable.
integer(c_intptr_t), INTENT(IN):: platform_id

! Helper variables to work with OpenCL API.
integer(c_int32_t)             :: err
integer(c_size_t)              :: zero_size = 0
integer(c_size_t)              :: temp_size
! For quering devices.
integer(c_int64_t)             :: device_type
integer(c_int32_t)             :: num_devices
integer(c_int)                 :: i
integer(c_intptr_t), allocatable, target :: device_ids(:)

! String arrays for holding platform details.
character, allocatable, target :: platform_profile(:)
character, allocatable, target :: platform_version(:)
character, allocatable, target :: platform_name(:)
character, allocatable, target :: platform_vendor(:)
character, allocatable, target :: platform_extensions(:)

! String array for holding device name.
character, allocatable, target :: device_name(:)
! Maximum compute units for device.
integer(c_int32_t), target     :: device_cu

! Profile.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_PROFILE, zero_size, C_NULL_PTR, temp_size)
allocate(platform_profile(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_PROFILE, temp_size, C_LOC(platform_profile), temp_size)
print *, 'Profile: ', platform_profile
deallocate(platform_profile)

! Version.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_VERSION, zero_size, C_NULL_PTR, temp_size)
allocate(platform_version(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_VERSION, temp_size, C_LOC(platform_version), temp_size)
print *, 'Version: ', platform_version
deallocate(platform_version)

! Name.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_NAME, zero_size, C_NULL_PTR, temp_size)
allocate(platform_name(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_NAME, temp_size, C_LOC(platform_name), temp_size)
print *, 'Name: ', platform_name
deallocate(platform_name)

! Vendor.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_VENDOR, zero_size, C_NULL_PTR, temp_size)
allocate(platform_vendor(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_VENDOR, temp_size, C_LOC(platform_vendor), temp_size)
print *, 'Vendor: ', platform_vendor
deallocate(platform_vendor)

! Extensions.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_EXTENSIONS, zero_size, C_NULL_PTR, temp_size)
allocate(platform_extensions(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_EXTENSIONS, temp_size, C_LOC(platform_extensions), temp_size)
print *, 'Extensions: ', platform_extensions
deallocate(platform_extensions)

!
! Print device information for this platform.
!
! Get device count.
!device_type = CL_DEVICE_TYPE_ALL
err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, 0, C_NULL_PTR, num_devices)

if (err /= CL_SUCCESS .or. num_devices < 1) then
  print *, 'No devices found: ', err
  return
end if

print *
print '(A, I2)', 'Num Devices: ', num_devices

! Allocate an array to hold device handles.
allocate(device_ids(num_devices))

! Get device IDs.
err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, num_devices, C_LOC(device_ids), num_devices)
if (err /= CL_SUCCESS) then
  print *, 'Error quering devices: ', err
  return
end if

! Loop over devices and print information.
do i = 1, num_devices
! Maximum compute units.
  temp_size = 4
  err = clGetDeviceInfo(device_ids(i), CL_DEVICE_MAX_COMPUTE_UNITS, temp_size, C_LOC(device_cu), temp_size)

! Name.
  err = clGetDeviceInfo(device_ids(i), CL_DEVICE_NAME, zero_size, C_NULL_PTR, temp_size)
  allocate(device_name(temp_size))
  err = clGetDeviceInfo(device_ids(i), CL_DEVICE_NAME, temp_size, C_LOC(device_name), temp_size)

! Print brief device details. Since this routine can be used by the user to determine GPU device IDs,
! we subtract 1 from the device ID to make sure that the CPU gets number 0...
  write (*, '(A,I2,A,I3,A)', advance='no') ' Device (#', i-1, ', Compute Units: ', device_cu, ') - '
  print *, device_name

  deallocate(device_name)
end do
end subroutine CLquery_platform_info


!--------------------------------------------------------------------------
!
! SUBROUTINE:CLread_source_file
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read an OpenCL source file and return the source properly formatted
!
!> @param sourcefile filename for the OpenCL source code
!> @param source c_str containing the source, NULL-terminated
!> @param slength source string length
!
!> @date 02/18/16  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine CLread_source_file(sourcefile, csource, slength)

use local
use error
use ISO_C_BINDING

IMPLICIT NONE

integer, parameter                      :: source_length = 50000

character(fnlen), INTENT(IN)            :: sourcefile
character(len=source_length, KIND=c_char),INTENT(OUT) :: csource
integer(c_size_t),INTENT(OUT)           :: slength


character(len=source_length),target     :: source
character(fnlen)                        :: fname
integer(kind=irg)                       :: irec, ierr 

! read the source file in the EMsoft/opencl folder
fname = trim(openclpathname)//trim(sourcefile)
fname = EMsoft_toNativePath(fname)
open(unit = dataunit, file = trim(fname), access='direct', status = 'old', &
     action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) call FatalError("CLread_source_file: ",'Cannot open file '//fname)

source = ''
irec = 1
do
  read(unit = dataunit, rec = irec, iostat = ierr) source(irec:irec)
  if (ierr /= 0) exit
  if(irec == source_length) call FatalError("CLread_source_file: ",'Error: CL source file is too big')
  irec = irec + 1
end do
close(unit=dataunit)

csource = trim(source)
csource(irec:irec) = C_NULL_CHAR
slength = irec

end subroutine CLread_source_file



end module CLsupport
