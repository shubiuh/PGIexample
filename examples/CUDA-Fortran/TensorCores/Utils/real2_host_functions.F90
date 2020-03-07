!
!     Copyright (c) 2019, NVIDIA CORPORATION.  All rights reserved.
!
! NVIDIA CORPORATION and its licensors retain all intellectual property
! and proprietary rights in and to this software, related documentation
! and any modifications thereto.  Any use, reproduction, disclosure or
! distribution of this software and related documentation without an express
! license agreement from NVIDIA CORPORATION is strictly prohibited.
!
!          THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT
!   WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT
!   NOT LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR
!   FITNESS FOR A PARTICULAR PURPOSE.
!
module real2_host_functions

  use pgi_acc_common

  interface assignment(=)
    module procedure hosthalf2float
    module procedure hostfloat2half
  end interface

  interface random_number
    module procedure rando_real2
  end interface

  interface __half2float
    module procedure __half2float1d
    module procedure __half2float2d
  end interface

  interface __float2half
    module procedure __float2half1d
    module procedure __float2half2d
  end interface

#if !defined OPENPOWER && !defined _WIN64
  interface
    integer function supportsHalfConvert() bind(C,name="__Cpuid_is_f16c")
    end function
  end interface
#endif

contains

! 1-D
function __half2float1d(a) result(res)
use pgi_acc_common, only : WMMAHALF
type(WMMAHALF) :: a(:)
real(4) :: res(size(a))
n = size(a)
if (n.eq.1) then
  call hosthalf2float(res(1),a(1))
else if (n.gt.1) then
  j = 1
#if !defined OPENPOWER && !defined _WIN64
  if (supportsHalfConvert() == 1) then
    if ((n-j+1).ge.8) then
      do i = j, n, 8
        call v8half2float(res(i),a(i))
      end do
      j = j + 8 * (n/8)
    end if
    if ((n-j+1).ge.4) then
      call v4half2float(res(j),a(j))
      j = j + 4
    end if
  end if
#endif
  do i = j, n
    call hosthalf2float(res(i),a(i))
  end do
end if
end function
  
function __float2half1d(a) result(res)
use pgi_acc_common, only : WMMAHALF
real(4) :: a(:)
type(WMMAHALF) :: res(size(a))
n = size(a)
if (n.eq.1) then
  call hostfloat2half(res(1),a(1))
else if (n.gt.1) then
  j = 1
#if !defined OPENPOWER && !defined _WIN64
  if (supportsHalfConvert() == 1) then
    if ((n-j+1).ge.8) then
      do i = j, n, 8
        call v8float2half(res(i),a(i))
      end do
      j = j + 8 * (n/8)
    end if
    if ((n-j+1).ge.4) then
      call v4float2half(res(j),a(j))
      j = j + 4
    end if
  end if
#endif
  do i = j, n
    call hostfloat2half(res(i),a(i))
  end do
end if
end function
  
! 2-D
function __half2float2d(a) result(res)
use pgi_acc_common, only : WMMAHALF
type(WMMAHALF) :: a(:,:)
real(4) :: res(size(a,dim=1),size(a,dim=2))
logical fastConvert
n1 = size(a,dim=1)
n2 = size(a,dim=2)
#if !defined OPENPOWER && !defined _WIN64
fastConvert = supportsHalfConvert() .eq. 1
#endif
do k = 1, n2
  if (n1.ge.1) then
    j = 1
#if !defined OPENPOWER && !defined _WIN64
    if (fastConvert) then
      if ((n1-j+1).ge.8) then
        do i = j, n1, 8
          call v8half2float(res(i,k),a(i,k))
        end do
        j = j + 8 * (n1/8)
      end if
      if ((n1-j+1).ge.4) then
        call v4half2float(res(j,k),a(j,k))
        j = j + 4
      end if
    end if
#endif
    do i = j, n1
      call hosthalf2float(res(i,k),a(i,k))
    end do
  end if
end do
end function

function __float2half2d(a) result(res)
use pgi_acc_common, only : WMMAHALF
real(4) :: a(:,:)
type(WMMAHALF) :: res(size(a,dim=1),size(a,dim=2))
logical fastConvert
n1 = size(a,dim=1)
n2 = size(a,dim=2)
#if !defined OPENPOWER && !defined _WIN64
fastConvert = supportsHalfConvert() .eq. 1
#endif
do k = 1, n2
  if (n1.ge.1) then
    j = 1
#if !defined OPENPOWER && !defined _WIN64
    if (fastConvert) then
      if ((n1-j+1).ge.8) then
        do i = j, n1, 8
          call v8float2half(res(i,k),a(i,k))
        end do
        j = j + 8 * (n1/8)
      end if
      if ((n1-j+1).ge.4) then
        call v4float2half(res(j,k),a(j,k))
        j = j + 4
      end if
    end if
#endif
    do i = j, n1
      call hostfloat2half(res(i,k),a(i,k))
    end do
  end if
end do
end function
  
subroutine rando_real2(a, scale, offset)
type(WMMAHALF), intent(out) :: a(:)
integer(2), dimension(size(a)) :: ia
real(4), optional :: scale, offset
real(4)  ax(100), scale1, offset1
integer(4)  ix(100)
equivalence(ax,ix)
n = size(a)
if (present(scale)) then
  scale1 = scale
else
  scale1 = 1.0
end if
if (present(offset)) then
  offset1 = offset
else
  offset1 = 0.0
end if
do i = 1, n, 100
  call random_number(ax)
  do j = i, min(n,i+99)
    k = j-i+1
    ax(k) = ax(k)*scale1+offset1
    ia(j) = ior(z'3c00',ishft(iand(z'7fe000',ix(k)),-13))
  end do
end do
a = transfer(ia,a)
return
end

subroutine hostfloat2half(yt,x)
use pgi_acc_common, only : WMMAHALF
use ieee_arithmetic
type(WMMAHALF), intent(out) :: yt
real(4), intent(in) :: x
integer(2) :: y
if (ieee_is_nan(x)) then
  y = z'7fff'
else if (x.eq.0.0) then
  itmp1 = transfer(x,itmp1)
  if (itmp1.lt.0) then
    y = z'8000' ! -0.0
  else
    y = z'0000' ! +0.0
  end if
else if (x.ge.65520.0) then  ! overflow
  y = z'7c00' ! +infinity
else if (x.le.-65520.0) then
  y = z'fc00' ! -infinity
else if (abs(x).lt.(2.0**-14.0)) then
  ! denormal
  itmp1 = transfer(x,itmp1)
  ! Use current rounding mode
  itmp2 = ieee_rint(abs(x)*(2.0**24.0))
  y = iand(ishft(itmp1,-16),z'08000') + itmp2
else
  itmp = transfer(x,itmp)
  itmp1 = ishft(itmp,1)
  itmp1 = ishft(itmp1,-1)
  iman = iand(itmp1,z'7fffff')
  iman = ieee_rint(real(iman)*(2.0**-13.0))
  iexp = (ishft(itmp1,-23) - 127) + 15
  y = iand(ishft(itmp,-16),z'08000') + ishft(iexp,10) + iman
  yt%X = y
endif
end subroutine

subroutine hosthalf2float(y,xt)
use pgi_acc_common, only : WMMAHALF
real(4), intent(out) :: y
type(WMMAHALF), intent(in) :: xt
integer(2) :: x
x = xt%X
itmp = iand(x,z'7fff')
iexp = ishft(itmp,-10)
isgn = iand(x,z'8000')
iman = iand(x,z'03ff')
if (iexp .eq. 0) then
  if (isgn .eq. 0) then
    y = real(iman) * 2.0**-24.0  ! Should be exact
  else
    y = -(real(iman) * 2.0**-24.0)
  end if
else
  if (iexp .eq. 31) then
    iexp = z'00ff'  ! Infinity or NaN
  else
    iexp = (iexp - 15) + 127
  end if
  itmp = ior(ior(ishft(isgn,16),ishft(iexp,23)),ishft(iman,13))
  y = transfer(itmp,y)
end if
end subroutine
end module real2_host_functions
