! !program AddTwoNumbers
! !   implicit none

!    ! Variable declarations
!    integer :: a, b, c

!    ! Add two numbers
!    a = 5
!    b = 10
!    c = a + b

!    ! Print result to screen
!    print*,'The sum is ', c

! end program


! program IntegerPitfall
!   implicit none

!   ! Variable declaration
!   real :: x
!   x = 1.0

!   ! Print results to screen
!   print*, 2.0/3.0

! end program


! program IfStatements
!    integer :: a
!    real :: b
!    a = -2
!    b = 0.2
!    if (a == 1) then
!       print* , 'a is equal to 1'
!    else if (a < 0) then
!       if (b > 0.1) then
!            print* , 'a is smaller than 0 and b is larger to 0.1'
!       end if
!    else
!        print*, 'a is larger or equal to 0'
!    end if
! end program

! program DoLoops
!    integer:: i
!    do i =  0, 10
!        if (o==5) then

!            exit
!        end if
!        print*, i
!    end do
! end program


! program PrintTemp
!  implicit none
!  real :: T = 21.5

!  ! Print using default formatting
!  print*, T

!  ! Print using a format specifier
!  print'(f5.2)', T

!  ! Print using default formatting
!  print*, 'Temp. is ', T

!  ! Print using a format specifier
!  print'(Af5.2)', 'Temp. is ', T

! end program


!program PrintNum
!
!  print*, 12.345
!
!end program

!program MainProgram
!   use Constants
!   implicit none
!   real :: area, radius
!
!   radius = 2
!   area = pi*radius**2
!   print*,'Area = ', area
!end program

! program MainProgram
!  use CircleModule
!  implicit none
!  real :: a

!  call area(a, 5.0)
!  print*,'Area = ', a
! end program


! program MainProgram
!  use AddModule
!  implicit none
!  integer :: b = 2
!  print*, 'b = ', b
!  call Sub(b)
!  print*, 'b = ', b
!  call Sub(b)
!  print*, 'b = ', b
! end program

!program MainProgram
!  use ArrayModule
!  real, dimension(3,3):: A
!  A(1,:) = (/ 1 ,1 ,1 ,1 /)
!  A(2,:) = (/ 2 ,2 ,2 ,2 /)
!  A(3,:) = (/ 3 ,3 ,3 ,3 /)
!  print'(4f10.3)',transpose(A)
!  print* ,' '
!  call subBad(A)
!  print* ,' '
!  call subNice(A)
!end program

! program MainProgram
!     integer :: a = 2
!     integer :: b = 3
!     print*, 'a = ', a
!     print*, 'b = ', b
!     print*, 'a * b = ', a*b
! end program