module Constants
   implicit none
   real :: pi = 3.1415
end module

module CircleModule
  implicit none
  contains

  subroutine area(a, r)
    real :: a, r, pi = 3.141
    a = pi*r**2
  end subroutine

end module

module AddModule
  implicit none
  contains
  subroutine Sub(b)
    implicit none
    integer :: b
    integer :: a
    a = 0
    a = a + 1
    b = b + a
  end subroutine
end module

module ArrayModule
  implicit none
  contains
  subroutine subBad(A)
    real, dimension(2,2):: A
    print'(4f10.3)',transpose(A)
  end subroutine

  subroutine subNice(A)
    real, dimension(:,:):: A
    print'(4f10.3)',transpose(A)
  end subroutine

end module
