module diagonal
   ! This module contains type diag, with many tools to manipulate them.
   implicit none
   private
   type :: diag
      double precision, allocatable, dimension(:) :: d
   end type
   !! Public type
   public :: diag
   !! Public functions/subroutines

   !! Public overload
contains
   function diag_mul_right(D,A) result(B)
      double precision, intent(in), dimension(:,:) :: A
      type(diag), intent(in) :: D
      double precision, dimension(size(A,1),size(A,2)) :: B
      ! local variables 
      integer M,I
      M=size(A,1)
      do I=1,M
         B(I,:)=D%d(I)*A(I,:)
      end do
   end function
end module diagonal
