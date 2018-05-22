module easy_mat_mod
   ! This module contains various type to execute operations in optimized way
   implicit none
   private
   type :: diag
      ! diagonal square matrix
      double precision, allocatable, dimension(:) :: el
   end type
   type :: split
      class(*), allocatable :: MM,NN ! Such that A=MM-NN
      type(diag) :: E 
   contains
      procedure :: split_size
   end type
   type :: msplit
      integer K ! number of splitting
      type(split), allocatable, dimension(:) :: S
   end type
   interface matmul
      module procedure diag_mul_mat
      module procedure diag_mul_vec
      module procedure mat_mul_diag
      module procedure diag_mul_diag
   end interface
   interface size
      module procedure diag_size
   end interface
   interface diag
      module procedure diag_merge_diag
      module procedure diag_merge_vec
      module procedure vec_merge_diag
   end interface
   !! Public class/type
   public :: diag
   public :: split
   public :: msplit
   !! Public functions/subroutines
   public :: norminf
   public :: inv
   public :: diag_full
   !! Public overload
   public :: matmul
   public :: size
contains
   !! Tools for diagonal type
   function diag_mul_mat(A,B) result(C)
      ! left product bay diagonal type
      type(diag) A
      double precision, intent(in), dimension(:,:) :: B
      double precision, dimension(size(B,1),size(B,2)) :: C
      ! local variables 
      integer M,I
      M=size(B,1)
      do I=1,M
         C(I,:)=A%el(I)*B(I,:)
      end do
   end function
   function diag_mul_vec(A,B) result(C)
      ! left product by diagonal type with vector
      type(diag) A
      double precision, intent(in), dimension(:) :: B
      double precision, dimension(size(B,1)) :: C
      ! local variables  
      C=A%el*B
   end function
   function mat_mul_diag(A,B) result(C)
      ! right product of diagonal with matrix
      double precision, intent(in), dimension(:,:) :: A
      type(diag), intent(in) :: B
      double precision, dimension(size(A,1),size(B)) :: C
      ! local variables 
      integer M,I
      M=size(A,1)
      do I=1,M
         C(:,I)=B%el(I)*A(:,I)
      end do
   end function
   function diag_mul_diag(A,B) result(C)
      type(diag), intent(in) :: A,B
      type(diag) :: C
      ! local variables 
      C=diag(A%el*B%el)
   end function
   pure function diag_size(D) result (M)
      type(diag), intent(in) :: D
      integer M
      M=size(D%el)
   end function
   pure function norminf(X) result (a)
       ! Take in input X and return the infinity norm of X
       ! X, double precision, dim(X)=m,
       implicit none
       double precision, intent(in), dimension(:) :: X
       double precision a
       a=maxval(abs(X))
   end function
   pure function inv(D) result(B)
      ! Invert matrix of type diagonal
      implicit none
      type(diag), intent(in) :: D
      type(diag) B
      integer I
      B=diag([(1/D%el(I),I=1,size(D))])
   end function
   pure function diag_merge_diag(A,B) result(C)
      type(diag), intent(in) :: A,B
      type(diag :: C
      ! local variable
      C=diag([A%el,B%el])
   end function
   pure function diag_merge_vec(A,B) result(C)
      type(diag), intent(in) :: A
      double precision, intent(in) :: B
      type(diag) :: C
      C=diag([A%el,B])
   end function
   pure function vec_merge_diag(A,B) result(C)
      double precision, intent(in) :: A
      type(diag), intent(in) :: B
      type(diag) :: C
      C=diag([A,B%el])
   end function
   pure function diag_full(D) result(F)
      ! return full matrix from diagonal type
      type(diag), intent(in) :: D
      double precision, dimension(size(D),size(D)) :: F
      ! local variables
      integer I
      F=0
      do I=1,size(D)
         F(I,I)=D%el(I)
      end do
   end function
   pure function split_size(S) result(J)
     class(split), intent(in) :: S
     integer J
     J=size(S%E)
   end function 
end module easy_mat_mod
