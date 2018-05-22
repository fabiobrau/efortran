module iter_sys
   !! In this module there are iterative metods for solving linear system
   use easy_mat_mod
   implicit none
   private
   double precision :: tol=1d-10
   integer :: max_iter=100000
   !! public overload
   
   !! public function/subroutines
   public :: jacob
   !! public variables
   public :: tol
   public :: max_iter
contains
   function  jacob(A,B,X_0) result(X) 
      ! Take in input A, B e return the solution of system AX=B
      ! A, double precision, dim(A)=MxM
      ! B, double precision, dim(B)=M
      ! X, double precision, dim(X)=M
      ! loc_err, double precision, scalar
      implicit none
      double precision, intent(in), dimension(:,:) :: A
      double precision, intent(in), dimension(:) :: B
      double precision, optional, intent(in), dimension(:) :: X_0
      double precision, dimension(size(A,2)) :: X
      ! local variables
      integer :: M, K, I
      double precision :: loc_err=1
      double precision, allocatable, dimension(:) ::Xold, Baux
      double precision, allocatable, dimension(:,:) :: N
      type(diag) :: D
      M=size(A,1)
      N=-A
      if (present(X_0)) then
         Xold=X_0
      else
         Xold=B;
      end if
      do K=1,M
         N(K,K)=0
      end do
      D=diag([(1/A(I,I),I=1,M)])
      Baux=matmul(D,B)
      do K=1,max_iter
         X=matmul(D,matmul(N,Xold))+ Baux
         loc_err=norminf(X-Xold)
         if (loc_err<tol) then
            exit
         end if
         Xold=X
      end do
    end function
end module
