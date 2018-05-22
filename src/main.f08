PROGRAM TEST
   use iter_sys
   use easy_mat_mod
   implicit none
   double precision, allocatable, dimension(:,:) :: A
   double precision, allocatable, dimension(:) :: B,X
   type(diag) :: C,D,E
   class(split), allocatable :: S
   double precision :: loc_err
   integer :: ENDD, L, LL, N0, I, J, K, M, N, P, R, II, MEM
   M=10
   max_iter=100
   allocate(A(M,M),B(M))!,X(M))
   do J=1,M
      A(:,J)=[(I+J,I=1,M)]
      A(J,J)=2*1e3
   end do 
   B=1
!   X=jacob(A,B)
!   write(*,*) 'vediamo se termina', norminf(matmul(A,X)-B)
   D=diag([(2*I,I=1,M)])
   C=D
   D=inv(C)
   E=diag(C,D)
   A=diag_full(D)
   allocate(S)
!   allocate(S%MM,source=C%el)
!   S%MM=C
!   S%NN=A
   S%E=D
!   E=matmul(C,D)
!   write(*,*) C%el
!   write(*,*) D%el
   write(*,*) E%el
END PROGRAM
