MODULE fastmul
        implicit none
        private
! make moltiplication more easy
        interface operator(*)
        !        module procedure MATMATMUL !CONFLITTO CON PRODOTTO TIPO .* DIO CANE
                module procedure MATVECMUL
                module procedure VECMATMUL
        end interface
        public :: operator(*)
contains
        PURE FUNCTION MATMATMUL(A,B) RESULT(C)
        ! ONLY FOR DOUBLE PRECISION MATRIX
        ! array(2) x array(2) multiplier
                implicit none
                double precision, dimension(:,:), intent(IN) :: A, B
                double precision, dimension(size(A,1),size(B,2)) :: C
                C=matmul(A,B)
        END FUNCTION MATMATMUL
        PURE FUNCTION MATVECMUL(A,B) RESULT(C)
        ! ONLY FOR DOUBLE PRECISION MATRIX
        ! array(2) x array(1) multiplier
                implicit none
                double precision, dimension(:,:), intent(IN) :: A
                double precision, dimension(:), intent(IN) :: B
                double precision, dimension(size(A,1)) :: C
                C=matmul(A,B)
        END FUNCTION MATVECMUL
        PURE FUNCTION VECMATMUL(A,B) RESULT(C)
        ! ONLY FOR DOUBLE PRECISION MATRIX
        ! array(1) x array(2) multiplier
                implicit none
                double precision, dimension(:), intent(IN) :: A
                double precision, dimension(:,:), intent(IN) :: B
                double precision, dimension(size(A)) :: C
                C=matmul(A,B)
        END FUNCTION VECMATMUL
END MODULE fastmul
