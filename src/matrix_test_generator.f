!File=src/matrix_test_generator.f
!Author=fabio
!Created=mar 08 mag 2018 14:41:03 CEST
!Last Modified=mar 08 mag 2018 14:41:03 CEST
program matrix_test_generator
    integer :: n=500
    double precision :: A
    function hilbert(n) result(B)
        ! build the hilbert matrix of size nxn
        implicit none
        integer, intent(IN) :: n
        double precision, dimension(n,n) :: B
        ! local variables
        integer :: i,j
        ! body of the function
        do i=1,n
        do j=1,n
            A(i,j)=1/(i+j-1)
        end do
        end do

        
    end function

end program matrix_test_generator

