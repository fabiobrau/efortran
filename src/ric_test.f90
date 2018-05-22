PROGRAM TEST
        use HMOD
        use LRA
        implicit none
        type(hm) :: HA
        double precision :: TOL
        double precision, allocatable, dimension(:,:) :: A, U, V
        double precision, allocatable, dimension(:) :: B
        integer :: K,M
        TOL=1e-15
        M=1000
        K=7
        A=matmul(randn(M,10),randn(10,M))
        HA=hm(A,10)

END PROGRAM
