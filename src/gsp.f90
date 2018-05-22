PROGRAM testgmp
        implicit none
!variabili 
        integer :: M=10, N=10, DMIN, I, J
        double precision :: TOL
        double precision, allocatable, dimension(:,:) :: A, R, AOLD
!inizailizzazione
        DMIN=min(M,N)
        TOL=1.d-20
        allocate(A(M,N), AOLD(M,N), R(DMIN,N))
        call random_number(A)
        AOLD=A
        call gsp(M,N,DMIN,A,R,TOL)
        write(*,*) norm2(A(:,DMIN))
        write(*,*) norm2(AOLD(:,DMIN) - matmul(A(:,DMIN),R)) 
        WRITE(*,'(10F8.2)') ((R(I,J),J=1,N),I=1,DMIN)
END PROGRAM
SUBROUTINE gsp(M,N,DMIN,A,R,TOL)
        implicit none
!variabili di input
        integer :: M, N, DMIN
        double precision ::  TOL
        double precision :: A(M,N), R(DMIN,N)

!variabili locali
        integer :: I, J
        double precision :: VTEMP(M)
        R=0
        VTEMP=A(:,1) ! primo vettore della base ortonormale
        R(1,1)=norm2(VTEMP)
        A(:,1)=VTEMP/R(1,1)
        if (abs(R(1,1)) .ge. TOL) then
                A(:,1)=A(:,1)/R(1,1)
        else
                stop -1
        end if
        do J=2,N 
                VTEMP=A(:,J)
                do I=1,J-1
                        R(I,J)=dot_product(reshape(A(:,I),(/M/)),VTEMP)
                        A(:,J)=A(:,J) - R(I,J)*A(:,I)
                end do
                R(J,J)=norm2(A(:,J))
                if (abs(R(J,J)) .ge. TOL) then
                        A(:,J)=A(:,J)/R(J,J)
                else
                        stop -J
                end if
        end do
END SUBROUTINE
