MODULE LRA
        implicit none
        private
        public :: DIAG         ! diagonal matrix from vector
        public :: SVD          ! Svd decomposition
        public :: QR           ! Ortonormalize in NON-compress form a bases
        public :: ON           ! ortonormalize in compress form a bases
        public :: GETQ         ! decompress a orto-normal base
        public :: QMULL        ! multuply on left a unitary houseolder matrix
        public :: QMULR        ! multiply on right a unitary householder matrix
contains
        PURE FUNCTION DIAG(S) RESULT(D)
                ! crrete diagonal matrix from vector S
                double precision, dimension(:), intent(IN) :: S
                double precision, dimension(size(S),size(S)) :: D
        ! local variables
                integer :: I, N
                D=0
                N=size(S)
                do I=1,N
                        D(I,I)=S(I)
                end do
        END FUNCTION DIAG
        SUBROUTINE SVD(A,U,B,V,K,TOL)
        ! K, integer               - bound of numerical rank
        ! TOL, double precision    - epsilon rank accuracy
        ! Z, double precision(*,*) - matrix to apply decomposition
        ! U, double precision(*,*) - left singolar vector
        ! V, double precision(*,*) - right singolar vector
        ! Apply svd decomposition to a real martix Z and check if
        ! numerical rank with tolerance TOL is bounded from K
                implicit none
        !variabili in input        
                integer, intent(INOUT) :: K
                double precision, intent(IN) :: TOL
                double precision, dimension(:,:), intent(IN) :: A
                double precision, allocatable, dimension(:,:),intent(OUT) :: U,V
                double precision, allocatable, dimension(:), intent(OUT) :: B

        !local scalars
                character :: JOBU, JOBVT
                logical :: contr=.true.
                integer :: LDA, LDU, LDVT, LWORK, INFO, DMIN, I, M, N, KB
        !local array
                double precision, dimension(size(A,1),size(A,2)) :: Z
                double precision, dimension(min(size(A,1),size(A,2)),size(A,2)) :: VT
                double precision, dimension(size(A,2),min(size(A,1),size(A,2))) :: VAUX
                double precision, dimension(size(A,1),min(size(A,1),size(A,2))) :: UAUX
                double precision, dimension(min(size(A,1),size(A,2))) :: BAUX
                double precision, allocatable, dimension(:) :: WORK
                !inizializzazione parametri
                M=size(A,1)
                N=size(A,2)
                DMIN=min(M,N)
                JOBU='S'
                JOBVT='S'
                LDA=M
                LDU=M
                LDVT=DMIN
                LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
                !controlli
                if (tol.le.0) then
                        write(*,*) 'TOL deve essere positiva'
                stop -2
                end if
                ALLOCATE(WORK(LWORK))
                Z=A
                CALL dgesvd(JOBU,JOBVT,M,N,Z,LDA,BAUX,UAUX,LDU,VT,LDVT,WORK,LWORK,INFO)
                ! calcolo rango numerico rispetto a tol
                KB=0
                do I=0,DMIN-1
                        if (BAUX(DMIN-I).ge. TOL) then
                                KB=DMIN -I
                                contr=.false. 
                                exit
                        end if
                end do
                if (contr) then ! mi accerto che tol piccolo abbastanza
                        KB=0
                        write(*,*) 'WARNING: 0 numerical rank with tol: ', TOL
                end if
                if (K < KB) then
                        write(*,*) 'WARNING:numerical rank:', KB, ', greater than bound K: ',K, ', with TOL ', TOL
                end if
                K=min(K,KB)
                U=UAUX(:,1:K)
                B=BAUX(1:K)
                VAUX=transpose(VT)
                V=VAUX(:,1:K)
                DEALLOCATE(WORK) !libero memoria
        END SUBROUTINE SVD
        FUNCTION QR(A) RESULT(U)
        ! A, double precision, dimension(:,:) - Bases (lunga e stretta!!)
        ! U, double precision, dimension(size(A)) - Ortonormal basis
        ! La function prende in input una base (non usare per rango basso)
        !  rappresentata dalle
        ! colonne di A, restituisce la base ortornomale ottenuta
        ! con il metodo econ-QR
        ! N.B non viene memorizzata la matrice R!!!
                implicit none
                double precision, dimension(:,:), intent(IN) :: A
        ! local variables
                double precision, dimension(size(A,1),size(A,2)) :: U
                double precision, dimension(size(A,1),size(A,1)) :: Q
                double precision, dimension(size(A,2)) :: TAU
                double precision, dimension(size(A,1),1) :: V
                double precision, dimension(size(A,2)) :: WORK
                integer :: M, N, LDA, INFO, I, J, LWORK
                M=size(A,1)
                N=size(A,2)
                if(M<N) then
                        write(*,*) 'A is not a bases'
                        stop -1
                end if
                LDA=M
                LWORK=N
                U=A
                call dgeqrf(M,N,U,LDA,TAU,WORK,LWORK,INFO)
!                write(*,'(13F6.2)') ((U(I,J),J=1,13),I=1,M)
                Q=0.d0
                do I=1,M
                        Q(I,I)=1.d0
                end do
                do I=1,N
                        V(1:I-1,1)=0.d0
                        V(I,1)=1.d0
                        V(I+1:M,1)=U(I+1:M,I)
                        Q = Q - matmul(matmul(Q,TAU(I)*V),transpose(V))
                end do
                U=Q(:,1:N)
        END FUNCTION QR
        FUNCTION ON(A) RESULT(U)
        ! A, double precision, dimension(:,:) - Bases
        ! U, double precision, dimension(size(A)) - ON BASES
        ! La function prende in input una base (non usare per rango basso)
        !  rappresentata dalle
        ! colonne di A, restituisce la base ortornomale ottenuta
        ! con il metodo econ-QR organizzata dentro U nel modo seguente
        ! U(J<I) ha per colonne i riflettor
        ! diag(U) i coeffcenti che servono per costruire Q
        ! N.B non viene memorizzata la matrice R!!!
! Variabili di input 
        double precision, intent(IN), dimension(:,:) :: A
        double precision, dimension(size(A,1),size(A,2)) :: U
!variabili locali
        integer :: M, N, LWORK, DMIN, LDA, INFO, I
        double precision, allocatable, dimension(:) :: WORK
        double precision,dimension(min(size(A,1),size(A,2))) :: TAU
!controllo
!dgeqrf( M, N, A, LDA, TAU, WORK, LWORK, INFO )
        M=size(A,1)
        N=size(A,2)
        DMIN=min(M,N)
        LDA=M
        LWORK=N
        allocate(WORK(LWORK))
        U=A
        call dgeqrf( M, N, U, LDA, TAU, WORK, LWORK, INFO)
        deallocate(WORK)
        do I=1,DMIN
                U(I,I)=TAU(I)
        end do
        END FUNCTION ON
        FUNCTION GETQ(U) RESULT(QOUT)
        ! decomprime la matrice unitaria U ottenuta con ON
                implicit none
                double precision, intent(IN), dimension(:,:) :: U
                double precision, dimension(size(U,1),size(U,2)) :: QOUT
        ! local variables
                integer :: M, N, I, DMIN
                double precision, dimension(size(U,1),1) :: V
                double precision, dimension(size(U,1),size(U,1)) :: Q
                M=size(U,1)
                N=size(U,2)
                DMIN=min(M,N)
                Q=0
                do I=1,M
                        Q(I,I)=1
                end do        
                do I=1,DMIN
                        V(1:I-1,1)=0
                        V(I,1)=1
                        V(I+1:M,1)=U(I+1:M,I)
                        Q= Q- matmul(matmul(U(I,I)*Q,V),transpose(V))
                end do
                QOUT=Q(:,1:N)
        END FUNCTION GETQ
        FUNCTION QMULL(U,A) result(B)
        ! A, double precision, dimension(:,:)
        ! U, ON(OUT)
        ! moltiplica a sinistra la matrice unitaria U 
        ! ottenuta con ON per una matrice A
                implicit none
                double precision, dimension(:,:), intent(IN) :: A, U
                double precision, dimension(size(U,1),size(A,2)) :: B
        ! local variables
                integer :: M, N, K, DMIN, I
                double precision, dimension(size(U,1),1) :: V
                M=size(U,1)
                N=size(U,2)
                K=size(A,2)
                DMIN=min(M,N)
                B=0
                B(1:K,:)=A(1:K,:)
                do I=DMIN,1,-1
                        V(1:I-1,1)=0
                        V(I,1)=1
                        V(I+1:M,1)=U(I+1:M,I)
                        B= B - matmul(U(I,I)*V,matmul(transpose(V),B))
                end do
        END FUNCTION QMULL
        FUNCTION QMULR(A,U) result(B)
        ! A, double precision, dimension(:,:)
        ! U, ON(OUT)
        ! moltiplica a destra la matrice unitaria U 
        ! ottenuta con ON per una matrice A
                implicit none
                double precision, dimension(:,:), intent(IN) :: A, U
                double precision, dimension(size(A,1),size(U,2)) :: B
        ! local variables
                integer :: M, N, DMIN, I
                double precision, dimension(size(U,1),1) :: V
                double precision, dimension(size(A,1),size(A,2)) :: BAUX
                M=size(U,1)
                N=size(U,2)
                DMIN=min(M,N)
                BAUX=A
                do I=1,DMIN
                        V(1:I-1,1)=0
                        V(I,1)=1
                        V(I+1:M,1)=U(I+1:M,I)
                        BAUX= BAUX - matmul(matmul(U(I,I)*BAUX,V),transpose(V))
                end do
                B=0
                B(:,1:N)=BAUX(:,1:N)
        END FUNCTION QMULR
END MODULE    
