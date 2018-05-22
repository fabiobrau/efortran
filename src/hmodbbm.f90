MODULE HMODBBM
    ! In this module there are fundamental functions to work with matrix
    ! in Hierachical Off Diagonal Low Rank form.
    ! We assume structured matrix the black box multiplier can be overload
    ! to perform operation.
    ! All strcuture are explicit, no recursion there.
    use LRA
    implicit none
    private
    type :: HODMAT
        double precision, allocatable, dimension(:,:) :: U, V, B
        ! U must be N x R
        ! V must be N x R
        ! B at level L must be 2**L x R
    end type HODMAT
    type :: HMAT
        integer, allocatable, dimension(:,:) :: T
        type(HODMAT), allocatable, dimension(:) :: OD
        double precision, allocatable, dimension(:,:) :: D
        ! T tree, dimension(2**(LL+1),2), with LL #levels
        ! OD Off Diagonal rappresentation with rank K, dimension(LL)
        ! D diagonal blocks, dimension(N ,NMAX) NMAX maximum length of leaves.
    end type HMAT
    interface transpose
        module procedure htrans
        module procedure hodtrans
    end interface transpose
    interface matmul
        module procedure hmull
    end interface matmul
    interface size
        module procedure sizehmat
    end interface size
    interface eye
        module procedure eyefull
        module procedure eyecomp
    end interface eye
    !! PUBLIC OVERLOAD
    public :: matmul
    public :: transpose
    public :: size
    !! PUBLIC TYPE
    public :: hmat      ! hierarchical rapresentation of a matrix HODLR
    public :: hodmat    ! hierach. off diag. rappresentation of a matrix HODLR
    !! PUBLIC FUNCTION/ SUBROUTINE
    public :: randnlow  ! create random matrix with low rank off diagonl blocks
    public :: eye       ! create identity in compress form
    public :: hodtrans  ! traspose compress truncated form (NO DIAGONAL BLOCKS)
    public :: hodmull   ! L multiply compress truncated form (NO DIAGONAL BLOCKS)
    public :: hmull     ! Lmultiply compress completed form (WITH DIAGONAL BLOCKS)        
    public :: tree      ! get a binary tree
    public :: randn     ! random array normal distribuition
    public :: sub       ! hierarchically sub-matrix
    public :: xdim      ! node dimension as intervall
    public :: ltronc    ! testing function
    public :: hodlrrs   ! hierarchically off diagonal low rank random sampling
    public :: hlrrs     ! hiechically completed form (WITH DIAGONAL BLOCKS)
    public :: level     ! get the level
contains
    FUNCTION SIZEHMAT(HA) RESULT(S)
        ! calculate storage memory cost
        implicit none
        type(HMAT) :: HA
        integer :: S
    ! local variables
        integer :: I,LL
        S=size(HA%D)
        LL=size(HA%OD)
        do I=1,LL
            S=S+size(HA%OD(I)%U)+size(HA%OD(I)%V)+size(HA%OD(I)%B)
        end do
    END FUNCTION SIZEHMAT
        FUNCTION RANDNLOW(K,LL,N0) RESULT(A)
                ! build a random matrix such as the off diagonal
                ! have rank lower than K
                implicit none
                integer, intent(IN) :: K, LL, N0
                double precision, dimension(2**LL*N0,2**LL*N0) :: A
        ! local variables
                double precision, dimension(N0*2**LL,K) :: U
                integer, dimension(2**(LL+1)-2,2) :: T
                integer :: N, I, J
                if (k>N0) then
                        write(*,*) 'K must be lower than N0'
                        stop -1
                end if
                T=tree(LL,N0)
                N=2**LL*N0
                U=randn(N,K)
                A=matmul(U,transpose(U))
                do J=1,2**LL
                        I=2**LL +J -2
                        A(T(I,1):T(I,2),T(I,1):T(I,2))=randn(N0,N0)
                end do
        END FUNCTION RANDNLOW  
        FUNCTION eyefull(N) RESULT(II)
        ! create full identity matrix
                implicit none
                integer :: N
                double precision, dimension(N,N) :: II
        ! local variables
                integer :: I
                II=0
                do I=1,N
                        II(I,I)=1.d0
                end do
        END FUNCTION eyefull
        FUNCTION eyecomp(LL,N0,R) RESULT(II)
        ! create compress HODLR form of identity matrix
        ! dimension it is N=2*LL*N0
        ! R is the numerical rank bound
                implicit none
                integer, intent(IN) :: LL, N0,R
                type(HMAT) :: II
        ! local variables
                type(HODMAT), dimension(LL) :: OD
                double precision, dimension(N0*2**LL,R) :: U, D
                double precision, dimension(2**LL,R) :: B
                double precision, dimension(R,R) :: ID
                integer, dimension(2**(LL+1)-2,2) :: T
                integer :: I, J, L, ENDD, h,k                
                T=tree(LL,N0)
                ID=eyefull(R)
                do L=1,LL
                        U=0.d0
                        do J=1,2**L
                                I=2**L-2+J
                                U(T(I,1):T(I,1)+R-1,:)=ID
                        end do
                        B=0.d0
                        OD(L)=hodmat(U,U,B(1:2**L,:))
                end do
                II=hmat(T,OD,U)

        END FUNCTION eyecomp
        FUNCTION tree(LL,N0) result(T)
        ! LL, scalar integer - tree depht (#levels)
        ! N0, scalar integer - leaf dimension
        ! T, integer, dimension(2**(LL+1)-2,2)- tree arranged in matrix form
        ! create a bilanced tree, all leaves are at the same level LL
        ! The tree T it's memorized in this way:
        ! The Ith row contains the extremal of the Ith node as interval.
                implicit none
                integer, intent(IN) :: N0, LL
                integer, dimension(2**(LL+1)-2,2) :: T
        ! local variables
                integer :: N, L, I, ML
                N=2**LL*N0
                do L=1,LL
                        ML=N/2**L
                        do I=1,2**L
                                T(2**L -2 + I,1)=1+(I-1)*ML
                                T(2**L -2 + I,2)=I*ML
                        end do
                end do
        END FUNCTION tree
        FUNCTION RANDN(N,K) RESULT(X)
        ! funzione di interfaccia per dlarnv del lapack,
        ! crea matrice NxK con distribuzione normale (0,1)
                implicit none
                integer, intent(IN) :: N,K
                double precision, dimension(N,K) :: X
        ! local variables
                integer, dimension(4) :: ISEED
                double precision, dimension(4) :: DSEED
                double precision, dimension(N*K) :: XVEC
                call random_number(DSEED)
                ISEED=floor(DSEED*4095)
                ISEED(4)=2*floor(DSEED(4)*2047) +1
                call dlarnv(3,ISEED,N*K,XVEC)
                X=reshape(XVEC,(/N,K/))
        END FUNCTION RANDN
        PURE FUNCTION XDIM(T,X) RESULT(NX)
        ! T, integer, dimension(:,2) - albero anche non bilanciato
        ! X, integer scalar          - posizione del nodo
        ! restituisce la lunghezza del intervallo corrispondente al nodo X
                implicit none
                integer, intent(IN) :: X
                integer, dimension(:,:), intent(IN) :: T
                integer :: NX
                NX = T(X,2) - T(X,1) +1
        END FUNCTION XDIM
        PURE FUNCTION SUB(A,X,Y,T) RESULT(B)
        ! A, real dimension(:,:) - matrice principale
        ! X,Y, integer scalar    - posizioni dei nodi sull'albero
        ! per convenzione X=0 significa A(:,Y)
        ! T, integer dimension(:,2)- binary tree (even not equidist.)
        ! Restituisce il sottoblocco della matrice A
        ! corrispondente ai nodi X Y nell'albero T
        ! ATTENZIONE: NON USARE A SINISTRA DI UN ASSEGNAMENTO
                implicit none
                integer, intent(IN) :: X ,Y
                integer, dimension(:,:), intent(IN) :: T
                double precision, dimension(:,:), intent(IN) :: A
                double precision, allocatable, dimension(:,:) :: B
        ! local variables
                integer :: NX ,NY
                if (X==0 .and. Y==0) then
                        NX=size(A,1)
                        NY=size(B,1)
                        allocate(B(NX,NY))
                        B=A
                else if (X==0) then
                        NX=size(A,1)
                        NY=xdim(T,Y)
                        allocate(B(NX,NY))
                        B=A(:,T(Y,1):T(Y,2))
                else if (Y==0) then
                        NX=xdim(T,X)
                        NY=size(A,2)
                        allocate(B(NX,NY))
                        B=A(T(X,1):T(X,2),:)
                else
                        NX=xdim(T,X)
                        NY=xdim(T,Y)
                        allocate(B(NX,NY))
                        B=A(T(X,1):T(X,2),T(Y,1):T(Y,2))
                end if
        END FUNCTION SUB
        FUNCTION LTRONC(A,T,L) result(B)
        ! It is a testing function.
        ! Costruisce la matrice troncata al livello L dell'albero T
                implicit none
                integer, intent(IN) :: L
                integer, dimension(:,:), intent(IN) :: T
                double precision, dimension(:,:), intent(IN) :: A
                double precision, dimension(size(A,1),size(A,2)) :: B
        ! local variables
                integer :: I, XLAST
                XLAST=2**(L+1)-2 ! considero ultimo nodo del livello L 
                B=0
                do I=1,XLAST-1,2
                        B(T(I,1):T(I,2),T(I+1,1):T(I+1,2))=sub(A,I,I+1,T)
                        B(T(I+1,1):T(I+1,2),T(I,1):T(I,2))=sub(A,I+1,I,T)
                end do
        END FUNCTION LTRONC
        SUBROUTINE HODLRRS(A,AZIP,K,P,LL,N0,T)
        ! LL, scalar, integer           - number of level
        ! N0, scalar, integer           - N/2**LL
        ! P, scalar, integer,           - oversamplig parameter
        ! K, scalar, integer            - rank off daigonal blocks bound
        ! T, integer, dimension(*,*)            - Tree (OPTIONAL)
        ! A, double precision, dimension(*,*)   - Matrix to compress
        ! AZIP, type(HODLR), dimension(LL)      - Compress matrix (OUT)
            
        ! Compress the matrix off diaganl blocks of A 
        ! ,using Random sampling, into HODLR matrix form
        ! out/in/opt variables
                implicit none
                integer, intent(IN) :: K, P, LL, N0
                double precision, intent(IN), dimension(:,:) :: A
                integer, dimension(:,:), intent(IN) :: T
                type(HODMAT), dimension(LL) , intent(OUT):: AZIP
       ! local variables
                double precision, dimension(N0*2**LL,K+P) :: O1,O2 ,Y1 ,Y2, Z1, Z2, U, V
                double precision, dimension(2**LL ,K+P) :: B
                type(HODMAT), dimension(LL) :: BZIP
                double precision :: TOL=1d-16
                double precision, dimension(K+P,K+P) :: UAUX
                integer :: I, J, L, N , R, NA, NB
                R=K+P
                N=size(A,1)
                do L=1,LL ! loop sui livelli 
                        ! costuisco le matrici casuali
                        ! Invece che richiamare tante volte randn
                        ! pongo a zero i blocchi pari per O1 e dipsari per O2
                        O1=0.d0
                        O2=0.d0
                        write(*,*) 'Processing level:',L
                        do J=1,2**L-1,2 !II,II+1 sibling iif I odd
                                I=2**L -2 +J
                                NA=xdim(T,I)
                                NB=xdim(T,I+1)
                                O1(T(I,1):T(I,2),:)=randn(NA,R)
                                O2(T(I+1,1):T(I+1,2),:)=randn(NB,R)
                        end do
                        ! costuisco le matrici campione
                        Y1=matmul(A,O2)! - matmul(ltronc(A,T,L),O2)
                        Y2=matmul(A,O1)! - matmul(ltronc(A,T,L),O1)
                        if (L>1) then ! uso compr. al livello prec
                        ! moltiplico A^(l) O1 e O2 in modo ottimizzato usando il livello prec
                                Y1=Y1 - hodmull(AZIP,O2,L-1,T)
                                Y2=Y2 - hodmull(AZIP,O1,L-1,T)
                        end if
                        ! aggiorno O1 O2 con la base ortonormale in forma esplicita
                        O1=0.d0
                        O2=0.d0
                        do J=1,2**L-1,2
                                I=2**L -2 +J
                                U(T(I,1):T(I,2),:)=qr(sub(Y1,I,0,T))
                                U(T(I+1,1):T(I+1,2),:)=qr(sub(Y2,I+1,0,T))
                                O1(T(I,1):T(I,2),:)=sub(U,I,0,T)
                                O2(T(I+1,1):T(I+1,2),:)=sub(U,I+1,0,T)
                        end do
                        Z1=matmul(transpose(A),O2)
                        Z2=matmul(transpose(A),O1)
                        if (L>1) then ! Completo costruzione di Z1 e Z2
                                BZIP(1:L-1)=transpose(AZIP(1:L-1)) ! transposing the compress type
                                Z1=Z1 - hodmull(BZIP,O2,L-1,T)
                                Z2=Z2 - hodmull(BZIP,O1,L-1,T)
                        end if
                        V=0.d0
                        ! calcolo la scomposizione SVD, aggiorno la forma compressa di A
                        do J=1,2**L-1,2
                                I=2**L -2 +J
                                write (*,*) 'processing node', J
                                ! blocco I+1,I
                                call svd(sub(Z1,I,0,T),V(T(I,1):T(I,2),:),UAUX,B(J+1,:),K,TOL)
                                U(T(I+1,1):T(I+1,2),:)=matmul(U(T(I+1,1):T(I+1,2),:),UAUX)
                                ! fratello destro
                                call svd(sub(Z2,I+1,0,T),V(T(I+1,1):T(I+1,2),:),UAUX,B(J,:),K,TOL)
                                U(T(I,1):T(I,2),:)=matmul(U(T(I,1):T(I,2),:),UAUX)
                        end do
                        !write(*,'(13F6.2)') B
                        AZIP(L)=hodmat(U,V,B(1:2**L,:))
                end do
        END SUBROUTINE HODLRRS
        SUBROUTINE  HLRRS(A,HA,P,K,LL,N0,T)
        ! Take in input a square matrix A and compress with random sampling
                implicit none
        ! input variables
                integer, intent(IN) :: LL, N0, K,P
                integer, dimension(:,:), intent(IN) :: T
                double precision, dimension(:,:), intent(IN) :: A
                type(HMAT), intent(INOUT) :: HA
        ! local variables
                integer :: NMAX,I, J, N, ENDD
                type(HODMAT), dimension(LL) :: ODAUX
                double precision, allocatable, dimension(:,:) :: D
                N=N0*2**LL
                ! massima dimensione delle foglie
                NMAX=0
                do J=1,2**LL
                        I=2**LL -2 +J
                        if (NMAX<xdim(T,I)) then
                                NMAX=xdim(T,I)
                        end if
                end do
                call hodlrrs(A,ODAUX,K,P,LL,N0,T)
                allocate(D(N,NMAX))
                do J=1,2**LL
                        I=2**LL-2+J
                        ENDD=xdim(T,I)
                        D(T(I,1):T(I,2),1:ENDD)=eye(ENDD)
                end do
                D=matmul(A,D) - hodmull(ODAUX,D,LL,T)
                do J=1,2**LL
                        I=2**LL-2+J
                        ENDD=xdim(T,I)
                        D(T(I,1):T(I,2),ENDD+1:NMAX)=0.d0
                end do
                HA=hmat(T,ODAUX,D)
        END SUBROUTINE HLRRS
        FUNCTION HODMULL(AZIP,X,LL,T) RESULT(Y)
        ! Esegue il prodotto fra AZIP*X fino al livello LL in AZIP
        ! La matrice AZIP compressa tramite svd in modo che
        ! AZIP(l)%U, dimension(N,K) vettori singolari sinistri blocchi al livello l
        ! AZIP(l)%V, dimension(N,K) vettori singolari destri blocchi livello l
        ! AZIP(l)%B(I,:), dimension(2**l,K) valori singolari blocco A_I,I+1 I dispari
        ! AZIP(l)%B(I,:), dimension(2**l,K) valori singolari blocco A_I,I-1 I pari
        ! T, integer, dimension(:,2) albero binario
                implicit none
                type(HODMAT), dimension(:), intent(IN) :: AZIP
                double precision, dimension(:,:), intent(IN) :: X
                integer, dimension(:,:), intent(IN) :: T
                integer, intent(IN) :: LL
                double precision, dimension(size(X,1),size(X,2)) :: Y
        ! local variables
                double precision, allocatable, dimension(:,:) :: SIG
                double precision, dimension(size(AZIP(1)%U,1),size(AZIP(1)%U,2)) :: U, V
                integer :: L, I, J, N, R, II
                N=size(AZIP(1)%U,1)
                R=size(AZIP(1)%U,2)
                allocate(SIG(R,R))
                Y=0.d0
                do L=1,LL
                        U=AZIP(L)%U
                        V=AZIP(L)%V
                        do J=1,2**L-1,2
                                I=2**L-2+J
                                SIG=0.d0 
                                do II=1,R
                                        SIG(II,II)=AZIP(L)%B(J,II)
                                end do
                                Y(T(I,1):T(I,2),:)=Y(T(I,1):T(I,2),:)&
                                +matmul(sub(U,I,0,T),matmul(SIG,matmul(&
                                transpose(sub(V,I+1,0,T)),X(T(I+1,1):T(I+1,2),:))))
                                SIG=0.d0
                                do II=1,R
                                        SIG(II,II)=AZIP(L)%B(J+1,II)
                                end do
                                Y(T(I+1,1):T(I+1,2),:)=Y(T(I+1,1):T(I+1,2),:)&
                                +matmul(sub(U,I+1,0,T),matmul(SIG,matmul(&
                                transpose(sub(V,I,0,T)),X(T(I,1):T(I,2),:))))
                        end do
                end do
        END FUNCTION HODMULL
        FUNCTION HMULL(HA,X) RESULT(Y)
        ! Esegue il prodotto fra la forma compressa tramite HODLRRS e un array double 
                implicit none
                type(HMAT), intent(IN) :: HA
                double precision, dimension(:,:), intent(IN) :: X
                double precision, dimension(size(X,1),size(X,2)) :: Y
        ! local variables
                integer :: I, J, ENDD, LL, LTREE
                LTREE=size(HA%T,1)
                LL=level(LTREE)               
                Y=HODMULL(HA%OD,X,LL,HA%T)
                !write(*,'(20F6.2)') ((Y(I,J),J=1,20),I=1,20)
                do J=1,2**LL
                        I=2**LL-2+J
                        ENDD=xdim(HA%T,I) ! Dimensione effettiva del blocco diagonale
                        Y(HA%T(I,1):HA%T(I,2),:)=sub(Y,I,0,HA%T)+matmul(&
                                HA%D(HA%T(I,1):HA%T(I,2),1:ENDD),sub(X,I,0,HA%T))
                end do
        END FUNCTION HMULL
        FUNCTION HODTRANS(AZIP) RESULT(BZIP)
        ! Transpose a HODLR matrix (without diagonal)
                type(HODMAT), dimension(:), intent(IN) :: AZIP
        ! local variables
                type(HODMAT), dimension(size(AZIP,1)) :: BZIP
                double precision, dimension(2**size(AZIP,1),size(AZIP(1)%B,2)) :: B
                integer :: I, LL, L, N, R
                N=size(AZIP(1)%U,1)
                R=size(AZIP(1)%U,2)
                LL=size(AZIP,1)
                do L=1,LL
                        B=0.d0
                        B((/((I),I=1,2**L-1,2)/),:)=AZIP(L)%B(2:2:2**L,:)
                        B((/((I),I=2,2**L,2)/),:)=AZIP(L)%B(1:2:2**L-1,:)
                        BZIP(L)=hodmat(AZIP(L)%V,AZIP(L)%U,B(1:2**L,:))
                end do
        END FUNCTION HODTRANS
        FUNCTION HTRANS(HA) result(HB)
        ! Transpose a compressed type matrix
                implicit none
                type(HMAT), intent(IN) :: HA
                type(HMAT) :: HB
        ! local variables
                integer :: I,L,LL, ENDD
                double precision, dimension(size(HA%D,1),size(HA%D,2)) :: D
                LL=size(HA%OD,1)
                do I=2**LL-1,2**(LL+1)-2
                        ENDD=xdim(HA%T,I) ! can have not same dimension
                        D(HA%T(I,1):HA%T(I,2),1:ENDD)=transpose(HA%D(HA%T(I,1):HA%T(I,2),1:ENDD))
                end do
                HB=hmat(HA%T,transpose(HA%OD),D)
        END FUNCTION
        PURE FUNCTION LEVEL(X) result(L)
        ! Data una posizione di un nodo in un albero binario, restituisce
        ! il livello del nodo corrispondente alla posizione x
                integer, intent(IN) :: X
                integer :: L
        ! local variables
                do L=1,32 ! massimo numero naturale 2**32
                        if ((X -(2*2**L - 2)).le. 0) then
                                exit
                        end if
                end do
        END FUNCTION LEVEL
END MODULE HMODBBM
