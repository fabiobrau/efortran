MODULE HMOD
        ! In this module there are fundamental functions to work with matrix
        ! in Hierachical Off Diagonal Low Rank form, the strucure is recursive
        use LRA
        implicit none
        private
        type :: hm
                type(hm), pointer :: A11, A22
                double precision, allocatable, dimension(:,:) :: U12, U21, V12, V21
                double precision, allocatable, dimension(:) :: B12, B21
                double precision, allocatable, dimension(:,:) :: F
        end type hm
        interface hm
                module procedure hm_construct
        end interface hm        
        !! PUBLIC TYPE/CONSTRUCTOR
        public :: hm
        !! PUBLIC FUNCTION/ SUBROUTINE
        public :: random_sampling
        public :: randn     ! random array normal distribuition
contains
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
        SUBROUTINE RANDOM_SAMPLING(A,U,B,V,K)
                ! comprime con la tecnica di random sampling una matrice
                implicit none
                double precision, dimension(:,:),intent(IN) :: A
                integer, intent(INOUT) :: K
                double precision, allocatable, dimension(:,:), intent(OUT) :: U, V
                double precision, allocatable, dimension(:), intent(OUT) :: B
        ! local variables
                double precision, allocatable, dimension(:,:):: O
                double precision, allocatable, dimension(:,:) :: UAUX
                double precision :: TOL
                integer :: R, M, N, P
                P=5
                N=size(A,2)
                R=min(P+K,N)
                TOL=1e-15
                O=randn(N,R)
                U=qr(matmul(A,O))
                O=matmul(transpose(A),U)
                call SVD(O,V,B,UAUX,K,TOL)
                U=matmul(U,UAUX)
        END SUBROUTINE
        RECURSIVE FUNCTION HM_CONSTRUCT(A,K) RESULT(HA)
        ! costruttor of type HODLR,
        ! A, double precision, matrix to compress
        ! K, integer, bound off-diagonal rank
                implicit none
                type(hm) :: HA
                double precision, dimension(:,:), intent(IN) :: A
                integer, intent(IN) :: K
        ! local variables
                double precision :: TOL
                double precision, allocatable, dimension(:,:) :: U, V
                double precision, allocatable :: B(:)
                integer :: N, MID, R
              
                N=size(A,1)
                MID =ceiling(real(N/2))
                if (K>MID) then
                        HA%F=A
                else
                        R=K
                        call random_sampling(A(MID+1:N,1:MID), U,B,V,R)
                        allocate(HA%U21(N,size(U,2)),HA%B21(size(B)),HA%V21(N,size(V,1)))
                        HA%U21=U
                        HA%V21=V
                        HA%B21=B
                        R=K
                        call random_sampling(A(1:MID,MID+1:N),U,B,V,R)
                        allocate(HA%U12(N,size(U,2)),HA%B12(size(B)),HA%V12(N,size(V,1)))
                        HA%U12=U
                        HA%V12=V
                        HA%B12=B
                        allocate(HA%A11,HA%A22)
                        HA%A11 = hm_construct(A(1:MID,1:MID),K)
                        HA%A22 = hm_construct(A(MID+1:N,MID+1:N),K)
                end if
        END FUNCTION HM_CONSTRUCT


END MODULE HMOD
