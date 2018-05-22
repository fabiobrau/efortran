MODULE HMOD
        implicit none
contains
SUBROUTINE HODLR(A,N,K,LL,AZIP,U,V) ! provisorio, tutto con svd
! N, scalar                 - N= 2**L * N0
! A, dimension NxN          - matrice da comprimere
! K, scalar                 - rank off-diagonal submatrix
! L, scalar                 - numero di livelli
! AZIP, dimension (K,K,N/N0)- approssimazioni dei blocchi, indicizzate per posizione
! U,V,  dimension (k*N*L,K) - basi ortonormali per scompattare

! Prende in input una matrice A, un limite al rango k e organizza in forma gerarchica
! in un tensore AZIP(1,j,k) con k posizione del nodo.
! le matrici U e V servono per decomprimere
        implicit none
! variabili di input/output
        integer, intent(IN) :: K, N, LL
        double precision, intent(IN) :: A(N,N)
        double precision, intent(OUT) :: AZIP(K,K,2**(LL+1)-2), U(N*LL,K), V(N*LL,K)
! variabili locali
        integer :: I, J, DMIN, ML, L, FROM, TO, POSA, POSUV
        double precision :: TOL 
        double precision, allocatable, dimension(:) :: S
        double precision, allocatable, dimension(:,:) :: AAUX, UAUX, VAUX
        
        TOL =1.e-10
        do L=1,LL ! L sta per livello
        write(*,*) 'Livello -> ', L 
                ML=N/2**L ! indica la dimensione di ogni matrice del livello corrente
                allocate(S(ML),AAUX(ML,ML),UAUX(ML,ML),VAUX(ML,ML))
                do J=1,2**L -1,2 ! J sta per fratello sinistro
                        write(*,*) 'Figlio J ->', J
                ! semplificano la notazione, descrivono posizione del blocco corrente
                        FROM = 1+(J-1)*ML 
                        TO = J*ML
                ! semplifica la notazione, indica il numero del blocco corrente
                        POSA = 2**L+(J-2)
                        POSUV= (L-1)*N + (J-1)*ML
                ! archivio il sottoblocco sopra diagonale
                        AAUX=A(FROM:TO,FROM+ML:TO+ML)
                        DMIN=ML ! serve per la funzione svd che sovrascrive ML
                        call svd(AAUX,UAUX,VAUX,S,ML,ML,DMIN,TOL)
                        do I=1,ML
                                AZIP(I,I,POSA) = S(I)
                        end do 
                        U(POSUV:POSUV+ML,:)=UAUX(:,1:K)
                        V(POSUV:POSUV+ML,:)=VAUX(:,1:K)
                ! ripeto la procedura per il blocco sottodiagonale
                ! semplifica la notazione, indica il numero del blocco corrente
                        POSA = 2**L+(J-2) + 1
                        POSUV= (L-1)*N + (J-1)*ML + ML
                ! archivio il sottoblocco sottodiagonale
                        AAUX=A(FROM+ML:TO+ML,FROM:TO)
                        DMIN=ML ! serve per la funzione svd che sovrascrive ML
                        call svd(AAUX,UAUX,VAUX,S,ML,ML,DMIN,TOL)
                        do I=1,ML
                                AZIP(I,I,POSA) = S(I)
                        end do 
                        U(POSUV:POSUV+ML,:)=UAUX(:,1:K)
                        V(POSUV:POSUV+ML,:)=VAUX(:,1:K)
                end do
                deallocate(S,AAUX,UAUX,VAUX)
        end do
        write(*,*) 'MATRICE NON COMPRESSA', storage_size(A)*sizeof(A)
        write(*,*) 'MATRICE COMPRESSA', storage_size(AZIP)*(sizeof(AZIP)+sizeof(U)+sizeof(V))
END SUBROUTINE
SUBROUTINE SVD(A,U,V,S,M,N,K,TOL)
! La subroutine è definita dai paramtri:
!       integer m,n,k
!       double precision, tol
!       double precision, dimension (m,n), A
!       double precision, dimension (m,k), U
!       double precision, dimension (n,k), V
!       double precision, dimension (k),   S
! prende in input k=min(m,n);
! calcola la svd di A e cerca il massimo k tale che il k-esimo sing.val
! sia maggiore di tol, sovrascrive A con A_k;
! sovrascrive k con il rango numerico ripetto a tol

        IMPLICIT NONE
!variabili in input        
        INTEGER :: M, N, K
        DOUBLE PRECISION :: TOL
        DOUBLE PRECISION :: A(M,N), U(M,K), S(K), V(N,K)

!scalari locali
        CHARACTER :: JOBU, JOBVT
        LOGICAL :: contr=.true.
        INTEGER :: LDA, LDU, LDVT, LWORK, INFO, DMIN, I, J
!array locali
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK
!        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SIG
        DOUBLE PRECISION :: VT(K,N)

!inizializzazione parametri
!controlli
        if (K .ne. min(M,N)) then
                write(*,*) 'K deve esere min(m,n)'
                stop -1
        end if
        if (tol.le.0) then
                write(*,*) 'TOL deve essere positiva'
                stop -2
        end if
        DMIN=K
        JOBU='S'
        JOBVT='S'
        LDA=M
        LDU=M
        LDVT=DMIN
        LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
        ALLOCATE(WORK(LWORK))

        CALL dgesvd(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO)

        V=TRANSPOSE(VT) !tutti i vettori singolari destri

! calcolo rango numerico rispetto a tol

!        do I=0,DMIN-1
!                if (S(DMIN-I).ge. TOL) then
!                        K=DMIN -I
!                        contr=.false. 
!                        exit
!                end if
!        end do
!        if (contr) then ! mi accerto che tol piccolo abbastanza
!                K=0;
!                A=0;
!                write(*,*) 'Matrice nulla con TOL assegnata, provare a diminuire TOL'
!                stop -3
!        end if
!        ALLOCATE(SIG(K,K))    
!        DO I=1,K    ! costruisco matrice valori singolari
!                DO J=1,K
!                        IF (I==J) THEN
!                                SIG(I,J)=S(I)
!                        ELSE
!                                SIG(I,J)=0
!                        ENDIF
!                END DO
!        END DO 

!sovrascrivo A con A_k
!        A=MATMUL(U(:,1:K),MATMUL(SIG,VT(1:K,:))) ! poco performante uncomment for testing

        DEALLOCATE(WORK)!,SIG) !libero memoria
END SUBROUTINE
END MODULE HMOD
PROGRAM test
        use HMOD
        integer :: K, N, LL
        double precision, allocatable, dimension(:,:) :: A, U, V
        double precision, allocatable, dimension(:,:,:) :: AZIP
        
        N=2**8 * 3
        LL=8
        K=3
        allocate(A(N,N), AZIP(K,K,2**(LL+1)-2), U(N*LL,K), V(N*LL,K))
!        A(1:3:N,10)=1
!        A(1:5:N,30)=1
!        A(1:7:N,100)=1
        call random_number(A)
        A=A/norm2(A)
        call HODLR(A,N,K,LL,AZIP,U,V)
END PROGRAM
