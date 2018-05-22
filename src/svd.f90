PROGRAM testsvd
        implicit none
        integer :: m=300, n=400, dmin, k, i, j
        double precision :: tol, err
        double precision, allocatable, dimension(:,:) :: U, V, SIG, A, Aold
        double precision, allocatable, dimension(:) :: S

dmin=min(m,n)
k=dmin
tol=1.d-10
allocate(A(m,n),Aold(m,n), U(m,dmin), V(n,dmin), SIG(dmin,dmin),S(dmin))

call random_number(A)
A=A/norm2(A)
!A(1:140,:)=A(141:280,:)
A(300,:)=tol**2 *A(300,:)
Aold=A

call svd(A,U,V,S,m,n,k,tol)

err=norm2(Aold - A)
WRITE(*,*) 'tolleranza= ', tol
WRITE(*,*) 'rango= ', k
WRITE(*,*) 'sigma_k+1', S(MIN(k+1,DMIN))

deallocate(A,Aold,U,V,SIG,S)
END PROGRAM

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
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SIG
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

do I=0,DMIN-1
        if (S(DMIN-I).ge. TOL) then
                K=DMIN -I
                contr=.false. 
                exit
        end if
end do

if (contr) then ! mi accerto che tol piccolo abbastanza
        K=0;
        A=0;
        write(*,*) 'Matrice nulla con TOL assegnata, provare a diminuire TOL'
        stop -3
end if
ALLOCATE(SIG(K,K))    
DO I=1,K    ! costruisco matrice valori singolari
        DO J=1,K
          IF (I==J) THEN
             SIG(I,J)=S(I)
          ELSE
             SIG(I,J)=0
          ENDIF
        END DO
END DO 

!sovrascrivo A con A_k
A=MATMUL(U(:,1:K),MATMUL(SIG,VT(1:K,:))) ! poco performante

DEALLOCATE(WORK,SIG) !libero memoria
END SUBROUTINE
