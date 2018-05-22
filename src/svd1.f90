PROGRAM testsvd
        implicit none
        integer :: m=30, n=40, dmin, k, i, j
        double precision :: tol
        double precision, allocatable, dimension(:,:) :: U, V, SIG, A, Aold
        double precision, allocatable, dimension(:) :: S

dmin=min(m,n)
k=16
allocate(A(m,n),Aold(m,n),U(m,k), V(n,k), SIG(k,k),S(k))
call random_number(A)
A(1:15,:)=A(16:30,:)
Aold=A

call svd(A,U,V,S,m,n,k,tol)


DO I=1,k
        DO J=1,k
          IF (I==J) THEN
             SIG(I,J)=S(I)
          ELSE
             SIG(I,J)=0
          ENDIF
        END DO
END DO

tol=norm2(Aold - A)
WRITE(*,*) 'errore '
WRITE(*,*) tol

WRITE(*,*) ' norma di U'
WRITE(*,*) norm2(U)

deallocate(A,Aold,U,V,SIG,S)
END PROGRAM

SUBROUTINE SVD(A,U,V,S,M,N,K,TOL)
        IMPLICIT NONE
!variabili in input        
        INTEGER :: M, N, K
        DOUBLE PRECISION :: TOL
        DOUBLE PRECISION :: A(M,N), U(M,K), S(K), V(N,K)
!subroutine dgesvd       (       
!        character       JOBU,
!        character       JOBVT,
!        integer         M,
!        integer         N,
!        double precision, dimension( lda, * )   A,
!        integer         LDA,
!        double precision, dimension( * )        S,
!        double precision, dimension( ldu, * )   U,
!        integer         LDU,
!        double precision, dimension( ldvt, * )  VT,
!        integer         LDVT,
!        double precision, dimension( * )        WORK,
!        integer         LWORK,
!        integer         INFO 
!        )               

!scalari locali
        CHARACTER :: JOBU, JOBVT
        INTEGER :: LDA, LDU, LDVT, LWORK, INFO, DMIN, I, J
!array locali
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: UL, VTL, VL, SIG
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK, SL

!inizializzazione parametri
DMIN=min(M,N)
JOBU='S'
JOBVT='S'
LDA=M
LDU=M
LDVT=DMIN
LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
ALLOCATE(WORK(LWORK),UL(M,DMIN),VTL(DMIN,N),VL(N,DMIN),SL(DMIN),SIG(K,K))
!controlli
IF (K.le.0 .OR. K>DMIN) THEN
        WRITE(*,*) "K DEVE ESSERE COMPRESO TRA 0 E DMIN"
        STOP -1
END IF

CALL dgesvd(JOBU,JOBVT,M,N,A,LDA,SL,UL,LDU,VTL,LDVT,WORK,LWORK,INFO)

VL=TRANSPOSE(VTL) !tutti i vettori singolari destri
U=UL(:,1:K) ! solo i primi k vettori singolari sinistri
V=VL(:,1:K) ! solo i primi k vettori singolari destri
S=SL(1:K)   ! solo i primi k valori singolari
DO I=1,k    ! costruisco matrice valori singolari
        DO J=1,k
          IF (I==J) THEN
             SIG(I,J)=S(I)
          ELSE
             SIG(I,J)=0
          ENDIF
        END DO
END DO
A=MATMUL(U,MATMUL(SIG,TRANSPOSE(V))) ! appross di A con rango k
DEALLOCATE(WORK,UL,VTL,VL,SL,SIG) ! libero memoria
END SUBROUTINE
