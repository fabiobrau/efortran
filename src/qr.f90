PROGRAM testqr
        implicit none
        integer :: m=1000, n=1000, K, DMIN, i, j
        integer, allocatable, dimension(:) :: JPVT
        double precision :: TOL
        double precision, allocatable, dimension(:,:) :: A, AOLD, Q, R


        DMIN=min(m,n)
        K=DMIN 
        TOL=1.d-10
        allocate(A(m,n), AOLD(m,n), Q(m,dmin), R(dmin,n), JPVT(n))
        call random_number(A)
        A=A/norm2(A)
        A(1:149,:)=A(150:298,:)
        A(300,:)=tol**2 *A(300,:)
        AOLD=A
        call qr(M,N,A,K,Q,R,JPVT,TOL)
        write(*,*) norm2(AOLD(:,JPVT) - MATMUL(Q,R))
        write(*,*) 'Rango = ', K
        deallocate(A, AOLD, Q, R, JPVT)
END PROGRAM
SUBROUTINE qr(M,N,A,K,Q,R,JPVT,TOL)
        implicit none
! La subroutine prende in input un set di vettori e sfruttando il QR
! con pivoting mando in output il set ortogonalizzato in Q 
! in modo che i primi k vettori siano lin.indipendenti
! K[in] deve essere inizalizzato come min(m,n)
! K[out] indice corrispondente al più piccolo R(I,I)<TOL, se non esiste K=DMIN
! A viene sovrascritta in modo che la parte triangolare superiore sia R
! il vettore tau non esce in output quindi non è possibile ricostruire
! le singole riflessioni poi magari si aggiusta.

! Variabili di input 
        integer :: M, N, K
        integer :: JPVT(N)
        double precision :: TOL
        double precision :: A(M,K), Q(M,K), R(K,N)
!variabili locali
        logical :: contr=.true.
        integer :: LDA, LWORK, INFO, DMIN, I, J
        double precision, allocatable, dimension(:) :: WORK
        double precision ::  ID(M,M), V(M,1), QL(M,M) 
        double precision :: TAU(K), RDIAG(K)
!controllo
if (K .ne. min(M,N)) then
        write(*,*) 'Inizializzare K come min(M,N)'
        stop -1
end if
!inizializzazione
DMIN=K
LDA=M
LWORK=3*N+1 !lascio a dgeqp3 la scelta ottimale
allocate(WORK(LWORK))
call dgeqp3(M,N,A,LDA,JPVT,TAU,WORK,LWORK,INFO)
!costruisco Q, R
do I=1,M ! costruisco matrice identica di taglia m
        do J=1,M
                if (J==I) then
                        ID(I,J)=1
                else
                        ID(I,J)=0
                end if
        end do
end do
R=0
do I=1,DMIN ! costruisco R
        do J=1,N
                if (I .le. J) then
                        R(I,J)=A(I,J)
                else
                        R(I,J)=0
                end if
        end do
end do
! Costruisco la matrice ortogonale quadrata
V(1,1)=1
V(2:M,1)=A(2:M,1)
QL=ID - matmul(TAU(1)*V,transpose(V))
do I=2,DMIN !costruisco Q
        V(1:I-1,1)=0
        V(I,1)=1
        if (I < DMIN) then
                V(I+1:M,1)=A(I+1:M,I)
        end if
        QL=QL - TAU(I)*matmul(matmul(QL,V),transpose(V)) ! sfrutta la struct.
end do
! costruisco Q di dimensione M X DMIN
if (DMIN==M) then
        Q=QL
else
        Q=QL(:,1:DMIN) ! N.B non è un appross. vedi teorema di Gu
end if

do I=0,DMIN-1
        if (abs(R(DMIN-I,DMIN-I)) > TOL) then !! attenzione da verificare
                K=DMIN -I
                contr=.false.
                exit
        end if
end do
! mi accerto che esista il minimo
if (contr) then
        K=0
        write(*,*) 'Matrice nulla con TOL asseganta, provare a diminuire TOL'
        stop -2
end if
deallocate(WORK)
END SUBROUTINE
