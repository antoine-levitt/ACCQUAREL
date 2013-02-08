MODULE constants
  DOUBLE PRECISION,PARAMETER :: PI=3.14159265358979323846D0
! speed of light in the vacuum in atomic units (for the relativistic case)
! Note : One has $c=\frac{e^2h_e}{\hbar\alpha}$, where $\alpha$ is the fine structure constant, $c$ is the speed of light in the vacuum, $e$ is the elementary charge, $\hbar$ is the reduced Planck constant and $k_e$ is the Coulomb constant. In Hartree atomic units, the numerical values of the electron mass, the elementary charge, the reduced Planck constant and the Coulomb constant are all unity by definition, so that $c=\alpha^{-1}$. The value chosen here is the one recommended in: P. J. Mohr, B. N. Taylor, and D. B. Newell, CODATA recommended values of the fundamental physical constants: 2006.
  DOUBLE PRECISION,PARAMETER :: SPEED_OF_LIGHT=137.035999967994D0
END MODULE

MODULE random
CONTAINS

SUBROUTINE INIT_RANDOM_SEED()
! Initialization of the seed of the pseudorandom number generator used by RANDOM_NUMBER based on the system's time (this routine must be called once).
  IMPLICIT NONE
  INTEGER :: I,N,CLOCK
  INTEGER,DIMENSION(:),ALLOCATABLE :: SEED

  CALL RANDOM_SEED(SIZE=N)
  ALLOCATE(SEED(N))

  CALL SYSTEM_CLOCK(COUNT=CLOCK)

  SEED=CLOCK+37*(/(I-1,I=1,N)/)
  CALL RANDOM_SEED(PUT=SEED)

  DEALLOCATE(SEED)
END SUBROUTINE INIT_RANDOM_SEED

FUNCTION GET_RANDOM(N) RESULT(RANDOM_ARRAY)
! Function that returns an array of random numbers of size N in (0, 1)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  REAL,DIMENSION(N) :: RANDOM_ARRAY

  CALL RANDOM_NUMBER(RANDOM_ARRAY)
END FUNCTION GET_RANDOM
END MODULE random

MODULE matrix_tools
INTERFACE PACK
  MODULE PROCEDURE PACK_symmetric,PACK_hermitian
END INTERFACE

INTERFACE UNPACK
  MODULE PROCEDURE UNPACK_symmetric,UNPACK_hermitian
END INTERFACE

INTERFACE ABA
  MODULE PROCEDURE ABA_symmetric,ABA_hermitian
END INTERFACE

INTERFACE ABCBA
  MODULE PROCEDURE ABCBA_symmetric,ABCBA_hermitian
END INTERFACE

INTERFACE ABC_CBA
  MODULE PROCEDURE ABC_CBA_symmetric,ABC_CBA_hermitian
END INTERFACE

INTERFACE BUILD_BLOCK_DIAGONAL
   MODULE PROCEDURE BUILD_BLOCK_DIAGONAL_symmetric
END INTERFACE BUILD_BLOCK_DIAGONAL

INTERFACE FINNERPRODUCT
  MODULE PROCEDURE FROBENIUSINNERPRODUCT_real,FROBENIUSINNERPRODUCT_complex
END INTERFACE

INTERFACE NORM
  MODULE PROCEDURE NORM_real,NORM_complex,NORM_symmetric,NORM_hermitian
END INTERFACE

INTERFACE INVERSE
  MODULE PROCEDURE INVERSE_real,INVERSE_complex,INVERSE_symmetric,INVERSE_hermitian
END INTERFACE

INTERFACE SQUARE_ROOT
  MODULE PROCEDURE SQUARE_ROOT_symmetric,SQUARE_ROOT_hermitian
END INTERFACE

INTERFACE EXPONENTIAL
   MODULE PROCEDURE EXPONENTIAL_real,EXPONENTIAL_complex
END INTERFACE EXPONENTIAL

INTERFACE TRACE
  MODULE PROCEDURE TRACE_symmetric,TRACE_hermitian
END INTERFACE

INTERFACE TRACEOFPRODUCT
  MODULE PROCEDURE TRACEOFPRODUCT_real,TRACEOFPRODUCT_complex,TRACEOFPRODUCT_symmetric,TRACEOFPRODUCT_hermitian
END INTERFACE

INTERFACE EIGENSOLVER
   MODULE PROCEDURE EIGENSOLVER_symmetric_prefactorized,EIGENSOLVER_hermitian_prefactorized
END INTERFACE

INTERFACE COMMUTATOR
   MODULE PROCEDURE COMMUTATOR_symmetric,COMMUTATOR_hermitian
END INTERFACE

INTERFACE PRINTMATRIX
   MODULE PROCEDURE PRINTMATRIX_symmetric,PRINTMATRIX_hermitian,PRINTMATRIX_complex,PRINTMATRIX_real
END INTERFACE

INTERFACE READMATRIX
   MODULE PROCEDURE READMATRIX_complex,READMATRIX_real
END INTERFACE READMATRIX

CONTAINS

! handling of symmetric and hermitian matrices stored in packed form.

FUNCTION PACK_symmetric(A,N) RESULT (PA)
! Function that stores the upper triangular part of a symmetric matrix in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PA

  INTEGER :: IJ,I,J

  IJ=0
  DO J=1,N
     DO I=1,J
        IJ=IJ+1
        PA(IJ)=(A(I,J)+A(J,I))/2.D0
     END DO
  END DO
END FUNCTION PACK_symmetric

FUNCTION PACK_hermitian(A,N) RESULT (PA)
! Function that stores the upper triangular part of a hermitian matrix in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N,N),INTENT(IN) :: A
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PA

  INTEGER :: IJ,I,J

  IJ=0
  DO J=1,N
     DO I=1,J
        IJ=IJ+1
        PA(IJ)=A(I,J)
        PA(IJ)=(A(I,J)+conjg(A(J,I)))/2.D0
     END DO
  END DO
END FUNCTION PACK_hermitian

FUNCTION UNPACK_symmetric(PA,N) RESULT (A)
! Function that unpacks a symmetric matrix which upper triangular part is stored in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE PRECISION,DIMENSION(N,N) :: A

  INTEGER :: I,J

  DO I=1,N
     J=I*(I-1)/2
     A(1:I,I)=PA(1+J:I+J)
     A(I,1:I-1)=PA(1+J:I-1+J)
  END DO
END FUNCTION UNPACK_symmetric

FUNCTION UNPACK_hermitian(PA,N) RESULT (A)
! Function that unpacks a hermitian matrix which upper triangular part is stored in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE COMPLEX,DIMENSION(N,N) :: A

  INTEGER :: I,J

  DO I=1,N
     J=I*(I-1)/2
     A(1:I,I)=PA(1+J:I+J)
     A(I,1:I-1)=CONJG(PA(1+J:I-1+J))
  END DO
END FUNCTION UNPACK_hermitian

FUNCTION ABA_symmetric(PA,PB,N) RESULT (PC)
! Function that computes the product ABA, where A and B are two symmetric matrices, which upper triangular parts are stored in packed format (the resulting matrix being also stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PC

  INTEGER :: I,J,K,L,IJ

  PC=(0.D0,0.D0)
  IJ=0
  DO J=1,N
     DO I=1,J
        IJ=IJ+1
        DO K=1,I
           DO L=1,K
              PC(IJ)=PC(IJ)+PA(K+(I-1)*I/2)*PB(L+(K-1)*K/2)*PA(L+(J-1)*J/2)
           END DO
           DO L=K+1,J
              PC(IJ)=PC(IJ)+PA(K+(I-1)*I/2)*PB(K+(L-1)*L/2)*PA(L+(J-1)*J/2)
           END DO
           DO L=J+1,N
              PC(IJ)=PC(IJ)+PA(K+(I-1)*I/2)*PB(K+(L-1)*L/2)*PA(J+(L-1)*L/2)
           END DO
        END DO
        DO K=I+1,J
           DO L=1,K
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(L+(K-1)*K/2)*PA(L+(J-1)*J/2)
           END DO
           DO L=K+1,J
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(K+(L-1)*L/2)*PA(L+(J-1)*J/2)
           END DO
           DO L=J+1,N
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(K+(L-1)*L/2)*PA(J+(L-1)*L/2)
           END DO
        END DO
        DO K=J+1,N
           DO L=1,J
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(L+(K-1)*K/2)*PA(L+(J-1)*J/2)
           END DO
           DO L=J+1,K
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(L+(K-1)*K/2)*PA(J+(L-1)*L/2)
           END DO
           DO L=K+1,N
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(K+(L-1)*L/2)*PA(J+(L-1)*L/2)
           END DO
        END DO
     END DO
  END DO
END FUNCTION ABA_symmetric

FUNCTION ABA_hermitian(PA,PB,N) RESULT (PC)
! Function that computes the product ABA, where A and B are two hermitian matrices, which upper triangular parts are stored in packed format (the resulting matrix being also stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PC

  INTEGER :: I,J,K,L,IJ

  PC=(0.D0,0.D0)
  IJ=0
  DO J=1,N
     DO I=1,J
        IJ=IJ+1
        DO K=1,I
           DO L=1,K
              PC(IJ)=PC(IJ)+CONJG(PA(K+(I-1)*I/2))*CONJG(PB(L+(K-1)*K/2))*PA(L+(J-1)*J/2)
           END DO
           DO L=K+1,J
              PC(IJ)=PC(IJ)+CONJG(PA(K+(I-1)*I/2))*PB(K+(L-1)*L/2)*PA(L+(J-1)*J/2)
           END DO
           DO L=J+1,N
              PC(IJ)=PC(IJ)+CONJG(PA(K+(I-1)*I/2))*PB(K+(L-1)*L/2)*CONJG(PA(J+(L-1)*L/2))
           END DO
        END DO
        DO K=I+1,J
           DO L=1,K
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*CONJG(PB(L+(K-1)*K/2))*PA(L+(J-1)*J/2)
           END DO
           DO L=K+1,J
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(K+(L-1)*L/2)*PA(L+(J-1)*J/2)
           END DO
           DO L=J+1,N
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(K+(L-1)*L/2)*CONJG(PA(J+(L-1)*L/2))
           END DO
        END DO
        DO K=J+1,N
           DO L=1,J
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*CONJG(PB(L+(K-1)*K/2))*PA(L+(J-1)*J/2)
           END DO
           DO L=J+1,K
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*CONJG(PB(L+(K-1)*K/2))*CONJG(PA(J+(L-1)*L/2))
           END DO
           DO L=K+1,N
              PC(IJ)=PC(IJ)+PA(I+(K-1)*K/2)*PB(K+(L-1)*L/2)*CONJG(PA(J+(L-1)*L/2))
           END DO
        END DO
     END DO
  END DO

  DO I=1,N
     PC(I*(I+1)/2) = REAL(PC(I*(I+1)/2))
  END DO
END FUNCTION ABA_hermitian

FUNCTION ABCBA_symmetric(PA,PB,PC,N) RESULT (PD)
! Function that computes the product ABCBA, where A, B, and C are three symmetric matrices, which upper triangular parts are stored in packed format (the resulting matrix being also stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB,PC
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PD

  DOUBLE PRECISION,DIMENSION(N,N) :: A,B,C

  A=UNPACK(PA,N) ; B=UNPACK(PB,N) ; C=UNPACK(PC,N)
  PD=PACK(MATMUL(A,MATMUL(B,MATMUL(C,MATMUL(B,A)))),N)
END FUNCTION ABCBA_symmetric

FUNCTION ABCBA_hermitian(PA,PB,PC,N) RESULT (PD)
! Function that computes the product ABCBA, where A, B, and C are three hermitian matrices, which upper triangular parts are stored in packed format (the resulting matrix being also stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB,PC
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PD

  PD=ABA(PA,ABA(PB,PC,N),N)
END FUNCTION ABCBA_hermitian

FUNCTION ABC_CBA_symmetric(PA,PB,PC,N) RESULT (PD)
! Function that computes the sum ABC+CBA, where A, B, and C are three symmetric matrices, which upper triangular parts are stored in packed format (the resulting matrix being also stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB,PC
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PD

  DOUBLE PRECISION,DIMENSION(N,N) :: A,B,C

  A=UNPACK(PA,N) ; B=UNPACK(PB,N) ; C=UNPACK(PC,N)
  PD=PACK(MATMUL(A,MATMUL(B,C))+MATMUL(C,MATMUL(B,A)),N)
END FUNCTION ABC_CBA_symmetric

FUNCTION ABC_CBA_hermitian(PA,PB,PC,N) RESULT (PD)
! Function that computes the sum ABC+CBA, where A, B, and C are three hermitian matrices, which upper triangular parts are stored in packed format (the resulting matrix being also stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB,PC
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PD

  DOUBLE COMPLEX,DIMENSION(N,N) :: A,B,C

  A=UNPACK(PA,N) ; B=UNPACK(PB,N) ; C=UNPACK(PC,N)
  PD=PACK(MATMUL(A,MATMUL(B,C))+MATMUL(C,MATMUL(B,A)),N)
END FUNCTION ABC_CBA_hermitian


SUBROUTINE BUILD_BLOCK_DIAGONAL_symmetric(PB,PA,N)
! Subroutine that forms the block-diagonal symmetric matrix B of order 2N from a symmetric matrix A of size N, both of which have their upper triangular parts are stored in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION, DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE PRECISION, DIMENSION(N*2*(N*2+1)/2),INTENT(OUT) :: PB

  DOUBLE PRECISION, DIMENSION(2*N,2*N) :: B

  B=0.D0
  B(1:N,1:N)=UNPACK(PA,N) ; B(N+1:2*N,N+1:2*N)=B(1:N,1:N)
  PB=PACK(B,2*N)
END SUBROUTINE BUILD_BLOCK_DIAGONAL_symmetric

! diverse linear algebra routines

FUNCTION INVERSE_real(A,N) RESULT(INVA)
! Function that computes the inverse of a square real matrix.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N,N) :: A
  DOUBLE PRECISION,DIMENSION(N,N) :: INVA

  INTEGER :: INFO
  INTEGER,DIMENSION(N) :: IPIV
  DOUBLE PRECISION,DIMENSION(N) :: WORK

  INVA=A
  CALL DGETRF(N,N,INVA,N,IPIV,INFO)
  IF (INFO/=0) GOTO 1
  CALL DGETRI(N,INVA,N,IPIV,WORK,N,INFO)
  IF (INFO/=0) GOTO 2
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DGETRF: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DGETRF: U(',INFO,',',INFO,') is exactly zero. The factorization &
    &has been completed, but the factor U is exactly singular, and division by zero will &
    &occur if it is used to solve a system of equations'
  END IF
  GO TO 3
2 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DGETRI: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DGETRI: U(',INFO,',',INFO,') is exactly zero; &
    &the matrix is singular and its inverse could not be computed'
  END IF
3 WRITE(*,*)'(called from function INVERSE)'
  STOP
END FUNCTION INVERSE_real

FUNCTION INVERSE_complex(A,N) RESULT(INVA)
! Function that computes the inverse of a square complex matrix.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N,N) :: A
  DOUBLE COMPLEX,DIMENSION(N,N) :: INVA

  INTEGER :: INFO
  INTEGER,DIMENSION(N) :: IPIV
  DOUBLE COMPLEX,DIMENSION(N) :: WORK

  INVA=A
  CALL ZGETRF(N,N,INVA,N,IPIV,INFO)
  IF (INFO/=0) GOTO 1
  CALL ZGETRI(N,INVA,N,IPIV,WORK,N,INFO)
  IF (INFO/=0) GOTO 2
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZGETRF: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZGETRF: U(',INFO,',',INFO,') is exactly zero. The factorization &
    &has been completed, but the factor U is exactly singular, and division by zero will &
    &occur if it is used to solve a system of equations'
  END IF
  GO TO 3
2 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZGETRI: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZGETRI: U(',INFO,',',INFO,') is exactly zero; &
    &the matrix is singular and its inverse could not be computed'
  END IF
3 WRITE(*,*)'(called from function INVERSE)'
  STOP
END FUNCTION INVERSE_complex

FUNCTION INVERSE_symmetric(PA,N) RESULT(PINVA)
! Function that computes the inverse of a symmetric matrix which upper triangular part is stored in packed format (its inverse being stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PINVA

  INTEGER :: INFO
  INTEGER,DIMENSION(N) :: IPIV
  DOUBLE PRECISION,DIMENSION(N) :: WORK
 
  PINVA=PA
  CALL DSPTRF('U',N,PINVA,IPIV,INFO)
  IF (INFO/=0) GOTO 1
  CALL DSPTRI('U',N,PINVA,IPIV,WORK,INFO)
  IF (INFO/=0) GOTO 2
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DSPTRF: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DSPTRF: D(',INFO,',',INFO,') is exactly zero. The factorization &
    &has been completed, but the block diagonal matrix D is exactly singular, and division &
    &by zero will occur if it is used to solve a system of equations'
  END IF
  GO TO 3
2 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DSPTRI: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DSPTRI: D(',INFO,',',INFO,') is exactly zero; &
    &the matrix is singular and its inverse could not be computed'
  END IF
3 WRITE(*,*)'(called from function INVERSE)'
  STOP
END FUNCTION INVERSE_symmetric

FUNCTION INVERSE_hermitian(PA,N) RESULT(PINVA)
! Function that computes the inverse of an hermitian matrix which upper triangular part is stored in packed format (its inverse being stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PINVA

  INTEGER :: INFO
  INTEGER,DIMENSION(N) :: IPIV
  DOUBLE COMPLEX,DIMENSION(N) :: WORK
 
  PINVA=PA
  CALL ZHPTRF('U',N,PINVA,IPIV,INFO)
  IF (INFO/=0) GOTO 1
  CALL ZHPTRI('U',N,PINVA,IPIV,WORK,INFO)
  IF (INFO/=0) GOTO 2
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZHPTRF: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZHPTRF: D(',INFO,',',INFO,') is exactly zero. The factorization &
    &has been completed, but the block diagonal matrix D is exactly singular, and division &
    &by zero will occur if it is used to solve a system of equations'
  END IF
  GO TO 3
2 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZHPTRI: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZHPTRI: D(',INFO,',',INFO,') is exactly zero; &
    &the matrix is singular and its inverse could not be computed'
  END IF
3 WRITE(*,*)'(called from function INVERSE)'
  STOP
END FUNCTION INVERSE_hermitian

FUNCTION SQUARE_ROOT_symmetric(PA,N) RESULT(PSQRA)
! Function that computes the square root of a symmetric, positive-definite matrix which upper triangular part is stored in packed format (its square root being stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PSQRA

  INTEGER :: INFO,I
  INTEGER,DIMENSION(N) :: IPIV
  DOUBLE PRECISION,DIMENSION(N) :: EIG
  DOUBLE PRECISION,DIMENSION(3*N) :: WORK
  DOUBLE PRECISION,DIMENSION(N,N) :: EIGVEC,M

  PSQRA=PA
  CALL DSPEV('V','U',N,PSQRA,EIG,EIGVEC,N,WORK,INFO)
  IF (INFO/=0) GOTO 1
  FORALL(I=1:N) M(:,I)=SQRT(EIG(I))*EIGVEC(:,I)
  PSQRA=PACK(MATMUL(M,TRANSPOSE(EIGVEC)),N)
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DSPEV: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DSPEV: the algorithm failed to converge; ',INFO, &
    &'off-diagonal elements of an intermediate tridiagonal form did not converge to zero'
  END IF
  WRITE(*,*)'(called from function SQUARE_ROOT)'
  STOP
END FUNCTION SQUARE_ROOT_symmetric

FUNCTION SQUARE_ROOT_hermitian(PA,N) RESULT(PSQRA)
! Function that computes the square root of an hermitian, positive-definite matrix which upper triangular part is stored in packed format (its square root being stored as such).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PSQRA

  INTEGER :: INFO,I
  INTEGER,DIMENSION(N) :: IPIV
  DOUBLE PRECISION,DIMENSION(N) :: EIG
  DOUBLE PRECISION,DIMENSION(3*N-2) :: RWORK
  DOUBLE COMPLEX,DIMENSION(2*N-1) :: WORK
  DOUBLE COMPLEX,DIMENSION(N,N) :: EIGVEC,M

  PSQRA=PA
  CALL ZHPEV('V','U',N,PSQRA,EIG,EIGVEC,N,WORK,RWORK,INFO)
  IF (INFO/=0) GOTO 1
  FORALL(I=1:N) M(:,I)=SQRT(EIG(I))*EIGVEC(:,I)
  PSQRA=PACK(MATMUL(M,TRANSPOSE(CONJG(EIGVEC))),N)
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZHPEV: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZHPEV: the algorithm failed to converge; ',INFO, &
    &'off-diagonal elements of an intermediate tridiagonal form did not converge to zero'
  END IF
  WRITE(*,*)'(called from function SQUARE_ROOT)'
  STOP
END FUNCTION SQUARE_ROOT_hermitian

FUNCTION EXPONENTIAL_real(T,A,N) result(EXPTA)
! Function that computes the matrix exponential exp(tA), where A is an N-by-N real matrix and t is a real scalar, using the Expokit software package (http://www.maths.uq.edu.au/expokit/).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,INTENT(IN) :: T
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A
  DOUBLE PRECISION,DIMENSION(N,N) :: EXPTA
  
  INTEGER :: IEXP,NS,IFLAG
  INTEGER,DIMENSION(N)  :: IWSP
  INTEGER,PARAMETER :: IDEG=6
  DOUBLE PRECISION,DIMENSION(4*N*N+IDEG+1) :: WSP

  CALL DGPADM(IDEG,N,T,A,N,WSP,SIZE(WSP,1),IWSP,IEXP,NS,IFLAG)
  IF (IFLAG/=0) GO TO 1
  EXPTA=RESHAPE(WSP(IEXP:IEXP+N*N-1),SHAPE(EXPTA))
  RETURN
1 WRITE(*,*)'Subroutine DGPADM: there is a problem'
  WRITE(*,*)'(called from function EXPONENTIAL)'
  STOP
END FUNCTION EXPONENTIAL_real

FUNCTION EXPONENTIAL_complex(T,A,N) result(EXPTA)
! Function that computes the matrix exponential exp(tA), where A is an N-by-N complex matrix and t is a real scalar, using the Expokit software package (http://www.maths.uq.edu.au/expokit/).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,INTENT(IN) :: T
  DOUBLE COMPLEX,DIMENSION(N,N),INTENT(IN) :: A
  DOUBLE COMPLEX,DIMENSION(N,N) :: EXPTA

  INTEGER :: IEXP,NS,IFLAG
  INTEGER,DIMENSION(N)  :: IWSP
  INTEGER,PARAMETER :: IDEG=6
  DOUBLE COMPLEX,DIMENSION(4*N*N+IDEG+1) :: WSP

  CALL ZGPADM(IDEG,N,T,A,N,WSP,SIZE(WSP,1),IWSP,IEXP,NS,IFLAG)
  IF (IFLAG/=0) GO TO 1
  EXPTA=RESHAPE(WSP(IEXP:IEXP+N*N-1),SHAPE(EXPTA))
  RETURN
1 WRITE(*,*)'Subroutine ZGPADM: there is a problem'
  WRITE(*,*)'(called from function EXPONENTIAL)'
  STOP
END FUNCTION EXPONENTIAL_complex

FUNCTION TRACE_real(A,N) RESULT (TRACE)
  ! Function that computes the trace of a symmetric matrix, which upper triangular part is stored in packed format.
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A
  DOUBLE PRECISION :: TRACE

  INTEGER :: I

  TRACE=0.D0
  DO I=1,N
     TRACE=TRACE+A(I,I)
  END DO
END FUNCTION TRACE_real

FUNCTION TRACE_symmetric(PA,N) RESULT (TRACE)
! Function that computes the trace of a symmetric matrix, which upper triangular part is stored in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE PRECISION :: TRACE

  INTEGER :: I

  TRACE=0.D0
  DO I=1,N
     TRACE=TRACE+PA((I+1)*I/2)
  END DO
END FUNCTION TRACE_symmetric

FUNCTION TRACE_hermitian(PA,N) RESULT (TRACE)
! Function that computes the trace of a hermitian matrix, which upper triangular part is stored in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE COMPLEX :: TRACE

  INTEGER :: I

  TRACE=(0.D0,0.D0)
  DO I=1,N
     TRACE=TRACE+PA((I+1)*I/2)
  END DO
END FUNCTION TRACE_hermitian

FUNCTION TRACEOFPRODUCT_real(A,B,N) RESULT (TRACE)
! Function that computes the trace of the product of two square matrices A and B.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A,B
  DOUBLE PRECISION :: TRACE

  INTEGER :: I,J

  TRACE=0.D0
  DO I=1,N
     DO J=1,N
        TRACE=TRACE+A(I,J)*B(J,I)
     END DO
  END DO
END FUNCTION TRACEOFPRODUCT_real

FUNCTION TRACEOFPRODUCT_complex(A,B,N) RESULT (TRACE)
! Function that computes the trace of the product of two square complex matrices A and B.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N,N),INTENT(IN) :: A,B
  DOUBLE COMPLEX :: TRACE

  INTEGER :: I,J

  TRACE=(0.D0,0.D0)
  DO I=1,N
     DO J=1,N
        TRACE=TRACE+A(I,J)*B(J,I)
     END DO
  END DO
END FUNCTION TRACEOFPRODUCT_complex

FUNCTION TRACEOFPRODUCT_symmetric(PA,PB,N) RESULT (TRACE)
! Function that computes the trace of the product of two symmetric matrices A and B, which upper triangular parts are stored in packed format. 
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB
  DOUBLE PRECISION :: TRACE

  INTEGER :: I,J,IJ,JI

  TRACE=0.D0
  DO J=1,N
     DO I=1,J
        IJ=I+(J-1)*J/2
        TRACE=TRACE+PA(IJ)*PB(IJ)
     END DO
     DO I=J+1,N
        JI=J+(I-1)*I/2
        TRACE=TRACE+PA(JI)*PB(JI)
     END DO
  END DO
END FUNCTION TRACEOFPRODUCT_symmetric

FUNCTION TRACEOFPRODUCT_hermitian(PA,PB,N) RESULT (TRACE)
! Function that computes the trace of the product of two hermitian matrices A and B, which upper triangular parts are stored in packed format.
  IMPLICIT NONE 
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB
  DOUBLE COMPLEX :: TRACE

  INTEGER :: I,J,IJ,JI

  TRACE=(0.D0,0.D0)
  DO J=1,N
     DO I=1,J
        IJ=I+(J-1)*J/2
        TRACE=TRACE+PA(IJ)*CONJG(PB(IJ))
     END DO
     DO I=J+1,N
        JI=J+(I-1)*I/2
        TRACE=TRACE+CONJG(PA(JI))*PB(JI)
     END DO
  END DO
END FUNCTION TRACEOFPRODUCT_hermitian

FUNCTION FROBENIUSINNERPRODUCT_real(A,B,N) RESULT (FIP)
! Function that computes the Frobenius inner product between two square real matrices (i.e. $<A,B>_F=\sum_{i,j=1}^n a_{ij}b_{ij}$).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A,B
  DOUBLE PRECISION :: FIP

  INTEGER :: I,J

  FIP=0.D0
  DO I=1,N
     DO J=1,N
        FIP=FIP+A(I,J)*B(I,J)
     END DO
  END DO
END FUNCTION FROBENIUSINNERPRODUCT_real

FUNCTION FROBENIUSINNERPRODUCT_complex(A,B,N) RESULT (FIP)
! Function that computes the Frobenius inner product between two square complex matrices (i.e. $<A,B>_F=\sum_{i,j=1}^n a_{ij}\overline{b_{ij}}$).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N,N),INTENT(IN) :: A,B
  DOUBLE COMPLEX :: FIP

  INTEGER :: I,J

  FIP=(0.D0,0.D0)
  DO I=1,N
     DO J=1,N
        FIP=FIP+A(I,J)*CONJG(B(I,J))
     END DO
  END DO
END FUNCTION FROBENIUSINNERPRODUCT_complex

SUBROUTINE NORM_check_norm(CHAR)
  CHARACTER(1),INTENT(IN) :: CHAR
  IF((CHAR /= 'F') .AND. &
       &(CHAR /= 'I') .AND. &
       &(CHAR /= '1') .AND. &
       &(CHAR /= 'M')) THEN
     WRITE(*,*) 'Invalid norm'
     STOP
  END IF
END SUBROUTINE NORM_check_norm

FUNCTION NORM_real(M,N,CHAR) RESULT (NORM)
! Function that computes the one norm, or the Frobenius norm, or the infinity norm of a real matrix.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: M
  CHARACTER(1),INTENT(IN) :: CHAR
  DOUBLE PRECISION :: NORM

  DOUBLE PRECISION,DIMENSION(N) :: WORK
  DOUBLE PRECISION,EXTERNAL :: DLANGE

  CALL NORM_check_norm(CHAR)
  NORM = DLANGE(CHAR,N,N,M,N,WORK)
END FUNCTION NORM_real

FUNCTION NORM_complex(M,N,CHAR) RESULT (NORM)
! Function that computes the one norm, or the Frobenius norm, or the infinity norm of a complex matrix.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N,N),INTENT(IN) :: M
  CHARACTER(1),INTENT(IN) :: CHAR
  DOUBLE PRECISION :: NORM

  DOUBLE PRECISION,DIMENSION(N) :: WORK
  DOUBLE PRECISION,EXTERNAL :: ZLANGE

  CALL NORM_check_norm(CHAR)
  NORM = ZLANGE(CHAR,N,N,M,N,WORK)
END FUNCTION NORM_complex

FUNCTION NORM_symmetric(PM,N,CHAR) RESULT (NORM)
! Function that returns the one norm, or the Frobenius norm, or the infinity norm of a real symmetric matrix, which upper triangular part is stored in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PM
  CHARACTER(1),INTENT(IN) :: CHAR
  DOUBLE PRECISION :: NORM

  DOUBLE PRECISION,DIMENSION(N) :: WORK
  DOUBLE PRECISION,EXTERNAL :: DLANSP

  CALL NORM_check_norm(CHAR)
  NORM = DLANSP(CHAR,'U',N,PM,WORK)
END FUNCTION NORM_symmetric

FUNCTION NORM_hermitian(PM,N,CHAR) RESULT (NORM)
! Function that returns the one norm, or the Frobenius norm, or the infinity norm of a hermitian matrix, which upper triangular part is stored in packed format.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PM
  CHARACTER(1),INTENT(IN) :: CHAR
  DOUBLE PRECISION :: NORM

  DOUBLE PRECISION,DIMENSION(N) :: WORK
  DOUBLE PRECISION,EXTERNAL :: ZLANHP

  CALL NORM_check_norm(CHAR)
  NORM = ZLANHP(CHAR,'U',N,PM,WORK)
END FUNCTION NORM_hermitian

SUBROUTINE EIGENSOLVER_symmetric_prefactorized(PA,PCFB,N,EIG,EIGVEC,INFO)
! Subroutine that computes all the eigenvalues and the eigenvectors of a real generalized symmetric-definite eigenproblem, of the form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric, their upper triangular part being stored in packed format, and B is also positive definite. It is also assumed that the Cholesky factorization of B has previously been computed and stored in packed format.
! Note: it is a simplification of LAPACK's DSPGV subroutine.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PCFB
  DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: EIG
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: EIGVEC
  INTEGER,INTENT(OUT) :: INFO

  INTEGER :: I,NEIG
  DOUBLE PRECISION,DIMENSION(3*N) :: WORK
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: AP,BP

  AP=PA ; BP=PCFB
! Transform problem to standard eigenvalue problem and solve
  CALL DSPGST(1,'U',N,AP,BP,INFO)
  IF (INFO/=0) GO TO 1
  CALL DSPEV('V','U',N,AP,EIG,EIGVEC,N,WORK,INFO)
  IF (INFO/=0) GO TO 2
! Backtransform eigenvectors to the original problem
  NEIG=N
  IF (INFO>0) NEIG=INFO-1
  DO I=1,NEIG
     CALL DTPSV('U','N','Non-unit',N,BP,EIGVEC(1,I),1)
  END DO
  RETURN
! MESSAGES
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DSPGST: the',-INFO,'-th argument had an illegal value'
  END IF
  GO TO 3
2 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DSPEV: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DSPEV: the algorithm failed to converge; ',INFO,'off-diagonal elements &
              &of an intermediate tridiagonal form did not converge to zero'
  END IF
3 WRITE(*,*)'(called from subroutine EIGENSOLVER)'
  RETURN
END SUBROUTINE EIGENSOLVER_symmetric_prefactorized

SUBROUTINE EIGENSOLVER_hermitian_prefactorized(PA,PCFB,N,EIG,EIGVEC,INFO)
! Subroutine that computes all the eigenvalues and the eigenvectors of a complex generalized hermitian-definite eigenproblem, of the form A*x=(lambda)*B*x. Here A and B are assumed to be hermitian, their upper triangular part being stored in packed format, and B is also positive definite. It is also assumed that the Cholesky factorization of B has previously been computed and stored in packed format.
! Note: it is a simplification of LAPACK's ZHPGV subroutine.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PCFB
  DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(N,N),INTENT(OUT) :: EIGVEC
  INTEGER,INTENT(OUT) :: INFO
  INTEGER :: I,NEIG
  DOUBLE PRECISION,DIMENSION(3*N-2) :: RWORK
  DOUBLE COMPLEX,DIMENSION(2*N-1) :: WORK
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: AP,BP

  AP=PA ; BP=PCFB
! Transform problem to standard eigenvalue problem and solve
  CALL ZHPGST(1,'U',N,AP,BP,INFO)
  IF (INFO/=0) GO TO 1
  CALL ZHPEV('V','U',N,AP,EIG,EIGVEC,N,WORK,RWORK,INFO)
  IF (INFO/=0) GO TO 2
! Backtransform eigenvectors to the original problem
  NEIG=N
  IF (INFO>0) NEIG=INFO-1
  DO I=1,NEIG
     CALL ZTPSV('U','N','Non-unit',N,BP,EIGVEC(1,I),1)
  END DO
  RETURN
! MESSAGES
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZHPGST: the',-INFO,'-th argument had an illegal value'
  END IF
  GO TO 3
2 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZHPEV: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZHPEV: the algorithm failed to converge; ',INFO,'off-diagonal elements &
               &of an intermediate tridiagonal form did not converge to zero'
  END IF
3 WRITE(*,*)'(called from subroutine EIGENSOLVER)'
  RETURN
END SUBROUTINE EIGENSOLVER_hermitian_prefactorized

FUNCTION COMMUTATOR_symmetric(PA,PB,PS,N) RESULT (C)
! Function that computes the "commutator" [A,B]=ABS-SBA in a discrete nonorthonormal basis, A and B being two symmetric matrices of size N (only the upper triangular part of the matrices is stored in packed format) and S being the overlap matrix of the basis (stored similarly).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB,PS
  DOUBLE PRECISION,DIMENSION(N,N) :: C

  DOUBLE PRECISION,DIMENSION(N,N) :: A,B,S

  A=UNPACK(PA,N) ; B=UNPACK(PB,N) ; S=UNPACK(PS,N)
  C=MATMUL(MATMUL(A,B),S)-MATMUL(S,MATMUL(B,A))
END FUNCTION COMMUTATOR_symmetric

FUNCTION COMMUTATOR_hermitian(PA,PB,PS,N) RESULT (C)
! Function that computes the "commutator" [A,B]=ABS-SBA in a discrete nonorthonormal basis, A and B being two hermitian matrices of size N (only the upper triangular part of the matrices is stored in packed format) and S being the overlap matrix of the basis (stored similarly).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA,PB,PS
  DOUBLE COMPLEX,DIMENSION(N,N) :: C

  DOUBLE COMPLEX,DIMENSION(N,N) :: A,B,S

  A=UNPACK(PA,N) ; B=UNPACK(PB,N) ; S=UNPACK(PS,N)
  C=MATMUL(MATMUL(A,B),S)-MATMUL(S,MATMUL(B,A))
END FUNCTION COMMUTATOR_hermitian

! input/output routines for matrices

SUBROUTINE READMATRIX_real(MAT,N,LOGUNIT)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,LOGUNIT
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: MAT

  INTEGER :: I
  DOUBLE PRECISION,DIMENSION(N) :: LINE

  DO I=1,N
     READ(LOGUNIT,*)LINE
     MAT(I,:)=LINE
  END DO
END SUBROUTINE READMATRIX_REAL

SUBROUTINE READMATRIX_complex(MAT,N,LOGUNIT_REAL,LOGUNIT_IMAG)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,LOGUNIT_REAL,LOGUNIT_IMAG
  DOUBLE COMPLEX,DIMENSION(N,N),INTENT(OUT) :: MAT

  DOUBLE PRECISION,DIMENSION(N,N) :: R,I

  CALL READMATRIX_real(R,N,LOGUNIT_REAL)
  CALL READMATRIX_real(I,N,LOGUNIT_IMAG)
  MAT=DCMPLX(R,I)
END SUBROUTINE READMATRIX_complex

SUBROUTINE PRINTMATRIX_real(MAT,N,LOGUNIT)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,LOGUNIT
  DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: MAT

  INTEGER :: I

  DO I=1,N
     WRITE(LOGUNIT,*)MAT(I,:)
  END DO
END SUBROUTINE PRINTMATRIX_real

SUBROUTINE PRINTMATRIX_complex(MAT,N,LOGUNIT_REAL,LOGUNIT_IMAG)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,LOGUNIT_REAL,LOGUNIT_IMAG
  DOUBLE COMPLEX,DIMENSION(N,N),INTENT(IN) :: MAT

  INTEGER :: I

  DO I=1,N
     IF (LOGUNIT_REAL==LOGUNIT_IMAG) THEN
        WRITE(LOGUNIT_REAL,*)MAT(I,:)
     ELSE
        WRITE(LOGUNIT_REAL,*)REAL(MAT(I,:))
        WRITE(LOGUNIT_IMAG,*)AIMAG(MAT(I,:))
     END IF
  END DO
END SUBROUTINE PRINTMATRIX_complex

SUBROUTINE PRINTMATRIX_symmetric(PMAT,N,LOGUNIT)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,LOGUNIT
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PMAT

  CALL PRINTMATRIX_real(UNPACK(PMAT,N),N,LOGUNIT)
END SUBROUTINE PRINTMATRIX_symmetric

SUBROUTINE PRINTMATRIX_hermitian(PMAT,N,LOGUNIT_REAL,LOGUNIT_IMAG)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,LOGUNIT_REAL,LOGUNIT_IMAG
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PMAT

  CALL PRINTMATRIX_complex(UNPACK(PMAT,N),N,LOGUNIT_REAL,LOGUNIT_IMAG)
END SUBROUTINE PRINTMATRIX_hermitian
END MODULE matrix_tools

MODULE mathematical_functions
INTERFACE DFACT
  MODULE PROCEDURE DOUBLE_FACTORIAL
END INTERFACE

CONTAINS

FUNCTION FACTORIAL(N) RESULT(FACT)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  INTEGER :: FACT

  INTEGER :: I

  IF (N<0) THEN
     STOP'Function FACTORIAL: the factorial is undefined for negative integers.'
  ELSE 
     FACT=1
     IF (N>1) THEN
        DO I=2,N
           FACT=FACT*I
        END DO
     END IF
  END IF
END FUNCTION FACTORIAL

FUNCTION DOUBLE_FACTORIAL(N) RESULT(DFACT)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  INTEGER :: DFACT

  INTEGER :: I

  IF (N<-1) THEN
     STOP'Function DOUBLE_FACTORIAL: the double factorial is undefined for negative integers lower than -1.'
  ELSE
     DFACT=1
     IF (N>1) THEN
        I=N
        DO WHILE (I>1)
           DFACT=DFACT*I
           I=I-2
        END DO
     END IF
  END IF
END FUNCTION DOUBLE_FACTORIAL
END MODULE

MODULE setup_tools
! name of the setup file
  CHARACTER(100) :: SETUP_FILE
CONTAINS

SUBROUTINE SETUP_FILENAME
! Subroutine that retrieves the name of the setup file

  CALL GETARG(1,SETUP_FILE)
  IF (SETUP_FILE=='') SETUP_FILE='setup'
END SUBROUTINE SETUP_FILENAME

SUBROUTINE LOOKFOR(NUNIT,SUBSTRING,INFO)
! Subroutine that looks for a given text string in an open unit.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NUNIT
  INTEGER,INTENT(OUT) :: INFO
  CHARACTER(*),INTENT(IN) :: SUBSTRING

  CHARACTER(80) :: STRING

  INFO=0
1 READ(100,'(a)',ERR=2,END=2) STRING
  IF (INDEX(STRING,SUBSTRING)==0) GOTO 1
  RETURN
2 INFO=1
END SUBROUTINE LOOKFOR


FUNCTION GETENV(STR) RESULT(RES)
  CHARACTER(*),INTENT(IN) :: STR
  CHARACTER(40) :: RESSTR
  DOUBLE PRECISION :: RES
  CALL GET_ENVIRONMENT_VARIABLE(STR,RESSTR)
  READ(RESSTR,*), RES
END FUNCTION GETENV

END MODULE
