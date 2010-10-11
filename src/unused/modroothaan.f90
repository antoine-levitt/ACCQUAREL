SUBROUTINE MODIFIEDROOTHAAN(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR)
! modified Roothaan algorithm (made in ACCQUAREL@CEREMADE)
  USE data_parameters ; USE basis_parameters ; USE common_functions ; USE matrices
  USE matrix_tools ; USE metric_relativistic ; USE scf_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(OUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(OUT) :: EIGVEC
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR

  INTEGER :: ITER,MXSUBITR,LOON,INFO,I
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PDM1,PPROJM
  LOGICAL :: NUMCONV

  INTERFACE RODA
    SUBROUTINE ODA_relativistic(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME,PPROJM)
      USE basis_parameters
      INTEGER,INTENT(IN) :: NBAST
      DOUBLE PRECISION,DIMENSION(NBAST),INTENT(INOUT) :: EIG
      DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(INOUT) :: EIGVEC
      DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
      TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
      DOUBLE PRECISION,INTENT(IN) :: TRSHLD
      INTEGER,INTENT(IN) :: MAXITR
      LOGICAL,INTENT(IN) :: RESUME
      DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PPROJM
  END SUBROUTINE
END INTERFACE

! INITIALIZATIONS AND PRELIMINARIES
  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2),PPROJM(1:NBAST*(NBAST+1)/2))

  MXSUBITR=20
  ITER=0
  PDM=(0.D0,0.D0)
  PTEFM=(0.D0,0.D0)
  ETOT1=0.D0

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix
  PFM=POEFM+PTEFM
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Computation of the projector
  CALL CHECKORB(EIG,NBAST,LOON)
  CALL FORMPROJ(PPROJM,EIGVEC,NBAST,LOON)
! Computation of the density matrix using the Optimal Damping Algorithm (adapted to a projected Dirac-Fock functional)
  CALL RODA(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MXSUBITR,.FALSE.,PPROJM)
!  CALL PROJECTEDROOTHAAN(EIG,EIGVEC,NBAST,POEFM,PHI,BILIST,BITYPE,BIVALUES,BINMBR,TRSHLD,MXSUBITR,PPROJM)
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)' '
  WRITE(*,*)'Total energy =',ETOT
! Numerical convergence check
  CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV)
  IF (NUMCONV) THEN
! Convergence reached
     GO TO 2
  ELSE IF (ITER==MAXITR) THEN
! Maximum number of iterations reached without convergence
     GO TO 3
  ELSE
! Convergence not reached, increment
     ETOT1=ETOT
     GO TO 1
  END IF

! MESSAGES
2 WRITE(*,*)'Subroutine MODIFIEDROOTHAAN: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
3 WRITE(*,*)'Subroutine MODIFIEDROOTHAAN: no convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
4 WRITE(*,*)'(called from subroutine MODIFIEDROOTHAAN)'
5 DEALLOCATE(PDM,PDM1,PTEFM,PFM,PPROJM)
  RETURN
END SUBROUTINE

SUBROUTINE ODA_relativistic(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME,PPROJM)
! Optimal Damping Algorithm (relativistic case)
! Reference: E. Cancès and C. Le Bris, Can we outperform the DIIS approach for electronic structure calculations?, Internat. J. Quantum Chem., 79(2), 2000.
! Note: the algorithm coded in the present routine has been adapted to the relativistic setting using a projected Dirac-based hamiltonian (hence the appearance of a projector matrix in the routine variables).
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_relativistic ; USE scf_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(INOUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(INOUT) :: EIGVEC
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PPROJM

  INTEGER :: ITER,LOON,INFO,I
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE PRECISION :: ALPHA,BETA,LAMBDA
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: PTEFM,PTTEFM,PPOEFM,PPTEFM,PPFM,PDM,PDM1,PTDM,PDMDIF
  LOGICAL :: NUMCONV

  WRITE(*,*)' '
  WRITE(*,*)'SCF subcycles using ODA'

! INITIALIZATIONS AND PRELIMINARIES
  ITER=0
  ETOT1=0.D0

  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2),PDMDIF(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PTTEFM(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PPOEFM(1:NBAST*(NBAST+1)/2),PPTEFM(1:NBAST*(NBAST+1)/2),PPFM(1:NBAST*(NBAST+1)/2))
! initialisation de l'algorithme en projetant la matrice densite obtenue a l'iteration precedente (D_n) sur l'orthogonal du noyau du projecteur "courant" P_{n+1}=(projecteur negatif associe a F(D_n))
! Note : cela necessite de s'assurer que le rang de la matrice projetee reste bien le meme (semble OK si ||P_n-P_{n+1}||>1) + re-orthonormalisation (gram-schmitt)
! Solution actuelle : on fabrique cette matrice a partir des NBE premiers vecteurs propres (dans le gap)
  CALL CHECKORB(EIG,NBAST,LOON)
  CALL FORMDM(PDM,EIGVEC,NBAST,LOON,LOON+NBE-1)
  PTDM=PDM
  PPOEFM=ABCBA(PS,PPROJM,POEFM,NBAST)

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'### SUBITER =',ITER

! Assembly and diagonalization of the projected Fock matrix associated to the pseudo-density matrix
  CALL BUILDTEFM(PTTEFM,NBAST,PHI,PTDM)
  PPFM=PPOEFM+ABCBA(PS,PPROJM,PTTEFM,NBAST)
  CALL EIGENSOLVER(PPFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix according to the aufbau principle
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  PPTEFM=ABCBA(PS,PPROJM,PTEFM,NBAST)
  ETOT=ENERGY(PPOEFM,PPTEFM,PTDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT
! Numerical convergence check
  CALL CHECKNUMCONV(PDM,PDM1,PPOEFM+PPTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV)
  IF (NUMCONV) THEN
! Convergence reached
     GO TO 2
  ELSE IF (ITER==MAXITR) THEN
! Maximum number of iterations reached without convergence
     GO TO 3
  ELSE
! Convergence not reached, increment
! Optimization step for the assembly of the pseudo-density matrix
     PDMDIF=PDM-PTDM
     BETA=REAL(TRACEOFPRODUCT(ABCBA(PS,PPROJM,POEFM+PTTEFM,NBAST),PDMDIF,NBAST))
     IF (BETA>0.D0) THEN
        STOP'Subroutine ODA: internal computation error (beta>0)'
     ELSE
        CALL BUILDTEFM(PTTEFM,NBAST,PHI,PDMDIF)
        ALPHA=REAL(TRACEOFPRODUCT(ABCBA(PS,PPROJM,PTTEFM,NBAST),PDMDIF,NBAST))
        IF (ALPHA>0.D0) THEN
           LAMBDA=-BETA/ALPHA
           IF (LAMBDA<1.D0) THEN
              WRITE(*,*)'lambda=',LAMBDA
              PTDM=PTDM+LAMBDA*PDMDIF
           ELSE
              WRITE(*,*)'lambda=1.'
              PTDM=PDM
           END IF
        ELSE IF (ALPHA<0.D0) THEN
           WRITE(*,*)'lambda=1.'
           PTDM=PDM
        ELSE
          STOP'Subroutine ODA: internal computation error (alpha=0).'
        END IF
     END IF
     ETOT1=ETOT
     GO TO 1
  END IF
! MESSAGES
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ODA: convergence after',ITER,'iteration(s).'
  GO TO 5
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ODA: no convergence after',ITER,'iteration(s).'
  GO TO 5
4 WRITE(*,*)'(called from subroutine ODA)'
5 DEALLOCATE(PDM,PDM1,PTDM,PDMDIF,PTEFM,PTTEFM,PPOEFM,PPTEFM,PPFM)
END SUBROUTINE
