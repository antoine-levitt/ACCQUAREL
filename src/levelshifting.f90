SUBROUTINE LEVELSHIFTING_relativistic(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Level-shifting algorithm (relativistic case)
! Reference: V. R. Saunders and I. H. Hillier, A "level-shifting" method for converging closed shell Hartree-Fock wave functions, Int. J. Quantum Chem., 7(4), 699-705, 1973.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_relativistic ; USE scf_tools ; USE setup_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(OUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(OUT) :: EIGVEC
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,LOON,INFO,I
  DOUBLE PRECISION :: SHIFT,ETOT,ETOT1
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PDM1
  LOGICAL :: NUMCONV

! INITIALIZATION AND PRELIMINARIES
! Reading of the shift parameter
  OPEN(100,FILE=SETUP_FILE,STATUS='OLD',ACTION='READ')
  CALL LOOKFOR(100,'LEVEL-SHIFTING ALGORITHM PARAMETERS',INFO)
  READ(100,'(/,f16.8)')SHIFT
  CLOSE(100)
  WRITE(*,*)'Shift parameter value =',SHIFT

  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))

  ITER=0
  PDM=(0.D0,0.D0)
  PTEFM=(0.D0,0.D0)
  ETOT1=0.D0
  OPEN(16,FILE='plots/shftenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/shftcrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/shftcrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,'(a)')' '
  WRITE(*,'(a,i3)')'# ITER = ',ITER

! Assembly and diagonalization of the Fock matrix
  PFM=POEFM+PTEFM-SHIFT*ABA(PS,PDM,NBAST)
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix according to the aufbau principle
  CALL CHECKORB(EIG,NBAST,LOON)
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,LOON,LOON+NBE-1)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT,'Ha'
! Numerical convergence check
  CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
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
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine LEVELSHIFTING: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine LEVELSHIFTING: no convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
4 WRITE(*,*)'(called from subroutine LEVELSHIFTING)'
5 DEALLOCATE(PDM,PDM1,PTEFM,PFM)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
END SUBROUTINE LEVELSHIFTING_relativistic

SUBROUTINE LEVELSHIFTING_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Level-shifting algorithm (restricted closed-shell Hartree-Fock formalism)
! Reference: V. R. Saunders and I. H. Hillier, A "level-shifting" method for converging closed shell Hartree-Fock wave functions, Int. J. Quantum Chem., 7(4), 699-705, 1973.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools ; USE setup_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(OUT) :: EIG
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(OUT) :: EIGVEC
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,INFO,I
  DOUBLE PRECISION :: SHIFT,ETOT,ETOT1
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PDM1
  LOGICAL :: NUMCONV

! INITIALIZATION AND PRELIMINARIES
! Reading of the shift parameter
  OPEN(100,FILE=SETUP_FILE,STATUS='OLD',ACTION='READ')
  CALL LOOKFOR(100,'LEVEL-SHIFTING ALGORITHM PARAMETERS',INFO)
  READ(100,'(/,f16.8)')SHIFT
  CLOSE(100)
  WRITE(*,*)'Shift parameter value =',SHIFT

  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))

  ITER=0
  PDM=0.D0
  PTEFM=0.D0
  ETOT1=0.D0
  OPEN(16,FILE='plots/shftenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/shftcrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/shftcrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix
  PFM=POEFM+PTEFM-SHIFT*ABA(PS,PDM,NBAST)
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix according to the aufbau principle
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE/2)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT,'Ha'
! Numerical convergence check
  CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
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
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine LEVELSHIFTING: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine LEVELSHIFTING: no convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
4 WRITE(*,*)'(called from subroutine LEVELSHIFTING)'
5 DEALLOCATE(PDM,PDM1,PTEFM,PFM)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
END SUBROUTINE LEVELSHIFTING_RHF

SUBROUTINE LEVELSHIFTING_RGHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Level-shifting algorithm (real general Hartree-Fock formalism)
! Reference: V. R. Saunders and I. H. Hillier, A "level-shifting" method for converging closed shell Hartree-Fock wave functions, Int. J. Quantum Chem., 7(4), 699-705, 1973.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools ; USE setup_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(OUT) :: EIG
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(OUT) :: EIGVEC
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,INFO,I
  DOUBLE PRECISION :: SHIFT,ETOT,ETOT1
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PDM1
  LOGICAL :: NUMCONV

! INITIALIZATION AND PRELIMINARIES
! Reading of the shift parameter
  OPEN(100,FILE=SETUP_FILE,STATUS='OLD',ACTION='READ')
  CALL LOOKFOR(100,'LEVEL-SHIFTING ALGORITHM PARAMETERS',INFO)
  READ(100,'(/,f16.8)')SHIFT
  CLOSE(100)
  WRITE(*,*)'Shift parameter value =',SHIFT

  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))

  ITER=0
  PDM=0.D0
  PTEFM=0.D0
  ETOT1=0.D0
  OPEN(16,FILE='plots/shftenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/shftcrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/shftcrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix
  PFM=POEFM+PTEFM-SHIFT*ABA(PS,PDM,NBAST)
  IF(.NOT.(RESUME .AND. ITER == 1)) THEN
     CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
     IF (INFO/=0) GO TO 4
  END IF
! Assembly of the density matrix according to the aufbau principle
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM_RGHF(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY_RGHF(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT,'Ha'
! Numerical convergence check
  CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
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
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine LEVELSHIFTING: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine LEVELSHIFTING: no convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
4 WRITE(*,*)'(called from subroutine LEVELSHIFTING)'
5 DEALLOCATE(PDM,PDM1,PTEFM,PFM)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
END SUBROUTINE LEVELSHIFTING_RGHF
