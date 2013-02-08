SUBROUTINE ODA_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Optimal Damping Algorithm (restricted closed-shell Hartree-Fock formalism)
! Reference: E. Cancès and C. Le Bris, Can we outperform the DIIS approach for electronic structure calculations?, Internat. J. Quantum Chem., 79(2), 82-90, 2000.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools
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
  DOUBLE PRECISION :: ETOT,ETOT1,ETTOT
  DOUBLE PRECISION :: ALPHA,BETA,LAMBDA
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PTTEFM,PFM,PDM,PDM1,PTDM,PDMDIF
  LOGICAL :: NUMCONV

! INITIALIZATIONS AND PRELIMINARIES
  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2),PDMDIF(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PTTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))

  ITER=0
  PDM=0.D0 ; PTDM=0.D0
  ETOT1=0.D0
  OPEN(16,FILE='plots/odaenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/odacrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/odacrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix associated to the pseudo-density matrix
  CALL BUILDTEFM(PTTEFM,NBAST,PHI,PTDM)
  PFM=POEFM+PTTEFM
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 5
! Assembly of the density matrix according to the aufbau principle
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE/2)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT
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
! Optimization step for the assembly of the pseudo-density matrix
     PDMDIF=PDM-PTDM
     BETA=TRACEOFPRODUCT(POEFM+PTTEFM,PDMDIF,NBAST)
     IF (BETA>0.D0) THEN
        WRITE(*,*)'Warning: internal computation error (beta>0).'
        GO TO 4
     ELSE
        CALL BUILDTEFM(PTTEFM,NBAST,PHI,PDMDIF)
        ALPHA=TRACEOFPRODUCT(PTTEFM,PDMDIF,NBAST)
        WRITE(*,*)'alpha =',ALPHA,'beta =',BETA
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
          WRITE(*,*)'Warning: internal computation error (alpha=0).'
          GO TO 4
        END IF
     END IF
     ETOT1=ETOT
     CALL BUILDTEFM(PTTEFM,NBAST,PHI,PTDM)
! Energy associated to the pseudo density matrix
     ETTOT=ENERGY(POEFM,PTTEFM,PTDM,NBAST)
     WRITE(*,*)'E(tilde{D}_n)=',ETTOT
     GO TO 1
  END IF
! MESSAGES
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ODA: convergence after',ITER,'iteration(s).'
  GO TO 6
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ODA: no convergence after',ITER,'iteration(s).'
  GO TO 6
4 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ODA: computation stopped after',ITER,'iteration(s).'
  GO TO 6
5 WRITE(*,*)'(called from subroutine ODA)'
  GO TO 7
6 OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
7 DEALLOCATE(PDM,PDM1,PTDM,PDMDIF,PTEFM,PTTEFM,PFM)
  CLOSE(16) ; CLOSE(17) ;CLOSE(18)
END SUBROUTINE ODA_RHF

SUBROUTINE ODA_RGHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Optimal Damping Algorithm (real general Hartree-Fock formalism)
! Reference: E. Cancès and C. Le Bris, Can we outperform the DIIS approach for electronic structure calculations?, Internat. J. Quantum Chem., 79(2),82-90,  2000.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools
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
  DOUBLE PRECISION :: ETOT,ETOT1,ETTOT
  DOUBLE PRECISION :: ALPHA,BETA,LAMBDA
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PTTEFM,PFM,PDM,PDM1,PTDM,PDMDIF
  LOGICAL :: NUMCONV

! INITIALIZATIONS AND PRELIMINARIES
  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2),PDMDIF(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PTTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))

  ITER=0
  PDM=0.D0 ; PTDM=0.D0
  ETOT1=0.D0
  OPEN(16,FILE='plots/odaenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/odacrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/odacrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix associated to the pseudo-density matrix
  CALL BUILDTEFM_RGHF(PTTEFM,NBAST,PHI,PTDM)
  PFM=POEFM+PTTEFM
  IF(.NOT.(RESUME .AND. ITER == 1)) THEN
     CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
     IF (INFO/=0) GO TO 5
  END IF
  
! Assembly of the density matrix according to the aufbau principle
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM_RGHF(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY_RGHF(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT
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
! Optimization step for the assembly of the pseudo-density matrix
     PDMDIF=PDM-PTDM
     BETA=TRACEOFPRODUCT(POEFM+PTTEFM,PDMDIF,NBAST)
     IF (BETA>0.D0) THEN
        WRITE(*,*)'Warning: internal computation error (beta>0).'
        write(*,*) BETA
        BETA=0
     END IF
     ! IF (BETA>0.D0) THEN
     !    WRITE(*,*)'Warning: internal computation error (beta>0).'
     !    write(*,*) BETA
     !    GO TO 4
     ! ELSE
        CALL BUILDTEFM_RGHF(PTTEFM,NBAST,PHI,PDMDIF)
        ALPHA=TRACEOFPRODUCT(PTTEFM,PDMDIF,NBAST)
        WRITE(*,*)'alpha =',ALPHA,'beta =',BETA
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
          WRITE(*,*)'Warning: internal computation error (alpha=0).'
          GO TO 4
        END IF
     ! END IF
     ETOT1=ETOT
     CALL BUILDTEFM_RGHF(PTTEFM,NBAST,PHI,PTDM)
! Energy associated to the pseudo density matrix
     ETTOT=ENERGY_RGHF(POEFM,PTTEFM,PTDM,NBAST)
     WRITE(*,*)'E(tilde{D}_n)=',ETTOT
     GO TO 1
  END IF
! MESSAGES
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ODA: convergence after',ITER,'iteration(s).'
  GO TO 6
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ODA: no convergence after',ITER,'iteration(s).'
  GO TO 6
4 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ODA: computation stopped after',ITER,'iteration(s).'
  GO TO 6
5 WRITE(*,*)'(called from subroutine ODA)'
  GO TO 7
6 OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
7 DEALLOCATE(PDM,PDM1,PTDM,PDMDIF,PTEFM,PTTEFM,PFM)
  CLOSE(16) ; CLOSE(17) ;CLOSE(18)
END SUBROUTINE ODA_RGHF
