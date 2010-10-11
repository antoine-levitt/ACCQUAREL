SUBROUTINE PENALIZEDROOTHAAN(EIG,EIGVEC,NBAST,NBAS,POEFM,PHI,TRSHLD,MAXITR)
! Penalized version (using total angular momentum operators) of Roothaan's algorithm (relativistic case)
! Reference: C. C. J. Roothaan, New developments in molecular orbital theory, Rev. Modern Phys., 23(2), 69-89, 1951.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_relativistic ; USE scf_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(OUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(OUT) :: EIGVEC
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR

  INTEGER :: ITER,LOON,INFO,I,J,IJ
  DOUBLE PRECISION :: ETOT,ETOT1,ETOTP,ETOTP1
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PDM1,PJ
  LOGICAL :: NUMCONV,NUMCONVP
  DOUBLE COMPLEX :: EPS
  
! INITIALIZATION AND PRELIMINARIES
  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2),PJ(1:NBAST*(NBAST+1)/2))

! penalization coefficient
  EPS=(1.D0,0.D0)
! assembly of the matrix of one the total angular momentum operators
  CALL BUILDTAMCM(PJ,PHI,NBAST,NBAS,3)

  ITER=0
  PDM=(0.D0,0.D0)
  PTEFM=(0.D0,0.D0)
  ETOT1=0.D0 ; ETOTP1=0.D0
  OPEN(16,FILE='plots/prootenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/prootcrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/prootcrit2.txt',STATUS='unknown',ACTION='write')

! PENALIZATION LOOP
1 CONTINUE
  WRITE(*,*)'penalization coefficient =',EPS
! SCF LOOP
2 CONTINUE
  ITER=ITER+1

  WRITE(*,'(a)')' '
  WRITE(*,'(a,i3)')'# ITER = ',ITER

! Assembly and diagonalization of the Fock matrix
  PFM=POEFM+PTEFM-EPS*PJ
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 5
! Assembly of the density matrix according to the aufbau principle
  CALL CHECKORB(EIG,NBAST,LOON)
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,LOON,LOON+NBE-1)
! Computation of the energy (wrt. the penalized hamiltonian) associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOTP=ENERGY(POEFM-EPS*PJ,PTEFM,PDM,NBAST)
  WRITE(*,*)'E_\eps(D_n)=',ETOTP
  WRITE(53,'(e22.14)')ETOTP
! Numerical convergence check for the SCF loop
  CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM-EPS*PJ,NBAST,ETOTP,ETOTP1,TRSHLD,NUMCONV)
  IF (NUMCONV) THEN
! Convergence reached for the SCF loop
! Numerical convergence check for the penalization loop
     ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
     WRITE(*,*)'E(D_n)=',ETOT
     WRITE(16,'(e22.14)')ETOT
     CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONVP)
     IF (NUMCONVP) THEN
! convergence reached for the penalization loop, exit
        GO TO 3
     ELSE
! the penalization coefficient is decreased, new SCF loop
        IF (ITER==MAXITR) GO TO 4
        ETOT1=ETOT
        EPS=0.9D0*EPS
        GO TO 1
     END IF
  ELSE IF (ITER==MAXITR) THEN
! Maximum number of iterations reached without convergence in the SCF loop
     GO TO 4
  ELSE
! Convergence not reached for the SCF loop, increment
     ETOTP1=ETOTP
     GO TO 2
  END IF

! MESSAGES
3 WRITE(*,*)'Subroutine PENROOTHAAN: convergence after',ITER,'iterations.'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
4 WRITE(*,*)'Subroutine PENROOTHAAN: no convergence after',ITER,'iterations.'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
5 WRITE(*,*)'(called from subroutine PENALIZEDROOTHAAN)'
6 DEALLOCATE(PDM,PDM1,PTEFM,PFM,PJ)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
  RETURN
END SUBROUTINE
