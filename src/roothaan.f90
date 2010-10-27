SUBROUTINE ROOTHAAN_relativistic(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Roothaan's algorithm (closed-shell Dirac-Hartree-Fock formalism).
! Reference: C. C. J. Roothaan, New developments in molecular orbital theory, Rev. Modern Phys., 23(2), 69-89, 1951.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_relativistic ; USE scf_tools
  USE graphics_tools
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(OUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(OUT) :: EIGVEC
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,LOON,INFO,I
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PDM1
  LOGICAL :: NUMCONV

! INITIALIZATION AND PRELIMINARIES
  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))

  ITER=0
  PDM=(0.D0,0.D0)
  PTEFM=(0.D0,0.D0)
  ETOT1=0.D0
  OPEN(16,FILE='plots/rootenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/rootcrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/rootcrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,'(a)')' '
  WRITE(*,'(a,i3)')'# ITER = ',ITER

! Assembly and diagonalization of the Fock matrix
  PFM=POEFM+PTEFM
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix according to the aufbau principle
  CALL CHECKORB(EIG,NBAST,LOON)
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,LOON,LOON+NBE-1)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT
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
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ROOTHAAN: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ROOTHAAN: no convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  CALL EXPORT_DENSITY(PDM,PHI,NBAST,-3.D0,3.D0,50,'density','cube')
  GO TO 5
4 WRITE(*,*)'(called from subroutine ROOTHAAN)'
5 DEALLOCATE(PDM,PDM1,PTEFM,PFM)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
END SUBROUTINE ROOTHAAN_relativistic

SUBROUTINE ROOTHAAN_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Roothaan's algorithm (restricted closed-shell Hartree-Fock formalism).
! Reference: C. C. J. Roothaan, New developments in molecular orbital theory, Rev. Modern Phys., 23(2), 69-89, 1951.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools
  USE graphics_tools
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(OUT) :: EIG
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(OUT) :: EIGVEC
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,INFO,I
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PDM1
  LOGICAL :: NUMCONV

! INITIALIZATION AND PRELIMINARIES
  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))

  ITER=0
  PDM=0.D0
  PTEFM=0.D0
  ETOT1=0.D0
  OPEN(16,FILE='plots/rootenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/rootcrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/rootcrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix
  PFM=POEFM+PTEFM
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix according to the aufbau principle
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE/2)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT
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
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ROOTHAAN: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ROOTHAAN: no convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
4 WRITE(*,*)'(called from subroutine ROOTHAAN)'
5 DEALLOCATE(PDM,PDM1,PTEFM,PFM)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
END SUBROUTINE ROOTHAAN_RHF

SUBROUTINE ROOTHAAN_UHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Roothaan's algorithm (unrestricted open-shell Hartree-Fock formalism).
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
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PEMS,PFMA,PFMB,PDMA,PDMB,PTDM,PSDM,PTDM1
  LOGICAL :: NUMCONV

! INITIALIZATION AND PRELIMINARIES
  ALLOCATE(PDMA(1:NBAST*(NBAST+1)/2),PDMB(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTDM(1:NBAST*(NBAST+1)/2),PSDM(1:NBAST*(NBAST+1)/2),PTDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PEMS(1:NBAST*(NBAST+1)/2),PFMA(1:NBAST*(NBAST+1)/2),PFMB(1:NBAST*(NBAST+1)/2))

  ITER=0
  PTDM=0.D0
  PFMA=POEFM ; PFMB=POEFM
  ETOT1=0.D0

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Diagonalization of the Fock matrix for $\alpha$ spin orbitals
  CALL EIGENSOLVER(PFMA,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix for $\alpha$ spin orbitals according to the aufbau principle
  CALL FORMDM(PDMA,EIGVEC,NBAST,1,NBEA)
! Diagonalization of the Fock matrix for $\beta$ spin orbitals
  CALL EIGENSOLVER(PFMB,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix for $\beta$ spin orbitals according to the aufbau principle
  CALL FORMDM(PDMB,EIGVEC,NBAST,1,NBEB)
! Assembly of the total and spin density matrices
  PTDM1=PTDM
  PTDM=PDMA+PDMB ; PSDM=PDMA-PDMB
! Computation of the energy associated to the density matrices
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PTDM)
  CALL BUILDEXCHANGE(PEMS,NBAST,PHI,PSDM)
  ETOT=ENERGY(POEFM,PTEFM,PEMS,PTDM,PSDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT
! Assembly of the Fock matrices
  PFMA=POEFM+0.5D0*(PTEFM-PEMS) ; PFMB=POEFM+0.5D0*(PTEFM+PEMS)
! Numerical convergence check
  CALL CHECKNUMCONV(PDMA,PDMB,PTDM1,PFMA,PFMB,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV)
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
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ROOTHAAN: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ROOTHAAN: no convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
4 WRITE(*,*)'(called from subroutine ROOTHAAN)'
5 DEALLOCATE(PDMA,PDMB,PTDM,PSDM,PTDM1,PTEFM,PEMS,PFMA,PFMB)
END SUBROUTINE ROOTHAAN_UHF

! NON WORKING PART!!!!!!!

SUBROUTINE ROOTHAAN_AOCOSDHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR)
! Roothaan's algorithm (average-of-configuration open-shell Dirac-Hartree-Fock formalism).
! References: C. C. J. Roothaan, Self-consistent field theory for open shells of electronic systems, Rev. Modern Phys., 32(2), 179-185, 1960.
! ref average-of-configuration?
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_relativistic ; USE scf_tools
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(OUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(OUT) :: EIGVEC
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR

  INTEGER :: ITER,LOON,INFO,I
  DOUBLE PRECISION :: f,a,alpha
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: PTEFMC,PTEFMO,PFMC,PFMO,PDMC,PDMO,PDMC1,PDMO1
  LOGICAL :: NUMCONV

! INITIALIZATION AND PRELIMINARIES
  ALLOCATE(PDMC(1:NBAST*(NBAST+1)/2),PDMO(1:NBAST*(NBAST+1)/2),PDMC1(1:NBAST*(NBAST+1)/2),PDMO1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFMC(1:NBAST*(NBAST+1)/2),PTEFMO(1:NBAST*(NBAST+1)/2),PFMC(1:NBAST*(NBAST+1)/2),PFMO(1:NBAST*(NBAST+1)/2))

  f=REAL(NBEOS)/REAL(NBOOS) ; a=NBOOS*(NBEOS-1.D0)/(NBEOS*(NBOOS-1.D0)) ; alpha=(1.D0-a)/(1.D0-f)

  ITER=0
  PDMC=(0.D0,0.D0) ; PDMO=(0.D0,0.D0)
  PTEFMC=(0.D0,0.D0) ; PTEFMO=(0.D0,0.D0)
  ETOT1=0.D0

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,'(a)')' '
  WRITE(*,'(a,i3)')'# ITER = ',ITER

! Assembly and diagonalization of the Fock matrix for closed-shell orbitals
  PFMC=POEFM+PTEFMC+PTEFMO+alpha*ABC_CBA(PS,PDMO,PTEFMO,NBAST)
  CALL EIGENSOLVER(PFMC,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix for closed-shell orbitals according to the aufbau principle
  CALL CHECKORB(EIG,NBAST,LOON)
  PDMC1=PDMC
  CALL FORMDM(PDMC,EIGVEC,NBAST,LOON,LOON+NBECS-1)
! Assembly and diagonalization of the Fock matrix for open-shell orbitals
  PFMO=POEFM+PTEFMC+a*PTEFMO+alpha*ABC_CBA(PS,PDMC,PTEFMO,NBAST)
  CALL EIGENSOLVER(PFMO,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
! Assembly of the density matrix for open-shell orbitals according to the aufbau principle
  CALL CHECKORB(EIG,NBAST,LOON)
  PDMO1=PDMO
  CALL FORMDM(PDMO,EIGVEC,NBAST,LOON+NBECS,LOON+NBECS+NBOOS-1)
  PDMO=f*PDMO
! Computation of the energy associated to the closed- and open-shell density matrices
  CALL BUILDTEFM(PTEFMC,NBAST,PHI,PDMC)
  CALL BUILDTEFM(PTEFMO,NBAST,PHI,PDMO)
  ETOT=ENERGY_AOCOSDHF(POEFM,PTEFMC,PTEFMO,PDMC,PDMO,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT
! Numerical convergence check
  PFMC=POEFM+PTEFMC+PTEFMO+alpha*ABC_CBA(PS,PDMO,PTEFMO,NBAST)
  PFMO=POEFM+PTEFMC+a*PTEFMO+alpha*ABC_CBA(PS,PDMC,PTEFMO,NBAST)
  CALL CHECKNUMCONV(PDMC,PDMO,PDMC1,PDMO1,PFMC,PFMO,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV)
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
2 WRITE(*,*)'Subroutine ROOTHAAN: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
3 WRITE(*,*)'Subroutine ROOTHAAN: no convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 5
4 WRITE(*,*)'(called from subroutine ROOTHAAN)'
5 DEALLOCATE(PDMC,PDMO,PDMC1,PDMO1,PTEFMC,PTEFMO,PFMC,PFMO)
END SUBROUTINE ROOTHAAN_AOCOSDHF
