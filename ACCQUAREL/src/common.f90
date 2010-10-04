MODULE debug
  ! This module contains a logical flag to perform various checks in the computation of ESA's THETA function.
  IMPLICIT NONE
  !  LOGICAL,PARAMETER :: THETA_CHECK=.TRUE.
  LOGICAL,PARAMETER :: THETA_CHECK=.FALSE.
END MODULE debug

MODULE metric_relativistic
  ! This module contains pointers assigned to the tensors relative to the metric (as we work in a nonorthogonal basis), which are the overlap matrix S, its inverse S^-{1}, its square root S^{1/2}, the inverse of its square root S^{-1/2} and the matrix of its Cholesky factorization (for the solution of the generalized eigenproblem), all stored in packed form.
  IMPLICIT NONE
  DOUBLE COMPLEX,POINTER,DIMENSION(:) :: PS,PIS,PSRS,PISRS,PCFS
END MODULE metric_relativistic

MODULE metric_nonrelativistic
  ! This module contains pointers assigned to the tensors relative to the metric (as we work in a nonorthogonal basis), which are the overlap matrix S, its inverse S^-{1}, its square root S^{1/2}, the inverse of its square root S^{-1/2} and the matrix of its Cholesky factorization (for the solution of the generalized eigenproblem), all stored in packed form.
  IMPLICIT NONE
  DOUBLE PRECISION,POINTER,DIMENSION(:) :: PS,PIS,PSRS,PISRS,PCFS
END MODULE metric_nonrelativistic

MODULE common_functions
  IMPLICIT NONE
  INTERFACE ENERGY
     MODULE PROCEDURE ENERGY_relativistic,ENERGY_AOCOSDHF,ENERGY_RHF,ENERGY_UHF
  END INTERFACE ENERGY

  INTERFACE ELECTRONIC_ENERGY
     MODULE PROCEDURE ELECTRONIC_ENERGY_relativistic,ELECTRONIC_ENERGY_AOCOSDHF,ELECTRONIC_ENERGY_RHF,ELECTRONIC_ENERGY_UHF
  END INTERFACE ELECTRONIC_ENERGY

CONTAINS

  FUNCTION ENERGY_relativistic(POEFM,PTEFM,PDM,N) RESULT(ETOT)
    ! Function that computes the Dirac-Fock total energy associated to a density matrix (whose upper triangular part is stored in packed form in PDM) of a given molecular system, POEFM and PTEFM containing respectively the core hamiltonian matrix and the two-electron part of the Fock matrix (both stored similarly).
    USE data_parameters
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
    DOUBLE PRECISION :: ETOT

    ETOT=ELECTRONIC_ENERGY_relativistic(POEFM,PTEFM,PDM,N)+INTERNUCLEAR_ENERGY
  END FUNCTION ENERGY_relativistic

  FUNCTION ENERGY_AOCOSDHF(POEFM,PTEFMC,PTEFMO,PDMC,PDMO,N) RESULT(ETOT)
    USE data_parameters
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFMC,PTEFMO,PDMC,PDMO
    DOUBLE PRECISION :: ETOT

    ETOT=ELECTRONIC_ENERGY_AOCOSDHF(POEFM,PTEFMC,PTEFMO,PDMC,PDMO,N)+INTERNUCLEAR_ENERGY
  END FUNCTION ENERGY_AOCOSDHF

  FUNCTION ENERGY_HF(POEFM,PTEFM,PDM,N) RESULT(ENERGY)
    ! Function that computes the Hartree-Fock energy associated to a density matrix (whose upper triangular part is stored in packed form in PDM) of a given system, POEFM and PTEFM containing respectively the (core hamiltonian) matrix and the (two-electron) part of the Fock matrix (both stored similarly).
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
    DOUBLE PRECISION :: ENERGY

    INTEGER :: I,J,IJ
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PEM

    PEM=POEFM+0.5D0*PTEFM

    ENERGY=0.D0
    IJ=0
    DO J=1,N
       DO I=1,J
          IJ=IJ+1
          IF (I.NE.J) THEN
             ENERGY=ENERGY+2.D0*PEM(IJ)*PDM(IJ)
          ELSE
             ENERGY=ENERGY+PEM(IJ)*PDM(IJ)
          END IF
       END DO
    END DO
  END FUNCTION ENERGY_HF

  FUNCTION ENERGY_RHF(POEFM,PTEFM,PDM,N) RESULT(ETOT)
    ! Function that computes the restricted closed-shell Hartree-Fock total energy associated to a density matrix (whose upper triangular part is stored in packed form in PDM) of a given molecular system, POEFM and PTEFM containing respectively the core hamiltonian matrix and the two-electron part of the Fock matrix (both stored similarly).
    USE data_parameters
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
    DOUBLE PRECISION :: ETOT

    ETOT=ELECTRONIC_ENERGY_RHF(POEFM,PTEFM,PDM,N)+INTERNUCLEAR_ENERGY
  END FUNCTION ENERGY_RHF

  FUNCTION ENERGY_UHF(POEFM,PTEFMA,PTEFMB,PDMA,PDMB,N) RESULT(ETOT)
    ! Function that computes the unrestricted open-shell Hartree-Fock total energy associated to the $\alpha$ and $\beta$ spin density matrices (whose upper triangular part is stored in packed form in PDM) of a given molecular system, POEFM and PTEFM containing respectively the core hamiltonian matrix and the two-electron part of the Fock matrix (both stored similarly).
    USE data_parameters
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFMA,PTEFMB,PDMA,PDMB
    DOUBLE PRECISION :: ETOT

    ETOT=ELECTRONIC_ENERGY_UHF(POEFM,PTEFMA,PTEFMB,PDMA,PDMB,N)+INTERNUCLEAR_ENERGY
  END FUNCTION ENERGY_UHF

  FUNCTION ELECTRONIC_ENERGY_relativistic(POEFM,PTEFM,PDM,N) RESULT(ENERGY)
    ! Function that computes the Dirac-Fock electronic energy associated to a density matrix (whose upper triangular part is stored in packed form in PDM) of a given molecular system, POEFM and PTEFM containing respectively the core hamiltonian matrix and the two-electron part of the Fock matrix (both stored similarly).
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
    DOUBLE PRECISION :: ENERGY

    INTEGER :: I,J,IJ
    DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PEM

    PEM=POEFM+0.5D0*PTEFM
    ENERGY=0.D0
    IJ=0
    DO J=1,N
       DO I=1,J
          IJ=IJ+1
          IF (I.NE.J) THEN
             ENERGY=ENERGY+REAL(PEM(IJ)*CONJG(PDM(IJ))+CONJG(PEM(IJ))*PDM(IJ))
          ELSE
             ENERGY=ENERGY+REAL(PEM(IJ)*PDM(IJ))
          END IF
       END DO
    END DO
  END FUNCTION ELECTRONIC_ENERGY_relativistic

  FUNCTION ELECTRONIC_ENERGY_AOCOSDHF(POEFM,PTEFMC,PTEFMO,PDMC,PDMO,N) RESULT(ENERGY)
    ! Function that computes average-of-configuration Dirac-Hartree-Fock electronic energy associated the closed- and open-shell density matrices (whose upper triangular parts are stored in packed form).
    USE data_parameters
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFMC,PTEFMO,PDMC,PDMO
    DOUBLE PRECISION :: ENERGY

    INTEGER :: I,J,IJ
    DOUBLE PRECISION :: a
    DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PEMC,PEMO

    a=REAL(NBOOS*(NBEOS-1))/REAL(NBEOS*(NBOOS-1))
    PEMC=POEFM+0.5D0*PTEFMC+PTEFMO ; PEMO=POEFM+0.5D0*a*PTEFMO
    ENERGY=0.D0
    IJ=0
    DO J=1,N
       DO I=1,J
          IJ=IJ+1
          IF (I.NE.J) THEN
             ENERGY=ENERGY+REAL(PEMC(IJ)*CONJG(PDMC(IJ))+CONJG(PEMC(IJ))*PDMC(IJ)  &
                  &                            +PEMO(IJ)*CONJG(PDMO(IJ))+CONJG(PEMO(IJ))*PDMO(IJ))
          ELSE
             ENERGY=ENERGY+REAL(PEMC(IJ)*PDMC(IJ)+PEMO(IJ)*PDMO(IJ))
          END IF
       END DO
    END DO
  END FUNCTION ELECTRONIC_ENERGY_AOCOSDHF

  FUNCTION ELECTRONIC_ENERGY_RHF(POEFM,PTEFM,PDM,N) RESULT(ENERGY)
    ! Function that computes the restricted closed-shell Hartree-Fock electronic energy associated to a density matrix (whose upper triangular part is stored in packed form in PDM) of a given molecular system, POEFM and PTEFM containing respectively the core hamiltonian matrix and the two-electron part of the Fock matrix (both stored similarly).
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
    DOUBLE PRECISION :: ENERGY

    INTEGER :: I,J,IJ
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PEM

    PEM=2.D0*POEFM+PTEFM
    ENERGY=0.D0
    IJ=0
    DO J=1,N
       DO I=1,J
          IJ=IJ+1
          IF (I.NE.J) THEN
             ENERGY=ENERGY+2.D0*PEM(IJ)*PDM(IJ)
          ELSE
             ENERGY=ENERGY+PEM(IJ)*PDM(IJ)
          END IF
       END DO
    END DO
  END FUNCTION ELECTRONIC_ENERGY_RHF

  FUNCTION ELECTRONIC_ENERGY_UHF(POEFM,PTEFMA,PTEFMB,PDMA,PDMB,N) RESULT(ENERGY)
    ! Function that computes the unrestricted open-shell Hartree-Fock electronic energy associated to the $\alpha$ and $\beta$ spin density matrices (whose upper triangular parts are stored in packed form).
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFMA,PTEFMB,PDMA,PDMB
    DOUBLE PRECISION :: ENERGY

    INTEGER :: I,J,IJ
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PEMA,PEMB

    PEMA=POEFM+0.5D0*PTEFMA ; PEMB=POEFM+0.5D0*PTEFMB
    ENERGY=0.D0
    IJ=0
    DO J=1,N
       DO I=1,J
          IJ=IJ+1
          IF (I.NE.J) THEN
             ENERGY=ENERGY+2.D0*PEMA(IJ)*PDMA(IJ)+PEMB(IJ)*PDMB(IJ)
          ELSE
             ENERGY=ENERGY+PEMA(IJ)*PDMA(IJ)+PEMB(IJ)*PDMB(IJ)
          END IF
       END DO
    END DO
  END FUNCTION ELECTRONIC_ENERGY_UHF

  FUNCTION FREEENERGY_HF(POEFM,PTEFM,PDM,N,TEMPERATURE,FUNC) RESULT(ENRG)
    ! Function that computes the free energy (of a Hartree(-Fock) model with temperature) associated to a density matrix (whose upper triangular part is stored in packed form in PDM) of a given system, POEFM and PTEFM containing respectively the (core hamiltonian) matrix and the (two-electron) part of the Fock matrix (both stored similarly).
    USE matrix_tools ; USE metric_nonrelativistic
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
    DOUBLE PRECISION,INTENT(IN) :: TEMPERATURE
    INTERFACE
       DOUBLE PRECISION FUNCTION FUNC(X)
         IMPLICIT NONE
         DOUBLE PRECISION,INTENT(IN) :: X
       END FUNCTION FUNC
    END INTERFACE
    DOUBLE PRECISION :: ENRG

    INTEGER :: I,J,IJ,INFO
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PTMP
    DOUBLE PRECISION,DIMENSION(N) :: W
    DOUBLE PRECISION,DIMENSION(3*N) :: WORK
    DOUBLE PRECISION,DIMENSION(N,N) :: Z,TMP

    PTMP=PDM
    CALL DSPEV('V','U',N,PTMP,W,Z,N,WORK,INFO)
    IF (INFO/=0) GO TO 1
    TMP=0.D0
    DO I=1,N
       TMP(I,I)=FUNC(W(I))
    END DO
    PTMP=PACK(MATMUL(Z,MATMUL(TMP,TRANSPOSE(Z))),N)

    ENRG=ENERGY_HF(POEFM,PTEFM,PDM,N)+TEMPERATURE*TRACEOFPRODUCT(PTMP,PS,N)
    RETURN
1   IF (INFO<0) THEN
       WRITE(*,*)'Subroutine DSPEV: the',-INFO,'-th argument had an illegal value'
    ELSE
       WRITE(*,*)'Subroutine DSPEV: the algorithm failed to converge; ',INFO,'off-diagonal elements &
            &of an intermediate tridiagonal form did not converge to zero'
    END IF
    WRITE(*,*)'(called from function FREEENERGY)'
    STOP
  END FUNCTION FREEENERGY_HF
END MODULE common_functions
