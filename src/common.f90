MODULE debug
! This module contains a logical flag to perform various checks in the computation of ESA's THETA function.
!  LOGICAL,PARAMETER :: THETA_CHECK=.TRUE.
  LOGICAL :: THETA_CHECK=.TRUE.
END MODULE

MODULE metric_relativistic
! This module contains pointers assigned to the tensors relative to the metric (as we work in a nonorthogonal basis), which are the overlap matrix S, its inverse S^-{1}, its square root S^{1/2}, the inverse of its square root S^{-1/2} and the matrix of its Cholesky factorization (for the solution of the generalized eigenproblem), all stored in packed form.  
   DOUBLE COMPLEX,POINTER,DIMENSION(:) :: PS,PIS,PSRS,PISRS,PCFS
END MODULE

MODULE metric_nonrelativistic
! This module contains pointers assigned to the tensors relative to the metric (as we work in a nonorthogonal basis), which are the overlap matrix S, its inverse S^-{1}, its square root S^{1/2}, the inverse of its square root S^{-1/2} and the matrix of its Cholesky factorization (for the solution of the generalized eigenproblem), all stored in packed form.  
   DOUBLE PRECISION,POINTER,DIMENSION(:) :: PS,PIS,PSRS,PISRS,PCFS
END MODULE

MODULE common_functions
INTERFACE ENERGY
  MODULE PROCEDURE ENERGY_relativistic,ENERGY_AOCOSDHF,ENERGY_RHF,ENERGY_UHF
END INTERFACE

INTERFACE ELECTRONIC_ENERGY
  MODULE PROCEDURE ELECTRONIC_ENERGY_relativistic,ELECTRONIC_ENERGY_AOCOSDHF,ELECTRONIC_ENERGY_RHF,ELECTRONIC_ENERGY_UHF
END INTERFACE

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

FUNCTION ENERGY_RHF(POEFM,PTEFM,PDM,N) RESULT(ETOT)
! Function that computes the restricted closed-shell Hartree-Fock total energy associated to a density matrix (whose upper triangular parts are stored in packed form) of a given molecular system, POEFM and PTEFM containing respectively the core hamiltonian matrix and the two-electron part of the Fock matrix (both stored similarly).
  USE data_parameters
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
  DOUBLE PRECISION :: ETOT

  ETOT=ELECTRONIC_ENERGY_RHF(POEFM,PTEFM,PDM,N)+INTERNUCLEAR_ENERGY
END FUNCTION ENERGY_RHF

FUNCTION ENERGY_UHF(POEFM,PTEFM,PEMS,PTDM,PSDM,N) RESULT(ETOT)
! Function that computes the unrestricted open-shell Hartree-Fock total energy associated to the total and spin density matrices (whose upper triangular part is stored in packed form in PDM) of a given molecular system, POEFM, PTEFM and PEMS containing respectively the core hamiltonian matrix, the two-electron part of the Fock matrix associated to the total density matrix and the exchange matrix associated to the spin density matrix (all stored similarly).
  USE data_parameters
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PEMS,PTDM,PSDM
  DOUBLE PRECISION :: ETOT

  ETOT=ELECTRONIC_ENERGY_UHF(POEFM,PTEFM,PEMS,PTDM,PSDM,N)+INTERNUCLEAR_ENERGY
END FUNCTION ENERGY_UHF

FUNCTION ENERGY_RGHF(POEFM,PTEFM,PDM,N) RESULT(ETOT)
! Function that computes the real general Hartree-Fock total energy associated to a density matrix (whose upper triangular parts are stored in packed form) of a given molecular system, POEFM and PTEFM containing respectively the core hamiltonian matrix and the two-electron part of the Fock matrix (both stored similarly).
  USE data_parameters
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
  DOUBLE PRECISION :: ETOT

  ETOT=ELECTRONIC_ENERGY_HF(POEFM,PTEFM,PDM,N)+INTERNUCLEAR_ENERGY
END FUNCTION ENERGY_RGHF

FUNCTION ENERGY_DIFF_RGHF(POEFM,PDM1,PDM2,PHI,N) RESULT(ETOT)
  ! computes the energy difference between PDM1 and PDM2.
  USE case_parameters ; USE data_parameters ; USE basis_parameters
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PDM1,PDM2
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PTEFM1,PTEFM2p1
  DOUBLE PRECISION :: ETOT
  TYPE(gaussianbasisfunction),DIMENSION(:) :: PHI

  CALL BUILDTEFM_RGHF(PTEFM1,N,PHI,PDM1)
  CALL BUILDTEFM_RGHF(PTEFM2p1,N,PHI,PDM1+PDM2)

  ETOT = ENERGY_RGHF(POEFM,PTEFM2p1,PDM2-PDM1,N)
  ! ETOT = ETOT + ENERGY_RGHF(0*POEFM,PTEFM2m1,PDM2,N)

  
  ETOT = ENERGY_RHF(POEFM,PTEFM2p1,PDM2-PDM1,N)
END FUNCTION ENERGY_DIFF_RGHF

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

FUNCTION ELECTRONIC_ENERGY_HF(POEFM,PTEFM,PDM,N) RESULT(ENERGY)
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
END FUNCTION ELECTRONIC_ENERGY_HF

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

FUNCTION ELECTRONIC_ENERGY_UHF(POEFM,PTEFM,PEMS,PTDM,PSDM,N) RESULT(ENERGY)
! Function that computes the unrestricted open-shell Hartree-Fock electronic energy associated to the total and spin density matrices (whose upper triangular parts are stored in packed form) of a given molecular system, POEFM, PTEFM and PEMS containing respectively the core hamiltonian matrix, the two-electron part of the Fock matrix associated to the total density matrix and the exchange matrix associated to the spin density matrix (all stored similarly).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PEMS,PTDM,PSDM
  DOUBLE PRECISION :: ENERGY

  INTEGER :: I,J,IJ
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PM

  PM=POEFM+0.25D0*PTEFM
  ENERGY=0.D0
  IJ=0
  DO J=1,N
     DO I=1,J
        IJ=IJ+1
        IF (I.NE.J) THEN
           ENERGY=ENERGY+2.D0*(PM(IJ)*PTDM(IJ)-0.25D0*PEMS(IJ)*PSDM(IJ))
        ELSE
           ENERGY=ENERGY+PM(IJ)*PTDM(IJ)-0.25D0*PEMS(IJ)*PSDM(IJ)
        END IF
     END DO
  END DO
END FUNCTION ELECTRONIC_ENERGY_UHF

FUNCTION FREEENERGY(POEFM,PTEFM,PDM,N,TEMPERATURE,FUNC) RESULT(ENRG)
! Function that computes the free energy of a Hartree model with temperature associated to a density matrix (whose upper triangular part is stored in packed form in PDM) of a given system, POEFM and PTEFM containing respectively the (core hamiltonian) matrix and the (two-electron) part of the Fock matrix (both stored similarly), FUNC being the entropy generating function defining the entropy functional.
  USE matrix_tools ; USE metric_nonrelativistic
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: POEFM,PTEFM,PDM
  DOUBLE PRECISION,INTENT(IN) :: TEMPERATURE
  INTERFACE
    DOUBLE PRECISION FUNCTION FUNC(X)
      DOUBLE PRECISION,INTENT(IN) :: X
    END FUNCTION
  END INTERFACE
  DOUBLE PRECISION :: ENRG
  
  INTEGER :: I,J,IJ,INFO
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PTMP
  DOUBLE PRECISION,DIMENSION(N) :: W
  DOUBLE PRECISION,DIMENSION(3*N) :: WORK
  DOUBLE PRECISION,DIMENSION(N,N) :: Z,TMP

  DOUBLE PRECISION :: ENTTERM

  PTMP=PDM
  CALL DSPEV('V','U',N,PTMP,W,Z,N,WORK,INFO)
  IF (INFO/=0) GO TO 1
  TMP=0.D0
  DO I=1,N
     TMP(I,I)=FUNC(W(I))
  END DO
  PTMP=PACK(MATMUL(Z,MATMUL(TMP,TRANSPOSE(Z))),N)

  ENRG=ELECTRONIC_ENERGY_HF(POEFM,PTEFM,PDM,N)+TEMPERATURE*TRACEOFPRODUCT(PTMP,PS,N)
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DSPEV: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DSPEV: the algorithm failed to converge; ',INFO,'off-diagonal elements &
       &of an intermediate tridiagonal form did not converge to zero'
  END IF
  WRITE(*,*)'(called from function FREEENERGY)'
  STOP
END FUNCTION FREEENERGY
END MODULE
