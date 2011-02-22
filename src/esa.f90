MODULE esa_mod
  USE basis_parameters
  INTEGER,POINTER :: NBAST_P
  DOUBLE COMPLEX,POINTER,DIMENSION(:) :: POEFM_P,PTDM_P,PDMDIF_P
  TYPE(twospinor),POINTER,DIMENSION(:) :: PHI_P
  CHARACTER,POINTER :: METHOD_P

CONTAINS

FUNCTION THETA(PDM,POEFM,N,PHI,METHOD) RESULT (PDMNEW)
! Function that computes the value of $\Theta(D)$ (with $D$ a hermitian matrix of size N, which upper triangular part is stored in packed form in PDM) with precision TOL (if the (geometrical) convergence of the fixed-point iterative sequence is attained), \Theta being the function defined in: A new definition of the Dirac-Fock ground state, Eric Séré, preprint (2009).
  USE case_parameters ; USE basis_parameters ; USE matrices
  USE matrix_tools ; USE metric_relativistic
  USE debug
  IMPLICIT NONE
  CHARACTER,INTENT(IN) :: METHOD
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PDM,POEFM
  TYPE(twospinor),DIMENSION(N),INTENT(IN) :: PHI
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PDMNEW

  INTEGER,PARAMETER :: ITERMAX=50
  DOUBLE PRECISION,PARAMETER :: TOL=1.D-8

  INTEGER :: ITER,INFO,LOON
  DOUBLE PRECISION :: RESIDUAL
  DOUBLE PRECISION,DIMENSION(N) :: EIG
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PTEFM,PFM,PDMOLD,PPROJM
  DOUBLE COMPLEX,DIMENSION(N,N) :: EIGVEC,DM,SF,TSF

  IF (METHOD=='N') THEN
     PDMNEW = PDM
     RETURN
  END IF

  ITER=0
  PDMOLD=PDM
  IF (THETA_CHECK) THEN
     WRITE(13,*)'tr(SD)=',REAL(TRACEOFPRODUCT(PS,PDM,N))
!     PTMP=PDM ; CALL ZHPEV('N','U',N,PTMP,EIG,EIGVEC,N,WORK,RWORK,INFO)
!     IF (INFO==0) WRITE(13,*)' Eigenvalues:',EIG
  END IF
! Fixed-point iterative sequence
  DO
     ITER=ITER+1
! Projection of the current density matrix
     CALL BUILDTEFM(PTEFM,N,PHI,PDMOLD)
     PFM=POEFM+PTEFM
     IF (METHOD=='D') THEN
! computation of the "positive" spectral projector associated to the Fock matrix based on the current density matrix via diagonalization
        CALL EIGENSOLVER(PFM,PCFS,N,EIG,EIGVEC,INFO)
        IF (INFO/=0) GO TO 1
        CALL FORMPROJ(PPROJM,EIGVEC,N,MINLOC(EIG,DIM=1,MASK=EIG>-C*C))
        PDMNEW=ABCBA(PPROJM,PS,PDMOLD,N)
     ELSEIF(METHOD=='S') THEN
! OR
! computation of the "positive" spectral projector via the sign function obtained by polynomial recursion
        DM=UNPACK(PDMOLD,N) ; SF=SIGN(PFM+C*C*PS,N)
        TSF=TRANSPOSE(CONJG(SF))
        PDMNEW=PACK((DM+MATMUL(DM,TSF)+MATMUL(SF,DM+MATMUL(DM,TSF)))/4.D0,N)
     END IF
     IF (THETA_CHECK) THEN
        WRITE(13,*)'Iter #',ITER
        WRITE(13,*)'tr(SPSDSP)=',REAL(TRACEOFPRODUCT(PS,PDMNEW,N))
!        PTMP=PDMNEW ; CALL ZHPEV('N','U',N,PTMP,EIG,EIGVEC,N,WORK,RWORK,INFO)
!        IF (INFO==0) WRITE(13,*)' Eigenvalues:',EIG
     END IF
     IF (THETA_CHECK) WRITE(13,*)'Residual =',NORM(ABA(PSRS,PDMNEW-PDMOLD,N),N,'I')
     RESIDUAL=NORM(PDMNEW-PDMOLD,N,'I')
     IF (RESIDUAL<TOL) THEN
        IF (THETA_CHECK) THEN
           WRITE(13,'(a,i3,a)')' Function THETA: convergence after ',ITER,' iteration(s).'
           WRITE(13,*)ITER,RESIDUAL
        END IF
        RETURN
     ELSE IF (ITER>=ITERMAX) THEN
        WRITE(*,*)'Warning in function THETA: no convergence after ',ITER,' iteration(s) (the residual is ',RESIDUAL,').'
        IF (THETA_CHECK) THEN
           WRITE(13,*)'Function THETA: no convergence after ',ITER,' iteration(s).'
           WRITE(13,*)'Residual =',RESIDUAL
        END IF
        RETURN
     END IF
     PDMOLD=PDMNEW
  END DO
1 WRITE(*,*)'(called from subroutine THETA)'
END FUNCTION THETA

FUNCTION SIGN(PA,N) RESULT (SA)
! Function that computes the sign function of the hermitian matrix of a selfadjoint operator, which upper triangular part is stored in packed form, using a polynomial recursion (PINVS contains the inverse of the overlap matrix, which upper triangular part is stored in packed form).
! Note: the result is an unpacked matrix due to subsequent use.
  USE matrix_tools ; USE metric_relativistic
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PA
  DOUBLE COMPLEX,DIMENSION(N,N) :: SA

  DOUBLE COMPLEX,DIMENSION(N,N) :: A,ISRS

  INTEGER,PARAMETER :: ITERMAX=50
  DOUBLE PRECISION,PARAMETER :: TOL=1.D-15
  INTEGER :: I,ITER

  A=UNPACK(PA,N) ; ISRS=UNPACK(PISRS,N)
  SA=MATMUL(UNPACK(PIS,N),A)/NORM(MATMUL(ISRS,MATMUL(A,ISRS)),N,'F')
  ITER=0
  DO
     ITER=ITER+1
     A=SA
     SA=(3.D0*SA-MATMUL(SA,MATMUL(SA,SA)))/2.D0
     IF (NORM(SA-A,N,'F')<TOL) THEN
        RETURN
     ELSE IF (ITER==ITERMAX) THEN
        WRITE(*,*)'Function SIGN: no convergence after ',ITER,' iteration(s).'
        STOP
     END IF
  END DO
END FUNCTION SIGN
END MODULE

FUNCTION DFE(LAMBDA) RESULT (ETOT)
  USE common_functions ; USE matrices ; USE matrix_tools ; USE esa_mod
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: LAMBDA
  DOUBLE PRECISION :: ETOT

  DOUBLE COMPLEX,DIMENSION(NBAST_P*(NBAST_P+1)/2) :: PDM,PTEFM,PTMP

  PDM=THETA(PTDM_P+LAMBDA*PDMDIF_P,POEFM_P,NBAST_P,PHI_P,METHOD_P)
  CALL BUILDTEFM(PTEFM,NBAST_P,PHI_P,PDM)
  ETOT=ENERGY(POEFM_P,PTEFM,PDM,NBAST_P)
END FUNCTION DFE

SUBROUTINE ESA(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! Eric Séré's Algorithm
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_relativistic ; USE scf_tools ; USE setup_tools
  USE optimization_tools ; USE esa_mod
  USE debug
  IMPLICIT NONE
  INTEGER,INTENT(IN),TARGET :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(INOUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(INOUT) :: EIGVEC
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN),TARGET :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN),TARGET :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  CHARACTER,TARGET :: METHOD
  INTEGER :: ITER,LOON,INFO,I,NUM
  DOUBLE PRECISION :: ALPHA,BETA,LAMBDA
  DOUBLE PRECISION :: DPDUMMY,TOL
  DOUBLE PRECISION :: ETOT,ETOT1,ETTOT,PDMDIFNORM
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: PTEFM,PTTEFM,PFM,PDM,PDM1
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE,TARGET :: PTDM,PDMDIF
  LOGICAL :: NUMCONV
  INTERFACE
    DOUBLE PRECISION FUNCTION DFE(LAMBDA)
      DOUBLE PRECISION,INTENT(IN) :: LAMBDA
    END FUNCTION
  END INTERFACE

! INITIALIZATIONS AND PRELIMINARIES
  OPEN(100,FILE=SETUP_FILE,STATUS='OLD',ACTION='READ')
  CALL LOOKFOR(100,'SERE''S ALGORITHM PARAMETERS',INFO)
  READ(100,'(/,a)')METHOD
  CLOSE(100)
  IF (METHOD=='D') THEN
     WRITE(*,*)'Function $Theta$ computed using diagonalization'
  ELSEIF(METHOD=='S') THEN
     WRITE(*,*)'Function $Theta$ computed using polynomial recursion'
  ELSE
     WRITE(*,*)'Function $Theta$ ignored'
  END IF
  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2),PDMDIF(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PTTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))
! Pointer assignments
  NBAST_P=>NBAST ; POEFM_P=>POEFM ; PTDM_P=>PTDM ; PDMDIF_P=>PDMDIF ; PHI_P=>PHI ; METHOD_P=>METHOD
  

  ITER=0
  TOL=1.D-4
  PDM=(0.D0,0.D0) ; PTDM=(0.D0,0.D0)
  ETOT1=0.D0
  OPEN(16,FILE='plots/esaenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/esacrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/esacrit2.txt',STATUS='unknown',ACTION='write')
  OPEN(19,FILE='plots/esaenrgyt.txt',STATUS='unknown',ACTION='write')

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
  CALL CHECKORB(EIG,NBAST,LOON)
  PDM1=PDM
  CALL FORMDM(PDM,EIGVEC,NBAST,LOON,LOON+NBE-1)
! Computation of the energy associated to the density matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT
! test ERIC
!  PFM=POEFM+PTEFM
!  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
!  CALL CHECKORB(EIG,NBAST,LOON)
! computation of the positive spectral projector associated to PFM
!  CALL FORMDM(PTMP,EIGVEC,NBAST,LOON,NBAST)
! infinity norm of the commutator between PFM and its positive spectral projector
!  WRITE(44,*)ITER,NORM(COMMUTATOR(POEFM+PTEFM,PTMP,PS,NBAST),NBAST,'I')
! fin test ERIC
! Computation of the pseudo density matrix
  IF (THETA_CHECK) WRITE(13,*)'theta(D-~D)'
  PDMDIF=THETA(PDM-PTDM,POEFM,NBAST,PHI,METHOD)
  PDMDIFNORM=NORM(ABA(PSRS,PDMDIF,NBAST),NBAST,'I')
  WRITE(*,*)'Infinity norm of the difference D_{n-1}-~D_{n-2}=',PDMDIFNORM
  IF (PDMDIFNORM<=TRSHLD) GO TO 2
  BETA=REAL(TRACEOFPRODUCT(POEFM+PTTEFM,PDMDIF,NBAST))
  write(*,*)'beta=',beta
  IF (BETA>0.D0) THEN
     WRITE(*,*)'Warning: internal computation error (beta>0).'
     GO TO 4
  ELSE
     CALL BUILDTEFM(PTTEFM,NBAST,PHI,PDMDIF)
     ALPHA=REAL(TRACEOFPRODUCT(PTTEFM,PDMDIF,NBAST))
     write(*,*)'alpha=',alpha
     IF (ALPHA>0.D0) THEN
        LAMBDA=-BETA/ALPHA
        IF (LAMBDA<1.D0) THEN
           WRITE(*,*)'lambda=',LAMBDA
! on peut raffiner avec la methode de la section doree (la fonction n'etant a priori pas quadratique), voir comment optimiser un peu cela...
!            WRITE(*,*)'refinement using golden section search (',0.D0,',',LAMBDA,',',1.D0,')'
!            DPDUMMY=GOLDEN(0.D0,LAMBDA,1.D0,DFE,TOL,LAMBDA)
           IF (THETA_CHECK) WRITE(13,*)'theta(~D+s(D-~D))'
           PTDM=THETA(PTDM+LAMBDA*PDMDIF,POEFM,NBAST,PHI,METHOD)
        ELSE
           WRITE(*,*)'lambda=1.'
           IF (THETA_CHECK) WRITE(13,*)'theta(D)'
           PTDM=THETA(PDM,POEFM,NBAST,PHI,METHOD)
        END IF
     ELSE IF (ALPHA<0.D0) THEN
        WRITE(*,*)'lambda=1.'
        IF (THETA_CHECK) WRITE(13,*)'theta(D)'
        PTDM=THETA(PDM,POEFM,NBAST,PHI,METHOD)
     ELSE
        WRITE(*,*)'Warning: internal computation error (alpha=0).'
        GO TO 4
     END IF
  END IF
! Trace of the pseudo density matrix
  WRITE(*,*)'tr(tilde{D}_n)=',REAL(TRACE(ABA(PSRS,PTDM,NBAST),NBAST))
! Energy associated to the pseudo density matrix
  CALL BUILDTEFM(PTTEFM,NBAST,PHI,PTDM)
  ETTOT=ENERGY(POEFM,PTTEFM,PTDM,NBAST)
  WRITE(19,*)ETTOT
  WRITE(*,*)'E(tilde{D}_n)=',ETTOT
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
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ESA: convergence after',ITER,'iteration(s).'
  GO TO 6
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ESA: no convergence after',ITER,'iteration(s).'
  GO TO 6
4 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine ESA: computation stopped after',ITER,'iteration(s).'
  GO TO 6
5 WRITE(*,*)'(called from subroutine ESA)'
  GO TO 7
6 OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,'(i4,e22.14)')I,EIG(I)
  END DO
  CLOSE(9)
7 NULLIFY(NBAST_P,POEFM_P,PTDM_P,PDMDIF_P,PHI_P)
  DEALLOCATE(PDM,PDM1,PTDM,PDMDIF,PTEFM,PTTEFM,PFM)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18) ; CLOSE(19)
END SUBROUTINE ESA
