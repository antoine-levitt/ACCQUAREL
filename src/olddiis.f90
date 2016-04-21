SUBROUTINE DIIS_relativistic(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! DIIS (Direct Inversion in the Iterative Subspace) algorithm (relativistic case)
! References:
! P. Pulay, Convergence acceleration of iterative sequences. The case of SCF iteration, Chem. Phys. Lett., 73(2), 393-398, 1980.
! P. Pulay, Improved SCF convergence acceleration, J. Comput. Chem., 3(4), 556-560, 1982.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_relativistic ; USE scf_tools ; USE setup_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(INOUT) :: EIG
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(INOUT) :: EIGVEC
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,LOON,INFO,I,J,IJ,K
  INTEGER :: MXSET,MSET,MM
  INTEGER,DIMENSION(:),ALLOCATABLE :: IPIV
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PTDM,PBM,DIISV
  DOUBLE COMPLEX,DIMENSION(:,:),ALLOCATABLE :: PDMSET,TMP,ISRS
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: ERRSET
  LOGICAL :: NUMCONV

! INITIALIZATIONS AND PRELIMINARIES
! Reading of the maximum dimension of the density matrix simplex
  OPEN(100,FILE=SETUP_FILE,STATUS='OLD',ACTION='READ')
  CALL LOOKFOR(100,'DIIS ALGORITHM PARAMETERS',INFO)
  READ(100,'(/,i2)')MXSET
  CLOSE(100)

  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))
  ALLOCATE(ERRSET(1:MXSET,1:NBAST,1:NBAST),PDMSET(1:MXSET,1:NBAST*(NBAST+1)/2))
  ALLOCATE(ISRS(1:NBAST,1:NBAST),TMP(1:NBAST,1:NBAST))

  ISRS=UNPACK(PISRS,NBAST)

  ITER=0
  MSET=0 ; PDMSET=(0.D0,0.D0)
  IF (RESUME) THEN
     CALL CHECKORB(EIG,NBAST,LOON)
     CALL FORMDM(PTDM,EIGVEC,NBAST,LOON,LOON+NBE-1)
  ELSE
     PTDM=(0.D0,0.D0)
     ETOT1=0.D0
  END IF
  OPEN(16,FILE='plots/diisenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/diiscrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/diiscrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PTDM)
  PFM=POEFM+PTEFM
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 5
! Assembly of the density matrix according to the aufbau principle
  CALL CHECKORB(EIG,NBAST,LOON)
  CALL FORMDM(PDM,EIGVEC,NBAST,LOON,LOON+NBE-1)
! Computation of the total energy
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'E(D_n)=',ETOT,'Ha'
! Numerical convergence check
  IF (ITER==1) THEN
     CALL CHECKNUMCONV(PDM,PDMSET(1,:),POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
  ELSE
     CALL CHECKNUMCONV(PDM,PDMSET(MSET,:),POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
  END IF
  IF (NUMCONV) THEN
! Convergence reached
     GO TO 2
  ELSE IF (ITER==MAXITR) THEN
! Maximum number of iterations reached without convergence
     GO TO 3
  ELSE
! Convergence not reached, increment
     ETOT1=ETOT
! Storage of the current density matrix, computation and storage of a new DIIS error vector
     IF (MSET==MXSET) THEN
! we have reached the maximum allowed dimension for the simplex
        PDMSET=EOSHIFT(PDMSET,SHIFT=1) ! elimination of the oldest stored density matrix
        ERRSET=EOSHIFT(ERRSET,SHIFT=1) ! elimination of the oldest stored error vector
     ELSE
        MSET=MSET+1 ! the current dimension of the simplex is increased by one
     END IF
     PDMSET(MSET,:)=PDM
     ERRSET(MSET,:,:)=COMMUTATOR(POEFM+PTEFM,PDMSET(MSET,:),PS,NBAST)
     IF (ITER==1) THEN
        PTDM=PDM
     ELSE
        WRITE(*,*)'Dimension of the density matrix simplex =',MSET
! Computation of the new pseudo-density matrix
! assembly and solving of the linear system associated to the DIIS equations
        MM=(MSET+1)*(MSET+2)/2
        ALLOCATE(PBM(1:MM))
        PBM=(0.D0,0.D0) ; PBM(MM-MSET:MM-1)=(-1.D0,0.D0)
        IJ=0
        DO J=1,MSET
           DO I=1,J
              IJ=IJ+1 ; PBM(IJ)=&
             &FINNERPRODUCT(MATMUL(ISRS,MATMUL(ERRSET(J,:,:),ISRS)),MATMUL(ISRS,MATMUL(ERRSET(I,:,:),ISRS)),NBAST)
           END DO
        END DO
        ALLOCATE(DIISV(1:MSET+1))
        DIISV(1:MSET)=(0.D0,0.D0) ; DIISV(MSET+1)=(-1.D0,0.D0)
        ALLOCATE(IPIV(1:MSET+1))
        CALL ZHPSV('U',MSET+1,1,PBM,IPIV,DIISV,MSET+1,INFO)
        DEALLOCATE(IPIV,PBM)
        IF (INFO/=0) GO TO 4
! TO DO: If the procedure fails, the data from the oldest iterations should be thrown out until the system of equations becomes solvable.
! assembly of the pseudo-density matrix
        PTDM=(0.D0,0.D0)
        DO I=1,MSET
           PTDM=PTDM+DIISV(I)*PDMSET(I,:)
        END DO
        DEALLOCATE(DIISV)
! Note that, sometimes, some of the coefficients are negative, so that the extrapolated density matrix does not belong to the convex set. However, there seems to be no benefit of constraining coefficients over a convex set within the DIIS scheme (see Kudin, Scuseria, Cancès 2002).
     END IF
     GO TO 1
  END IF
! MESSAGES
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: no convergence after',MAXITR,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
4 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZHPSV: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZHPSV: the factorization has been completed, but the block diagonal matrix is &
    &exactly singular, so the solution could not be computed'
  END IF
  GO TO 5
5 WRITE(*,*)'(called from subroutine DIIS)'
6 DEALLOCATE(PDM,PTDM,PDMSET,PTEFM,PFM,TMP,ISRS,ERRSET)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
END SUBROUTINE DIIS_relativistic

SUBROUTINE DIIS_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! DIIS (Direct Inversion in the Iterative Subspace) algorithm (restricted closed-shell Hartree-Fock formalism)
! Note: the present implementation proposes a modification (active only when the simplex dimension is set to 0) of the usual algorithm in the form of possible restarts according to an approximate linear dependence condition for the error vectors. An adaptive version of the restart condition is used whenever the parameter used to check for restart is nonpositive.
! References:
! P. Pulay, Convergence acceleration of iterative sequences. The case of SCF iteration, Chem. Phys. Lett., 73(2), 393-398, 1980.
! P. Pulay, Improved SCF convergence acceleration, J. Comput. Chem., 3(4), 556-560, 1982.
! H. F. Walker, P. Ni, Anderson acceleration for fixed-point iterations, SIAM J. Numer. Anal., 49(4), 1715-1735, 2011.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools ; USE setup_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(INOUT) :: EIG
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(INOUT) :: EIGVEC
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,I,J,IJ,K,MXSET,MSET,MM
  INTEGER :: INFO
  INTEGER,DIMENSION(:),ALLOCATABLE :: IPIV,IWORK
  DOUBLE PRECISION :: RCP,ETOT,ETOT1
  DOUBLE PRECISION :: RCOND
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PTDM,PBM,DIISV,TMPVEC
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: FPBM,WORK,FERR,BERR
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: SUMERR
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: ERRSET,DIFERRSET,PDMSET,TMP
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: B,X
  LOGICAL :: ADAPTRC,NUMCONV,EXPERT

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: ISRS

! INITIALIZATIONS AND PRELIMINARIES
  EXPERT=.FALSE.
! Reading of the maximum dimension of the density matrix simplex
  OPEN(100,FILE=SETUP_FILE,STATUS='OLD',ACTION='READ')
  CALL LOOKFOR(100,'DIIS ALGORITHM PARAMETERS',INFO)
  READ(100,'(/,i2)')MXSET
  IF (MXSET==0) THEN
! DIIS algorithm with restarts
     READ(100,'(f16.8)')RCP
     IF (RCP<=0.) THEN
        ADAPTRC=.TRUE. ! adaptive restart condition
        WRITE(*,'(a)')' (variant with restarts and adaptive condition)'
     ELSE
        ADAPTRC=.FALSE.
        WRITE(*,'(a)')' (variant with restarts)'
     END IF
  END IF
  CLOSE(100)

  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))
  IF (MXSET/=0) THEN
     ALLOCATE(ERRSET(1:MXSET,1:NBAST*NBAST),PDMSET(1:MXSET,1:NBAST*(NBAST+1)/2))
     ALLOCATE(SUMERR(1:NBAST*NBAST))
  ELSE
     ALLOCATE(TMPVEC(1:NBAST*NBAST))
! CHECK THIS
     ALLOCATE(ERRSET(1:200,1:NBAST*NBAST),PDMSET(1:200,1:NBAST*(NBAST+1)/2),DIFERRSET(1:199,1:NBAST*NBAST))
  END IF
  ALLOCATE(TMP(1:NBAST,1:NBAST))

  ITER=0
  MSET=0 ; PDMSET=0.D0
  PTDM=0.D0
  ETOT1=0.D0
  OPEN(16,FILE='plots/diisenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/diiscrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/diiscrit2.txt',STATUS='unknown',ACTION='write')
 
! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PTDM)
  PFM=POEFM+PTEFM
  ISRS=UNPACK(PISRS,NBAST)
  WRITE(*,*)'Frobenius norm of the commutator [F(~D_n),~D_n] =',&
 &NORM(MATMUL(ISRS,MATMUL(COMMUTATOR(PFM,PTDM,PS,NBAST),ISRS)),NBAST,'F')
  IF(.NOT.(RESUME .AND. ITER == 1)) THEN
     CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
     IF (INFO/=0) GO TO 5
  END IF
! Assembly of the density matrix according to the aufbau principle
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE/2)
! Computation of the total energy
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'Total energy =',ETOT,'Ha'
! Numerical convergence check
  IF (ITER==1) THEN
     CALL CHECKNUMCONV(PDM,PDMSET(1,:),POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
  ELSE
     CALL CHECKNUMCONV(PDM,PDMSET(MSET,:),POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
  END IF
  IF (NUMCONV) THEN
! Convergence reached
     GO TO 2
  ELSE IF (ITER==MAXITR) THEN
! Maximum number of iterations reached without convergence
     GO TO 3
  ELSE
! Convergence not reached, increment
     ETOT1=ETOT
! Storage of the current density matrix, computation and storage of a new DIIS error vector
     IF ((MXSET/=0).AND.(MSET==MXSET)) THEN
! we have reached the maximum allowed dimension for the simplex
        PDMSET=EOSHIFT(PDMSET,SHIFT=1) ! elimination of the oldest stored density matrix
        ERRSET=EOSHIFT(ERRSET,SHIFT=1) ! elimination of the oldest stored error vector
     ELSE
        MSET=MSET+1 ! the current dimension of the simplex is increased by one
     END IF
     PDMSET(MSET,:)=PDM
     TMP=COMMUTATOR(POEFM+PTEFM,PDM,PS,NBAST)
     IJ=0
     DO I=1,NBAST
        DO J=1,NBAST
           IJ=IJ+1
           ERRSET(MSET,IJ)=TMP(I,J)
        END DO
     END DO
     IF ((MXSET==0).AND.(MSET>1)) THEN
! restart condition
        IF (ADAPTRC) THEN
! adaptive case
           RCP=NORM2(ERRSET(1,:))**(1./MSET)
           WRITE(*,*)'tau=',RCP
           IF (RCP>1) RCP=1.D0
        END IF
!        DIFERRSET(MSET-1,:)=ERRSET(MSET,:)-ERRSET(1,:)
        DIFERRSET(MSET-1,:)=ERRSET(MSET,:)-ERRSET(MSET-1,:)
! compute a orthonormal basis using modified Gram--Schmidt's process
        IF (MSET==2) THEN
           DIFERRSET(MSET-1,:)=DIFERRSET(MSET-1,:)/NORM2(DIFERRSET(MSET-1,:))
        ELSE
           TMPVEC=DIFERRSET(MSET-1,:)
           DO I=1,MSET-2
              TMPVEC=TMPVEC-DOT_PRODUCT(TMPVEC,DIFERRSET(I,:))*DIFERRSET(I,:)
           END DO
           IF (RCP*NORM2(DIFERRSET(MSET-1,:))>NORM2(TMPVEC)) THEN
! restart
              WRITE(*,*)'Restart at step',ITER
              PDMSET(1,:)=PDMSET(MSET,:)
              ERRSET(1,:)=ERRSET(MSET,:)
              MSET=1
           ELSE
! no restart
              DIFERRSET(MSET-1,:)=TMPVEC/NORM2(TMPVEC)
           END IF
        END IF
     END IF
     WRITE(*,*)'Dimension of the density matrix simplex =',MSET
     IF (MSET==1) THEN
        PTDM=PDM
     ELSE
! Computation of the new pseudo-density matrix
! assembly and solving of the linear system associated to the DIIS equations
        MM=(MSET+1)*(MSET+2)/2
        ALLOCATE(PBM(1:MM))
        PBM=0.D0 ; PBM(MM-MSET:MM-1)=-1.D0
        IJ=0
        DO J=1,MSET
           DO I=1,J
              IJ=IJ+1 ; PBM(IJ)=PBM(IJ)+DOT_PRODUCT(ERRSET(J,:),ERRSET(I,:))
           END DO
        END DO
        ALLOCATE(DIISV(1:MSET+1))
        DIISV(1:MSET)=0.D0 ; DIISV(MSET+1)=-1.D0
        IF (EXPERT) THEN
           ALLOCATE(FPBM(1:MM),IPIV(1:MSET+1),B(MSET+1,1),X(MSET+1,1),FERR(1),BERR(1),WORK(3*MSET+3),IWORK(MSET+1))
           B(:,1)=DIISV
           CALL DSPSVX('N','U',MSET+1,1,PBM,FPBM,IPIV,B,MSET+1,X,MSET+1,RCOND,FERR,BERR,WORK,IWORK,INFO)
           IF (INFO/=0) GO TO 4
           WRITE(*,*)'Estimate of the reciprocal condition number of the matrix of the normal equations =',RCOND
           WRITE(*,*)'Estimated forward error bound =',FERR
           DIISV=X(:,1)
           DEALLOCATE(FPBM,IPIV,B,X,FERR,BERR,WORK,IWORK)
        ELSE
           ALLOCATE(IPIV(1:MSET+1))
           CALL DSPSV('U',MSET+1,1,PBM,IPIV,DIISV,MSET+1,INFO)
           IF (INFO/=0) GO TO 4
           DEALLOCATE(IPIV)
        END IF
        DEALLOCATE(PBM)
        IF (INFO/=0) GO TO 4
! TO DO: If this procedure fails with the algorithm without restarts, the data from the oldest iterations should be thrown out until the system of equations becomes solvable.
        WRITE(*,*)'c_i = ',DIISV(1:MSET)
! assembly of the pseudo-density matrix
        PTDM=0.D0; SUMERR=0.D0
        DO I=1,MSET
           PTDM=PTDM+DIISV(I)*PDMSET(I,:)
           SUMERR=SUMERR+DIISV(I)*ERRSET(I,:)
        END DO
        WRITE(*,*)'value of minimized criterion=',NORM2(SUMERR)
        DEALLOCATE(DIISV)
! Note: sometimes, some of the coefficients are negative, so that the extrapolated density matrix does not belong to the convex set. However, there seems to be no benefit of constraining coefficients over a convex set within the DIIS scheme (see Kudin, Scuseria, Cancès 2002).
     END IF
     GO TO 1
  END IF
! MESSAGES
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: no convergence after',MAXITR,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
4 IF (INFO<0) THEN
     IF (EXPERT) THEN
        WRITE(*,*)'Subroutine DSPSVX: the',-INFO,'-th argument had an illegal value.'
     ELSE
        WRITE(*,*)'Subroutine DSPSV: the',-INFO,'-th argument had an illegal value.'
     END IF
  ELSE
     IF (EXPERT) THEN
        IF (INFO<=MSET+1) THEN
           WRITE(*,*)'Subroutine DSPSVX: the factorization has been completed, but the block diagonal matrix is &
          &exactly singular, so the solution and error bounds could not be computed.'
        ELSE IF (INFO==MSET+2) THEN
           WRITE(*,*)'Subroutine DSPSVX: the block diagonal matrix is nonsingular, but the estimate of the &
          &reciprocal of its condition number is less than machine precision. The factorization has been &
          &completed, but the matrix is singular to working precision, so the solution and error bounds have &
          &not been computed.'
        END IF
     ELSE
        WRITE(*,*)'Subroutine DSPSV: the factorization has been completed, but the block diagonal matrix is &
       &exactly singular, so the solution could not be computed.'
     END IF
  END IF
  GO TO 5
5 WRITE(*,*)'(called from subroutine DIIS)'
6 DEALLOCATE(PDM,PTDM,PDMSET,PTEFM,PFM,TMP,ERRSET)
  DEALLOCATE(SUMERR)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
END SUBROUTINE DIIS_RHF

SUBROUTINE DIIS_RGHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! DIIS (Direct Inversion in the Iterative Subspace) algorithm (real general Hartree-Fock formalism)
! References:
! P. Pulay, Convergence acceleration of iterative sequences. The case of SCF iteration, Chem. Phys. Lett., 73(2), 393-398, 1980.
! P. Pulay, Improved SCF convergence acceleration, J. Comput. Chem., 3(4), 556-560, 1982.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools ; USE setup_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(INOUT) :: EIG
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(INOUT) :: EIGVEC
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,INFO,I,J,IJ,K
  INTEGER :: MXSET,MSET,MM
  INTEGER,DIMENSION(:),ALLOCATABLE :: IPIV
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PTDM,PBM,DIISV
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: ERRSET,PDMSET,TMP
  LOGICAL :: NUMCONV

! INITIALIZATIONS AND PRELIMINARIES
! Reading of the maximum dimension of the density matrix simplex
  OPEN(100,FILE=SETUP_FILE,STATUS='OLD',ACTION='READ')
  CALL LOOKFOR(100,'DIIS ALGORITHM PARAMETERS',INFO)
  READ(100,'(/,i2)')MXSET
  CLOSE(100)

  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))
  ALLOCATE(ERRSET(1:MXSET,1:NBAST*NBAST),PDMSET(1:MXSET,1:NBAST*(NBAST+1)/2),TMP(1:NBAST,1:NBAST))

  ITER=0
  MSET=0 ; PDMSET=0.D0
  PTDM=0.D0
  ETOT1=0.D0
  OPEN(16,FILE='plots/diisenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/diiscrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/diiscrit2.txt',STATUS='unknown',ACTION='write')

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Assembly and diagonalization of the Fock matrix
  CALL BUILDTEFM_RGHF(PTEFM,NBAST,PHI,PTDM)
  PFM=POEFM+PTEFM
  IF(.NOT.(RESUME .AND. ITER == 1)) THEN
     CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
     IF (INFO/=0) GO TO 5
  END IF
! Assembly of the density matrix according to the aufbau principle
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE)
! Computation of the total energy
  CALL BUILDTEFM_RGHF(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY_RGHF(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'Total energy =',ETOT,'Ha'
! Numerical convergence check
  IF (ITER==1) THEN
     CALL CHECKNUMCONV(PDM,PDMSET(1,:),POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
  ELSE
     CALL CHECKNUMCONV(PDM,PDMSET(MSET,:),POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
  END IF
  IF (NUMCONV) THEN
! Convergence reached
     GO TO 2
  ELSE IF (ITER==MAXITR) THEN
! Maximum number of iterations reached without convergence
     GO TO 3
  ELSE
! Convergence not reached, increment
     ETOT1=ETOT
! Storage of the current density matrix, computation and storage of a new DIIS error vector
     IF (MSET==MXSET) THEN
! we have reached the maximum allowed dimension for the simplex
        PDMSET=EOSHIFT(PDMSET,SHIFT=1) ! elimination of the oldest stored density matrix
        ERRSET=EOSHIFT(ERRSET,SHIFT=1) ! elimination of the oldest stored error vector
     ELSE
        MSET=MSET+1 ! the current dimension of the simplex is increased by one
     END IF
     PDMSET(MSET,:)=PDM
     TMP=COMMUTATOR(POEFM+PTEFM,PDMSET(MSET,:),PS,NBAST)
     IJ=0
     DO I=1,NBAST
        DO J=1,NBAST
           IJ=IJ+1
           ERRSET(MSET,IJ)=TMP(I,J)
        END DO
     END DO
     IF (ITER==1) THEN
        PTDM=PDM
     ELSE
        WRITE(*,*)'Dimension of the density matrix simplex =',MSET
! Computation of the new pseudo-density matrix
! assembly and solving of the linear system associated to the DIIS equations
        MM=(MSET+1)*(MSET+2)/2
        ALLOCATE(PBM(1:MM))
        PBM=0.D0 ; PBM(MM-MSET:MM-1)=-1.D0
        IJ=0
        DO J=1,MSET
           DO I=1,J
              IJ=IJ+1 ; PBM(IJ)=PBM(IJ)+DOT_PRODUCT(ERRSET(J,:),ERRSET(I,:))
           END DO
        END DO
        ALLOCATE(DIISV(1:MSET+1))
        DIISV(1:MSET)=0.D0 ; DIISV(MSET+1)=-1.D0
        ALLOCATE(IPIV(1:MSET+1))
        CALL DSPSV('U',MSET+1,1,PBM,IPIV,DIISV,MSET+1,INFO)
        DEALLOCATE(IPIV,PBM)
        IF (INFO/=0) GO TO 4
! TO DO: If the procedure fails, the data from the oldest iterations should be thrown out until the system of equations becomes solvable.
! assembly of the pseudo-density matrix
        PTDM=0.D0
        DO I=1,MSET
           PTDM=PTDM+DIISV(I)*PDMSET(I,:)
        END DO
        DEALLOCATE(DIISV)
! Note: sometimes, some of the coefficients are negative, so that the extrapolated density matrix does not belong to the convex set. However, there seems to be no benefit of constraining coefficients over a convex set within the DIIS scheme (see Kudin, Scuseria, Cancès 2002).
     END IF
     GO TO 1
  END IF
! MESSAGES
2 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: convergence after',ITER,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: no convergence after',MAXITR,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
4 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZHPSV: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZHPSV: the factorization has been completed, but the block diagonal matrix is &
    &exactly singular, so the solution could not be computed'
  END IF
  GO TO 5
5 WRITE(*,*)'(called from subroutine DIIS)'
6 DEALLOCATE(PDM,PTDM,PDMSET,PTEFM,PFM,TMP,ERRSET)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18)
END SUBROUTINE DIIS_RGHF
