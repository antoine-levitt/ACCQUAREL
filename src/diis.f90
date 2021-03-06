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
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PTFM,PDM,PTDM,PBM,PROJDIFF
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: FPBM,TAU,WORK,FERR,BERR
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: SUMERR,RHS,COEF
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: PDMSET,PFMSET,ERRSET,DIFERRSET,TMP
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: QRM
  LOGICAL :: ADAPTRC,NUMCONV

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: ISRS

! INITIALIZATIONS AND PRELIMINARIES
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
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2),PTFM(1:NBAST*(NBAST+1)/2))
  IF (MXSET/=0) THEN
! NOTE: it is not needed to store the density matrices (nor assemble the pseudo density matrix)
     ALLOCATE(PDMSET(1:NBAST*(NBAST+1)/2,1:MXSET),PFMSET(1:NBAST*(NBAST+1)/2,1:MXSET))
     ALLOCATE(ERRSET(1:NBAST*NBAST,1:MXSET),DIFERRSET(1:NBAST*NBAST,1:MXSET-1))
     ALLOCATE(QRM(1:NBAST*NBAST,1:MXSET-1))
  ELSE
! RETHINK THESE CHOICES...
     ALLOCATE(PDMSET(1:NBAST*(NBAST+1)/2,1:100),PFMSET(1:NBAST*(NBAST+1)/2,1:100))
     ALLOCATE(ERRSET(1:NBAST*NBAST,1:100),DIFERRSET(1:NBAST*NBAST,1:99))
     ALLOCATE(QRM(1:NBAST*NBAST,1:99))
     ALLOCATE(PROJDIFF(1:NBAST*NBAST))
  END IF
  ALLOCATE(RHS(1:NBAST*NBAST),TMP(1:NBAST,1:NBAST))
! NEEDED?
  ALLOCATE(SUMERR(1:NBAST*NBAST))

  ITER=0
  MSET=0 ; PDMSET=0.D0; PFMSET=0.D0
  PTDM=0.D0
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PTDM)
  PTFM=POEFM+PTEFM
  ETOT1=0.D0
  OPEN(16,FILE='plots/diisenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/diiscrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/diiscrit2.txt',STATUS='unknown',ACTION='write')
 
  ISRS=UNPACK(PISRS,NBAST)

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER

! Diagonalization of the Fock matrix associated with the pseudo-density matrix
  WRITE(*,*)'Frobenius norm of the commutator [F(~D_n),~D_n] =',&
 &NORM(MATMUL(ISRS,MATMUL(COMMUTATOR(PTFM,PTDM,PS,NBAST),ISRS)),NBAST,'F')
  IF (.NOT.(RESUME.AND.(ITER==1))) THEN
     CALL EIGENSOLVER(PTFM,PCFS,NBAST,EIG,EIGVEC,INFO)
     IF (INFO/=0) GO TO 20
  END IF
! Assembly of the density matrix according to the Aufbau principle
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE/2)
! Computation of the total energy
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'Total energy =',ETOT,'Ha'
! Assembly of the Fock matrix associated with the density matrix
  PFM=POEFM+PTEFM
! Numerical convergence check
  IF (ITER==1) THEN
     CALL CHECKNUMCONV(PDM,PDMSET(:,1),PFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
  ELSE
     CALL CHECKNUMCONV(PDM,PDMSET(:,MSET),PFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
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
! Storage of the current density and Fock matrice, computation and storage of a new DIIS error vector
     IF ((MXSET/=0).AND.(MSET==MXSET)) THEN
! we have reached the maximum allowed dimension for the simplex
        PDMSET=EOSHIFT(PDMSET,SHIFT=1,DIM=2) ! elimination of the oldest density matrix stored
        PFMSET=EOSHIFT(PFMSET,SHIFT=1,DIM=2) ! elimination of the oldest Fock matrix stored
        ERRSET=EOSHIFT(ERRSET,SHIFT=1,DIM=2) ! elimination of the oldest error vector stored
        DIFERRSET=EOSHIFT(DIFERRSET,SHIFT=1,DIM=2) ! elimination of the oldest difference stored
     ELSE
        MSET=MSET+1 ! the current dimension of the simplex is increased by one
     END IF
     PDMSET(:,MSET)=PDM
     PFMSET(:,MSET)=PFM
     TMP=COMMUTATOR(PFM,PDM,PS,NBAST)
     IJ=0
     DO I=1,NBAST
        DO J=1,NBAST
           IJ=IJ+1
           ERRSET(IJ,MSET)=TMP(I,J)
        END DO
     END DO
     IF (MSET>1) THEN
! computation and storage of a new difference
        DIFERRSET(:,MSET-1)=ERRSET(:,MSET)-ERRSET(:,MSET-1)
        IF (MXSET==0) THEN
           IF (ADAPTRC) THEN
! adaptation the parameter of the restart condition
              RCP=(0.00001*NORM2(ERRSET(:,1)))**(1./MSET)
              WRITE(*,*)'tau=',RCP
              IF (RCP>1) RCP=1.D0
           END IF
           IF (MSET>2) THEN
! computation of the orthogonal projection of the latest difference on the linear subspace 
! generated by the previous ones by using the QR factorization computed during the previous iteration
              PROJDIFF=DIFERRSET(:,MSET-1)
              ALLOCATE(WORK(1:1))
              CALL DORMQR('L','T',NBAST*NBAST,1,MSET-2,QRM(:,1:MSET-2),NBAST*NBAST,TAU(1:MSET-2),PROJDIFF,NBAST*NBAST,WORK,1,INFO)
              CALL DORMQR('L','N',MSET-2,1,MSET-2,QRM(:,1:MSET-2),NBAST*NBAST,TAU(1:MSET-2),PROJDIFF,NBAST*NBAST,WORK,1,INFO)
              DEALLOCATE(TAU,WORK)
! check the restart condition
              IF (RCP*NORM2(DIFERRSET(:,MSET-1))>NORM2(PROJDIFF)) THEN
                 WRITE(*,*)'Restart at step',ITER
                 PDMSET(:,1)=PDMSET(:,MSET)
                 PFMSET(:,1)=PFMSET(:,MSET)
                 ERRSET(:,1)=ERRSET(:,MSET)
                 MSET=1
              END IF
           END IF
        END IF
     END IF
     WRITE(*,*)'Dimension of the density matrix simplex =',MSET
     IF (MSET==1) THEN
        PTDM=PDM
        PTFM=PFM
     ELSE IF (MSET>=2) THEN
! Computation (using a QR factorization) of the DIIS coefficients by solving the unsconstrained form of the DIIS least-squares problem
! TO DO: If this procedure fails with the algorithm without restarts, the data from the oldest iterations should be thrown out until the system of equations becomes solvable.
! computation the QR factorization
        ALLOCATE(TAU(1:MSET-1),WORK(1:MSET-1))
        QRM(:,1:MSET-1)=DIFERRSET(:,1:MSET-1)
        CALL DGEQRF(NBAST*NBAST,MSET-1,QRM(:,1:MSET-1),NBAST*NBAST,TAU,WORK,MSET-1,INFO)
        DEALLOCATE(WORK)
        IF (INFO/=0) GO TO 4
! solution of the linear system
        RHS=ERRSET(:,MSET)
        ALLOCATE(WORK(1:1))
        CALL DORMQR('L','T',NBAST*NBAST,1,MSET-1,QRM(:,1:MSET-1),NBAST*NBAST,TAU,RHS,NBAST*NBAST,WORK,1,INFO)
        DEALLOCATE(WORK)
        IF (MXSET/=0) DEALLOCATE(TAU)
        IF (INFO/=0) GO TO 5
        CALL DTRTRS('U','N','N',MSET-1,1,QRM(1:MSET-1,1:MSET-1),MSET-1,RHS(1:MSET-1),MSET-1,INFO)
        IF (INFO/=0) GO TO 6
! computation of the proper coefficients
        ALLOCATE(COEF(1:MSET))
        COEF(1)=RHS(1)
        DO I=2,MSET-1
           COEF(I)=RHS(I)-RHS(I-1)
        END DO
        COEF(MSET)=1-RHS(MSET-1)
! NOTE: sometimes, some of the coefficients are negative, so that the extrapolated density matrix does not belong to the convex set. However, there seems to be no benefit of constraining coefficients over a convex set within the DIIS scheme (see Kudin, Scuseria, Cancès 2002).
        WRITE(*,*)'c_i = ',COEF(1:MSET)
! assembly of the pseudo-density matrix and associated Fock matrix
        PTDM=0.D0; PTFM=0.D0; SUMERR=0.D0
        DO I=1,MSET
           PTDM=PTDM+COEF(I)*PDMSET(:,I)
           PTFM=PTFM+COEF(I)*PFMSET(:,I)
           SUMERR=SUMERR+COEF(I)*ERRSET(:,I)
        END DO
        WRITE(*,*)'value of minimized criterion =',NORM2(SUMERR)
        DEALLOCATE(COEF)
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
  GO TO 21
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: no convergence after',MAXITR,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 21
4 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DGEQRF: the',-INFO,'-th argument had an illegal value'
  END IF
  GO TO 20
5 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DORGQR: the',-INFO,'-th argument had an illegal value'
  END IF
  GO TO 20
6 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DTRTRS: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DTRTRS: the',INFO,'-th diagonal element of A is zero, &
               &indicating that the matrix is singular and the solutions have not been computed.'
  END IF
20 WRITE(*,*)'(called from subroutine DIIS)'
21 DEALLOCATE(PDM,PTDM,PDMSET,PFMSET,PTEFM,PFM,PTFM,TMP)
  DEALLOCATE(ERRSET,DIFERRSET,SUMERR)
  DEALLOCATE(QRM,RHS)
  IF (MXSET==0) DEALLOCATE(PROJDIFF,TAU)
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
