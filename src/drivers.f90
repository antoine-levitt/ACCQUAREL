SUBROUTINE DRIVER_relativistic
  !$ USE omp_lib
  USE case_parameters ; USE basis_parameters ; USE scf_parameters
  USE basis ; USE integrals ; USE matrices ; USE matrix_tools
  USE metric_relativistic ; USE scf_algorithms
  INTEGER :: NBAST
  TYPE(twospinor),DIMENSION(:),ALLOCATABLE :: PHI
  INTEGER,DIMENSION(:),ALLOCATABLE :: NBAS
  TYPE(gaussianbasisfunction),DIMENSION(:),ALLOCATABLE :: GBF
  INTEGER,DIMENSION(2) :: NGBF
  INTEGER :: TYPENMBR(3),INFO,I
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: EIG
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: POEFM
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE,TARGET :: POM,PIOM,PSROM,PISROM,PCFOM
  DOUBLE COMPLEX,DIMENSION(:,:),ALLOCATABLE :: EIGVEC
  LOGICAL :: RESUME

! POUR L'INSTANT
  RESUME=.FALSE.

! Computation of the discretization basis
  CALL FORMBASIS(PHI,NBAS,GBF,NGBF)
  NBAST=SUM(NBAS)
! Computations of the tensors relative to the metric
! - computation and assembly of the overlap matrix
  WRITE(*,'(/,a)')'* Computation and assembly of the overlap matrix'
  ALLOCATE(POM(1:NBAST*(NBAST+1)/2))
  CALL BUILDOM(POM,PHI,NBAST,NBAS) ; PS=>POM
! - computation of the inverse of the overlap matrix
  ALLOCATE(PIOM(1:NBAST*(NBAST+1)/2))
  PIOM=INVERSE(POM,NBAST) ; PIS=>PIOM
! - computation of the square root of the overlap matrix
  ALLOCATE(PSROM(1:NBAST*(NBAST+1)/2))
  PSROM=SQUARE_ROOT(POM,NBAST) ; PSRS=>PSROM
! - computation of the inverse of the square root of the overlap matrix
  ALLOCATE(PISROM(1:NBAST*(NBAST+1)/2))
  PISROM=INVERSE(PSROM,NBAST) ; PISRS=>PISROM
! - computation of the Cholesky factorization of the overlap matrix
  ALLOCATE(PCFOM(1:NBAST*(NBAST+1)/2))
  PCFOM=POM
  CALL ZPPTRF('U',NBAST,PCFOM,INFO) ; PCFS=>PCFOM
  IF (INFO/=0) GOTO 1
! Computation and assembly of the core hamiltonian matrix
  WRITE(*,'(a)')'* Computation and assembly of the core hamiltonian matrix'
  ALLOCATE(POEFM(1:NBAST*(NBAST+1)/2))
  CALL BUILDOEFM(POEFM,PHI,NBAST,NBAS)
! Creation of the list of nonzero bielectronic integrals
  WRITE(*,'(a)')'* Creation of the list of nonzero bielectronic integrals'
  CALL BUILDBILIST(PHI,NBAS,BINMBR,TYPENMBR)
  IF (DIRECT) THEN
     IF (.NOT.USEDISK) THEN
! storage of the list of nonzero bielectronic integrals in memory
        ALLOCATE(BILIST(1:BINMBR,1:4))
        OPEN(LUNIT,form='UNFORMATTED')
        DO I=1,BINMBR
           READ(LUNIT)BILIST(I,:)
        END DO
        CLOSE(LUNIT,STATUS='DELETE')
     END IF
  ELSE IF (.NOT.DIRECT) THEN
! Precomputation of the bielectronic integrals is preferred to "on the fly" computation
! Precomputation of the bielectronic integrals involving gaussian basis functions
     WRITE(*,'(a)')'* Computation and storage in memory of the bielectronic integrals of GBF'
     CALL PRECOMPUTEGBFCOULOMBVALUES(GBF,NGBF)
     IF (SEMIDIRECT) THEN
        IF (.NOT.USEDISK) THEN
! storage the list and the type of nonzero bielectronic integrals (in order to use the precomputed GBF bielectronic integrals) in memory
           ALLOCATE(BILIST(1:BINMBR,1:4),BITYPE(1:BINMBR))
           OPEN(LUNIT,form='UNFORMATTED')
           DO I=1,BINMBR
              READ(LUNIT)BILIST(I,:),BITYPE(I)
           END DO
           CLOSE(LUNIT,STATUS='DELETE')
        END IF
     ELSE
        WRITE(*,'(a)')'* Computation of the bielectronic integrals of 2-spinors basis functions'
        IF (USEDISK) THEN
! storage of the list and values of nonzero bielectronic integrals on disk
           ALLOCATE(BILIST(1:1,1:4),BITYPE(1:1))
           OPEN(LUNIT,form='UNFORMATTED') ; OPEN(BIUNIT,form='UNFORMATTED')
           DO I=1,BINMBR
              READ(LUNIT)BILIST(1,:),BITYPE(1)
              WRITE(BIUNIT)BILIST(1,:),COULOMBVALUE(PHI(BILIST(1,1)),PHI(BILIST(1,2)),PHI(BILIST(1,3)),PHI(BILIST(1,4)),BITYPE(1))
              IF (BINMBR>=10.AND.MODULO(I,BINMBR/10)==0) WRITE(*,'(a1,i3,a6)')' ',CEILING(100.*I/BINMBR),'% done'
           END DO
           CLOSE(LUNIT,STATUS='DELETE') ; CLOSE(BIUNIT)
           DEALLOCATE(BILIST,BITYPE)
        ELSE
! storage of the list and values of nonzero bielectronic integrals in memory
           ALLOCATE(BILIST(1:BINMBR,1:4),BITYPE(1:1),CBIVALUES(1:BINMBR))
           OPEN(LUNIT,form='UNFORMATTED')
           DO I=1,BINMBR
              READ(LUNIT)BILIST(I,:),BITYPE(1)
              CBIVALUES(I)=COULOMBVALUE(PHI(BILIST(I,1)),PHI(BILIST(I,2)),PHI(BILIST(I,3)),PHI(BILIST(I,4)),BITYPE(1))
              IF (BINMBR>=10.AND.MODULO(I,BINMBR/10)==0) WRITE(*,'(a1,i3,a6)')' ',CEILING(100.*I/BINMBR),'% done'
           END DO
           CLOSE(LUNIT,STATUS='DELETE')
           DEALLOCATE(BITYPE)
        END IF
        CALL DEALLOCATE_INTEGRALS
     END IF 
  END IF
!
! SCF CYCLES
!
  ALLOCATE(EIG(1:NBAST),EIGVEC(1:NBAST,1:NBAST))
  DO I=1,NBALG
     SELECT CASE (ALG(I))
        CASE (1)
        WRITE(*,'(/,a)')' Roothaan''s algorithm'
        SELECT CASE (MODEL)
        CASE (1)
        CALL ROOTHAAN(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
        CASE (2)
        CALL ROOTHAAN_AOCOSDHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR)
        END SELECT
        CASE (2)
        WRITE(*,'(/,a)')' level-shifting algorithm'
        CALL LEVELSHIFTING(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
        CASE (3)
        WRITE(*,'(/,a)')' DIIS algorithm'
        CALL DIIS(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
        CASE (5)
        WRITE(*,'(/,a)')' Eric Sere''s algorithm'
        CALL ESA(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
     END SELECT
  END DO
  IF (DIRECT) THEN
     IF (USEDISK) THEN
        OPEN(LUNIT,form='UNFORMATTED') ; CLOSE(LUNIT,STATUS='DELETE')
     ELSE
        DEALLOCATE(BILIST)
     END IF
  ELSE
     IF (SEMIDIRECT) THEN
        IF (USEDISK) THEN
           OPEN(LUNIT,form='UNFORMATTED') ; CLOSE(LUNIT,STATUS='DELETE')
        ELSE
           DEALLOCATE(BILIST,BITYPE)
        END IF
        CALL DEALLOCATE_INTEGRALS
     ELSE
        IF (USEDISK) THEN
           OPEN(BIUNIT,form='UNFORMATTED') ; CLOSE(BIUNIT,STATUS='DELETE')
        ELSE
           DEALLOCATE(BILIST,CBIVALUES)
        END IF
     END IF
  END IF
  DEALLOCATE(NBAS,PHI,EIG,EIGVEC,POEFM,POM,PIOM,PSROM,PISROM)
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZPPTRF: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZPPTRF: the leading minor of order',INFO,'is not positive definite, &
    &and the factorization could not be completed'
  END IF
  WRITE(*,*)'(called from subroutine DRIVER)'
  STOP
END SUBROUTINE DRIVER_relativistic

SUBROUTINE DRIVER_nonrelativistic
  !$ USE omp_lib
  USE case_parameters ; USE basis_parameters ; USE scf_parameters
  USE basis ; USE integrals ; USE matrices ; USE matrix_tools ; USE metric_nonrelativistic
  USE scf_algorithms
  INTEGER :: NBAST
  TYPE(gaussianbasisfunction),DIMENSION(:),ALLOCATABLE :: PHI
  INTEGER,DIMENSION(:),ALLOCATABLE :: NBAS
  INTEGER :: INFO,I
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: POEFM,EIG
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,TARGET :: POM,PIOM,PSROM,PISROM,PCFOM
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: EIGVEC
  LOGICAL :: RESUME

! Computation of the discretization basis
  CALL FORMBASIS(PHI,NBAS)
  NBAST=SUM(NBAS)
! Computations of the tensors relative to the metric
! - computation and assembly of the overlap matrix
  WRITE(*,'(/,a)')'* Computation and assembly of the overlap matrix'
  ALLOCATE(POM(1:NBAST*(NBAST+1)/2))
  CALL BUILDOM(POM,PHI,NBAST) ; PS=>POM
! - computation of the inverse of the overlap matrix
  ALLOCATE(PIOM(1:NBAST*(NBAST+1)/2))
  PIOM=INVERSE(POM,NBAST) ; PIS=>PIOM
! - computation of the square root of the overlap matrix
  ALLOCATE(PSROM(1:NBAST*(NBAST+1)/2))
  PSROM=SQUARE_ROOT(POM,NBAST) ; PSRS=>PSROM
! - computation of the inverse of the square root of the overlap matrix
  ALLOCATE(PISROM(1:NBAST*(NBAST+1)/2))
  PISROM=INVERSE(PSROM,NBAST) ; PISRS=>PISROM
! - computation of the Cholesky factorization of the overlap matrix
  ALLOCATE(PCFOM(1:NBAST*(NBAST+1)/2))
  PCFOM=POM
  CALL DPPTRF('U',NBAST,PCFOM,INFO) ; PCFS=>PCFOM
  IF (INFO/=0) GOTO 1
! Computation and assembly of the matrix of the free hamiltonian
  WRITE(*,'(a)')'* Computation and assembly of the core hamiltonian matrix'
  ALLOCATE(POEFM(1:NBAST*(NBAST+1)/2))
  CALL BUILDOEFM(POEFM,PHI,NBAST)
! Creation of the list of the nonzero bielectronic integrals
  WRITE(*,'(a)')'* Creation of the list of nonzero bielectronic integrals'
  CALL BUILDBILIST(PHI,NBAS,BINMBR)
  IF (DIRECT) THEN
! Computation of the bielectronic integrals will be done "on the fly"
! storage of the list of nonzero bielectronic integrals in memory
     ALLOCATE(BILIST(1:BINMBR,1:4))
     OPEN(LUNIT,form='UNFORMATTED')
     DO I=1,BINMBR
        READ(LUNIT)BILIST(I,:)
     END DO
     CLOSE(LUNIT,STATUS='DELETE')
  ELSE
! Precomputation of the bielectronic integrals
     WRITE(*,'(a)')'* Computation of the bielectronic integrals of GBF basis functions'
     IF (USEDISK) THEN
        ALLOCATE(BILIST(1:1,1:4))
        OPEN(LUNIT,form='UNFORMATTED') ; OPEN(BIUNIT,form='UNFORMATTED')
	!$OMP PARALLEL DO SCHEDULE(STATIC, 1)
        DO I=1,BINMBR
           READ(LUNIT)BILIST(1,:)
           WRITE(BIUNIT)BILIST(1,:),COULOMBVALUE(PHI(BILIST(1,1)),PHI(BILIST(1,2)),PHI(BILIST(1,3)),PHI(BILIST(1,4)))
           IF (BINMBR>=10.AND.MODULO(I,BINMBR/10)==0) WRITE(*,'(a1,i3,a6)')' ',CEILING(100.*I/BINMBR),'% done'
        END DO
	!$OMP END PARALLEL DO
        DEALLOCATE(BILIST)
        CLOSE(LUNIT,STATUS='DELETE') ; CLOSE(BIUNIT)
     ELSE
        ALLOCATE(BILIST(1:BINMBR,1:4),RBIVALUES(1:BINMBR))
        OPEN(LUNIT,form='UNFORMATTED')
	!$OMP PARALLEL DO SCHEDULE(STATIC, 1)
        DO I=1,BINMBR
           READ(LUNIT)BILIST(I,:)
           RBIVALUES(I)=COULOMBVALUE(PHI(BILIST(I,1)),PHI(BILIST(I,2)),PHI(BILIST(I,3)),PHI(BILIST(I,4)))
           IF (BINMBR>=10.AND.MODULO(I,BINMBR/10)==0) WRITE(*,'(a1,i3,a6)')' ',CEILING(100.*I/BINMBR),'% done'
        END DO
	!$OMP END PARALLEL DO
        CLOSE(LUNIT,STATUS='DELETE')
     END IF
  END IF
!
! SCF CYCLES
!
  ALLOCATE(EIG(1:NBAST),EIGVEC(1:NBAST,1:NBAST))
  DO I=1,NBALG
     SELECT CASE (ALG(I))
        CASE (1)
        WRITE(*,'(/,a)')' Roothaan''s algorithm'
        SELECT CASE (MODEL)
           CASE (1)
           CALL ROOTHAAN_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (2)
           CALL ROOTHAAN_UHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (3)
           WRITE(*,*)' Not implemented yet!'
        END SELECT
        CASE (2)
        WRITE(*,'(/,a)')' level-shifting algorithm'
        SELECT CASE (MODEL)
           CASE (1)
           CALL LEVELSHIFTING(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (2)
           WRITE(*,*)' Not implemented yet!'
           CASE (3)
           WRITE(*,*)' Not implemented yet!'
        END SELECT
        CASE (3)
        WRITE(*,'(/,a)')' DIIS algorithm'
        SELECT CASE (MODEL)
           CASE (1)
           CALL DIIS(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (2)
           WRITE(*,*)' Not implemented yet!'
           CASE (3)
           WRITE(*,*)' Not implemented yet!'
        END SELECT
        CASE (4)
        WRITE(*,'(/,a)')' Optimal damping algorithm (ODA)'
        SELECT CASE (MODEL)
           CASE (1)
           CALL ODA(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (2)
           WRITE(*,*)' Not implemented yet!'
           CASE (3)
           WRITE(*,*)' Not implemented yet!'
        END SELECT
     END SELECT
  END DO
  IF (DIRECT) THEN
     DEALLOCATE(BILIST)
  ELSE
     IF (USEDISK) THEN
        OPEN(BIUNIT,form='UNFORMATTED') ; CLOSE(BIUNIT,STATUS='DELETE')
     ELSE
        DEALLOCATE(BILIST,RBIVALUES)
     END IF
  END IF
  DEALLOCATE(NBAS,PHI,EIG,EIGVEC,POEFM,POM,PIOM,PSROM,PISROM)
  RETURN
1 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZPPTRF: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZPPTRF: the leading minor of order',INFO,'is not positive definite, &
    &and the factorization could not be completed'
  END IF
  WRITE(*,*)'(called from subroutine DRIVER)'
  STOP
END SUBROUTINE DRIVER_nonrelativistic

SUBROUTINE DRIVER_boson_star
! Preliminary driver for the boson star model (G. Aki, J. Dolbeault: a Hartree model with temperature for boson stars)
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE scf_parameters
  USE basis ; USE integrals ; USE matrices ; USE matrix_tools ; USE common_functions
  USE metric_nonrelativistic ; USE scf_tools ; USE rootfinding_tools ; USE constants
  INTEGER :: NBAST,ITER,I,J,INFO
  INTEGER,TARGET :: RANK
  INTEGER,DIMENSION(:),ALLOCATABLE :: NBAS
  DOUBLE PRECISION :: MU,ETOT,ETOT1,SHIFT
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: POEFM,PTEFM,PFM,PDM,PDM1,PTMP,LAMBDA
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTTEFM,PTDM,PDMDIF
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,TARGET :: POM,PIOM,PSROM,PISROM,PCFOM,EIG
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: EIGVEC
  TYPE(gaussianbasisfunction),DIMENSION(:),ALLOCATABLE :: PHI
  LOGICAL :: NUMCONV

! boucle sur la temperature
  INTEGER :: K

! test
!  DOUBLE PRECISION :: RCOND
!  INTEGER,DIMENSION(:),ALLOCATABLE :: IPIV,IWORK
!  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PA,WORK
  DOUBLE PRECISION,PARAMETER :: SCALING=32.15530615624011D0

! Computation of the discretization basis
  CALL FORMBASIS(PHI,NBAS)
  NBAST=NBAS(1)
! Computations of the tensors relative to the metric
! - computation and assembly of the overlap matrix
  WRITE(*,'(/,a)')'* Computation and assembly of the overlap matrix'
  ALLOCATE(POM(1:NBAST*(NBAST+1)/2))
  CALL BUILDOM(POM,PHI,NBAST) ; PS=>POM
! tests
!  ALLOCATE(PA(1:NBAST*(NBAST+1)/2),WORK(1:2*NBAST),IPIV(NBAST),IWORK(NBAST))
!  PA=POM
!  CALL DSPTRF('U',NBAST,PA,IPIV,INFO)
!  CALL DSPCON('U',NBAST,PA,IPIV,NORM(POM,NBAST,'1'),RCOND,WORK,IWORK,INFO)
!  DEALLOCATE(PA,WORK,IPIV,IWORK)
!  write(*,*)' estimate of the reciprocal cond. numb. of this matrix =',RCOND
! tests
! - computation of the inverse of the overlap matrix
  ALLOCATE(PIOM(1:NBAST*(NBAST+1)/2))
  PIOM=INVERSE(POM,NBAST) !; PIS=>PIOM
! - computation of the square root of the overlap matrix
  ALLOCATE(PSROM(1:NBAST*(NBAST+1)/2))
  PSROM=SQUARE_ROOT(POM,NBAST) ; PSRS=>PSROM
! - computation of the inverse of the square root of the overlap matrix
  ALLOCATE(PISROM(1:NBAST*(NBAST+1)/2))
  PISROM=INVERSE(PSROM,NBAST) ; PISRS=>PISROM
! - computation of the Cholesky factorization of the overlap matrix
  ALLOCATE(PCFOM(1:NBAST*(NBAST+1)/2))
  PCFOM=POM
  CALL DPPTRF('U',NBAST,PCFOM,INFO) !; PCFS=>PCFOM
  IF (INFO/=0) GO TO 5
! Computation and assembly of the matrix of the free hamiltonian
  WRITE(*,'(a)')'* Computation and assembly of the hamiltonian matrix'
  ALLOCATE(POEFM(1:NBAST*(NBAST+1)/2))
  CALL BUILDKPFM(POEFM,PHI,NBAST)
! Creation of the list of nonzero bielectronic integrals
  WRITE(*,'(a)')'* Creation of the list of nonzero bielectronic integrals'
  CALL BUILDBILIST(PHI,NBAS,BINMBR)
! Computation of the Coulomb integrals
  WRITE(*,'(a)')'* Computation of the bielectronic integrals'
  IF (USEDISK) THEN
     ALLOCATE(BILIST(1:1,1:4))
     OPEN(LUNIT,form='UNFORMATTED') ; OPEN(BIUNIT,form='UNFORMATTED')
     !$OMP PARALLEL DO SCHEDULE(STATIC, 1)
     DO I=1,BINMBR
        READ(LUNIT)BILIST(1,:)
        WRITE(BIUNIT)BILIST(1,:),COULOMBVALUE(PHI(BILIST(1,1)),PHI(BILIST(1,2)),PHI(BILIST(1,3)),PHI(BILIST(1,4)))
        IF (BINMBR>=10.AND.MODULO(I,BINMBR/10)==0) WRITE(*,'(a1,i3,a6)')' ',CEILING(100.*I/BINMBR),'% done'
     END DO
     !$OMP END PARALLEL DO
     DEALLOCATE(BILIST)
     CLOSE(LUNIT,STATUS='DELETE') ; CLOSE(BIUNIT)
  ELSE
     ALLOCATE(BILIST(1:BINMBR,1:4),RBIVALUES(1:BINMBR))
     OPEN(LUNIT,form='UNFORMATTED')
     !$OMP PARALLEL DO SCHEDULE(STATIC, 1)
     DO I=1,BINMBR
        READ(LUNIT)BILIST(I,:)
        RBIVALUES(I)=COULOMBVALUE(PHI(BILIST(I,1)),PHI(BILIST(I,2)),PHI(BILIST(I,3)),PHI(BILIST(I,4)))
        IF (BINMBR>=10.AND.MODULO(I,BINMBR/10)==0) WRITE(*,'(a1,i3,a6)')' ',CEILING(100.*I/BINMBR),'% done'
     END DO
     !$OMP END PARALLEL DO
     CLOSE(LUNIT,STATUS='DELETE')
  END IF
!
! SCF CYCLES
!
  ALLOCATE(EIG(1:NBAST),EIGVEC(1:NBAST,1:NBAST))
  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PDM1(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2),PTMP(1:NBAST*(NBAST+1)/2))
  ALLOCATE(LAMBDA(1:NBAST))
  ALLOCATE(PTTEFM(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2),PDMDIF(1:NBAST*(NBAST+1)/2))

! Initialisation
  CALL EIGENSOLVER(POEFM,PCFOM,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 6
  OPEN(33)
  DO I=1,NBAST
     write(33,*)EIG(I)
  END DO
  CLOSE(33)
  RANK=1 ; CALL FORMDM(PDM,EIGVEC,NBAST,1,1) ; PDM=MASS*PDM
  PDM1=0.D0
  LAMBDA=0.D0 ;
!  IF (TEMPERATURE==0) LAMBDA(1)=MASS
  LAMBDA(1)=MASS
  ETOT1=0.D0
  CALL BUILDCOULOMB(PTEFM,NBAST,PHI,PDM)
  PTEFM=KAPPA*MASS*PTEFM/(4.D0*PI)

  DO K=0,1000
! TEMPERATURE LOOP
     TEMPERATURE=DBLE(K)*0.0001
  ITER=0
! ROOTHAAN'S ALGORITHM LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER
! Assembly and diagonalization of the Fock matrix associated to the density matrix
  PFM=POEFM+PTEFM
  CALL EIGENSOLVER(PFM,PCFOM,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 6
  IF (TEMPERATURE>0.D0) THEN
! Computation (using the bisection method) of the lagrangian multiplier $\mu$ associated to the constraint on the trace of the density matrix
     MU_I=>EIG ; RANK_P=>RANK
     MU=RTBIS(FUNCFORMU,EIG(1),EIG(NBAST),1.D-16)
     WRITE(*,*)'mu=',MU
     WRITE(*,*)'Residual f(mu)=',FUNCFORMU(MU)
! Computation of the updated occupation numbers
     DO I=1,NBAST
        IF (MU-EIG(I)>0.D0) THEN
           LAMBDA(I)=RECIP_DENTFUNC((MU-EIG(I))/TEMPERATURE)
           RANK=I
        ELSE
           LAMBDA(I)=0.D0
        END IF
!        WRITE(*,*)'lambda_',I,'=',LAMBDA(I)
     END DO
     write(*,*)'Rank(gamma)=',RANK,', sumlambda_i=',SUM(LAMBDA)
  END IF
! Assembly of the density matrix
  PDM1=PDM ; PDM=0.D0
  DO I=1,RANK
     CALL DSPR('U',NBAST,LAMBDA(I),EIGVEC(:,I),1,PDM)
  END DO
! Computation of the energy associated to the density matrix
  CALL BUILDCOULOMB(PTEFM,NBAST,PHI,PDM)
  PTEFM=KAPPA*MASS*PTEFM/(4.D0*PI)
  ETOT=ENERGY_HF(POEFM,PTEFM,PDM,NBAST)
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
2 WRITE(*,*)' ' ; WRITE(*,*)'Convergence after',ITER,'iterations.'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  WRITE(33,*)TEMPERATURE,ETOT,MU,RANK
! plot of the eigenfunction associated a pure state
  IF (RANK==1) THEN
     WRITE(*,*)' '
     WRITE(*,*)'The minimizer is achieved by a pure state.'
     OPEN(9,FILE='gnuplot.batch',STATUS='UNKNOWN',ACTION='WRITE')
     DO J=1,NBAST
        WRITE(9,*)'psi',J-1,'(x)=',SCALING**2,'*(\'
        DO I=1,NBAST
           IF (EIGVEC(I,J)>=0.D0) THEN
              WRITE(9,*)'+',EIGVEC(I,J),'*exp(-',FIRST_TERM*COMMON_RATIO**(I-1),'*(',SCALING,'*x)**2)\'
           ELSE
              WRITE(9,*)EIGVEC(I,J),'*exp(-',FIRST_TERM*COMMON_RATIO**(I-1),'*(',SCALING,'*x)**2)\'
           END IF
        END DO
        WRITE(9,*)')'
     END DO
     CLOSE(9)
  ELSE IF (RANK>1) THEN
     WRITE(*,*)' '
     WRITE(*,*)'The minimizer is achieved by a mixed state.'
  END IF
  GO TO 4
3 WRITE(*,*)' ' ; WRITE(*,*)'No convergence after',ITER,'iterations.'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
!  GO TO 4
4 CONTINUE
  END DO
  DEALLOCATE(LAMBDA,EIG,EIGVEC,POEFM,PTEFM,PFM,PDM,PDM1)
  DEALLOCATE(PTTEFM,PTDM,PDMDIF)
  IF (USEDISK) THEN
     OPEN(BIUNIT,form='UNFORMATTED') ; CLOSE(BIUNIT,STATUS='DELETE')
  ELSE
     DEALLOCATE(BILIST,RBIVALUES)
  END IF
  GO TO 7
  RETURN
5 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DPPTRF: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DPPTRF: the leading minor of order',INFO,'is not positive definite, &
    &and the factorization could not be completed'
  END IF
6 WRITE(*,*)'(called from subroutine DRIVER_boson_star)'
7 DEALLOCATE(POM,PIOM,PSROM,PISROM,NBAS,PHI)
END SUBROUTINE DRIVER_boson_star
