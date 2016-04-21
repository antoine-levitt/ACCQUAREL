SUBROUTINE DRIVER_relativistic
  USE case_parameters ; USE basis_parameters ; USE scf_parameters
  USE basis ; USE integrals ; USE matrices ; USE matrix_tools
  USE metric_relativistic ; USE scf_algorithms
  IMPLICIT NONE
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
! Condition number
  WRITE(*,'(a,f16.2)')'* Condition number of the overlap matrix:', NORM(PS,NBAST,'1')*NORM(PIS,NBAST,'1')
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
        OPEN(LUNIT,access='STREAM')
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
           OPEN(LUNIT,access='STREAM')
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
           OPEN(LUNIT,access='STREAM') ; OPEN(BIUNIT,access='STREAM')
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
           OPEN(LUNIT,access='STREAM')
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
  IF(RESUME) THEN
     OPEN(100,FILE='eig.txt')
     READ(100,*) EIG
     OPEN(101,FILE='eigvec_r.txt')
     OPEN(102,FILE='eigvec_i.txt')
     CALL READMATRIX(EIGVEC,NBAST,101,102)
     CLOSE(100)
     CLOSE(101)
     CLOSE(102)
  END IF
  DO I=1,NBALG
     SELECT CASE (ALG(I))
        CASE (1)
        WRITE(*,'(/,a)')' Roothaan''s algorithm'
        SELECT CASE (MODEL)
           CASE (1)
           CALL ROOTHAAN(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
!           CALL ROOTHAAN_test(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR)
           CASE (2)
           CALL ROOTHAAN_AOCOSDHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR)
           CASE (3)
           CALL ROOTHAAN(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
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
        CASE (6)
        WRITE(*,'(/,a)')' Roothaan/Gradient algorithm'
        CALL GRADIENT(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
     END SELECT
  END DO
  IF (DIRECT) THEN
     IF (USEDISK) THEN
        OPEN(LUNIT,access='STREAM') ; CLOSE(LUNIT,STATUS='DELETE')
     ELSE
        DEALLOCATE(BILIST)
     END IF
  ELSE
     IF (SEMIDIRECT) THEN
        IF (USEDISK) THEN
           OPEN(LUNIT,access='STREAM') ; CLOSE(LUNIT,STATUS='DELETE')
        ELSE
           DEALLOCATE(BILIST,BITYPE)
        END IF
        CALL DEALLOCATE_INTEGRALS
     ELSE
        IF (USEDISK) THEN
           OPEN(BIUNIT,access='STREAM') ; CLOSE(BIUNIT,STATUS='DELETE')
        ELSE
           DEALLOCATE(BILIST,CBIVALUES)
        END IF
     END IF
  END IF
  
  OPEN(100,FILE='eig.txt')
  WRITE(100,*) EIG
  OPEN(101,FILE='eigvec_r.txt')
  OPEN(102,FILE='eigvec_i.txt')
  CALL PRINTMATRIX(EIGVEC,NBAST,101,102)
  CLOSE(100)
  CLOSE(101)
  CLOSE(102)
  
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
  USE case_parameters ; USE basis_parameters ; USE scf_parameters
  USE basis ; USE integrals ; USE matrices ; USE matrix_tools ; USE metric_nonrelativistic
  USE scf_algorithms
  IMPLICIT NONE
  INTEGER :: NBAST
  TYPE(gaussianbasisfunction),DIMENSION(:),ALLOCATABLE :: PHI
  INTEGER,DIMENSION(:),ALLOCATABLE :: NBAS
  INTEGER :: INFO,I
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: POEFM,EIG
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,TARGET :: POM,PIOM,PSROM,PISROM,PCFOM
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: EIGVEC

! Computation of the discretization basis
  IF(MODEL == 4) THEN
     CALL FORMBASIS_RGHF(PHI,NBAS)
  ELSE
     CALL FORMBASIS(PHI,NBAS)
  END IF
  NBAST=SUM(NBAS)
! Computations of the tensors relative to the metric
! - computation and assembly of the overlap matrix
  WRITE(*,'(/,a)')'* Computation and assembly of the overlap matrix'
  ALLOCATE(POM(1:NBAST*(NBAST+1)/2))
  IF(MODEL == 4) THEN
     CALL BUILDOM_RGHF(POM,PHI,NBAST)
  ELSE
     CALL BUILDOM(POM,PHI,NBAST)
  END IF
  PS=>POM
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
! Condition number
  WRITE(*,'(a,f16.2)')'* Condition number of the overlap matrix:', NORM(PS,NBAST,'1')*NORM(PIS,NBAST,'1')
! Computation and assembly of the matrix of the free hamiltonian
  WRITE(*,'(a)')'* Computation and assembly of the core hamiltonian matrix'
  ALLOCATE(POEFM(1:NBAST*(NBAST+1)/2))
  IF(MODEL == 4) THEN
     CALL BUILDOEFM_RGHF(POEFM,PHI,NBAST)
  ELSE
     CALL BUILDOEFM(POEFM,PHI,NBAST)
  END IF
! Creation of the list of the nonzero bielectronic integrals
  WRITE(*,'(a)')'* Creation of the list of nonzero bielectronic integrals'
  CALL BUILDBILIST(PHI,NBAST,BINMBR)
  IF (DIRECT) THEN
! Computation of the bielectronic integrals will be done "on the fly"
! storage of the list of nonzero bielectronic integrals in memory
     ALLOCATE(BILIST(1:BINMBR,1:4))
     OPEN(LUNIT,access='STREAM')
     DO I=1,BINMBR
        READ(LUNIT)BILIST(I,:)
     END DO
     CLOSE(LUNIT,STATUS='DELETE')
  ELSE
! Precomputation of the bielectronic integrals
     WRITE(*,'(a)')'* Computation of the bielectronic integrals of GBF basis functions'
     IF (USEDISK) THEN
        ALLOCATE(BILIST(1:1,1:4))
        OPEN(LUNIT,access='STREAM') ; OPEN(BIUNIT,access='STREAM')
        DO I=1,BINMBR
           READ(LUNIT)BILIST(1,:)
           WRITE(BIUNIT)BILIST(1,:),COULOMBVALUE(PHI(BILIST(1,1)),PHI(BILIST(1,2)),PHI(BILIST(1,3)),PHI(BILIST(1,4)))
           IF (BINMBR>=10.AND.MODULO(I,BINMBR/10)==0) WRITE(*,'(a1,i3,a6)')' ',CEILING(100.*I/BINMBR),'% done'
        END DO
        DEALLOCATE(BILIST)
        CLOSE(LUNIT,STATUS='DELETE') ; CLOSE(BIUNIT)
     ELSE
        ALLOCATE(BILIST(1:BINMBR,1:4),RBIVALUES(1:BINMBR))
        OPEN(LUNIT,access='STREAM')
	!$OMP PARALLEL DO SCHEDULE(STATIC,1)
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
  IF(RESUME) THEN
     OPEN(100,FILE='eig.txt')
     READ(100,*) EIG
     OPEN(101,FILE='eigvec.txt')
     CALL READMATRIX(EIGVEC,NBAST,101)
     CLOSE(100)
     CLOSE(101)
  END IF
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
           CASE (4)
           CALL ROOTHAAN_RGHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
        END SELECT
        CASE (2)
        WRITE(*,'(/,a)')' level-shifting algorithm'
        SELECT CASE (MODEL)
           CASE (1)
           CALL LEVELSHIFTING_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (2)
           WRITE(*,*)' Not implemented yet!'
           CASE (3)
           WRITE(*,*)' Not implemented yet!'
           CASE (4)
           CALL LEVELSHIFTING_RGHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
        END SELECT
        CASE (3)
        WRITE(*,'(/,a)')' DIIS algorithm'
        SELECT CASE (MODEL)
           CASE (1)
           CALL DIIS_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (2)
           WRITE(*,*)' Not implemented yet!'
           CASE (3)
           WRITE(*,*)' Not implemented yet!'
           CASE (4)
           CALL DIIS_RGHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
        END SELECT
        CASE (4)
        WRITE(*,'(/,a)')' Optimal damping algorithm (ODA)'
        SELECT CASE (MODEL)
           CASE (1)
           CALL ODA_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (2)
           WRITE(*,*)' Not implemented yet!'
           CASE (3)
           WRITE(*,*)' Not implemented yet!'
           CASE (4)
           CALL ODA_RGHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
        END SELECT
        CASE(6)
        WRITE(*,*)' Roothaan/Gradient algorithm'
        SELECT CASE (MODEL)
           CASE (1)
           CALL GRADIENT(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
           CASE (2)
           WRITE(*,*)' Not implemented yet!'
           CASE (3)
           WRITE(*,*)' Not implemented yet!'
           CASE (4)
           CALL GRADIENT_RGHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
        END SELECT

     END SELECT
  END DO

  OPEN(100,FILE='eig.txt')
  WRITE(100,*) EIG
  OPEN(101,FILE='eigvec.txt')
  CALL PRINTMATRIX(EIGVEC,NBAST,101)
  CLOSE(100)
  CLOSE(101)
  
  IF (DIRECT) THEN
     DEALLOCATE(BILIST)
  ELSE
     IF (USEDISK) THEN
        OPEN(BIUNIT,access='STREAM') ; CLOSE(BIUNIT,STATUS='DELETE')
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
! Preliminary driver for a boson star model (G. Aki, J. Dolbeault: a Hartree model with temperature for boson stars)
! Reference: G. Aki, J. Dolbeault, and C. Sparber, Thermal effects in gravitational Hartree systems, Ann. H. Poincaré, 12(6), 1055-1079, 2011.
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE scf_parameters
  USE basis ; USE integrals ; USE matrices ; USE matrix_tools ; USE common_functions
  USE metric_nonrelativistic ; USE scf_tools ; USE rootfinding_tools ; USE constants
  USE gnufor2
  IMPLICIT NONE
  INTEGER :: ITER,RITER,I,J,INFO
  INTEGER,TARGET :: RANK,NBAST
  INTEGER,DIMENSION(:),ALLOCATABLE :: NBAS
  DOUBLE PRECISION :: MU,ETOT,ETOT1,SHIFT
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: POEFM,PTEFM,PFM,PDM,PDM1,PTMP,LAMBDA,OLDLAMBDA
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTTEFM,PTDM,PDMDIF
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE,TARGET :: POM,PIOM,PSROM,PISROM,PCFOM,EIG
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: EIGVEC
  TYPE(gaussianbasisfunction),DIMENSION(:),ALLOCATABLE :: PHI
  LOGICAL :: NUMCONV
! graphical plots
  DOUBLE PRECISION :: ITERATIONARRAY(MAXITR),ENERGYARRAY(MAXITR)
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: TEMPERATUREARRAY,MUARRAY,RANKARRAY

! loop on temperature
  INTEGER :: TEMPMAXITR,K

! test
!  DOUBLE PRECISION :: RCOND
!  INTEGER,DIMENSION(:),ALLOCATABLE :: IPIV,IWORK
!  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PA,WORK
  DOUBLE PRECISION :: LAMBDADIFF
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
  IF (INFO/=0) GO TO 3
! Computation and assembly of the matrix of the free hamiltonian
  WRITE(*,'(a)')'* Computation and assembly of the hamiltonian matrix'
  ALLOCATE(POEFM(1:NBAST*(NBAST+1)/2))
  CALL BUILDKPFM(POEFM,PHI,NBAST)
! Creation of the list of nonzero bielectronic integrals
  WRITE(*,'(a)')'* Creation of the list of nonzero bielectronic integrals'
  CALL BUILDBILIST(PHI,NBAST,BINMBR)
! Computation of the Coulomb integrals
  WRITE(*,'(a)')'* Computation of the bielectronic integrals'
  IF (USEDISK) THEN
     ALLOCATE(BILIST(1:1,1:4))
     OPEN(LUNIT,ACCESS='STREAM') ; OPEN(BIUNIT,ACCESS='STREAM')
     DO I=1,BINMBR
        READ(LUNIT)BILIST(1,:)
        WRITE(BIUNIT)BILIST(1,:),COULOMBVALUE(PHI(BILIST(1,1)),PHI(BILIST(1,2)),PHI(BILIST(1,3)),PHI(BILIST(1,4)))
        IF (BINMBR>=10.AND.MODULO(I,BINMBR/10)==0) WRITE(*,'(a1,i3,a6)')' ',CEILING(100.*I/BINMBR),'% done'
     END DO
     DEALLOCATE(BILIST)
     CLOSE(LUNIT,STATUS='DELETE') ; CLOSE(BIUNIT)
  ELSE
     ALLOCATE(BILIST(1:BINMBR,1:4),RBIVALUES(1:BINMBR))
     OPEN(LUNIT,ACCESS='STREAM')
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
  ALLOCATE(LAMBDA(1:NBAST),OLDLAMBDA(1:NBAST))
  ALLOCATE(PTTEFM(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2),PDMDIF(1:NBAST*(NBAST+1)/2))

! Initialisation
  CALL EIGENSOLVER(POEFM,PCFOM,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 4
!  OPEN(33)
!  DO I=1,NBAST
!     write(33,*)EIG(I)
!  END DO
!  CLOSE(33)
  RANK=1 ; CALL FORMDM(PDM,EIGVEC,NBAST,1,1) ; PDM=MASS*PDM
  PDM1=0.D0
  LAMBDA=0.D0 ;
!  IF (TEMPERATURE==0) LAMBDA(1)=MASS
  LAMBDA(1)=MASS
  ETOT1=0.D0
  CALL BUILDCOULOMB(PTEFM,NBAST,PHI,PDM)
  PTEFM=KAPPA*MASS*PTEFM/(4.D0*PI)

  ITERATIONARRAY=(/ (I,I=1,MAXITR) /)
  TEMPMAXITR=50
  ALLOCATE(TEMPERATUREARRAY(0:TEMPMAXITR),MUARRAY(0:TEMPMAXITR),RANKARRAY(0:TEMPMAXITR))
  DO K=0,TEMPMAXITR
! TEMPERATURE LOOP
     TEMPERATURE=DBLE(K)*1.D-4
     TEMPERATUREARRAY(K)=TEMPERATURE
     WRITE(*,*)' ' ; WRITE(*,*)'*** TEMPERATURE=',TEMPERATURE ; WRITE(*,*)' '
     ITER=0
1    ITER=ITER+1
     WRITE(*,*)' '
     WRITE(*,*)'# ITER =',ITER
     IF (TEMPERATURE>0.D0) THEN
! Computation of the Lagrange multiplier $\mu$ associated to the constraint on the trace of the density matrix
        MU_I=>EIG ; RANK_P=>NBAST
        MU=RTBIS(FUNCFORMU,EIG(1),EIG(NBAST),1.D-16)
        WRITE(*,*)'mu=',MU,' (residual = ',FUNCFORMU(MU),')'
! Computation of the updated occupation numbers
        OLDLAMBDA=LAMBDA ; LAMBDA=0.D0
        DO I=1,NBAST
           IF (MU-EIG(I)>0.D0) THEN
              LAMBDA(I)=RECIP_DENTFUNC((MU-EIG(I))/TEMPERATURE)
              RANK=I
              WRITE(*,*)'lambda_',I,'=',LAMBDA(I)
           END IF
        END DO
        write(*,*)'rank = ',RANK
     END IF
     RITER=0
! ROOTHAAN'S ALGORITHM LOOP
10   CONTINUE
     RITER=RITER+1
! Assembly and diagonalization of the Fock matrix associated to the density matrix
     PFM=POEFM+PTEFM
     CALL EIGENSOLVER(PFM,PCFOM,NBAST,EIG,EIGVEC,INFO)
     IF (INFO/=0) GO TO 4
! Assembly of the density matrix
     PDM1=PDM ; PDM=0.D0
     IF (TEMPERATURE>0.D0) THEN
        DO I=1,RANK
           CALL DSPR('U',NBAST,LAMBDA(I),EIGVEC(:,I),1,PDM)
        END DO
     ELSE
        CALL FORMDM(PDM,EIGVEC,NBAST,1,1) ; PDM=MASS*PDM
     END IF
! Computation of the energy associated to the density matrix
     CALL BUILDCOULOMB(PTEFM,NBAST,PHI,PDM)
     PTEFM=KAPPA*MASS*PTEFM/(4.D0*PI)
! Hartree energy
!     ETOT=ELECTRONIC_ENERGY_HF(POEFM,PTEFM,PDM,NBAST)
! Total free energy
     ETOT=FREEENERGY(POEFM,PTEFM,PDM,NBAST,TEMPERATURE,ENTROPY_FUNCTION)
     ENERGYARRAY(RITER)=ETOT
!     WRITE(*,*)'E(D_n)=',ETOT
! Numerical convergence check
     CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.FALSE.)
     IF (NUMCONV) THEN
! Convergence reached
!        CALL PLOT(ITERATIONARRAY(1:RITER),ENERGYARRAY(1:RITER),persist='no',terminal='png')
        CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
        GO TO 20
     ELSE IF (RITER==MAXITR) THEN
! Maximum number of iterations reached without convergence
        CALL PLOT(ITERATIONARRAY,ENERGYARRAY,persist='no',terminal='png')
        CALL CHECKNUMCONV(PDM,PDM1,POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV,.TRUE.)
        GO TO 30
     ELSE
! Convergence not reached, increment
        ETOT1=ETOT
        GO TO 10
     END IF
! MESSAGES
20   WRITE(*,*)' ' ; WRITE(*,*)'Roothaan algorithm: convergence after',RITER,'iterations.'
     WRITE(*,*)'E(D)=',ETOT
     IF (TEMPERATURE==0.D0) THEN
        WRITE(*,*)' '
        WRITE(*,*)'The predicted critical temperature is T_c =',(EIG(2)-EIG(1))/DENTFUNC(MASS)
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
     END IF
     GO TO 40
30   WRITE(*,*)' ' ; WRITE(*,*)'Roothaan algorithm: no convergence after',RITER,'iterations.'
40   CONTINUE
! END OF ROOTHAAN'S ALGORITHM LOOP
     IF (TEMPERATURE>0.D0) THEN
        LAMBDADIFF=MAXVAL(ABS(LAMBDA-OLDLAMBDA))
        WRITE(*,*)LAMBDADIFF
        IF (LAMBDADIFF<=TRSHLD) THEN
           WRITE(*,*)' ' ; WRITE(*,*)'Convergence after',ITER,'iterations.'
           IF (RANK==1) THEN
              WRITE(*,*)' '
              WRITE(*,*)'The minimizer is achieved by a pure state.'
           ELSE IF (RANK>1) THEN
              WRITE(*,*)' '
              WRITE(*,*)'The minimizer is achieved by a mixed state.'
           END IF
           GO TO 2
        ELSE IF (ITER==30) THEN
! Maximum number of iterations reached without convergence
           WRITE(*,*)' ' ; WRITE(*,*)'No convergence after',ITER,'iterations.'
           GO TO 2
        ELSE
! Convergence not reached, increment
           GO TO 1
        END IF
2       MUARRAY(K)=MU; RANKARRAY(K)=RANK
        WRITE(33,*)TEMPERATURE,ETOT,MU,RANK,ITER
     ELSE
        RANKARRAY(K)=RANK
     END IF
  END DO
  CALL PLOT(TEMPERATUREARRAY(1:TEMPMAXITR),MUARRAY(1:TEMPMAXITR),persist='no',terminal='png',filename='mu.png')
  CALL PLOT(TEMPERATUREARRAY,RANKARRAY,persist='no',terminal='png',filename='rank.png')
  DEALLOCATE(LAMBDA,OLDLAMBDA,EIG,EIGVEC,POEFM,PTEFM,PFM,PDM,PDM1)
  DEALLOCATE(PTTEFM,PTDM,PDMDIF)
  DEALLOCATE(TEMPERATUREARRAY,MUARRAY,RANKARRAY)
  IF (USEDISK) THEN
     OPEN(BIUNIT,ACCESS='STREAM') ; CLOSE(BIUNIT,STATUS='DELETE')
  ELSE
     DEALLOCATE(BILIST,RBIVALUES)
  END IF
  GO TO 5
  RETURN
3 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine DPPTRF: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine DPPTRF: the leading minor of order',INFO,'is not positive definite, &
    &and the factorization could not be completed'
  END IF
4 WRITE(*,*)'(called from subroutine DRIVER_boson_star)'
5 DEALLOCATE(POM,PIOM,PSROM,PISROM,NBAS,PHI)
END SUBROUTINE DRIVER_boson_star
