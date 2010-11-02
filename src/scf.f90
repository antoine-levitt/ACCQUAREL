MODULE scf_tools
INTERFACE CHECKNUMCONV
  MODULE PROCEDURE CHECKNUMCONV_relativistic,CHECKNUMCONV_AOCOSDHF,CHECKNUMCONV_RHF,CHECKNUMCONV_UHF
END INTERFACE

CONTAINS

SUBROUTINE CHECKORB(EIG,N,LOON)
! Subroutine that determines the number of the lowest and highest occupied electronic orbitals and checks if they are both in the spectral gap (in the relavistic case).
  USE case_parameters ; USE data_parameters
  DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: EIG
  INTEGER,INTENT(IN) :: N
  INTEGER,INTENT(OUT) :: LOON

  INTEGER :: HGEN
  DOUBLE PRECISION :: AC

  AC=SCALING_FACTOR*C

! Determination of the number of the lowest occupied orbital (i.e., the one relative to the first eigenvalue associated to an electronic state in the gap)
  LOON=MINLOC(EIG,DIM=1,MASK=EIG>-AC*AC)
  IF (LOON.EQ.0) THEN
     STOP'Subroutine CHECKORB: no eigenvalue associated to an electronic state.'
  ELSE
     WRITE(*,'(a,i3,a,i3,a)')' Number of the lowest occupied electronic orbital = ',LOON,'(/',N,')'
     IF (N-LOON.LT.NBE) THEN
        WRITE(*,'(a)')' Subroutine CHECKORB: there are not enough eigenvalues associated to electronic states (',N-LOON,').'
        STOP
     END IF
! Determination of the number of the highest orbital relative to an eigenvalue associated to an electronic state in the gap
     HGEN=MAXLOC(EIG,DIM=1,MASK=EIG<0.D0)
     WRITE(*,'(a,i3)')' Number of eigenvalues associated to electronic states in the gap = ',HGEN-LOON+1
     IF (HGEN-LOON+1.LT.NBE) THEN
        WRITE(*,'(a,i2,a)')' Warning: there are less than ',NBE,' eigenvalues associated to electronic states in the gap.'
     END IF
  END IF
END SUBROUTINE CHECKORB

SUBROUTINE CHECKNUMCONV_relativistic(PDMN,PDMO,PFM,N,ETOTN,ETOTO,TRSHLD,NUMCONV)
! Subroutine that checks several numerical convergence criteria for the SCF solutions of Hartree-Fock equations (restricted closed-shell Hartree-Fock and closed-shell Dirac-Hartree-Fock formalisms).
  USE matrix_tools ; USE metric_relativistic
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PDMN,PDMO,PFM
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,INTENT(IN) :: ETOTN,ETOTO,TRSHLD
  LOGICAL,INTENT(OUT) :: NUMCONV

  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PDIFF
  DOUBLE COMPLEX,DIMENSION(N,N) :: CMT,ISRS
  DOUBLE PRECISION :: FNDFDN,FNCMT
  DOUBLE PRECISION,DIMENSION(N) :: WORK
  LOGICAL :: CONVD,CONVC

  CONVD=.FALSE. ; CONVC=.FALSE.

  PDIFF=ABA(PSRS,PDMN-PDMO,N)
  WRITE(*,*)'Infinity norm of the difference D_n-D_{n-1} =',NORM(PDIFF,N,'I')
  FNDFDN=NORM(PDIFF,N,'F')
  WRITE(*,*)'Frobenius norm of the difference D_n-D_{n-1} =',FNDFDN
  WRITE(17,'(e22.14)')FNDFDN
  IF (FNDFDN<=TRSHLD) CONVD=.TRUE.

  ISRS=UNPACK(PISRS,N)
  CMT=MATMUL(ISRS,MATMUL(COMMUTATOR(PFM,PDMN,PS,N),ISRS))
  WRITE(*,*)'Infinity norm of the commutator [F(D_n),D_n] =',NORM(CMT,N,'I')
  FNCMT=NORM(CMT,N,'F')
  WRITE(*,*)'Frobenius norm of the commutator [F(D_n),D_n] =',FNCMT
  WRITE(18,'(e22.14)')FNCMT
  IF (FNCMT<=TRSHLD) CONVC=.TRUE.
! This criterion is not used to assert convergence
  WRITE(*,*)'Difference of the energies E_n-E_{n-1} =',ETOTN-ETOTO
  WRITE(16,'(e22.14)')ETOTN

  IF (CONVD.AND.CONVC) THEN
     NUMCONV=.TRUE.
  ELSE
     NUMCONV=.FALSE.
  END IF
END SUBROUTINE CHECKNUMCONV_relativistic

SUBROUTINE CHECKNUMCONV_RHF(PDMN,PDMO,PFM,N,ETOTN,ETOTO,TRSHLD,NUMCONV)
! Subroutine that checks several numerical convergence criteria for the SCF solutions of Hartree-Fock equations (restricted closed-shell Hartree-Fock formalism).
  USE matrix_tools ; USE metric_nonrelativistic
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PDMN,PDMO,PFM
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,INTENT(IN) :: ETOTN,ETOTO,TRSHLD
  LOGICAL,INTENT(OUT) :: NUMCONV

  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PDIFF
  DOUBLE PRECISION,DIMENSION(N,N) :: CMT,ISRS
  DOUBLE PRECISION :: FNDFDN,FNCMT
  DOUBLE PRECISION,DIMENSION(N) :: WORK
  LOGICAL :: CONVD,CONVC

  CONVD=.FALSE. ; CONVC=.FALSE.

  PDIFF=ABA(PSRS,PDMN-PDMO,N)
  WRITE(*,*)'Infinity norm of the difference D_n-D_{n-1} =',NORM(PDIFF,N,'I')
  FNDFDN=NORM(PDIFF,N,'F')
  WRITE(*,*)'Frobenius norm of the difference D_n-D_{n-1} =',FNDFDN
  WRITE(17,'(e22.14)')FNDFDN
  IF (FNDFDN<=TRSHLD) CONVD=.TRUE.

  ISRS=UNPACK(PISRS,N)
  CMT=MATMUL(ISRS,MATMUL(COMMUTATOR(PFM,PDMN,PS,N),ISRS))
  WRITE(*,*)'Infinity norm of the commutator [F(D_n),D_n] =',NORM(CMT,N,'I')
  FNCMT=NORM(CMT,N,'F')
  WRITE(*,*)'Frobenius norm of the commutator [F(D_n),D_n] =',FNCMT
  WRITE(18,'(e22.14)')FNCMT
  IF (FNCMT<=TRSHLD) CONVC=.TRUE.
! This criterion is not used to assert convergence
  WRITE(*,*)'Difference of the energies E_n-E_{n-1} =',ETOTN-ETOTO
  WRITE(16,'(e22.14)')ETOTN

  IF (CONVD.AND.CONVC) THEN
     NUMCONV=.TRUE.
  ELSE
     NUMCONV=.FALSE.
  END IF
END SUBROUTINE CHECKNUMCONV_RHF

SUBROUTINE CHECKNUMCONV_UHF(PDMA,PDMB,PTDMO,PFMA,PFMB,N,ETOTN,ETOTO,TRSHLD,NUMCONV)
! Subroutine that checks several numerical convergence criteria for the SCF solutions of Hartree-Fock type equations (unrestricted open-shell Hartree-Fock formalism).
  USE matrix_tools ; USE metric_nonrelativistic
  IMPLICIT NONE
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PDMA,PDMB,PTDMO,PFMA,PFMB
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,INTENT(IN) :: ETOTN,ETOTO,TRSHLD
  LOGICAL,INTENT(OUT) :: NUMCONV

  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PDIFF
  DOUBLE PRECISION,DIMENSION(N,N) :: CMT,ISRS
  DOUBLE PRECISION :: FNDFDN,FNCMTA,FNCMTB
  LOGICAL :: CONVD,CONVC

  CONVD=.FALSE. ; CONVC=.FALSE.

  PDIFF=ABA(PSRS,PDMA+PDMB-PTDMO,N)
  WRITE(*,*)'Infinity norm of the difference D_n-D_{n-1} =',NORM(PDIFF,N,'I')
  FNDFDN=NORM(PDIFF,N,'F')
  WRITE(*,*)'Frobenius norm of the difference D_n-D_{n-1} =',FNDFDN
  IF (FNDFDN<=TRSHLD) CONVD=.TRUE.

  ISRS=UNPACK(PISRS,N)
  CMT=MATMUL(ISRS,MATMUL(COMMUTATOR(PFMA,PDMA,PS,N),ISRS))
  WRITE(*,*)'Infinity norm of the commutator [F(D_n^a),D_n^a] =',NORM(CMT,N,'I')
  FNCMTA=NORM(CMT,N,'F')
  WRITE(*,*)'Frobenius norm of the commutator [F(D_n^a),D_n^a] =',FNCMTA
  CMT=MATMUL(ISRS,MATMUL(COMMUTATOR(PFMB,PDMB,PS,N),ISRS))
  WRITE(*,*)'Infinity norm of the commutator [F(D_n^b),D_n^b] =',NORM(CMT,N,'I')
  FNCMTB=NORM(CMT,N,'F')
  WRITE(*,*)'Frobenius norm of the commutator [F(D_n^b),D_n^b] =',FNCMTB
  IF ((FNCMTA<=TRSHLD).AND.(FNCMTB<=TRSHLD)) CONVC=.TRUE.
! This criterion is not used to assert convergence
  WRITE(*,*)'Difference of the energies E_n-E_{n-1} =',ETOTN-ETOTO

  IF (CONVD.AND.CONVC) THEN
     NUMCONV=.TRUE.
  ELSE
     NUMCONV=.FALSE.
  END IF
END SUBROUTINE CHECKNUMCONV_UHF

SUBROUTINE CHECKNUMCONV_AOCOSDHF(PDMCN,PDMON,PDMCO,PDMOO,PFMC,PFMO,N,ETOTN,ETOTO,TRSHLD,NUMCONV)
! Subroutine that checks several numerical convergence criteria for the SCF solutions of Hartree-Fock type equations (average-of-configuration open-shell Dirac-Hartree-Fock formalism).
  USE matrix_tools ; USE metric_relativistic
  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PDMCN,PDMON,PDMCO,PDMOO,PFMC,PFMO
  INTEGER,INTENT(IN) :: N
  DOUBLE PRECISION,INTENT(IN) :: ETOTN,ETOTO,TRSHLD
  LOGICAL,INTENT(OUT) :: NUMCONV

  DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PDIFF
  DOUBLE COMPLEX,DIMENSION(N,N) :: CMT,ISRS
  DOUBLE PRECISION :: FNDFDNC,FNDFDNO,FNCMTC,FNCMTO
  LOGICAL :: CONVD,CONVC

  CONVD=.FALSE. ; CONVC=.FALSE.

  PDIFF=PDMCN-PDMCO
  FNDFDNC=NORM(PDIFF,N,'F')
  WRITE(*,*)'Infinity norm of the difference D_n^c-DC_{n-1}^c =',NORM(PDIFF,N,'I')
  WRITE(*,*)'Frobenius norm of the difference D_n^c-DC_{n-1}^c =',FNDFDNC
  PDIFF=PDMON-PDMOO
  FNDFDNO=NORM(PDIFF,N,'F')
  WRITE(*,*)'Infinity norm of the difference D_n^o-DC_{n-1}^o =',NORM(PDIFF,N,'I')
  WRITE(*,*)'Frobenius norm of the difference D_n^o-DC_{n-1}^o =',FNDFDNO
  IF ((FNDFDNC<=TRSHLD).AND.(FNDFDNO<=TRSHLD)) CONVD=.TRUE.

  ISRS=UNPACK(PISRS,N)
  CMT=MATMUL(ISRS,MATMUL(COMMUTATOR(PFMC,PDMCN,PS,N),ISRS))
  FNCMTC=NORM(CMT,N,'F')
  WRITE(*,*)'Infinity norm of the commutator [F(D_n^c),D_n^c] =',NORM(CMT,N,'I')
  WRITE(*,*)'Frobenius norm of the commutator [F(D_n^c),D_n^c] =',FNCMTC
  CMT=MATMUL(ISRS,MATMUL(COMMUTATOR(PFMO,PDMON,PS,N),ISRS))
  FNCMTO=NORM(CMT,N,'F')
  WRITE(*,*)'Infinity norm of the commutator [F(D_n^o),D_n^o] =',NORM(CMT,N,'I')
  WRITE(*,*)'Frobenius norm of the commutator [F(D_n^o),D_n^o] =',FNCMTO
  IF ((FNCMTC<=TRSHLD).AND.(FNCMTO<=TRSHLD)) CONVC=.TRUE.
! This criterion is not used to assert convergence
  WRITE(*,*)'Difference of the energies E_n-E_{n-1} =',ETOTN-ETOTO

  IF (CONVD.AND.CONVC) THEN
     NUMCONV=.TRUE.
  ELSE
     NUMCONV=.FALSE.
  END IF
END SUBROUTINE CHECKNUMCONV_AOCOSDHF
END MODULE

MODULE graphics_tools

INTERFACE DENSITY_POINTWISE_VALUE
  MODULE PROCEDURE DENSITY_POINTWISE_VALUE_relativistic,DENSITY_POINTWISE_VALUE_nonrelativistic
END INTERFACE

INTERFACE EXPORT_DENSITY
  MODULE PROCEDURE EXPORT_DENSITY_relativistic,EXPORT_DENSITY_nonrelativistic
END INTERFACE

CONTAINS

FUNCTION DENSITY_POINTWISE_VALUE_relativistic(PDM,PHI,NBAST,POINT) RESULT(VALUE)
! Function that computes the value of the electronic density associated to a given density matrix (only the upper triangular part of this matrix is stored in packed format) at a given point of space.
  USE basis_parameters ; USE matrix_tools
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(3),INTENT(IN) :: POINT
  DOUBLE PRECISION :: VALUE

  INTEGER :: I
  DOUBLE COMPLEX,DIMENSION(NBAST,2) :: POINTWISE_VALUES
  DOUBLE COMPLEX,DIMENSION(2,2) :: MATRIX

  DO I=1,NBAST
     POINTWISE_VALUES(I,:)=TWOSPINOR_POINTWISE_VALUE(PHI(I),POINT)
  END DO
  MATRIX=MATMUL(TRANSPOSE(CONJG(POINTWISE_VALUES)),MATMUL(UNPACK(PDM,NBAST),POINTWISE_VALUES))
  VALUE=REAL(MATRIX(1,1)+MATRIX(2,2))
END FUNCTION DENSITY_POINTWISE_VALUE_relativistic

FUNCTION DENSITY_POINTWISE_VALUE_nonrelativistic(PDM,PHI,NBAST,POINT) RESULT(VALUE)
! Function that computes the value of the electronic density associated to a given density matrix (only the upper triangular part of this matrix is stored in packed format) at a given point of space.
  USE basis_parameters ; USE matrix_tools
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(3),INTENT(IN) :: POINT
  DOUBLE PRECISION :: VALUE

  INTEGER :: I
  REAL(KIND=C_DOUBLE),DIMENSION(NBAST) :: POINTWISE_VALUES

  DO I=1,NBAST
     POINTWISE_VALUES(I)=GBF_POINTWISE_VALUE(PHI(I),POINT)
  END DO
  VALUE=DOT_PRODUCT(POINTWISE_VALUES,MATMUL(UNPACK(PDM,NBAST),POINTWISE_VALUES))
END FUNCTION DENSITY_POINTWISE_VALUE_nonrelativistic

SUBROUTINE EXPORT_DENSITY_relativistic(PDM,PHI,NBAST,RMIN,RMAX,NPOINTS,FILENAME,FILEFORMAT)
  USE basis_parameters ; USE data_parameters ; USE matrices
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2) :: PDM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,INTENT(IN) :: NBAST,NPOINTS
  DOUBLE PRECISION,INTENT(IN) :: RMIN,RMAX
  CHARACTER(*),INTENT(IN) :: FILENAME
  CHARACTER(*),INTENT(IN) :: FILEFORMAT

  INTEGER :: I,J,K
  DOUBLE PRECISION :: GRID_STEPSIZE
  DOUBLE PRECISION,DIMENSION(3) :: X

  IF (NPOINTS/=1) THEN
     GRID_STEPSIZE=(RMAX-RMIN)/(NPOINTS-1)
  ELSE
     STOP
  END IF
  IF ((FILEFORMAT=='matlab').OR.(FILEFORMAT=='MATLAB')) THEN
     OPEN(UNIT=42,FILE='plots/'//FILENAME)
     DO I=0,NPOINTS-1
        DO J=0,NPOINTS-1
           DO K=0,NPOINTS-1
              X(1)=RMIN+REAL(I)*GRID_STEPSIZE
              X(2)=RMIN+REAL(J)*GRID_STEPSIZE
              X(3)=RMIN+REAL(K)*GRID_STEPSIZE
              WRITE(42,*)X(:),DENSITY_POINTWISE_VALUE(PDM,PHI,NBAST,X)
           END DO
        END DO
     END DO
     CLOSE(42)
  ELSE IF ((FILEFORMAT=='cube').OR.(FILEFORMAT=='CUBE')) THEN
     OPEN(UNIT=42,FILE='plots/'//FILENAME//'.cube')
     WRITE(42,*)'CUBE FILE.'
     WRITE(42,*)'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
     WRITE(42,*)NBN,RMIN,RMIN,RMIN
     WRITE(42,*)NPOINTS,(RMAX-RMIN)/NPOINTS,0.D0,0.D0
     WRITE(42,*)NPOINTS,0.D0,(RMAX-RMIN)/NPOINTS,0.D0
     WRITE(42,*)NPOINTS,0.D0,0.D0,(RMAX-RMIN)/NPOINTS
     DO I=1,NBN
        WRITE(42,*)Z(I),0.D0,CENTER(:,I)
     END DO
     DO I=0,NPOINTS-1
        DO J=0,NPOINTS-1
           DO K=0,NPOINTS-1
              X(1)=RMIN+REAL(I)*GRID_STEPSIZE
              X(2)=RMIN+REAL(J)*GRID_STEPSIZE
              X(3)=RMIN+REAL(K)*GRID_STEPSIZE
              WRITE(42,*)DENSITY_POINTWISE_VALUE(PDM,PHI,NBAST,X)
           END DO
        END DO
     END DO
     CLOSE(42)
  END IF
END SUBROUTINE EXPORT_DENSITY_relativistic

SUBROUTINE EXPORT_DENSITY_nonrelativistic(PDM,PHI,NBAST,RMIN,RMAX,NPOINTS,FILENAME,FILEFORMAT)
  USE basis_parameters ; USE data_parameters ; USE matrices
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2) :: PDM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,INTENT(IN) :: NBAST,NPOINTS
  DOUBLE PRECISION,INTENT(IN) :: RMIN,RMAX
  CHARACTER(*),INTENT(IN) :: FILENAME
  CHARACTER(*),INTENT(IN) :: FILEFORMAT

  INTEGER :: I,J,K
  DOUBLE PRECISION :: GRID_STEPSIZE
  DOUBLE PRECISION,DIMENSION(3) :: X

  IF (NPOINTS/=1) THEN
     GRID_STEPSIZE=(RMAX-RMIN)/(NPOINTS-1)
  ELSE
     STOP
  END IF
  IF ((FILEFORMAT=='matlab').OR.(FILEFORMAT=='MATLAB')) THEN
     OPEN(UNIT=42,FILE='plots/'//FILENAME)
     DO I=0,NPOINTS-1
        DO J=0,NPOINTS-1
           DO K=0,NPOINTS-1
              X(1)=RMIN+REAL(I)*GRID_STEPSIZE
              X(2)=RMIN+REAL(J)*GRID_STEPSIZE
              X(3)=RMIN+REAL(K)*GRID_STEPSIZE
              WRITE(42,*)X(:),DENSITY_POINTWISE_VALUE(PDM,PHI,NBAST,X)
           END DO
        END DO
     END DO
     CLOSE(42)
  ELSE IF ((FILEFORMAT=='cube').OR.(FILEFORMAT=='CUBE')) THEN
     OPEN(UNIT=42,FILE='plots/'//FILENAME//'.cube')
     WRITE(42,*)'CUBE FILE.'
     WRITE(42,*)'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
     WRITE(42,*)NBN,RMIN,RMIN,RMIN
     WRITE(42,*)NPOINTS,(RMAX-RMIN)/NPOINTS,0.D0,0.D0
     WRITE(42,*)NPOINTS,0.D0,(RMAX-RMIN)/NPOINTS,0.D0
     WRITE(42,*)NPOINTS,0.D0,0.D0,(RMAX-RMIN)/NPOINTS
     DO I=1,NBN
        WRITE(42,*)Z(I),0.D0,CENTER(:,I)
     END DO
     DO I=0,NPOINTS-1
        DO J=0,NPOINTS-1
           DO K=0,NPOINTS-1
              X(1)=RMIN+REAL(I)*GRID_STEPSIZE
              X(2)=RMIN+REAL(J)*GRID_STEPSIZE
              X(3)=RMIN+REAL(K)*GRID_STEPSIZE
              WRITE(42,*)DENSITY_POINTWISE_VALUE(PDM,PHI,NBAST,X)
           END DO
        END DO
     END DO
     CLOSE(42)
  END IF
END SUBROUTINE EXPORT_DENSITY_nonrelativistic
END MODULE

MODULE output
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ALL_EIG
  INTERFACE OUTPUT_ITER
     MODULE PROCEDURE OUTPUT_ITER_nonrelativistic, OUTPUT_ITER_relativistic
  END INTERFACE OUTPUT_ITER
  INTERFACE OUTPUT_FINALIZE
     MODULE PROCEDURE OUTPUT_FINALIZE_nonrelativistic, OUTPUT_FINALIZE_relativistic
  END INTERFACE OUTPUT_FINALIZE

  CONTAINS
  
  SUBROUTINE OUTPUT_ITER_nonrelativistic(ITER,MAXITR,PDM,PHI,NBAST,EIG,EIGVEC,ETOT)
    USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
    USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools

    INTEGER :: NBAST,ITER,MAXITR
    DOUBLE PRECISION :: ETOT
    DOUBLE PRECISION,DIMENSION(NBAST) :: EIG
    DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: EIGVEC
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2) :: PDM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST) :: PHI
    
    IF(ITER == 1) THEN
       ! Init
       ALLOCATE(ALL_EIG(NBAST,MAXITR))
    END IF
    ALL_EIG(:,ITER) = EIG
  END SUBROUTINE OUTPUT_ITER_nonrelativistic


  
  SUBROUTINE OUTPUT_FINALIZE_nonrelativistic(ITER,MAXITR,PDM,PHI,NBAST,EIG,EIGVEC,ETOT)
    USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
    USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools ; USE graphics_tools

    INTEGER :: NBAST,ITER,MAXITR,I
    DOUBLE PRECISION :: ETOT
    DOUBLE PRECISION,DIMENSION(NBAST) :: EIG
    DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: EIGVEC
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2) :: PDM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST) :: PHI
    
    OPEN(42,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
    DO I=1,ITER
       WRITE(42,*) ALL_EIG(:,I)
    END DO
    CLOSE(42)

    CALL EXPORT_DENSITY_nonrelativistic(PDM,PHI,NBAST,-12.D0,12.D0,20,'density','matlab')
  END SUBROUTINE OUTPUT_FINALIZE_nonrelativistic

  
  SUBROUTINE OUTPUT_ITER_relativistic(ITER,MAXITR,PDM,PHI,NBAST,EIG,EIGVEC,ETOT)
    USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
    USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools

    INTEGER :: NBAST,ITER,MAXITR
    DOUBLE PRECISION :: ETOT
    DOUBLE PRECISION,DIMENSION(NBAST) :: EIG
    DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: EIGVEC
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2) :: PDM
    TYPE(twospinor),DIMENSION(NBAST) :: PHI
    
    IF(ITER == 1) THEN
       ! Init
       ALLOCATE(ALL_EIG(NBAST,MAXITR))
    END IF
    ALL_EIG(:,ITER) = EIG
  END SUBROUTINE OUTPUT_ITER_relativistic
  SUBROUTINE OUTPUT_FINALIZE_relativistic(ITER,MAXITR,PDM,PHI,NBAST,EIG,EIGVEC,ETOT)
    USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
    USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools ; USE graphics_tools

    INTEGER :: NBAST,ITER,MAXITR,I
    DOUBLE PRECISION :: ETOT
    DOUBLE PRECISION,DIMENSION(NBAST) :: EIG
    DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: EIGVEC
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2) :: PDM
    TYPE(twospinor),DIMENSION(NBAST) :: PHI
    
    OPEN(42,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
    DO I=1,ITER
       WRITE(42,*) ALL_EIG(:,I)
    END DO
    CLOSE(42)

    CALL EXPORT_DENSITY_relativistic(PDM,PHI,NBAST,-12.D0,12.D0,20,'density','matlab')
  END SUBROUTINE OUTPUT_FINALIZE_relativistic
END MODULE output
