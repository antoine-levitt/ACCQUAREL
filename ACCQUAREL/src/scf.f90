MODULE scf_tools
  INTERFACE CHECKNUMCONV
     MODULE PROCEDURE CHECKNUMCONV_relativistic,CHECKNUMCONV_AOCOSDHF,CHECKNUMCONV_RHF,CHECKNUMCONV_UHF
  END INTERFACE CHECKNUMCONV

CONTAINS
  SUBROUTINE CHECKORB(EIG,N,LOON)
    ! Subroutine that determines the number of the lowest and highest occupied electronic orbitals and checks if they are both in the spectral gap (in the relavistic case).
    USE case_parameters ; USE data_parameters
    IMPLICIT NONE
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
    IMPLICIT NONE
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
    IMPLICIT NONE
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

  SUBROUTINE CHECKNUMCONV_AOCOSDHF(PDMCN,PDMON,PDMCO,PDMOO,PFMC,PFMO,N,ETOTN,ETOTO,TRSHLD,NUMCONV)
    ! Subroutine that checks several numerical convergence criteria for the SCF solutions of Hartree-Fock type equations (average-of-configuration open-shell Dirac-Hartree-Fock formalism).
    USE matrix_tools ; USE metric_relativistic
    IMPLICIT NONE
    DOUBLE COMPLEX,DIMENSION(N*(N+1)/2),INTENT(IN) :: PDMCN,PDMON,PDMCO,PDMOO,PFMC,PFMO
    INTEGER,INTENT(IN) :: N
    DOUBLE PRECISION,INTENT(IN) :: ETOTN,ETOTO,TRSHLD
    LOGICAL,INTENT(OUT) :: NUMCONV

    DOUBLE COMPLEX,DIMENSION(N*(N+1)/2) :: PDIFF
    DOUBLE COMPLEX,DIMENSION(N,N) :: PTMPB
    DOUBLE PRECISION :: MXDFDNC,MXDFDNO,MXCMC,MXCMO
    LOGICAL :: CONVD,CONVC

    CONVD=.FALSE. ; CONVC=.FALSE.

    PDIFF=PDMCN-PDMCO
    MXDFDNC=NORM(PDIFF,N,'I')
    WRITE(*,*)'Infinity norm of the difference DC_n-DC_{n-1} =',MXDFDNC
    WRITE(*,*)'Frobenius norm of the difference DC_n-DC_{n-1} =',NORM(PDIFF,N,'F')
    PDIFF=PDMON-PDMOO
    MXDFDNO=NORM(PDIFF,N,'I')
    WRITE(*,*)'Infinity norm of the difference DO_n-DO_{n-1} =',MXDFDNO
    WRITE(*,*)'Frobenius norm of the difference DO_n-DO_{n-1} =',NORM(PDIFF,N,'F')
    IF ((MXDFDNC<=TRSHLD).AND.(MXDFDNO<=TRSHLD)) CONVD=.TRUE.

    PTMPB=COMMUTATOR(PFMC,PDMCN,PS,N)
    MXCMC=NORM(PTMPB,N,'I')
    WRITE(*,*)'Infinity norm of the commutator [F(DC_n),DC_n] =',MXCMC
    WRITE(*,*)'Frobenius norm of the commutator [F(DC_n),DC_n] =',NORM(PTMPB,N,'F')
    MXCMO=NORM(PTMPB,N,'I')
    WRITE(*,*)'Infinity norm of the commutator [F(DO_n),DO_n] =',MXCMO
    WRITE(*,*)'Frobenius norm of the commutator [F(DO_n),DO_n] =',NORM(PTMPB,N,'F')
    IF ((MXCMC<=TRSHLD).AND.(MXCMO<=TRSHLD)) CONVC=.TRUE.
    ! This criterion is not used to assert convergence
    WRITE(*,*)'Difference of the energies E_n-E_{n-1} =',ETOTN-ETOTO

    IF (CONVD.AND.CONVC) THEN
       NUMCONV=.TRUE.
    ELSE
       NUMCONV=.FALSE.
    END IF
  END SUBROUTINE CHECKNUMCONV_AOCOSDHF

  SUBROUTINE CHECKNUMCONV_UHF(PDMN,PDMO,N,ETOTN,ETOTO,TRSHLD,NUMCONV)
    ! Subroutine that checks several numerical convergence criteria for the SCF solutions of Hartree-Fock type equations (unrestricted open-shell Hartree-Fock formalism).
    USE matrix_tools ; USE metric_nonrelativistic
    IMPLICIT NONE
    DOUBLE PRECISION,DIMENSION(N*(N+1)/2),INTENT(IN) :: PDMN,PDMO
    INTEGER,INTENT(IN) :: N
    DOUBLE PRECISION,INTENT(IN) :: ETOTN,ETOTO,TRSHLD
    LOGICAL,INTENT(OUT) :: NUMCONV

    DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PDIFF
    DOUBLE PRECISION :: MXDFDN
    LOGICAL :: CONVD

    CONVD=.FALSE.

    PDIFF=ABA(PSRS,PDMN-PDMO,N)
    MXDFDN=NORM(PDIFF,N,'I')
    WRITE(*,*)'Infinity norm of the difference D_n-D_{n-1} =',MXDFDN
    WRITE(*,*)'Frobenius norm of the difference D_n-D_{n-1} =',NORM(PDIFF,N,'F')
    IF (MXDFDN<=TRSHLD) CONVD=.TRUE.
    ! This criterion is not used to assert convergence
    WRITE(*,*)'Difference of the energies E_n-E_{n-1} =',ETOTN-ETOTO

    IF (CONVD) THEN
       NUMCONV=.TRUE.
    ELSE
       NUMCONV=.FALSE.
    END IF
  END SUBROUTINE CHECKNUMCONV_UHF
END MODULE scf_tools
