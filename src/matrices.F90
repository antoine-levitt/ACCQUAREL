SUBROUTINE FORMDM_relativistic(PDM,EIGVEC,NBAST,LOON,HOON)
! Assembly of the density matrix from selected eigenvectors associated to (occupied) electronic orbitals (only the upper triangular part of the matrix is stored in packed format).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST,LOON,HOON
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PDM
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(IN) :: EIGVEC

  INTEGER :: I

  PDM=(0.D0,0.D0)
  DO I=LOON,HOON
     CALL ZHPR('U',NBAST,1.D0,EIGVEC(:,I),1,PDM)
  END DO
END SUBROUTINE FORMDM_relativistic

SUBROUTINE FORMDM_nonrelativistic_nonorthogonal(PDM,EIGVEC,NBAST,LOON,HOON)
! Assembly of the density matrix from selected eigenvectors associated to (occupied) electronic orbitals (only the upper triangular part of the matrix is stored in packed format). Same as FORMDM_nonrelativistic, but does not expect orthogonal eigenvectors
  USE metric_nonrelativistic ; USE matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST,LOON,HOON
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PDM
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(IN) :: EIGVEC
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: EIGVEC_GS,S

  INTEGER :: I,J

  EIGVEC_GS=EIGVEC
  S=UNPACK(PS,NBAST)

  PDM=0.D0
  DO I=LOON,HOON
     DO J=LOON,I-1
        EIGVEC_GS(:,I)=EIGVEC_GS(:,I) - dot_product(EIGVEC_GS(:,J),MATMUL(S,EIGVEC_GS(:,I))) * EIGVEC_GS(:,J)
     END DO
     EIGVEC_GS(:,I) = EIGVEC_GS(:,I) / SQRT(dot_product(EIGVEC_GS(:,I),MATMUL(S,EIGVEC_GS(:,I))))
     CALL DSPR('U',NBAST,1.D0,EIGVEC(:,I),1,PDM)
  END DO
END SUBROUTINE FORMDM_nonrelativistic_nonorthogonal


SUBROUTINE FORMDM_nonrelativistic(PDM,EIGVEC,NBAST,LOON,HOON)
! Assembly of the density matrix from selected eigenvectors associated to (occupied) electronic orbitals (only the upper triangular part of the matrix is stored in packed format).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST,LOON,HOON
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PDM
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(IN) :: EIGVEC

  INTEGER :: I

  PDM=0.D0
  DO I=LOON,HOON
     CALL DSPR('U',NBAST,1.D0,EIGVEC(:,I),1,PDM)
  END DO
END SUBROUTINE FORMDM_nonrelativistic

SUBROUTINE FORMPROJ(PPROJM,EIGVEC,NBAST,LOON)
! Assembly of the matrix of the projector on the "positive" space (i.e., the electronic states) associated to a Dirac-Fock Hamiltonian (only the upper triangular part of the matrix is stored in packed format).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST,LOON
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PPROJM
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(IN) :: EIGVEC

  INTEGER :: I

  PPROJM=(0.D0,0.D0)
  DO I=0,NBAST-LOON
     CALL ZHPR('U',NBAST,1.D0,EIGVEC(:,LOON+I),1,PPROJM)
  END DO
END SUBROUTINE FORMPROJ

SUBROUTINE BUILDOM_relativistic(POM,PHI,NBAST,NBAS)
! Computation and assembly of the overlap matrix between basis functions, i.e., the Gram matrix of the basis with respect to the $L^2(\mathbb{R}^3,\mathbb{C}^4)$ inner product (only the upper triangular part of the matrix is stored in packed format).
  USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS

  INTEGER :: I,J,K,L,M
  DOUBLE COMPLEX :: VALUE

  POM=(0.D0,0.D0)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE+PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))         &
 &                           *OVERLAPVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L))
              END DO
           END DO
        END DO
        POM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE+PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))         &
 &                           *OVERLAPVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L))
              END DO
           END DO
        END DO
        POM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
END SUBROUTINE BUILDOM_relativistic

SUBROUTINE BUILDOM_nonrelativistic(POM,PHI,NBAST)
! Computation and assembly of the overlap matrix between basis functions, i.e. the Gram matrix of the basis with respacet to the $L^2(\mathbb{R}^3)$ inner product (only the upper triangular part of the matrix is stored in packed format), for the RHF, ROHF or UHF models.
  USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI

  INTEGER :: I,J

  DO J=1,NBAST
     DO I=1,J
        POM(I+(J-1)*J/2)=OVERLAPVALUE(PHI(I),PHI(J))
     END DO
  END DO
END SUBROUTINE BUILDOM_nonrelativistic

SUBROUTINE BUILDOM_RGHF(POM,PHI,NBAST)
! Computation and assembly of the block-diagonal overlap matrix between basis functions, i.e. the Gram matrix of the basis with respacet to the $L^2(\mathbb{R}^3)$ inner product (only the upper triangular part of the matrix is stored in packed format), in the real general Hartree-Fock formalism.
  USE basis_parameters ; USE integrals ; USE matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POM
  DOUBLE PRECISION,DIMENSION(NBAST/2*(NBAST/2+1)/2) :: POM_nonrelativistic
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI

  CALL BUILDOM_nonrelativistic(POM_nonrelativistic,PHI(1:NBAST/2),NBAST/2)
  CALL BUILD_BLOCK_DIAGONAL(POM,POM_nonrelativistic,NBAST/2)
END SUBROUTINE BUILDOM_RGHF

SUBROUTINE BUILDKPFM_nonrelativistic(PKPFM,PHI,NBAST)
! Computation and assembly of the kinetic part of the Fock matrix (only the upper triangular part of the matrix is stored in packed format).
  USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PKPFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI

  INTEGER :: I,J,N
  DOUBLE PRECISION :: VALUE

  DO J=1,NBAST
     DO I=1,J
        PKPFM(I+(J-1)*J/2)=KINETICVALUE(PHI(I),PHI(J))/2.D0
     END DO
  END DO
END SUBROUTINE BUILDKPFM_nonrelativistic

SUBROUTINE BUILDOEFM_relativistic(POEFM,PHI,NBAST,NBAS)
! Computation and assembly of the monoelectronic part of the Fock matrix (only the upper triangular part of the matrix is stored in packed form)
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS

  INTEGER :: I,J,K,L,M,N
  DOUBLE COMPLEX :: TMP,VALUE
  
  POEFM=(0.D0,0.D0)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 DO N=1,NBN
                    VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                            &
 &                              *Z(N)*POTENTIALVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),CENTER(:,N))
                 END DO
              END DO
           END DO
        END DO
        POEFM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  IF(MODEL==3) THEN
! CGHF model: Schrodinger kinetic energy
     DO J=1,NBAS(1)
        DO I=1,J
           VALUE=(0.D0,0.D0)
           DO K=1,2
              DO L=1,PHI(I)%nbrofcontractions(K)
                 DO M=1,PHI(J)%nbrofcontractions(K)
                    VALUE=VALUE+PHI(I)%coefficients(K,M)*PHI(J)%coefficients(K,L)                     &
 &                              *KINETICVALUE(PHI(I)%contractions(K,M),PHI(J)%contractions(K,L))/2.D0
                 END DO
              END DO
           END DO
           POEFM(I+(J-1)*J/2)=POEFM(I+(J-1)*J/2)+VALUE
        END DO
     END DO
  ELSE
     DO J=NBAS(1)+1,SUM(NBAS)
        DO I=1,NBAS(1)
           VALUE=(0.D0,0.D0)
           DO L=1,PHI(I)%nbrofcontractions(1)
              DO M=1,PHI(J)%nbrofcontractions(1)
                 VALUE=VALUE-C*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(1,L))                    &
 &                           *DCMPLX(0.D0,DERIVVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(1,L),3))
              END DO
           END DO
           DO L=1,PHI(I)%nbrofcontractions(1)
              DO M=1,PHI(J)%nbrofcontractions(2)
                 VALUE=VALUE-C*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))               &
 &                           *DCMPLX(DERIVVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L),2), &
 &                                   DERIVVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L),1))
              END DO
           END DO
           DO L=1,PHI(I)%nbrofcontractions(2)
              DO M=1,PHI(J)%nbrofcontractions(1)
                 VALUE=VALUE-C*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))                &
 &                           *DCMPLX(-DERIVVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L),2), &
 &                                   DERIVVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L),1))
              END DO
           END DO
           DO L=1,PHI(I)%nbrofcontractions(2)
              DO M=1,PHI(J)%nbrofcontractions(2)
                 VALUE=VALUE+C*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(2,L))                    &
 &                           *DCMPLX(0.D0,DERIVVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(2,L),3))
              END DO 
           END DO
           POEFM(I+(J-1)*J/2)=VALUE
        END DO
     END DO
     DO J=NBAS(1)+1,SUM(NBAS)
        DO I=NBAS(1)+1,J
           VALUE=(0.D0,0.D0)
           DO K=1,2
              DO L=1,PHI(I)%nbrofcontractions(K) 
                 DO M=1,PHI(J)%nbrofcontractions(K)
                    TMP=2.D0*C*C*OVERLAPVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L))
                    DO N=1,NBN
                       TMP=TMP+Z(N)*POTENTIALVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),CENTER(:,N))
                    END DO
                    VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))*TMP
                 END DO
              END DO
           END DO
           POEFM(I+(J-1)*J/2)=VALUE
        END DO
     END DO
  END IF
END SUBROUTINE BUILDOEFM_relativistic

SUBROUTINE BUILDOEFM_nonrelativistic(POEFM,PHI,NBAST)
! Computation and assembly of the monoelectronic part of the Fock matrix (only the upper triangular part of the matrix is stored in packed format).
  USE data_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI

  INTEGER :: I,J,N
  DOUBLE PRECISION :: VALUE

  POEFM=0.D0
  DO J=1,NBAST
     DO I=1,J
        VALUE=KINETICVALUE(PHI(I),PHI(J))/2.D0
        DO N=1,NBN
           VALUE=VALUE-Z(N)*POTENTIALVALUE(PHI(I),PHI(J),CENTER(:,N))
        END DO
        POEFM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
END SUBROUTINE BUILDOEFM_nonrelativistic

SUBROUTINE BUILDOEFM_RGHF(POEFM,PHI,NBAST)
! Computation and assembly of the monoelectronic part of the block-diagonal Fock matrix in the real general Hartree-Fock formalism (only the upper triangular part of the matrix is stored in packed format).
  USE data_parameters ; USE basis_parameters ; USE integrals ; USE matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POEFM
  DOUBLE PRECISION,DIMENSION(NBAST/2*(NBAST/2+1)/2) :: POEFM_nonrelativistic
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI

  CALL BUILDOEFM_nonrelativistic(POEFM_nonrelativistic,PHI(1:NBAST/2),NBAST/2)
  CALL BUILD_BLOCK_DIAGONAL(POEFM,POEFM_nonrelativistic,NBAST/2)
END SUBROUTINE BUILDOEFM_RGHF

SUBROUTINE BUILDTEFM_relativistic(PTEFM,NBAST,PHI,PDM)
! Computation and assembly of the bielectronic part of the Fock matrix associated to a given density matrix using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals ; use matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  INTEGER(smallint) :: I,J,K,L
  INTEGER :: N
  CHARACTER(2) :: CLASS
  DOUBLE COMPLEX :: INTGRL
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: TEFM,DM

  TEFM=(0.D0,0.D0)
  DM=UNPACK(PDM,NBAST)

  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) OPEN(LUNIT,access='STREAM')
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) OPEN(BIUNIT,access='STREAM')
  DO N=1,BINMBR
     IF (DIRECT) THEN
! the value of the bielectronic integral is computed "on the fly"
        IF (USEDISK) THEN
           READ(LUNIT)I,J,K,L
        ELSE
           I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
        END IF
        INTGRL=COULOMBVALUE(PHI(I),PHI(J),PHI(K),PHI(L))
     ELSE
        IF (SEMIDIRECT) THEN
! the value of the bielectronic integral is computed "on the fly", but using the precomputed values of the involved CGTO bielectronic integrals
           IF (USEDISK) THEN
              READ(LUNIT)I,J,K,L,CLASS
           ELSE
              I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
              CLASS=BITYPE(N)
           END IF
           INTGRL=COULOMBVALUE(PHI(I),PHI(J),PHI(K),PHI(L),CLASS)
        ELSE
           IF (USEDISK) THEN
! the value of the bielectronic integral is read on disk
              READ(BIUNIT)I,J,K,L,INTGRL
           ELSE
              I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
              INTGRL=CBIVALUES(N)
           END IF
        END IF
     END IF
!     IF ((I==J).AND.(I==K).AND.(I==L)) THEN
!       NO CONTRIBUTION
!     ELSE IF ((I==J).AND.(I/=K).AND.(I==L)) THEN
!        TEFM(I,I)=TEFM(I,I)+INTGRL*(DM(K,I)-DM(I,K))
!     ELSE IF ((I==J).AND.(I==K).AND.(I/=L)) THEN
!        TEFM(I,I)=TEFM(I,I)+INTGRL*(DM(I,L)-DM(L,I))
!     ELSE IF ((I/=J).AND.(I==K).AND.(I==L)) THEN
!        TEFM(I,I)=TEFM(I,I)+INTGRL*(DM(I,J)-DM(J,I))
!     ELSE
     IF ((I/=J).AND.(I==K).AND.(J==L)) THEN
        TEFM(I,J)=TEFM(I,J)+INTGRL*(DM(I,J)-DM(J,I))
!     ELSE IF ((I/=J).AND.(J==K).AND.(J==L)) THEN
!        TEFM(J,J)=TEFM(J,J)+INTGRL*(DM(I,J)-DM(J,I))
!     ELSE IF ((I/=J).AND.(I==K).AND.(I/=L).AND.(J/=L)) THEN
!        TEFM(I,J)=TEFM(I,J)+INTGRL*(DM(I,L)-DM(L,I))
!        TEFM(I,L)=TEFM(I,L)+INTGRL*(DM(I,J)-DM(J,I))
!     ELSE IF ((I/=J).AND.(I/=K).AND.(J/=K).AND.(J==L)) THEN
!        TEFM(I,J)=TEFM(I,J)+INTGRL*(DM(K,J)-DM(J,K))
!        TEFM(K,J)=TEFM(K,J)+INTGRL*(DM(I,J)-DM(J,I))
     ELSE
        TEFM(I,J)=TEFM(I,J)+INTGRL*DM(K,L)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DM(J,K)
        TEFM(K,L)=TEFM(K,L)+INTGRL*DM(I,J)
        TEFM(K,J)=TEFM(K,J)-INTGRL*DM(L,I)
     END IF
  END DO
  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) CLOSE(LUNIT)
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) CLOSE(BIUNIT)
  PTEFM=PACK(TEFM,NBAST)
END SUBROUTINE BUILDTEFM_relativistic

SUBROUTINE BUILDTEFM_RHF(PTEFM,NBAST,PHI,PDM)
! Computation and assembly of the two-electron part of the Fock matrix associated to a given density matrix in the restricted closed-shell Hartree-Fock formalism, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
! Note: G(D)=2J(D)-K(D), with J(D) the Coulomb term and K(D) the exchange term.
  USE scf_parameters ; USE basis_parameters ; USE integrals ; use matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: TEFM,DM
  INTEGER(smallint) :: I,J,K,L
  INTEGER :: N
  DOUBLE PRECISION :: INTGRL

  TEFM=0.D0
  DM=UNPACK(PDM,NBAST)
  IF (.NOT.DIRECT.AND.USEDISK) OPEN(BIUNIT,access='STREAM')
  DO N=1,BINMBR
     IF (DIRECT) THEN
! the values of the bielectronic integrals are computed "on the fly"
        I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
        INTGRL=COULOMBVALUE(PHI(I),PHI(J),PHI(K),PHI(L))
     ELSE
        IF (USEDISK) THEN
! the list and values of the bielectronic integrals are read on disk
           READ(BIUNIT)I,J,K,L,INTGRL
        ELSE
! the list and values of the bielectronic integrals are read in memory
           I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
           INTGRL=RBIVALUES(N)
        END IF
     END IF
! 1 value for the 4 indices
     IF ((I==J).AND.(J==K).AND.(K==L)) THEN
        TEFM(I,I)=TEFM(I,I)+INTGRL*DM(I,I)
! 2 distinct values for the 4 indices
     ELSE IF ((I>J).AND.(J==K).AND.(K==L)) THEN
        TEFM(I,J)=TEFM(I,J)+INTGRL*DM(J,J)
        TEFM(J,I)=TEFM(J,I)+INTGRL*DM(J,J)
        TEFM(J,J)=TEFM(J,J)+INTGRL*(DM(I,J)+DM(J,I))
     ELSE IF ((I==J).AND.(J==K).AND.(K>L)) THEN
        TEFM(L,I)=TEFM(L,I)+INTGRL*DM(I,I)
        TEFM(I,L)=TEFM(I,L)+INTGRL*DM(I,I)
        TEFM(I,I)=TEFM(I,I)+INTGRL*(DM(L,I)+DM(I,L))
     ELSE IF ((I==J).AND.(J>K).AND.(K==L)) THEN
        TEFM(I,I)=TEFM(I,I)+2.D0*INTGRL*DM(K,K)
        TEFM(K,K)=TEFM(K,K)+2.D0*INTGRL*DM(I,I)
        TEFM(I,K)=TEFM(I,K)-INTGRL*DM(I,K)
        TEFM(K,I)=TEFM(K,I)-INTGRL*DM(K,I)
     ELSE IF ((I==K).AND.(K>J).AND.(J==L)) THEN
        TEFM(I,J)=TEFM(I,J)+2.D0*INTGRL*DM(I,J)
        TEFM(J,I)=TEFM(J,I)+2.D0*INTGRL*DM(J,I)
        TEFM(I,J)=TEFM(I,J)+INTGRL*DM(J,I)
        TEFM(J,I)=TEFM(J,I)+INTGRL*DM(I,J)
        TEFM(I,I)=TEFM(I,I)-INTGRL*DM(J,J)
        TEFM(J,J)=TEFM(J,J)-INTGRL*DM(I,I)
! 3 distinct values for the 4 indices
     ELSE IF ((I==K).AND.(K>J).AND.(J>L)) THEN
        TEFM(I,J)=TEFM(I,J)+2.D0*INTGRL*(DM(I,L)+DM(L,I))
        TEFM(J,I)=TEFM(J,I)+2.D0*INTGRL*(DM(I,L)+DM(L,I))
        TEFM(I,L)=TEFM(I,L)+2.D0*INTGRL*(DM(I,J)+DM(J,I))
        TEFM(L,I)=TEFM(L,I)+2.D0*INTGRL*(DM(I,J)+DM(J,I))
        TEFM(I,I)=TEFM(I,I)-INTGRL*(DM(J,L)+DM(L,J))
        TEFM(L,I)=TEFM(L,I)-INTGRL*DM(I,J)
        TEFM(I,J)=TEFM(I,J)-INTGRL*DM(L,I)
        TEFM(L,J)=TEFM(L,J)-INTGRL*DM(I,I)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DM(J,I)
        TEFM(J,I)=TEFM(J,I)-INTGRL*DM(I,L)
        TEFM(J,L)=TEFM(J,L)-INTGRL*DM(I,I)
     ELSE IF ((I>J).AND.(J==K).AND.(K>L)) THEN
        TEFM(J,I)=TEFM(J,I)+2.D0*INTGRL*(DM(J,L)+DM(L,J))
        TEFM(I,J)=TEFM(I,J)+2.D0*INTGRL*(DM(J,L)+DM(L,J))
        TEFM(J,L)=TEFM(J,L)+2.D0*INTGRL*(DM(J,I)+DM(I,J))
        TEFM(L,J)=TEFM(L,J)+2.D0*INTGRL*(DM(J,I)+DM(I,J))
        TEFM(J,J)=TEFM(J,J)-INTGRL*(DM(I,L)+DM(L,I))
        TEFM(L,J)=TEFM(L,J)-INTGRL*DM(J,I)
        TEFM(J,I)=TEFM(J,I)-INTGRL*DM(L,J)
        TEFM(L,I)=TEFM(L,I)-INTGRL*DM(J,J)
        TEFM(J,L)=TEFM(J,L)-INTGRL*DM(I,J)
        TEFM(I,J)=TEFM(I,J)-INTGRL*DM(J,L)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DM(J,J)
     ELSE IF ((I>K).AND.(K>J).AND.(J==L)) THEN
        TEFM(J,I)=TEFM(J,I)+2.D0*INTGRL*(DM(J,K)+DM(K,J))
        TEFM(I,J)=TEFM(I,J)+2.D0*INTGRL*(DM(J,K)+DM(K,J))
        TEFM(J,K)=TEFM(J,K)+2.D0*INTGRL*(DM(J,I)+DM(I,J))
        TEFM(K,J)=TEFM(K,J)+2.D0*INTGRL*(DM(J,I)+DM(I,J))
        TEFM(J,J)=TEFM(J,J)-INTGRL*(DM(I,K)+DM(K,I))
        TEFM(K,J)=TEFM(K,J)-INTGRL*DM(J,I)
        TEFM(J,I)=TEFM(J,I)-INTGRL*DM(K,J)
        TEFM(K,I)=TEFM(K,I)-INTGRL*DM(J,J)
        TEFM(J,K)=TEFM(J,K)-INTGRL*DM(I,J) 
        TEFM(I,J)=TEFM(I,J)-INTGRL*DM(J,K)
        TEFM(I,K)=TEFM(I,K)-INTGRL*DM(J,J)
     ELSE IF ((I>J).AND.(I>K).AND.(K==L)) THEN
        TEFM(I,J)=TEFM(I,J)+2.D0*INTGRL*DM(K,K)
        TEFM(J,I)=TEFM(J,I)+2.D0*INTGRL*DM(K,K)
        TEFM(K,K)=TEFM(K,K)+2.D0*INTGRL*(DM(I,J)+DM(J,I))
        TEFM(K,I)=TEFM(K,I)-INTGRL*DM(K,J)
        TEFM(K,J)=TEFM(K,J)-INTGRL*DM(K,I)
        TEFM(I,K)=TEFM(I,K)-INTGRL*DM(J,K)
        TEFM(J,K)=TEFM(J,K)-INTGRL*DM(I,K)
     ELSE IF ((I==J).AND.(J>K).AND.(K>L)) THEN
        TEFM(K,L)=TEFM(K,L)+2.D0*INTGRL*DM(I,I)
        TEFM(L,K)=TEFM(L,K)+2.D0*INTGRL*DM(I,I)
        TEFM(I,I)=TEFM(I,I)+2.D0*INTGRL*(DM(K,L)+DM(L,K))
        TEFM(I,K)=TEFM(I,K)-INTGRL*DM(I,L)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DM(I,K)
        TEFM(K,I)=TEFM(K,I)-INTGRL*DM(L,I)
        TEFM(L,I)=TEFM(L,I)-INTGRL*DM(K,I)
! 4 distinct values for the 4 indices
     ELSE IF (    ((I>J).AND.(J>K).AND.(K>L)) &
              .OR.((I>K).AND.(K>J).AND.(J>L)) &
              .OR.((I>K).AND.(K>L).AND.(L>J))) THEN
        TEFM(I,J)=TEFM(I,J)+2.D0*INTGRL*(DM(K,L)+DM(L,K))
        TEFM(J,I)=TEFM(J,I)+2.D0*INTGRL*(DM(K,L)+DM(L,K))
        TEFM(K,L)=TEFM(K,L)+2.D0*INTGRL*(DM(I,J)+DM(J,I))
        TEFM(L,K)=TEFM(L,K)+2.D0*INTGRL*(DM(I,J)+DM(J,I))
        TEFM(K,I)=TEFM(K,I)-INTGRL*DM(L,J)
        TEFM(L,I)=TEFM(L,I)-INTGRL*DM(K,J)
        TEFM(K,J)=TEFM(K,J)-INTGRL*DM(L,I)
        TEFM(L,J)=TEFM(L,J)-INTGRL*DM(K,I)
        TEFM(I,K)=TEFM(I,K)-INTGRL*DM(J,L)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DM(J,K)
        TEFM(J,K)=TEFM(J,K)-INTGRL*DM(I,L)
        TEFM(J,L)=TEFM(J,L)-INTGRL*DM(I,K)
     END IF
  END DO
  IF (.NOT.DIRECT.AND.USEDISK) CLOSE(BIUNIT)
  PTEFM=PACK(TEFM,NBAST)
END SUBROUTINE BUILDTEFM_RHF

SUBROUTINE BUILDTEFM_RGHF(PTEFM,NBAST,PHI,PDM)
! Computation and assembly of the two-electron part of the Fock matrix associated to a given density matrix in the real general Hartree-Fock formalism, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
! Note: G(D)=J(D)-K(D), with J(D) the Coulomb term and K(D) the exchange term.
  USE scf_parameters ; USE basis_parameters ; USE integrals ; use matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2) :: PCM,PEM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  CALL BUILDCOULOMB_nonrelativistic(PCM,NBAST,PHI,PDM)
  CALL BUILDEXCHANGE_nonrelativistic(PEM,NBAST,PHI,PDM)
  PTEFM=PCM-PEM
END SUBROUTINE BUILDTEFM_RGHF

SUBROUTINE BUILDCOULOMB_relativistic(PCM,NBAST,PHI,PDM)
! Computation and assembly of the Coulomb term in the Fock matrix associated to a given density matrix, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals ; USE matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PCM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  INTEGER(smallint) :: I,J,K,L
  INTEGER :: N
  CHARACTER(2) :: CLASS
  DOUBLE COMPLEX :: INTGRL
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: CM,DM

  CM=(0.D0,0.D0)
  DM=UNPACK(PDM,NBAST)
  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) OPEN(LUNIT,access='STREAM')
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) OPEN(BIUNIT,access='STREAM')
  DO N=1,BINMBR
     IF (DIRECT) THEN
! the value of the bielectronic integral is computed "on the fly"
        IF (USEDISK) THEN
           READ(LUNIT)I,J,K,L
        ELSE
           I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
        END IF
        INTGRL=COULOMBVALUE(PHI(I),PHI(J),PHI(K),PHI(L))
     ELSE
        IF (SEMIDIRECT) THEN
! the value of the bielectronic integral is computed "on the fly", but using the precomputed values of the involved CGTO bielectronic integrals
           IF (USEDISK) THEN
              READ(LUNIT)I,J,K,L,CLASS
           ELSE
              I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
              CLASS=BITYPE(N)
           END IF
           INTGRL=COULOMBVALUE(PHI(I),PHI(J),PHI(K),PHI(L),CLASS)
        ELSE
           IF (USEDISK) THEN
! the value of the bielectronic integral is read on disk
              READ(BIUNIT)I,J,K,L,INTGRL
           ELSE
              I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
              INTGRL=CBIVALUES(N)
           END IF
        END IF
     END IF
     IF ((I/=J).AND.(I==K).AND.(J==L)) THEN
        CM(I,J)=CM(I,J)+INTGRL*DM(I,J)
     ELSE
        CM(I,J)=CM(I,J)+INTGRL*DM(K,L)
        CM(K,L)=CM(K,L)+INTGRL*DM(I,J)
     END IF
  END DO
  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) CLOSE(LUNIT)
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) CLOSE(BIUNIT)
  PCM=PACK(CM,NBAST)
END SUBROUTINE BUILDCOULOMB_relativistic

SUBROUTINE BUILDCOULOMB_nonrelativistic(PCM,NBAST,PHI,PDM)
! Computation and assembly of the Coulomb term in the Fock matrix associated to a given density matrix, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals ; use matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PCM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: CM,DM
  INTEGER(smallint) :: I,J,K,L
  INTEGER :: N
  DOUBLE PRECISION :: INTGRL

  CM=0.D0
  DM=UNPACK(PDM,NBAST)
#define ACTION(I,J,K,L) CM(I,J)=CM(I,J)+INTGRL*DM(K,L)
#include "forall.f90"
#undef ACTION
  PCM=PACK(CM,NBAST)
END SUBROUTINE BUILDCOULOMB_nonrelativistic

SUBROUTINE BUILDEXCHANGE_relativistic(PEM,NBAST,PHI,PDM)
! Computation and assembly of the exchange term in the Fock matrix associated to a given density matrix, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals ; USE matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PEM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  INTEGER(smallint) :: I,J,K,L
  INTEGER :: N
  CHARACTER(2) :: CLASS
  DOUBLE COMPLEX :: INTGRL
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: EM,DM

  EM=(0.D0,0.D0)
  DM=UNPACK(PDM,NBAST)
  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) OPEN(LUNIT,access='STREAM')
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) OPEN(BIUNIT,access='STREAM')
  DO N=1,BINMBR
     IF (DIRECT) THEN
! the value of the bielectronic integral is computed "on the fly"
        IF (USEDISK) THEN
           READ(LUNIT)I,J,K,L
        ELSE
           I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
        END IF
        INTGRL=COULOMBVALUE(PHI(I),PHI(J),PHI(K),PHI(L))
     ELSE
        IF (SEMIDIRECT) THEN
! the value of the bielectronic integral is computed "on the fly", but using the precomputed values of the involved CGTO bielectronic integrals
           IF (USEDISK) THEN
              READ(LUNIT)I,J,K,L,CLASS
           ELSE
              I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
              CLASS=BITYPE(N)
           END IF
           INTGRL=COULOMBVALUE(PHI(I),PHI(J),PHI(K),PHI(L),CLASS)
        ELSE
           IF (USEDISK) THEN
! the value of the bielectronic integral is read on disk
              READ(BIUNIT)I,J,K,L,INTGRL
           ELSE
              I=BILIST(N,1) ; J=BILIST(N,2) ; K=BILIST(N,3) ; L=BILIST(N,4)
              INTGRL=CBIVALUES(N)
           END IF
        END IF
     END IF
     IF ((I/=J).AND.(I==K).AND.(J==L)) THEN
        EM(I,J)=EM(I,J)+INTGRL*DM(J,I)
     ELSE
        EM(I,L)=EM(I,L)+INTGRL*DM(J,K)
        EM(K,J)=EM(K,J)+INTGRL*DM(L,I)
     END IF
  END DO
  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) CLOSE(LUNIT)
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) CLOSE(BIUNIT)
  PEM=PACK(EM,NBAST)
END SUBROUTINE BUILDEXCHANGE_relativistic

SUBROUTINE BUILDEXCHANGE_nonrelativistic(PEM,NBAST,PHI,PDM)
! Computation and assembly of the exchange term in the Fock matrix associated to a given density matrix, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals ; USE matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PEM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: EM,DM
  INTEGER(smallint) :: I,J,K,L
  INTEGER :: N
  DOUBLE PRECISION :: INTGRL

  EM=0.D0
  DM=UNPACK(PDM,NBAST)
#define ACTION(I,J,K,L) EM(I,K)=EM(I,K)+INTGRL*DM(L,J)
#include "forall.f90"
#undef ACTION
  PEM=PACK(EM,NBAST)
END SUBROUTINE BUILDEXCHANGE_nonrelativistic

SUBROUTINE BUILDSAMCM(PSAMCM,PHI,NBAST,NBAS,COMPONENT)
! Computation and assembly of the matrix associated to one of the three components of the spin angular momentum operator S=-i/4\alpha^\alpha (only the upper triangular part of the matrix is stored in packed format).
  USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PSAMCM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS
  INTEGER :: COMPONENT

  INTEGER :: I,J,K,L,M
  DOUBLE COMPLEX :: VALUE

  PSAMCM=(0.D0,0.D0)
  SELECT CASE (COMPONENT)
  CASE (1)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))    &
 &                        *OVERLAPVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))    &
 &                        *OVERLAPVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L))
           END DO
        END DO
        PSAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))    &
 &                        *OVERLAPVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))    &
 &                        *OVERLAPVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L))   
           END DO
        END DO
        PSAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  CASE (2)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE-.5D0*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))                 &
 &                        *DCMPLX(0.D0,OVERLAPVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L)))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))                 &
 &                        *DCMPLX(0.D0,OVERLAPVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L)))
           END DO
        END DO
        PSAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE-.5D0*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))                 &
 &                        *DCMPLX(0.D0,OVERLAPVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L)))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))                 &
 &                        *DCMPLX(0.D0,OVERLAPVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L)))
           END DO
        END DO
        PSAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  CASE (3)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                         &
 &                           *.5D0*(-1.D0)**K*OVERLAPVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L))
              END DO
           END DO
        END DO
        PSAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                         &
 &                           *.5D0*(-1.D0)**K*OVERLAPVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L))
              END DO
           END DO
        END DO
        PSAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  END SELECT
END SUBROUTINE BUILDSAMCM

SUBROUTINE BUILDOAMCM(POAMCM,PHI,NBAST,NBAS,COMPONENT)
! Computation and assembly of the matrix associated to one of the three components of the orbital angular momentum operator L=x^p (only the upper triangular part of the matrix is stored in packed format).
  USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POAMCM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS
  INTEGER :: COMPONENT

  INTEGER :: I,J,K,L,M
  DOUBLE COMPLEX :: VALUE

  POAMCM=(0.D0,0.D0)
  SELECT CASE (COMPONENT)
  CASE (1)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                          &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),3,2)   &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),2,3))
              END DO
           END DO
        END DO
        POAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                          &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),3,2)   &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),2,3))
              END DO
           END DO
        END DO
        POAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  CASE (2)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                          &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),1,3)   &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),3,1))
              END DO
           END DO
        END DO
        POAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                          &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),1,3)   &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),3,1))
              END DO
           END DO
        END DO
        POAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  CASE (3)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                          &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),2,1)   &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),1,2))
              END DO
           END DO
        END DO
        POAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                          &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),2,1)   &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),1,2))
              END DO
           END DO
        END DO
        POAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  END SELECT
END SUBROUTINE BUILDOAMCM

SUBROUTINE BUILDTAMCM(PTAMCM,PHI,NBAST,NBAS,COMPONENT)
! Computation and assembly of the matrix associated to one of the three components of the total angular momentum operator J=L+S (only the upper triangular part of the matrix is stored in packed format).
  USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTAMCM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS
  INTEGER :: COMPONENT

  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2) :: PSAMCM,POAMCM

  CALL BUILDSAMCM(PSAMCM,PHI,NBAST,NBAS,COMPONENT)
  CALL BUILDOAMCM(POAMCM,PHI,NBAST,NBAS,COMPONENT)
  PTAMCM=PSAMCM+POAMCM
END SUBROUTINE BUILDTAMCM

SUBROUTINE BUILDA(PA,NBAST,NBO,PHI,EIG,EIGVEC)
! Computation and assembly of the matrix A defined by equation (8c) in the reference below.
! Reference: G. B. Bacskay, A quadratically convergent Hartree-Fock (QC-SCF) method. Application to closed shell systems, Chem. Phys., 61(3), 385-404, 1981.
  USE scf_parameters ; USE basis_parameters ; USE integrals ; USE matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST,NBO
  DOUBLE PRECISION,DIMENSION(NBO*(NBAST-NBO)*(NBO*(NBAST-NBO)+1)/2),INTENT(OUT) :: PA
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(IN) :: EIG
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(IN) :: EIGVEC

  INTEGER(smallint) :: I,J,K,L
  INTEGER :: N,NBV,OI,OJ,VA,VB,IA,JB
  DOUBLE PRECISION :: INTGRL
  DOUBLE PRECISION,DIMENSION(NBO*(NBAST-NBO),NBO*(NBAST-NBO)) :: A

  NBV=NBAST-NBO
  A=0.D0
  DO OI=1,NBO
     DO VA=1,NBV
        IA=NBV*(OI-1)+VA
        A(IA,IA)=EIG(NBO+VA)-EIG(OI)
     END DO
  END DO
#define ACTION(I,J,K,L) \
  DO OI=1,NBO ;\
     DO VA=1,NBV ;\
        DO OJ=1,NBO ;\
           DO VB=1,NBV ;\
              IA=NBV*(OI-1)+VA ; JB=NBV*(OJ-1)+VB ;\
              A(IA,JB)=A(IA,JB)+INTGRL*EIGVEC(I,NBO+VA)*EIGVEC(K,OJ)*(EIGVEC(J,OI)*EIGVEC(L,NBO+VB)-EIGVEC(J,NBO+VB)*EIGVEC(L,OI)) ;\
           END DO ;\
        END DO ;\
     END DO ;\
  END DO
#include "forall.f90"
#undef ACTION
  PA=PACK(A,(NBAST-NBO)*NBO)
END SUBROUTINE BUILDA

SUBROUTINE BUILDB(PB,NBAST,NBO,PHI,EIG,EIGVEC)
! Computation and assembly of the matrix B defined by equation (8d) in the reference below.
! Reference: G. B. Bacskay, A quadratically convergent Hartree-Fock (QC-SCF) method. Application to closed shell systems, Chem. Phys., 61(3), 385-404, 1981.
  USE scf_parameters ; USE basis_parameters ; USE integrals ; USE matrix_tools
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST,NBO
  DOUBLE PRECISION,DIMENSION(NBO*(NBAST-NBO)*(NBO*(NBAST-NBO)+1)/2),INTENT(OUT) :: PB
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(IN) :: EIG
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(IN) :: EIGVEC

  INTEGER(smallint) :: I,J,K,L
  INTEGER :: N,NBV,OI,OJ,VA,VB,IA,JB
  DOUBLE PRECISION :: INTGRL
  DOUBLE PRECISION,DIMENSION(NBO*(NBAST-NBO),NBO*(NBAST-NBO)) :: B

  NBV=NBAST-NBO
  B=0.D0
#define ACTION(I,J,K,L) \
  DO OI=1,NBO ;\
     DO VA=1,NBV ;\
        DO OJ=1,NBO ;\
           DO VB=1,NBV ;\
              IA=NBV*(OI-1)+VA ; JB=NBV*(OJ-1)+VB ;\
              B(IA,VB)=B(IA,JB)+INTGRL*EIGVEC(I,OI)*EIGVEC(K,OJ)*(EIGVEC(J,NBO+VA)*EIGVEC(L,NBO+VB)-EIGVEC(J,NBO+VB)*EIGVEC(L,NBO+VA)) ;\
           END DO ;\
        END DO ;\
     END DO ;\
  END DO
#include "forall.f90"
#undef ACTION
  PB=PACK(B,NBO*(NBAST-NBO))
END SUBROUTINE BUILDB

MODULE matrices
INTERFACE FORMDM
  SUBROUTINE FORMDM_relativistic(PDM,EIGVEC,NBAST,LOON,HOON)
    INTEGER,INTENT(IN) :: NBAST,LOON,HOON
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PDM
    DOUBLE COMPLEX,DIMENSION(NBAST,NBAST),INTENT(IN) :: EIGVEC
  END SUBROUTINE

  SUBROUTINE FORMDM_nonrelativistic(PDM,EIGVEC,NBAST,LOON,HOON)
    INTEGER,INTENT(IN) :: NBAST,LOON,HOON
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PDM
    DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(IN) :: EIGVEC
  END SUBROUTINE
END INTERFACE

INTERFACE BUILDOM
  SUBROUTINE BUILDOM_relativistic(POM,PHI,NBAST,NBAS)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POM
    TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
    INTEGER,DIMENSION(2),INTENT(IN) :: NBAS
  END SUBROUTINE

  SUBROUTINE BUILDOM_nonrelativistic(POM,PHI,NBAST)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  END SUBROUTINE
END INTERFACE

INTERFACE BUILDKPFM
  SUBROUTINE BUILDKPFM_nonrelativistic(PKPFM,PHI,NBAST)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PKPFM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  END SUBROUTINE
END INTERFACE

INTERFACE BUILDOEFM
  SUBROUTINE BUILDOEFM_relativistic(POEFM,PHI,NBAST,NBAS)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POEFM
    TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
    INTEGER,DIMENSION(2),INTENT(IN) :: NBAS
  END SUBROUTINE

  SUBROUTINE BUILDOEFM_nonrelativistic(POEFM,PHI,NBAST)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POEFM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  END SUBROUTINE
END INTERFACE

INTERFACE BUILDTEFM
  SUBROUTINE BUILDTEFM_relativistic(PTEFM,NBAST,PHI,PDM)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
    TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM
  END SUBROUTINE

  SUBROUTINE BUILDTEFM_RHF(PTEFM,NBAST,PHI,PDM)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM
  END SUBROUTINE
END INTERFACE

INTERFACE BUILDCOULOMB
  SUBROUTINE BUILDCOULOMB_relativistic(PCM,NBAST,PHI,PDM)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PCM
    TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM
  END SUBROUTINE

  SUBROUTINE BUILDCOULOMB_nonrelativistic(PCM,NBAST,PHI,PDM)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PCM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM
  END SUBROUTINE
END INTERFACE

INTERFACE BUILDEXCHANGE
  SUBROUTINE BUILDEXCHANGE_relativistic(PEM,NBAST,PHI,PDM)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PEM
    TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
    DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM
  END SUBROUTINE

  SUBROUTINE BUILDEXCHANGE_nonrelativistic(PEM,NBAST,PHI,PDM)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PEM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM
  END SUBROUTINE
END INTERFACE
END MODULE
