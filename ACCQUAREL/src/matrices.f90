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
END SUBROUTINE

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
END SUBROUTINE

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
END SUBROUTINE

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
END SUBROUTINE

SUBROUTINE BUILDOM_nonrelativistic(POM,PHI,NBAST)
! Computation and assembly of the overlap matrix between basis functions, i.e. the Gram matrix of the basis with respacet to the $L^2(\mathbb{R}^3)$ inner product (only the upper triangular part of the matrix is stored in packed format).
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
END SUBROUTINE

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
END SUBROUTINE

SUBROUTINE BUILDOEFM_relativistic(POEFM,PHI,NBAST,NBAS)
! Computation and assembly of the monoelectronic part of the Fock matrix (only the upper triangular part of the matrix is stored in packed form)
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: POEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS

  INTEGER :: I,J,K,L,M,N
  DOUBLE PRECISION :: AC
  DOUBLE COMPLEX :: TMP,VALUE

  AC=SCALING_FACTOR*C

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
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=1,NBAS(1)
        VALUE=(0.D0,0.D0)
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE-AC*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(1,L))                   &
 &                        *DCMPLX(0.D0,DERIVVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(1,L),3))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE-AC*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))               &
 &                        *DCMPLX(DERIVVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L),2),  &
 &                                DERIVVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L),1))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE-AC*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))                &
 &                        *DCMPLX(-DERIVVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L),2),  &
 &                                DERIVVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L),1))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE+AC*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(2,L))                   &
 &                        *DCMPLX(0.D0,DERIVVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(2,L),3))
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
                 TMP=2.D0*AC*AC*OVERLAPVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L))
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
END SUBROUTINE

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
END SUBROUTINE

SUBROUTINE BUILDTEFM_relativistic(PTEFM,NBAST,PHI,PDM)
! Computation and assembly of the bielectronic part of the Fock matrix associated to a given density matrix using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  INTEGER :: I,J,K,L,N
  CHARACTER(2) :: CLASS
  DOUBLE COMPLEX :: INTGRL
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: TEFM,DM

  TEFM=(0.D0,0.D0)
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        DM(I,J)=PDM(N)
        IF (I/=J) DM(J,I)=CONJG(DM(I,J))
     END DO
  END DO

  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) OPEN(LUNIT,form='UNFORMATTED')
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) OPEN(BIUNIT,form='UNFORMATTED')
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
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        PTEFM(N)=TEFM(I,J)
     END DO
  END DO
END SUBROUTINE

SUBROUTINE BUILDTEFM_RHF(PTEFM,NBAST,PHI,PDM)
! Computation and assembly of the two-electron part of the Fock matrix associated to a given density matrix in the restricted closed-shell Hartree-Fock formalism, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
! Note: G(D)=2J(D)-K(D), with J(D) the Coulomb term and K(D) the exchange term.
  USE scf_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: TEFM,DM
  INTEGER :: I,J,K,L,N
  DOUBLE PRECISION :: INTGRL

  TEFM=0.D0
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        DM(I,J)=PDM(N)
        DM(J,I)=DM(I,J)
     END DO
  END DO
  IF (.NOT.DIRECT.AND.USEDISK) OPEN(BIUNIT,form='UNFORMATTED')
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
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        PTEFM(N)=TEFM(I,J)
     END DO
  END DO
END SUBROUTINE

SUBROUTINE BUILDTEFM_UHF(PTEFM,NBAST,PHI,PDMA,PDMB)
! Computation and assembly of the two-electron part of one of the two Fock matrices in the unrestricted open-shell Hartree-Fock formalism, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
! Note: G(D_a)=J(D_a)-K(D_a)+J(D_b), with J(D_a) the Coulomb term and K(D_a) the exchange term associated to D_a, and J(D_b) the Coulomb term associated to D_b.
  USE scf_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDMA,PDMB

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: TEFM,DM,DMA,DMB
  INTEGER :: I,J,K,L,N
  DOUBLE PRECISION :: INTGRL

  TEFM=0.D0
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        DMA(I,J)=PDMA(N) ; DMB(I,J)=PDMB(N)
        DMA(J,I)=DMA(I,J) ; DMB(J,I)=DMB(I,J)
     END DO
  END DO
  DM=DMA+DMB
  IF (.NOT.DIRECT.AND.USEDISK) OPEN(BIUNIT,form='UNFORMATTED')
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
        TEFM(I,I)=TEFM(I,I)+INTGRL*DMB(I,I)
! 2 distinct values for the 4 indices C
     ELSE IF ((I>J).AND.(J==K).AND.(K==L)) THEN
        TEFM(I,J)=TEFM(I,J)+INTGRL*DMB(J,J)
        TEFM(J,I)=TEFM(J,I)+INTGRL*DMB(J,J)
        TEFM(J,J)=TEFM(J,J)+INTGRL*(DMB(I,J)+DMB(J,I))
     ELSE IF ((I==J).AND.(J==K).AND.(K>L)) THEN
        TEFM(L,I)=TEFM(L,I)+INTGRL*DMB(I,I)
        TEFM(I,L)=TEFM(I,L)+INTGRL*DMB(I,I)
        TEFM(I,I)=TEFM(I,I)+INTGRL*(DMB(L,I)+DMB(I,L))
     ELSE IF ((I==J).AND.(J>K).AND.(K==L)) THEN
        TEFM(I,I)=TEFM(I,I)+INTGRL*DM(K,K)
        TEFM(K,K)=TEFM(K,K)+INTGRL*DM(I,I)
        TEFM(I,K)=TEFM(I,K)-INTGRL*DMA(I,K)
        TEFM(K,I)=TEFM(K,I)-INTGRL*DMA(K,I)
     ELSE IF ((I==K).AND.(K>J).AND.(J==L)) THEN
        TEFM(I,J)=TEFM(I,J)+INTGRL*DM(I,J)
        TEFM(J,I)=TEFM(J,I)+INTGRL*DM(J,I)
        TEFM(I,J)=TEFM(I,J)+INTGRL*DMA(J,I)
        TEFM(J,I)=TEFM(J,I)+INTGRL*DMA(I,J)
        TEFM(I,I)=TEFM(I,I)-INTGRL*DMA(J,J)
        TEFM(J,J)=TEFM(J,J)-INTGRL*DMA(I,I)
! 3 distinct values for the 4 indices
     ELSE IF ((I==K).AND.(K>J).AND.(J>L)) THEN
        TEFM(I,J)=TEFM(I,J)+INTGRL*(DM(I,L)+DM(L,I))
        TEFM(J,I)=TEFM(J,I)+INTGRL*(DM(I,L)+DM(L,I))
        TEFM(I,L)=TEFM(I,L)+INTGRL*(DM(I,J)+DM(J,I))
        TEFM(L,I)=TEFM(L,I)+INTGRL*(DM(I,J)+DM(J,I))
        TEFM(I,I)=TEFM(I,I)-INTGRL*(DMA(J,L)+DMA(L,J))
        TEFM(L,I)=TEFM(L,I)-INTGRL*DMA(I,J)
        TEFM(I,J)=TEFM(I,J)-INTGRL*DMA(L,I)
        TEFM(L,J)=TEFM(L,J)-INTGRL*DMA(I,I)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DMA(J,I)
        TEFM(J,I)=TEFM(J,I)-INTGRL*DMA(I,L)
        TEFM(J,L)=TEFM(J,L)-INTGRL*DMA(I,I)
     ELSE IF ((I>J).AND.(J==K).AND.(K>L)) THEN
        TEFM(J,I)=TEFM(J,I)+INTGRL*(DM(J,L)+DM(L,J))
        TEFM(I,J)=TEFM(I,J)+INTGRL*(DM(J,L)+DM(L,J))
        TEFM(J,L)=TEFM(J,L)+INTGRL*(DM(J,I)+DM(I,J))
        TEFM(L,J)=TEFM(L,J)+INTGRL*(DM(J,I)+DM(I,J))
        TEFM(J,J)=TEFM(J,J)-INTGRL*(DMA(I,L)+DMA(L,I))
        TEFM(L,J)=TEFM(L,J)-INTGRL*DMA(J,I)
        TEFM(J,I)=TEFM(J,I)-INTGRL*DMA(L,J)
        TEFM(L,I)=TEFM(L,I)-INTGRL*DMA(J,J)
        TEFM(J,L)=TEFM(J,L)-INTGRL*DMA(I,J)
        TEFM(I,J)=TEFM(I,J)-INTGRL*DMA(J,L)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DMA(J,J)
     ELSE IF ((I>K).AND.(K>J).AND.(J==L)) THEN
        TEFM(J,I)=TEFM(J,I)+INTGRL*(DM(J,K)+DM(K,J))
        TEFM(I,J)=TEFM(I,J)+INTGRL*(DM(J,K)+DM(K,J))
        TEFM(J,K)=TEFM(J,K)+INTGRL*(DM(J,I)+DM(I,J))
        TEFM(K,J)=TEFM(K,J)+INTGRL*(DM(J,I)+DM(I,J))
        TEFM(J,J)=TEFM(J,J)-INTGRL*(DMA(I,K)+DMA(K,I))
        TEFM(K,J)=TEFM(K,J)-INTGRL*DMA(J,I)
        TEFM(J,I)=TEFM(J,I)-INTGRL*DMA(K,J)
        TEFM(K,I)=TEFM(K,I)-INTGRL*DMA(J,J)
        TEFM(J,K)=TEFM(J,K)-INTGRL*DMA(I,J)
        TEFM(I,J)=TEFM(I,J)-INTGRL*DMA(J,K)
        TEFM(I,K)=TEFM(I,K)-INTGRL*DMA(J,J)
     ELSE IF ((I>J).AND.(I>K).AND.(K==L)) THEN
        TEFM(I,J)=TEFM(I,J)+INTGRL*DM(K,K)
        TEFM(J,I)=TEFM(J,I)+INTGRL*DM(K,K)
        TEFM(K,K)=TEFM(K,K)+INTGRL*(DM(I,J)+DM(J,I))
        TEFM(K,I)=TEFM(K,I)-INTGRL*DMA(K,J)
        TEFM(K,J)=TEFM(K,J)-INTGRL*DMA(K,I)
        TEFM(I,K)=TEFM(I,K)-INTGRL*DMA(J,K)
        TEFM(J,K)=TEFM(J,K)-INTGRL*DMA(I,K)
     ELSE IF ((I==J).AND.(J>K).AND.(K>L)) THEN
        TEFM(K,L)=TEFM(K,L)+INTGRL*DM(I,I)
        TEFM(L,K)=TEFM(L,K)+INTGRL*DM(I,I)
        TEFM(I,I)=TEFM(I,I)+INTGRL*(DM(K,L)+DM(L,K))
        TEFM(I,K)=TEFM(I,K)-INTGRL*DMA(I,L)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DMA(I,K)
        TEFM(K,I)=TEFM(K,I)-INTGRL*DMA(L,I)
        TEFM(L,I)=TEFM(L,I)-INTGRL*DMA(K,I)
! 4 distinct values for the 4 indices
     ELSE IF (    ((I>J).AND.(J>K).AND.(K>L)) &
              .OR.((I>K).AND.(K>J).AND.(J>L)) &
              .OR.((I>K).AND.(K>L).AND.(L>J))) THEN
        TEFM(I,J)=TEFM(I,J)+INTGRL*(DM(K,L)+DM(L,K))
        TEFM(J,I)=TEFM(J,I)+INTGRL*(DM(K,L)+DM(L,K))
        TEFM(K,L)=TEFM(K,L)+INTGRL*(DM(I,J)+DM(J,I))
        TEFM(L,K)=TEFM(L,K)+INTGRL*(DM(I,J)+DM(J,I))
        TEFM(K,I)=TEFM(K,I)-INTGRL*DMA(L,J)
        TEFM(L,I)=TEFM(L,I)-INTGRL*DMA(K,J)
        TEFM(K,J)=TEFM(K,J)-INTGRL*DMA(L,I)
        TEFM(L,J)=TEFM(L,J)-INTGRL*DMA(K,I)
        TEFM(I,K)=TEFM(I,K)-INTGRL*DMA(J,L)
        TEFM(I,L)=TEFM(I,L)-INTGRL*DMA(J,K)
        TEFM(J,K)=TEFM(J,K)-INTGRL*DMA(I,L)
        TEFM(J,L)=TEFM(J,L)-INTGRL*DMA(I,K)
     END IF
  END DO
  IF (.NOT.DIRECT.AND.USEDISK) CLOSE(BIUNIT)
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        PTEFM(N)=TEFM(I,J)
     END DO
  END DO
END SUBROUTINE

SUBROUTINE BUILDCOULOMB_relativistic(PCM,NBAST,PHI,PDM)
! Computation and assembly of the Coulomb term in the Fock matrix associated to a given density matrix, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PCM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  INTEGER :: I,J,K,L,N
  CHARACTER(2) :: CLASS
  DOUBLE COMPLEX :: INTGRL
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: CM,DM

  CM=(0.D0,0.D0)
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        DM(I,J)=PDM(N)
        IF (I/=J) DM(J,I)=CONJG(DM(I,J))
     END DO
  END DO

  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) OPEN(LUNIT,form='UNFORMATTED')
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) OPEN(BIUNIT,form='UNFORMATTED')
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
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        PCM(N)=CM(I,J)
     END DO
  END DO
END SUBROUTINE

SUBROUTINE BUILDCOULOMB_nonrelativistic(PCM,NBAST,PHI,PDM)
! Computation and assembly of the Coulomb term in the Fock matrix associated to a given density matrix, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PCM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: CM,DM
  INTEGER :: I,J,K,L,N
  DOUBLE PRECISION :: INTGRL

  CM=0.D0
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        DM(I,J)=PDM(N)
        DM(J,I)=DM(I,J)
     END DO
  END DO
  IF (.NOT.DIRECT.AND.USEDISK) OPEN(BIUNIT,form='UNFORMATTED')
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
        CM(I,I)=CM(I,I)+INTGRL*DM(I,I)
! 2 distinct values for the 4 indices
     ELSE IF ((I>J).AND.(J==K).AND.(K==L)) THEN
        CM(I,J)=CM(I,J)+INTGRL*DM(J,J)
        CM(J,I)=CM(J,I)+INTGRL*DM(J,J)
        CM(J,J)=CM(J,J)+INTGRL*(DM(I,J)+DM(J,I))
     ELSE IF ((I==J).AND.(J==K).AND.(K>L)) THEN
        CM(L,I)=CM(L,I)+INTGRL*DM(I,I)
        CM(I,L)=CM(I,L)+INTGRL*DM(I,I)
        CM(I,I)=CM(I,I)+INTGRL*(DM(L,I)+DM(I,L))
     ELSE IF ((I==J).AND.(J>K).AND.(K==L)) THEN
        CM(I,I)=CM(I,I)+INTGRL*DM(K,K)
        CM(K,K)=CM(K,K)+INTGRL*DM(I,I)
     ELSE IF ((I==K).AND.(K>J).AND.(J==L)) THEN
        CM(I,J)=CM(I,J)+INTGRL*(DM(I,J)+DM(J,I))
        CM(J,I)=CM(J,I)+INTGRL*(DM(J,I)+DM(I,J))
! 3 distinct values for the 4 indices
     ELSE IF ((I==K).AND.(K>J).AND.(J>L)) THEN
        CM(I,J)=CM(I,J)+INTGRL*(DM(I,L)+DM(L,I))
        CM(J,I)=CM(J,I)+INTGRL*(DM(I,L)+DM(L,I))
        CM(I,L)=CM(I,L)+INTGRL*(DM(I,J)+DM(J,I))
        CM(L,I)=CM(L,I)+INTGRL*(DM(I,J)+DM(J,I))
     ELSE IF ((I>J).AND.(J==K).AND.(K>L)) THEN
        CM(J,I)=CM(J,I)+INTGRL*(DM(J,L)+DM(L,J))
        CM(I,J)=CM(I,J)+INTGRL*(DM(J,L)+DM(L,J))
        CM(J,L)=CM(J,L)+INTGRL*(DM(J,I)+DM(I,J))
        CM(L,J)=CM(L,J)+INTGRL*(DM(J,I)+DM(I,J))
     ELSE IF ((I>K).AND.(K>J).AND.(J==L)) THEN
        CM(J,I)=CM(J,I)+INTGRL*(DM(J,K)+DM(K,J))
        CM(I,J)=CM(I,J)+INTGRL*(DM(J,K)+DM(K,J))
        CM(J,K)=CM(J,K)+INTGRL*(DM(J,I)+DM(I,J))
        CM(K,J)=CM(K,J)+INTGRL*(DM(J,I)+DM(I,J))
     ELSE IF ((I>J).AND.(I>K).AND.(K==L)) THEN
        CM(I,J)=CM(I,J)+INTGRL*DM(K,K)
        CM(J,I)=CM(J,I)+INTGRL*DM(K,K)
        CM(K,K)=CM(K,K)+INTGRL*(DM(I,J)+DM(J,I))
     ELSE IF ((I==J).AND.(J>K).AND.(K>L)) THEN
        CM(K,L)=CM(K,L)+INTGRL*DM(I,I)
        CM(L,K)=CM(L,K)+INTGRL*DM(I,I)
        CM(I,I)=CM(I,I)+INTGRL*(DM(K,L)+DM(L,K))
! 4 distinct values for the 4 indices
     ELSE IF (    ((I>J).AND.(J>K).AND.(K>L)) &
              .OR.((I>K).AND.(K>J).AND.(J>L)) &
              .OR.((I>K).AND.(K>L).AND.(L>J))) THEN
        CM(I,J)=CM(I,J)+INTGRL*(DM(K,L)+DM(L,K))
        CM(J,I)=CM(J,I)+INTGRL*(DM(K,L)+DM(L,K))
        CM(K,L)=CM(K,L)+INTGRL*(DM(I,J)+DM(J,I))
        CM(L,K)=CM(L,K)+INTGRL*(DM(I,J)+DM(J,I))
     END IF
  END DO
  IF (.NOT.DIRECT.AND.USEDISK) CLOSE(BIUNIT)
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        PCM(N)=CM(I,J)
     END DO
  END DO
END SUBROUTINE

SUBROUTINE BUILDEXCHANGE_relativistic(PEM,NBAST,PHI,PDM)
! Computation and assembly of the exchange term in the Fock matrix associated to a given density matrix, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PEM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  INTEGER :: I,J,K,L,N
  CHARACTER(2) :: CLASS
  DOUBLE COMPLEX :: INTGRL
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: EM,DM

  EM=(0.D0,0.D0)
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        DM(I,J)=PDM(N)
        IF (I/=J) DM(J,I)=CONJG(DM(I,J))
     END DO
  END DO

  IF ((DIRECT.OR.SEMIDIRECT).AND.USEDISK) OPEN(LUNIT,form='UNFORMATTED')
  IF ((.NOT.DIRECT.AND..NOT.SEMIDIRECT).AND.USEDISK) OPEN(BIUNIT,form='UNFORMATTED')
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
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        PEM(N)=EM(I,J)
     END DO
  END DO
END SUBROUTINE

SUBROUTINE BUILDEXCHANGE_nonrelativistic(PEM,NBAST,PHI,PDM)
! Computation and assembly of the exchange term in the Fock matrix associated to a given density matrix, using a list of the nonzero integrals (only the upper triangular part of the matrix is stored in packed format).
  USE scf_parameters ; USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PEM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDM

  DOUBLE PRECISION,DIMENSION(NBAST,NBAST) :: EM,DM
  INTEGER :: I,J,K,L,N
  DOUBLE PRECISION :: INTGRL

  EM=0.D0
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        DM(I,J)=PDM(N)
        DM(J,I)=DM(I,J)
     END DO
  END DO
  IF (.NOT.DIRECT.AND.USEDISK) OPEN(BIUNIT,form='UNFORMATTED')
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
        EM(I,I)=EM(I,I)+INTGRL*DM(I,I)
! 2 distinct values for the 4 indices
     ELSE IF ((I>J).AND.(J==K).AND.(K==L)) THEN
        EM(I,J)=EM(I,J)+INTGRL*DM(J,J)
        EM(J,I)=EM(J,I)+INTGRL*DM(J,J)
        EM(J,J)=EM(J,J)+INTGRL*(DM(I,J)+DM(J,I))
     ELSE IF ((I==J).AND.(J==K).AND.(K>L)) THEN
        EM(L,I)=EM(L,I)+INTGRL*DM(I,I)
        EM(I,L)=EM(I,L)+INTGRL*DM(I,I)
        EM(I,I)=EM(I,I)+INTGRL*(DM(L,I)+DM(I,L))
     ELSE IF ((I==J).AND.(J>K).AND.(K==L)) THEN
        EM(I,K)=EM(I,K)+INTGRL*DM(I,K)
        EM(K,I)=EM(K,I)+INTGRL*DM(K,I)
     ELSE IF ((I==K).AND.(K>J).AND.(J==L)) THEN
        EM(I,I)=EM(I,I)+INTGRL*DM(J,J)
        EM(J,J)=EM(J,J)+INTGRL*DM(I,I)
        EM(I,J)=EM(I,J)+INTGRL*DM(J,I)
        EM(J,I)=EM(J,I)+INTGRL*DM(I,J)
! 3 distinct values for the 4 indices
     ELSE IF ((I==K).AND.(K>J).AND.(J>L)) THEN
        EM(I,I)=EM(I,I)+INTGRL*(DM(J,L)+DM(L,J))
        EM(L,I)=EM(L,I)+INTGRL*DM(I,J)
        EM(I,J)=EM(I,J)+INTGRL*DM(L,I)
        EM(L,J)=EM(L,J)+INTGRL*DM(I,I)
        EM(I,L)=EM(I,L)+INTGRL*DM(J,I)
        EM(J,I)=EM(J,I)+INTGRL*DM(I,L)
        EM(J,L)=EM(J,L)+INTGRL*DM(I,I)
     ELSE IF ((I>J).AND.(J==K).AND.(K>L)) THEN
        EM(J,J)=EM(J,J)+INTGRL*(DM(I,L)+DM(L,I))
        EM(L,J)=EM(L,J)+INTGRL*DM(J,I)
        EM(J,I)=EM(J,I)+INTGRL*DM(L,J)
        EM(L,I)=EM(L,I)+INTGRL*DM(J,J)
        EM(J,L)=EM(J,L)+INTGRL*DM(I,J)
        EM(I,J)=EM(I,J)+INTGRL*DM(J,L)
        EM(I,L)=EM(I,L)+INTGRL*DM(J,J)
     ELSE IF ((I>K).AND.(K>J).AND.(J==L)) THEN
        EM(J,J)=EM(J,J)+INTGRL*(DM(I,K)+DM(K,I))
        EM(K,J)=EM(K,J)+INTGRL*DM(J,I)
        EM(J,I)=EM(J,I)+INTGRL*DM(K,J)
        EM(K,I)=EM(K,I)+INTGRL*DM(J,J)
        EM(J,K)=EM(J,K)+INTGRL*DM(I,J)
        EM(I,J)=EM(I,J)+INTGRL*DM(J,K)
        EM(I,K)=EM(I,K)+INTGRL*DM(J,J)
     ELSE IF ((I>J).AND.(I>K).AND.(K==L)) THEN
        EM(K,I)=EM(K,I)+INTGRL*DM(K,J)
        EM(K,J)=EM(K,J)+INTGRL*DM(K,I)
        EM(I,K)=EM(I,K)+INTGRL*DM(J,K)
        EM(J,K)=EM(J,K)+INTGRL*DM(I,K)
     ELSE IF ((I==J).AND.(J>K).AND.(K>L)) THEN
        EM(I,K)=EM(I,K)+INTGRL*DM(I,L)
        EM(I,L)=EM(I,L)+INTGRL*DM(I,K)
        EM(K,I)=EM(K,I)+INTGRL*DM(L,I)
        EM(L,I)=EM(L,I)+INTGRL*DM(K,I)
! 4 distinct values for the 4 indices
     ELSE IF (    ((I>J).AND.(J>K).AND.(K>L)) &
              .OR.((I>K).AND.(K>J).AND.(J>L)) &
              .OR.((I>K).AND.(K>L).AND.(L>J))) THEN
        EM(K,I)=EM(K,I)+INTGRL*DM(L,J)
        EM(L,I)=EM(L,I)+INTGRL*DM(K,J)
        EM(K,J)=EM(K,J)+INTGRL*DM(L,I)
        EM(L,J)=EM(L,J)+INTGRL*DM(K,I)
        EM(I,K)=EM(I,K)+INTGRL*DM(J,L)
        EM(I,L)=EM(I,L)+INTGRL*DM(J,K)
        EM(J,K)=EM(J,K)+INTGRL*DM(I,L)
        EM(J,L)=EM(J,L)+INTGRL*DM(I,K)
     END IF
  END DO
  IF (.NOT.DIRECT.AND.USEDISK) CLOSE(BIUNIT)
  N=0
  DO J=1,NBAST
     DO I=1,J
        N=N+1
        PEM(N)=EM(I,J)
     END DO
  END DO
END SUBROUTINE

SUBROUTINE BUILDTAMCM(PTAMCM,PHI,NBAST,NBAS,COMPONENT)
! Computation and assembly of the matrix associated to one of the three components of the total angular momentum operator J (only the upper triangular part of the matrix is stored in packed format).
  USE basis_parameters ; USE integrals
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTAMCM
  TYPE(twospinor),DIMENSION(NBAST),INTENT(IN) :: PHI
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS
  INTEGER :: COMPONENT

  INTEGER :: I,J,K,L,M
  DOUBLE COMPLEX :: VALUE

  PTAMCM=(0.D0,0.D0)
  SELECT CASE (COMPONENT)
  CASE (1)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                         &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),3,2)  &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),2,3))
              END DO
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))   &
 &                        *OVERLAPVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))   &
 &                        *OVERLAPVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L))
           END DO
        END DO
        PTAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                         &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),3,2)  &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),2,3))
              END DO
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))   &
 &                        *OVERLAPVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))   &
 &                        *OVERLAPVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L))   
           END DO
        END DO
        PTAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  CASE (2)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                         &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),1,3)  &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),3,1))
              END DO
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE-.5D0*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))                &
 &                        *DCMPLX(0.D0,OVERLAPVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L)))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))                &
 &                        *DCMPLX(0.D0,OVERLAPVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L)))
           END DO
        END DO
        PTAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                         &
 &                           *DCMPLX(0.D0,XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),1,3)  &
 &                                        -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),3,1))
              END DO
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(1)
           DO M=1,PHI(J)%nbrofcontractions(2)
              VALUE=VALUE-.5D0*PHI(J)%coefficients(2,M)*CONJG(PHI(I)%coefficients(1,L))                &
 &                        *DCMPLX(0.D0,OVERLAPVALUE(PHI(J)%contractions(2,M),PHI(I)%contractions(1,L)))
           END DO
        END DO
        DO L=1,PHI(I)%nbrofcontractions(2)
           DO M=1,PHI(J)%nbrofcontractions(1)
              VALUE=VALUE+.5D0*PHI(J)%coefficients(1,M)*CONJG(PHI(I)%coefficients(2,L))                &
 &                        *DCMPLX(0.D0,OVERLAPVALUE(PHI(J)%contractions(1,M),PHI(I)%contractions(2,L)))
           END DO
        END DO
        PTAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  CASE (3)
  DO J=1,NBAS(1)
     DO I=1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                                 &
 &                           *DCMPLX(.5D0*(-1.D0)**K*OVERLAPVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L)), &
 &                                   XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),2,1)               &
 &                                   -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),1,2))
              END DO
           END DO
        END DO
        PTAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  DO J=NBAS(1)+1,SUM(NBAS)
     DO I=NBAS(1)+1,J
        VALUE=(0.D0,0.D0)
        DO K=1,2
           DO L=1,PHI(I)%nbrofcontractions(K)
              DO M=1,PHI(J)%nbrofcontractions(K)
                 VALUE=VALUE-PHI(J)%coefficients(K,M)*CONJG(PHI(I)%coefficients(K,L))                                 &
 &                           *DCMPLX(.5D0*(-1.D0)**K*OVERLAPVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L)), &
 &                                   XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),2,1)               &
 &                                   -XDERIVVALUE(PHI(J)%contractions(K,M),PHI(I)%contractions(K,L),1,2))
              END DO
           END DO
        END DO
        PTAMCM(I+(J-1)*J/2)=VALUE
     END DO
  END DO
  END SELECT
END SUBROUTINE

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

  SUBROUTINE BUILDTEFM_UHF(PTEFM,NBAST,PHI,PDMA,PDMB)
    USE basis_parameters
    INTEGER,INTENT(IN) :: NBAST
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(OUT) :: PTEFM
    TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
    DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: PDMA,PDMB
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
