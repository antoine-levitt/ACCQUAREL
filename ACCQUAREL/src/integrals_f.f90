MODULE integrals
! Note: all the CGTO integrals are computed by the A.S.P.I.C. code (written in C++ by F. Lodier, see http://www.ann.jussieu.fr/A.S.P.I.C/).
  IMPLICIT NONE
! number of a priori nonzero bielectronic integrals
  INTEGER :: BINMBR
! arrays for the list, values (real/complex for CGTO/2-spinor basis functions in the non-relativistic/relativistic case) and "class" (relativistic case only) of bielectronic integrals (when stored in memory)
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: BILIST
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: RBIVALUES
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: CBIVALUES
  CHARACTER(2),DIMENSION(:),ALLOCATABLE :: BITYPE
! arrays for the values of precomputed CGTO bielectronic integrals used to compute more efficiently (thanks to the use of symmetries) the 2-spinor bielectronic integrals in the relativistic case
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: LLIJKL,LLIKJL,LLILJK,SLIJKL,SSIJKL,SSIKJL,SSILJK
  INTEGER,DIMENSION(2) :: NBF
! unit number for the list of nonzero bielectronic integrals (when stored on disk)
  INTEGER,PARAMETER :: LUNIT=10
! unit number for both the list and values of nonzero bielectronic integrals (when stored on disk)
  INTEGER,PARAMETER :: BIUNIT=11

INTERFACE
  REAL(KIND=C_DOUBLE) FUNCTION unrolloverlap(nbrofprimitives_a,center_a,exponents_a,coefficients_a,monomialdegree_a, &
 &                       nbrofprimitives_b,center_b,exponents_b,coefficients_b,monomialdegree_b)                     &
 & BIND(C,NAME="unrolloverlap")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(KIND=C_INT),VALUE :: nbrofprimitives_a,nbrofprimitives_b
    REAL(KIND=C_DOUBLE),DIMENSION(3) :: center_a,center_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: exponents_a,exponents_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: coefficients_a,coefficients_b
    INTEGER(KIND=C_INT),DIMENSION(3) :: monomialdegree_a,monomialdegree_b
  END FUNCTION
END INTERFACE

INTERFACE
  REAL(KIND=C_DOUBLE) FUNCTION unrollkinetic(nbrofprimitives_a,center_a,exponents_a,coefficients_a,monomialdegree_a, &
 &                                           nbrofprimitives_b,center_b,exponents_b,coefficients_b,monomialdegree_b) &
 & BIND(C,NAME="unrollkinetic")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(KIND=C_INT),VALUE :: nbrofprimitives_a,nbrofprimitives_b
    REAL(KIND=C_DOUBLE),DIMENSION(3) :: center_a,center_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: exponents_a,exponents_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: coefficients_a,coefficients_b
    INTEGER(KIND=C_INT),DIMENSION(3) :: monomialdegree_a,monomialdegree_b
  END FUNCTION
END INTERFACE

INTERFACE
  REAL(KIND=C_DOUBLE) FUNCTION unrollderiv(nbrofprimitives_a,center_a,exponents_a,coefficients_a,monomialdegree_a, &
 &                                         nbrofprimitives_b,center_b,exponents_b,coefficients_b,monomialdegree_b, &
 &                                         dimension)                                                              &
 & BIND(C,NAME="unrollderiv")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(KIND=C_INT),VALUE :: nbrofprimitives_a,nbrofprimitives_b
    REAL(KIND=C_DOUBLE),DIMENSION(3) :: center_a,center_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: exponents_a,exponents_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: coefficients_a,coefficients_b
    INTEGER(KIND=C_INT),DIMENSION(3) :: monomialdegree_a,monomialdegree_b
    INTEGER(KIND=C_INT),VALUE :: dimension
  END FUNCTION
END INTERFACE

INTERFACE
  REAL(KIND=C_DOUBLE) FUNCTION unrollpotential(nbrofprimitives_a,center_a,exponents_a,coefficients_a,monomialdegree_a, &
 &                                             nbrofprimitives_b,center_b,exponents_b,coefficients_b,monomialdegree_b, &
 &                                             center)                                                                 &
 & BIND(C,NAME="unrollpotential")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(KIND=C_INT),VALUE :: nbrofprimitives_a,nbrofprimitives_b
    REAL(KIND=C_DOUBLE),DIMENSION(3) :: center_a,center_b,center
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: exponents_a,exponents_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: coefficients_a,coefficients_b
    INTEGER(KIND=C_INT),DIMENSION(3) :: monomialdegree_a,monomialdegree_b
  END FUNCTION
END INTERFACE

INTERFACE
  REAL(KIND=C_DOUBLE) FUNCTION unrollxderiv(nbrofprimitives_a,center_a,exponents_a,coefficients_a,monomialdegree_a, &
 &                                          nbrofprimitives_b,center_b,exponents_b,coefficients_b,monomialdegree_b, &
 &                                          dimension1,dimension2)                                                  &
 & BIND(C,NAME="unrollxderiv")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(KIND=C_INT),VALUE :: nbrofprimitives_a,nbrofprimitives_b
    REAL(KIND=C_DOUBLE),DIMENSION(3) :: center_a,center_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: exponents_a,exponents_b
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: coefficients_a,coefficients_b
    INTEGER(KIND=C_INT),DIMENSION(3) :: monomialdegree_a,monomialdegree_b
    INTEGER(KIND=C_INT),VALUE :: dimension1,dimension2
  END FUNCTION
END INTERFACE

INTERFACE
  REAL(KIND=C_DOUBLE) FUNCTION unrollcoulomb(nbrofprimitives_a,center_a,exponents_a,coefficients_a,monomialdegree_a, &
 &                                           nbrofprimitives_b,center_b,exponents_b,coefficients_b,monomialdegree_b, &
 &                                           nbrofprimitives_c,center_c,exponents_c,coefficients_c,monomialdegree_c, &
 &                                           nbrofprimitives_d,center_d,exponents_d,coefficients_d,monomialdegree_d) &
 & BIND(C,NAME="unrollcoulomb")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(KIND=C_INT),VALUE :: nbrofprimitives_a,nbrofprimitives_b,nbrofprimitives_c,nbrofprimitives_d
    REAL(KIND=C_DOUBLE),DIMENSION(3) :: center_a,center_b,center_c,center_d
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: exponents_a,exponents_b,exponents_c,exponents_d
    REAL(KIND=C_DOUBLE),DIMENSION(6) :: coefficients_a,coefficients_b,coefficients_c,coefficients_d
    INTEGER(KIND=C_INT),DIMENSION(3) :: monomialdegree_a,monomialdegree_b,monomialdegree_c,monomialdegree_d
  END FUNCTION
END INTERFACE

INTERFACE COULOMBVALUE
  MODULE PROCEDURE COULOMBVALUE_relativistic,COULOMBVALUE_nonrelativistic,COULOMBVALUE_precomputed
END INTERFACE

INTERFACE BUILDBILIST
  MODULE PROCEDURE BUILDBILIST_relativistic,BUILDBILIST_nonrelativistic
END INTERFACE

CONTAINS

FUNCTION POINTVALUE(PHI, X) RESULT(VALUE)
! Function that computes the value of PHI at point X
  USE iso_c_binding ; USE basis_parameters
  TYPE(gaussianbasisfunction),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(3), INTENT(IN) :: X
  REAL(KIND=C_DOUBLE) :: VALUE, MONOMIAL
  INTEGER :: I,D

  VALUE = PRODUCT((X - PHI%center)**PHI%monomialdegree) &
       & * DOT_PRODUCT(PHI%coefficients(1:PHI%nbrofexponents), EXP(-PHI%exponents(1:PHI%nbrofexponents) * SUM((X - PHI%center)**2)))

END FUNCTION POINTVALUE
  
FUNCTION OVERLAPVALUE(PHI_A,PHI_B) RESULT (VALUE)
! Function that computes the value of the integral over R^3 of the product of two GTO basis functions.
  USE iso_c_binding ; USE basis_parameters
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),INTENT(IN) :: PHI_A,PHI_B
  REAL(KIND=C_DOUBLE) :: VALUE

  VALUE=unrolloverlap(PHI_A%nbrofexponents,PHI_A%center,PHI_A%exponents,PHI_A%coefficients,PHI_A%monomialdegree, &
 &                    PHI_B%nbrofexponents,PHI_B%center,PHI_B%exponents,PHI_B%coefficients,PHI_B%monomialdegree)
END FUNCTION

FUNCTION KINETICVALUE(PHI_A,PHI_B) RESULT (VALUE)
! Function that computes the value of the integral over R^3 of the scalar product between the gradients of two GTO basis functions.
  USE iso_c_binding ; USE basis_parameters
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),INTENT(IN) :: PHI_A,PHI_B
  REAL(KIND=C_DOUBLE) :: VALUE

  VALUE=unrollkinetic(PHI_A%nbrofexponents,PHI_A%center,PHI_A%exponents,PHI_A%coefficients,PHI_A%monomialdegree, &
 &                    PHI_B%nbrofexponents,PHI_B%center,PHI_B%exponents,PHI_B%coefficients,PHI_B%monomialdegree)
END FUNCTION

FUNCTION DERIVVALUE(PHI_A,PHI_B,DIMENSION) RESULT (VALUE)
! Function that computes the value of the integral over R^3 of the product of the partial derivative (with respect to a space variable) of a basis function with another basis function (this kind of integrals appear in the variational formulation involving the Dirac operator).
! Note: if DIMENSION = 1 (respectively 2, 3) then partial derivative with respect to x (respectively y, z).
  USE iso_c_binding ; USE basis_parameters
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),INTENT(IN) :: PHI_A,PHI_B
  INTEGER(KIND=C_INT),INTENT(IN) :: DIMENSION
  REAL(KIND=C_DOUBLE) :: VALUE

  VALUE=unrollderiv(PHI_A%nbrofexponents,PHI_A%center,PHI_A%exponents,PHI_A%coefficients,PHI_A%monomialdegree, &
 &                  PHI_B%nbrofexponents,PHI_B%center,PHI_B%exponents,PHI_B%coefficients,PHI_B%monomialdegree, &
 &                  DIMENSION-1)
END FUNCTION

FUNCTION POTENTIALVALUE(PHI_A,PHI_B,CENTER) RESULT (VALUE)
! Function that computes the value of the integral over R^3 of the product between two GTO basis functions and a coulombic potential centered on a given point.
  USE iso_c_binding ; USE basis_parameters
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),INTENT(IN) :: PHI_A,PHI_B
  REAL(KIND=C_DOUBLE),DIMENSION(3),INTENT(IN) :: CENTER
  REAL(KIND=C_DOUBLE) :: VALUE

  VALUE=unrollpotential(PHI_A%nbrofexponents,PHI_A%center,PHI_A%exponents,PHI_A%coefficients,PHI_A%monomialdegree, &
 &                      PHI_B%nbrofexponents,PHI_B%center,PHI_B%exponents,PHI_B%coefficients,PHI_B%monomialdegree, &
 &                      CENTER)
END FUNCTION

FUNCTION XDERIVVALUE(PHI_A,PHI_B,DIMENSION1,DIMENSION2) RESULT (VALUE)
! Function that computes the value of the integral over R^3 of the product of the partial derivative (with respect to the space variable x, y or z) of a GTO basis function with another GTO basis function, times x, y or z (this kind of integral appears in variational formulations involving the J operator).
! Notes: - if DIMENSION1 = 1 (respectively 2, 3) then partial derivative with respect to x (respectively y, z).
!        - if DIMENSION2 = 1 (respectively 2, 3) then the product is multiplied by x (respectively y, z).
  USE iso_c_binding ; USE basis_parameters
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),INTENT(IN) :: PHI_A,PHI_B
  INTEGER(KIND=C_INT),INTENT(IN) :: DIMENSION1,DIMENSION2
  REAL(KIND=C_DOUBLE) :: VALUE

  VALUE=unrollxderiv(PHI_A%nbrofexponents,PHI_A%center,PHI_A%exponents,PHI_A%coefficients,PHI_A%monomialdegree, &
 &                   PHI_B%nbrofexponents,PHI_B%center,PHI_B%exponents,PHI_B%coefficients,PHI_B%monomialdegree, &
 &                   DIMENSION1-1,DIMENSION2-1)
END FUNCTION

FUNCTION COULOMBVALUE_nonrelativistic(PHI_A,PHI_B,PHI_C,PHI_D) RESULT (VALUE)
! Function that computes the value of the bielectronic integral between four GTO basis functions.
  USE iso_c_binding ; USE basis_parameters
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),INTENT(IN) :: PHI_A,PHI_B,PHI_C,PHI_D
  REAL(KIND=C_DOUBLE) :: VALUE

  VALUE=unrollcoulomb(PHI_A%nbrofexponents,PHI_A%center,PHI_A%exponents,PHI_A%coefficients,PHI_A%monomialdegree, &
 &                    PHI_B%nbrofexponents,PHI_B%center,PHI_B%exponents,PHI_B%coefficients,PHI_B%monomialdegree, &
 &                    PHI_C%nbrofexponents,PHI_C%center,PHI_C%exponents,PHI_C%coefficients,PHI_C%monomialdegree, &
 &                    PHI_D%nbrofexponents,PHI_D%center,PHI_D%exponents,PHI_D%coefficients,PHI_D%monomialdegree)
END FUNCTION

SUBROUTINE BUILDBILIST_nonrelativistic(PHI,NBAS,LISTSIZE)
! Subroutine that generates the list (without redundancy as symmetries are taken into account) of the bielectronic integrals with nonzero value.
! Reference: R. Ahlrichs, Methods for efficient evaluation of integrals for gaussian type basis sets, Theoret. Chim. Acta, 33, 157-167, 1974.
  USE case_parameters ; USE basis_parameters
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),DIMENSION(:),INTENT(IN) :: PHI
  INTEGER,DIMENSION(1),INTENT(IN) :: NBAS
  INTEGER,INTENT(OUT) :: LISTSIZE

  INTEGER :: I,J,K,L
  INTEGER,DIMENSION(3) :: GLOBALMONOMIALDEGREE
  LOGICAL :: SC

  OPEN(LUNIT,form='UNFORMATTED')
! determination of the number of elements (i.e., integer quadruples) that compose the list
  LISTSIZE=0
  DO I=1,NBAS(1) ; DO J=1,I ; DO K=1,J ; DO L=1,K
     SC=((PHI(I)%center_id==PHI(J)%center_id).AND.(PHI(J)%center_id==PHI(K)%center_id).AND.(PHI(K)%center_id==PHI(L)%center_id))
     GLOBALMONOMIALDEGREE=PHI(I)%monomialdegree+PHI(J)%monomialdegree+PHI(K)%monomialdegree+PHI(L)%monomialdegree
! parity check on the product of the monomials if the four functions share the same center
     IF ((SC.AND.ALL(MOD(GLOBALMONOMIALDEGREE,2)==0)).OR.(.NOT.SC)) THEN
        LISTSIZE=LISTSIZE+1 ; WRITE(LUNIT)I,J,K,L
        IF (K<J) THEN
           LISTSIZE=LISTSIZE+1 ; WRITE(LUNIT)I,K,J,L
        END IF
        IF ((J<I).AND.(L<K)) THEN
           LISTSIZE=LISTSIZE+1 ; WRITE(LUNIT)I,L,J,K
        END IF
     END IF
  END DO ; END DO ; END DO ; END DO
  CLOSE(LUNIT)
  WRITE(*,*)' Number of CGTO bielectronic integrals to be computed:',LISTSIZE
END SUBROUTINE

FUNCTION COULOMBVALUE_relativistic(PHI_A,PHI_B,PHI_C,PHI_D) RESULT (VALUE)
! Function that computes the value of the bielectronic integral between four 2-spinor basis functions.
  USE basis_parameters
  IMPLICIT NONE
  TYPE(twospinor),INTENT(IN) :: PHI_A,PHI_B,PHI_C,PHI_D
  DOUBLE COMPLEX :: VALUE

  INTEGER :: I,IA,IB,J,JC,JD

  VALUE=(0.D0,0.D0)
  DO I=1,2
     DO IA=1,PHI_A%nbrofcontractions(I) ; DO IB=1,PHI_B%nbrofcontractions(I)
        DO J=1,2
           DO JC=1,PHI_C%nbrofcontractions(J) ; DO JD=1,PHI_D%nbrofcontractions(J)
              VALUE=VALUE+PHI_A%coefficients(I,IA)*CONJG(PHI_B%coefficients(I,IB))                                  &
 &                        *PHI_C%coefficients(J,JC)*CONJG(PHI_D%coefficients(J,JD))                                 &
 &                        *COULOMBVALUE(PHI_A%contractions(I,IA),PHI_B%contractions(I,IB),PHI_C%contractions(J,JC), &
 &                                      PHI_D%contractions(J,JD))
           END DO ; END DO
        END DO
     END DO ; END DO
  END DO
END FUNCTION

FUNCTION COULOMBVALUE_precomputed(PHI_A,PHI_B,PHI_C,PHI_D,CLASS) RESULT (VALUE)
! Function that computes the value of the bielectronic integral between four 2-spinor basis functions from lists containing the precomputed values of the bielectronic integrals between scalar GTO functions.
  USE basis_parameters
  IMPLICIT NONE
  TYPE(twospinor),INTENT(IN) :: PHI_A,PHI_B,PHI_C,PHI_D
  CHARACTER(2),INTENT(IN) :: CLASS
  DOUBLE COMPLEX :: VALUE

  INTEGER :: I,IA,IB,J,JC,JD

  VALUE=(0.D0,0.D0)
  DO I=1,2
     DO IA=1,PHI_A%nbrofcontractions(I) ; DO IB=1,PHI_B%nbrofcontractions(I)
        DO J=1,2
           DO JC=1,PHI_C%nbrofcontractions(J) ; DO JD=1,PHI_D%nbrofcontractions(J)
              VALUE=VALUE+PHI_A%coefficients(I,IA)*CONJG(PHI_B%coefficients(I,IB))                                      &
 &                        *PHI_C%coefficients(J,JC)*CONJG(PHI_D%coefficients(J,JD))                                     &
 &                        *PRECOMPUTEDCOULOMBVALUE(PHI_A%contidx(I,IA),PHI_B%contidx(I,IB),PHI_C%contidx(J,JC),         &
 &                                                 PHI_D%contidx(J,JD),CLASS)
           END DO ; END DO
        END DO
     END DO ; END DO
  END DO
END FUNCTION

SUBROUTINE BUILDBILIST_relativistic(PHI,NBAS,LISTSIZE,SUBSIZE)
! Subroutine that generates the list (more or less without redundancy since the a priori symmetries for complex 2-spinor functions are taken into account) of the bielectronic integrals with nonzero value.
  USE case_parameters ; USE basis_parameters ; USE scf_parameters
  IMPLICIT NONE
  TYPE(twospinor),DIMENSION(:),INTENT(IN) :: PHI
  INTEGER,DIMENSION(2),INTENT(IN) :: NBAS
  INTEGER,INTENT(OUT) :: LISTSIZE,SUBSIZE(3)

  INTEGER :: I,J,K,L,I1,I2,I3,I4,I5,I6
  INTEGER,DIMENSION(3) :: GLOBALMONOMIALDEGREE
  LOGICAL :: SC

  OPEN(LUNIT,form='UNFORMATTED')
! determination of the number of elements (i.e., integer quadruples) that compose the list
  SUBSIZE=0
! LLLL-type integrals
  DO I=1,NBAS(1) ; DO J=1,NBAS(1) ; DO K=1,NBAS(1) ; DO L=1,NBAS(1)
     IF (L+K*NBAS(1)<=J+I*NBAS(1)) THEN
        DO I1=1,2
           DO I2=1,PHI(I)%nbrofcontractions(I1)
              DO I3=1,PHI(J)%nbrofcontractions(I1)
                 DO I4=1,2
                    DO I5=1,PHI(K)%nbrofcontractions(I4)
                       DO I6=1,PHI(L)%nbrofcontractions(I4)
                          SC=((PHI(I)%contractions(I1,I2)%center_id==PHI(J)%contractions(I1,I3)%center_id)      &
 &                            .AND.(PHI(J)%contractions(I1,I3)%center_id==PHI(K)%contractions(I4,I5)%center_id) &
 &                            .AND.(PHI(K)%contractions(I4,I5)%center_id==PHI(L)%contractions(I4,I6)%center_id))
                          GLOBALMONOMIALDEGREE= PHI(I)%contractions(I1,I2)%monomialdegree &
 &                                             +PHI(J)%contractions(I1,I3)%monomialdegree &
 &                                             +PHI(K)%contractions(I4,I5)%monomialdegree &
 &                                             +PHI(L)%contractions(I4,I6)%monomialdegree
                          IF ((SC.AND.ALL(MOD(GLOBALMONOMIALDEGREE,2)==0)).OR.(.NOT.SC)) THEN
                             SUBSIZE(1)=SUBSIZE(1)+1
                             WRITE(LUNIT)I,J,K,L,'LL'
                             GOTO 1
                          END IF
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
1       CONTINUE
     END IF
  END DO ; END DO ; END DO ; END DO
! SSLL-type integrals
  DO I=NBAS(1)+1,SUM(NBAS) ; DO J=NBAS(1)+1,SUM(NBAS) ; DO K=1,NBAS(1) ; DO L=1,NBAS(1)
     DO I1=1,2
        DO I2=1,PHI(I)%nbrofcontractions(I1)
           DO I3=1,PHI(J)%nbrofcontractions(I1)
              DO I4=1,2
                 DO I5=1,PHI(K)%nbrofcontractions(I4)
                    DO I6=1,PHI(L)%nbrofcontractions(I4)
                       SC=((PHI(I)%contractions(I1,I2)%center_id==PHI(J)%contractions(I1,I3)%center_id)      &
 &                         .AND.(PHI(J)%contractions(I1,I3)%center_id==PHI(K)%contractions(I4,I5)%center_id) &
 &                         .AND.(PHI(K)%contractions(I4,I5)%center_id==PHI(L)%contractions(I4,I6)%center_id))
                       GLOBALMONOMIALDEGREE= PHI(I)%contractions(I1,I2)%monomialdegree &
 &                                          +PHI(J)%contractions(I1,I3)%monomialdegree &
 &                                          +PHI(K)%contractions(I4,I5)%monomialdegree &
 &                                          +PHI(L)%contractions(I4,I6)%monomialdegree
                       IF ((SC.AND.ALL(MOD(GLOBALMONOMIALDEGREE,2)==0)).OR.(.NOT.SC)) THEN
                          WRITE(LUNIT)I,J,K,L,'SL'
                          SUBSIZE(2)=SUBSIZE(2)+1
                          GOTO 2
                       END IF
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
2    CONTINUE
  END DO ; END DO ; END DO ; END DO
  IF (SSINTEGRALS) THEN
! SSSS-type integrals
     DO I=NBAS(1)+1,SUM(NBAS) ; DO J=NBAS(1)+1,SUM(NBAS) ; DO K=NBAS(1)+1,SUM(NBAS) ; DO L=NBAS(1)+1,SUM(NBAS)
        IF (L+K*NBAS(2)<=J+I*NBAS(2)) THEN
           DO I1=1,2
              DO I2=1,PHI(I)%nbrofcontractions(I1)
                 DO I3=1,PHI(J)%nbrofcontractions(I1)
                    DO I4=1,2
                       DO I5=1,PHI(K)%nbrofcontractions(I4)
                          DO I6=1,PHI(L)%nbrofcontractions(I4)
                             SC=((PHI(I)%contractions(I1,I2)%center_id==PHI(J)%contractions(I1,I3)%center_id)      &
 &                               .AND.(PHI(J)%contractions(I1,I3)%center_id==PHI(K)%contractions(I4,I5)%center_id) &
 &                               .AND.(PHI(K)%contractions(I4,I5)%center_id==PHI(L)%contractions(I4,I6)%center_id))
                             GLOBALMONOMIALDEGREE= PHI(I)%contractions(I1,I2)%monomialdegree &
 &                                                +PHI(J)%contractions(I1,I3)%monomialdegree &
 &                                                +PHI(K)%contractions(I4,I5)%monomialdegree &
 &                                                +PHI(L)%contractions(I4,I6)%monomialdegree
                             IF ((SC.AND.ALL(MOD(GLOBALMONOMIALDEGREE,2)==0)).OR.(.NOT.SC)) THEN
                                WRITE(LUNIT)I,J,K,L,'SS'
                                SUBSIZE(3)=SUBSIZE(3)+1
                                GOTO 3
                             END IF
                          END DO
                       END DO
                    END DO
                 END DO
              END DO
           END DO
3          CONTINUE
        END IF
     END DO ; END DO ; END DO ; END DO
  END IF
  LISTSIZE=SUM(SUBSIZE)
  CLOSE(LUNIT)
  WRITE(*,*)' Number of 2-spinor-type orbital bielectronic integrals to be computed:',LISTSIZE
END SUBROUTINE

SUBROUTINE PRECOMPUTECGTOCOULOMBVALUES(CGTO,NCGTO)
! Routine that computes the values of the bielectronic integrals over a cartesian gaussian-type orbital basis, taking into account the eightfold permutational symmetry of the integrals (see R. Ahlrichs, Methods for efficient evaluation of integrals for gaussian type basis sets, Theoret. Chim. Acta, 33, 157-167, 1974). These values are next used to compute more efficiently the bielectronic integrals over a cartesian 2-spinor-type orbital basis in the relativistic case (see the GETPRECOMPUTEDCOULOMBVALUE function).
  USE basis_parameters ; USE scf_parameters
  IMPLICIT NONE
  INTEGER,DIMENSION(2),INTENT(IN) :: NCGTO  
  TYPE(gaussianbasisfunction),DIMENSION(SUM(NCGTO)),INTENT(IN) :: CGTO

  INTEGER :: I,J,K,L,M,N,O
  INTEGER,DIMENSION(3) :: GLOBALMONOMIALDEGREE
  LOGICAL :: SC

  NBF=NCGTO
! (LL|LL) integrals
  WRITE(*,*)'- Computing LL integrals'
  ALLOCATE(LLIJKL(1:NCGTO(1)*(NCGTO(1)+1)*(NCGTO(1)**2+5*NCGTO(1)+6)/24), &
 &         LLIKJL(1:NCGTO(1)*(NCGTO(1)+1)*(NCGTO(1)**2+NCGTO(1)-2)/24),   &
 &         LLILJK(1:NCGTO(1)*(NCGTO(1)+1)*(NCGTO(1)**2-3*NCGTO(1)+2)/24))
  M=0 ; N=0 ; O=0
  DO I=1,NCGTO(1) ; DO J=1,I ; DO K=1,J ; DO L=1,K
     SC=((CGTO(I)%center_id==CGTO(J)%center_id).AND.(CGTO(J)%center_id==CGTO(K)%center_id) &
 &                                             .AND.(CGTO(K)%center_id==CGTO(L)%center_id))
     GLOBALMONOMIALDEGREE=CGTO(I)%monomialdegree+CGTO(J)%monomialdegree+CGTO(K)%monomialdegree+CGTO(L)%monomialdegree
! parity check (one center case)
     IF ((SC.AND.ALL(MOD(GLOBALMONOMIALDEGREE,2)==0)).OR.(.NOT.SC)) THEN
        M=M+1 ; LLIJKL(M)=COULOMBVALUE(CGTO(I),CGTO(J),CGTO(K),CGTO(L))
        IF (K<J) THEN
           N=N+1 ; LLIKJL(N)=COULOMBVALUE(CGTO(I),CGTO(K),CGTO(J),CGTO(L))
        END IF
        IF ((J<I).AND.(L<K)) THEN
           O=O+1 ; LLILJK(O)=COULOMBVALUE(CGTO(I),CGTO(L),CGTO(J),CGTO(K))
        END IF
     ELSE
        M=M+1 ; LLIJKL(M)=(0.D0,0.D0)
        IF (K<J) THEN
           N=N+1 ; LLIKJL(N)=(0.D0,0.D0)
        END IF
        IF ((J<I).AND.(L<K)) THEN
           O=O+1 ; LLILJK(O)=(0.D0,0.D0)
        END IF
     END IF
  END DO ; END DO ; END DO ; END DO
! (SS|LL) integrals
  WRITE(*,*)'- Computing SL integrals'
  ALLOCATE(SLIJKL(1:NCGTO(1)*(NCGTO(1)+1)*NCGTO(2)*(NCGTO(2)+1)/4))
  N=0
  DO I=NCGTO(1)+1,SUM(NCGTO) ; DO J=NCGTO(1)+1,I ; DO K=1,NCGTO(1) ; DO L=1,K
     SC=((CGTO(I)%center_id==CGTO(J)%center_id).AND.(CGTO(J)%center_id==CGTO(K)%center_id) &
 &                                             .AND.(CGTO(K)%center_id==CGTO(L)%center_id))
     GLOBALMONOMIALDEGREE=CGTO(I)%monomialdegree+CGTO(J)%monomialdegree+CGTO(K)%monomialdegree+CGTO(L)%monomialdegree
! parity check (one center case)
     IF ((SC.AND.ALL(MOD(GLOBALMONOMIALDEGREE,2)==0)).OR.(.NOT.SC)) THEN
        N=N+1 ; SLIJKL(N)=COULOMBVALUE(CGTO(I),CGTO(J),CGTO(K),CGTO(L))
     ELSE
        N=N+1 ; SLIJKL(N)=(0.D0,0.D0)
     END IF
  END DO ; END DO ; END DO ; END DO
  IF (SSINTEGRALS) THEN
! (SS|SS) integrals
     WRITE(*,*)'- Computing SS integrals'
     ALLOCATE(SSIJKL(1:NCGTO(2)*(NCGTO(2)+1)*(NCGTO(2)**2+5*NCGTO(2)+6)/24), &
 &            SSIKJL(1:NCGTO(2)*(NCGTO(2)+1)*(NCGTO(2)**2+NCGTO(2)-2)/24),   &
 &            SSILJK(1:NCGTO(2)*(NCGTO(2)+1)*(NCGTO(2)**2-3*NCGTO(2)+2)/24))
     M=0 ; N=0 ; O=0
     DO I=NCGTO(1)+1,SUM(NCGTO) ; DO J=NCGTO(1)+1,I ; DO K=NCGTO(1)+1,J ; DO L=NCGTO(1)+1,K
        SC=((CGTO(I)%center_id==CGTO(J)%center_id).AND.(CGTO(J)%center_id==CGTO(K)%center_id) &
 &                                                .AND.(CGTO(K)%center_id==CGTO(L)%center_id))
        GLOBALMONOMIALDEGREE=CGTO(I)%monomialdegree+CGTO(J)%monomialdegree+CGTO(K)%monomialdegree+CGTO(L)%monomialdegree
! parity check (one center case)
        IF ((SC.AND.ALL(MOD(GLOBALMONOMIALDEGREE,2)==0)).OR.(.NOT.SC)) THEN
           M=M+1 ; SSIJKL(M)=COULOMBVALUE(CGTO(I),CGTO(J),CGTO(K),CGTO(L))
           IF (K<J) THEN
              N=N+1 ; SSIKJL(N)=COULOMBVALUE(CGTO(I),CGTO(K),CGTO(J),CGTO(L))
           END IF
           IF ((J<I).AND.(L<K)) THEN
              O=O+1 ; SSILJK(O)=COULOMBVALUE(CGTO(I),CGTO(L),CGTO(J),CGTO(K))
           END IF
        ELSE
           M=M+1 ; SSIJKL(M)=(0.D0,0.D0)
           IF (K<J) THEN
              N=N+1 ; SSIKJL(N)=(0.D0,0.D0)
           END IF
           IF ((J<I).AND.(L<K)) THEN
              O=O+1 ; SSILJK(O)=(0.D0,0.D0)
           END IF
        END IF
     END DO ; END DO ; END DO ; END DO
  END IF
END SUBROUTINE

FUNCTION PRECOMPUTEDCOULOMBVALUE(I,J,K,L,CLASS) RESULT(VALUE)
! Functions that returns the value of a precomputed bielectronic integral of class LL, SL or SS between four (real) cartesian gaussian-type orbital scalar basis functions stored in a list taking into account the eightfold permutational symmetry of the integrals (see R. Ahlrichs, Methods for efficient evaluation of integrals for gaussian type basis sets, Theoret. Chim. Acta, 33, 157-167, 1974).
! note: this function is called for the computation of bielectronic integrals over a (complex) 2-spinor, cartesian gaussian-type orbital basis, which does not naturally possess as many symmetries as a real scalar basis.
  USE basis_parameters
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: I,J,K,L
  CHARACTER(2),INTENT(IN) :: CLASS
  DOUBLE COMPLEX :: VALUE

  INTEGER :: IDX
  INTEGER,DIMENSION(4) :: N,TMP

! Preliminaries
  N=(/I,J,K,L/)
  TMP=N ; IF (N(1)<N(2)) N(1:2)=(/TMP(2),TMP(1)/)
  TMP=N ; IF (N(3)<N(4)) N(3:4)=(/TMP(4),TMP(3)/)
  TMP=N ; IF (N(1)<N(3)) N=(/TMP(3),TMP(4),TMP(1),TMP(2)/)
  IF (CLASS.EQ.'LL') THEN
! integral (LL|LL) between four "upper 2-spinor" scalar CGTO functions
     IF (N(3)<=N(2)) THEN
! integral of type (IJ|KL) 
        IDX=N(4)+N(3)*(N(3)-1)/2+(N(2)+1)*N(2)*(N(2)-1)/6+(N(1)-1)*N(1)*(N(1)+1)*(N(1)+2)/24
        VALUE=LLIJKL(IDX)
     ELSE IF ((N(1)==N(3)).AND.(N(1)==N(4))) THEN
! integral of type (IJ|KL)
        IDX=N(2)+N(1)*(N(1)-1)/2+(N(1)+1)*N(1)*(N(1)-1)/6+(N(1)-1)*N(1)*(N(1)+1)*(N(1)+2)/24
        VALUE=LLIJKL(IDX)
     ELSE IF (N(4)<=N(2)) THEN
! integral of type (IK|JL) 
        IDX=N(4)+N(2)*(N(2)-1)/2+N(3)*(N(3)-1)*(N(3)-2)/6+N(1)*(N(1)-1)*(N(1)*(N(1)-1)-2)/24
        VALUE=LLIKJL(IDX)
     ELSE IF (N(1)==N(3)) THEN
! integral of type (IK|JL)
        IDX=N(2)+N(4)*(N(4)-1)/2+N(1)*(N(1)-1)*(N(1)-2)/6+N(1)*(N(1)-1)*(N(1)*(N(1)-1)-2)/24
        VALUE=LLIKJL(IDX)
     ELSE
! integral of type (IL|JK)
        IDX=N(2)+(N(4)-1)*(N(4)-2)/2+N(3)*(N(3)-1)*(N(3)-2)/6+N(1)*(N(1)-1)*((N(1)-1)*(N(1)-4)+2)/24
        VALUE=LLILJK(IDX)
     END IF
  ELSE IF (CLASS=='SL') THEN
! integral (SS|LL) between two "lower 2-spinor" and two "upper 2-spinor" scalar CGTO functions
     N(1:2)=N(1:2)-(/NBF(1),NBF(1)/)
     IDX=N(4)+N(3)*(N(3)-1)/2+(N(2)+N(1)*(N(1)-1)/2-1)*NBF(1)*(NBF(1)+1)/2
     VALUE=SLIJKL(IDX)
  ELSE IF (CLASS=='SS') THEN
     N=N-(/NBF(1),NBF(1),NBF(1),NBF(1)/)
! integral (SS|SS) between four "lower 2-spinor" scalar CGTO functions
     IF (N(3)<=N(2)) THEN
! integral of type (IJ|KL) 
        IDX=N(4)+N(3)*(N(3)-1)/2+(N(2)+1)*N(2)*(N(2)-1)/6+(N(1)-1)*N(1)*(N(1)+1)*(N(1)+2)/24
        VALUE=SSIJKL(IDX)
     ELSE IF ((N(1)==N(3)).AND.(N(1)==N(4))) THEN
! integral of type (IJ|KL)
        IDX=N(2)+N(1)*(N(1)-1)/2+(N(1)+1)*N(1)*(N(1)-1)/6+(N(1)-1)*N(1)*(N(1)+1)*(N(1)+2)/24
        VALUE=SSIJKL(IDX)
     ELSE IF (N(4)<=N(2)) THEN
! integral of type (IK|JL) 
        IDX=N(4)+N(2)*(N(2)-1)/2+N(3)*(N(3)-1)*(N(3)-2)/6+N(1)*(N(1)-1)*(N(1)*(N(1)-1)-2)/24
        VALUE=SSIKJL(IDX)
     ELSE IF (N(1)==N(3)) THEN
! integral of type (IK|JL)
        IDX=N(2)+N(4)*(N(4)-1)/2+N(1)*(N(1)-1)*(N(1)-2)/6+N(1)*(N(1)-1)*(N(1)*(N(1)-1)-2)/24
        VALUE=SSIKJL(IDX)
     ELSE
! integral of type (IL|JK)
        IDX=N(2)+(N(4)-1)*(N(4)-2)/2+N(3)*(N(3)-1)*(N(3)-2)/6+N(1)*(N(1)-1)*((N(1)-1)*(N(1)-4)+2)/24
        VALUE=SSILJK(IDX)
     END IF
  END IF  
END FUNCTION

SUBROUTINE DEALLOCATE_INTEGRALS
  USE scf_parameters
! Routine that deallocate the arrays containing the values of the CGTO bielectronic integrals.
  IMPLICIT NONE
  DEALLOCATE(LLIJKL,LLIKJL,LLILJK,SLIJKL)
  IF (SSINTEGRALS) DEALLOCATE(SSIJKL,SSIKJL,SSILJK)
END SUBROUTINE


FUNCTION COMPUTE_DENSITY(X, PDM, PHI, N) RESULT(VALUE)
! Subroutine that computes the density defined by PDM in the PHI basis at point X
  USE basis_parameters ; USE matrix_tools
  IMPLICIT NONE
  integer N, I
  DOUBLE PRECISION,DIMENSION(N*(N+1)/2) :: PDM
  TYPE(gaussianbasisfunction),DIMENSION(N),INTENT(IN) :: PHI
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: X
  real, dimension(N) :: BASISVALUE
  real value

  DO I=1, N
     BASISVALUE(I) = POINTVALUE(PHI(I), X)
  END DO

  VALUE = DOT_PRODUCT(BASISVALUE, MATMUL(UNPACK(PDM, N), BASISVALUE))
  
END FUNCTION

END MODULE
