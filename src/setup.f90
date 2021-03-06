MODULE case_parameters
! relativistic or non-relativistic case flag
  LOGICAL :: RELATIVISTIC
! apparent speed of light (can be defined by the user via a scaling factor) 
  DOUBLE PRECISION :: C
! approximation index (1 is Hartree, 2 is Hartree-Fock) 
  INTEGER :: APPROXIMATION
! model/formalism index (see below)
  INTEGER :: MODEL
CONTAINS

SUBROUTINE SETUP_CASE
  USE constants ; USE setup_tools
  IMPLICIT NONE
  DOUBLE PRECISION :: SCALING_FACTOR
  CHARACTER(2) :: CHAR
  INTEGER :: INFO

  OPEN(100,FILE=SETUP_FILE,STATUS='old',ACTION='read')
  CALL LOOKFOR(100,'## CASE',INFO)
  IF (INFO/=0) STOP' The computational case is not given.'
  READ(100,'(a2)') CHAR
  IF (CHAR=='RE') THEN
     RELATIVISTIC=.TRUE.
     WRITE(*,'(a)')' Relativistic case'
     REWIND(100)
     CALL LOOKFOR(100,'## SPEED OF LIGHT SCALING FACTOR',INFO)
     IF (INFO==0) THEN
        READ(100,*) SCALING_FACTOR
     ELSE
        SCALING_FACTOR=1.D0
     END IF
     C=SCALING_FACTOR*SPEED_OF_LIGHT
     WRITE(*,'(a,f20.12)')' Speed of light c = ',C
  ELSE IF (CHAR=='NO') THEN
     RELATIVISTIC=.FALSE.
     WRITE(*,'(a)')' Non-relativistic case'
  ELSE
     WRITE(*,*)'Subroutine SETUP_CASE: error!'
     STOP
  END IF
  CLOSE(100)
END SUBROUTINE SETUP_CASE

SUBROUTINE SETUP_APPROXIMATION
  USE setup_tools
  IMPLICIT NONE
  CHARACTER(12) :: CHAR
  INTEGER :: INFO

  OPEN(100,FILE=SETUP_FILE,STATUS='old',ACTION='read')
  CALL LOOKFOR(100,'## APPROXIMATION',INFO)
  IF (INFO/=0) STOP' The approximation type is not given.'
  READ(100,'(a12)') CHAR
  IF (LEN_TRIM(CHAR)==7) THEN
     APPROXIMATION=1
     WRITE(*,'(a)')' Hartree approximation (Hartree product wave function)'
     IF (RELATIVISTIC) STOP' This option is incompatible with the previous one (for the time being)' 
  ELSE IF (LEN_TRIM(CHAR)>7) THEN
     APPROXIMATION=2
     WRITE(*,'(a)')' Hartree-Fock approximation (Slater determinant wave function)'
  ELSE
     WRITE(*,*)'Subroutine SETUP_APPROXIMATION: error!'
     STOP
  END IF
  CLOSE(100)
END SUBROUTINE SETUP_APPROXIMATION

SUBROUTINE SETUP_FORMALISM
  USE setup_tools
  IMPLICIT NONE
  CHARACTER(3) :: CHAR
  INTEGER :: INFO

  OPEN(100,FILE=SETUP_FILE,STATUS='old',ACTION='read')
  CALL LOOKFOR(100,'## MODEL',INFO)
  IF (INFO/=0) STOP' The formalism is not given.'
  READ(100,'(a3)') CHAR
  IF (RELATIVISTIC) THEN
     IF (APPROXIMATION==2) THEN
        IF (CHAR=='CSD') THEN
           MODEL=1
           WRITE(*,'(a)')' Closed-shell formalism'
        ELSE IF (CHAR=='OSD') THEN
           MODEL=2
           WRITE(*,'(a)')' Average-of-configuration open-shell formalism'
        ELSE IF (CHAR=='CGH') THEN
           MODEL=3
           WRITE(*,'(a)')' Complex General Hartree-Fock'
        ELSE
           GO TO 1
        END IF
     END IF
  ELSE
     IF (APPROXIMATION==1) THEN
        IF (CHAR=='BOS') THEN
!           MODEL=1
           WRITE(*,'(a)')' Boson star model (preliminary)'
        ELSE
           GO TO 1
        END IF
     ELSE IF (APPROXIMATION==2) THEN
        IF (CHAR=='RHF') THEN
! Restricted Hartree-Fock (RHF) formalism (doubly occupied orbitals)
           MODEL=1
           WRITE(*,'(a)')' (Real) Restricted Hartree-Fock (RHF) formalism'
        ELSE IF (CHAR=='UHF') THEN
! Unrestricted Hartree-Fock (UHF) formalism (different orbitals for different spins (DODS method))
! Reference: J. A. Pople and R. K. Nesbet, Self-consistent orbitals for radicals, J. Chem. Phys., 22(3), 571-572, 1954.
           MODEL=2
           WRITE(*,'(a)')' (Real) Unrestricted (open-shell) Hartree-Fock (UHF) formalism'
        ELSE IF (CHAR=='ROH') THEN
! Restricted Open-shell Hartree-Fock (ROHF) formalism
! Reference: C. C. J. Roothaan, Self-consistent field theory for open shells of electronic systems, Rev. Modern Phys., 32(2), 179-185, 1960.
           MODEL=3
           WRITE(*,'(a)')' (Real) Restricted Open-shell Hartree-Fock (ROHF) formalism'
           WRITE(*,*)'Option not implemented yet!'
           STOP
        ELSE IF (CHAR=='RGH') THEN
! Real General Hartree-Fock (RGHF) formalism
           MODEL=4
           WRITE(*,'(a)')' Real General Hartree-Fock formalism'
!        ELSE IF (CHAR=='CGH') THEN
! Complex General Hartree-Fock (CGHF) formalism
!           MODEL=5
!           WRITE(*,'(a)')' Complex General Hartree-Fock formalism'
        ELSE
           GO TO 1
        END IF
     END IF
  END IF
  CLOSE(100)
  RETURN
1 WRITE(*,*)'Subroutine SETUP_FORMALISM: no known formalism given!'
  STOP
END SUBROUTINE SETUP_FORMALISM
END MODULE

MODULE data_parameters
  USE iso_c_binding
! ** DATA FOR ATOMIC OR MOLECULAR SYSTEMS **
! number of nuclei in the molecular system (at most 10)
  INTEGER :: NBN
! atomic numbers of the nuclei
  INTEGER,DIMENSION(10) :: Z
! positions of the nucleii (in atomic units)
  DOUBLE PRECISION,DIMENSION(3,10) :: CENTER
! flags for the point group symmetries of the molecular system
! - plane symmetries
  LOGICAL :: SYM_SX,SYM_SY,SYM_SZ
! total number of electrons in the molecular system
  INTEGER :: NBE
! total number of electrons in the closed shells (open-shell DHF formalism)
  INTEGER :: NBECS
! number of open shell electrons and number of open shell orbitals (open-shell DHF formalism)
  INTEGER :: NBEOS,NBOOS
! respective numbers of electrons of $\alpha$ and $\beta$ spin (UHF and ROHF formalisms)
  INTEGER :: NBEA,NBEB
! internuclear repulsion energy (in atomic units (Hartrees))
  DOUBLE PRECISION :: INTERNUCLEAR_ENERGY
! ** DATA FOR BOSON STAR MODEL **
  DOUBLE PRECISION,PARAMETER :: KAPPA=-1.D0
! mass
  DOUBLE PRECISION :: MASS
! temperature
  DOUBLE PRECISION :: TEMPERATURE
! exponent for the entropy generating function
  DOUBLE PRECISION :: MB
!
  INTEGER,POINTER :: RANK_P
  DOUBLE PRECISION,POINTER,DIMENSION(:) :: MU_I

CONTAINS

SUBROUTINE SETUP_DATA
  USE case_parameters ; USE setup_tools
  IMPLICIT NONE
  CHARACTER(30) :: NAME
  INTEGER :: INFO,N

  OPEN(100,FILE=SETUP_FILE,STATUS='old',ACTION='read')
  IF (.NOT.RELATIVISTIC.AND.APPROXIMATION==1) THEN
     CALL LOOKFOR(100,'## DESCRIPTION OF THE BOSON STAR',INFO)
     IF (INFO/=0) STOP' The description of the boson star is not given.'
     WRITE(*,'(a)')' --- **** ---'
     READ(100,*) NAME
     READ(100,*) MASS
     WRITE(*,'(a,f5.3)')' * Mass = ',MASS
     READ(100,*) TEMPERATURE
     WRITE(*,'(a,f5.3)')' * Temperature = ',TEMPERATURE
     READ(100,*) MB
     WRITE(*,'(a,f5.3)')' * Exponent for the entropy generating function = ',MB
     WRITE(*,'(a)')' --- **** ---'
  ELSE
     CALL LOOKFOR(100,'## DESCRIPTION OF THE MOLECULAR SYSTEM',INFO)
     IF (INFO/=0) STOP' The description of the molecular system is not given.'
     WRITE(*,'(a)')' --- Molecular system ---'
     READ(100,'(1/,a)')NAME
     WRITE(*,'(a,a)')' ** NAME: ',NAME
     READ(100,'(i2)') NBN
     DO N=1,NBN
        WRITE(*,'(a,i2)')' * Nucleus #',N
        READ(100,*) Z(N),CENTER(:,N)
        WRITE(*,'(a,i3,a,a,a)')'   Charge Z=',Z(N),' (element: ',IDENTIFYZ(Z(N)),')'
        WRITE(*,'(a,3(f16.8))')'   Center position = ',CENTER(:,N)
     END DO
     CALL INTERNUCLEAR_REPULSION_ENERGY
     WRITE(*,'(a,f16.8,a)')' * Internuclear repulsion energy of the system = ',INTERNUCLEAR_ENERGY,' Ha'
     READ(100,'(i3)') NBE
     WRITE(*,'(a,i3)')' * Total number of electrons in the system = ',NBE
     IF (RELATIVISTIC) THEN
        IF (MODEL==2) THEN
           READ(100,'(3(i3))')NBECS,NBEOS,NBOOS
           WRITE(*,'(a,i3)')'   - number of closed-shell electrons = ',NBECS
           WRITE(*,'(a,i3)')'   - number of open-shell electrons = ',NBEOS
           WRITE(*,'(a,i3)')'   - number of open-shell orbitals = ',NBOOS
           IF (NBE/=NBECS+NBEOS) STOP' Problem with the total number of electrons'
           IF (NBOOS<=NBEOS) STOP' Problem with the number of open-shell orbitals!'
        END IF
     ELSE
        IF (MODEL==1) THEN
           IF (MODULO(NBE,2)/=0) STOP' Problem: the number of electrons must be even!'
        ELSE IF (MODEL==2) THEN
           READ(100,'(2(i3))')NBEA,NBEB
           WRITE(*,'(a,i3)')'   - number of electrons of $\alpha$ spin = ',NBEA
           WRITE(*,'(a,i3)')'   - number of electrons of $\beta$ spin = ',NBEB
           IF (NBE/=NBEA+NBEB) STOP' Problem with the total number of electrons!'
        END IF
     END IF
     REWIND(100)
     CALL LOOKFOR(100,'SYMMETRY SX',INFO)
     IF (INFO==0) THEN
        WRITE(*,'(a)')' The system possesses a symmetry wrt the x-plane.'
        SYM_SX=.TRUE.
     ELSE
        SYM_SX=.FALSE.
     END IF
     REWIND(100)
     CALL LOOKFOR(100,'SYMMETRY SY',INFO)
     IF (INFO==0) THEN
        WRITE(*,'(a)')' The system possesses a symmetry wrt the y-plane.'
        SYM_SY=.TRUE.
     ELSE
        SYM_SY=.FALSE.
     END IF
     REWIND(100)
     CALL LOOKFOR(100,'SYMMETRY SZ',INFO)
     IF (INFO==0) THEN
        WRITE(*,'(a)')' The system possesses a symmetry wrt the z-plane.'
        SYM_SZ=.TRUE.
     ELSE
        SYM_SZ=.FALSE.
     END IF
     WRITE(*,'(a)')' --------- **** ---------'
  END IF
  CLOSE(100)
END SUBROUTINE SETUP_DATA

FUNCTION IDENTIFYZ(Z) RESULT (SYMBOL)
! Function returning the symbol of a chemical element given its atomic number Z.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: Z
  CHARACTER(2) :: SYMBOL

  IF (Z>104) THEN
     WRITE(*,*)'Function IDENTIFYZ: unknown chemical element!'
     STOP
  END IF
! List of symbols from Hydrogen up to Rutherfordium.
  SELECT CASE (Z)
     CASE (1)  ; SYMBOL='H'
     CASE (2)  ; SYMBOL='He'
     CASE (3)  ; SYMBOL='Li'
     CASE (4)  ; SYMBOL='Be'
     CASE (5)  ; SYMBOL='B'
     CASE (6)  ; SYMBOL='C'
     CASE (7)  ; SYMBOL='N'
     CASE (8)  ; SYMBOL='O'
     CASE (9)  ; SYMBOL='F'
     CASE (10) ; SYMBOL='Ne'
     CASE (11) ; SYMBOL='Na'
     CASE (12) ; SYMBOL='Mg'
     CASE (13) ; SYMBOL='Al'
     CASE (14) ; SYMBOL='Si'
     CASE (15) ; SYMBOL='P'
     CASE (16) ; SYMBOL='S'
     CASE (17) ; SYMBOL='Cl'
     CASE (18) ; SYMBOL='Ar'
     CASE (19) ; SYMBOL='K'
     CASE (20) ; SYMBOL='Ca'
     CASE (21) ; SYMBOL='Sc'
     CASE (22) ; SYMBOL='Ti'
     CASE (23) ; SYMBOL='V'
     CASE (24) ; SYMBOL='Cr'
     CASE (25) ; SYMBOL='Mn'
     CASE (26) ; SYMBOL='Fe'
     CASE (27) ; SYMBOL='Co'
     CASE (28) ; SYMBOL='Ni'
     CASE (29) ; SYMBOL='Cu'
     CASE (30) ; SYMBOL='Zn'
     CASE (31) ; SYMBOL='Ga'
     CASE (32) ; SYMBOL='Ge'
     CASE (33) ; SYMBOL='As'
     CASE (34) ; SYMBOL='Se'
     CASE (35) ; SYMBOL='Br'
     CASE (36) ; SYMBOL='Kr'
     CASE (37) ; SYMBOL='Rb'
     CASE (38) ; SYMBOL='Sr'
     CASE (39) ; SYMBOL='Y'
     CASE (40) ; SYMBOL='Zr'
     CASE (41) ; SYMBOL='Nb'
     CASE (42) ; SYMBOL='Mo'
     CASE (43) ; SYMBOL='Tc'
     CASE (44) ; SYMBOL='Ru'
     CASE (45) ; SYMBOL='Rh'
     CASE (46) ; SYMBOL='Pd'
     CASE (47) ; SYMBOL='Ag'
     CASE (48) ; SYMBOL='Cd'
     CASE (49) ; SYMBOL='In'
     CASE (50) ; SYMBOL='Sn'
     CASE (51) ; SYMBOL='Sb'
     CASE (52) ; SYMBOL='Te'
     CASE (53) ; SYMBOL='I'
     CASE (54) ; SYMBOL='Xe'
     CASE (55) ; SYMBOL='Cs'
     CASE (56) ; SYMBOL='Ba'
     CASE (57) ; SYMBOL='La'
! Lanthanide elements (lanthanoids)
     CASE (58) ; SYMBOL='Ce'
     CASE (59) ; SYMBOL='Pr'
     CASE (60) ; SYMBOL='Nd'
     CASE (61) ; SYMBOL='Pm'
     CASE (62) ; SYMBOL='Sm'
     CASE (63) ; SYMBOL='Eu'
     CASE (64) ; SYMBOL='Gd'
     CASE (65) ; SYMBOL='Tb'
     CASE (66) ; SYMBOL='Dy'
     CASE (67) ; SYMBOL='Ho'
     CASE (68) ; SYMBOL='Er'
     CASE (69) ; SYMBOL='Tm'
     CASE (70) ; SYMBOL='Yb'
     CASE (71) ; SYMBOL='Lu'

     CASE (72) ; SYMBOL='Hf'
     CASE (73) ; SYMBOL='Ta'
     CASE (74) ; SYMBOL='W'
     CASE (75) ; SYMBOL='Re'
     CASE (76) ; SYMBOL='Os'
     CASE (77) ; SYMBOL='Ir'
     CASE (78) ; SYMBOL='Pt'
     CASE (79) ; SYMBOL='Au'
     CASE (80) ; SYMBOL='Hg'
     CASE (81) ; SYMBOL='Tl'
     CASE (82) ; SYMBOL='Pb'
     CASE (83) ; SYMBOL='Bi'
     CASE (84) ; SYMBOL='Po'
     CASE (85) ; SYMBOL='As'
     CASE (86) ; SYMBOL='Rn'
     CASE (87) ; SYMBOL='Fr'
     CASE (88) ; SYMBOL='Ra'
     CASE (89) ; SYMBOL='Ac'
! Actinide elements (actinoids)
     CASE (90) ; SYMBOL='Th'
     CASE (91) ; SYMBOL='Pa'
     CASE (92) ; SYMBOL='U'
     CASE (93) ; SYMBOL='Np'
     CASE (94) ; SYMBOL='Pu'
     CASE (95) ; SYMBOL='Am'
     CASE (96) ; SYMBOL='Cm'
     CASE (97) ; SYMBOL='Bk'
     CASE (98) ; SYMBOL='Cf'
     CASE (99) ; SYMBOL='Es'
     CASE (100); SYMBOL='Fm'
     CASE (101); SYMBOL='Md'
     CASE (102); SYMBOL='No'
     CASE (103); SYMBOL='Lr'

     CASE (104); SYMBOL='Rf'
  END SELECT
END FUNCTION IDENTIFYZ

SUBROUTINE INTERNUCLEAR_REPULSION_ENERGY
! Function that computes the internuclear repulsion energy for the given specific geometry of the molecular system.
  IMPLICIT NONE
  INTEGER :: I,J
  DOUBLE PRECISION,DIMENSION(3) :: DIFF

  INTERNUCLEAR_ENERGY=0.D0
  DO I=1,NBN
     DO J=I+1,NBN
        DIFF=CENTER(:,I)-CENTER(:,J)
        INTERNUCLEAR_ENERGY=INTERNUCLEAR_ENERGY+Z(I)*Z(J)/SQRT(DOT_PRODUCT(DIFF,DIFF))
     END DO
  END DO
END SUBROUTINE INTERNUCLEAR_REPULSION_ENERGY

! Various functions for the Hartree model with temperature (Hwtemp model)

FUNCTION POSITIVE_PART(X) RESULT(FUNC)
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: X
  DOUBLE PRECISION :: FUNC

  IF (X<0.D0) THEN
     FUNC=0.D0
  ELSE
     FUNC=X
  END IF
END FUNCTION POSITIVE_PART

FUNCTION ENTROPY_FUNCTION(X) RESULT(FUNC)
! entropy generating function beta
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: X
  DOUBLE PRECISION :: FUNC

  FUNC=POSITIVE_PART(X)**MB/MB
END FUNCTION ENTROPY_FUNCTION

FUNCTION DENTFUNC(X) RESULT(FUNC)
! derivative of the entropy generating function beta
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: X
  DOUBLE PRECISION :: FUNC

  FUNC=POSITIVE_PART(X)**(MB-1.D0)
END FUNCTION DENTFUNC

FUNCTION RECIP_DENTFUNC(X) RESULT(FUNC)
! reciprocal function of the derivative of the entropy generating function beta
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: X
  DOUBLE PRECISION :: FUNC

  IF (X<0.D0) THEN
     STOP'beta is not a bijection on R_-'
  ELSE
     FUNC=X**(1.D0/(MB-1.D0))
  END IF
END FUNCTION RECIP_DENTFUNC

FUNCTION DRECIP_DENTFUNC(X) RESULT(FUNC)
! derivative of the reciprocal function of the derivative of the entropy generating function beta
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: X
  DOUBLE PRECISION :: FUNC

  IF (X<0.D0) THEN
     STOP'beta is not a bijection on R_-'
  ELSE IF (X==0.D0) THEN
     STOP'beta is not differentiable at the origin'
  ELSE
     FUNC=X**((2.D0-MB)/(MB-1.D0))/(MB-1.D0)
  END IF
END FUNCTION

FUNCTION FUNCFORMU(X) RESULT(FUNC)
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: X
  DOUBLE PRECISION :: FUNC

  INTEGER :: I
  DOUBLE PRECISION :: Y

  FUNC=-MASS
  DO I=1,RANK_P
     Y=(X-MU_I(I))/TEMPERATURE
     IF (Y>0.D0) FUNC=FUNC+RECIP_DENTFUNC(Y)
  END DO
END FUNCTION FUNCFORMU

SUBROUTINE RDENTFUNCD(X,FVAL,FDERIV)
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: X
  DOUBLE PRECISION,INTENT(OUT) :: FVAL,FDERIV

  INTEGER :: I
  DOUBLE PRECISION :: Y

  FVAL=-MASS ; FDERIV=0.D0
  DO I=1,RANK_P
     Y=(X-MU_I(I))/TEMPERATURE
     IF (Y>0.D0) THEN
        FVAL=FVAL+RECIP_DENTFUNC(Y)
        FDERIV=FDERIV+DRECIP_DENTFUNC(Y)
     END IF
  END DO
END SUBROUTINE RDENTFUNCD
END MODULE

MODULE basis_parameters
  USE iso_c_binding
! flag for the choice of the basis set type (either existing in the library or even-tempered)
  LOGICAL :: LIBRARY
! PARAMETERS FOR A GIVEN BASIS SET
  CHARACTER(26) :: BASISFILE
  INTEGER,PARAMETER :: MAQN=4,MNOP=38,MNOC=38,MNOGBF=4
! Note: MAQN is the maximum number of different cartesian GBF function types (= maximum angular quantum number + 1) allowed, MNOP is the maximum number of primitives (of different exponents) allowed in any of these types, MNOC is the maximum number of contractions allowed in any of these types, MNOGBF is the maximum number of different GBF allowed in each component of a 2-spinor basis function (necessary for the lower 2-spinor basis due to the use of the Restricted Kinetic Balance scheme). MAQN, MNOP and MNOC depend on the basis that is used, MNOGBF depends on MAQN through the RKB scheme.

! PARAMETERS FOR AN EVEN-TEMPERED BASIS SET
  INTEGER :: NUMBER_OF_TERMS
  DOUBLE PRECISION :: FIRST_TERM,COMMON_RATIO

! Various flags for the contraction of the primitives and the Kinetic Balance scheme
  LOGICAL :: UNCONT
  LOGICAL :: KINBAL,UKB

TYPE gaussianbasisfunction
! Definition of a contracted cartesian gaussian type "orbital" (CGTO) basis function.
! nbrofexponents: the number of different gaussian primitives present in the contraction
! center: coordinates (x,y,z) of the center of the basis function
! center_id: number of the nucleus relative to the center of the basis function in the list of the nuclei forming the molecular system (used for checking the parity of the bielectronic integrands when the four basis functions share the same center)
! exponents: array containing the exponent of each of the gaussian primitives present in the contraction
! coefficients: array containing the coefficient of each of the gaussian primitives present in the contraction
! monomialdegrees: array containing the degrees (n_x,n_y,n_z) of the monomial common to each of the gaussian primitives
! Note: the maximum number of terms in a contraction is set to 6 (see the basis for Cr in 6-31G for instance).
  INTEGER(KIND=C_INT) :: nbrofexponents
  REAL(KIND=C_DOUBLE),DIMENSION(3) :: center
  INTEGER :: center_id
  REAL(KIND=C_DOUBLE),DIMENSION(6) :: exponents
  REAL(KIND=C_DOUBLE),DIMENSION(6) :: coefficients
  INTEGER(KIND=C_INT),DIMENSION(3) :: monomialdegree
END TYPE gaussianbasisfunction

TYPE twospinor
! Definition of a Pauli 2-spinor type basis function using gaussian basis functions.
! nbrofcontractions: array containing the number of different contractions (<=MNOGT0) present in each of the components of the 2-spinor
! contractions: array containing the contractions present in each of the components of the 2-spinor
! contidx : array containing the indices of the gaussian primitives appearing in the contractions with respect to a secondary array of gaussian primitives (used for precomputation purposes)
! coefficients: array containing the complex coefficient of each of the contractions present in each of the components of the 2-spinor
! Note: if one of the components of the 2-spinor is zero then the corresponding nbrofcontractions is set to 0
  INTEGER,DIMENSION(2) :: nbrofcontractions
  TYPE(gaussianbasisfunction),DIMENSION(2,MNOGBF) :: contractions
  INTEGER,DIMENSION(2,MNOGBF) :: contidx
  DOUBLE COMPLEX,DIMENSION(2,MNOGBF) :: coefficients
END TYPE twospinor

CONTAINS

SUBROUTINE SETUP_BASIS
  USE case_parameters ; USE setup_tools
  IMPLICIT NONE
  CHARACTER(LEN=4) :: CHAR
  CHARACTER(LEN=26) :: BASISNAME
  INTEGER :: INFO

  OPEN(100,FILE=SETUP_FILE,STATUS='old',ACTION='read')
  CALL LOOKFOR(100,'## BASIS DEFINITION',INFO)
  IF (INFO/=0) STOP' The basis definition is not given.'
  READ(100,'(3/,a)') BASISNAME
  IF (BASISNAME(1:6)=='BASIS ') THEN
! The basis set is an existing one in the basis library
     LIBRARY=.TRUE.
     BASISFILE='basis/'//BASISNAME(7:)
     READ(100,'(a4)') CHAR
     IF (CHAR=='UNCO') THEN
        UNCONT=.TRUE.
        WRITE(*,'(a,a,a)')' Basis set: ',BASISNAME,' (uncontracted)'
     ELSE
        UNCONT=.FALSE.
        WRITE(*,'(a,a,a)')' Basis set: ',BASISNAME,' (contracted)'
     END IF
  ELSE IF (BASISNAME(1:4)=='EVEN') THEN
! The basis set is an even-tempered one
     LIBRARY=.FALSE.
     IF (RELATIVISTIC.OR.MODEL>1) STOP' Option not implemented (even-tempered basis set)'
     WRITE(*,'(a)')' Even-tempered basis set (preliminary support)'
     READ(100,'(i4)') NUMBER_OF_TERMS
     WRITE(*,'(a,i4)')'  * number of exponents = ',NUMBER_OF_TERMS
     READ(100,*) FIRST_TERM
     WRITE(*,'(a,f16.8)')'  * first term of the geometric series = ',FIRST_TERM
     READ(100,*) COMMON_RATIO
     WRITE(*,'(a,f16.8)')'  * common ratio of the geometric series = ',COMMON_RATIO
  ELSE
     STOP' Unknown basis set type'
  END IF
  IF (RELATIVISTIC) THEN
     READ(100,'(a2)') CHAR
     IF (CHAR=='KI') THEN
        KINBAL=.TRUE.
        READ(100,'(a4)') CHAR
        IF (CHAR=='REST') THEN
           UKB=.FALSE.
           WRITE(*,'(a)')' Restricted kinetic balance'
        ELSE
           UKB=.TRUE.
           WRITE(*,'(a)')' (impaired) Unrestricted kinetic balance'
        END IF
     ELSE
        KINBAL=.FALSE.
        WRITE(*,'(a)')' No kinetic balance'
     END IF
  END IF
  CLOSE(100)
END SUBROUTINE SETUP_BASIS

FUNCTION GBF_POINTWISE_VALUE(GBF,POINT) RESULT(VALUE)
! Function that computes the value of a gaussian basis function at a given point of space.
  USE iso_c_binding
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),INTENT(IN) :: GBF
  DOUBLE PRECISION,DIMENSION(3),INTENT(IN) :: POINT
  REAL(KIND=C_DOUBLE) :: VALUE

  VALUE=PRODUCT((POINT-GBF%center)**GBF%monomialdegree)         &
      & *DOT_PRODUCT(GBF%coefficients(1:GBF%nbrofexponents),    &
      &              EXP(-GBF%exponents(1:GBF%nbrofexponents)*SUM((POINT-GBF%center)**2)))
END FUNCTION GBF_POINTWISE_VALUE

SUBROUTINE PRINTGBF(PHI,NUNIT)
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),INTENT(IN) :: PHI
  INTEGER,INTENT(IN) :: NUNIT

  INTEGER :: I

  WRITE(NUNIT,*)' number of exponents:',PHI%nbrofexponents
  WRITE(NUNIT,*)' center:',PHI%center
  WRITE(NUNIT,*)' common monomial:',PHI%monomialdegree
  DO I=1,PHI%nbrofexponents
     WRITE(NUNIT,*)'  gaussian primitive #',I
     WRITE(NUNIT,*)'    exponent:',PHI%exponents(I)
     WRITE(NUNIT,*)'    coefficient:',PHI%coefficients(I)
  END DO
END SUBROUTINE PRINTGBF

FUNCTION TWOSPINOR_POINTWISE_VALUE(PHI,POINT) RESULT(VALUE)
! Function that computes the value of a Pauli 2-spinor basis function at a given point of space.
  USE iso_c_binding
  IMPLICIT NONE
  TYPE(twospinor),INTENT(IN) :: PHI
  DOUBLE PRECISION,DIMENSION(3),INTENT(IN) :: POINT
  DOUBLE COMPLEX,DIMENSION(2) :: VALUE

  INTEGER :: I,J
  
  DO I=1,2
     VALUE(I)=(0.D0,0.D0)
     DO J=1,PHI%nbrofcontractions(I)
        VALUE(I)=VALUE(I)+PHI%coefficients(I,J)*GBF_POINTWISE_VALUE(PHI%contractions(I,J),POINT)
     END DO
  END DO
END FUNCTION TWOSPINOR_POINTWISE_VALUE

SUBROUTINE PRINT2SPINOR(PHI,NUNIT)
  IMPLICIT NONE
  TYPE(twospinor),INTENT(IN) :: PHI
  INTEGER,INTENT(IN) :: NUNIT

  INTEGER :: I,J,K

  DO K=1,2
     WRITE(NUNIT,*)'component #',K
     IF (PHI%nbrofcontractions(K)==0) THEN
        WRITE(NUNIT,*)' no contraction'
     ELSE
        WRITE(NUNIT,*)' number of contractions:',PHI%nbrofcontractions(K)
        DO I=1,PHI%nbrofcontractions(K)
           WRITE(NUNIT,*)'  contraction #',I
           WRITE(NUNIT,*)'    coefficient:',PHI%coefficients(K,I)
           WRITE(NUNIT,*)'    number of gaussian primitives:',PHI%contractions(K,I)%nbrofexponents
           WRITE(NUNIT,*)'    common monomial:',PHI%contractions(K,I)%monomialdegree
           WRITE(NUNIT,*)'    center:',PHI%contractions(K,I)%center
           DO J=1,PHI%contractions(K,I)%nbrofexponents
              WRITE(NUNIT,*)'      gaussian primitive #',J
              WRITE(NUNIT,*)'        exponent:',PHI%contractions(K,I)%exponents(J)
              WRITE(NUNIT,*)'        coefficient:',PHI%contractions(K,I)%coefficients(J)
           END DO
        END DO
     END IF
  END DO
END SUBROUTINE PRINT2SPINOR
END MODULE

MODULE scf_parameters
! number of different SCF algorithms to be used
  INTEGER :: NBALG
! SCF algorithm index (1: Roothaan, 2: level-shifting, 3: DIIS, 4: ODA (non-relativistic case only), 5: \'Eric S\'er\'e's (relativistic case only))
  INTEGER,DIMENSION(5) :: ALG
! threshold for numerical convergence
  DOUBLE PRECISION :: TRSHLD
! maximum number of iterations allowed
  INTEGER :: MAXITR
! flag for the direct computation of the bielectronic integrals
  LOGICAL :: DIRECT
! flag for the "semi-direct" computation of the bielectronic integrals (relativistic case only)
! Note: the case is considered as a (DIRECT==.FALSE.) subcase: GBF bielectronic integrals are precomputed and kept in memory, 2-spinor bielectronic integrals being computed "directly" using these values afterwards.
  LOGICAL :: SEMIDIRECT
! flag for the storage of the computed bielectronic integrals (and/or their list) on disk or in memory
  LOGICAL :: USEDISK
! flag for the use of the SS-bielectronic integrals
  LOGICAL :: SSINTEGRALS
! flag for the use of the LS-bielectronic integrals
  LOGICAL :: SLINTEGRALS
! flag to resume a previous computation (using EIG and EIGVEC)
  LOGICAL :: RESUME
CONTAINS

SUBROUTINE SETUP_SCF
  !$ USE omp_lib
  USE case_parameters ; USE setup_tools
  IMPLICIT NONE
  CHARACTER :: METHOD
  CHARACTER(4) :: CHAR
  INTEGER :: I,INFO,DIISMXSET,STAT,NUMBER_OF_THREADS
  DOUBLE PRECISION :: LSSHIFT,DIISRCP

  OPEN(100,FILE=SETUP_FILE,STATUS='old',ACTION='read')
  CALL LOOKFOR(100,'## SCF PARAMETERS',INFO)
  IF (INFO/=0) STOP' The SCF parameters are not given.'
  READ(100,'(7/,i1)') NBALG
  DO I=1,NBALG
     READ(100,'(i1)') ALG(I)
     IF (RELATIVISTIC.AND.(ALG(I)==4)) THEN
        STOP' The Optimal Damping Algorithm is intended for the non-relativistic case only.'
     ELSE IF ((.NOT.RELATIVISTIC).AND.(ALG(I)==5)) THEN
        STOP' ES''s algorithm is intended for the relativistic case only.'
     END IF
  END DO
  READ(100,*) TRSHLD
  WRITE(*,*)'Threshold =',TRSHLD
  READ(100,'(i5)') MAXITR
  WRITE(*,*)'Maximum number of iterations =',MAXITR
  READ(100,'(a3)') CHAR
  IF (RELATIVISTIC) THEN
     IF (CHAR=='DIR') THEN
        DIRECT=.TRUE.
        WRITE(*,'(a)')' Direct computation of the bielectronic integrals'
! check if the list of the bielectronic integrals is stored on disk or in memory
        READ(100,'(a4)') CHAR
        IF (CHAR=='DISK') THEN
           USEDISK=.TRUE.
        ELSE
           USEDISK=.FALSE.
        END IF
     ELSE IF (CHAR=='NOT') THEN
        DIRECT=.FALSE. ; SEMIDIRECT=.FALSE.
        READ(100,'(a4)') CHAR
! the list of the bielectronic integrals is stored in the same way as the integrals are (on disk or in memory)
        IF (CHAR=='DISK') THEN
           USEDISK=.TRUE.
           WRITE(*,'(a)')' Computed bielectronic integrals stored on disk'
        ELSE
           USEDISK=.FALSE.
           WRITE(*,'(a)')' Computed bielectronic integrals stored in memory'
        END IF
     ELSE IF (CHAR=='SEM') THEN
        DIRECT=.FALSE. ; SEMIDIRECT=.TRUE.
        WRITE(*,'(a)')' "Semi-direct" computation of the bielectronic integrals'
! check if the list of the bielectronic integrals is stored on disk or in memory
        READ(100,'(a4)') CHAR
        IF (CHAR=='DISK') THEN
           USEDISK=.TRUE.
        ELSE
           USEDISK=.FALSE.
        END IF
     ELSE
       WRITE(*,*)'Subroutine SETUP_SCF: unknown type of computation for bielectronic integrals.'
       STOP
     END IF
     IF (MODEL==3) THEN
! CGHF is a subcase of the relativistic case (for the moment...)
        SSINTEGRALS=.FALSE.
        SLINTEGRALS=.FALSE.
     ELSE
        REWIND(100)
        CALL LOOKFOR(100,'NOSS',INFO)
        IF (INFO==0) THEN
           SSINTEGRALS=.FALSE.
           WRITE(*,'(a)')' (SS-integrals are not used in the computation)'
        ELSE
           SSINTEGRALS=.TRUE.
        END IF
        REWIND(100)
        CALL LOOKFOR(100,'NOSL',INFO)
        IF (INFO==0) THEN
           SLINTEGRALS=.FALSE.
           WRITE(*,'(a)')' (SL-integrals are not used in the computation)'
        ELSE
           SLINTEGRALS=.TRUE.
        END IF
     END IF
  ELSE
     IF (CHAR=='DIR') THEN
        DIRECT=.TRUE.
        WRITE(*,'(a)')' Direct computation of the bielectronic integrals'
! the list of the bielectronic integrals is stored in memory
     ELSE IF (CHAR=='NOT') THEN
        DIRECT=.FALSE.
! the list of the bielectronic integrals is stored in the same way as the integrals are (on disk or in memory)
        READ(100,'(a4)') CHAR
        IF (CHAR=='DISK') THEN
           USEDISK=.TRUE.
           WRITE(*,'(a)')' Computed bielectronic integrals stored on disk'
        ELSE
           USEDISK=.FALSE.
           WRITE(*,'(a)')' Computed bielectronic integrals stored in memory'
        END IF
     ELSE
       WRITE(*,*)'Subroutine SETUP_SCF: unknown type of computation for bielectronic integrals.'
       STOP
     END IF
  END IF
! additional verifications on the parameters of the algorithms
  DO I=1,NBALG
     IF (ALG(I)==2) THEN
        REWIND(100)
        CALL LOOKFOR(100,'LEVEL-SHIFTING ALGORITHM PARAMETERS',INFO)
        IF (INFO/=0) GO TO 20
        READ(100,'(/,f16.8)',ERR=20,END=20) LSSHIFT
        IF (LSSHIFT<=0.) GO TO 21
     ELSE IF (ALG(I)==3) THEN
        REWIND(100)
        CALL LOOKFOR(100,'DIIS ALGORITHM PARAMETERS',INFO)
        IF (INFO/=0) GO TO 30
        READ(100,'(/,i2)',ERR=30,END=30) DIISMXSET
        IF (DIISMXSET==0) THEN
! The DIIS algorithm with restarts is used: read the parameter used to check for restart
           READ(100,'(f16.8)',ERR=32,END=32) DIISRCP
           IF (DIISRCP>=1.) GO TO 33
        ELSE IF (DIISMXSET<2) THEN
           GO TO 31
        END IF
     ELSE IF (ALG(I)==5) THEN
        REWIND(100)
        CALL LOOKFOR(100,'SERE''S ALGORITHM PARAMETERS',INFO)
        IF (INFO/=0) GO TO 50
        READ(100,'(/,a)',ERR=50,END=50)METHOD
        IF ((METHOD/='D').AND.(METHOD/='S').AND.(METHOD/='N')) GO TO 50
     END IF
  END DO
  
  RESUME=.FALSE.
  CALL LOOKFOR(100,'RESUME',INFO)
  IF(INFO == 0) THEN
     READ(100,'(a)')CHAR
     IF(CHAR == 'YES') THEN
        RESUME = .TRUE.
     END IF
  END IF
! determination of the number of threads to be used by OpenMP
  !$ CALL LOOKFOR(100,'PARALLELIZATION PARAMETERS',INFO)
  !$ IF (INFO==0) THEN
  !$    READ(100,'(/,i3)')NUMBER_OF_THREADS
  !$    CALL OMP_SET_NUM_THREADS(NUMBER_OF_THREADS)
  !$ END IF
  !$ WRITE(*,'(a,i2,a)') ' The number of thread(s) to be used is ',OMP_GET_MAX_THREADS(),'.'
  CLOSE(100)
  RETURN
! MESSAGES
20 STOP'No valid shift parameter has been given for the level-shifting algorithm.'
21 STOP'The shift parameter for the level-shifting algorithm must be positive.'
30 STOP'No valid maximum simplex dimension has been given for the DIIS algorithm.'
31 STOP'The maximum simplex dimension for the DIIS algorithm must be at least equal to two.'
32 STOP'No valid parameter to check for restart was given for the DIIS algorithm with restarts.'
33 STOP'The parameter to check for restart for the DIIS algorithm with restarts must be smaller than 1.'
50 STOP'No valid method given for the computation of $Theta$ in ES''s algorithm.'
END SUBROUTINE SETUP_SCF
END MODULE
