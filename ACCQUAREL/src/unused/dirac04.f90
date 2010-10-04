! routines relatives a l'utilisation du code DIRAC04
MODULE dirac04
CONTAINS

SUBROUTINE FORMOM(POM,NBAST)
! Lecture de la matrice d'overlap issue de DIRAC04 et stockage sous forme
! compacte triangulaire superieure
  IMPLICIT NONE
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2) :: POM
  INTEGER :: NBAST,N,I,J,K,SHIFTK
  DOUBLE PRECISION :: VALUE

  N=NBAST/2
  POM=(0.D0,0.D0)
  OPEN(UNIT=40,STATUS='OLD',ACTION='READ')
  DO I=1,N*(N+1)/2
     READ(40,'(E22.14)')VALUE
     POM(I)=DCMPLX(VALUE,0.D0)
  END DO
  CLOSE(40)

  K=0 ; SHIFTK=N*(N+1)/2
  DO J=1,N
     SHIFTK=SHIFTK+N
     DO I=1,J
        K=K+1 ; SHIFTK=SHIFTK+1
        POM(SHIFTK)=POM(K)
     END DO
  END DO
END SUBROUTINE

SUBROUTINE FORMOEFM(POEFM,NBAST)
! Lecture de la matrice de Fock uni-electronique issue de DIRAC04 et stockage
! sous forme compacte triangulaire superieure
  IMPLICIT NONE
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2) :: POEFM
  INTEGER :: NBAST,N,I,J,K,IDUMMY
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: QOEFM

  N=NBAST/2
! La matrice est sous forme quaternionique dans DIRAC04 : on lit quatre matrices reelles de taille N*N
  ALLOCATE(QOEFM(1:4,1:N,1:N))
  OPEN(UNIT=30,STATUS='OLD',ACTION='READ')
  DO K=1,4
     READ(30,'(I4)')IDUMMY
     DO J=1,N
        DO I=1,N
           READ(30,'(E22.14)')QOEFM(K,I,J)
        END DO
     END DO
  END DO
  CLOSE(30)
  POEFM=(0.D0,0.D0)
  DO J=1,N
     DO I=1,J
        POEFM(I+(J-1)*J/2)=DCMPLX(QOEFM(1,I,J),0.D0)+DCMPLX(0.D0,QOEFM(2,I,J))
        POEFM(I+(J+N-1)*(J+N)/2)=DCMPLX(0.D0,QOEFM(4,I,J))+DCMPLX(QOEFM(3,I,J),0.D0)
        POEFM(J+(I+N-1)*(I+N)/2)=DCMPLX(0.D0,QOEFM(4,J,I))+DCMPLX(QOEFM(3,J,I),0.D0)
        POEFM(I+N+(J+N-1)*(J+N)/2)=DCMPLX(QOEFM(1,I,J),0.D0)-DCMPLX(0.D0,QOEFM(2,I,J))
     END DO
  END DO
END SUBROUTINE

SUBROUTINE FORMTEFM(PTEFM,NBAST,PDM,LUNIT)
! Construction de la matrice de Fock bi-electronique (a partir de la matrice de densite
! et des integrales doubles issues de DIRAC04) et stockage sous forme compacte
! triangulaire superieure
  IMPLICIT NONE
  DOUBLE COMPLEX,DIMENSION(NBAST*(NBAST+1)/2) :: PTEFM,PDM
  DOUBLE COMPLEX,DIMENSION(NBAST,NBAST) :: TEFM,DM
  INTEGER :: NBAST,LUNIT,I,J,K,L,IJ,IOS
  DOUBLE PRECISION :: VALUE

  IJ=0
  DO J=1,NBAST
     DO I=1,J
        IJ=IJ+1
        DM(I,J)=PDM(IJ)
        IF (I.NE.J) DM(J,I)=CONJG(DM(I,J))
     END DO
  END DO
  TEFM=(0.D0,0.D0)
  OPEN(UNIT=LUNIT,STATUS='OLD',ACTION='READ')
  IOS=0
  DO WHILE (.TRUE.)
     READ(LUNIT,'(4I4,2X,E22.14)',IOSTAT=IOS)I,J,K,L,VALUE
     IF (IOS.NE.0) GO TO 1
! 1 valeur pour les 4 indices
     IF ((I.EQ.J).AND.(J.EQ.K).AND.(K.EQ.L)) THEN
        CALL DOIIII(TEFM,DM,NBAST,VALUE,I)
! 2 valeurs distinctes pour les 4 indices
     ELSE IF ((I.GT.J).AND.(J.EQ.K).AND.(K.EQ.L)) THEN
        CALL DOIJJJ(TEFM,DM,NBAST,VALUE,I,J)
     ELSE IF ((I.EQ.J).AND.(J.EQ.K).AND.(K.GT.L)) THEN
        CALL DOIJJJ(TEFM,DM,NBAST,VALUE,L,I)
     ELSE IF ((I.EQ.J).AND.(J.GT.K).AND.(K.EQ.L)) THEN
        CALL DOIIKK(TEFM,DM,NBAST,VALUE,I,K)
     ELSE IF ((I.EQ.K).AND.(K.GT.J).AND.(J.EQ.L)) THEN
        CALL DOIJIJ(TEFM,DM,NBAST,VALUE,I,J)
! 3 valeurs distinctes pour les 4 indices
     ELSE IF ((I.EQ.K).AND.(K.GT.J).AND.(J.GT.L)) THEN
        CALL DOIJIL(TEFM,DM,NBAST,VALUE,I,J,L)
     ELSE IF ((I.GT.J).AND.(J.EQ.K).AND.(K.GT.L)) THEN
        CALL DOIJIL(TEFM,DM,NBAST,VALUE,J,I,L)
     ELSE IF ((I.GT.K).AND.(K.GT.J).AND.(J.EQ.L)) THEN
        CALL DOIJIL(TEFM,DM,NBAST,VALUE,J,I,K)
     ELSE IF ((I.GT.J).AND.(I.GT.K).AND.(K.EQ.L)) THEN
        CALL DOIJKK(TEFM,DM,NBAST,VALUE,I,J,K)
     ELSE IF ((I.EQ.J).AND.(J.GT.K).AND.(K.GT.L)) THEN
        CALL DOIJKK(TEFM,DM,NBAST,VALUE,K,L,I)
! 4 valeurs distinctes pour les 4 indices
     ELSE IF (    ((I.GT.J).AND.(J.GT.K).AND.(K.GT.L)) &
              .OR.((I.GT.K).AND.(K.GT.J).AND.(J.GT.L)) &
              .OR.((I.GT.K).AND.(K.GT.L).AND.(L.GT.J))) THEN
        CALL DOIJKL(TEFM,DM,NBAST,VALUE,I,J,K,L)
     END IF
  END DO
1 CLOSE(LUNIT)
  IJ=0
  DO J=1,NBAST
     DO I=1,J
        IJ=IJ+1
        PTEFM(IJ)=TEFM(I,J)
     END DO
  END DO
END SUBROUTINE
END MODULE
