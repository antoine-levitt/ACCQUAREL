MODULE basis
  USE iso_c_binding
  INTEGER,DIMENSION(4),PARAMETER :: NGPMUSC=(/1,3,6,10/)
  INTEGER,DIMENSION(4),PARAMETER :: NGPMLSC=(/3,7,13,21/)

INTERFACE FORMBASIS
  MODULE PROCEDURE FORMBASIS_relativistic,FORMBASIS_nonrelativistic
END INTERFACE
CONTAINS

SUBROUTINE FORMBASIS_relativistic(PHI,NBAS,GBF,NGBF)
! Subroutine that builds the 2-spinor basis functions and related contracted cartesian Gaussian-type orbital functions for the molecular system considered, from the coefficients read in a file.
  USE basis_parameters ; USE data_parameters ; USE case_parameters ; USE mathematical_functions ; USE constants
  IMPLICIT NONE
  TYPE(twospinor),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: PHI
  INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: NBAS
  TYPE(gaussianbasisfunction),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: GBF
  INTEGER,DIMENSION(2),INTENT(OUT) :: NGBF

  INTEGER,DIMENSION(NBN) :: HAQN
  INTEGER,DIMENSION(MAQN,NBN) :: NOP,NOC
  DOUBLE PRECISION :: NCOEF
  DOUBLE PRECISION,DIMENSION(MNOP,MAQN,NBN) :: ALPHA
  DOUBLE PRECISION,DIMENSION(MNOP,MNOC,MAQN,NBN) :: CCOEF
  INTEGER,DIMENSION(2,NBN) :: NBASN,NGBFN
  INTEGER :: I,J,K,L,M,IUA,IUB,ILA,ILB,IP,IPIC,ICOEF,RNOC,NBOPIC,IUGBF,ILGBF,IDX

  ALLOCATE(NBAS(1:2))
  DO M=1,NBN
     CALL READBASIS(Z(M),HAQN(M),NOP(:,M),NOC(:,M),ALPHA(:,:,M),CCOEF(:,:,:,M),NBASN(:,M),NGBFN(:,M))
  END DO
  NBAS=SUM(NBASN,DIM=2) ; NGBF=SUM(NGBFN,DIM=2)
  ALLOCATE(PHI(1:SUM(NBAS)),GBF(1:SUM(NGBF)))
! Initialization of the indices for the upper 2-spinor basis functions and related gaussian basis functions
  IUA=0 ; IUB=NBAS(1)/2 ; IUGBF=0
! Indices for the lower 2-spinor basis functions and related gaussian basis functions
  ILA=NBAS(1) ; ILB=NBAS(1)+NBAS(2)/2 ; ILGBF=NGBF(1)
  DO M=1,NBN
! Loop on the nuclei of the molecular system
     DO I=1,HAQN(M)
        IF (UNCONT) THEN
! uncontracted basis
           DO J=1,NOP(I,M)
              IF (KINBAL) THEN
! Construction of the gaussian primitives related to the lower 2-spinor basis functions when a kinetic balance scheme is used
! note: the coefficients in the contractions are derived from the RKB scheme and put to 1 if the UKB scheme is used.
                 DO K=1,NGPMLSC(I)
                    ILGBF=ILGBF+1
                    GBF(ILGBF)%center=CENTER(:,M)
                    GBF(ILGBF)%center_id=M
                    GBF(ILGBF)%nbrofexponents=1
                    GBF(ILGBF)%exponents(1)=ALPHA(J,I,M)
! note: see how the UKB scheme is to be interpreted before modifying and uncommenting.
!                    IF (((I==2).AND.(K==7)).OR.((I==3).AND.(K>=11)).OR.((I==4).AND.(K>=16))) THEN
                       GBF(ILGBF)%coefficients(1)=1.D0 
!                    ELSE
!                       GBF(ILGBF)%coefficients(1)=2.D0*ALPHA(J,I,M)
!                    END IF
                    GBF(ILGBF)%monomialdegree=MDLSC(I,K)
                    IF (UKB) THEN
! Construction of the lower 2-spinor basis functions (unrestricted kinetic balance scheme)
! (nota bene: this scheme is not really functional due to linear depencies which have yet to be eliminated)
                       GBF(ILGBF)%coefficients(1)=1.D0
                       ILA=ILA+1 ; ILB=ILB+1 
                       PHI(ILA)%nbrofcontractions=(/1,0/) ; PHI(ILB)%nbrofcontractions=(/0,1/) 
                       PHI(ILA)%contractions(1,1)=GBF(ILGBF) ; PHI(ILB)%contractions(2,1)=GBF(ILGBF)
                       PHI(ILA)%contidx(1,1)=ILGBF ; PHI(ILB)%contidx(2,1)=ILGBF
                       PHI(ILA)%coefficients(1,1)=(1.D0,0.D0) ; PHI(ILB)%coefficients(2,1)=(1.D0,0.D0)
                    END IF
                 END DO
              END IF
              DO K=1,NGPMUSC(I)
! Construction of the uncontracted GBF functions related to the upper 2-spinor basis functions
                 IUGBF=IUGBF+1
                 GBF(IUGBF)%center=CENTER(:,M)
                 GBF(IUGBF)%center_id=M
                 GBF(IUGBF)%nbrofexponents=1
                 GBF(IUGBF)%exponents(1)=ALPHA(J,I,M)
                 GBF(IUGBF)%coefficients(1)=1.D0
                 GBF(IUGBF)%monomialdegree=MDUSC(I,K)
! Construction of the upper 2-spinor basis functions
                 IUA=IUA+1 ; IUB=IUB+1
                 PHI(IUA)%nbrofcontractions=(/1,0/) ; PHI(IUB)%nbrofcontractions=(/0,1/)
                 PHI(IUA)%contractions(1,1)=GBF(IUGBF) ; PHI(IUB)%contractions(2,1)=GBF(IUGBF)
                 PHI(IUA)%contidx(1,1)=IUGBF ; PHI(IUB)%contidx(2,1)=IUGBF
                 PHI(IUA)%coefficients(1,1)=(1.D0,0.D0) ; PHI(IUB)%coefficients(2,1)=(1.D0,0.D0)
! Construction of the lower 2-spinor basis functions (restricted kinetic balance scheme)
                 ILA=ILA+1 ; ILB=ILB+1
                 IF (KINBAL.AND..NOT.UKB) GO TO 3
1             END DO
           END DO
        ELSE
! a priori contracted basis
           IF (KINBAL.AND.UKB) GO TO 4
           ICOEF=1
           IP=1
           IF (NOC(I,M)==0) THEN
              RNOC=NOP(I,M)
           ELSE
              RNOC=NOC(I,M)
           END IF
           DO J=1,RNOC
              NBOPIC=0
              DO
                 IF (CCOEF(ICOEF,J,I,M)==0.D0) EXIT
                 NBOPIC=NBOPIC+1 
                 IF (ICOEF==NOP(I,M)) EXIT
                 ICOEF=ICOEF+1
              END DO
              IF (KINBAL) THEN
! Construction of the gaussian basis functions related to the lower 2-spinor basis functions when a kinetic balance scheme is used
                 DO K=1,NGPMLSC(I)
                    ILGBF=ILGBF+1
                    GBF(ILGBF)%center=CENTER(:,M)
                    GBF(ILGBF)%center_id=M
                    GBF(ILGBF)%nbrofexponents=NBOPIC
                    IPIC=0
                    DO L=IP,IP+NBOPIC-1
                       IPIC=IPIC+1
                       GBF(ILGBF)%exponents(IPIC)=ALPHA(L,I,M)
! part of the normalization coefficient depending only on the exponent and the monomial total degree
                       NCOEF=SQRT((2.D0*ALPHA(L,I,M)/PI)**1.5D0*(4.D0*ALPHA(L,I,M))**(I-1))
! nota bene: the multiplication of the contraction coefficient by the exponent, when needed, is done here.
                       IF (((I==2).AND.(K==7)).OR.((I==3).AND.(K>=11)).OR.((I==4).AND.(K>=16))) THEN
                          GBF(ILGBF)%coefficients(IPIC)=NCOEF*CCOEF(L,J,I,M)
                       ELSE
                          GBF(ILGBF)%coefficients(IPIC)=ALPHA(L,I,M)*NCOEF*CCOEF(L,J,I,M)
                       END IF
                       GBF(ILGBF)%monomialdegree=MDLSC(I,K)
                    END DO
                 END DO
              END IF
              DO K=1,NGPMUSC(I)
! Construction of the gaussian basis functions related to the upper 2-spinor basis functions
                 IUGBF=IUGBF+1
                 GBF(IUGBF)%center=CENTER(:,M)
                 GBF(IUGBF)%center_id=M
                 GBF(IUGBF)%nbrofexponents=NBOPIC
                 GBF(IUGBF)%monomialdegree=MDUSC(I,K)
                 IPIC=0
                 DO L=IP,IP+NBOPIC-1
                    IPIC=IPIC+1
                    GBF(IUGBF)%exponents(IPIC)=ALPHA(L,I,M)
! normalization coefficient
                    NCOEF=SQRT((2.D0*ALPHA(L,I,M)/PI)**1.5D0*(4.D0*ALPHA(L,I,M))**(I-1)/(DFACT(2*GBF(IUGBF)%monomialdegree(1)-1) &
 &                             *DFACT(2*GBF(IUGBF)%monomialdegree(2)-1)*DFACT(2*GBF(IUGBF)%monomialdegree(3)-1)))
                    GBF(IUGBF)%coefficients(IPIC)=NCOEF*CCOEF(L,J,I,M)
                 END DO
! Construction of the upper 2-spinor basis functions
                 IUA=IUA+1 ; IUB=IUB+1
                 PHI(IUA)%nbrofcontractions=(/1,0/) ; PHI(IUB)%nbrofcontractions=(/0,1/)
                 PHI(IUA)%contractions(1,1)=GBF(IUGBF) ; PHI(IUB)%contractions(2,1)=GBF(IUGBF)
                 PHI(IUA)%contidx(1,1)=IUGBF ; PHI(IUB)%contidx(2,1)=IUGBF
                 PHI(IUA)%coefficients(1,1)=(1.D0,0.D0) ; PHI(IUB)%coefficients(2,1)=(1.D0,0.D0)
! Construction of the lower 2-spinor basis functions (restricted kinetic balance scheme)
                 ILA=ILA+1 ; ILB=ILB+1
                 IF (KINBAL.AND..NOT.UKB) GO TO 3
2             END DO
              IP=IP+NBOPIC
           END DO
        END IF
     END DO
  END DO
  IF (.NOT.KINBAL) THEN
! If no kinetic balance scheme is used, basis functions for the lower 2-spinor are exactly the same as the ones for the upper 2-spinor
     DO I=1,NGBF(2)
! HIS IS NOT OPTIMAL IN TERMS OF COMPUTATION...
        GBF(NGBF(1)+I)=GBF(I)
     END DO
     DO I=1,NBAS(2)
        PHI(NBAS(1)+I)=PHI(I)
     END DO
  END IF
  WRITE(*,'(a,i4)')' Total number of basis functions =',SUM(NBAS)
  RETURN
3 IDX=ILGBF-NGPMLSC(I)
  SELECT CASE (I)
     CASE (1)
        PHI(ILA)%nbrofcontractions=(/1,2/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+3) ; PHI(ILA)%contractions(2,1:2)=(/GBF(IDX+1),GBF(IDX+2)/)
        PHI(ILA)%contidx(1,1)=IDX+3 ; PHI(ILA)%contidx(2,1:2)=(/IDX+1,IDX+2/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:2)=(/(0.D0,2.D0),(-2.D0,0.D0)/)
! NOTE: the part of the normalization coefficient depending only the monomial degrees is equal to 1.
        PHI(ILA)%coefficients(:,:)=C*PHI(ILA)%coefficients(:,:)
     CASE (2)
     SELECT CASE (K)
        CASE (1)
        PHI(ILA)%nbrofcontractions=(/1,3/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+3) ; PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+1),GBF(IDX+2),GBF(IDX+7)/)
        PHI(ILA)%contidx(1,1)=IDX+3 ; PHI(ILA)%contidx(2,1:3)=(/IDX+1,IDX+2,IDX+7/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(0.D0,-1.D0)/)
! Note: the part of the normalization coefficient depending only the monomial degrees is equal to 1.
        PHI(ILA)%coefficients(:,:)=C*PHI(ILA)%coefficients(:,:)
        CASE (2)
        PHI(ILA)%nbrofcontractions=(/1,3/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+5) ; PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+2),GBF(IDX+4),GBF(IDX+7)/)
        PHI(ILA)%contidx(1,1)=IDX+5 ; PHI(ILA)%contidx(2,1:3)=(/IDX+2,IDX+4,IDX+7/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(1.D0,0.D0)/)
! Note: the part of the normalization coefficient depending only the monomial degrees is equal to 1.
        PHI(ILA)%coefficients(:,:)=C*PHI(ILA)%coefficients(:,:)
        CASE (3)
        PHI(ILA)%nbrofcontractions=(/2,2/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+6),GBF(IDX+7)/) ; PHI(ILA)%contractions(2,1:2)=(/GBF(IDX+3),GBF(IDX+5)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+6,IDX+7/) ; PHI(ILA)%contidx(2,1:2)=(/IDX+3,IDX+5/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-1.D0)/) ; PHI(ILA)%coefficients(2,1:2)=(/(0.D0,2.D0),(-2.D0,0.D0)/)
! Note: the part of the normalization coefficient depending only the monomial degrees is equal to 1.
        PHI(ILA)%coefficients(:,:)=C*PHI(ILA)%coefficients(:,:)
     END SELECT
     CASE (3)
     SELECT CASE (K)
        CASE (1)
        PHI(ILA)%nbrofcontractions=(/1,3/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+3) ; PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+1),GBF(IDX+2),GBF(IDX+11)/)
        PHI(ILA)%contidx(1,1)=IDX+3 ; PHI(ILA)%contidx(2,1:3)=(/IDX+1,IDX+2,IDX+11/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(0.D0,-2.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (2)
        PHI(ILA)%nbrofcontractions=(/1,4/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+6) ; PHI(ILA)%contractions(2,1:4)=(/GBF(IDX+2),GBF(IDX+4),GBF(IDX+11),GBF(IDX+12)/)
        PHI(ILA)%contidx(1,1)=IDX+6 ; PHI(ILA)%contidx(2,1:4)=(/IDX+2,IDX+4,IDX+11,IDX+12/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:4)=(/(0.D0,2.D0),(-2.D0,0.D0),(1.D0,0.D0),(0.D0,-1.D0)/)
! Note: the part of the normalization coefficient depending only the monomial degrees is equal to 1.
        PHI(ILA)%coefficients(:,:)=C*PHI(ILA)%coefficients(:,:)
        CASE (3)
        PHI(ILA)%nbrofcontractions=(/2,3/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+5),GBF(IDX+11)/)
        PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+3),GBF(IDX+6),GBF(IDX+13)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+5,IDX+11/) ; PHI(ILA)%contidx(2,1:3)=(/IDX+3,IDX+6,IDX+13/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-1.D0)/)
        PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(0.D0,-1.D0)/)
! Note: the part of the normalization coefficient depending only the monomial degrees is equal to 1.
        PHI(ILA)%coefficients(:,:)=C*PHI(ILA)%coefficients(:,:)
        CASE (4)
        PHI(ILA)%nbrofcontractions=(/1,3/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+8) ; PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+4),GBF(IDX+7),GBF(IDX+12)/)
        PHI(ILA)%contidx(1,1)=IDX+8 ; PHI(ILA)%contidx(2,1:3)=(/IDX+4,IDX+7,IDX+12/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(2.D0,0.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (5)
        PHI(ILA)%nbrofcontractions=(/2,3/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+9),GBF(IDX+12)/)
        PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+6),GBF(IDX+8),GBF(IDX+13)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+9,IDX+12/) ; PHI(ILA)%contidx(2,1:3)=(/IDX+6,IDX+8,IDX+13/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-1.D0)/)
        PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(1.D0,0.D0)/)
! Note: the part of the normalization coefficient depending only the monomial degrees is equal to 1.
        PHI(ILA)%coefficients(:,:)=C*PHI(ILA)%coefficients(:,:)
        CASE (6)
        PHI(ILA)%nbrofcontractions=(/2,2/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+10),GBF(IDX+13)/) ; PHI(ILA)%contractions(2,1:2)=(/GBF(IDX+5),GBF(IDX+9)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+10,IDX+13/) ; PHI(ILA)%contidx(2,1:2)=(/IDX+5,IDX+9/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-2.D0)/) ; PHI(ILA)%coefficients(2,1:2)=(/(0.D0,2.D0),(-2.D0,0.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
     END SELECT
     CASE (4)
     SELECT CASE (K)
        CASE (1)
        PHI(ILA)%nbrofcontractions=(/1,3/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+3) ; PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+1),GBF(IDX+2),GBF(IDX+16)/)
        PHI(ILA)%contidx(1,1)=IDX+3 ; PHI(ILA)%contidx(2,1:3)=(/IDX+1,IDX+2,IDX+16/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(0.D0,-3.D0)/)
        NCOEF=15**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (2)
        PHI(ILA)%nbrofcontractions=(/1,4/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+5)
        PHI(ILA)%contractions(2,1:4)=(/GBF(IDX+2),GBF(IDX+4),GBF(IDX+16),GBF(IDX+17)/)
        PHI(ILA)%contidx(1,1)=IDX+5 ; PHI(ILA)%contidx(2,1:4)=(/IDX+2,IDX+4,IDX+16,IDX+17/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0)
        PHI(ILA)%coefficients(2,1:4)=(/(0.D0,2.D0),(-2.D0,0.D0),(1.D0,0.D0),(0.D0,-2.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (3)
        PHI(ILA)%nbrofcontractions=(/2,3/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+6),GBF(IDX+16)/)
        PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+3),GBF(IDX+5),GBF(IDX+18)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+6,IDX+16/) ; PHI(ILA)%contidx(2,1:3)=(/IDX+3,IDX+5,IDX+18/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-1.D0)/)
        PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(0.D0,-2.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (4)
        PHI(ILA)%nbrofcontractions=(/1,4/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+8) ; PHI(ILA)%contractions(2,1:4)=(/GBF(IDX+4),GBF(IDX+7),GBF(IDX+17),GBF(IDX+19)/)
        PHI(ILA)%contidx(1,1)=IDX+8 ; PHI(ILA)%contidx(2,1:4)=(/IDX+4,IDX+7,IDX+17,IDX+19/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:4)=(/(0.D0,2.D0),(-2.D0,0.D0),(2.D0,0.D0),(0.D0,-1.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (5)
        PHI(ILA)%nbrofcontractions=(/2,3/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+10),GBF(IDX+18)/)
        PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+6),GBF(IDX+9),GBF(IDX+21)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+10,IDX+18/) ; PHI(ILA)%contidx(2,1:3)=(/IDX+6,IDX+9,IDX+21/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-2.D0)/)
        PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(0.D0,-1.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (6)
        PHI(ILA)%nbrofcontractions=(/2,4/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+9),GBF(IDX+17)/)
        PHI(ILA)%contractions(2,1:4)=(/GBF(IDX+5),GBF(IDX+8),GBF(IDX+18),GBF(IDX+20)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+9,IDX+17/) ; PHI(ILA)%contidx(2,1:4)=(/IDX+5,IDX+8,IDX+18,IDX+20/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-1.D0)/)
        PHI(ILA)%coefficients(2,1:4)=(/(0.D0,2.D0),(-2.D0,0.D0),(1.D0,0.D0),(0.D0,-1.D0)/)
! Note: the part of the normalization coefficient depending only the monomial degrees is equal to 1.
        PHI(ILA)%coefficients(:,:)=C*PHI(ILA)%coefficients(:,:)
        CASE (7)
        PHI(ILA)%nbrofcontractions=(/1,3/)
        PHI(ILA)%contractions(1,1)=GBF(IDX+12) ; PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+10),GBF(IDX+11),GBF(IDX+19)/)
        PHI(ILA)%contidx(1,1)=IDX+12 ; PHI(ILA)%contidx(2,1:3)=(/IDX+10,IDX+11,IDX+19/)
        PHI(ILA)%coefficients(1,1)=(0.D0,2.D0) ; PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(3.D0,0.D0)/)
        NCOEF=15**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (8)
        PHI(ILA)%nbrofcontractions=(/2,3/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+13),GBF(IDX+19)/)
        PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+8),GBF(IDX+12),GBF(IDX+20)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+13,IDX+19/) ; PHI(ILA)%contidx(2,1:3)=(/IDX+8,IDX+12,IDX+20/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-1.D0)/)
        PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(2.D0,0.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (9)
        PHI(ILA)%nbrofcontractions=(/2,3/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+14),GBF(IDX+20)/)
        PHI(ILA)%contractions(2,1:3)=(/GBF(IDX+9),GBF(IDX+13),GBF(IDX+21)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+14,IDX+20/) ; PHI(ILA)%contidx(2,1:3)=(/IDX+9,IDX+13,IDX+21/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-2.D0)/)
        PHI(ILA)%coefficients(2,1:3)=(/(0.D0,2.D0),(-2.D0,0.D0),(1.D0,0.D0)/)
        NCOEF=3**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
        CASE (10)
        PHI(ILA)%nbrofcontractions=(/2,2/)
        PHI(ILA)%contractions(1,1:2)=(/GBF(IDX+15),GBF(IDX+21)/) ; PHI(ILA)%contractions(2,1:2)=(/GBF(IDX+10),GBF(IDX+14)/)
        PHI(ILA)%contidx(1,1:2)=(/IDX+15,IDX+21/) ; PHI(ILA)%contidx(2,1:2)=(/IDX+10,IDX+14/)
        PHI(ILA)%coefficients(1,1:2)=(/(0.D0,2.D0),(0.D0,-3.D0)/) ; PHI(ILA)%coefficients(2,1:2)=(/(0.D0,2.D0),(-2.D0,0.D0)/)
        NCOEF=15**(-0.5D0)
        PHI(ILA)%coefficients(:,:)=C*NCOEF*PHI(ILA)%coefficients(:,:)
     END SELECT
  END SELECT
  PHI(ILB)%nbrofcontractions(1)=PHI(ILA)%nbrofcontractions(2)
  PHI(ILB)%contractions(1,:)=PHI(ILA)%contractions(2,:)
  PHI(ILB)%contidx(1,:)=PHI(ILA)%contidx(2,:)
  PHI(ILB)%coefficients(1,:)=-CONJG(PHI(ILA)%coefficients(2,:))
  PHI(ILB)%nbrofcontractions(2)=PHI(ILA)%nbrofcontractions(1)
  PHI(ILB)%contractions(2,:)=PHI(ILA)%contractions(1,:)
  PHI(ILB)%contidx(2,:)=PHI(ILA)%contidx(1,:)
  PHI(ILB)%coefficients(2,:)=-PHI(ILA)%coefficients(1,:)
  IF (UNCONT) THEN
     GO TO 1
  ELSE
     GO TO 2
  END IF
4 IF (KINBAL.AND.UKB) STOP 'Option not implemented yet!'
END SUBROUTINE FORMBASIS_relativistic

SUBROUTINE FORMBASIS_nonrelativistic(PHI,NBAS)
! Subroutine that builds the gaussian basis functions for the molecular system considered, from coefficients either read in a file (basis set from the existing literature) or depending on some given parameters (even-tempered basis set).
  USE basis_parameters ; USE data_parameters ; USE mathematical_functions ; USE constants
  IMPLICIT NONE
  TYPE(gaussianbasisfunction),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: PHI
  INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: NBAS

  INTEGER :: I,J,K,L,M,IBF,IP,IPIC,ICOEF,NBOPIC
  INTEGER,DIMENSION(NBN) :: HAQN
  INTEGER,DIMENSION(1,NBN) :: NBASN
  INTEGER,DIMENSION(MAQN,NBN) :: NOP,NOC
  DOUBLE PRECISION :: NCOEF
  DOUBLE PRECISION,DIMENSION(MNOP,MAQN,NBN) :: ALPHA
  DOUBLE PRECISION,DIMENSION(MNOP,MNOC,MAQN,NBN) :: CCOEF

  ALLOCATE(NBAS(1:1))
  IF (LIBRARY) THEN
     DO M=1,NBN
        CALL READBASIS(Z(M),HAQN(M),NOP(:,M),NOC(:,M),ALPHA(:,:,M),CCOEF(:,:,:,M),NBASN(:,M))
     END DO
     NBAS=SUM(NBASN,DIM=2)
     ALLOCATE(PHI(1:SUM(NBAS)))
     IBF=0
     DO M=1,NBN
! Loop on the nuclei of the molecular system
        DO I=1,HAQN(M)
           IF (UNCONT) THEN
              DO J=1,NOP(I,M)
                 DO K=1,NGPMUSC(I)
                    IBF=IBF+1
                    PHI(IBF)%center=CENTER(:,M)
                    PHI(IBF)%center_id=M
                    PHI(IBF)%nbrofexponents=1
                    PHI(IBF)%exponents(1)=ALPHA(J,I,M)
                    PHI(IBF)%coefficients(1)=1.D0
                    PHI(IBF)%monomialdegree=MDUSC(I,K)
                 END DO
              END DO
           ELSE
              ICOEF=1
              IP=1
              DO J=1,NOC(I,M)
                 NBOPIC=0
                 DO
                    IF (CCOEF(ICOEF,J,I,M)==0.D0) EXIT
                    NBOPIC=NBOPIC+1
                    IF (ICOEF==NOP(I,M)) EXIT
                    ICOEF=ICOEF+1
                 END DO
                 DO K=1,NGPMUSC(I)
                    IBF=IBF+1
                    PHI(IBF)%center=CENTER(:,M)
                    PHI(IBF)%center_id=M
                    PHI(IBF)%nbrofexponents=NBOPIC
                    PHI(IBF)%monomialdegree=MDUSC(I,K)
                    IPIC=0
                    DO L=IP,IP+NBOPIC-1
                       IPIC=IPIC+1
                       PHI(IBF)%exponents(IPIC)=ALPHA(L,I,M)
! normalization coefficient
                       NCOEF=SQRT((2.D0*ALPHA(L,I,M)/PI)**1.5D0*(4.D0*ALPHA(L,I,M))**(I-1)/(DFACT(2*PHI(IBF)%monomialdegree(1)-1) &
 &                                *DFACT(2*PHI(IBF)%monomialdegree(2)-1)*DFACT(2*PHI(IBF)%monomialdegree(3)-1)))
                       PHI(IBF)%coefficients(IPIC)=NCOEF*CCOEF(L,J,I,M)
                    END DO
                 END DO
                 IP=IP+NBOPIC
              END DO
           END IF
        END DO
     END DO
  ELSE
     NBAS=NUMBER_OF_TERMS
     open(22)
     ALLOCATE(PHI(1:SUM(NBAS)))
     DO I=1,SUM(NBAS)
        PHI(I)%center=(/0.D0,0.D0,0.D0/)
        PHI(I)%nbrofexponents=1
        PHI(I)%exponents(1)=FIRST_TERM*COMMON_RATIO**(I-1)
        write(22,*)FIRST_TERM*COMMON_RATIO**(I-1)
        PHI(I)%coefficients(1)=1.D0
        PHI(I)%monomialdegree=(/0,0,0/)
     END DO
     close(22)
  END IF
  WRITE(*,'(a,i4)')' Total number of basis functions =',SUM(NBAS)
END SUBROUTINE FORMBASIS_nonrelativistic

SUBROUTINE READBASIS(Z,HAQN,NOP,NOC,ALPHA,CCOEF,NBAS,NGBF)
! Subroutine that reads (for a given chemical element) the coefficients of the functions of a chosen basis in the DALTON library and computes the number of basis functions (and related cartesian gaussian-type orbital functions) for the components of the upper and lower 2-spinors with respect to the indicated options (contracted basis or not, use of the "kinetic balance" process or not, etc...)
! Note: the name (and path) of the file containing the basis name is stored in the variable BASISFILE from the module basis_parameters.
! Z : atomic number of the element under consideration
! MAQN : maximum angular quantum number + 1 (= number of function types) allowed (s=1, p=2, d=3, etc...)
! MNOP : maximum number of primitives of any given type allowed
! MNOC : maximum number of contractions of any given type allowed
! HAQN : highest angular quantum number of the element basis functions + 1
! NOP : array containing the number of primitives of each type of basis function
! NOC : array containing the number of contractions for each type of basis function
! ALPHA : array containing the exponents of the primitives for each type of basis functions
! CCOEF : array containing the contraction coefficients for each type of basis functions
  USE case_parameters ; USE basis_parameters
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: Z
  INTEGER,INTENT(OUT) :: HAQN
  INTEGER,DIMENSION(MAQN),INTENT(OUT) :: NOP,NOC
  DOUBLE PRECISION,DIMENSION(MNOP,MAQN),INTENT(OUT) :: ALPHA
  DOUBLE PRECISION,DIMENSION(MNOP,MNOC,MAQN),INTENT(OUT) :: CCOEF
  INTEGER,DIMENSION(:),INTENT(OUT) :: NBAS
  INTEGER,DIMENSION(2),INTENT(OUT),OPTIONAL :: NGBF

  INTEGER :: LUBAS,PRINUM,CONNUM,IDUMMY,I
  LOGICAL :: LEX,CHKNELM

  LUBAS=20
  INQUIRE(FILE=BASISFILE,EXIST=LEX)
  IF (LEX) THEN
! Open the file containing the basis
     OPEN(UNIT=LUBAS,FILE=BASISFILE,ACTION='READ')
! Search for the right atomic element in this file
     CALL FINDELM(Z,BASISFILE,LUBAS)
     NOP=0 ; NOC=0 ; HAQN=0
! Find the number of primitives and contractions of each type and read them
1    CALL READPCN(CHKNELM,PRINUM,CONNUM,IDUMMY,BASISFILE,LUBAS)
     IF (UNCONT) CONNUM=-CONNUM
     IF (.NOT.CHKNELM) THEN
        HAQN=HAQN+1
        IF (HAQN>MAQN) GO TO 3
        NOP(HAQN)=PRINUM ; NOC(HAQN)=CONNUM
! Fill ALPHA and CCOEF by reading the corresponding exponents and contraction coefficients
        CALL READECC(MNOP,MNOC,NOP(HAQN),NOC(HAQN),ALPHA(:,HAQN),CCOEF(:,:,HAQN),BASISFILE,LUBAS)
        GO TO 1
     END IF
     CLOSE(LUBAS)
  ELSE
     GO TO 2
  END IF
! Computation of the number of basis functions
  IF (RELATIVISTIC) THEN
! relativistic case
     NBAS=0 ; NGBF=0
     DO I=1,HAQN
        IF (UNCONT) THEN
! uncontracted basis
           NBAS(1)=NBAS(1)+NOP(I)*NGPMUSC(I) 
           NGBF(1)=NGBF(1)+NOP(I)*NGPMUSC(I)
           IF (KINBAL) THEN
              NGBF(2)=NGBF(2)+NOP(I)*NGPMLSC(I)
              IF (UKB) NBAS(2)=NBAS(2)+NOP(I)*NGPMLSC(I)
           END IF
        ELSE
! a priori contracted basis
           IF (NOC(I)/=0) THEN
              NBAS(1)=NBAS(1)+NOC(I)*NGPMUSC(I)
              NGBF(1)=NGBF(1)+NOC(I)*NGPMUSC(I)
           ELSE
              NBAS(1)=NBAS(1)+NOP(I)*NGPMUSC(I)
              NGBF(1)=NGBF(1)+NOP(I)*NGPMUSC(I)
           END IF
           IF (KINBAL) THEN
              IF (NOC(I)/=0) THEN
                 NGBF(2)=NGBF(2)+NOC(I)*NGPMLSC(I)
                 IF (UKB) NBAS(2)=NBAS(2)+NOC(I)*NGPMLSC(I)
              ELSE
                 NGBF(2)=NGBF(2)+NOP(I)*NGPMLSC(I)
                 IF (UKB) NBAS(2)=NBAS(2)+NOP(I)*NGPMLSC(I)
              END IF
           END IF
        END IF
     END DO
     IF (KINBAL.AND..NOT.UKB) NBAS(2)=NBAS(1)
     IF (.NOT.KINBAL) THEN
        NBAS(2)=NBAS(1) ; NGBF(2)=NGBF(1)
     END IF
     NBAS=2*NBAS
  ELSE
! non-relativistic case
     NBAS(1)=0
     DO I=1,HAQN
        IF (UNCONT) THEN
! uncontracted basis
           NBAS(1)=NBAS(1)+NOP(I)*NGPMUSC(I)
        ELSE
! a priori contracted basis
           IF (NOC(I)/=0) THEN
              NBAS(1)=NBAS(1)+NOC(I)*NGPMUSC(I)
           ELSE
              NBAS(1)=NBAS(1)+NOP(I)*NGPMUSC(I)
           END IF
        END IF
     END DO
  END IF
  RETURN
! MESSAGES
2 WRITE(*,'(a)')' Subroutine READBASIS: the file containing the basis could not be found.'
  GO TO 4
3 WRITE(*,'(a)')' Subroutine READBASIS: this type of gaussian basis function is not supported!'
4 STOP
END SUBROUTINE READBASIS

SUBROUTINE FINDELM(Z,BASNAM,LUBAS)
  IMPLICIT NONE
  INTEGER :: Z,LUBAS
  CHARACTER(*) :: BASNAM
  CHARACTER(20) :: STRING
  CHARACTER :: CHDUMMY
  INTEGER :: IDUMMY,IOERR

1 READ(LUBAS,'(a20)',IOSTAT=IOERR,ERR=2,END=3)STRING
  READ(STRING,'(a1)')CHDUMMY
  IF ((CHDUMMY=='a').OR.(CHDUMMY=='A')) THEN
     READ(STRING,'(BN,a1,i4)')CHDUMMY,IDUMMY
     IF (Z==IDUMMY) THEN
       RETURN
     ELSE
        GO TO 1
     END IF
  ELSE
     GO TO 1
  END IF

! MESSAGES
2 WRITE(*,'(a,a,a,i2,a)')' Subroutine FINDELM: an error ocurred during the reading of the file containing the basis ', &
 & BASNAM,' (IOSTAT code=',IOERR,') '
  STOP
3 WRITE(*,'(a,i2,a,a,a)')' Subroutine FINDELM: the basis functions relative to atomic element number ',Z, &
 & ' could not be found in the basis ',BASNAM,'.'
  STOP
END SUBROUTINE FINDELM

SUBROUTINE READPCN(CHKNELM,PRINUM,CONNUM,INTISG,BASNAM,LUBAS)
  IMPLICIT NONE
  LOGICAL :: CHKNELM
  INTEGER :: PRINUM,CONNUM,INTISG,LUBAS
  CHARACTER(*) :: BASNAM
  LOGICAL :: BLANK
  CHARACTER(88) :: STRING
  CHARACTER :: CHDUMMY
  INTEGER :: IDUMMY,IOERR

  CHKNELM=.FALSE.
1 READ(LUBAS,'(a88)',IOSTAT=IOERR,ERR=2,END=3)STRING
  READ(STRING,'(a1)')CHDUMMY
  IF (CHDUMMY==' ') THEN
     CALL CHBLANK(BLANK,STRING)
     IF (BLANK) THEN
        GO TO 1
     ELSE
        READ (STRING,'(BN,3i5)')PRINUM,CONNUM,INTISG
        RETURN
     END IF
  ELSE IF ((CHDUMMY=='a').OR.(CHDUMMY=='A')) THEN
     CHKNELM=.TRUE.
     RETURN
  ELSE
     GO TO 1
  END IF
2 WRITE(*,'(a,a,a,i2,a)')' Subroutine READPCN: an error ocurred during the reading of the file containing the basis ', &
 & BASNAM,' (IOSTAT code=',IOERR,') '
  STOP
3 CHKNELM=.TRUE.
  RETURN
END SUBROUTINE READPCN

SUBROUTINE READECC(MNOP,MNOC,NOP,NOC,ALPHA,CCOEF,BASNAM,LUBAS)
  IMPLICIT NONE
  INTEGER :: MNOP,MNOC,NOP,NOC,LUBAS
  DOUBLE PRECISION,DIMENSION(MNOP) :: ALPHA
  DOUBLE PRECISION,DIMENSION(MNOP,MNOC) :: CCOEF
  CHARACTER(*) :: BASNAM
  INTEGER :: I,J,K,L,IOERR,NUMLIN,NUMRCC,NUM
  CHARACTER :: CHDUMMY
  CHARACTER(88) :: STRING
  LOGICAL :: BLANK
! NUMCCA is the maximum number of contraction coefficients on the first line of a block. If NOC is greater than NUMCCA, the reading of the coefficients will continue on the following line, which contains NUMCCB coefficients at most.
  INTEGER,PARAMETER :: NUMCCA=6,NUMCCB=7

  ALPHA=0.D0 ; CCOEF=0.D0
  J=0
1 IF (J<NOP) THEN
     READ(LUBAS,'(a88)',IOSTAT=IOERR,ERR=3)STRING
     READ(STRING,'(a1)')CHDUMMY
     IF (CHDUMMY==' ') THEN
        CALL CHBLANK(BLANK,STRING)
        IF (.NOT.BLANK) THEN
! This line contains an exponent and some contraction coefficients
           J=J+1
           L=ABS(NOC)
           CALL CONLIN(L,NUMLIN)
           IF (NOC<=0) THEN
! We want an uncontracted basis set: we just read the exponent and skip the continuation lines
              READ(STRING,'(f16.9)')ALPHA(J)
              CCOEF(J,J)=1.D0
              DO I=2,NUMLIN
                 READ(LUBAS,*)
              END DO
              GO TO 1
           ELSE
! Get the right format for the read statement
              IF (NUMLIN==1) THEN
                 L=NOC
              ELSE
                 L=NUMCCA
              END IF
! Read the exponent and contraction coefficients on the first line
              READ(STRING,'(f16.9,6f12.9)')ALPHA(J),(CCOEF(J,I),I=1,L)
! If there are more lines with contraction coefficients, read them
              DO I=2,NUMLIN
                 NUMRCC=NUMCCA+(I-1)*NUMCCB
! Get the right format for the read statement
                 NUM=MIN(NUMRCC,NOC)
2                READ(LUBAS,'(a88)',IOSTAT=IOERR,ERR=3)STRING
                 READ (STRING,'(a1)')CHDUMMY
                 IF (CHDUMMY==' ') THEN
                    CALL CHBLANK(BLANK,STRING)
                    IF (.NOT.BLANK) THEN
! This line contains contraction coefficients
                       READ(STRING,'(7f12.9)')(CCOEF(J,K),K=NUMCCA+(I-2)*NUMCCB+1,NUM)
                    END IF
                 ELSE
! This line is a blank one, read the next one
                    GO TO 2
                 END IF
              END DO
              GO TO 1
           END IF
        ELSE
! This line is a blank one, read the next one
           GO TO 1
        END IF  
     ELSE
! Skip this line
       GO TO 1
     END IF
  END IF
  RETURN
3 WRITE(*,'(a,a,a,i2,a)')' Subroutine READECC: an error ocurred during the reading of the file containing the basis ', &
 & BASNAM,' (IOSTAT code=',IOERR,') '
  STOP
END SUBROUTINE READECC

SUBROUTINE CHBLANK(BLANK,STRING)
! Subroutine that checks whether a line is a blank one or not.
  IMPLICIT NONE
  LOGICAL :: BLANK
  CHARACTER*(*) :: STRING
  INTEGER :: J

  BLANK=.TRUE.
  DO J=1,LEN(STRING)
     BLANK=BLANK.AND.(STRING(J:J)==' ')
  END DO
  RETURN
END SUBROUTINE CHBLANK

SUBROUTINE CONLIN(NOC,NUMLIN)
! Subroutine that finds out on how many lines the contraction coefficients are written in a given block (knowing there at most NUMCCA coefficients on the first line and NUMCCB on the following ones), and returns that number.
  IMPLICIT NONE
  INTEGER :: NOC,NUMLIN
  INTEGER,PARAMETER :: NUMCCA=6,NUMCCB=7

  IF (NOC<=NUMCCA) THEN
     NUMLIN=1
  ELSE IF (MOD(NOC-NUMCCA,NUMCCB)==0) THEN
     NUMLIN=(NOC-NUMCCA)/NUMCCB+1
  ELSE
     NUMLIN=(NOC-NUMCCA)/NUMCCB+2
  END IF
END SUBROUTINE CONLIN

FUNCTION MDUSC(I,J) RESULT (MONOMIALDEGREE)
! Function that returns the monomial degree of each cartesian gaussian primitive function of shell type "I" (where "I=1" means "s", "I=2" means "p", "I=3" means "d", etc...) appearing in the components of the upper 2-spinor basis functions.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: I,J
  INTEGER,DIMENSION(3) :: MONOMIALDEGREE

  SELECT CASE (I)
     CASE (1)
! shell type s (1 associated monomial)
     MONOMIALDEGREE=(/0,0,0/)
     CASE (2)
! shell type p (3 associated monomials)
     SELECT CASE (J)
        CASE (1) ; MONOMIALDEGREE=(/1,0,0/)
        CASE (2) ; MONOMIALDEGREE=(/0,1,0/)
        CASE (3) ; MONOMIALDEGREE=(/0,0,1/)
     END SELECT
     CASE (3)
! shell type d (6 associated monomials)
     SELECT CASE (J)
        CASE (1) ; MONOMIALDEGREE=(/2,0,0/)
        CASE (2) ; MONOMIALDEGREE=(/1,1,0/)
        CASE (3) ; MONOMIALDEGREE=(/1,0,1/)
        CASE (4) ; MONOMIALDEGREE=(/0,2,0/)
        CASE (5) ; MONOMIALDEGREE=(/0,1,1/)
        CASE (6) ; MONOMIALDEGREE=(/0,0,2/)
     END SELECT
     CASE (4)
! shell type f (10 associated monomials)
     SELECT CASE (J)
        CASE (1)  ; MONOMIALDEGREE=(/3,0,0/)
        CASE (2)  ; MONOMIALDEGREE=(/2,1,0/)
        CASE (3)  ; MONOMIALDEGREE=(/2,0,1/)
        CASE (4)  ; MONOMIALDEGREE=(/1,2,0/)
        CASE (5)  ; MONOMIALDEGREE=(/1,0,2/)
        CASE (6)  ; MONOMIALDEGREE=(/1,1,1/)
        CASE (7)  ; MONOMIALDEGREE=(/0,3,0/)
        CASE (8)  ; MONOMIALDEGREE=(/0,2,1/)
        CASE (9)  ; MONOMIALDEGREE=(/0,1,2/)
        CASE (10) ; MONOMIALDEGREE=(/0,0,3/)
     END SELECT
  END SELECT
END FUNCTION MDUSC

FUNCTION MDLSC(I,J) RESULT (MONOMIALDEGREE)
! Function that returns the monomial degree of each cartesian gaussian primitive function appearing in the components of the lower 2-spinor basis functions when a kinetic balance scheme is applied on the upper 2-spinor basis functions.
! Note: the integer I refers to the shell type of the gaussian basis functions appearing in the components of the upper 2-spinor basis function to which the (R/U)KB scheme is applied (see the function MDUSC above).
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: I,J
  INTEGER,DIMENSION(3) :: MONOMIALDEGREE

  SELECT CASE (I)
     CASE (1)
! shell type s (3 associated monomials)
     SELECT CASE (J)
        CASE (1) ; MONOMIALDEGREE=(/1,0,0/)
        CASE (2) ; MONOMIALDEGREE=(/0,1,0/)
        CASE (3) ; MONOMIALDEGREE=(/0,0,1/)
     END SELECT
     CASE (2)
! shell type p (7 associated monomials)
     SELECT CASE (J)
        CASE (1) ; MONOMIALDEGREE=(/2,0,0/)
        CASE (2) ; MONOMIALDEGREE=(/1,1,0/)
        CASE (3) ; MONOMIALDEGREE=(/1,0,1/)
        CASE (4) ; MONOMIALDEGREE=(/0,2,0/)
        CASE (5) ; MONOMIALDEGREE=(/0,1,1/)
        CASE (6) ; MONOMIALDEGREE=(/0,0,2/)
        CASE (7) ; MONOMIALDEGREE=(/0,0,0/)
     END SELECT
     CASE (3)
! shell type d (13 associated monomials)
     SELECT CASE (J)
        CASE (1)  ; MONOMIALDEGREE=(/3,0,0/)
        CASE (2)  ; MONOMIALDEGREE=(/2,1,0/)
        CASE (3)  ; MONOMIALDEGREE=(/2,0,1/)
        CASE (4)  ; MONOMIALDEGREE=(/1,2,0/)
        CASE (5)  ; MONOMIALDEGREE=(/1,0,2/)
        CASE (6)  ; MONOMIALDEGREE=(/1,1,1/)
        CASE (7)  ; MONOMIALDEGREE=(/0,3,0/)
        CASE (8)  ; MONOMIALDEGREE=(/0,2,1/)
        CASE (9)  ; MONOMIALDEGREE=(/0,1,2/)
        CASE (10) ; MONOMIALDEGREE=(/0,0,3/)
        CASE (11) ; MONOMIALDEGREE=(/1,0,0/)
        CASE (12) ; MONOMIALDEGREE=(/0,1,0/)
        CASE (13) ; MONOMIALDEGREE=(/0,0,1/)
     END SELECT
     CASE (4)
! shell type f (21 associated monomials)
     SELECT CASE (J)
        CASE (1)  ; MONOMIALDEGREE=(/4,0,0/)
        CASE (2)  ; MONOMIALDEGREE=(/3,1,0/)
        CASE (3)  ; MONOMIALDEGREE=(/3,0,1/)
        CASE (4)  ; MONOMIALDEGREE=(/2,2,0/)
        CASE (5)  ; MONOMIALDEGREE=(/2,1,1/)
        CASE (6)  ; MONOMIALDEGREE=(/2,0,2/)
        CASE (7)  ; MONOMIALDEGREE=(/1,3,0/)
        CASE (8)  ; MONOMIALDEGREE=(/1,2,1/)
        CASE (9)  ; MONOMIALDEGREE=(/1,1,2/)
        CASE (10) ; MONOMIALDEGREE=(/1,0,3/)
        CASE (11) ; MONOMIALDEGREE=(/0,4,0/)
        CASE (12) ; MONOMIALDEGREE=(/0,3,1/)
        CASE (13) ; MONOMIALDEGREE=(/0,2,2/)
        CASE (14) ; MONOMIALDEGREE=(/0,1,3/)
        CASE (15) ; MONOMIALDEGREE=(/0,0,4/)
        CASE (16) ; MONOMIALDEGREE=(/2,0,0/)
        CASE (17) ; MONOMIALDEGREE=(/1,1,0/)
        CASE (18) ; MONOMIALDEGREE=(/1,0,1/)
        CASE (19) ; MONOMIALDEGREE=(/0,2,0/)
        CASE (20) ; MONOMIALDEGREE=(/0,1,1/)
        CASE (21) ; MONOMIALDEGREE=(/0,0,2/)
     END SELECT
  END SELECT
END FUNCTION MDLSC
END MODULE
