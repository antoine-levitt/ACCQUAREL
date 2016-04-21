SUBROUTINE DIIS_RHF(EIG,EIGVEC,NBAST,POEFM,PHI,TRSHLD,MAXITR,RESUME)
! DIIS (Direct Inversion in the Iterative Subspace) algorithm (restricted closed-shell Hartree-Fock formalism)
! Reference: P. Pulay, Convergence acceleration of iterative sequences. The case of SCF iteration, Chem. Phys. Lett., 73(2), 393-398, 1980.
! DIIS with hard restart when linear dependencies encountered
  USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE common_functions
  USE matrices ; USE matrix_tools ; USE metric_nonrelativistic ; USE scf_tools ; USE setup_tools
  INTEGER,INTENT(IN) :: NBAST
  DOUBLE PRECISION,DIMENSION(NBAST),INTENT(INOUT) :: EIG
  DOUBLE PRECISION,DIMENSION(NBAST,NBAST),INTENT(INOUT) :: EIGVEC
  DOUBLE PRECISION,DIMENSION(NBAST*(NBAST+1)/2),INTENT(IN) :: POEFM
  TYPE(gaussianbasisfunction),DIMENSION(NBAST),INTENT(IN) :: PHI
  DOUBLE PRECISION,INTENT(IN) :: TRSHLD
  INTEGER,INTENT(IN) :: MAXITR
  LOGICAL,INTENT(IN) :: RESUME

  INTEGER :: ITER,INFO,I,J,IJ,K
  INTEGER :: MXSET,MSET,MM
  DOUBLE PRECISION :: MXDIV,Y1,Y2
  INTEGER,DIMENSION(:),ALLOCATABLE :: IPIV
  DOUBLE PRECISION :: ETOT,ETOT1
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: PTEFM,PFM,PDM,PTDM,PBM,DIISV
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: ERRSET,PDMSET,TMP, YRTGSET !YRTGSET contains y_n-i orthogonal to the previous y_n-j 
  LOGICAL :: NUMCONV
  DOUBLE PRECISION :: TSTART,TEND
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: DMEND,DMI 
  CHARACTER(30) :: file_name,lastfile_name

! INITIALIZATIONS AND PRELIMINARIES
! Reading of the maximum dimension of the density matrix simplex
  OPEN(100,FILE=SETUP_FILE,STATUS='OLD',ACTION='READ')
  CALL LOOKFOR(100,'DIIS ALGORITHM PARAMETERS',INFO)
  READ(100,'(/,i2)')MXSET
  READ(100,'(/,f16.8)')MXDIV !parameter to check the independence of vectors of Gram matrix
  CLOSE(100)
  WRITE(*,*)'Lindepcheck =',MXDIV

  ALLOCATE(PDM(1:NBAST*(NBAST+1)/2),PTDM(1:NBAST*(NBAST+1)/2))
  ALLOCATE(PTEFM(1:NBAST*(NBAST+1)/2),PFM(1:NBAST*(NBAST+1)/2))
  ALLOCATE(ERRSET(1:MXSET,1:NBAST*NBAST),PDMSET(1:MXSET,1:NBAST*(NBAST+1)/2),TMP(1:NBAST,1:NBAST),YRTGSET(1:MXSET,1:NBAST*NBAST))

  ITER=0
  MSET=0 ; PDMSET=0.D0
  PTDM=0.D0
  ETOT1=0.D0
  OPEN(16,FILE='plots/diisenrgy.txt',STATUS='unknown',ACTION='write')
  OPEN(17,FILE='plots/diiscrit1.txt',STATUS='unknown',ACTION='write')
  OPEN(18,FILE='plots/diiscrit2.txt',STATUS='unknown',ACTION='write')
  OPEN(19,FILE='plots/diisrstrt.txt',STATUS='unknown',ACTION='write')
  CALL CPU_TIME(TSTART)

! LOOP
1 CONTINUE
  ITER=ITER+1
  WRITE(*,*)' '
  WRITE(*,*)'# ITER =',ITER
  WRITE(file_name,fmt='(A14,I0,A4)')'plots/DM_DIIS_',ITER,'.txt'
  OPEN(20,FILE=file_name,STATUS='unknown')

! Assembly and diagonalization of the Fock matrix
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PTDM)
  PFM=POEFM+PTEFM
  CALL EIGENSOLVER(PFM,PCFS,NBAST,EIG,EIGVEC,INFO)
  IF (INFO/=0) GO TO 5
! Assembly of the density matrix according to the aufbau principle
  CALL FORMDM(PDM,EIGVEC,NBAST,1,NBE/2)
! Computation of the total energy
  CALL BUILDTEFM(PTEFM,NBAST,PHI,PDM)
  ETOT=ENERGY(POEFM,PTEFM,PDM,NBAST)
  WRITE(*,*)'Total energy =',ETOT
! Numerical convergence check
  IF (ITER==1) THEN
     CALL CHECKNUMCONV(PDM,PDMSET(1,:),POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV)
  ELSE
     CALL CHECKNUMCONV(PDM,PDMSET(MSET,:),POEFM+PTEFM,NBAST,ETOT,ETOT1,TRSHLD,NUMCONV)
  END IF
  IF (NUMCONV) THEN
! Convergence reached
     CALL PRINTMATRIX(PDM,NBAST,20)
     GO TO 2
  ELSE IF (ITER==MAXITR) THEN
! Maximum number of iterations reached without convergence
     GO TO 3
  ELSE
! Convergence not reached, increment
     ETOT1=ETOT
! Storage of the current density matrix, computation and storage of a new DIIS error vector
     IF (MSET==MXSET) THEN
! we have reached the maximum allowed dimension for the simplex
        PDMSET=0.D0 ! elimination of all the stored density matrix except the last computed
        ERRSET=0.D0 ! elimination of all the stored error vector except the last computed
	YRTGSET=0.D0
	MSET=1
	WRITE(19,'(a,i3,a)') 'Restart at step',ITER,' max dimension of simplex reached'
     ELSE
        MSET=MSET+1 ! the current dimension of the simplex is increased by one
     END IF
     PDMSET(MSET,:)=PDM
     TMP=COMMUTATOR(POEFM+PTEFM,PDMSET(MSET,:),PS,NBAST)
     IJ=0
     DO I=1,NBAST
        DO J=1,NBAST
           IJ=IJ+1
           ERRSET(MSET,IJ)=TMP(I,J)
        END DO
     END DO
! check the linear independence of the Gram matrix
     IF (MSET==1 .OR. ITER==1) THEN
        PTDM=PDM
     ELSE IF (MSET==2) THEN
	YRTGSET(1,:)=ERRSET(2,:)-ERRSET(1,:)
WRITE(*,*)'Dimension of the density matrix simplex =',MSET
! Computation of the new pseudo-density matrix
! assembly and solving of the linear system associated to the DIIS equations
        MM=(MSET+1)*(MSET+2)/2
        ALLOCATE(PBM(1:MM))
        PBM=0.D0 ; PBM(MM-MSET:MM-1)=-1.D0
        IJ=0
        DO J=1,MSET
           DO I=1,J
              IJ=IJ+1 ; PBM(IJ)=PBM(IJ)+DOT_PRODUCT(ERRSET(J,:),ERRSET(I,:))
           END DO
        END DO
        ALLOCATE(DIISV(1:MSET+1))
        DIISV(1:MSET)=0.D0 ; DIISV(MSET+1)=-1.D0
        ALLOCATE(IPIV(1:MSET+1))
        CALL DSPSV('U',MSET+1,1,PBM,IPIV,DIISV,MSET+1,INFO)
        DEALLOCATE(IPIV,PBM)
        IF (INFO/=0) GO TO 4
! assembly of the pseudo-density matrix
        PTDM=0.D0
        DO I=1,MSET
           PTDM=PTDM+DIISV(I)*PDMSET(I,:)
        END DO
        DEALLOCATE(DIISV)
     ELSE
!general formula
	I=2
	YRTGSET(MSET-1,:)=ERRSET(MSET,:)-ERRSET(MSET-1,:)
	DO I=1,MSET-2
        YRTGSET(MSET-1,:)=YRTGSET(MSET-1,:)-DOT_PRODUCT(YRTGSET(I,:),ERRSET(MSET,:)-ERRSET(MSET-1,:))&
&/(NORM2(YRTGSET(I,:))**2)*YRTGSET(I,:)
	END DO
!	     IF (YCHK .AND. (MXDIV**2)*DOT_PRODUCT(ERRSET(MSET,:)-ERRSET(MSET-1,:),ERRSET(MSET,:)-ERRSET(MSET-1,:)) > &
!&DOT_PRODUCT(YRTGSET(MSET-1,:)-ERRSET(MSET,:)+ERRSET(MSET-1,:),YRTGSET(MSET-1,:)-ERRSET(MSET,:)+ERRSET(MSET-1,:))) THEN
		! check the linear independence of the Gram matrix
WRITE(*,*)'Y NORM =',NORM2(ERRSET(MSET,:)-ERRSET(MSET-1,:))
WRITE(*,*)'Y-QY NORM =',NORM2(YRTGSET(MSET-1,:))
	     IF (MXDIV*NORM2(ERRSET(MSET,:)-ERRSET(MSET-1,:)) > NORM2(YRTGSET(MSET-1,:))) THEN
		WRITE(19,'(a,i2,a)')'Restart at step ',ITER,' linear dependencies encountered'
		YRTGSET=0.D0
		PTDM=PDM
		PDMSET=EOSHIFT(PDMSET,MSET-1) ! elimination of all the stored density matrix except the last computed
        	ERRSET=EOSHIFT(ERRSET,MSET-1) ! elimination of all the stored error vector except the last computed
		MSET=1
	     ELSE
        	WRITE(*,*)'Dimension of the density matrix simplex =',MSET
		! Computation of the new pseudo-density matrix
		! assembly and solving of the linear system associated to the DIIS equations
        	MM=(MSET+1)*(MSET+2)/2
        	ALLOCATE(PBM(1:MM))
        	PBM=0.D0 ; PBM(MM-MSET:MM-1)=-1.D0
        	IJ=0
        	DO J=1,MSET
        	   DO I=1,J
              		IJ=IJ+1 ; PBM(IJ)=PBM(IJ)+DOT_PRODUCT(ERRSET(J,:),ERRSET(I,:))
        	   END DO
        	END DO
        	ALLOCATE(DIISV(1:MSET+1))
        	DIISV(1:MSET)=0.D0 ; DIISV(MSET+1)=-1.D0
        	ALLOCATE(IPIV(1:MSET+1))
        	CALL DSPSV('U',MSET+1,1,PBM,IPIV,DIISV,MSET+1,INFO)
        	DEALLOCATE(IPIV,PBM)
        	IF (INFO/=0) GO TO 4
		! assembly of the pseudo-density matrix
        	PTDM=0.D0
        	DO I=1,MSET
          	PTDM=PTDM+DIISV(I)*PDMSET(I,:)
        	END DO
        	DEALLOCATE(DIISV)
	    END IF
	END IF
	CALL PRINTMATRIX(PDM,NBAST,20)
     	GO TO 1
  END IF
! MESSAGES
2 CALL CPU_TIME(TEND)
  WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: convergence after',ITER,'iteration(s) and',TEND-TSTART,'s.'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
! Computation of ||D_i - D_*||
  WRITE(lastfile_name,fmt='(A14,I0,A4)')'plots/DM_DIIS_',ITER,'.txt'
  OPEN(21,FILE=lastfile_name,STATUS='old',ACTION='READ')
  OPEN(22,FILE='plots/DM_DIIS_norm.txt',STATUS='unknown',ACTION='WRITE')
  ALLOCATE(DMEND(1:NBAST*(NBAST+1)/2))
  CALL READMATRIX(21,NBAST,DMEND)
  DO I=1,ITER
     WRITE(file_name,fmt='(A14,I0,A4)')'plots/DM_DIIS_',I,'.txt'
     OPEN(23,FILE=file_name,STATUS='old',ACTION='READ')
     ALLOCATE(DMI(1:NBAST*(NBAST+1)/2))
     CALL READMATRIX(23,NBAST,DMI)
     WRITE(22,'(f20.12)')NORM(DMI-DMEND,NBAST,'F')
     DEALLOCATE(DMI)
  END DO
  DEALLOCATE(DMEND)
  CLOSE(21) ; CLOSE(22) ;
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
3 WRITE(*,*)' ' ; WRITE(*,*)'Subroutine DIIS: no convergence after',MAXITR,'iteration(s).'
  OPEN(9,FILE='eigenvalues.txt',STATUS='UNKNOWN',ACTION='WRITE')
  DO I=1,NBAST
     WRITE(9,*)I,EIG(I)
  END DO
  CLOSE(9)
  GO TO 6
4 IF (INFO<0) THEN
     WRITE(*,*)'Subroutine ZHPSV: the',-INFO,'-th argument had an illegal value'
  ELSE
     WRITE(*,*)'Subroutine ZHPSV: the factorization has been completed, but the block diagonal matrix is &
    &exactly singular, so the solution could not be computed'
  END IF
  GO TO 5
5 WRITE(*,*)'(called from subroutine DIIS)'
6 DEALLOCATE(PDM,PTDM,PDMSET,PTEFM,PFM,TMP,ERRSET,YRTGSET)
  CLOSE(16) ; CLOSE(17) ; CLOSE(18) ; CLOSE(19) ; 
END SUBROUTINE DIIS_RHF
