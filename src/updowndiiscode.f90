     ELSE IF (MSET==2) THEN
! Computation of the new pseudo-density matrix by solving the unsconstrained form
! of the least-squares problem defining the DIIS using the QR factorisation with
! updating (and downdating)
! first QR factorization
        ALLOCATE(TAU(1:1),WORK(1:1))
        QRM(:,1)=DIFERRSET(:,1)
        CALL DGEQRF(NBAST*NBAST,1,QRM(:,1),NBAST*NBAST,TAU,WORK,1,INFO)
        DEALLOCATE(WORK)
        IF (INFO/=0) GO TO 4
! form Q
        QM(:,1)=QRM(:,1)
        ALLOCATE(WORK(1:NBAST*NBAST))
        CALL DORGQR(NBAST*NBAST,NBAST*NBAST,1,QM,NBAST*NBAST,TAU,WORK,NBAST*NBAST,INFO)
        DEALLOCATE(TAU,WORK)
        IF (INFO/=0) GO TO 5
     ELSE IF ((MSET>2).AND.((MSET<MXSET).OR.(MXSET==0))) THEN
! update of the QR factorization
        QRM(:,1:MSET-1)=MATMUL(TRANSPOSE(QM),DIFERRSET(:,1:MSET-1))
! update R
        ALLOCATE(TAU(1:1),WORK(1:1))
        CALL ADDCOLS(NBAST*NBAST,MSET-1,QRM(:,1:MSET-1),NBAST*NBAST,MSET-1,1,TAU,WORK,1,INFO)
        DEALLOCATE(WORK)
        IF (INFO/=0) GO TO 6
!  write(33,*)QRM(:,1:MSET-1)
! update Q
        ALLOCATE(WORK(1:2*MSET-2))
        CALL ADDCOLSQ(NBAST*NBAST,MSET-1,QRM(:,1:MSET-1),NBAST*NBAST,QM,NBAST*NBAST,MSET-1,1,TAU,WORK,INFO)
  write(*,*)'toto'
        DEALLOCATE(TAU,WORK)
        IF (INFO/=0) GO TO 7
     ELSE IF ((MSET>2).AND.(MSET==MXSET)) THEN
! downdate of the QR factorization
        QRM(:,1:MSET-2)=MATMUL(TRANSPOSE(QM),DIFERRSET(:,1:MSET-2))
! downdate R
        ALLOCATE(TAU(1:MSET-2),WORK(1:2))
        CALL DELCOLS(NBAST*NBAST,MSET-2,QRM(:,1:MSET-2),NBAST*NBAST,1,1,TAU,WORK,INFO)
        IF (INFO/=0) GO TO 8
! downdate Q
        CALL DELCOLSQ(NBAST*NBAST,MSET-2,QRM(:,1:MSET-2),NBAST*NBAST,QM,NBAST*NBAST,1,1,TAU,WORK,INFO)
        DEALLOCATE(TAU,WORK)
        IF (INFO/=0) GO TO 9
! update of the QR factorization
        QRM(:,1:MSET-1)=MATMUL(TRANSPOSE(QM),DIFERRSET(:,1:MSET-1))
! update R
        ALLOCATE(TAU(1:1),WORK(1:1))
        CALL ADDCOLS(NBAST*NBAST,MSET-1,QRM(:,1:MSET-1),NBAST*NBAST,MSET-1,1,TAU,WORK,1,INFO)
        DEALLOCATE(WORK)
        IF (INFO/=0) GO TO 6
! update Q
        ALLOCATE(WORK(1:2*MSET-2))
        CALL ADDCOLSQ(NBAST*NBAST,MSET-1,QRM(:,1:MSET-1),NBAST*NBAST,QM,NBAST*NBAST,MSET-1,1,TAU,WORK,INFO)
        DEALLOCATE(TAU,WORK)
        IF (INFO/=0) GO TO 7
     END IF
! solve the linear system
     RHS=MATMUL(TRANSPOSE(QM),ERRSET(:,MSET))
     CALL DTRTRS('U','N','N',MSET-1,1,QRM,MSET-1,RHS,MSET-1,INFO)

! routines for QR decomposition updating and downdating
! Source: Updating the QR factorization and the least squares problem
! Sven Hammarling and Craig Lucas
! MIMS EPrint: 2008.111
! November 2008

SUBROUTINE DELCOLS(M,N,A,LDA,K,P,TAU,WORK,INFO)
IMPLICIT NONE
!
! Craig Lucas,University of Manchester
! March,2004
!
! .. Scalar Arguments ..
INTEGER INFO,K,LDA,M,N,P
! ..
! .. Array Arguments ..
DOUBLE PRECISION A(LDA,*),TAU(*),WORK(*)
! ..
!
! Purpose
! =======
!
! Given a real m by (n+p) matrix, B, and the QR factorization
! B=Q_B * R_B, DELCOLS computes the QR factorization
! C=Q * R where C is the matrix B with p columns deleted
! from the kth column onwards.
!
! The input to this routine is Q_B' * C
!
! Arguments
! =========
!
! M (input) INTEGER
! The number of rows of the matrix C. M >= 0.
!
! N (input) INTEGER
! The number of columns of the matrix C. N >= 0.
!
! A (input/output) DOUBLE PRECISION array, dimension (LDA,N)
! On entry,the matrix Q_B' * C. The elements in columns
! 1:K-1 are not referenced.
! On exit, the elements on and above the diagonal contain
! the n by n upper triangular part of the matrix R. The
! elements below the diagonal in columns k:n,together with
! TAU represent the orthogonal matrix Q as a product of
! elementary reflectors (see Further Details).
!
! LDA (input) INTEGER
! The leading dimension of the array A. LDA >= max(1,M).
!
! K (input) INTEGER
! The position of the first column deleted from B.
! 0 < K <= N+P.
!
! P (input) INTEGER
! The number of columns deleted from B. P > 0.
!
! TAU (output) DOUBLE PRECISION array, dimension(N-K+1)
! The scalar factors of the elementary reflectors
! (see Further Details).
!
! WORK (workspace) DOUBLE PRECISION array, dimension (P+1)
! Work space.
!
! INFO (output) INTEGER
! =0: successful exit
! < 0: if INFO=-I, the I-th argument had an illegal value.
!
! Further Details
! ===============
!
! The matrix Q is represented as a product of Q_B and elementary
! reflectors
!
! Q=Q_B * H(k) * H(k+1) *...* H(last),last=min(m-1,n).
!
! Each H(j) has the form
!
! H(j)=I - tau*v*v'
!
! where tau is a real scalar, and v is a real vector with
! v(1:j-1)=0,v(j)=1, v(j+1:j+lenh-1), lenh=min(p+1,m-j+1),
! stored on exit in A(j+1:j+lenh-1,j) and v(j+lenh:m)=0, tau is
! stored in TAU(j).
!
! The matrix Q can be formed with DELCOLSQ
!
! =====================================================================
!
! .. Parameters ..
DOUBLE PRECISION ONE
PARAMETER (ONE=1.0D+0)
! ..
! .. Local Scalars ..
DOUBLE PRECISION AJJ
INTEGER J,LAST,LENH
! ..
! .. External Subroutines ..
EXTERNAL DLARF,DLARFG,XERBLA
! ..
!
! Test the input parameters
!
INFO=0
IF (M.LT.0) THEN
   INFO=-1
ELSE IF (N.LT.0) THEN
   INFO=-2
ELSE IF (LDA.LT.MAX(1,M)) THEN
   INFO=-4
ELSE IF (K.GT.N+P.OR.K.LE.0) THEN
   INFO=-5
ELSE IF (P.LE.0) THEN
   INFO=-6
END IF
IF (INFO.NE.0) THEN
   CALL XERBLA('DELCOLS',-INFO)
   RETURN
END IF
LAST=MIN(M-1,N)
DO J=K,LAST
!
! Generate elementary reflector H(J) to annihilate the nonzero entries below A(J,J)
!
   LENH=MIN(P+1,M-J+1)
   CALL DLARFG(LENH,A(J,J),A(J+1,J),1,TAU(J-K+1))
   IF (J.LT.N) THEN
!
! Apply H(J) to trailing matrix from left
!
      AJJ=A(J,J)
      A(J,J)=ONE
      CALL DLARF('L',LENH,N-J,A(J,J),1,TAU(J-K+1),A(J,J+1),LDA,WORK)
      A(J,J)=AJJ
   END IF
END DO
END SUBROUTINE DELCOLS

SUBROUTINE DELCOLSQ(M,N,A,LDA,Q,LDQ,K,P,TAU,WORK,INFO)
IMPLICIT NONE
!
! Craig Lucas, University of Manchester
! March, 2004
!
! .. Scalar Arguments ..
INTEGER INFO,K,LDA,LDQ,M,N,P
! ..
! .. Array Arguments ..
DOUBLE PRECISION A(LDA,*),Q(LDQ,*),TAU(*),WORK(*)
! ..
!
! Purpose
! =======
!
! DELCOLSQ generates an m by m real matrix Q with orthogonal columns,
! which is defined as the product of Q_B and elementary reflectors
!
! Q = Q_B * H(k) * H(k+1) *...* H(last), last = min( m-1, n ) .
!
! where the H(j) are as returned by DELCOLSQ, such that C = Q * R and
! C is the matrix B = Q_B * R_B, with p columns deleted from the
! kth column onwards.
!
! Arguments
! =========
!
! M (input) INTEGER
! The number of rows of the matrix A. M >= 0.
!
! N (input) INTEGER
! The number of columns of the matrix A. N >= 0.
!
! A (input) DOUBLE PRECISION array, dimension (LDA,N)
! On entry, the elements below the diagonal in columns k:n
! must contain the vector which defines the elementary
! reflector H(J) as returned by DELCOLS.
!
! LDA (input) INTEGER
! The leading dimension of the array A. LDA >= max(1,M).
!
! Q (input/output) DOUBLE PRECISION array, dimension (LDA,N)
! On entry, the matrix Q_B.
! On exit, the matrix Q.
!
! LDQ (input) INTEGER
! The leading dimension of the array Q. LDQ >= M.
!
! K (input) INTEGER
! The position of the first column deleted from B.
! 0 < K <= N+P.
!
! P (input) INTEGER
! The number of columns deleted from B. P > 0.
!
! TAU (input) DOUBLE PRECISION array, dimension(N-K+1)
! TAU(J) must contain the scalar factor of the elementary
! reflector H(J), as returned by DELCOLS.
!
! WORK DOUBLE PRECISION array, dimension (P+1)
! Work space.
!
! INFO (output) INTEGER
! = 0: successful exit
! < 0: if INFO = -I, the I-th argument had an illegal value.
!
! =====================================================================
!
! .. Parameters ..
DOUBLE PRECISION ONE
PARAMETER (ONE=1.0D+0)
! ..
! .. Local Scalars ..
DOUBLE PRECISION AJJ
INTEGER J,LAST,LENH
! ..
! .. External Subroutines ..
EXTERNAL DLARF,XERBLA
! ..
!
! Test the input parameters
!
INFO=0
IF (M.LT.0) THEN
   INFO=-1
ELSE IF (N.LT.0) THEN
   INFO=-2
ELSE IF (LDA.LT.MAX(1,M)) THEN
   INFO=-4
ELSE IF (K.GT.N+P.OR.K.LE.0) THEN
   INFO=-5
ELSE IF (P.LE.0) THEN
   INFO=-6
END IF
IF (INFO.NE.0) THEN
   CALL XERBLA('DELCOLSQ',-INFO)
   RETURN
END IF
LAST=MIN(M-1,N)
DO J=K,LAST
   LENH=MIN(P+1,M-J+1)
!
! Apply H(J) from right
!
   AJJ=A(J,J)
   A(J,J)=ONE
   CALL DLARF('R',M,LENH,A(J,J),1,TAU(J-K+1),Q(1,J),LDQ,WORK)
   A(J,J)=AJJ
END DO
END SUBROUTINE DELCOLSQ

SUBROUTINE ADDCOLS(M,N,A,LDA,K,P,TAU,WORK,LWORK,INFO)
IMPLICIT NONE
!
! Craig Lucas, University of Manchester
! March, 2004
!
! .. Scalar Arguments ..
INTEGER INFO,K,LDA,LWORK,M,N,P
! ..
! .. Array Arguments ..
DOUBLE PRECISION A(LDA,*),TAU(*),WORK(*)
! ..
!
! Purpose
! =======
!
! Given a real m by (n-p) matrix, B, and the QR factorization
! B = Q_B * R_B, ADDCOLS computes the QR factorization
! C = Q * R where C is the matrix B with p columns added
! in the kth column onwards.
!
! The input to this routine is Q_B' * C
!
! Arguments
! =========
!
! M (input) INTEGER
! The number of rows of the matrix C. M >= 0.
!
! N (input) INTEGER
! The number of columns of the matrix C. N >= 0.
!
! A (input/output) DOUBLE PRECISION array, dimension (LDA,N)
! On entry, the matrix Q_B' * C. The elements in columns
! 1:K-1 are not referenced.
! On exit, the elements on and above the diagonal contain
! the n by n upper triangular part of the matrix R. The
! elements below the diagonal in columns K:N, together with
! TAU represent the orthogonal matrix Q as a product of
! elementary reflectors and Givens rotations.
! (see Further Details).
!
! LDA (input) INTEGER
! The leading dimension of the array A. LDA >= max(1,M).
!
! K (input) INTEGER
! The position of the first column added to B.
! 0 < K <= N-P+1.
!
! P (input) INTEGER
! The number of columns added to B. P > 0.
!
! TAU (output) DOUBLE PRECISION array, dimension(P)
! The scalar factors of the elementary reflectors
! (see Further Details).
!
! WORK (workspace) DOUBLE PRECISION array, dimension ( LWORK )
! Work space.
!
! LWORK (input) INTEGER
! The dimension of the array WORK. LWORK >= P.
! For optimal performance LWORK >= P*NB, where NB is the
! optimal block size.
!
! INFO (output) INTEGER
! = 0: successful exit
! < 0: if INFO = -I, the I-th argument had an illegal value.
!
! Further Details
! ===============
!
! The matrix Q is represented as a product of Q_B, elementary
! reflectors and Givens rotations
!
! Q = Q_B * H(k) * H(k+1) *...* H(k+p-1) * G(k+p-1,k+p) *...
! *G(k,k+1) * G(k+p,k+p+1) *...* G(k+2p-2,k+2p-1)
!
! Each H(j) has the form
!
! H(j) = I - tau*v*v'
!
! where tau is a real scalar, and v is a real vector with
! v(1:n-p-j+1) = 0, v(j) = 1, and v(j+1:m) stored on exit in
! A(j+1:m,j), tau is stored in TAU(j).
!
! Each G(i,j) has the form
!
! i-1 i
! [ I ]
! [ c -s ] i-1
! G(i,j) = [ s c ] i
! [ I ]
!
! and zero A(i,j), where c and s are encoded in scalar and
! stored in A(i,j) and
!
! IF A(i,j) = 1, c = 0, s = 1
! ELSE IF | A(i,j) | < 1, s = A(i,j), c = sqrt(1-s**2)
! ELSE c = 1 / A(i,j), s = sqrt(1-c**2)
!
! The matrix Q can be formed with ADDCOLSQ
!
! =====================================================================
!
! .. Local Scalars ..
DOUBLE PRECISION C,S
INTEGER I,INC,ISTART,J,JSTOP,UPLEN
! ..
! .. External Subroutines ..
EXTERNAL DGEQRF,DLASR,DROT,DROTG,XERBLA
! ..
!
! Test the input parameters
!
INFO=0
IF (M.LT.0) THEN
   INFO=-1
ELSE IF (N.LT.0) THEN
   INFO=-2
ELSE IF (LDA.LT.MAX(1,M)) THEN
   INFO=-4
ELSE IF (K.GT.N-P+1.OR.K.LE.0) THEN
   INFO=-5
ELSE IF (P.LE.0) THEN
   INFO=-6
END IF
IF (INFO.NE.0) THEN
   CALL XERBLA('ADDCOLS',-INFO)
   RETURN
END IF
!
! Do a QR factorization on rows below N-P, if there is more than one
!
IF (M.GT.N-P+1) THEN
!
! Level 3 QR factorization
!
   CALL DGEQRF(M-N+P,P,A(N-P+1,K),LDA,TAU,WORK,LWORK,INFO)
END IF
!
! If K not equal to number of columns in B and not <= M-1 then
! there is some elimination by Givens to do
!
IF (K+P-1.NE.N.AND.K.LE.M-1) THEN
!
! Zero out the rest with Givens
! Allow for M < N
!
   JSTOP=MIN(P+K-1,M-1)
   DO J=K,JSTOP
!
! Allow for M < N
!
      ISTART=MIN(N-P+J-K+1,M)
      UPLEN=N-K-P-ISTART+J+1
      INC=ISTART-J
      DO I=ISTART,J+1,-1
!
! Recall DROTG updates A( I-1, J ) and stores C and S encoded as scalar in A( I, J )
!
         CALL DROTG(A(I-1,J),A(I,J),C,S)
         WORK(INC)=C
         WORK(N+INC)=S
!
! Update nonzero rows of R
! Do the next two line this way round because A( I-1, N-UPLEN+1 ) gets updated
!
         A(I,N-UPLEN)=-S*A(I-1,N-UPLEN)
         A(I-1,N-UPLEN)=C*A(I-1,N-UPLEN)
         CALL DROT(UPLEN,A(I-1,N-UPLEN+1),LDA,A(I,N-UPLEN+1),LDA,C,S)
         UPLEN=UPLEN+1
         INC=INC-1
      END DO
!
! Update inserted columns in one go
! Max number of rotations is N-1, we've allowed N
!
      IF (J.LT.P+K-1) THEN
         CALL DLASR('L','V','B',ISTART-J+1,K+P-1-J,WORK(1),WORK(N+1),A(J,J+1),LDA)
      END IF
   END DO
END IF
END SUBROUTINE ADDCOLS

SUBROUTINE ADDCOLSQ(M,N,A,LDA,Q,LDQ,K,P,TAU,WORK,INFO)
IMPLICIT NONE
!
! Craig Lucas, University of Manchester
! March, 2004
!
! .. Scalar Arguments ..
INTEGER INFO,K,LDA,LDQ,M,N,P
! ..
! .. Array Arguments ..
DOUBLE PRECISION A(LDA,*),Q(LDQ,*),TAU(*),WORK(*)
! ..
!
! Purpose
! =======
!
! ADDCOLSQ generates an m by m real matrix Q with orthogonal columns,
! which is defined as the product of Q_B, elementary reflectors and
! Givens rotations
!
! Q = Q_B * H(k) * H(k+1) *...* H(k+p-1) * G(k+p-1,k+p) *...
! *G(k,k+1) * G(k+p,k+p+1) *...* G(k+2p-2,k+2p-1)
!
! where the H(j) and G(i,j) are as returned by ADDCOLS, such that
! C = Q * R and C is the matrix B = Q_B * R_B, with p columns added
! from the kth column onwards.
!
! Arguments
! =========
!
! M (input) INTEGER
! The number of rows of the matrix A. M >= 0.
!
! N (input) INTEGER
! The number of columns of the matrix A. N >= 0.
!
! A (input) DOUBLE PRECISION array, dimension (LDA,N)
! On entry, the elements below the diagonal in columns
! K:K+P-1 (if M > M-P+1) must contain the vector which defines
! the elementary reflector H(J). The elements above these
! vectors and below the diagonal store the scalars such that
! the Givens rotations can be constructed, as returned by
! ADDCOLS.
!
! LDA (input) INTEGER
! The leading dimension of the array A. LDA >= max(1,M).
!
! Q (input/output) DOUBLE PRECISION array, dimension (LDA,N)
! On entry, the matrix Q_B.
! On exit, the matrix Q.
!
! LDQ (input) INTEGER
! The leading dimension of the array Q. LDQ >= M.
!
! K (input) INTEGER
! The postion of first column added to B.
! 0 < K <= N-P+1.
!
! P (input) INTEGER
! The number columns added. P > 0.
!
! TAU (output) DOUBLE PRECISION array, dimension(N-K+1)
! The scalar factors of the elementary reflectors.
!
! WORK (workspace) DOUBLE PRECISION array, dimension (2*N)
! Work space.
!
! INFO (output) INTEGER
! = 0: successful exit
! < 0: if INFO = -I, the I-th argument had an illegal value
!
! =====================================================================
!
! .. Parameters ..
DOUBLE PRECISION ONE,ZERO
PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
! ..
! .. Local Scalars ..
DOUBLE PRECISION DTEMP
INTEGER COL,I,INC,ISTART,J,JSTOP
! ..
! .. External Subroutines ..
EXTERNAL DLARF,DLASR,XERBLA
! ..
!
! Test the input parameters
!
INFO=0
IF (M.LT.0) THEN
   INFO=-1
ELSE IF (N.LT.0) THEN
   INFO=-2
ELSE IF (LDA.LT.MAX(1,M)) THEN
   INFO=-4
ELSE IF (K.GT.N-P+1.OR.K.LE.0) THEN
   INFO=-5
ELSE IF(P.LE.0) THEN
   INFO=-6
END IF
IF (INFO.NE.0) THEN
   CALL XERBLA('ADDCLQ',-INFO)
   RETURN
END IF
!
! We did a QR factorization on rows below N-P+1
!
IF (M.GT.N-P+1) THEN
   COL=N-P+1
   DO J=K,K+P-1
      DTEMP=A(COL,J)
      A(COL,J)=ONE
!
! If N+P > M-N we have only factored the first M-N columns
!
      IF (M-COL+1.LE.0) GO TO 10
      CALL DLARF('R',M,M-COL+1,A(COL,J),1,TAU(J-K+1),Q(1,COL),LDQ,WORK)
      A(COL,J)=DTEMP
      COL=COL+1
10    CONTINUE
   END DO
END IF
!
! If K not equal to number of columns in B then there was some elimination by Givens
!
IF (K+P-1.LT.N.AND.K.LE.M-1) THEN
!
! Allow for M < N, i.e DO P wide unless hit the bottom first
!
   JSTOP=MIN(P+K-1,M-1)
   DO J=K,JSTOP
      ISTART=MIN(N-P+J-K+1,M)
      INC=ISTART-J
!
! Compute vectors of C and S for rotations
!
      DO I=ISTART,J+1,-1
         IF (A(I,J).EQ.ONE) THEN
            WORK(INC)=ZERO
            WORK(N+INC)=ONE
         ELSE IF (ABS(A(I,J)).LT.ONE) THEN
            WORK(N+INC)=A(I,J)
            WORK(INC)=SQRT((1-A(I,J)**2))
         ELSE
            WORK(INC)=ONE/A(I,J)
            WORK(N+INC)=SQRT((1-WORK(INC)**2))
         END IF
         INC=INC-1
      END DO
!
! Apply rotations to the Jth column from the right
!
      CALL DLASR('R','V','b',M,ISTART-I+1,WORK(1),WORK(N+1),Q(1,I),LDQ)
   END DO
END IF
END SUBROUTINE ADDCOLSQ
