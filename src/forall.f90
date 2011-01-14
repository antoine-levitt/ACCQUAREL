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
        ACTION(I,I,I,I)
! 2 distinct values for the 4 indices
     ELSE IF ((I>J).AND.(J==K).AND.(K==L)) THEN
        ACTION(I,J,J,J)
        ACTION(J,J,I,J)
        ACTION(J,J,J,I)
        ACTION(J,I,J,J)
     ELSE IF ((I==J).AND.(J==K).AND.(K>L)) THEN
        ACTION(L,I,I,I)
        ACTION(I,I,L,I)
        ACTION(I,I,I,L)
        ACTION(I,L,I,I)
     ELSE IF ((I==J).AND.(J>K).AND.(K==L)) THEN
        ACTION(I,I,K,K)
        ACTION(K,K,I,I)
     ELSE IF ((I==K).AND.(K>J).AND.(J==L)) THEN
        ACTION(I,J,I,J)
        ACTION(J,I,J,I)
        ACTION(J,I,I,J)
        ACTION(I,J,J,I)
! 3 distinct values for the 4 indices
     ELSE IF ((I==K).AND.(K>J).AND.(J>L)) THEN
        ACTION(I,J,I,L)
        ACTION(J,I,I,L)
        ACTION(I,J,L,I)
        ACTION(J,I,L,I)
        
        ACTION(I,L,I,J)
        ACTION(L,I,I,J)
        ACTION(I,L,J,I)
        ACTION(L,I,J,I)
     ELSE IF ((I>J).AND.(J==K).AND.(K>L)) THEN
        ACTION(I,J,J,L)
        ACTION(J,I,J,L)
        ACTION(I,J,L,J)
        ACTION(J,I,L,J)
        
        ACTION(J,L,I,J)
        ACTION(L,J,I,J)
        ACTION(J,L,J,I)
        ACTION(L,J,J,I)
     ELSE IF ((I>K).AND.(K>J).AND.(J==L)) THEN
        ACTION(I,J,K,J)
        ACTION(J,I,K,J)
        ACTION(I,J,J,K)
        ACTION(J,I,J,K)

        ACTION(K,J,I,J)
        ACTION(J,K,I,J)
        ACTION(K,J,J,I)
        ACTION(J,K,J,I)
     ELSE IF ((I>J).AND.(I>K).AND.(K==L)) THEN
        ACTION(I,J,K,K)
        ACTION(J,I,K,K)
        
        ACTION(K,K,I,J)
        ACTION(K,K,J,I)
     ELSE IF ((I==J).AND.(J>K).AND.(K>L)) THEN
        ACTION(I,I,K,L)
        ACTION(I,I,L,K)
        
        ACTION(K,L,I,I)
        ACTION(L,K,I,I)
! 4 distinct values for the 4 indices
     ELSE IF (    ((I>J).AND.(J>K).AND.(K>L)) &
              .OR.((I>K).AND.(K>J).AND.(J>L)) &
              .OR.((I>K).AND.(K>L).AND.(L>J))) THEN
        ACTION(I,J,K,L)
        ACTION(J,I,K,L)
        ACTION(I,J,L,K)
        ACTION(J,I,L,K)
        
        ACTION(K,L,I,J)
        ACTION(L,K,I,J)
        ACTION(K,L,J,I)
        ACTION(L,K,J,I)
     END IF
  END DO
  IF (.NOT.DIRECT.AND.USEDISK) CLOSE(BIUNIT)






 
