C> @file
C>
C> Compute gradient in spectral space.
C> @author IREDELL @date 92-10-31

C> Computes the horizontal vector gradient of a scalar field
c> in spectral space.
C>
C> Subprogram speps() should be called already.
C>
C> If l is the zonal wavenumber, n is the total wavenumber,
c> eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) and a is earth radius,
c> then the zonal gradient of q(l,n) is simply i*l/a*q(l,n)
c> while the meridional gradient of q(l,n) is computed as
c> eps(l,n+1)*(n+2)/a*q(l,n+1)-eps(l,n+1)*(n-1)/a*q(l,n-1).
C>
C> Extra terms are computed over top of the spectral domain.
C>
C> Advantage is taken of the fact that eps(l,l)=0
c> in order to vectorize over the entire spectral domain.
C>
C> @param I SPECTRAL DOMAIN SHAPE (0 FOR TRIANGULAR, 1 FOR RHOMBOIDAL)
C> @param M SPECTRAL TRUNCATION
C> @param ENN1 ((M+1)*((I+1)*M+2)/2) N*(N+1)/A**2
C> @param ELONN1 ((M+1)*((I+1)*M+2)/2) L/(N*(N+1))*A
C> @param EON ((M+1)*((I+1)*M+2)/2) EPSILON/N*A
C> @param EONTOP (M+1) EPSILON/N*A OVER TOP
C> @param Q ((M+1)*((I+1)*M+2)) SCALAR FIELD
C> @param QDX ((M+1)*((I+1)*M+2)) ZONAL GRADIENT (TIMES COSLAT)
C> @param QDY ((M+1)*((I+1)*M+2)) MERID GRADIENT (TIMES COSLAT)
C> @param QDYTOP (2*(M+1)) MERID GRADIENT (TIMES COSLAT) OVER TOP
      SUBROUTINE SPGRADQ(I,M,ENN1,ELONN1,EON,EONTOP,Q,QDX,QDY,QDYTOP)

      REAL ENN1((M+1)*((I+1)*M+2)/2),ELONN1((M+1)*((I+1)*M+2)/2)
      REAL EON((M+1)*((I+1)*M+2)/2),EONTOP(M+1)
      REAL Q((M+1)*((I+1)*M+2))
      REAL QDX((M+1)*((I+1)*M+2)),QDY((M+1)*((I+1)*M+2))
      REAL QDYTOP(2*(M+1))
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TAKE ZONAL AND MERIDIONAL GRADIENTS
      K=1
      QDX(2*K-1)=0.
      QDX(2*K)=0.
      QDY(2*K-1)=EON(K+1)*ENN1(K+1)*Q(2*K+1)
      QDY(2*K)=EON(K+1)*ENN1(K+1)*Q(2*K+2)
      DO K=2,(M+1)*((I+1)*M+2)/2-1
        QDX(2*K-1)=-ELONN1(K)*ENN1(K)*Q(2*K)
        QDX(2*K)=ELONN1(K)*ENN1(K)*Q(2*K-1)
        QDY(2*K-1)=EON(K+1)*ENN1(K+1)*Q(2*K+1)-EON(K)*ENN1(K-1)*Q(2*K-3)
        QDY(2*K)=EON(K+1)*ENN1(K+1)*Q(2*K+2)-EON(K)*ENN1(K-1)*Q(2*K-2)
      ENDDO
      K=(M+1)*((I+1)*M+2)/2
      QDX(2*K-1)=-ELONN1(K)*ENN1(K)*Q(2*K)
      QDX(2*K)=ELONN1(K)*ENN1(K)*Q(2*K-1)
      QDY(2*K-1)=-EON(K)*ENN1(K-1)*Q(2*K-3)
      QDY(2*K)=-EON(K)*ENN1(K-1)*Q(2*K-2)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TAKE MERIDIONAL GRADIENT OVER TOP
      DO L=0,M
        K=L*(2*M+(I-1)*(L-1))/2+I*L+M+1
        QDYTOP(2*L+1)=-EONTOP(L+1)*ENN1(K)*Q(2*K-1)
        QDYTOP(2*L+2)=-EONTOP(L+1)*ENN1(K)*Q(2*K)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
