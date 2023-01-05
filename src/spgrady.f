C> @file
C> @brief Compute y-gradient in spectral space.
C> @author IREDELL @date 92-10-31

C> Computes the horizontal vector y-gradient of a scalar field
c> in spectral space.
C>
C> Subprogram speps should be called already.
C>
C> if L is the zonal wavenumber, N is the total wavenumber,
C> EPS(L,N)=SQRT((N**2-L**2)/(4*N**2-1)) and A is Earth radius,
C> then the meridional gradient of Q(L,N) is computed as
C> EPS(L,N+1)*(N+2)/A*Q(L,N+1)-EPS(L,N+1)*(N-1)/A*Q(L,N-1).
C>
C> Extra terms are computed over top of the spectral domain.
C>
C> advantage is taken of the fact that EPS(L,L)=0
C> in order to vectorize over the entire spectral domain.
C>
C> @param I spectral domain shape
c> (0 for triangular, 1 for rhomboidal)
C> @param M spectral truncation
C> @param ENN1 N*(N+1)/A**2
C> @param EON EPSILON/N*A
C> @param EONTOP EPSILON/N*A over top
C> @param Q scalar field
C> @param QDY merid gradient (times coslat)
C> @param QDYTOP merid gradient (times coslat) over top
C>
C> @author IREDELL @date 92-10-31
      SUBROUTINE SPGRADY(I,M,ENN1,EON,EONTOP,Q,QDY,QDYTOP)

      REAL ENN1((M+1)*((I+1)*M+2)/2)
      REAL EON((M+1)*((I+1)*M+2)/2),EONTOP(M+1)
      REAL Q((M+1)*((I+1)*M+2))
      REAL QDY((M+1)*((I+1)*M+2))
      REAL QDYTOP(2*(M+1))

C  TAKE MERIDIONAL GRADIENT
      K=1
      QDY(2*K-1)=EON(K+1)*ENN1(K+1)*Q(2*K+1)
      QDY(2*K)=EON(K+1)*ENN1(K+1)*Q(2*K+2)
      DO K=2,(M+1)*((I+1)*M+2)/2-1
        QDY(2*K-1)=EON(K+1)*ENN1(K+1)*Q(2*K+1)-EON(K)*ENN1(K-1)*Q(2*K-3)
        QDY(2*K)=EON(K+1)*ENN1(K+1)*Q(2*K+2)-EON(K)*ENN1(K-1)*Q(2*K-2)
      ENDDO
      K=(M+1)*((I+1)*M+2)/2
      QDY(2*K-1)=-EON(K)*ENN1(K-1)*Q(2*K-3)
      QDY(2*K)=-EON(K)*ENN1(K-1)*Q(2*K-2)

C  TAKE MERIDIONAL GRADIENT OVER TOP
      DO L=0,M
        K=L*(2*M+(I-1)*(L-1))/2+I*L+M+1
        QDYTOP(2*L+1)=-EONTOP(L+1)*ENN1(K)*Q(2*K-1)
        QDYTOP(2*L+2)=-EONTOP(L+1)*ENN1(K)*Q(2*K)
      ENDDO

      RETURN
      END
