C> @file
C> @brief Perform multiple fast fourier transforms.
C> @author IREDELL @date 96-02-20

C> This subprogram performs multiple fast fourier transforms
C> between complex amplitudes in fourier space and real values
C> in cyclic physical space.
C>
C> Subprogram spfft must be invoked first with idir=0
C> to initialize trigonemetric data. Use subprogram spfft1
C> to perform an fft without previous initialization.
C> This version invokes the ibm essl fft.
C>
C> The restrictions on imax are that it must be a multiple
C> of 1 to 25 factors of two, up to 2 factors of three,
C> and up to 1 factor of five, seven and eleven.
C>
C> If IDIR=0, then W and G need not contain any valid data.
C> the other parameters must be supplied and cannot change
C> in succeeding calls until the next time it is called with IDIR=0.
C>
C> This subprogram is not thread-safe when IDIR=0. On the other hand,
C> when IDIR is not zero, it can be called from a threaded region.
C>
C> @param IMAX NUMBER OF VALUES IN THE CYCLIC PHYSICAL SPACE
C> (SEE LIMITATIONS ON IMAX IN REMARKS BELOW.)
C> @param INCW FIRST DIMENSION OF THE COMPLEX AMPLITUDE ARRAY
C> (INCW >= IMAX/2+1)
C> @param INCG FIRST DIMENSION OF THE REAL VALUE ARRAY
C> (INCG >= IMAX)
C> @param KMAX NUMBER OF TRANSFORMS TO PERFORM
C> @param[out] W COMPLEX AMPLITUDES IF IDIR>0
C> @param[out] G REAL VALUES IF IDIR<0
C> @param IDIR DIRECTION FLAG
C> - IDIR=0 TO INITIALIZE INTERNAL TRIGONOMETRIC DATA
C> - IDIR>0 TO TRANSFORM FROM FOURIER TO PHYSICAL SPACE
C> - IDIR<0 TO TRANSFORM FROM PHYSICAL TO FOURIER SPACE
C>
C> @author IREDELL @date 96-02-20
      SUBROUTINE SPFFT(IMAX,INCW,INCG,KMAX,W,G,IDIR)

        IMPLICIT NONE
        INTEGER,INTENT(IN):: IMAX,INCW,INCG,KMAX,IDIR
        COMPLEX,INTENT(INOUT):: W(INCW,KMAX)
        REAL,INTENT(INOUT):: G(INCG,KMAX)
        INTEGER,SAVE:: NAUX1=0
        REAL,SAVE,ALLOCATABLE:: AUX1CR(:),AUX1RC(:)
        INTEGER:: NAUX2
        REAL:: AUX2(20000+INT(0.57*IMAX))

        NAUX2=20000+INT(0.57*IMAX)

C  INITIALIZATION.
C  ALLOCATE AND FILL AUXILIARY ARRAYS WITH TRIGONOMETRIC DATA
        SELECT CASE(IDIR)
          CASE(0)
            IF(NAUX1.GT.0) DEALLOCATE(AUX1CR,AUX1RC)
            NAUX1=25000+INT(0.82*IMAX)
            ALLOCATE(AUX1CR(NAUX1),AUX1RC(NAUX1))
            CALL SCRFT(1,W,INCW,G,INCG,IMAX,KMAX,-1,1.,
     &                 AUX1CR,NAUX1,AUX2,NAUX2,0.,0)
            CALL SRCFT(1,G,INCG,W,INCW,IMAX,KMAX,+1,1./IMAX,
     &                 AUX1RC,NAUX1,AUX2,NAUX2,0.,0)

C  FOURIER TO PHYSICAL TRANSFORM.
          CASE(1:)
            CALL SCRFT(0,W,INCW,G,INCG,IMAX,KMAX,-1,1.,
     &                 AUX1CR,NAUX1,AUX2,NAUX2,0.,0)

C  PHYSICAL TO FOURIER TRANSFORM.
          CASE(:-1)
            CALL SRCFT(0,G,INCG,W,INCW,IMAX,KMAX,+1,1./IMAX,
     &                 AUX1RC,NAUX1,AUX2,NAUX2,0.,0)
        END SELECT
      END SUBROUTINE
