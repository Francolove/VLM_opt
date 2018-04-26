MODULE Subroutines
!Questo modulo contiene le subroutines usate nel programma (Eccetto quella di lettura dei dati)
USE InFiles
IMPLICIT NONE

CONTAINS
   SUBROUTINE cross(a,b,c)
      REAL, DIMENSION(3), INTENT(IN) :: a, b
      REAL, DIMENSION(3), INTENT(OUT) :: c
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   END SUBROUTINE cross

   SUBROUTINE VortexSegment(xa, xb, x, gamma, v)
      REAL, DIMENSION(3), INTENT(IN) :: xa, xb, x
      REAL, INTENT(IN) :: gamma
      REAL, DIMENSION(3), INTENT(OUT) :: v

      REAL, DIMENSION(3) :: r1, r2, rab, e, r1rab
      REAL :: cosT1, cosT2, h

      r1 = x - xa; r2 = x - xb; rab = xb - xa
      CALL cross(r1, rab, r1rab)
      CALL cross(r1, r2, e)
      cosT1 = DOT_PRODUCT(r1, rab) / (NORM2(r1)*NORM2(rab))
      cosT2 = DOT_PRODUCT(r2, rab) / (NORM2(r2)*NORM2(rab))
      h = NORM2(r1rab) / NORM2(rab)

      IF (NORM2(r1rab) < 1E-5) THEN
         v = 0.
      ELSE
         v = gamma/(16.*ATAN(1.) * h) * (cosT1 - cosT2) * e / NORM2(e)
      ENDIF
   END SUBROUTINE VortexSegment

   SUBROUTINE VortexRing(xa, xb, xc, xd, x, gamma, v)
      REAL, DIMENSION(3), INTENT(IN) :: xa, xb, xc, xd, x
      REAL, INTENT(IN) :: gamma
      REAL, DIMENSION(3), INTENT(OUT) :: v

      REAL, DIMENSION(3) :: v1, v2, v3, v4

      CALL VortexSegment(xa, xb, x, gamma, v1)
      CALL VortexSegment(xb, xc, x, gamma, v2)
      CALL VortexSegment(xc, xd, x, gamma, v3)
      CALL VortexSegment(xd, xa, x, gamma, v4)
      v = v1 + v2 + v3 + v4
   END SUBROUTINE VortexRing

   SUBROUTINE Velocita(Vortices, Wake, x, m, n, v)
      IMPLICIT NONE
      TYPE(ptrVortici), DIMENSION(:), POINTER,INTENT(IN) :: Vortices, Wake
      REAL, DIMENSION(3), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: m, n
      REAL, DIMENSION(3), INTENT(OUT) :: v

      INTEGER :: j
      REAL, DIMENSION(3) :: xa, xb, xc, xd, &
                            u, u_, uw, uw_
      u = 0.; u_ = 0.; uw = 0.; uw_ = 0.; v = 0.
      DO j=1,m*n
         xa = Vortices(j)%ptrVortex%Vertici(1,1,:); xd = Vortices(j)%ptrVortex%Vertici(1,2,:)
         xb = Vortices(j)%ptrVortex%Vertici(2,1,:); xc = Vortices(j)%ptrVortex%Vertici(2,2,:)
         CALL VortexRing(xa, xb, xc, xd, x, Vortices(j)%ptrVortex%Gamma, u)
         xa(2) = -xa(2); xb(2) = -xb(2); xc(2) = -xc(2); xd(2) = -xd(2)
         CALL VortexRing(xa, xb, xc, xd, x, Vortices(j)%ptrVortex%Gamma, u_)
         IF (j >= n*(m-1)+1) THEN
            xa = Wake(j)%ptrVortex%Vertici(1,1,:); xd = Wake(j)%ptrVortex%Vertici(1,2,:)
            xb = Wake(j)%ptrVortex%Vertici(2,1,:); xc = Wake(j)%ptrVortex%Vertici(2,2,:)
            CALL VortexRing(xa, xb, xc, xd, x, Wake(j)%ptrVortex%Gamma, uw)
            xa(2) = -xa(2); xb(2) = -xb(2); xc(2) = -xc(2); xd(2) = -xd(2)
            CALL VortexRing(xa, xb, xc, xd, x, Wake(j)%ptrVortex%Gamma, uw_)
         ELSE
            uw = 0.
            uw_ = 0.
         ENDIF
         v = v + u - u_ + uw - uw_
         u = 0.; u_ = 0.; uw = 0.; uw_ = 0.
      ENDDO

   END SUBROUTINE Velocita

   SUBROUTINE VortexVertices(ptr0, ptr1, m, n)
      REAL, DIMENSION(:,:,:), POINTER :: ptr0, ptr1
      INTEGER, INTENT(IN) :: m, n

      INTEGER :: ptrStat, i

      ALLOCATE(ptr1(1:m+1, 1:n+1, 3), STAT=ptrStat)
      IF (ptrStat == 0) THEN
         DO i = 1, m

            ptr1(i,:,1) = ptr0(i,:,1) + 0.25*(ptr0(i+1,:,1)-ptr0(i,:,1))
            ptr1(i,:,2) = ptr0(i,:,2) + (ptr0(i+1,:,2)-ptr0(i,:,2)) / (ptr0(i+1,:,1)-ptr0(i,:,1))* &
                         (ptr1(i,:,1)-ptr0(i,:,1))
            ptr1(i,:,3) = ptr0(i,:,3) + (ptr0(i+1,:,3)-ptr0(i,:,3)) / (ptr0(i+1,:,1)-ptr0(i,:,1))* &
                         (ptr1(i,:,1)-ptr0(i,:,1))
         ENDDO
         ptr1(m+1,:,1) = ptr0(m+1,:,1) + 0.25*(ptr0(m+1,:,1)-ptr0(m,:,1))
         ptr1(m+1,:,2) = ptr0(m+1,:,2)
         ptr1(m+1,:,3) = ptr0(m+1,:,3)
      ELSE
         WRITE(*,*) 'Unable to allocate the pointer'
         STOP
      ENDIF
   END SUBROUTINE VortexVertices

   SUBROUTINE PanelPctrl(ptr0, ptr1, m, n)
      REAL, DIMENSION(:,:,:), POINTER :: ptr0, ptr1
      INTEGER, INTENT(IN) :: m, n

      INTEGER :: ptrStat, i, j
      REAL, DIMENSION(:,:,:), POINTER :: ptr2

      ALLOCATE(ptr1(1:m, 1:n, 3), ptr2(1:m, 1:n+1, 3), STAT=ptrStat)
      IF (ptrStat == 0) THEN
         DO i = 1, m
            ptr2(i,:,1) = ptr0(i,:,1) + 0.75*(ptr0(i+1,:,1)-ptr0(i,:,1))
            ptr2(i,:,2) = ptr0(i,:,2) + (ptr0(i+1,:,2)-ptr0(i,:,2)) / (ptr0(i+1,:,1)-ptr0(i,:,1))* &
                         (ptr2(i,:,1)-ptr0(i,:,1))
            ptr2(i,:,3) = ptr0(i,:,3) + (ptr0(i+1,:,3)-ptr0(i,:,3)) / (ptr0(i+1,:,1)-ptr0(i,:,1))* &
                         (ptr2(i,:,1)-ptr0(i,:,1))
         ENDDO
         DO j = 1, n
            ptr1(:,j,1) = (ptr2(:,j,1) + ptr2(:,j+1,1))/2
            ptr1(:,j,2) = (ptr2(:,j,2) + ptr2(:,j+1,2))/2
            ptr1(:,j,3) = (ptr2(:,j,3) + ptr2(:,j+1,3))/2
         ENDDO
         DEALLOCATE(ptr2)
      ELSE
         WRITE(*,*) 'Unable to allocate the pointer'
         STOP
      ENDIF
   END SUBROUTINE PanelPctrl

   SUBROUTINE GeneraPannelli(ptr0, Panels, m, n)

      REAL, DIMENSION(:,:,:), POINTER :: ptr0
      TYPE(ptrPannelli), DIMENSION(:), POINTER :: Panels
      INTEGER, INTENT(IN) :: m, n

      INTEGER :: i, j, k, ptrStat
      REAL, DIMENSION(:,:,:), POINTER :: ptr1

      CALL PanelPctrl(ptr0, ptr1, m, n)

      ALLOCATE(Panels(1:m*n), STAT=ptrStat)

      IF (ptrStat == 0) THEN
         k = 1
         DO i=1,m
            DO j=1,n
               ALLOCATE(Panels(k)%ptrPanel)
               CALL cross( (ptr0(i+1,j+1,:) - ptr0(i,j,:)), (ptr0(i,j+1,:) - ptr0(i+1,j,:)), &
                            Panels(k)%ptrPanel%N )
               Panels(k)%ptrPanel%N = Panels(k)%ptrPanel%N / NORM2(Panels(k)%ptrPanel%N)
               Panels(k)%ptrPanel%Pctrl = ptr1(i,j,:)
               k = k + 1
            ENDDO
         ENDDO
      ELSE
         WRITE(*,*) 'Unable to allocate the pointer'
         STOP
      ENDIF

   END SUBROUTINE GeneraPannelli

   SUBROUTINE GeneraVortici(ptr0, Vortices, m, n)
      REAL, DIMENSION(:,:,:), POINTER :: ptr0
      TYPE(ptrVortici), DIMENSION(:), POINTER :: Vortices
      INTEGER, INTENT(IN) :: m, n

      INTEGER :: i, j, k, ptrStat
      REAL, DIMENSION(:,:,:), POINTER :: ptr1

      CALL VortexVertices(ptr0, ptr1, m, n)

      ALLOCATE(Vortices(1:m*n), STAT=ptrStat)

      IF (ptrStat == 0) THEN
         k = 1
         DO i=1,m
            DO j=1,n
               ALLOCATE(Vortices(k)%ptrVortex)
               Vortices(k)%ptrVortex%vertici(1,1,:) = ptr1(i,j,:)
               Vortices(k)%ptrVortex%vertici(1,2,:) = ptr1(i,j+1,:)
               Vortices(k)%ptrVortex%vertici(2,2,:) = ptr1(i+1,j+1,:)
               Vortices(k)%ptrVortex%vertici(2,1,:) = ptr1(i+1,j,:)
               k = k + 1
            ENDDO
         ENDDO
      ELSE
         WRITE(*,*) 'Unable to allocate the pointer'
         STOP
      ENDIF
   END SUBROUTINE GeneraVortici

!   SUBROUTINE CoeffVel(Vortices, Wake, m, n, U, alpha, v1, v2, v3, v4, &
 !                      p1, p2, p3, p4)! SemiSpan, AR, Cf, Cl, Cdi)
!      TYPE(ptrVortici), DIMENSION(:), POINTER :: Vortices, Wake
!      INTEGER, INTENT(IN) :: m, n
!      REAL, DIMENSION(3), INTENT(IN) :: U
!      REAL, INTENT(IN) :: alpha!, SemiSpan, AR
!      REAL, DIMENSION(m*n,3), INTENT(OUT) :: v1, v2, v3, v4, p1, p2, p3, p4!Cf
!      !REAL, INTENT(OUT) :: Cl, Cdi

!      INTEGER :: i
!      REAL, DIMENSION(3) :: f, u1, u2, u3, u4, &
!                            x1, x2, x3, x4, f1, &
!                            f2, f3, f4

      !f = 0.; f1 = 0.; f2 = 0.; f3 = 0.; f4 = 0.
!      u1 = 0.; u2 = 0.; u3 = 0.; u4 = 0.
!      x1 = 0.; x2 = 0.; x3 = 0.; x4 = 0.
!      DO i=1,m*n
!         x1 = (Vortices(i)%ptrVortex%Vertici(1,1,:) + Vortices(i)%ptrVortex%Vertici(2,1,:)) / 2
!         x2 = (Vortices(i)%ptrVortex%Vertici(2,1,:) + Vortices(i)%ptrVortex%Vertici(2,2,:)) / 2
!         x3 = (Vortices(i)%ptrVortex%Vertici(2,2,:) + Vortices(i)%ptrVortex%Vertici(1,2,:)) / 2
!         x4 = (Vortices(i)%ptrVortex%Vertici(1,2,:) + Vortices(i)%ptrVortex%Vertici(1,1,:)) / 2
!         CALL Velocita(Vortices, Wake, x1, m, n, u1)
!         CALL Velocita(Vortices, Wake, x2, m, n, u2)
!         CALL Velocita(Vortices, Wake, x3, m, n, u3)
!         CALL Velocita(Vortices, Wake, x4, m, n, u4)

!         v1(i,:) = u1; v2(i,:)=u2; v3(i,:)=u3; v4(i,:)=u4
         !p1(i,:)=x1; p2(i,:)=x2; p3(i,:)=x3; p4(i,:)=x4
!         u1 = 0.; u2 = 0.; u3 = 0.; u4 = 0.
!         x1 = 0.; x2 = 0.; x3 = 0.; x4 = 0.
!         x1 = Vortices(i)%ptrVortex%Vertici(2,1,:) - Vortices(i)%ptrVortex%Vertici(1,1,:)
!         x2 = Vortices(i)%ptrVortex%Vertici(2,2,:) - Vortices(i)%ptrVortex%Vertici(2,1,:)
!         x3 = Vortices(i)%ptrVortex%Vertici(1,2,:) - Vortices(i)%ptrVortex%Vertici(2,2,:)
!         x4 = Vortices(i)%ptrVortex%Vertici(1,1,:) - Vortices(i)%ptrVortex%Vertici(1,2,:)
!         p1(i,:)=x1; p2(i,:)=x2; p3(i,:)=x3; p4(i,:)=x4
         !CALL cross(U+u1, x1, f1); CALL cross(U+u2, x2, f2)
         !CALL cross(U+u3, x3, f3); CALL cross(U+u4, x4, f4)

         !f = f + Vortices(i)%ptrVortex%Gamma * (f1 + f2 + f3 + f4)
         !f1 = 0.; f2 = 0.; f3 = 0.; f4 = 0.
!         x1 = 0.; x2 = 0.; x3 = 0.; x4 = 0.
!      ENDDO
      !Cf = AR/(NORM2(U)*SemiSpan)**2 * f
      !Cl = Cf(3) * COS(alpha) - Cf(1) * SIN(alpha)
      !Cdi = Cf(1) * COS(alpha) + Cf(3) * SIN(alpha)
!   END SUBROUTINE CoeffVel

SUBROUTINE AerCoeff(Vortices, Wake, m, n, U, alpha, SemiSpan, AR, Cf, Cl, Cdi)
      TYPE(ptrVortici), DIMENSION(:), POINTER :: Vortices, Wake
      INTEGER, INTENT(IN) :: m, n
      REAL, DIMENSION(3), INTENT(IN) :: U
      REAL, INTENT(IN) :: alpha, SemiSpan, AR
      REAL, DIMENSION(3), INTENT(OUT) :: Cf
      REAL, INTENT(OUT) :: Cl, Cdi

      INTEGER :: i
      REAL, DIMENSION(3) :: f, u1, u2, u3, u4, &
                            x1, x2, x3, x4, f1, &
                            f2, f3, f4

      f = 0.; f1 = 0.; f2 = 0.; f3 = 0.; f4 = 0.
      u1 = 0.; u2 = 0.; u3 = 0.; u4 = 0.
      x1 = 0.; x2 = 0.; x3 = 0.; x4 = 0.
      cl = 0
!      OPEN (UNIT=800, FILE='m_pt', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
!      OPEN (UNIT=900, FILE='v_m_pt', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)

      DO i=1,m*n
         if (i <= m*n-n) THEN
            x1 = (Vortices(i)%ptrVortex%Vertici(1,1,:) + Vortices(i)%ptrVortex%Vertici(2,1,:)) / 2
            x2 = (Vortices(i)%ptrVortex%Vertici(2,1,:) + Vortices(i)%ptrVortex%Vertici(2,2,:)) / 2
            x3 = (Vortices(i)%ptrVortex%Vertici(2,2,:) + Vortices(i)%ptrVortex%Vertici(1,2,:)) / 2
            x4 = (Vortices(i)%ptrVortex%Vertici(1,2,:) + Vortices(i)%ptrVortex%Vertici(1,1,:)) / 2

!            write(800,*) x1
!            write(800,*) x2
!            write(800,*) x3
!            write(800,*) x4


            CALL Velocita(Vortices, Wake, x1, m, n, u1)
            CALL Velocita(Vortices, Wake, x2, m, n, u2)
            CALL Velocita(Vortices, Wake, x3, m, n, u3)
            CALL Velocita(Vortices, Wake, x4, m, n, u4)

!            write(900,*) u1+U
!            write(900,*) u2+U
!            write(900,*) u3+U
!            write(900,*) u4+U

            x1 = Vortices(i)%ptrVortex%Vertici(2,1,:) - Vortices(i)%ptrVortex%Vertici(1,1,:)
            x2 = Vortices(i)%ptrVortex%Vertici(2,2,:) - Vortices(i)%ptrVortex%Vertici(2,1,:)
            x3 = Vortices(i)%ptrVortex%Vertici(1,2,:) - Vortices(i)%ptrVortex%Vertici(2,2,:)
            x4 = Vortices(i)%ptrVortex%Vertici(1,1,:) - Vortices(i)%ptrVortex%Vertici(1,2,:)

            CALL cross(U+u1, x1, f1); CALL cross(U+u2, x2, f2)
            CALL cross(U+u3, x3, f3); CALL cross(U+u4, x4, f4)

            f = f + Vortices(i)%ptrVortex%Gamma * (f1 + f2 + f3 + f4)

         ELSE
            x1 = (Vortices(i)%ptrVortex%Vertici(1,1,:) + Vortices(i)%ptrVortex%Vertici(2,1,:)) / 2
            x3 = (Vortices(i)%ptrVortex%Vertici(2,2,:) + Vortices(i)%ptrVortex%Vertici(1,2,:)) / 2
            x4 = (Vortices(i)%ptrVortex%Vertici(1,2,:) + Vortices(i)%ptrVortex%Vertici(1,1,:)) / 2

!            write(800,*) x1
!            write(800,*) x3
!            write(800,*) x4

            CALL Velocita(Vortices, Wake, x1, m, n, u1)
            CALL Velocita(Vortices, Wake, x3, m, n, u3)
            CALL Velocita(Vortices, Wake, x4, m, n, u4)

!            write(900,*) u1+U
!            write(900,*) u3+U
!            write(900,*) u4+U

            x1 = Vortices(i)%ptrVortex%Vertici(2,1,:) - Vortices(i)%ptrVortex%Vertici(1,1,:)
            x3 = Vortices(i)%ptrVortex%Vertici(1,2,:) - Vortices(i)%ptrVortex%Vertici(2,2,:)
            x4 = Vortices(i)%ptrVortex%Vertici(1,1,:) - Vortices(i)%ptrVortex%Vertici(1,2,:)

            CALL cross(U+u1, x1, f1);
            CALL cross(U+u3, x3, f3); CALL cross(U+u4, x4, f4)

            f = f + Vortices(i)%ptrVortex%Gamma * (f1 + f3 + f4)
         ENDIF
         f1 = 0.; f2 = 0.; f3 = 0.; f4 = 0.
         x1 = 0.; x2 = 0.; x3 = 0.; x4 = 0.
      ENDDO
      Cf = AR/(NORM2(U)**2 *SemiSpan**2) * f
      Cl = Cf(3) * COS(alpha) - Cf(1) * SIN(alpha)
      Cdi = Cf(1) * COS(alpha) + Cf(3) * SIN(alpha)
   END SUBROUTINE AerCoeff

END MODULE Subroutines
