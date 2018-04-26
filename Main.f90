PROGRAM Main

USE Subroutines
USE InFiles
IMPLICIT NONE

REAL, DIMENSION(:,:,:), POINTER :: ptrNodes
type(ptrPannelli), DIMENSION(:), POINTER :: Pannelli
type(ptrVortici), DIMENSION(:), POINTER :: Anelli, Scia
REAL, DIMENSION(:), POINTER :: u, u_, uw, uw_, x, xa, xb, xc, xd
REAL, DIMENSION(:,:), POINTER :: A
REAL, DIMENSION(:), POINTER :: RHS, pivot
INTEGER :: info
REAL, DIMENSION(:,:), POINTER :: v1,v2,v3,v4, p1,p2,p3,p4

!Leggo Data e le coordinate dei vertici dei pannelli (queste le salvo in Nodes(:,:,1..3))
CALL SaveData

!Trasformo alpha in radianti e trovo Uinf
alpha = alpha * 4*ATAN(1.)/180
Uinf = 1E-5 * Re/b * (/COS(alpha), 0., SIN(alpha)/)

!Genero le strutture di pannelli e vortici ad anello
ptrNodes => Nodes
CALL GeneraPannelli(ptrNodes, Pannelli, m, n)
CALL GeneraVortici(ptrNodes, Anelli, m, n)
!Genero la scia a mano
ALLOCATE(Scia((m*n-n + 1):m*n))  !Gli indici di Scia non partono da 1!!!!!!

DO j=(m*n-n + 1), m*n
   ALLOCATE(Scia(j)%ptrVortex)
   Scia(j)%ptrVortex%Vertici(1,1,:) = Anelli(j)%ptrVortex%Vertici(2,1,:)
   Scia(j)%ptrVortex%Vertici(1,2,:) = Anelli(j)%ptrVortex%Vertici(2,2,:)
   Scia(j)%ptrVortex%Vertici(2,2,:) = Anelli(j)%ptrVortex%Vertici(2,2,:) + Uinf * 10.
   Scia(j)%ptrVortex%Vertici(2,1,:) = Anelli(j)%ptrVortex%Vertici(2,1,:) + Uinf * 10.
ENDDO

!Echo
OPEN (UNIT=100, FILE='Pc', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
!OPEN (UNIT=110, FILE='N', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
OPEN (UNIT=120, FILE='Xvor', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
!OPEN (UNIT=130, FILE='Xwake', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
DO i=1,m*n
      WRITE(100,*) Pannelli(i)%ptrPanel%Pctrl
!      WRITE(110,*) Pannelli(i,j)%ptrPanel%N
      WRITE(120,*) Anelli(i)%ptrVortex%Vertici(1,1,:)!, Anelli(i,j)%ptrVortex%Vertici(1,2,:)
!      IF (i==m) WRITE(130,*) Scia(j)%ptrVortex%Vertici(2,1,:), Scia(j)%ptrVortex%Vertici(2,2,:)
ENDDO

!Costruisco la matrice dei coefficienti d'influenza e il termine noto
ALLOCATE(A(1:m*n, 1:m*n), RHS(1:m*n), x(3), &
         u(3), u_(3), uw(3), uw_(3), xa(3), &
         xb(3), xc(3), xd(3))
u = 0.; u_ = 0.; uw = 0.; uw_ = 0.
DO i=1, m*n
   x = Pannelli(i)%ptrPanel%Pctrl
   DO j=1, m*n
      xa = Anelli(j)%ptrVortex%Vertici(1,1,:); xd = Anelli(j)%ptrVortex%Vertici(1,2,:)
      xb = Anelli(j)%ptrVortex%Vertici(2,1,:); xc = Anelli(j)%ptrVortex%Vertici(2,2,:)
      CALL VortexRing(xa, xb, xc, xd, x, 1., u)
      xa(2) = -xa(2); xb(2) = -xb(2); xc(2) = -xc(2); xd(2) = -xd(2)
      CALL VortexRing(xa, xb, xc, xd, x, 1., u_)
      IF (j >= n*(m-1)+1) THEN
         xa = Scia(j)%ptrVortex%Vertici(1,1,:); xd = Scia(j)%ptrVortex%Vertici(1,2,:)
         xb = Scia(j)%ptrVortex%Vertici(2,1,:); xc = Scia(j)%ptrVortex%Vertici(2,2,:)
         CALL VortexRing(xa, xb, xc, xd, x, 1., uw)
         xa(2) = -xa(2); xb(2) = -xb(2); xc(2) = -xc(2); xd(2) = -xd(2)
         CALL VortexRing(xa, xb, xc, xd, x, 1., uw_)
      ELSE
         uw = 0.
         uw_ = 0.
      ENDIF
      u = u - u_ + uw - uw_
      A(i,j) = DOT_PRODUCT(u, Pannelli(i)%ptrPanel%N)
      u = 0.; u_ = 0.; uw = 0.; uw_ = 0.
   ENDDO
   RHS(i) = -DOT_PRODUCT(Uinf, Pannelli(i)%ptrPanel%N)
ENDDO
DEALLOCATE(u_, uw, uw_, xa, xb, xc, xd, x)

!Risolvo sistema lineare
allocate(pivot(1:m*n))
CALL SGESV(m*n, 1, A, m*n,pivot, RHS, m*n, info)
write(*,*) 'info=', info

!Distribuisco la circolazione nei vari vortici
DO i=1,m*n
   Anelli(i)%ptrVortex%Gamma = RHS(i)
   IF (i >= n*(m-1)+1) THEN
      Scia(i)%ptrVortex%Gamma = RHS(i)
   ENDIF
ENDDO

!Calcolo Velocità in vari punti
OPEN (UNIT=140, FILE='Vel', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
OPEN (UNIT=150, FILE='VelVor', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
u = 0.
DO i=1,m*n
   CALL Velocita(Anelli, Scia, Pannelli(i)%ptrPanel%Pctrl, m, n, u)
   WRITE(140,*) u + Uinf
   CALL Velocita(Anelli, Scia, Anelli(i)%ptrVortex%Vertici(1,1,:), m, n, u)
   WRITE(150,*) u + Uinf
   u = 0.
ENDDO
DEALLOCATE(u, A, RHS)

!Calcolo coefficienti aerodinamici
CALL AerCoeff(Anelli, Scia, m, n, Uinf, alpha, b, AR, Cf, Cl, Cdi)
OPEN (UNIT=160, FILE='Coeff', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
WRITE(*,*) Cf
WRITE(160,*) Cl
WRITE(160,*) Cdi

!ALLOCATE(v1(m*n,3), v2(m*n,3), v3(m*n,3), v4(m*n,3), p1(m*n,3), p2(m*n,3), p3(m*n,3), p4(m*n,3))
!CALL CoeffVel(Anelli, Scia, m, n, Uinf, alpha, v1, v2, v3, v4, p1, p2, p3, p4)! b, AR, Cf, Cl, Cdi)
!OPEN (UNIT=161, FILE='v1', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
!OPEN (UNIT=162, FILE='v2', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
!OPEN (UNIT=163, FILE='v3', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
!OPEN (UNIT=164, FILE='v4', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=OpenStat)
!DO i=1,m*n

!   write(161,*) p1(i,:), Uinf + v1(i,:)
!   write(162,*) p2(i,:), Uinf + v2(i,:)
!   write(163,*) p3(i,:), Uinf + v3(i,:)
!   write(164,*) p4(i,:), Uinf + v4(i,:)
!ENDDO

END PROGRAM Main
