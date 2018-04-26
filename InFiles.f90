MODULE InFiles
!Modulo contenente la dichiarazione delle variabili globali e di
!DERIVED data TYPEs.
!Legge i dati dai file e inizializza le variabili.

IMPLICIT NONE
SAVE

TYPE :: panel
   REAL, DIMENSION(3) :: N !Normale al pannello
   REAL, DIMENSION(3) :: Pctrl !Punto di controllo
   REAL :: Area   !Area del pannello
END TYPE panel

TYPE :: vortex
   REAL, DIMENSION(2,2,3) :: vertici   !Array contenente i vertici dell'anello
   REAL :: Gamma
END TYPE vortex

TYPE :: ptrPannelli
   TYPE(panel), POINTER :: ptrPanel !Pointer di tipo panel
END TYPE ptrPannelli

TYPE :: ptrVortici
   TYPE(vortex), POINTER :: ptrVortex !Pointer di tipo vortex
END TYPE ptrVortici

INTEGER :: m, n, Re
REAL :: b, alpha, Cl, Cdi, AR
REAL, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: Nodes   !Contiene le coordinate dei pannelli
REAL, DIMENSION(3) :: Uinf, Cf
!REAL, DIMENSION(m*n, m*n) :: A
!REAL, DIMENSION(m*n) :: RHS

!Variabili super stupide
INTEGER :: OpenStat, AllocStat, i, j


CONTAINS

!Reding stuff
   SUBROUTINE SaveData
      !Lettura file Data.dat
      OPEN (UNIT=10, FILE='Data.dat', STATUS='UNKNOWN', ACTION='READ', IOSTAT=OpenStat)
      IF (OpenStat == 0) THEN
         READ(10,100, IOSTAT=OpenStat) m, n, b, AR, alpha, Re
         100 FORMAT (2x, I4, /, 2x, I4, /, 2x, F4.2, /, 3x, F4.2, /, 6x, F10.6, /, 3x, I8)
         !Echo
         WRITE(*,*) 'DATA'
         WRITE(*,*) '===================================================================='
         WRITE(*,110) m, n, m*n
         110 FORMAT (1x, 'Chord panels = ', I4, 2x, '(Half)Span panels = ', I4, 2x, &
                     'Total panels = ', I6)
         WRITE(*,120) alpha, b, AR, Re
         120 FORMAT (1x, 'Alpha = ', F6.2,'[deg]', /, '(Half)Span = ', F10.6, '[m]', /, &
                    'Aspect Ratio = ', F4.2, /, 'Re = ', I8)
         
         ALLOCATE(Nodes(1:m+1, 1:n+1, 3), STAT=AllocStat)
         IF (Allocstat /= 0 ) WRITE(*,*) 'Problema allocazione Nodes'
      ELSE
         WRITE(*,*) 'Riscontrato problema nella lettura del file Data.dat!'
      ENDIF
      !Lettura file X.dat
      OPEN (UNIT=20, FILE='X.dat', STATUS='UNKNOWN', ACTION='READ', IOSTAT=OpenStat)
      IF (OpenStat == 0) THEN
         i = 1
         DO
         READ(20,*, IOSTAT=OpenStat) Nodes(i,:,1)
         IF (OpenStat /= 0) EXIT
         i = i + 1
         ENDDO
      ELSE
         WRITE(*,*) 'Riscontrato problema nella lettura del file X.dat!'
      ENDIF
      !Letura file Y.dat
      OPEN (UNIT=30, FILE='Y.dat', STATUS='UNKNOWN', ACTION='READ', IOSTAT=OpenStat)
      IF (OpenStat == 0) THEN
         i = 1
         DO
         READ(30,*, IOSTAT=OpenStat) Nodes(i,:,2)
         IF (OpenStat /= 0) EXIT
         i = i + 1
         ENDDO
      ELSE
         WRITE(*,*) 'Riscontrato problema nella lettura del file Y.dat!'
      ENDIF
      !Lettura file Z.dat
      OPEN (UNIT=40, FILE='Z.dat', STATUS='UNKNOWN', ACTION='READ', IOSTAT=OpenStat)
      IF (OpenStat == 0) THEN
         i = 1
         DO
         READ(40,*, IOSTAT=OpenStat) Nodes(i,:,3)
         IF (OpenStat /= 0) EXIT
         i = i + 1
         ENDDO
      ELSE
         WRITE(*,*) 'Riscontrato problema nella lettura del file Z.dat!'
      ENDIF             

   END SUBROUTINE SaveData



END MODULE InFiles
