include 'mpif.h'
DOUBLE PRECISION :: xmin=-2, xmax=1, ymin=-1.5, ymax=1.5, h=0.001, x, y
COMPLEX*16 :: f, f0
INTEGER :: hx, hy, htot, i, j, k, nmax=255
INTEGER :: nprocs, rank, ierr
INTEGER, DIMENSION(:,:), ALLOCATABLE :: map, maps, n, ns

hx=((xmax-xmin)/h+1)
hy=((ymax-ymin)/h+1)
htot=hx*hy

ALLOCATE(map(hx, hy))
ALLOCATE(maps(hx, hy))
ALLOCATE(n(hx,hy))
ALLOCATE(ns(hx,hy))

map=0
maps=0
n=0
ns=0

CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

DO i=rank+1,hx,nprocs
  x=(i-1)*h+xmin
  DO j=1,hy
    y=(j-1)*h+ymin
    f0=CMPLX(x,y)
    f=0
    n(i,j)=nmax
    DO k=1, nmax
      f = f*f + f0
      IF (ABS(f) .GT. 2) THEN
        map(i,j)=1
        n(i,j)=k
        EXIT
      END IF 
    END DO
  END DO
END DO

CALL MPI_ALLREDUCE(map,maps,hx*hy,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(n,ns,hx*hy,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

IF (rank .EQ. 0) THEN
  OPEN(99,FILE='map.dat',STATUS='replace')
  DO i=1,hx
    x=(i-1)*h+xmin
    DO j=1,hy
      y=(j-1)*h+ymin
      WRITE(99,*) x, y, ns(i,j)
    END DO
  END DO 
END IF

CALL MPI_FINALIZE(ierr)
STOP

END



