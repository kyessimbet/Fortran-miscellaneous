"""
This program simulates the classical dynamics of N gravitating particles using the velocity Verlet algorithm.


"""



  use omp_lib
  implicit none

  ! Parameters
  integer :: k, i, j, a, b, q
  double precision :: ti, tf, dt, e, er, t
  double precision :: f1x, f1y, f1z, f1, r2, r
  double precision, parameter :: eps2 = 1.0d-12   ! softening to avoid singularities

  ! Arrays
  double precision, allocatable :: v1x(:,:), v1y(:,:), v1z(:,:)
  double precision, allocatable :: v2x(:,:), v2y(:,:), v2z(:,:)
  double precision, allocatable :: x(:,:), y(:,:), z(:,:)
  double precision, allocatable :: f2x(:), f2y(:), f2z(:), m(:)

  ! I/O
  integer :: iu_in, iu_out, idx

  ! Simulation setup
  ti = 0.0d0
  tf = 5.0d0
  dt = 0.001d0
  q  = int((tf - ti)/dt) + 1

  ! Read number of particles
  read(*,*) k

  ! Allocate
  allocate(v1x(q-1,k), v1y(q-1,k), v1z(q-1,k))
  allocate(v2x(q,k),   v2y(q,k),   v2z(q,k))
  allocate(x(q,k),     y(q,k),     z(q,k))
  allocate(f2x(k),     f2y(k),     f2z(k), m(k))

  ! Initialize to zero
  v1x = 0.0d0; v1y = 0.0d0; v1z = 0.0d0
  v2x = 0.0d0; v2y = 0.0d0; v2z = 0.0d0
  x   = 0.0d0; y   = 0.0d0; z   = 0.0d0
  f2x = 0.0d0; f2y = 0.0d0; f2z = 0.0d0
  m   = 0.0d0

  ! Read initial data (skip header line)
  iu_in = 100
  open(unit=iu_in, file="init.dat", form="formatted", status="old", action="read")
  read(iu_in, *)
  do i = 1, k
     read(iu_in, *) idx, m(i), x(1,i), y(1,i), z(1,i), v2x(1,i), v2y(1,i), v2z(1,i)
  end do
  close(iu_in)

  ! Initial forces
  do i = 1, k
     f1x = 0.0d0; f1y = 0.0d0; f1z = 0.0d0
     do j = 1, k
        if (i /= j) then
           r2 = (x(1,i)-x(1,j))*(x(1,i)-x(1,j)) + &
                 (y(1,i)-y(1,j))*(y(1,i)-y(1,j)) + &
                 (z(1,i)-z(1,j))*(z(1,i)-z(1,j)) + eps2
           r  = sqrt(r2)
           f1 = (m(i)*m(j))/r2
           f1x = f1x + f1*(x(1,j)-x(1,i))/r
           f1y = f1y + f1*(y(1,j)-y(1,i))/r
           f1z = f1z + f1*(z(1,j)-z(1,i))/r
        end if
     end do
     f2x(i) = f1x
     f2y(i) = f1y
     f2z(i) = f1z
  end do

  ! Initial energy
  er = 0.0d0
  do i = 1, k
     er = er + 0.5d0*m(i)*( v2x(1,i)*v2x(1,i) + v2y(1,i)*v2y(1,i) + v2z(1,i)*v2z(1,i) )
     do a = i+1, k
        r2 = (x(1,a)-x(1,i))*(x(1,a)-x(1,i)) + &
              (y(1,a)-y(1,i))*(y(1,a)-y(1,i)) + &
              (z(1,a)-z(1,i))*(z(1,a)-z(1,i)) + eps2
        er = er - m(i)*m(a)/sqrt(r2)
     end do
  end do

  t = ti

  !$omp parallel default(shared) private(j,i,a,b,f1x,f1y,f1z,r2,r,f1)
  !$omp do schedule(static)
  do j = 1, q-1
     ! Velocity half-step and position update
     do i = 1, k
        v1x(j,i) = v2x(j,i) + 0.5d0*f2x(i)*dt/m(i)
        v1y(j,i) = v2y(j,i) + 0.5d0*f2y(i)*dt/m(i)
        v1z(j,i) = v2z(j,i) + 0.5d0*f2z(i)*dt/m(i)
        x(j+1,i) = x(j,i) + v1x(j,i)*dt
        y(j+1,i) = y(j,i) + v1y(j,i)*dt
        z(j+1,i) = z(j,i) + v1z(j,i)*dt
     end do

     ! New forces at j+1
     do a = 1, k
        f1x = 0.0d0; f1y = 0.0d0; f1z = 0.0d0
        do b = 1, k
           if (a /= b) then
              r2 = (x(j+1,a)-x(j+1,b))*(x(j+1,a)-x(j+1,b)) + &
                    (y(j+1,a)-y(j+1,b))*(y(j+1,a)-y(j+1,b)) + &
                    (z(j+1,a)-z(j+1,b))*(z(j+1,a)-z(j+1,b)) + eps2
              r  = sqrt(r2)
              f1 = (m(a)*m(b))/r2
              f1x = f1x + f1*(x(j+1,b)-x(j+1,a))/r
              f1y = f1y + f1*(y(j+1,b)-y(j+1,a))/r
              f1z = f1z + f1*(z(j+1,b)-z(j+1,a))/r
           end if
        end do
        f2x(a) = f1x
        f2y(a) = f1y
        f2z(a) = f1z
     end do

     ! Velocity full-step
     do i = 1, k
        v2x(j+1,i) = v1x(j,i) + 0.5d0*f2x(i)*dt/m(i)
        v2y(j+1,i) = v1y(j,i) + 0.5d0*f2y(i)*dt/m(i)
        v2z(j+1,i) = v1z(j,i) + 0.5d0*f2z(i)*dt/m(i)
     end do

     ! Energy at step j+1
     e = 0.0d0
     !$omp parallel do reduction(+:e) private(a,r2)
     do i = 1, k
        e = e + 0.5d0*m(i)*( v2x(j+1,i)*v2x(j+1,i) + v2y(j+1,i)*v2y(j+1,i) + v2z(j+1,i)*v2z(j+1,i) )
        do a = i+1, k
           r2 = (x(j+1,a)-x(j+1,i))*(x(j+1,a)-x(j+1,i)) + &
                 (y(j+1,a)-y(j+1,i))*(y(j+1,a)-y(j+1,i)) + &
                 (z(j+1,a)-z(j+1,i))*(z(j+1,a)-z(j+1,i)) + eps2
           e = e - m(i)*m(a)/sqrt(r2)
        end do
     end do

     print *, j, dt*j, e, (e-er)/er
     er = e
     t  = t + dt
  end do
  !$omp end do
  !$omp end parallel

  ! Write final state in init.dat-like format
  iu_out = 200
  open(unit=iu_out, file="final.dat", form="formatted", status="replace", action="write")
  write(iu_out,'(A)') '# final.dat'
  do i = 1, k
     write(iu_out,'(I6,1X,ES12.5,3(1X,ES12.5),3(1X,ES12.5))') i, m(i), x(q,i), y(q,i), z(q,i), &
                                                              v2x(q,i), v2y(q,i), v2z(q,i)
  end do
  close(iu_out)

end program nbody_verlet

