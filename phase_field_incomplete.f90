program phase_field
  implicit none

  ! Parameters
  real(8), parameter :: pi = 4.d0 * atan(1.d0)
  real(8), parameter :: M = 1.0d0
  real(8), parameter :: alpha = -1.0d0
  real(8), parameter :: beta = 1.0d0
  real(8), parameter :: gamma = 1.0d0
  integer, parameter :: Lx = 256, Ly = 256
  real(8), parameter :: deltat = 0.005d0
  integer, parameter :: ts = 10000, output_frq = 100
  real(8), parameter :: phi_ini = 0.0d0

  ! Variables
  real(8) :: phi(0:Lx-1, 0:Ly-1, 2)
  real(8) :: mu(0:Lx-1, 0:Ly-1)
  real(8) :: lap_phi, lap_mu
  real(8) :: ram_1, ram_2
  integer :: ix, iy, t, i_output
  integer :: ixp, ixm, iyp, iym
  integer :: layerp, layern
  integer :: seedsize
  integer, allocatable :: seed(:)
  character(90) :: term

  layerp = 1
  layern = 2

  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do ix = 1, seedsize
    call system_clock(count=seed(ix))
  end do
  call random_seed(put=seed(:))

  ! Initialize phi
  do iy = 0, Ly - 1
    do ix = 0, Lx - 1
      call random_number(ram_1)
      call random_number(ram_2)
      phi(ix, iy, layerp) = sqrt(-2.d0 * log(ram_1)) * cos(2.d0 * pi * ram_2) + phi_ini
    end do
  end do

  ! Output initial state
  open(21, file='phi_0000.d')
  do ix = 0, Lx - 1
    do iy = 0, Ly - 1
      write(21, '(2i4, f12.5)') ix, iy, phi(ix, iy, layerp)
    end do
  end do
  close(21)

  ! Time evolution
  i_output = 0
  do t = 1, ts
    if (mod(t, output_frq) == 0) write(*,*) t, '/', ts

    ! Calculate mu
    do ix = 0, Lx - 1
      do iy = 0, Ly - 1
        ixp = mod(ix + 1, Lx)
        ixm = mod(ix - 1 + Lx, Lx)
        iyp = mod(iy + 1, Ly)
        iym = mod(iy - 1 + Ly, Ly)

        lap_phi = phi(ixp, iy, layerp) + phi(ixm, iy, layerp) + &
                  phi(ix, iyp, layerp) + phi(ix, iym, layerp) - &
                  4.d0 * phi(ix, iy, layerp)

        mu(ix, iy) = alpha * phi(ix, iy, layerp) + &
                     beta * phi(ix, iy, layerp)**3 - &
                     gamma * lap_phi
      end do
    end do

    ! Update phi
    do ix = 0, Lx - 1
      do iy = 0, Ly - 1
        ixp = mod(ix + 1, Lx)
        ixm = mod(ix - 1 + Lx, Lx)
        iyp = mod(iy + 1, Ly)
        iym = mod(iy - 1 + Ly, Ly)

        lap_mu = mu(ixp, iy) + mu(ixm, iy) + mu(ix, iyp) + mu(ix, iym) - &
                 4.d0 * mu(ix, iy)

        phi(ix, iy, layern) = phi(ix, iy, layerp) + deltat * M * lap_mu
      end do
    end do

    ! Output results
    if (mod(t, output_frq) == 0) then
      i_output = i_output + 1
      write(term, '("phi_", i4.4, ".d")') i_output
      open(21, file=term)
      do ix = 0, Lx - 1
        do iy = 0, Ly - 1
          write(21, '(2i4, f12.5)') ix, iy, phi(ix, iy, layern)
        end do
      end do
      close(21)
    end if

    ! Swap layers
    layerp = 3 - layerp
    layern = 3 - layern
  end do
end program phase_field
