program phase_field_case
  implicit none

  ! === パラメータ選択用 ===
  integer, parameter :: case_id = 3   ! 1: 混合, 2: スピノーダル, 3: 核形成
  real(8), parameter :: pi = 4.d0 * atan(1.d0) 
  real(8) :: alpha, phi_ini
  character(20) :: prefix

  ! === 共通パラメータ ===
  real(8), parameter :: beta = 1.0d0, gamma = 1.0d0
  real(8), parameter :: M = 1.0d0
  integer, parameter :: Lx = 256, Ly = 256
  real(8), parameter :: deltat = 0.005d0
  integer, parameter :: ts = 10000, output_frq = 100

  ! === 初期化変数・格納配列 ===
  real(8) :: phi(0:Lx-1, 0:Ly-1, 2), mu(0:Lx-1, 0:Ly-1)
  integer :: ix, iy, ixp, ixm, iyp, iym
  integer :: layerp, layern, t, i_output
  real(8) :: ram_1, ram_2, lap_phi, lap_mu
  character(90) :: term
  integer :: seedsize
  integer, allocatable :: seed(:)

  ! === ケースごとの設定 ===
  select case (case_id)
  case (1) ! 混合状態
    alpha = 1.0d0
    phi_ini = 0.0d0
    prefix = 'case1_phi_'
  case (2) ! スピノーダル分解
    alpha = -1.0d0
    phi_ini = 0.0d0
    prefix = 'case2_phi_'
  case (3) ! 核形成・成長
    alpha = -1.0d0
    phi_ini = 0.6d0
    prefix = 'case3_phi_'
  end select

  ! === ランダム初期化 ===
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do ix = 1, seedsize
    call system_clock(count=seed(ix))
  end do
  call random_seed(put=seed(:))

  layerp = 1
  layern = 2

  ! === 初期濃度設定（平均 phi_ini, 分散1.0）===
  do iy = 0, Ly - 1
    do ix = 0, Lx - 1
      call random_number(ram_1)
      call random_number(ram_2)
      phi(ix, iy, layerp) = sqrt(-2.d0 * log(ram_1)) * cos(2.d0 * pi * ram_2) + phi_ini
    end do
  end do

  ! === 初期状態を出力 ===
  write(term, '(a, a)') prefix, '0000.d'
  open(21, file=term)
  do ix = 0, Lx - 1
    do iy = 0, Ly - 1
      write(21, '(2i4, f12.5)') ix, iy, phi(ix, iy, layerp)
    end do
  end do
  close(21)

  ! === 時間発展 ===
  i_output = 0
  do t = 1, ts
    if (mod(t, output_frq) == 0) write(*,*) t, '/', ts

    ! μの計算
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
                     beta  * phi(ix, iy, layerp)**3 - &
                     gamma * lap_phi
      end do
    end do

    ! φの更新
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

    ! 出力
    if (mod(t, output_frq) == 0) then
      i_output = i_output + 1
      write(term, '(a, i4.4, a)') prefix, i_output, '.d'
      open(21, file=term)
      do ix = 0, Lx - 1
        do iy = 0, Ly - 1
          write(21, '(2i4, f12.5)') ix, iy, phi(ix, iy, layern)
        end do
      end do
      close(21)
    end if

    ! レイヤー入れ替え
    layerp = 3 - layerp
    layern = 3 - layern
  end do

end program phase_field_case
