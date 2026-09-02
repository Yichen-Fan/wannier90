program settings_growth
  use w90_library
  use mpi
  implicit none

  type(lib_common_type) :: w90
  integer :: ierr, stdout, stderr
  integer :: mp_grid(3) = [1, 1, 1]
  integer :: distk(1) = [0]
  integer :: num_wann = 2
  real(kind=8) :: unit_cell(3, 3), kpoints(3, 1), atoms_frac(3, 2)
  character(len=1) :: symbols(2) = ['H', 'H']

  call w90_get_fortran_stdout(stdout)
  call w90_get_fortran_stderr(stderr)

  call MPI_Init(ierr)

  unit_cell = reshape([3.0d0, 0.0d0, 0.0d0, &
                       0.0d0, 3.0d0, 0.0d0, &
                       0.0d0, 0.0d0, 3.0d0], [3, 3])
  kpoints(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  atoms_frac(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  atoms_frac(:, 2) = [0.5d0, 0.0d0, 0.0d0]

  ! Match the Python C-interface call order and reach the initial settings
  ! capacity exactly. Growing the settings array must preserve the repeated
  ! one-letter symbols so that H:s expands to one projection on each H atom.
  call w90_set_comm(w90, MPI_COMM_WORLD)
  call w90_set_option(w90, 'kpoints', kpoints)
  call w90_set_option(w90, 'mp_grid', mp_grid)
  call w90_set_option(w90, 'distk', distk)
  call w90_set_option(w90, 'total_bands', num_wann)
  call w90_set_option(w90, 'num_wann', num_wann)
  call w90_set_option(w90, 'unit_cell_cart', unit_cell)
  call w90_set_option(w90, 'symbols', symbols)
  call w90_set_option(w90, 'atoms_frac', atoms_frac)
  call w90_set_option(w90, 'gamma_only', .false.)
  call w90_set_option(w90, 'spinors', .false.)
  call w90_set_option(w90, 'use_bloch_phases', .false.)
  call w90_set_option(w90, 'projections', 'H:s')
  call w90_set_option(w90, 'iprint', 0)
  call w90_set_option(w90, 'num_iter', 10)
  call w90_set_option(w90, 'conv_window', 3)
  call w90_set_option(w90, 'conv_tol', 1.0d-8)
  call w90_set_option(w90, 'num_cg_steps', 5)
  call w90_set_option(w90, 'num_dump_cycles', 100)
  call w90_set_option(w90, 'num_print_cycles', 1)
  call w90_set_option(w90, 'dis_num_iter', 200)

  call w90_input_setopt(w90, 'wannier90-growth', stdout, stderr, ierr)

  if (ierr /= 0) then
    write (*, '(a,i0)') 'FAIL: settings-growth w90_input_setopt returned ierr = ', ierr
    call MPI_Finalize(ierr)
    stop 1
  end if

  call MPI_Finalize(ierr)

end program settings_growth
