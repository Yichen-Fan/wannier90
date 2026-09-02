!-*- mode: F90 -*-!
!------------------------------------------------------------!
! Copyright (C) 2026 Wannier Developer Group                 !
!                                                            !
! This library is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU Lesser General !
! Public License as published by the Free Software           !
! Foundation; either version 2.1 of the License, or (at your  !
! option) any later version.                                 !
!                                                            !
! This library is distributed in the hope that it will be    !
! useful,but WITHOUT ANY WARRANTY; without even the implied  !
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR    !
! PURPOSE.  See the GNU Lesser General Public License for    !
! more details.                                              !
!                                                            !
! You should have received a copy of the GNU Lesser General  !
! Public License along with this library; if not, see        !
! <https://www.gnu.org/licenses/>.                           !
!------------------------------------------------------------!

program lbfgs_unit

  use, intrinsic :: iso_fortran_env, only: error_unit
#ifdef MPI
#ifndef MPIH
  use mpi, only: MPI_COMM_WORLD, MPI_INTEGER, MPI_SUM, MPI_SUCCESS, &
                 MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, MPI_Finalize, MPI_Init
#endif
#endif
  use w90_comms, only: w90_comm_type
  use w90_constants, only: dp
  use w90_error, only: w90_error_type
  use w90_lbfgs, only: lbfgs_apply, lbfgs_history_capacity, lbfgs_history_count, &
                       lbfgs_initialize, lbfgs_is_initialized, lbfgs_push, lbfgs_reset, &
                       lbfgs_state_type

  implicit none

#ifdef MPIH
  include 'mpif.h'
#endif

  type(w90_comm_type) :: comm
  integer :: failures, rank, total_failures
#ifdef MPI
  integer :: ierr, num_ranks
#endif

  failures = 0
#ifdef MPI
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Init failed in L-BFGS unit test'
  comm%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Comm_rank failed in L-BFGS unit test'
  call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Comm_size failed in L-BFGS unit test'
#else
  rank = 0
#endif

  call test_identity_and_one_pair(comm, rank, failures)
  call test_damping_and_skipped_pair(comm, rank, failures)
  call test_circular_history_and_reset(comm, rank, failures)
#ifdef MPI
  call test_distributed_products(comm, rank, num_ranks, failures)
  call test_zero_local_vectors(comm, rank, num_ranks, failures)

  call MPI_Allreduce(failures, total_failures, 1, MPI_INTEGER, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Allreduce failed in L-BFGS unit test'
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Finalize failed in L-BFGS unit test'
#else
  total_failures = failures
#endif

  if (total_failures /= 0) then
    if (rank == 0) write (error_unit, '(a,i0)') 'L-BFGS unit-test failures: ', total_failures
    error stop 1
  end if
  if (rank == 0) write (*, '(a)') 'L-BFGS unit tests passed'

contains

  subroutine test_identity_and_one_pair(comm, rank, failures)
    type(w90_comm_type), intent(in) :: comm
    integer, intent(in) :: rank
    integer, intent(inout) :: failures

    type(lbfgs_state_type) :: state
    type(w90_error_type), allocatable :: error
    complex(kind=dp) :: actual(1, 1, 2), expected(1, 1, 2)
    complex(kind=dp) :: step(1, 1, 2), vector(1, 1, 2), y(1, 1, 2)
    logical :: accepted, damped

    call lbfgs_initialize(state, 1, 1, 2, 3, error, comm)
    if (.not. expect_no_error(error, 'initialize identity history', rank, failures)) return
    call expect_true(lbfgs_is_initialized(state), 'history is initialized', rank, failures)
    call expect_integer(lbfgs_history_capacity(state), 3, 'history capacity', rank, failures)
    call expect_integer(lbfgs_history_count(state), 0, 'empty history count', rank, failures)

    vector(1, 1, 1) = cmplx(3.0_dp, 2.0_dp, kind=dp)
    vector(1, 1, 2) = cmplx(-4.0_dp, 0.5_dp, kind=dp)
    call lbfgs_apply(state, vector, actual, error, comm)
    if (.not. expect_no_error(error, 'apply empty history', rank, failures)) return
    call expect_close(actual, vector, 'empty history is identity', rank, failures)

    step = cmplx(0.0_dp, 0.0_dp, kind=dp)
    y = cmplx(0.0_dp, 0.0_dp, kind=dp)
    step(1, 1, 1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
    y(1, 1, 1) = cmplx(2.0_dp, 0.0_dp, kind=dp)
    call lbfgs_push(state, step, y, 0.1_dp, accepted, damped, error, comm)
    if (.not. expect_no_error(error, 'push exact pair', rank, failures)) return
    call expect_true(accepted, 'well-conditioned pair is accepted', rank, failures)
    call expect_false(damped, 'well-conditioned pair is not damped', rank, failures)
    call expect_integer(lbfgs_history_count(state), 1, 'count after exact pair', rank, failures)

    vector(1, 1, 1) = cmplx(3.0_dp, 0.0_dp, kind=dp)
    vector(1, 1, 2) = cmplx(4.0_dp, 0.0_dp, kind=dp)
    expected(1, 1, 1) = cmplx(1.5_dp, 0.0_dp, kind=dp)
    expected(1, 1, 2) = cmplx(4.0_dp, 0.0_dp, kind=dp)
    call lbfgs_apply(state, vector, actual, error, comm)
    if (.not. expect_no_error(error, 'apply exact pair', rank, failures)) return
    call expect_close(actual, expected, 'one-pair two-loop product', rank, failures)
  end subroutine test_identity_and_one_pair

  subroutine test_damping_and_skipped_pair(comm, rank, failures)
    type(w90_comm_type), intent(in) :: comm
    integer, intent(in) :: rank
    integer, intent(inout) :: failures

    type(lbfgs_state_type) :: state
    type(w90_error_type), allocatable :: error
    complex(kind=dp) :: actual(1, 1, 1), expected(1, 1, 1)
    complex(kind=dp) :: step(1, 1, 1), vector(1, 1, 1), y(1, 1, 1)
    logical :: accepted, damped

    call lbfgs_initialize(state, 1, 1, 1, 2, error, comm)
    if (.not. expect_no_error(error, 'initialize damping history', rank, failures)) return

    step = cmplx(1.0_dp, 0.0_dp, kind=dp)
    y = cmplx(-1.0_dp, 0.0_dp, kind=dp)
    call lbfgs_push(state, step, y, 0.2_dp, accepted, damped, error, comm)
    if (.not. expect_no_error(error, 'push damped pair', rank, failures)) return
    call expect_true(accepted, 'negative-curvature pair is damped and accepted', rank, failures)
    call expect_true(damped, 'negative-curvature pair reports damping', rank, failures)

    vector = cmplx(1.0_dp, 0.0_dp, kind=dp)
    expected = cmplx(5.0_dp, 0.0_dp, kind=dp)
    call lbfgs_apply(state, vector, actual, error, comm)
    if (.not. expect_no_error(error, 'apply damped pair', rank, failures)) return
    call expect_close(actual, expected, 'damped pair reaches curvature floor', rank, failures)

    step = cmplx(0.0_dp, 0.0_dp, kind=dp)
    y = cmplx(1.0_dp, 0.0_dp, kind=dp)
    call lbfgs_push(state, step, y, 0.2_dp, accepted, damped, error, comm)
    if (.not. expect_no_error(error, 'push zero step', rank, failures)) return
    call expect_false(accepted, 'zero step is skipped', rank, failures)
    call expect_false(damped, 'skipped pair is not reported as damped', rank, failures)
    call expect_integer(lbfgs_history_count(state), 1, 'skipped pair preserves count', rank, failures)
  end subroutine test_damping_and_skipped_pair

  subroutine test_circular_history_and_reset(comm, rank, failures)
    type(w90_comm_type), intent(in) :: comm
    integer, intent(in) :: rank
    integer, intent(inout) :: failures

    type(lbfgs_state_type) :: state
    type(w90_error_type), allocatable :: error
    complex(kind=dp) :: actual(1, 1, 3), expected(1, 1, 3)
    complex(kind=dp) :: step(1, 1, 3), vector(1, 1, 3), y(1, 1, 3)
    integer :: pair_index
    logical :: accepted, damped

    call lbfgs_initialize(state, 1, 1, 3, 2, error, comm)
    if (.not. expect_no_error(error, 'initialize circular history', rank, failures)) return

    do pair_index = 1, 3
      step = cmplx(0.0_dp, 0.0_dp, kind=dp)
      y = cmplx(0.0_dp, 0.0_dp, kind=dp)
      step(1, 1, pair_index) = cmplx(1.0_dp, 0.0_dp, kind=dp)
      y(1, 1, pair_index) = cmplx(2.0_dp**pair_index, 0.0_dp, kind=dp)
      call lbfgs_push(state, step, y, 0.1_dp, accepted, damped, error, comm)
      if (.not. expect_no_error(error, 'push circular pair', rank, failures)) return
      call expect_true(accepted, 'circular pair is accepted', rank, failures)
      call expect_false(damped, 'circular pair is not damped', rank, failures)
    end do
    call expect_integer(lbfgs_history_count(state), 2, 'circular history is capped', rank, failures)

    vector(1, 1, 1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
    vector(1, 1, 2) = cmplx(4.0_dp, 0.0_dp, kind=dp)
    vector(1, 1, 3) = cmplx(8.0_dp, 0.0_dp, kind=dp)
    expected = cmplx(1.0_dp, 0.0_dp, kind=dp)
    call lbfgs_apply(state, vector, actual, error, comm)
    if (.not. expect_no_error(error, 'apply circular history', rank, failures)) return
    call expect_close(actual, expected, 'circular history retains newest pairs', rank, failures)

    call lbfgs_reset(state)
    call expect_integer(lbfgs_history_count(state), 0, 'reset empties history', rank, failures)
    call expect_integer(lbfgs_history_capacity(state), 2, 'reset preserves capacity', rank, failures)
    call lbfgs_apply(state, vector, actual, error, comm)
    if (.not. expect_no_error(error, 'apply reset history', rank, failures)) return
    call expect_close(actual, vector, 'reset restores identity product', rank, failures)
  end subroutine test_circular_history_and_reset

#ifdef MPI
  subroutine test_distributed_products(comm, rank, num_ranks, failures)
    type(w90_comm_type), intent(in) :: comm
    integer, intent(in) :: rank, num_ranks
    integer, intent(inout) :: failures

    type(lbfgs_state_type) :: state
    type(w90_error_type), allocatable :: error
    complex(kind=dp) :: actual(1, 1, 1), expected(1, 1, 1)
    complex(kind=dp) :: step(1, 1, 1), vector(1, 1, 1), y(1, 1, 1)
    logical :: accepted, damped

    if (num_ranks /= 2) then
      if (rank == 0) write (*, '(a)') 'Skipping distributed-product check: it requires two MPI ranks'
      return
    end if

    call lbfgs_initialize(state, 1, 1, 1, 2, error, comm)
    if (.not. expect_no_error(error, 'initialize distributed history', rank, failures)) return
    step = cmplx(1.0_dp, 0.0_dp, kind=dp)
    if (rank == 0) then
      y = cmplx(2.0_dp, 0.0_dp, kind=dp)
      vector = cmplx(3.0_dp, 0.0_dp, kind=dp)
      expected = cmplx(16.0_dp/9.0_dp, 0.0_dp, kind=dp)
    else
      y = cmplx(4.0_dp, 0.0_dp, kind=dp)
      vector = cmplx(5.0_dp, 0.0_dp, kind=dp)
      expected = cmplx(10.0_dp/9.0_dp, 0.0_dp, kind=dp)
    end if
    call lbfgs_push(state, step, y, 0.1_dp, accepted, damped, error, comm)
    if (.not. expect_no_error(error, 'push distributed pair', rank, failures)) return
    call expect_true(accepted, 'distributed pair is accepted', rank, failures)
    call expect_false(damped, 'distributed pair is not damped', rank, failures)
    call lbfgs_apply(state, vector, actual, error, comm)
    if (.not. expect_no_error(error, 'apply distributed pair', rank, failures)) return
    call expect_close(actual, expected, 'two-loop recursion uses global products', rank, failures)
  end subroutine test_distributed_products

  subroutine test_zero_local_vectors(comm, rank, num_ranks, failures)
    type(w90_comm_type), intent(in) :: comm
    integer, intent(in) :: rank, num_ranks
    integer, intent(inout) :: failures

    type(lbfgs_state_type) :: state
    type(w90_error_type), allocatable :: error
    complex(kind=dp), allocatable :: actual(:, :, :), expected(:, :, :)
    complex(kind=dp), allocatable :: step(:, :, :), vector(:, :, :), y(:, :, :)
    integer :: num_vectors_local
    logical :: accepted, damped

    if (num_ranks /= 2) return
    if (rank == 0) then
      num_vectors_local = 1
    else
      num_vectors_local = 0
    end if
    allocate (actual(1, 1, num_vectors_local), expected(1, 1, num_vectors_local), &
              step(1, 1, num_vectors_local), vector(1, 1, num_vectors_local), &
              y(1, 1, num_vectors_local))

    call lbfgs_initialize(state, 1, 1, num_vectors_local, 2, error, comm)
    if (.not. expect_no_error(error, 'initialize zero-local history', rank, failures)) return
    if (rank == 0) then
      step = cmplx(2.0_dp, 0.0_dp, kind=dp)
      y = cmplx(8.0_dp, 0.0_dp, kind=dp)
      vector = cmplx(6.0_dp, 0.0_dp, kind=dp)
      expected = cmplx(1.5_dp, 0.0_dp, kind=dp)
    end if
    call lbfgs_push(state, step, y, 0.1_dp, accepted, damped, error, comm)
    if (.not. expect_no_error(error, 'push zero-local pair', rank, failures)) return
    call expect_true(accepted, 'zero-local rank shares accepted update', rank, failures)
    call expect_false(damped, 'zero-local update is not damped', rank, failures)
    call expect_integer(lbfgs_history_count(state), 1, 'zero-local rank advances history', rank, failures)
    call lbfgs_apply(state, vector, actual, error, comm)
    if (.not. expect_no_error(error, 'apply zero-local pair', rank, failures)) return
    call expect_close(actual, expected, 'zero-local rank participates in recursion', rank, failures)
  end subroutine test_zero_local_vectors
#endif

  logical function expect_no_error(error, label, rank, failures)
    type(w90_error_type), allocatable, intent(inout) :: error
    character(len=*), intent(in) :: label
    integer, intent(in) :: rank
    integer, intent(inout) :: failures

    expect_no_error = .not. allocated(error)
    if (.not. expect_no_error) then
      failures = failures + 1
      write (error_unit, '(a,i0,3a)') 'rank ', rank, ': ', trim(label), ': '//trim(error%message)
      deallocate (error)
    end if
  end function expect_no_error

  subroutine expect_close(actual, expected, label, rank, failures)
    complex(kind=dp), intent(in) :: actual(:, :, :), expected(:, :, :)
    character(len=*), intent(in) :: label
    integer, intent(in) :: rank
    integer, intent(inout) :: failures

    real(kind=dp), parameter :: tolerance = 256.0_dp*epsilon(1.0_dp)
    real(kind=dp) :: error_norm, scale

    if (any(shape(actual) /= shape(expected))) then
      failures = failures + 1
      write (error_unit, '(a,i0,3a)') 'rank ', rank, ': ', trim(label), ': shape mismatch'
      return
    end if
    if (size(actual) == 0) return
    error_norm = maxval(abs(actual - expected))
    scale = max(1.0_dp, maxval(abs(expected)))
    if (error_norm > tolerance*scale) then
      failures = failures + 1
      write (error_unit, '(a,i0,3a,es12.4)') &
        'rank ', rank, ': ', trim(label), ': maximum scaled error = ', error_norm/scale
    end if
  end subroutine expect_close

  subroutine expect_integer(actual, expected, label, rank, failures)
    integer, intent(in) :: actual, expected, rank
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    if (actual /= expected) then
      failures = failures + 1
      write (error_unit, '(a,i0,3a,i0,a,i0)') &
        'rank ', rank, ': ', trim(label), ': got ', actual, ', expected ', expected
    end if
  end subroutine expect_integer

  subroutine expect_true(condition, label, rank, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(in) :: rank
    integer, intent(inout) :: failures

    if (.not. condition) then
      failures = failures + 1
      write (error_unit, '(a,i0,3a)') 'rank ', rank, ': ', trim(label), ': expected true'
    end if
  end subroutine expect_true

  subroutine expect_false(condition, label, rank, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(in) :: rank
    integer, intent(inout) :: failures

    if (condition) then
      failures = failures + 1
      write (error_unit, '(a,i0,3a)') 'rank ', rank, ': ', trim(label), ': expected false'
    end if
  end subroutine expect_false

end program lbfgs_unit
