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
!                                                            !
! The webpage of the Wannier90 code is                       !
! <https://www.wannier.org>.                                 !
!                                                            !
! The Wannier90 code is hosted on GitHub                     !
! <https://github.com/wannier-developers/wannier90>          !
!------------------------------------------------------------!
!                                                            !
!  w90_lbfgs: limited-memory inverse-Hessian products        !
!                                                            !
!------------------------------------------------------------!

module w90_lbfgs
  !! Distributed L-BFGS history for dense complex tangent vectors.
  !!
  !! Vectors are rank-three arrays, normally anti-Hermitian matrices at
  !! the k-points owned by the calling MPI rank. Inner products are real
  !! Frobenius products, summed over all ranks. The module operates on
  !! conventional gradients: a secant pair is
  !!
  !!   s_k = x_(k+1) - x_k,  y_k = g_(k+1) - g_k,
  !!
  !! expressed in a common tangent-space coordinate system. Transport,
  !! manifold updates, line searches, and the sign used for a search
  !! direction belong to the caller. In particular, [[lbfgs_apply]] returns
  !! H_k v; a minimizer normally forms its direction as -H_k g_k.

  use w90_constants, only: dp
  use w90_comms, only: comms_allreduce, w90_comm_type
  use w90_error, only: set_error_alloc, set_error_fatal, w90_error_type
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private

  type, public :: lbfgs_state_type
    !! Circular limited-memory store of distributed secant pairs.
    !!
    !! The array storage is local to each MPI rank. All ranks must call the
    !! public collective routines in the same order and must use identical
    !! history settings. The local third dimension may differ between ranks.
    private
    logical :: initialized = .false.
    integer :: history_size = 0
    integer :: num_stored = 0
    integer :: next_slot = 1
    complex(kind=dp), allocatable :: step(:, :, :, :)
    complex(kind=dp), allocatable :: gradient_change(:, :, :, :)
    real(kind=dp), allocatable :: inverse_curvature(:)
  end type lbfgs_state_type

  public :: lbfgs_apply
  public :: lbfgs_history_capacity
  public :: lbfgs_history_count
  public :: lbfgs_initialize
  public :: lbfgs_is_initialized
  public :: lbfgs_push
  public :: lbfgs_reset

contains

  !===========================================================
  subroutine lbfgs_initialize(state, matrix_size_1, matrix_size_2, num_vectors_local, &
                              history_size, error, comm)
    !===========================================================
    !! Allocate an empty L-BFGS history for a distributed tangent vector.
    !!
    !! The first two dimensions describe a square matrix and the third is the
    !! number of vectors owned by this rank. A zero local third dimension is
    !! valid. `history_size` must be positive. Reinitializing a state discards
    !! its previous allocation and all stored secant pairs.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(out) :: state
    integer, intent(in) :: matrix_size_1
    integer, intent(in) :: matrix_size_2
    integer, intent(in) :: num_vectors_local
    integer, intent(in) :: history_size
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    real(kind=dp) :: allocation_failed
    real(kind=dp) :: distributed_setting(6)
    real(kind=dp) :: invalid_argument(4)
    integer :: ierr

    ! Synchronize validation before any rank returns. The local vector count
    ! may legitimately differ, so argument failures are not assumed uniform.
    invalid_argument = 0.0_dp
    if (matrix_size_1 /= matrix_size_2) invalid_argument(1) = 1.0_dp
    if (matrix_size_1 < 1 .or. matrix_size_2 < 1) invalid_argument(2) = 1.0_dp
    if (num_vectors_local < 0) invalid_argument(3) = 1.0_dp
    if (history_size < 1) invalid_argument(4) = 1.0_dp
    call comms_allreduce(invalid_argument(1), size(invalid_argument), 'MAX', error, comm)
    if (allocated(error)) return

    if (invalid_argument(1) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS tangent matrices must be square', comm)
      return
    end if
    if (invalid_argument(2) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS tangent matrices must have positive size', comm)
      return
    end if
    if (invalid_argument(3) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS local vector count must be non-negative', comm)
      return
    end if
    if (invalid_argument(4) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS history size must be positive', comm)
      return
    end if

    ! The local vector count may differ; matrix and history dimensions may not.
    distributed_setting = [real(matrix_size_1, dp), -real(matrix_size_1, dp), &
                           real(matrix_size_2, dp), -real(matrix_size_2, dp), &
                           real(history_size, dp), -real(history_size, dp)]
    call comms_allreduce(distributed_setting(1), size(distributed_setting), 'MAX', error, comm)
    if (allocated(error)) return
    if (distributed_setting(1) /= -distributed_setting(2) .or. &
        distributed_setting(3) /= -distributed_setting(4) .or. &
        distributed_setting(5) /= -distributed_setting(6)) then
      call set_error_fatal(error, 'L-BFGS dimensions are inconsistent across ranks', comm)
      return
    end if

    allocate (state%step(matrix_size_1, matrix_size_2, num_vectors_local, history_size), &
              state%gradient_change(matrix_size_1, matrix_size_2, num_vectors_local, history_size), &
              state%inverse_curvature(history_size), stat=ierr)
    allocation_failed = 0.0_dp
    if (ierr /= 0) allocation_failed = 1.0_dp
    call comms_allreduce(allocation_failed, 1, 'MAX', error, comm)
    if (allocated(error)) return
    if (allocation_failed > 0.5_dp) then
      call set_error_alloc(error, 'Error allocating L-BFGS history', comm)
      return
    end if

    state%history_size = history_size
    state%initialized = .true.
    call lbfgs_reset(state)

  end subroutine lbfgs_initialize

  !===========================================================
  subroutine lbfgs_reset(state)
    !===========================================================
    !! Empty a history without releasing its allocated storage.
    !!
    !! Calling this routine on an uninitialized history is harmless.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(inout) :: state

    state%num_stored = 0
    state%next_slot = 1
    if (allocated(state%inverse_curvature)) state%inverse_curvature = 0.0_dp

  end subroutine lbfgs_reset

  !===========================================================
  subroutine lbfgs_push(state, step, gradient_change, curvature_tolerance, accepted, damped, &
                        error, comm)
    !===========================================================
    !! Validate and append a distributed secant pair.
    !!
    !! `step` is s_k and `gradient_change` is y_k in a common body-coordinate
    !! tangent space. A valid pair has positive curvature. A pair below the
    !! relative curvature threshold is shifted along s to reach the target.
    !! The tolerance is a cosine-like relative threshold and must lie strictly
    !! between zero and one.
    !! A pair that is zero, non-finite, or still fails the roundoff-level
    !! positive-curvature check is skipped. A skipped pair does not evict an
    !! existing entry from a full circular history.
    !!
    !! `accepted` reports whether a pair was stored and `damped` reports
    !! whether the stored gradient change was shifted. Both are false for a
    !! nonfatal skipped update and when `error` is set.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(inout) :: state
    complex(kind=dp), intent(in) :: step(:, :, :)
    complex(kind=dp), intent(in) :: gradient_change(:, :, :)
    real(kind=dp), intent(in) :: curvature_tolerance
    logical, intent(out) :: accepted
    logical, intent(out) :: damped
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    complex(kind=dp), allocatable :: stored_gradient_change(:, :, :)
    real(kind=dp) :: allocation_failed
    real(kind=dp) :: curvature, curvature_floor, damping_shift
    real(kind=dp) :: inverse_curvature
    real(kind=dp) :: step_norm, step_norm2
    real(kind=dp) :: change_norm, change_norm2
    real(kind=dp) :: distributed_setting(8)
    real(kind=dp) :: invalid_argument(3)
    real(kind=dp) :: products(3)
    real(kind=dp) :: stored_products(2)
    integer :: ierr
    integer :: slot
    logical :: pair_damped
    logical :: storage_is_valid

    accepted = .false.
    damped = .false.
    pair_damped = .false.

    ! Shape and initialization failures can be rank-local. Synchronize all
    ! validation results before returning so later collectives remain aligned.
    invalid_argument = 0.0_dp
    storage_is_valid = lbfgs_storage_is_valid(state)
    if (.not. storage_is_valid) invalid_argument(1) = 1.0_dp
    if (storage_is_valid) then
      if (.not. lbfgs_shape_matches(state, step) .or. &
          .not. lbfgs_shape_matches(state, gradient_change)) invalid_argument(2) = 1.0_dp
    end if
    if (.not. ieee_is_finite(curvature_tolerance) .or. &
        curvature_tolerance <= 0.0_dp .or. curvature_tolerance >= 1.0_dp) then
      invalid_argument(3) = 1.0_dp
    end if
    call comms_allreduce(invalid_argument(1), size(invalid_argument), 'MAX', error, comm)
    if (allocated(error)) return

    if (invalid_argument(1) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS history is not initialized', comm)
      return
    end if
    if (invalid_argument(2) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS secant vector has incompatible shape', comm)
      return
    end if
    if (invalid_argument(3) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS curvature tolerance must lie strictly between zero and one', comm)
      return
    end if

    distributed_setting = [real(state%history_size, dp), -real(state%history_size, dp), &
                           real(state%num_stored, dp), -real(state%num_stored, dp), &
                           real(state%next_slot, dp), -real(state%next_slot, dp), &
                           curvature_tolerance, -curvature_tolerance]
    call comms_allreduce(distributed_setting(1), size(distributed_setting), 'MAX', error, comm)
    if (allocated(error)) return
    if (distributed_setting(1) /= -distributed_setting(2) .or. &
        distributed_setting(3) /= -distributed_setting(4) .or. &
        distributed_setting(5) /= -distributed_setting(6) .or. &
        distributed_setting(7) /= -distributed_setting(8)) then
      call set_error_fatal(error, 'L-BFGS state or curvature tolerance is inconsistent across ranks', comm)
      return
    end if

    products(1) = lbfgs_local_inner_product(step, step)
    products(2) = lbfgs_local_inner_product(gradient_change, gradient_change)
    products(3) = lbfgs_local_inner_product(step, gradient_change)
    call comms_allreduce(products(1), 3, 'SUM', error, comm)
    if (allocated(error)) return

    step_norm2 = products(1)
    change_norm2 = products(2)
    curvature = products(3)
    if (.not. ieee_is_finite(step_norm2) .or. .not. ieee_is_finite(change_norm2) .or. &
        .not. ieee_is_finite(curvature) .or. step_norm2 <= tiny(1.0_dp) .or. &
        change_norm2 <= tiny(1.0_dp)) then
      return
    end if

    step_norm = sqrt(max(step_norm2, 0.0_dp))
    change_norm = sqrt(max(change_norm2, 0.0_dp))
    curvature_floor = curvature_tolerance*step_norm*change_norm
    damping_shift = 0.0_dp

    if (curvature < curvature_floor) then
      damping_shift = (curvature_floor - curvature)/step_norm2
      if (.not. ieee_is_finite(damping_shift)) return

      ! Stage the exact vector that would be stored. In particular, do not
      ! overwrite the oldest circular slot until the recomputed pair passes.
      allocate (stored_gradient_change(size(gradient_change, 1), size(gradient_change, 2), &
                                       size(gradient_change, 3)), stat=ierr)
      allocation_failed = 0.0_dp
      if (ierr /= 0) allocation_failed = 1.0_dp
      call comms_allreduce(allocation_failed, 1, 'MAX', error, comm)
      if (allocated(error)) return
      if (allocation_failed > 0.5_dp) then
        call set_error_alloc(error, 'Error allocating damped L-BFGS secant vector', comm)
        return
      end if

      stored_gradient_change = gradient_change + damping_shift*step
      stored_products(1) = lbfgs_local_inner_product(stored_gradient_change, stored_gradient_change)
      stored_products(2) = lbfgs_local_inner_product(step, stored_gradient_change)
      call comms_allreduce(stored_products(1), size(stored_products), 'SUM', error, comm)
      if (allocated(error)) return

      change_norm2 = stored_products(1)
      curvature = stored_products(2)
      if (.not. ieee_is_finite(change_norm2) .or. .not. ieee_is_finite(curvature) .or. &
          change_norm2 <= tiny(1.0_dp)) return
      change_norm = sqrt(max(change_norm2, 0.0_dp))
      pair_damped = .true.
    end if

    if (.not. ieee_is_finite(curvature) .or. &
        curvature <= epsilon(1.0_dp)*step_norm*max(change_norm, tiny(1.0_dp))) then
      return
    end if

    inverse_curvature = 1.0_dp/curvature
    if (.not. ieee_is_finite(inverse_curvature)) return

    slot = state%next_slot
    state%step(:, :, :, slot) = step
    if (pair_damped) then
      state%gradient_change(:, :, :, slot) = stored_gradient_change
    else
      state%gradient_change(:, :, :, slot) = gradient_change
    end if
    state%inverse_curvature(slot) = inverse_curvature
    state%num_stored = min(state%num_stored + 1, state%history_size)
    state%next_slot = modulo(slot, state%history_size) + 1
    accepted = .true.
    damped = pair_damped

  end subroutine lbfgs_push

  !===========================================================
  subroutine lbfgs_apply(state, vector, inverse_product, error, comm)
    !===========================================================
    !! Apply the L-BFGS inverse-Hessian approximation to `vector`.
    !!
    !! The standard two-loop recursion uses all stored pairs in chronological
    !! order. The initial inverse Hessian is the identity. The returned
    !! `inverse_product` is H_k `vector`; this routine never changes its sign.
    !!
    !! With an empty history the result is simply a copy of `vector`.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(in) :: state
    complex(kind=dp), intent(in) :: vector(:, :, :)
    complex(kind=dp), intent(out) :: inverse_product(:, :, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    real(kind=dp), allocatable :: alpha(:)
    real(kind=dp) :: allocation_failed
    real(kind=dp) :: beta, product
    real(kind=dp) :: distributed_setting(6)
    real(kind=dp) :: invalid_argument(2)
    integer :: chronological_index, ierr, oldest_slot, slot
    logical :: storage_is_valid

    invalid_argument = 0.0_dp
    storage_is_valid = lbfgs_storage_is_valid(state)
    if (.not. storage_is_valid) invalid_argument(1) = 1.0_dp
    if (storage_is_valid) then
      if (.not. lbfgs_shape_matches(state, vector) .or. &
          size(inverse_product, 1) /= size(state%step, 1) .or. &
          size(inverse_product, 2) /= size(state%step, 2) .or. &
          size(inverse_product, 3) /= size(state%step, 3)) invalid_argument(2) = 1.0_dp
    end if
    call comms_allreduce(invalid_argument(1), size(invalid_argument), 'MAX', error, comm)
    if (allocated(error)) return

    if (invalid_argument(1) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS history is not initialized', comm)
      return
    end if
    if (invalid_argument(2) > 0.5_dp) then
      call set_error_fatal(error, 'L-BFGS product vector has incompatible shape', comm)
      return
    end if

    distributed_setting = [real(state%history_size, dp), -real(state%history_size, dp), &
                           real(state%num_stored, dp), -real(state%num_stored, dp), &
                           real(state%next_slot, dp), -real(state%next_slot, dp)]
    call comms_allreduce(distributed_setting(1), size(distributed_setting), 'MAX', error, comm)
    if (allocated(error)) return
    if (distributed_setting(1) /= -distributed_setting(2) .or. &
        distributed_setting(3) /= -distributed_setting(4) .or. &
        distributed_setting(5) /= -distributed_setting(6)) then
      call set_error_fatal(error, 'L-BFGS state is inconsistent across ranks', comm)
      return
    end if

    inverse_product = vector
    if (state%num_stored == 0) return

    allocate (alpha(state%num_stored), stat=ierr)
    allocation_failed = 0.0_dp
    if (ierr /= 0) allocation_failed = 1.0_dp
    call comms_allreduce(allocation_failed, 1, 'MAX', error, comm)
    if (allocated(error)) return
    if (allocation_failed > 0.5_dp) then
      call set_error_alloc(error, 'Error allocating L-BFGS recursion workspace', comm)
      return
    end if

    oldest_slot = modulo(state%next_slot - state%num_stored - 1, state%history_size) + 1

    ! Traverse newest to oldest in the first recursion.
    do chronological_index = state%num_stored, 1, -1
      slot = modulo(oldest_slot + chronological_index - 2, state%history_size) + 1
      product = lbfgs_local_inner_product(state%step(:, :, :, slot), inverse_product)
      call comms_allreduce(product, 1, 'SUM', error, comm)
      if (allocated(error)) return
      alpha(chronological_index) = state%inverse_curvature(slot)*product
      inverse_product = inverse_product - alpha(chronological_index)*state%gradient_change(:, :, :, slot)
    end do

    ! Traverse oldest to newest in the second recursion.
    do chronological_index = 1, state%num_stored
      slot = modulo(oldest_slot + chronological_index - 2, state%history_size) + 1
      product = lbfgs_local_inner_product(state%gradient_change(:, :, :, slot), inverse_product)
      call comms_allreduce(product, 1, 'SUM', error, comm)
      if (allocated(error)) return
      beta = state%inverse_curvature(slot)*product
      inverse_product = inverse_product + &
                        (alpha(chronological_index) - beta)*state%step(:, :, :, slot)
    end do

  end subroutine lbfgs_apply

  !===========================================================
  pure integer function lbfgs_history_count(state)
    !===========================================================
    !! Return the number of secant pairs currently stored.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(in) :: state

    lbfgs_history_count = state%num_stored

  end function lbfgs_history_count

  !===========================================================
  pure integer function lbfgs_history_capacity(state)
    !===========================================================
    !! Return the configured maximum number of secant pairs.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(in) :: state

    lbfgs_history_capacity = state%history_size

  end function lbfgs_history_capacity

  !===========================================================
  pure logical function lbfgs_is_initialized(state)
    !===========================================================
    !! Report whether the history owns a complete allocation.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(in) :: state

    lbfgs_is_initialized = lbfgs_storage_is_valid(state)

  end function lbfgs_is_initialized

  !===========================================================
  pure logical function lbfgs_storage_is_valid(state)
    !===========================================================
    !! Check internal allocation and counter invariants.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(in) :: state

    lbfgs_storage_is_valid = state%initialized .and. &
                             state%history_size > 0 .and. &
                             state%num_stored >= 0 .and. &
                             state%num_stored <= state%history_size .and. &
                             state%next_slot >= 1 .and. &
                             state%next_slot <= state%history_size .and. &
                             allocated(state%step) .and. &
                             allocated(state%gradient_change) .and. &
                             allocated(state%inverse_curvature)

  end function lbfgs_storage_is_valid

  !===========================================================
  pure logical function lbfgs_shape_matches(state, vector)
    !===========================================================
    !! Check the local dimensions of a tangent vector.
    !===========================================================

    implicit none

    type(lbfgs_state_type), intent(in) :: state
    complex(kind=dp), intent(in) :: vector(:, :, :)

    lbfgs_shape_matches = size(vector, 1) == size(state%step, 1) .and. &
                          size(vector, 2) == size(state%step, 2) .and. &
                          size(vector, 3) == size(state%step, 3)

  end function lbfgs_shape_matches

  !===========================================================
  pure real(kind=dp) function lbfgs_local_inner_product(left, right)
    !===========================================================
    !! Return the local real Frobenius product Re Tr(left^H right).
    !===========================================================

    implicit none

    complex(kind=dp), intent(in) :: left(:, :, :)
    complex(kind=dp), intent(in) :: right(:, :, :)

    lbfgs_local_inner_product = sum(real(conjg(left)*right, kind=dp))

  end function lbfgs_local_inner_product

end module w90_lbfgs
