module array_operations_mod
   !! This module provides supplementary functions and subroutines to handle
   !! matrix allocation, invertion, appending etc.
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   implicit none
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   interface allocate_1d
      pure module subroutine allocate_1d_int32(array_, size_)
         !! allocate a 1d array and fill it with 0s (integer version)
         integer(int32), intent(inout), allocatable :: array_(:)
         integer(int32), intent(in)                 :: size_
      end subroutine allocate_1d_int32
      pure module subroutine allocate_1d_sp(array_, size_)
         !! allocate a 1d array and fill it with 0s (single precision version)
         real(sp), intent(inout), allocatable :: array_(:)
         integer(int32), intent(in)           :: size_
      end subroutine allocate_1d_sp
      pure module subroutine allocate_1d_dp(array_, size_)
         !! allocate a 1d array and fill it with 0s (double precision version)
         real(dp), intent(inout), allocatable :: array_(:)
         integer(int32), intent(in)           :: size_
      end subroutine allocate_1d_dp
   end interface
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   interface allocate_2d
      pure module subroutine allocate_2d_int32(array_, size1_, size2_)
         !! allocate a 2d array and fill it with 0s (integer version)
         integer(int32), intent(inout), allocatable :: array_(:,:)
         integer(int32), intent(in)                 :: size1_, size2_
      end subroutine allocate_2d_int32
      pure module subroutine allocate_2d_sp(array_, size1_, size2_)
         !! allocate a 2d array and fill it with 0s (single precision version)
         real(sp), intent(inout), allocatable :: array_(:,:)
         integer(int32), intent(in)           :: size1_, size2_
      end subroutine allocate_2d_sp
      pure module subroutine allocate_2d_dp(array_, size1_, size2_)
         !! allocate a 2d array and fill it with 0s (double precision version)
         real(dp), intent(inout), allocatable :: array_(:,:)
         integer(int32), intent(in)           :: size1_, size2_
      end subroutine allocate_2d_dp
   end interface
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   interface allocate_3d
      pure module subroutine allocate_3d_int32(array_, size1_, size2_, size3_)
         !! allocate a 3d array and fill it with 0s (integer version)
         integer(int32), intent(inout), allocatable :: array_(:,:,:)
         integer(int32), intent(in)                 :: size1_, size2_, size3_
      end subroutine allocate_3d_int32
      pure module subroutine allocate_3d_sp(array_, size1_, size2_, size3_)
         !! allocate a 3d array and fill it with 0s (single precision version)
         real(sp), intent(inout), allocatable :: array_(:,:,:)
         integer(int32), intent(in)           :: size1_, size2_, size3_
      end subroutine allocate_3d_sp
      pure module subroutine allocate_3d_dp(array_, size1_, size2_, size3_)
         !! allocate a 3d array and fill it with 0s (double precision version)
         real(dp), intent(inout), allocatable :: array_(:,:,:)
         integer(int32), intent(in)           :: size1_, size2_, size3_
      end subroutine allocate_3d_dp
   end interface
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   interface append
      pure module subroutine append_int32(array_, element_)
         !! append element to an array (integer version)
         integer(int32), intent(inout), allocatable :: array_(:)
         integer(int32), intent(in)                 :: element_
      end subroutine append_int32
      pure module subroutine append_sp(array_, element_)
         !! append element to an array (single precision version)
         real(sp), intent(inout), allocatable :: array_(:)
         real(sp), intent(in)                 :: element_
      end subroutine append_sp
      pure module subroutine append_dp(array_, element_)
         !! append element to an array (double precision version)
         real(dp), intent(inout), allocatable :: array_(:)
         real(dp), intent(in)                 :: element_
      end subroutine append_dp
   end interface
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   interface invert_symmetric_matrix
      module subroutine invert_symmetric_matrix_int32(matrix_)
         !! invert a symmetric matrix using DSYTRI method (integer version)
         integer(int32), intent(inout) :: matrix_(:,:)
      end subroutine invert_symmetric_matrix_int32
      module subroutine invert_symmetric_matrix_sp(matrix_)
         !! invert a symmetric matrix using DSYTRI method
         !! (single precision version)
         real(sp), intent(inout) :: matrix_(:,:)
      end subroutine invert_symmetric_matrix_sp
      module subroutine invert_symmetric_matrix_dp(matrix_)
         !! invert a symmetric matrix using DSYTRI method
         !! (double precision version)
         real(dp), intent(inout) :: matrix_(:,:)
      end subroutine invert_symmetric_matrix_dp
   end interface
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   interface fill_symmetric_matrix
      module subroutine fill_symmetric_matrix_int32(matrix_, upper_lower_)
         !! fill the upper/lower triangle of a symmetric matrix (integer version)
         integer(int32), intent(inout) :: matrix_(:,:)
         character(len = 1), intent(in) :: upper_lower_
      end subroutine fill_symmetric_matrix_int32
      module subroutine fill_symmetric_matrix_sp(matrix_, upper_lower_)
         !! fill the upper/lower triangle of a symmetric matrix
         !! (single precision version)
         real(sp), intent(inout) :: matrix_(:,:)
         character(len = 1), intent(in) :: upper_lower_
      end subroutine fill_symmetric_matrix_sp
      module subroutine fill_symmetric_matrix_dp(matrix_, upper_lower_)
         !! fill the upper/lower triangle of a symmetric matrix
         !!  (double precision version)
         real(dp), intent(inout) :: matrix_(:,:)
         character(len = 1), intent(in) :: upper_lower_
      end subroutine fill_symmetric_matrix_dp
   end interface
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   interface add_scalar_to_diagonal
      module subroutine add_scalar_to_diagonal_int32(matrix_, scalar_)
         !! add a scalar value to the matrix diagonal (integer version)
         integer(int32), intent(inout) :: matrix_(:,:)
         integer(int32), intent(in) :: scalar_
      end subroutine add_scalar_to_diagonal_int32
      module subroutine add_scalar_to_diagonal_sp(matrix_, scalar_)
         !! add a scalar value to the matrix diagonal
         !! (single precision version)
         real(sp), intent(inout) :: matrix_(:,:)
         real(sp), intent(in) :: scalar_
      end subroutine add_scalar_to_diagonal_sp
      module subroutine add_scalar_to_diagonal_dp(matrix_, scalar_)
         !!! add a scalar value to the matrix diagonal
         !! (double precision version)
         real(dp), intent(inout) :: matrix_(:,:)
         real(dp), intent(in) :: scalar_
      end subroutine add_scalar_to_diagonal_dp
   end interface
   !---------------------------------------------------------------------------!
end module array_operations_mod
