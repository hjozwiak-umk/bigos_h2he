submodule (array_operations_mod) array_operations_add_scalar_to_diagonal_submod
   !! a submodule for add_scalar_to_diagonal subroutines
implicit none

   contains

      module subroutine add_scalar_to_diagonal_int32(matrix_, scalar_)
          !! add a scalar value to the matrix diagonal  (integer).
          integer(int32), intent(inout) :: matrix_(:,:)
          integer(int32), intent(in) :: scalar_
          !---------------------------------------------------------------------!
          integer(int32) :: i_, size_
          !---------------------------------------------------------------------!
          size_ = size(matrix_, dim = 1)
          !---------------------------------------------------------------------!
          do i_ = 1, size_
            matrix_(i_, i_) = matrix_(i_, i_) + scalar_
          enddo
          !---------------------------------------------------------------------!
      end subroutine add_scalar_to_diagonal_int32

      module subroutine add_scalar_to_diagonal_sp(matrix_, scalar_)
          !! add a scalar value to the matrix diagonal (single precision).
          real(sp), intent(inout) :: matrix_(:,:)
          real(sp), intent(in) :: scalar_
          !---------------------------------------------------------------------!
          integer(int32) :: i_, size_
          !---------------------------------------------------------------------!
          size_ = size(matrix_, dim = 1)
          !---------------------------------------------------------------------!
          do i_ = 1, size_
            matrix_(i_, i_) = matrix_(i_, i_) + scalar_
          enddo
         !---------------------------------------------------------------------!
      end subroutine add_scalar_to_diagonal_sp

      module subroutine add_scalar_to_diagonal_dp(matrix_, scalar_)
         !! add a scalar value to the matrix diagonal 
         real(dp), intent(inout) :: matrix_(:,:)
         real(dp), intent(in) :: scalar_
         !---------------------------------------------------------------------!
         integer(int32) :: i_, size_
         !---------------------------------------------------------------------!
         size_ = size(matrix_, dim = 1)
         !---------------------------------------------------------------------!
         do i_ = 1, size_
           matrix_(i_, i_) = matrix_(i_, i_) + scalar_
         enddo
         !---------------------------------------------------------------------!
      end subroutine add_scalar_to_diagonal_dp

end submodule array_operations_add_scalar_to_diagonal_submod
