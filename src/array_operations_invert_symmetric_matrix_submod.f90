submodule (array_operations_mod) array_operations_invert_symmetric_matrix_submod
   !! a submodule for append subroutines
   use utility_functions_mod, only: integer_to_character, write_message,       &
      write_error
   implicit none
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      module subroutine invert_symmetric_matrix_int32(matrix_)
         !! invert a symmetric matrix using DSYTRI method (integer version)
         integer(int32), intent(inout) :: matrix_(:,:)
         !---------------------------------------------------------------------!
         call write_error("Integer version of Lapack inversion procedures " // &
            "does not exist: convert the integer arrays to real")
         !---------------------------------------------------------------------!
      end subroutine invert_symmetric_matrix_int32
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      module subroutine invert_symmetric_matrix_sp(matrix_)
         !! invert a symmetric matrix using SSYTRI method
         !! (single precision version)
         real(sp), intent(inout) :: matrix_(:,:)
         !---------------------------------------------------------------------!
         integer(int32) :: size_1_, size_2_, size_, lwork_, nb_, ok_, ILAENV
         integer(int32), allocatable :: pivot_(:)
         real(sp), allocatable :: work_(:)
         !---------------------------------------------------------------------!
         size_1_ = size(matrix_, dim = 1)
         size_2_ = size(matrix_, dim = 2)
         if (size_1_ .eq. size_2_) then
            size_ = size_1_
         else
            call write_message("Error in invert_symmetric_matrix_sp: size "    &
               // "in dim = 1 ("//trim(adjustl(integer_to_character(size_1_))) &
               // ") is different than in dim = 2 (" //                        &
               trim(adjustl(integer_to_character(size_2_))) // ")")
            call write_error("Adapt this subroutine to rectangle matrices")
         endif
         !---------------------------------------------------------------------!
         if (allocated(pivot_)) deallocate(pivot_)
         allocate(pivot_(size_))
         !---------------------------------------------------------------------!
         nb_ = ILAENV(1,'SSYTRF','L',size_,-1,-1,-1)
         lwork_ = nb_ * size_
         if (allocated(work_)) deallocate(work_)
         allocate(work_(lwork_))
         work_ = 0
         !---------------------------------------------------------------------!
         call SSYTRF('L',size_,matrix_,size_,pivot_,work_,lwork_,ok_)
         if (ok_ /= 0) then
            call write_error("SSYTRF failed with status: " //                  &
               trim(adjustl(integer_to_character(ok_))) )
         endif
         !---------------------------------------------------------------------!
         call SSYTRI('L',size_,matrix_,size_,pivot_,work_,ok_)
         if (ok_ /= 0) then
            call write_error("SSYTRI failed with status: " //                  &
               trim(adjustl(integer_to_character(ok_))) )
         endif
         !---------------------------------------------------------------------!
      end subroutine invert_symmetric_matrix_sp
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      module subroutine invert_symmetric_matrix_dp(matrix_)
         !! invert a symmetric matrix using DSYTRI method (double precision version)
         real(dp), intent(inout) :: matrix_(:,:)
         !---------------------------------------------------------------------!
         integer(int32) :: size_1_, size_2_, size_, lwork_, nb_, ok_, ILAENV, i_size_, i_size_2_
         integer(int32), allocatable :: pivot_(:)
         real(dp), allocatable :: work_(:)
         !---------------------------------------------------------------------!
         size_1_ = size(matrix_, dim = 1)
         size_2_ = size(matrix_, dim = 2)
         if (size_1_ .eq. size_2_) then
            size_ = size_1_
         else
            call write_message("Error in invert_symmetric_matrix_dp: size "    &
               // "in dim = 1 ("//trim(adjustl(integer_to_character(size_1_))) &
               // ") is different than in dim = 2 (" //                        &
               trim(adjustl(integer_to_character(size_2_))) // ")")
            call write_error("Adapt this subroutine to rectangle matrices")
         endif
         !---------------------------------------------------------------------!
         if (allocated(pivot_)) deallocate(pivot_)
         allocate(pivot_(size_))
         !---------------------------------------------------------------------!
         nb_ = ILAENV(1,'DSYTRF','L',size_,-1,-1,-1)
         lwork_ = nb_ * size_
         if (allocated(work_)) deallocate(work_)
         allocate(work_(lwork_))
         work_ = 0
         !---------------------------------------------------------------------!
         call DSYTRF('L',size_,matrix_,size_,pivot_,work_,lwork_,ok_)
         if (ok_ /= 0) then
            call write_error("DSYTRF failed with status: " //                  &
               trim(adjustl(integer_to_character(ok_))) )
         endif
         !---------------------------------------------------------------------!
         call DSYTRI('L',size_,matrix_,size_,pivot_,work_,ok_)
         if (ok_ /= 0) then
            call write_error("DSYTRI failed with status: " //                  &
               trim(adjustl(integer_to_character(ok_))) )
         endif
         !---------------------------------------------------------------------!
      end subroutine invert_symmetric_matrix_dp
   !---------------------------------------------------------------------------!
end submodule array_operations_invert_symmetric_matrix_submod
