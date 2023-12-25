submodule (additional_mod) fill_symmetric_matrix_mod
   !! a submodule for append subroutines
   use utility_functions_mod, only: to_lowercase
implicit none

   contains

      module subroutine fill_symmetric_matrix_int32(matrix_, upper_lower_)
          !! Fill the upper/lower triangle of a symmetric matrix (integer).
          integer, intent(inout) :: matrix_(:,:)
          character(len = 1), intent(in) :: upper_lower_
          !---------------------------------------------------------------------!
          integer :: i_, size_1_, size_2_, size_
          !---------------------------------------------------------------------!
          size_1_ = size(matrix_, dim = 1)
          size_2_ = size(matrix_, dim = 2)
          if (size_1_ .eq. size_2_) then
             size_ = size_1_
          else
             print *, "Error in fill_symmetric_matrix_int: size in dim = 1 (",size_1_,&
               ") is different than in dim = 2 (", size_2_, ")"
             print *, "Adapt this subroutine to rectangle matrices"
             stop
          endif
          !---------------------------------------------------------------------!
          select case(to_lowercase(upper_lower_))
              case('l')
                  do i_ = 1, size_ - 1
                      matrix_(i_ + 1:size_, i_) = matrix_(i_, i_ + 1:size_)
                  enddo
              case('u')
                  do i_ = 1, size_ - 1
                      matrix_(i_, i_ + 1:size_) = matrix_(i_ + 1:size_, i_)
                  enddo
              case default
                  write(*,*) "Invalid argument to fill_symmetric_matrix_int32 subroutine"
                  write(*,*) "upper_lower_:", upper_lower_
                  write(*,*) "'u' or 'l' expected"
                  stop
          end select
          !---------------------------------------------------------------------!
      end subroutine fill_symmetric_matrix_int32

      module subroutine fill_symmetric_matrix_sp(matrix_, upper_lower_)
         !! Fill the upper/lower triangle of a symmetric matrix (single precision).
         real(sp), intent(inout) :: matrix_(:,:)
         character(len = 1), intent(in) :: upper_lower_
         !---------------------------------------------------------------------!
         integer :: i_, size_1_, size_2_, size_
         !---------------------------------------------------------------------!
         size_1_ = size(matrix_, dim = 1)
         size_2_ = size(matrix_, dim = 2)
         if (size_1_ .eq. size_2_) then
            size_ = size_1_
         else
            print *, "Error in fill_symmetric_matrix_sp: size in dim = 1 (",size_1_,&
               ") is different than in dim = 2 (", size_2_, ")"
            print *, "Adapt this subroutine to rectangle matrices"
            stop
         endif
         !---------------------------------------------------------------------!
         select case(to_lowercase(upper_lower_))
             case('l')
                 do i_ = 1, size_ - 1
                     matrix_(i_ + 1:size_, i_) = matrix_(i_, i_ + 1:size_)
                 enddo
             case('u')
                 do i_ = 1, size_ - 1
                     matrix_(i_, i_ + 1:size_) = matrix_(i_ + 1:size_, i_)
                 enddo
             case default
                 write(*,*) "Invalid argument to fill_symmetric_matrix_sp subroutine"
                 write(*,*) "upper_lower_:", upper_lower_
                 write(*,*) "'u' or 'l' expected"
                 stop
         end select
         !---------------------------------------------------------------------!
      end subroutine fill_symmetric_matrix_sp



      module subroutine fill_symmetric_matrix_dp(matrix_, upper_lower_)
         !! fill the upper/lower triangle of a symmetric matrix
         real(dp), intent(inout) :: matrix_(:,:)
         character(len = 1), intent(in) :: upper_lower_
         !---------------------------------------------------------------------!
         integer :: i_, size_1_, size_2_, size_
         !---------------------------------------------------------------------!
         size_1_ = size(matrix_, dim = 1)
         size_2_ = size(matrix_, dim = 2)
         if (size_1_ .eq. size_2_) then
            size_ = size_1_
         else
            print *, "Error in fill_symmetric_matrix_dp: size in dim = 1 (",size_1_,&
               ") is different than in dim = 2 (", size_2_, ")"
            print *, "Adapt this subroutine to rectangle matrices"
            stop
         endif
         !---------------------------------------------------------------------!
         select case(to_lowercase(upper_lower_))
            case('l')
             do i_ = 1, size_ - 1
                matrix_(i_ + 1 : size_, i_) = matrix_(i_, i_ + 1 : size_)
             enddo
            case('u')
             do i_ = 1, size_ - 1
                matrix_(i_, i_ + 1: size_) = matrix_(i_ + 1 : size_, i_)
             enddo
            case default
               write(*,*) "Invalid argument of the fill_symetric_matrix_dp subroutine"
               write(*,*) "upper_lower_", upper_lower_
               write(*,*) "'u' or 'l' expected"
               stop
         end select
         !---------------------------------------------------------------------!
      end subroutine fill_symmetric_matrix_dp

end submodule fill_symmetric_matrix_mod
