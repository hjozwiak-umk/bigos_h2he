submodule (array_operations_mod) array_operations_append_submod
   !! a submodule for append subroutines
implicit none

contains

      pure module subroutine append_int32(array_, element_)
         !! append element to an array (intger version)
         integer(int32), allocatable, intent(inout) :: array_(:)
         integer(int32), intent(in)              :: element_
         !---------------------------------------------------------------------!
         integer(int32), allocatable :: tmp_array_(:)
         !---------------------------------------------------------------------!
         if (allocated(array_)) then
            allocate(tmp_array_(size(array_)+1))
            tmp_array_(:size(array_)) = array_(:)
            tmp_array_(size(array_)+1)    = element_
            deallocate(array_)
            allocate(array_(size(tmp_array_)))
            call move_alloc(tmp_array_, array_)
         else
            allocate(array_(1))
            array_(1) = element_
         endif
         
      end subroutine append_int32
      pure module subroutine append_sp(array_, element_)
         !! append element to an array (single precision version)
         real(sp), allocatable, intent(inout)  :: array_(:)
         real(sp), intent(in)               :: element_
         !---------------------------------------------------------------------!
         real(sp), allocatable :: tmp_array_(:)
         !---------------------------------------------------------------------!
         if (allocated(array_)) then
            allocate(tmp_array_(size(array_)+1))
            tmp_array_(:size(array_)) = array_(:)
            tmp_array_(size(array_)+1)    = element_
            deallocate(array_)
            allocate(array_(size(tmp_array_)))
            call move_alloc(tmp_array_, array_)
         else
            allocate(array_(1))
            array_(1) = element_
         endif
         
      end subroutine append_sp
      pure module subroutine append_dp(array_, element_)
         !! append element to an array (double precision version)
         real(dp), allocatable, intent(inout)  :: array_(:)
         real(dp), intent(in)               :: element_
         !---------------------------------------------------------------------!
         real(dp), allocatable :: tmp_array_(:)
         !---------------------------------------------------------------------!
         if (allocated(array_)) then
            allocate(tmp_array_(size(array_)+1))
            tmp_array_(:size(array_)) = array_(:)
            tmp_array_(size(array_)+1)    = element_
            deallocate(array_)
            allocate(array_(size(tmp_array_)))
            call move_alloc(tmp_array_, array_)
         else
            allocate(array_(1))
            array_(1) = element_
         endif

      end subroutine append_dp

end submodule array_operations_append_submod
