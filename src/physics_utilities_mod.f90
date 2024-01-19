module physics_utilities_mod
   !! This module provides helper functions: "units_conversion", "total_energy()"
   !! "wavevector_squared_from_energy", and functions that count and save
   !! open levels in the rovibrational basis
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_error
   use array_operations_mod, only: allocate_1d
   use data_mod
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   contains
   !---------------------------------------------------------------------------!
      subroutine units_conversion
         !! converts all physical quantities to atomic units
         !---------------------------------------------------------------------!
         integer(int32) :: level_index_1_
         !---------------------------------------------------------------------!
         reduced_mass = reduced_mass*amutoau
         energy=energy/hartreetocm
         vdepth=vdepth/hartreetocm
         !---------------------------------------------------------------------!
         do level_index_1_=1,nlevel
            elevel(level_index_1_)=elevel(level_index_1_)/hartreetocm
         enddo
         !---------------------------------------------------------------------!
         units_converted = .true.
         !---------------------------------------------------------------------!
      end subroutine units_conversion
!------------------------------------------------------------------------------!
      function total_energy() result(etot_)
         !! returns the total energy
         !---------------------------------------------------------------------!
         real(dp) ::  etot_
         !---------------------------------------------------------------------!
         if (relative_energy_flag.eq.0) then
            etot_ = energy
         else if (relative_energy_flag.eq.1) then
            etot_ = energy+elevel(initial)
         endif
         !---------------------------------------------------------------------!
      end function total_energy
!------------------------------------------------------------------------------!
      function wavevector_squared_from_energy(energy_) result(k_)
         !! returns the squared wavevector, \\(k_{a}^{2}\\),
         !! given the energy of a given state, \\(E_{a}\\);
         !! calls etot() function; atomic units in the whole function
         !! \\( k_{a} = \sqrt(2 \mu (E_{tot} - E_{a}) \\)
         !! since it uses reduced_mass and total_energy(), the function checks
         !! if units are already converted
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: energy_
            !! energy of a given state, \\( E_{a} \\), in a.u.
         real(dp) :: k_
            !! wavevector, \\(k_{a}\\), in a.u.
         !---------------------------------------------------------------------!
         if (units_converted) then
            k_ = 2*reduced_mass*(total_energy() - energy_)
         else
            call write_error("wavevector_squared_from_energy called but units are not " //&
               "converted yet")
         endif
         !---------------------------------------------------------------------!
      end function wavevector_squared_from_energy
!------------------------------------------------------------------------------!
      function is_open(energy_) result(is_open_)
         !! checks if a channel/level is energetically accessible (open)
         !! by comparing energy with total_energy
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: energy_
            !! level/channel energy
         logical :: is_open_
         !---------------------------------------------------------------------!
         is_open_ = ( energy_ <= total_energy() )
         !---------------------------------------------------------------------!
      end function is_open
!------------------------------------------------------------------------------!
      function count_open_basis_levels() result(open_)
         !! counts the energetically accessible levels in the basis
         !---------------------------------------------------------------------!
         integer(int32) :: open_, level_index_1_
         !---------------------------------------------------------------------!
         open_ = 0
         do level_index_1_ = 1, nlevel
            if (is_open(elevel(level_index_1_))) open_ = open_ + 1
         enddo
         !---------------------------------------------------------------------!
      end function count_open_basis_levels
!------------------------------------------------------------------------------!
      subroutine save_open_basis_levels(number_of_open_basis_levels,           &
         open_basis_levels, basis_wavevectors)
         !! saves indices to open levels in the basis and corresponding
         !! wavevectors (in A^2)
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_open_basis_levels
            !! number of energetically accessible levels in the basis
         integer(int32), intent(inout), allocatable :: open_basis_levels(:)
            !! array holding indices to energetically accessible levels in the basis
         real(dp), intent(inout), allocatable :: basis_wavevectors(:)
            !! array holding wavevectors calculated w.r.t energetically accessible levels in the basis
         !---------------------------------------------------------------------!
         integer(int32) :: count_, level_index_1_
         !---------------------------------------------------------------------!
         call allocate_1d(open_basis_levels, number_of_open_basis_levels)
         call allocate_1d(basis_wavevectors, number_of_open_basis_levels)
         !---------------------------------------------------------------------!
         count_ = 0
         do level_index_1_ = 1, nlevel
            if (is_open(elevel(level_index_1_))) then
               count_ = count_ + 1
               open_basis_levels(count_) = level_index_1_
               basis_wavevectors(count_) =                                     &
                  sqrt( wavevector_squared_from_energy(elevel(level_index_1_)) ) / bohrtoangstrom 
            endif
         enddo
         !---------------------------------------------------------------------!
      end subroutine save_open_basis_levels
!------------------------------------------------------------------------------!
end module physics_utilities_mod
