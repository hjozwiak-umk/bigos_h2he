module input_validation
   !! This module provides subroutines validating read values
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_error, write_message,&
      incorrect_value, integer_to_character
   use data_mod
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   contains
   !---------------------------------------------------------------------------!
      subroutine check_namelist_input
         !! Check variables read from namelist "input"
         !---------------------------------------------------------------------!
         if (reduced_mass.lt.0) then
            call incorrect_value("reduced_mass", reduced_mass, 5)
         endif

         if ((relative_energy_flag.ne.0).and.(relative_energy_flag.ne.1)) then
            call incorrect_value("relative_energy_flag", relative_energy_flag, 5)
         endif

         if (energy.lt.0) then
            call incorrect_value("energy", energy, 5)
         endif

         if (rmin.le.0) then
            call incorrect_value("rmin", rmin, 5)
         endif

         if (rmax.le.0) then
            call incorrect_value("rmax", rmax, 5)
         endif

         if (rmax.lt.rmin) then
            call incorrect_value("rmax/rmin", rmax / rmin, 5)
         endif

         if (steps.le.0.d0) then
            call incorrect_value("steps", steps, 5)
         endif

         if (vdepth.lt.0.d0) then
            call incorrect_value("vdepth", vdepth, 5)
         endif

         if (jtotmin.lt.0) then
            call incorrect_value("jtotmin", jtotmin, 5)
         endif

         if (jtotmax.lt.0) then

            if (consecutive_blocks_threshold.lt.0) then
               call write_message("JTOTMAX < 0:")
               call incorrect_value("consecutive_blocks_threshold", consecutive_blocks_threshold, 5)
            endif

            if (elastic_xs_threshold.lt.0) then
               call write_message("JTOTMAX < 0:")
               call incorrect_value("elastic_xs_threshold", elastic_xs_threshold, 5)
            endif

            if (inelastic_xs_threshold.lt.0) then
               call write_message("JTOTMAX < 0:")
               call incorrect_value("inelastic_xs_threshold", inelastic_xs_threshold, 5)
            endif

         else
         
            if (jtotmax.lt.jtotmin) then
               call write_message("jtotmax is smaller than jtotmin")
               call incorrect_value("jtotmax/jtotmin",                         &
                  real(jtotmax/jtotmin, dp), 5)
            endif

         endif

         if (nlevel.le.0) then
            call incorrect_value("nlevel", nlevel, 5)
         endif

         if (relative_energy_flag.eq.1) then

            if (initial.le.0) then
               call write_message("relative_energy_flag = 1:")
               call incorrect_value("initial", initial, 5)
            endif

            if (initial.gt.nlevel) then
               call write_message("relative_energy_flag = 1:")
               call write_message("nlevel = " //                               &
                  trim(adjustl(integer_to_character(nlevel))))
               call incorrect_value("initial > nlevel", initial, 5)
            endif

         endif

         if (nr.le.0) then
            call incorrect_value("nr", nr, 5)
         endif

         if (nterms.le.0) then
            call incorrect_value("nterms", nterms, 5)
         endif

         if (total_number_of_coupling_terms.le.0) then
            call incorrect_value("total_number_of_coupling_terms", total_number_of_coupling_terms, 5)
         endif

         if (n_skip_lines.lt.0) then
            call incorrect_value("n_skip_lines", n_skip_lines, 5)
         endif

         if ((iunits.ne.0).and.(iunits.ne.1)) then
            call incorrect_value("iunits", iunits, 5)
         endif

         inquire(file = potentialfile, exist = pes_file_exists)
         if (pes_file_exists.eqv..false.) then
            call write_error(trim(adjustl(potentialfile)) // " does not exist")
         endif

         if (prntlvl.lt.0) then
            call incorrect_value("prntlvl", prntlvl, 5)
         endif
         !---------------------------------------------------------------------!
      end subroutine check_namelist_input
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine check_namelist_basis
         !! Check variables read from namelist "basis"
         !---------------------------------------------------------------------!
         integer(int32) :: level_index_
         !---------------------------------------------------------------------!
         do level_index_ = 1, nlevel
            if (v1array(level_index_).lt.0) then
               call incorrect_value("v1array(" //                              &
                  integer_to_character(level_index_) // ")", v1array(level_index_), 5)
            endif

            if (j1array(level_index_).lt.0) then
               call incorrect_value("j1array(" //                              &
                  integer_to_character(level_index_) // ")", j1array(level_index_), 5)
            endif

            if (elevel(level_index_).lt.0.0_dp) then
               call incorrect_value("elevel(" //                              &
                  integer_to_character(level_index_) // ")", elevel(level_index_), 5)
            endif
         enddo

      end subroutine check_namelist_basis
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine check_namelist_potential
         !! Check variables read from namelist "potential"
         !---------------------------------------------------------------------!
         integer(int32) :: legendre_index_, column_index_
         !---------------------------------------------------------------------!
         do legendre_index_ = 1, nterms
            if (l1tab(legendre_index_).lt.0) then
               call incorrect_value("l1tab(" //                                &
                  integer_to_character(legendre_index_) // ")",                &
                  l1tab(legendre_index_), 5)
            endif
         enddo

         do column_index_ = 1, total_number_of_coupling_terms
            if (v1pes(column_index_).lt.0) then
               call incorrect_value("v1pes(" //                                &
                  integer_to_character(column_index_) // ")",                  &
                  v1pes(column_index_), 5)
            endif

            if (j1pes(column_index_).lt.0) then
               call incorrect_value("j1pes(" //                                &
                  integer_to_character(column_index_) // ")",                  &
                  j1pes(column_index_), 5)
            endif

            if (v1ppes(column_index_).lt.0) then
               call incorrect_value("vp1pes(" //                               &
                  integer_to_character(column_index_) // ")",                  &
                  v1ppes(column_index_), 5)
            endif

            if (j1ppes(column_index_).lt.0) then
               call incorrect_value("j1ppes(" //                               &
                  integer_to_character(column_index_) // ")",                  &
                  j1ppes(column_index_), 5)
            endif
         enddo

      end subroutine check_namelist_potential
   !------------------------------------------------------------------------------!
end module input_validation
