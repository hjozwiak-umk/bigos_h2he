module statetostateXS
   !! This modules contains the subroutine that calculates the state-to-state XS
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use utility_functions_mod, only: write_message, time_count_summary
   implicit none
   contains
   !---------------------------------------------------------------------------!
      subroutine CROSSSECTION(jj, nopen, number_of_channels,                   &
         number_of_open_basis_levels, open_basis_levels,                       &
         open_basis_wavevectors, srmatrix, simatrix, channels_level_indices,   &
         channels_l_values, xs_array)
         !! calculate the state-to-state cross-section
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: jj
            !! total angular momentum
         integer(int32), intent(in) :: nopen
            !! number of open channels
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: number_of_open_basis_levels
            !! number of all possible state-to-state XS
         integer(int32), intent(in) :: open_basis_levels(number_of_open_basis_levels)
            !! holds indices to the basis arrays which correspond to open channels 
         real(dp), intent(in) :: open_basis_wavevectors(number_of_open_basis_levels)
            !! holds wavenumbers k_{i}
         real(dp), intent(in) :: srmatrix(nopen, nopen), simatrix(nopen, nopen)
            !! real and imaginary parts of the S-matrix
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_l_values(number_of_channels)
            !! holds all values of l
         real(dp), intent(inout) :: xs_array(number_of_open_basis_levels*number_of_open_basis_levels)
            !! array holding all XSs
         !---------------------------------------------------------------------!
         integer(int32) :: jinit, jfin, njoccur, njpoccur, v1lev, j1lev,       &
            v1plev, j1plev, ltmp, lptmp, iinit, ifin, ioccur1, ioccur2, ii,    &
            il, ilp
         real(dp) :: waveinit, xssum, telr, teli, telsq, time_start, time_finish, xs_time
         integer(int32), allocatable :: jblockarr(:), jpblockarr(:)
         !---------------------------------------------------------------------!
         call CPU_TIME(time_start)
         xs_array = 0
         !---------------------------------------------------------------------!
         do iinit = 1, number_of_open_basis_levels
            jinit = open_basis_levels(iinit)
            v1lev = v1array(jinit)
            j1lev = j1array(jinit)
            waveinit = open_basis_wavevectors(iinit)
            do ifin = 1, number_of_open_basis_levels
               jfin = open_basis_levels(ifin)
               v1plev = v1array(jfin)
               j1plev = j1array(jfin)
               njoccur = 0
               njpoccur = 0
               do ii = 1, number_of_channels
                  if ((j1array(channels_level_indices(ii)).eq.j1lev)&
                     .and.(v1array(channels_level_indices(ii)).eq.v1lev)) then
                     njoccur = njoccur+1
                  endif
                  if (j1array(channels_level_indices(ii)).eq.j1plev&
                     .and.(v1array(channels_level_indices(ii)).eq.v1plev)) then
                     njpoccur = njpoccur+1
                  endif
               enddo

               call allocate_1d(jblockarr,njoccur)
               call allocate_1d(jpblockarr,njpoccur)

               ioccur1 = 0
               ioccur2 = 0
               do ii = 1,number_of_channels
                  if ((j1array(channels_level_indices(ii)).eq.j1lev)&
                     .and.(v1array(channels_level_indices(ii)).eq.v1lev)) then
                     ioccur1 = ioccur1+1
                     jblockarr(ioccur1) = ii
                  endif
                  if (j1array(channels_level_indices(ii)).eq.j1plev&
                     .and.(v1array(channels_level_indices(ii)).eq.v1plev)) then
                     ioccur2 = ioccur2+1
                     jpblockarr(ioccur2) = ii
                  endif
               enddo

               xssum = 0.d0

               do il = 1,njoccur
                  ltmp = channels_l_values(jblockarr(il))
                  do ilp = 1,njpoccur
                     lptmp = channels_l_values(jpblockarr(ilp))
                     if ((j1lev.eq.j1plev).and.(v1lev.eq.v1plev)) then
                        if (ltmp.eq.lptmp) then
                           telr = 1.d0-srmatrix(jblockarr(il),&
                                                jpblockarr(ilp))
                           teli = -simatrix(jblockarr(il),jpblockarr(ilp))
                           telsq = telr**2.d0+teli**2.d0
                           xssum = xssum+telsq
                        else
                           telr = -srmatrix(jblockarr(il),jpblockarr(ilp))
                           teli = -simatrix(jblockarr(il),jpblockarr(ilp))
                           telsq = telr**2.d0+teli**2.d0
                           xssum = xssum+telsq
                        endif
                     else
                        telr = -srmatrix(jblockarr(il),jpblockarr(ilp))
                        teli = -simatrix(jblockarr(il),jpblockarr(ilp))
                        telsq = telr**2.d0+teli**2.d0
                        xssum = xssum+telsq
                     endif
                  enddo
               enddo

               xssum = xssum*(2.d0*jj+1.d0)

               xs_array((iinit-1)*number_of_open_basis_levels+ifin) = xssum*pi/&
                 ((2.d0*j1lev+1.d0)*waveinit**2.d0)
            enddo
         enddo
         !---------------------------------------------------------------------!
         CALL CPU_TIME(time_finish)

         if (prntlvl.ge.2) call time_count_summary(time_start, time_finish,    &
            xs_time, "Cross-sections calculations completed in ")

      end
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
      subroutine print_largest_partial_xs(jj, maxXSdiag, maxXSoff, jinddiag, jindoff1,  &
         jindoff2, number_of_open_basis_levels, open_basis_levels)
         !! print the largest partial elastic and inelastic state-to-state XS
         !! in given J-block
         !---------------------------------------------------------------------------!
         integer(int32), intent(in) :: jj
            !! total angular momentum
         real(dp), intent(in) :: maxXSdiag
            !! the largest partial elastic state-to-state XS in this J-block
         real(dp), intent(in) :: maxXSoff
            !! the largest partial inelastic state-to-state XS in this J-block
         integer(int32), intent(in) :: jinddiag
            !! index pointing indirectly to quantum numbers associated with
            !! the largest partial elastic state-to-state XS in this J-block
         integer(int32), intent(in) :: jindoff1, jindoff2
            !! indices pointing indirectly to quantum numbers associated with
            !! the largest partial inelastic state-to-state XS in this J-block
         integer(int32), intent(in) :: number_of_open_basis_levels
            !! number of all possible state-to-state XS (size of open_basis_levels array)
         integer(int32), intent(in) :: open_basis_levels(number_of_open_basis_levels)
            !! holds indices to the basis arrays which correspond to open channels 
         !---------------------------------------------------------------------------!
         character(len=200) :: header_line, line
         !---------------------------------------------------------------------------!
         if ((prntlvl.eq.1).or.(prntlvl.eq.2)) then
            call write_message("Largest partial elastic state-to-state for JTOT = " // &
               integer_to_character(jj) // ": " // float_to_character(maxXSdiag))
            call write_message("Largest partial inelastic state-to-state for JTOT = " // &
               integer_to_character(jj) // ": " // float_to_character(maxXSoff))
            call write_message(repeat(" ", 43) // "***")
         else if (prntlvl.ge.3) then
            !------------------------------------------------------------------------!
            call write_message("Largest partial elastic state-to-state for JTOT = " // &
               integer_to_character(jj) )
            write(header_line, "(2x,a4,2x,a4,2x,a2,2x,a4,2x,a4,16x,a2)") "v1_f",     &
               "j1_f", "<-", "v1_i", "j1_i", "XS"
            call write_message(header_line)
            write(line, "(2X,I4,2X,I4,6X,I4,2X,I4,2X,E16.8)")                        &
               v1array(open_basis_levels(jinddiag)), j1array(open_basis_levels(jinddiag)),                   &
               v1array(open_basis_levels(jinddiag)), j1array(open_basis_levels(jinddiag)), maxXSdiag
            call write_message(line)
            !------------------------------------------------------------------------!
            call write_message("Largest partial inelastic state-to-state for JTOT = " // &
               integer_to_character(jj) )
            call write_message(header_line)
            write(line, "(2X,I4,2X,I4,6X,I4,2X,I4,2X,E16.8)")                        &
               v1array(open_basis_levels(jindoff2)), j1array(open_basis_levels(jindoff2)),                   &
               v1array(open_basis_levels(jindoff1)), j1array(open_basis_levels(jindoff1)), maxXSoff
            call write_message(line)
            call write_message(repeat(" ", 43) // "***")
            !------------------------------------------------------------------------!
         endif
      end subroutine print_largest_partial_xs
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
      subroutine check_dtol_otol(maxXSdiag, maxXSoff, ncacdiag, ncacoff, terminate)
         !! check if the dtol/otol condition on partial XS is already fulfilled
         !---------------------------------------------------------------------------!
         real(dp), intent(in)    :: maxXSdiag, maxXSoff
            !! largest elastic and inelastic XS
         integer(int32), intent(inout) :: ncacdiag, ncacoff
            !! number of consecutive blocks for which dtol/otol condition is already
            !! fulfilled. Gets incremented within the subroutine
         logical, intent(inout) :: terminate
            !! if .true. the dtol/otol condition is fulfilled;
            !! J-tot loop is terminated
         !---------------------------------------------------------------------------!
         integer(int32) :: icount, icount2
         logical :: diagcontr, offcontr
         !---------------------------------------------------------------------------!
         terminate = .false.
         !---------------------------------------------------------------------------!
         ! diagcontr and offcontr check if dtol and otol are fulfilled, respectively
         !---------------------------------------------------------------------------!
         diagcontr = .false.
         offcontr = .false.
         !---------------------------------------------------------------------------!
         if (maxXSdiag.le.dtol) diagcontr = .true.
         if (maxXSoff.le.otol) offcontr = .true.
         !---------------------------------------------------------------------------!
         if (diagcontr) then
            ncacdiag = ncacdiag + 1
         else
            ncacdiag = 0
         endif
         !------------------------------------------------------------------------------!
         if (offcontr) then
            ncacoff = ncacoff + 1
         else
            ncacoff = 0
         endif
         !------------------------------------------------------------------------------!
         ! Finish the calculations if both contributions are smaller than the limits:   !
         !------------------------------------------------------------------------------!
         if ((ncacdiag.ge.ncac).and.(ncacoff.ge.ncac)) terminate = .true.
         !------------------------------------------------------------------------------!
      end subroutine check_dtol_otol
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
end module statetostateXS
