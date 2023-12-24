module io_mod
   !! this module provides following functions and subroutines:
   !! 1. input_file - reads the input file prepared by the user
   !! 2. input_check - checks the variables supplied in the input file
   !! 3. input_summary - summary of the input variables
   !! 5. etotal (function) - returns the total energy of the system
   !! 6. wavenumberekin (function) - returns the wavenumber (in 1/A^{2})
   !! 7. units_conversion - converts all physical quantities to atomic units
   !! 8. count_available_xs (function) - counts energetically accessible
   !!    levels in the basis
   !! 9. jaccess (function) - returns jopen and waveopen - arrays needed for
   !!    calculations of the state-to-state XS
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use supplementary_mod, only: file_io_status, write_error, write_message,    &
      incorrect_value, integer_to_character, float_to_character
   use additional_mod, only: allocate_1d, allocate_2d, allocate_3d
!------------------------------------------------------------------------------!
   implicit none
!------------------------------------------------------------------------------!
   character(len = 80) :: label, potentialfile, smatrixfile, partialfile
   integer(int32) :: ietoterel, jtotmin, jtotmax, jtotstep,                    &
      iabsdr, steps, ncac, initial, nlevel, nr, nterms, ncol, ncoupl, totalcol,&
      n_skip_lines, iunits, ipart, prntlvl, saveeigenvectors, nmodlevels
   real(dp) :: reducedmass, energy, rmin, rmax, dr, vdepth, dtol, otol
   real(dp), parameter :: amutoau = 1822.8884862d0
   real(dp), parameter :: bohrtoangstrom = 0.5291772109d0
   real(dp), parameter :: hartreetocm = 219474.631363d0
   real(dp), parameter :: pi = dacos(-1.d0)
   integer(int32), allocatable :: v1array(:), j1array(:), l1tab(:), l2tab(:),  &
      lltab(:), v1pes(:), j1pes(:), v1ppes(:), j1ppes(:), reduced_v1pes(:),    &
      reduced_j1pes(:), reduced_v1ppes(:), reduced_j1ppes(:)
   real(dp), allocatable :: elevel(:), rmat(:), b(:), c(:), d(:), vmat(:,:),   &
      bmat(:,:), cmat(:,:), dmat(:,:), read_vmat3D(:,:,:), vmat3D(:,:,:),      &
      bmat3D(:,:,:), cmat3D(:,:,:), dmat3D(:,:,:)
   logical :: pes_file_exists
   !------------------------------------------------------------------------------!
   contains

      subroutine read_input_file
         !! reads the input file prepared by the user using NAMELIST feature
         !! the code uses 3 namelists: input, basis and potential
         !------------------------------------------------------------------------!
         character(len = 200):: err_message
         integer(int32) :: ilevel, iilevel, icoupl, icol, il, io_status
         !------------------------------------------------------------------------!
         namelist / INPUT / label, reducedmass, ietoterel, energy,                &
            jtotmin, jtotmax, jtotstep, rmin, rmax, iabsdr, dr, steps, vdepth,    &
            ncac, dtol, otol, nlevel, initial, nr, nterms, ncol, totalcol,        &
            n_skip_lines, iunits, potentialfile, smatrixfile, ipart, partialfile, &
            prntlvl
         namelist / BASIS / v1array, j1array, elevel
         namelist / POTENTIAL / l1tab, l2tab, lltab, v1pes, j1pes, v1ppes, j1ppes
!------------------------------------------------------------------------------!
! Pre-declaration of variables                                                 !
!------------------------------------------------------------------------------!
! The most important variables (if they are not specified, the code stops):    !
!------------------------------------------------------------------------------!
         reducedmass = -1.0_dp
         energy = -1.0_dp
         rmin = -1.0_dp
         rmax = -1.0_dp
         nlevel = -1
         initial = -1
         nr = -1
         nterms = -1
!------------------------------------------------------------------------------!
! Additional variables (the code runs with the pre-declared values):           !
!------------------------------------------------------------------------------!
         ietoterel = 0
         jtotmin = 0
         jtotmax = -1
         jtotstep = 1
         iabsdr = 0
         dr = -1.0_dp
         steps = 10
         vdepth = 0.0_dp
         ncac = 1
         dtol = 0.1_dp
         otol = 0.1_dp
         ncol = 1
         totalcol = 1
         n_skip_lines = 0
         iunits = 0
         potentialfile = 'RadialTerms.dat'
         smatrixfile = 'SmatrixFile.dat'
         ipart = 0
         partialfile = 'PartialFile.dat'
         prntlvl = 2
!------------------------------------------------------------------------------!
         open(unit=5, action='read', form='formatted', access='sequential',    &
            status = 'old', iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'o')
!------------------------------------------------------------------------------!
! Read the input namelist:                                                     !
!------------------------------------------------------------------------------!
         read(unit=5, nml=INPUT, iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'r')
!------------------------------------------------------------------------------!
! Check if the variables from input namelist are supplied correctly:           !
!------------------------------------------------------------------------------!
         call input_check(1)

         if (jtotmax.eq.-1) jtotmax = 999999

         call allocate_1d(v1array,nlevel)
         call allocate_1d(j1array,nlevel)
         call allocate_1d(elevel,nlevel)

         call allocate_1d(l1tab,nterms)
         call allocate_1d(v1pes,totalcol)
         call allocate_1d(v1ppes,totalcol)
         call allocate_1d(j1pes,totalcol)
         call allocate_1d(j1ppes,totalcol)
!------------------------------------------------------------------------------!
! Read the basis namelist & check if the values were supplied correctly:       !
!------------------------------------------------------------------------------!
         read(unit=5, nml=BASIS, iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'r')
         call input_check(2)
!------------------------------------------------------------------------------!
! If itype = 2/4 the code reads all the totalcol coupling terms, but some of   !
! them will not be used in the calculations. Here, the code prepares           !
! the arrays of ncoupl size, that will hold only the necessary terms           !
!------------------------------------------------------------------------------!
         ncoupl = nlevel * (nlevel + 1) / 2

         call allocate_1d(reduced_j1pes,ncoupl)
         call allocate_1d(reduced_j1ppes,ncoupl)
         call allocate_1d(reduced_v1pes,ncoupl)
         call allocate_1d(reduced_v1ppes,ncoupl)

         icoupl = 0

         do ilevel = 1, nlevel
            do iilevel = ilevel, nlevel
               icoupl = icoupl + 1
               reduced_v1pes(icoupl)  = v1array(ilevel)
               reduced_j1pes(icoupl)  = j1array(ilevel)
               reduced_v1ppes(icoupl) = v1array(iilevel)
               reduced_j1ppes(icoupl) = j1array(iilevel)
            enddo
         enddo
!------------------------------------------------------------------------------!
! Read the potential namelist & check if the values were supplied correctly:   !
!------------------------------------------------------------------------------!
         read(unit=5, nml=POTENTIAL, iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'r')
         call input_check(3)

         close(5)
!------------------------------------------------------------------------------!
! Prepare the arrays that are needed for interpolation of the coupling terms:  !
!------------------------------------------------------------------------------!
         call allocate_1d(rmat,nr)
         call allocate_1d(b,nr)
         call allocate_1d(c,nr)
         call allocate_1d(d,nr)

         call allocate_3d(read_vmat3D,nr,nterms,totalcol)
         call allocate_3d(vmat3D,nr,nterms,ncoupl)
         call allocate_3d(bmat3D,nr,nterms,ncoupl)
         call allocate_3d(cmat3D,nr,nterms,ncoupl)
         call allocate_3d(dmat3D,nr,nterms,ncoupl)
!------------------------------------------------------------------------------!
! Summarize the input parameters:                                              !
!------------------------------------------------------------------------------!
         call input_summary

      end subroutine read_input_file
!------------------------------------------------------------------------------!
   subroutine input_check(nmlistind)
      !! checks if the supplied input parameters are correct
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: nmlistind
         !! nmlistind = 1: namelist INPUT
         !! nmlistind = 2: namelist BASIS
         !! nmlistind = 3: namelist POTENTIAL
      !------------------------------------------------------------------------!
      integer(int32) :: ilevel, il, icol
      !------------------------------------------------------------------------!
      if (nmlistind.eq.1) then
      !------------------------------------------------------------------------!
      ! Namelist input:
      !------------------------------------------------------------------------!
         if (reducedmass.lt.0) then
            call incorrect_value("reducedmass", reducedmass, 5)
         endif

         if ((ietoterel.ne.0).and.(ietoterel.ne.1)) then
            call incorrect_value("ietoterel", ietoterel, 5)
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

         if ((iabsdr.ne.0).and.(iabsdr.ne.1)) then
            call incorrect_value("iabsdr", iabsdr, 5)
         endif

         if ((iabsdr.eq.1).and.(dr.lt.0.d0)) then
            call incorrect_value("dr", dr, 5)
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

            if (ncac.lt.0) then
               call write_message("JTOTMAX < 0:")
               call incorrect_value("ncac", ncac, 5)
            endif

            if (dtol.lt.0) then
               call write_message("JTOTMAX < 0:")
               call incorrect_value("dtol", dtol, 5)
            endif

            if (otol.lt.0) then
               call write_message("JTOTMAX < 0:")
               call incorrect_value("otol", otol, 5)
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

         if (ietoterel.eq.1) then

            if (initial.le.0) then
               call write_message("ietoterel = 1:")
               call incorrect_value("initial", initial, 5)
            endif

            if (initial.gt.nlevel) then
               call write_message("ietoterel = 1:")
               call write_message("nlevel = " // integer_to_character(nlevel))
               call incorrect_value("initial > nlevel", initial, 5)
            endif

         endif

         if (nr.le.0) then
            call incorrect_value("nr", nr, 5)
         endif

         if (nterms.le.0) then
            call incorrect_value("nterms", nterms, 5)
         endif

         if (totalcol.le.0) then
            call incorrect_value("totalcol", totalcol, 5)
         endif

         if (n_skip_lines.lt.0) then
            call incorrect_value("n_skip_lines", n_skip_lines, 5)
         endif

         if ((iunits.ne.0).and.(iunits.ne.1)) then
            call incorrect_value("iunits", iunits, 5)
         endif

         inquire(file = potentialfile, exist = pes_file_exists)
         if (pes_file_exists.eqv..false.) then
            call write_error(trim(adjustl(potentialfile)) // "does not exist")
         endif

         if ((ipart.ne.0).and.(ipart.ne.1)) then
            call incorrect_value("ipart", ipart, 5)
         endif

         if (prntlvl.lt.0) then
            call incorrect_value("prntlvl", prntlvl, 5)
         endif
      !------------------------------------------------------------------------!
      else if (nmlistind.eq.2) then
      !------------------------------------------------------------------------!
      ! Namelist basis:
      !------------------------------------------------------------------------!
         do ilevel = 1, nlevel
            if (v1array(ilevel).lt.0) then
               call incorrect_value("v1array(" // integer_to_character(ilevel) // ")", v1array(ilevel), 5)
            endif

            if (j1array(ilevel).lt.0) then
               call incorrect_value("j1array(" // integer_to_character(ilevel) // ")", j1array(ilevel), 5)
            endif

            if (elevel(ilevel).lt.0.0_dp) then
               call incorrect_value("elevel(" // integer_to_character(ilevel) // ")", elevel(ilevel), 5)
            endif
         enddo
      !------------------------------------------------------------------------!
      else if (nmlistind.eq.3) then
      !------------------------------------------------------------------------!
      ! Namelist potential:
      !------------------------------------------------------------------------!
         do il = 1, nterms
            if (l1tab(il).lt.0) then
               call incorrect_value("l1tab(" // integer_to_character(il) // ")", l1tab(il), 5)
            endif
         enddo

         do icol = 1, totalcol
            if (v1pes(icol).lt.0) then
               call incorrect_value("v1pes(" // integer_to_character(icol) // ")", v1pes(icol), 5)
            endif

            if (j1pes(icol).lt.0) then
               call incorrect_value("j1pes(" // integer_to_character(icol) // ")", j1pes(icol), 5)
            endif

            if (v1ppes(icol).lt.0) then
               call incorrect_value("vp1pes(" // integer_to_character(icol) // ")", v1ppes(icol), 5)
            endif

            if (j1ppes(icol).lt.0) then
               call incorrect_value("j1ppes(" // integer_to_character(icol) // ")", j1ppes(icol), 5)
            endif
         enddo

      endif
      !------------------------------------------------------------------------!
   end subroutine input_check
!------------------------------------------------------------------------------!
   subroutine input_summary
      !! summarize the input parameters for the current run
      !------------------------------------------------------------------------!
      integer(int32) :: ilevel
      !------------------------------------------------------------------------!
      call write_message("User-supplied label: " // label)
      call write_message("Reduced mass: " //                                   &
         trim(adjustl(float_to_character(reducedmass, "(F10.4)"))) // " a.m.u.")
      call write_message("*** Energy levels in the basis set: ***")
      
      call write_message("   v       j            Energy (cm^{-1})")
      do ilevel = 1,nlevel
         write(*,"(I4,4X,I4,16X,F12.4)") v1array(ilevel), j1array(ilevel), elevel(ilevel)
      enddo

      write(*,"(44X,A3)") "***"

      if (jtotmax.ne.999999) then
         call write_message("The equations will be solved" //                  &
            "for total angular momentum J from " //                            &
            trim(adjustl(integer_to_character(jtotmin))) // " to "             &
            // trim(adjustl(integer_to_character(jtotmax))) // "with step" //  &
            trim(adjustl(integer_to_character(jtotstep))))
      else
         call write_message("The loop over JTOT will be performed from " //    &
            trim(adjustl(integer_to_character(jtotmin))) // " with step " //   &
            trim(adjustl(integer_to_character(jtotstep))) // " until " //      &
            trim(adjustl(integer_to_character(ncac))) //                       &
            " consecutive JTOT-blocks contribute less than " //                &
            trim(adjustl(float_to_character(dtol, "(E10.4)"))) //              &
            " A^2 to the elastic XS and less than " //                         &
            trim(adjustl(float_to_character(otol, "(E10.4)"))) //              &
            " A^2 to the inelastic XS")
      endif

      if (ietoterel.eq.0) then
         call write_message("The calculations will be performed for the total energy equal to "&
            // trim(adjustl(float_to_character(ETOTAL(), "(F10.4)")))//" cm-1")
      else if(ietoterel.eq.1) then
         call write_message("Relative kinetic energy of the colliding system: " //&
            trim(adjustl(float_to_character(energy, "(F10.4)"))) // " cm-1")
         call write_message("and the corresponding wavenumber: " // &
            trim(adjustl(float_to_character(WAVENUMBEREKIN(), "(F10.4)"))) //  &
            " 1/Ang")
         call write_message("The kinetic energy is calculated with respect to the" //&
            " v = " // trim(adjustl(integer_to_character(v1array(initial)))) //&
            " j = " // trim(adjustl(integer_to_character(j1array(initial)))) //&
            " level in the basis set with the rotational energy " //           &
            trim(adjustl(float_to_character(elevel(initial), "(F10.4)")))      &
            // " cm-1.")
         call write_message("This gives the total energy equal to " // &
            trim(adjustl(float_to_character(ETOTAL(), "(F10.4)"))) // " cm-1")
      endif

      if (ipart.eq.1) then
         call write_message("Partial cross sections will be saved into " // partialfile )
      endif

      call write_message("S-matrix elements will be saved into " // smatrixfile )
      !------------------------------------------------------------------------!
   end subroutine input_summary
!------------------------------------------------------------------------------!
   function ETOTAL() result(etot_)
      !! returns the total energy
      !------------------------------------------------------------------------!
      real(dp) ::  etot_
      !------------------------------------------------------------------------!
      if (ietoterel.eq.0) then
         etot_ = energy
      else if (ietoterel.eq.1) then
         etot_ = energy+elevel(initial)
      endif
      !------------------------------------------------------------------------!
   end function
!------------------------------------------------------------------------------!
   function WAVENUMBEREKIN() result(k_)
      !! returns the wavenumber (in 1/A^{2})
      !------------------------------------------------------------------------!
      real(dp) :: k_
      !------------------------------------------------------------------------!
      if (ietoterel.eq.1) then
         k_ = dsqrt(2*reducedmass*amutoau*&
          (ETOTAL()-elevel(initial))/hartreetocm)/bohrtoangstrom
      else
         call write_error("Wavenumberekin called when ietoterel = 0")
      endif
      !------------------------------------------------------------------------!
   end function
!------------------------------------------------------------------------------!
   subroutine units_conversion
      !! converts all physical quantities to atomic units
      !---------------------------------------------------------------------!
      integer(int32) :: ilevel
      !---------------------------------------------------------------------!
      reducedmass = reducedmass*amutoau
      energy=energy/hartreetocm
      vdepth=vdepth/hartreetocm
      !---------------------------------------------------------------------!
      do ilevel=1,nlevel
         elevel(ilevel)=elevel(ilevel)/hartreetocm
      enddo
      !---------------------------------------------------------------------!
   end subroutine units_conversion
!------------------------------------------------------------------------------!
      function count_open_basis_levels() result(open_)
         !! counts the energetically accessible levels in the basis
         !---------------------------------------------------------------------!
         integer(int32) :: open_, ilevel
         !---------------------------------------------------------------------!
         open_ = 0
         do ilevel = 1, nlevel
            if (ETOTAL()-elevel(ilevel).ge.0) then
               open_ = open_ + 1
            endif
         enddo
         !---------------------------------------------------------------------!
      end function count_open_basis_levels
!------------------------------------------------------------------------------!
      subroutine save_open_basis_levels(number_of_open_basis_levels,           &
         open_basis_levels, open_basis_wavevectors)
         !! saves indices to open levels in the basis and corresponding
         !! wavenumbers
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_open_basis_levels
            !! number of energetically accessible levels in the basis
         integer(int32), intent(inout) :: open_basis_levels(number_of_open_basis_levels)
            !! array holding indices to energetically accessible levels in the basis
         real(dp), intent(inout) :: open_basis_wavevectors(number_of_open_basis_levels)
            !! array holding wavevectors calculated w.r.t energetically accessible levels in the basis
         !---------------------------------------------------------------------!
         integer(int32) :: count_, ilevel
         !---------------------------------------------------------------------!
         count_ = 0
         do ilevel = 1, nlevel
            if (ETOTAL() - elevel(ilevel) >= 0) then
               count_ = count_ + 1
               open_basis_levels(count_) = ilevel
               open_basis_wavevectors(count_) = dsqrt((2 * reducedmass         &
                  * (ETOTAL() - elevel(ilevel)))) / bohrtoangstrom
            endif
         enddo
         !---------------------------------------------------------------------!
      end subroutine save_open_basis_levels
!------------------------------------------------------------------------------!
end module io_mod
