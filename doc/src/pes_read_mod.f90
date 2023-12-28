module POTENTIAL
   !! this module provides functions for reading and interpolating radial coupling
   !! coefficients of the PES:
   !! 1. potential_read: reads radial coupling terms from external file
   !! 2. potential_reduction: reduces the number of radial terms that couple
   !!    different rovibrational levels to a smaller number of terms that
   !!    describe the couplings between the levels provided in the basis
   !! 3. potential_interpolation: cubic spline interpolation of the coupling terms
   !!--------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use utility_functions_mod, only: file_io_status, write_error, write_message,&
      time_count_summary
   use math_functions_mod, only: spline
   implicit none
   contains
   !---------------------------------------------------------------------------!
      subroutine potential_read
         !! reads the radial coupling terms from the external file.
         !! The file is assumed to be formatted as follows:
         !! R, V_{l1,v1,j1,v1',j1'}(R) 
         !! The read radial coupling terms are kept in vmat/read_vmat3D
         !---------------------------------------------------------------------!
         character(len = 200) :: err_message
         integer(int32) :: nrtmp, l1, iskip_, il, ir, icol, io_status
         !---------------------------------------------------------------------!
         open (8,file=trim(potentialfile),form='formatted',status='old',       &
            iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'o')
         !---------------------------------------------------------------------!
         ! Skip the informative lines at the beginning                         
         !---------------------------------------------------------------------!
         do iskip_ = 1, n_skip_lines
            read(8,*) 
         enddo

         do il = 1, nterms
            !------------------------------------------------------------------!
            ! Check if current value of l1 is in line with l1tab           
            !------------------------------------------------------------------!
            read (8,*) l1
            if (l1.ne.l1tab(il)) then
               close(8)
               close(11)
               if (ipart.eq.1) close(12)
               call write_error("potential_read: l1 = " // integer_to_character(l1) //&
                  " differs from expected balue in l1tab (" // integer_to_character(il) // &
                  ") = " // integer_to_character(l1tab(il)))
            endif
            do ir=1,nr
               read(8,*) rmat(ir), (read_vmat3D(ir,il,icol),&
                          icol = 1, totalcol)
            enddo
         enddo

         close(8, iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'c')
         !---------------------------------------------------------------------!
         ! Check if supplied radial terms cover a sufficient range of R
         !---------------------------------------------------------------------!
         if (rmin.lt.rmat(1)) then
            close(11)
            if (ipart.eq.1) close(12)
            call write_error("rmin value provided by the user (" //            &
               float_to_character(rmin, "(F10.4)") // ") is smaller than " //  &
               "rmin supplied in " // trim(adjustl(potentialfile)) // " ( "//  &
               float_to_character(rmat(1), "(F10.4)") // ")")
         else if (rmax.gt.rmat(nr)) then
            close(11)
            if (ipart.eq.1) close(12)
            call write_error("rmax value provided by the user (" //            &
               float_to_character(rmax, "(F10.4)") // ") is larger than " //   &
               "rmax supplied in " // trim(adjustl(potentialfile)) // " ( "//  &
               float_to_character(rmat(nr), "(F10.4)") // ")")
         endif
         !---------------------------------------------------------------------!
      end subroutine
!------------------------------------------------------------------------------!
      subroutine potential_reduction
         !! reducesof the read_vmat3D matrix to retain only the necessary coupling terms
         !! which are kept in reduced_*. The number of necessary coupling terms (ncoupl)
         !! is already determined in io_mod module.
         !! If totalcol = ncoupl, the procedure is ignored
         !---------------------------------------------------------------------!
         integer(int32) :: il, icoupl, ir, icol
         !---------------------------------------------------------------------!
         ! Print the original set of quantum numbers
         !---------------------------------------------------------------------!
         if (prntlvl.ge.2) then
            call write_message("*** Reducing the number of the radial coupling terms... ***")
            if (prntlvl.ge.3) then
               call write_message("Original set of quantum numbers:")
               call write_message("       v1  j1  v1`  j1`")
               do icol = 1, totalcol
                  write(*,"(5X,2(2X,I2),2(2X,I2))") v1pes(icol), j1pes(icol),  &
                   v1ppes(icol), j1ppes(icol)
               enddo
            endif
         endif
         !---------------------------------------------------------------------!
         ! Pick only the necessary terms:
         !---------------------------------------------------------------------!
         do il=1,nterms
            do ir=1,nr
               do icoupl=1,ncoupl
                  do icol = 1, totalcol
                     if ((reduced_j1pes(icoupl) == j1pes(icol)).and.&
                         (reduced_j1ppes(icoupl) == j1ppes(icol)).and.&
                         (reduced_v1pes(icoupl) == v1pes(icol)).and.&
                         (reduced_v1ppes(icoupl) == v1ppes(icol))) then
                        vmat3D(ir,il,icoupl) = read_vmat3D(ir,il,icol)
                     endif
                  enddo
               enddo
            enddo
         enddo
         !---------------------------------------------------------------------!
         ! Print the reduced set of quantum numbers
         !---------------------------------------------------------------------!
         if (prntlvl.ge.2) then
            if (prntlvl.ge.3) then
               call write_message("Reduced set of quantum numbers:")
               call write_message("       v1  j1  v1`  j1`")
               do icoupl = 1, ncoupl
                  write(*,"(5X,2(2X,I2),2(2X,I2))") reduced_v1pes(icoupl),     &
                     reduced_j1pes(icoupl), reduced_v1ppes(icoupl), reduced_j1ppes(icoupl)
               enddo
            endif
            call write_message("*** Reduced " //                               &
               trim(adjustl(integer_to_character(totalcol))) //                &
               " radial terms to " //                                          &
               trim(adjustl(integer_to_character(ncoupl))) // " ***")
         endif
         !---------------------------------------------------------------------!
         ! The original matrix is no longer needed
         !---------------------------------------------------------------------!
         deallocate(read_vmat3D)
         !---------------------------------------------------------------------!
      end subroutine      
!------------------------------------------------------------------------------!
      subroutine potential_interpolation
         !! intertpolates the necessary radial coupling terms using cubic
         !! spline functions. The coefficients are kept in bmat/cmat/dmat matrices
         !---------------------------------------------------------------------!
         integer(int32) :: il, icoupl, ir
         !---------------------------------------------------------------------!
         do il=1,nterms
            do icoupl=1,ncoupl
               call SPLINE(nr,rmat,vmat3D(:,il,icoupl),b,c,d)
               do ir=1,nr
                  bmat3D(ir,il,icoupl) = b(ir)
                  cmat3D(ir,il,icoupl) = c(ir)
                  dmat3D(ir,il,icoupl) = d(ir)
               enddo
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine
!------------------------------------------------------------------------------!
end module POTENTIAL
