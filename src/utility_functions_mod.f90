module utility_functions_mod
   !! utility_functions_mod contains functions which handle writing 
   !! messages/errors/warnings on screen, formatting headers, summary of the 
   !! calculations and a few other supporting functions.
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64, &
                                            output_unit
   !---------------------------------------------------------------------------!
   implicit none
   private
   public :: write_header, write_message, write_warning, write_error, time_count_summary,    &
             alloc_status, file_io_status, incorrect_value, to_lowercase,      &
             integer_to_character, float_to_character
   !---------------------------------------------------------------------------!
   character(len=*), parameter :: letters   =                                  &
                          "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
   character(len=*), parameter :: uppercase = letters(1:26)
   character(len=*), parameter :: lowercase = letters(27:)
   !---------------------------------------------------------------------------!
   interface incorrect_value
      !! interface for the following message:
      !! ``incorrect value encountered:
      !!   variable_name = variable_value``
      module procedure incorrect_value_ch
         !! for character variables
      module procedure incorrect_value_int32
         !! for integer variables
      module procedure incorrect_value_sp
         !! for single precision variables
      module procedure incorrect_value_dp
         !! for double precision variables
   end interface
   !---------------------------------------------------------------------------!
   contains
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine write_message(message_, unit_)
         !! writes a message on a chosen unit
         !---------------------------------------------------------------------!
         character(len = *), intent(in)       :: message_
            !! a message to be written
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         if (present(unit_)) then
            write(unit_, '(a)') trim(message_)
         else
            write(output_unit, '(a)') trim(message_)
         endif
         !---------------------------------------------------------------------!
      end subroutine write_message
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine write_warning(message_, unit_)
         !! writes a warning message on a chosen unit
         !---------------------------------------------------------------------!
         character(len = *), intent(in)       :: message_
            !! a message to be written
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         call write_message('Warning: '//trim(message_), unit_)
         !---------------------------------------------------------------------!
      end subroutine write_warning
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine write_error(message_, unit_)
         !! writes an error message on a chosen unit
         !---------------------------------------------------------------------!
         character(len = *), intent(in)       :: message_
            !! a message to be written
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         call write_message('Error: '//trim(message_), unit_)
         stop
         !---------------------------------------------------------------------!
      end subroutine write_error
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine write_header(header_type, opt_integer_)
         !! writes headers on screen
         character(len = *), intent(in) :: header_type
            !! specifies the type of the header: 'main', 'input_read',
            !! 'input_check', 'input_summary', 'initialization', 'check_norm',
            !! 'save_basis', 'save_pes', 'radial_terms', 'save_radial_terms',
            !! 'reconstruction'
         integer(int32), optional, intent(in) :: opt_integer_
            !! optional integer used in case "block" to pass jtot value
         !---------------------------------------------------------------------!
         character(len = 100) :: header_star, header_str
         character(len = 10)  :: tmp_str_
         integer(int32) :: len_str_
         !---------------------------------------------------------------------!
         select case(trim(header_type))
            case('main')
               write(header_star, fmt = "(a90)") repeat("*", 90)
               call write_message(header_star)
               call write_message(header_star)
               write(header_str, fmt = '(a,25x,a43,20x,a)')                    &
                       '*','BIGOS quantum scattering package, vs. 0.01.','*'
               call write_message(header_str)
               write(header_str, fmt = '(a,36x,a19,33x,a)')                    &
                       '*', 'the SCATTERING code','*'
               call write_message(header_str)
               write(header_str, fmt = '(a,29x,a31,28x,a)')                    &
                       '*', 'adjusted for H2-He calculations','*'
               call write_message(header_str)
               write(header_str, fmt = '(a,37x,a17,34x,a)')                    &
                       '*', 'by Hubert Jozwiak','*'
               call write_message(header_str)
               write(header_str, fmt = '(a,40x,a11,37x,a)')                    &
                       '*', '20/12/2023 ','*'
               call write_message(header_str)
               call write_message(header_star)
            case('block')
               call write_message(repeat('*', 90))
               if (present(opt_integer_)) then
                  write(tmp_str_, "(i10)") opt_integer_
                  len_str_ = len_trim(tmp_str_)
                  write(*, '("*", A, "JTOT = ", A, A, "*")')                   &
                     repeat(' ', 40 - len_str_), tmp_str_, repeat(' ', 41)
                  call write_message(repeat('*', 90))
               else
                  call write_error("**** JTOT value not provided in " //       &
                     "write_header_block ****")
               endif
            case default
               call incorrect_value('header_type (write_header)', header_type)
         end select
         !---------------------------------------------------------------------!
      end subroutine write_header
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine time_count_summary(start_, stop_, time_, message_)
         !! print the message about the time it took to complete a single task
         !---------------------------------------------------------------------!
         real(dp), intent(in)  :: start_
            !! initial time
         real(dp), intent(in)  :: stop_
            !! final time
         real(dp), intent(out) :: time_
            !! stop_ - start_
         character(len = *), optional, intent(in) :: message_
            !! (optional) a message to print instead of a default
            !! "Completed in ... s"
         !---------------------------------------------------------------------!
         character(len = 12)  :: default_message = 'Completed in'
         character(len = 100) :: time_msg
         !---------------------------------------------------------------------!
         time_ = stop_ - start_
         if (present(message_)) then
            write(time_msg, fmt='(a,x,a,es11.4,x,a)')                          &
                    '--', trim(message_), time_, 's'
         else
            write(time_msg, fmt='(a,x,a,es11.4,x,a)')                          &
                    '--', default_message, time_, 's'
         endif

         call write_message(time_msg)
         !---------------------------------------------------------------------!
      end subroutine time_count_summary
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine alloc_status(istat_, message_, op_, unit_)
         !! check the status after allocation
         !---------------------------------------------------------------------!
         integer(int32)                       :: istat_
            !! result of stat=istat in (de)allocate
         character(len = *), intent(in)       :: message_
            !! a message to be written
         character(len = 1), intent(in)       :: op_
            !! 'a' for allocation, 'd' for deallocation
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         character(len = :), allocatable :: add_prefix_
         !---------------------------------------------------------------------!
         add_prefix_ = ''
         if (istat_ /= 0) then
            select case(op_)
               case('a')
                  add_prefix_ = 'memory allocation: '//trim(message_)
               case('d')
                  add_prefix_ = 'memory deallocation: '//trim(message_)
               case default
                  call write_error                                             &
                         ('Incorrect op_ argument in alloc_status subroutine ('&
                          //trim(op_)//')')
            end select
            call write_error(add_prefix_, unit_)
         endif
         !---------------------------------------------------------------------!
      end subroutine alloc_status
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine file_io_status(istat_, iomsg_, channel_, op_, unit_)
         !! check the status during various io operations on files
         !---------------------------------------------------------------------!
         integer(int32)                       :: istat_
            !! result of iostat in open/read/write/close
         character(len = *), intent(in)       :: iomsg_
            !! result of iomsg in open/read/write/close
         integer(int32), intent(in)           :: channel_
            !! name of the file
         character(len = 1), intent(in)       :: op_
            !! 'o' for opening of the file, 'r' for reading, 'w' for writing,
            !! 'c' for closing
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         character(len = :), allocatable :: add_prefix_
         !---------------------------------------------------------------------!
         add_prefix_ = ''
         if (istat_ /= 0) then
            select case(op_)
               case('o')
                  add_prefix_ = 'opening file on channel: '//                  &
                     integer_to_character(channel_)
               case('r')
                  add_prefix_ = 'reading file on channel: '//                  &
                     integer_to_character(channel_)
               case('w')
                  add_prefix_ = 'writing to file on channel: '//               &
                     integer_to_character(channel_)
               case('c')
                  add_prefix_ = 'closing file on channel: '//                  &
                     integer_to_character(channel_)
               case default
                  call write_error                                             &
                       ('Incorrect op_ argument in file_io_status subroutine ('&
                          //trim(op_)//')')
            end select
            call write_error(trim(add_prefix_) // " with message: " // iomsg_, unit_)
         endif
         !---------------------------------------------------------------------!
      end subroutine file_io_status
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine incorrect_value_ch(name_, value_, unit_)
         !! ``incorrect value encountered:
         !!   variable_name = variable_value``
         !---------------------------------------------------------------------!
         character(len = *), intent(in)       :: name_
            !! name of the variable
         character(len = *), intent(in)       :: value_
            !! value of the variable
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         character(len = :), allocatable  :: add_prefix_
         !---------------------------------------------------------------------!
         add_prefix_ = 'Incorrect value encountered:  '//trim(name_)//' = '//  &
                       trim(value_)
         call write_error(add_prefix_, unit_)
         !---------------------------------------------------------------------!
      end subroutine incorrect_value_ch
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine incorrect_value_int32(name_, value_, unit_)
         !! ``incorrect value encountered:
         !!   variable_name = variable_value``
         !---------------------------------------------------------------------!
         character(len = *), intent(in)       :: name_
            !! name of the variable
         integer(int32), intent(in)           :: value_
            !! value of the variable
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         character(len = 20)              :: tmp_
         character(len = :),  allocatable :: add_prefix_
         !---------------------------------------------------------------------!
         write(tmp_, '(i5)') value_
         add_prefix_ = 'Incorrect value encountered:  '//trim(name_)//' = '//  &
                       trim(tmp_)
         call write_error(add_prefix_, unit_)
         !---------------------------------------------------------------------!
      end subroutine incorrect_value_int32
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine incorrect_value_sp(name_, value_, unit_)
         !! ``incorrect value encountered:
         !!   variable_name = variable_value``
         !---------------------------------------------------------------------!
         character(len = *), intent(in)       :: name_
            !! name of the variable
         real(sp), intent(in)                 :: value_
            !! value of the variable
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         character(len = 16), allocatable :: tmp_
         character(len = :), allocatable  :: add_prefix_
         !---------------------------------------------------------------------!
         write(tmp_, '(e16.8)') value_

         add_prefix_ = 'Incorrect value encountered:  '//trim(name_)//' = '//  &
                       trim(tmp_)
         call write_error(add_prefix_, unit_)
         !---------------------------------------------------------------------!
      end subroutine incorrect_value_sp
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine incorrect_value_dp(name_, value_, unit_)
         !! ``incorrect value encountered:
         !!   variable_name = variable_value``
         !---------------------------------------------------------------------!
         character(len = *), intent(in)       :: name_
            !! name of the variable
         real(dp), intent(in)                 :: value_
            !! value of the variable
         integer(int32), optional, intent(in) :: unit_
            !! optional, unit where the message will be written
         !---------------------------------------------------------------------!
         character(len = 16), allocatable :: tmp_
         character(len = :), allocatable  :: add_prefix_
         !---------------------------------------------------------------------!
         write(tmp_, '(e16.8)') value_

         add_prefix_ = 'Incorrect value encountered:  '//trim(name_)//' = '//  &
                       trim(tmp_)
         call write_error(add_prefix_, unit_)
         !---------------------------------------------------------------------!
      end subroutine incorrect_value_dp
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      function to_lowercase(str) result(low_str)
         !! forces lowercase on given string
         !---------------------------------------------------------------------!      
         character(len=*), intent(in) :: str
            !! input string
         character(len=len(str))      :: low_str
            !! output (lowercase) string
         !---------------------------------------------------------------------!            
         integer(int32)               :: i
         !---------------------------------------------------------------------!
         do i = 1, len(str)
            low_str(i:i) = char_to_lowercase(str(i:i))
         enddo
         !---------------------------------------------------------------------!
      end function to_lowercase
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      function char_to_lowercase(s) result(l_s)
         !! forces lowercase on a single character
         !---------------------------------------------------------------------!
         character(len=1), intent(in) :: s
            !! input character
         character(len=1)             :: l_s
            !! output (lowercase) character
         !---------------------------------------------------------------------!
         integer(int32)               :: indx
         !---------------------------------------------------------------------!
         indx = index(uppercase, s)

         if (indx > 0) then
            l_s = lowercase(indx:indx)
         else
            l_s = s
         endif
         !---------------------------------------------------------------------!
      end function char_to_lowercase
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      function integer_to_character(i, format_string) result(res)
         !! transfers integer to a character
         !---------------------------------------------------------------------!
         integer, intent(in) :: i
            !! input integer
         character(len=*), intent(in), optional :: format_string
            !! Optional format string.
         character(len=32) :: res
            !! output character
         !---------------------------------------------------------------------!
         character(len=32) :: default_format, user_format
         !---------------------------------------------------------------------!
         ! Deafult format
         !---------------------------------------------------------------------!
         default_format =  '(i0)'
         !---------------------------------------------------------------------!
         if (present(format_string)) then
            user_format = trim(format_string)
         else
            user_format = default_format
         endif
         !---------------------------------------------------------------------!
         write (res, user_format) i
         res = adjustl(res)
         !---------------------------------------------------------------------!
      end function integer_to_character
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      function float_to_character(f, format_string) result(res)
         !! Converts a floating-point number to a character string.
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: f
            !! input floating-point number
         character(len=*), intent(in), optional :: format_string
            !! Optional format string.
         character(len=64) :: res
            !! Output character string.
         !---------------------------------------------------------------------!
         character(len=32) :: default_format, user_format
         !---------------------------------------------------------------------!
         ! Default format: 6 decimal places
         !---------------------------------------------------------------------!
         default_format = '(F0.6)'  
         !---------------------------------------------------------------------!
         if (present(format_string)) then
            user_format = trim(format_string)
         else
            user_format = default_format
         endif
         !---------------------------------------------------------------------!
         write(res, user_format) f
         res = adjustl(res)
         !---------------------------------------------------------------------!
      end function float_to_character
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
end module utility_functions_mod
